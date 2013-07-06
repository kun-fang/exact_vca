module hubbard_cluster
	use hop_mod
	use basis
	use diagnal
	use sparse_real
	implicit none
	
	type hubbard_cluster_type
		character(len=20)::name
		type(cluster),pointer::hop
		real(kind=8)::U,miu,M,Tep,tx,ty,delta
		integer,pointer,dimension(:,:)::struct
		integer::orbit,dim
		type(basis_type),pointer::bas
		type(sparse_real_type),pointer::nakedH
	end type
	
	public::cluster_build,cluster_clean,cluster_solver
	public::cluster_get_n_orbit,cluster_get_dim,cluster_show_basis,cluster_dress

	private::naked_H_build,buildH,LowCase
	!***********************!
	!   cluster structure   !
	!***********************!
	private::square,square_aniso,betts8,ladder8,ladder6,site10
	
		
	contains
!--------------------------public---------------------------------

	function cluster_build(w,miu) result(p)
		implicit none
		type(hubbard_cluster_type),pointer::p
		class(cluster),pointer::cl
		character(len=20)::rc
		integer::i,j,k,it
		real(kind=8)::w,miu
		allocate(p)
		open(unit=7,file='cluster.input',status='old',action='read')
		read(7,*) p%name
		read(7,*) rc,miu
		read(7,*) rc,p%U
		read(7,*) rc,p%tx
		read(7,*) rc,p%ty
		read(7,*) rc,p%miu
		read(7,*) rc,p%Tep
		read(7,*) rc,p%M
		read(7,*) rc,p%delta
		read(7,*) rc,w
		close(7)
		rc=trim(adjustl(p%name))//'.conf'
		allocate(p%hop)
		cl=>p%hop
		call cluster_init(cl,trim(adjustl(rc)))
		p%orbit=p%hop%site
		p%bas=>basis_init(p%orbit)
		p%dim=basis_get_n_basis(p%bas)
		p%nakedH=>naked_H_build(p%bas,p%hop)
	end function

	subroutine cluster_clean(clust)
		implicit none
		type(hubbard_cluster_type),pointer,intent(inout)::clust
		integer::i
		if(.not.associated(clust)) return
		if(associated(clust%nakedH)) call sparse_real_clean(clust%nakedH)
		if(associated(clust%bas)) call basis_clean(clust%bas)
		deallocate(clust%hop%coordinate)
		deallocate(clust%hop%vector)
		deallocate(clust%hop%hop)
		deallocate(clust%hop)
		deallocate(clust)
	end subroutine

	subroutine cluster_solver(cluster,e,X)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(kind=8),pointer,dimension(:),intent(out)::e
		real(kind=8),pointer,dimension(:,:),intent(out)::X
		integer,pointer,dimension(:)::bk
		real(kind=8),pointer,dimension(:)::d
		real(kind=8),pointer,dimension(:,:)::H,A
		integer::n,block,i,j,k,ai,af,m
		n=cluster_get_dim(cluster)
		block=basis_get_n_block(cluster%bas)
		allocate(e(n))
		allocate(X(n,n))
		allocate(bk(block+1))
		X(1:n,1:n)=0.0
		do i=1,block
			bk(i)=basis_get_block(i,cluster%bas)
		end do
		bk(block+1)=n+1
		do i=1,block
			ai=bk(i)
			af=bk(i+1)-1
			m=af-ai+1
			H=>cluster_dress(cluster,ai,af,cluster%bas)
			d=>e(ai:af)
			A=>X(ai:af,ai:af)
			call diagnal_exact_real(m,H,d,A)
			deallocate(H)
		end do
		deallocate(bk)
	end subroutine

	function cluster_dress(p,ai,af,bas) result(pt)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::p
		integer,intent(in)::ai,af
		type(basis_type),pointer,intent(in)::bas
		real(kind=8),pointer,dimension(:,:)::pt
		integer,pointer,dimension(:,:)::b
		real(kind=8)::sum,q,x
		integer::i,j,k,n,l,ti,tj,step,elec,site
		pt=>sparse_real_clone_matrix(p%nakedH,ai,af)
		step=ai-1
		n=af-ai+1
		site=basis_get_n_orbit(bas)
		do i=1,n
			ti=i+step
			elec=basis_elec(ti,bas)
			k=0
			x=pt(i,i)
			if(abs(p%M)>1.d-10) then
				b=>basis_get_basis(ti,bas)
				k=0
				do j=1,site
					k=k-(-1)**j*(b(j,1)-b(j,2))
				end do
				deallocate(b)
			end if
			x=p%U*x-p%miu*elec-p%M*k
			pt(i,i)=x
			do j=1,n
				if(i==j) cycle
				if(j>i) then
					pt(i,j)=0.0
				end if
				tj=j+step
				x=pt(i,j)
				if(abs(x)<1.d-6) cycle
				if(x<0) then
					l=-1
				else
					l=1
				end if
				k=int(abs(x))
				if(k==1) x=-p%tx*l
				if(k==2) x=-p%ty*l
				pt(i,j)=x
			end do
		end do
	end function

	function cluster_get_n_orbit(cluster) result(x)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		integer::x
		if(.not.associated(cluster)) return
		x=cluster%orbit
	end function
	
	function cluster_get_dim(cluster) result(x)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		integer::x
		if(.not.associated(cluster)) return
		x=cluster%dim
	end function
	
	function cluster_show_basis(cluster) result(p)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(basis_type),pointer::p
		if(.not.associated(cluster)) then
			nullify(p)
			return
		end if
		p=>cluster%bas
	end function

!--------------------------private--------------------------------

	function naked_H_build(b,cl) result(H)
		implicit none
		type(basis_type),pointer,intent(in)::b
		type(cluster),pointer::cl
		integer,pointer,dimension(:,:)::hop
		type(sparse_real_type),pointer::H
		integer::i,j,orbit,dim
		orbit=basis_get_n_orbit(b)
		dim=basis_get_n_basis(b)
		H=>sparse_real_init(dim)
		allocate(hop(orbit,orbit))
		hop=cl%hop(1:orbit,1:orbit,0)
		call buildH(cl,dim,orbit,H,hop,b)
		deallocate(hop)
	end function

	subroutine buildH(cl,dim,site,H,t,b)
		implicit none
		type(cluster),pointer::cl
		integer,intent(in)::dim,site
		type(sparse_real_type),pointer,intent(inout)::H
		integer,pointer,dimension(:,:),intent(in)::t
		type(basis_type),pointer,intent(in)::b
		real(kind=8),dimension(2)::d
		integer::i,j,k1,k,l,r,s,alpha,beta,sp,init,final
		integer::n_block
		real(kind=8)::x
		
		n_block=basis_get_n_block(b)
		do l=1,n_block
			init=basis_get_block(l,b)
			if(l==n_block) then
				final=dim
			else
				final=basis_get_block(l+1,b)-1
			end if
			do r=init,final
				x=basis_double_occupy(r,b)*1.0
				call sparse_real_add_element(H,r,r,x)
				do sp=1,2
					do alpha=1,site
						i=basis_c(b,alpha,sp,r)
						if(i==0) cycle
						k1=i/abs(i)
						i=abs(i)
						do beta=1,site
							if(t(alpha,beta)==0) cycle
							s=basis_cplus(b,beta,sp,i)
							if(s==0) cycle
							k=k1*s/abs(s)
							s=abs(s)
							if(t(alpha,beta)==1) then
								d=cl%coordinate(1:2,alpha)-cl%coordinate(1:2,beta)
								if(abs(d(1))>abs(d(2))) then
									x=k*1.0
								else
									x=k*2.0
								end if
							end if
							call sparse_real_add_element(H,s,r,x)
						end do
					end do
				end do
			end do
		end do
	end subroutine

	function LowCase(a) result(b)
		implicit none
		character(len=20)::a,b
		integer::i,j,l
		b=a
		do i=1,20
			j=iachar(a(i:i))
			if(j>=65.and.j<=90) then
				b(i:i)=achar(j+32)
			end if
		end do
	end function

	!******************* cluster structure *********************


	subroutine square(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=4
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,3)=2
		hop(1,4)=1
		hop(2,1)=1
		hop(2,3)=1
		hop(2,4)=2
		hop(3,1)=2
		hop(3,2)=1
		hop(3,4)=1
		hop(4,1)=1
		hop(4,2)=2
		hop(4,3)=1
	end subroutine

	subroutine square_aniso(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=4
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,4)=2
		hop(2,1)=1
		hop(2,3)=2
		hop(3,2)=2
		hop(3,4)=1
		hop(4,1)=2
		hop(4,3)=1
	end subroutine
	
	subroutine betts8(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=8
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,4)=1
		hop(1,6)=1
		hop(1,8)=1
		hop(1,3)=2
		hop(1,7)=2
		hop(2,1)=1
		hop(2,3)=1
		hop(2,5)=1
		hop(2,7)=1
		hop(2,4)=2
		hop(2,8)=2
		hop(3,2)=1
		hop(3,4)=1
		hop(3,6)=1
		hop(3,8)=1
		hop(3,1)=2
		hop(3,5)=2
		hop(4,1)=1
		hop(4,3)=1
		hop(4,5)=1
		hop(4,7)=1
		hop(4,2)=2
		hop(4,6)=2
		hop(5,2)=1
		hop(5,4)=1
		hop(5,6)=1
		hop(5,8)=1
		hop(5,3)=2
		hop(5,7)=2
		hop(6,1)=1
		hop(6,3)=1
		hop(6,5)=1
		hop(6,7)=1
		hop(6,4)=2
		hop(6,8)=2
		hop(7,2)=1
		hop(7,4)=1
		hop(7,6)=1
		hop(7,8)=1
		hop(7,1)=2
		hop(7,5)=2
		hop(8,1)=1
		hop(8,3)=1
		hop(8,5)=1
		hop(8,7)=1
		hop(8,2)=2
		hop(8,6)=2
	end subroutine
	
	subroutine ladder8(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=8
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,4)=1
		hop(1,6)=1
		hop(1,8)=1
		hop(1,3)=2
		hop(1,7)=2
		hop(2,1)=1
		hop(2,3)=1
		hop(2,5)=1
		hop(2,7)=1
		hop(2,4)=2
		hop(2,8)=2
		hop(3,2)=1
		hop(3,4)=1
		hop(3,6)=1
		hop(3,8)=1
		hop(3,1)=2
		hop(3,5)=2
		hop(4,1)=1
		hop(4,3)=1
		hop(4,5)=1
		hop(4,7)=1
		hop(4,2)=2
		hop(4,6)=2
		hop(5,2)=1
		hop(5,4)=1
		hop(5,6)=1
		hop(5,8)=1
		hop(5,3)=2
		hop(5,7)=2
		hop(6,1)=1
		hop(6,3)=1
		hop(6,5)=1
		hop(6,7)=1
		hop(6,4)=2
		hop(6,8)=2
		hop(7,2)=1
		hop(7,4)=1
		hop(7,6)=1
		hop(7,8)=1
		hop(7,1)=2
		hop(7,5)=2
		hop(8,1)=1
		hop(8,3)=1
		hop(8,5)=1
		hop(8,7)=1
		hop(8,2)=2
		hop(8,6)=2
	end subroutine
	
	subroutine ladder6(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=6
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,4)=1
		hop(1,6)=1
		hop(1,8)=1
		hop(1,3)=2
		hop(1,7)=2
		hop(2,1)=1
		hop(2,3)=1
		hop(2,5)=1
		hop(2,7)=1
		hop(2,4)=2
		hop(2,8)=2
		hop(3,2)=1
		hop(3,4)=1
		hop(3,6)=1
		hop(3,8)=1
		hop(3,1)=2
		hop(3,5)=2
		hop(4,1)=1
		hop(4,3)=1
		hop(4,5)=1
		hop(4,7)=1
		hop(4,2)=2
		hop(4,6)=2
		hop(5,2)=1
		hop(5,4)=1
		hop(5,6)=1
		hop(5,8)=1
		hop(5,3)=2
		hop(5,7)=2
		hop(6,1)=1
		hop(6,3)=1
		hop(6,5)=1
		hop(6,7)=1
		hop(6,4)=2
		hop(6,8)=2
	end subroutine
	
	subroutine site10(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=10
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,4)=1
		hop(1,8)=1
		hop(1,10)=1
		hop(1,3)=2
		hop(1,5)=2
		hop(1,7)=2
		hop(1,9)=2
		hop(2,1)=1
		hop(2,3)=1
		hop(2,5)=1
		hop(2,9)=1
		hop(2,4)=2
		hop(2,6)=2
		hop(2,8)=2
		hop(2,10)=2
		hop(3,2)=1
		hop(3,4)=1
		hop(3,6)=1
		hop(3,10)=1
		hop(3,1)=2
		hop(3,5)=2
		hop(3,7)=2
		hop(3,9)=2
		hop(4,1)=1
		hop(4,3)=1
		hop(4,5)=1
		hop(4,7)=1
		hop(4,2)=2
		hop(4,6)=2
		hop(4,8)=2
		hop(4,10)=2
		hop(5,2)=1
		hop(5,4)=1
		hop(5,6)=1
		hop(5,8)=1
		hop(5,1)=2
		hop(5,3)=2
		hop(5,7)=2
		hop(5,9)=2
		hop(6,3)=1
		hop(6,5)=1
		hop(6,7)=1
		hop(6,9)=1
		hop(6,2)=2
		hop(6,4)=2
		hop(6,8)=2
		hop(6,10)=2
		hop(7,4)=1
		hop(7,6)=1
		hop(7,8)=1
		hop(7,10)=1
		hop(7,1)=2
		hop(7,3)=2
		hop(7,5)=2
		hop(7,9)=2
		hop(8,1)=1
		hop(8,5)=1
		hop(8,7)=1
		hop(8,9)=1
		hop(8,2)=2
		hop(8,4)=2
		hop(8,6)=2
		hop(8,10)=2
		hop(9,2)=1
		hop(9,6)=1
		hop(9,8)=1
		hop(9,10)=1
		hop(9,1)=2
		hop(9,3)=2
		hop(9,5)=2
		hop(9,7)=2
		hop(10,1)=1
		hop(10,3)=1
		hop(10,7)=1
		hop(10,9)=1
		hop(10,2)=2
		hop(10,4)=2
		hop(10,6)=2
		hop(10,8)=2
	end subroutine

	subroutine site12(hop,site)
		implicit none
		integer,intent(out)::site
		integer,pointer,dimension(:,:),intent(out)::hop
		integer::i,j
		site=12
		allocate(hop(site,site))
		do i=1,site
			do j=1,site
				hop(j,i)=0
			end do
		end do
		hop(1,2)=1
		hop(1,4)=1
		hop(1,5)=1
		hop(1,9)=1
		hop(1,6)=2
		hop(1,8)=2
		hop(1,10)=2
		hop(1,12)=2
		hop(2,1)=1
		hop(2,3)=1
		hop(2,6)=1
		hop(2,10)=1
		hop(2,5)=2
		hop(2,7)=2
		hop(2,9)=2
		hop(2,11)=2
		hop(3,2)=1
		hop(3,4)=1
		hop(3,7)=1
		hop(3,11)=1
		hop(3,6)=2
		hop(3,8)=2
		hop(3,10)=2
		hop(3,12)=2
		hop(4,1)=1
		hop(4,3)=1
		hop(4,8)=1
		hop(4,12)=1
		hop(4,5)=2
		hop(4,7)=2
		hop(4,9)=2
		hop(4,11)=2
		hop(5,1)=1
		hop(5,6)=1
		hop(5,8)=1
		hop(5,9)=1
		hop(5,2)=2
		hop(5,4)=2
		hop(5,10)=2
		hop(5,12)=2
		hop(6,2)=1
		hop(6,5)=1
		hop(6,7)=1
		hop(6,10)=1
		hop(6,1)=2
		hop(6,3)=2
		hop(6,9)=2
		hop(6,11)=2
		hop(7,3)=1
		hop(7,6)=1
		hop(7,8)=1
		hop(7,11)=1
		hop(7,3)=2
		hop(7,6)=2
		hop(7,8)=2
		hop(7,11)=2
		hop(8,4)=1
		hop(8,5)=1
		hop(8,7)=1
		hop(8,12)=1
		hop(8,1)=2
		hop(8,3)=2
		hop(8,9)=2
		hop(8,11)=2
		hop(9,1)=1
		hop(9,5)=1
		hop(9,10)=1
		hop(9,12)=1
		hop(9,2)=2
		hop(9,4)=2
		hop(9,6)=2
		hop(9,8)=2
		hop(10,2)=1
		hop(10,6)=1
		hop(10,9)=1
		hop(10,11)=1
		hop(10,1)=2
		hop(10,3)=2
		hop(10,5)=2
		hop(10,7)=2
		hop(11,3)=1
		hop(11,7)=1
		hop(11,10)=1
		hop(11,12)=1
		hop(11,2)=2
		hop(11,4)=2
		hop(11,6)=2
		hop(11,8)=2
		hop(12,4)=1
		hop(12,8)=1
		hop(12,9)=1
		hop(12,11)=1
		hop(12,1)=2
		hop(12,3)=2
		hop(12,5)=2
		hop(12,7)=2
	end subroutine


end module hubbard_cluster
