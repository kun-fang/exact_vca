module hubbard_cluster
	use hop_mod
	use basis
	use diagnal
	use sparse_real
	implicit none
	
	type hubbard_cluster_type
		type(hop),pointer::hop
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
	!private::square,square_aniso,betts8,ladder8,ladder6,site10
	
		
	contains
!--------------------------public---------------------------------

	function cluster_build() result(p)
		implicit none
		type(hubbard_cluster_type),pointer::p
		allocate(p)
		p%hop=>cluster_init()
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
		if(associated(clust%hop)) call cluster_delete(clust%hop)
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
			if(abs(p%hop%M)>1.d-10) then
				b=>basis_get_basis(ti,bas)
				k=0
				do j=1,site
					k=k-(-1)**j*(b(j,1)-b(j,2))
				end do
				deallocate(b)
			end if
			x=p%hop%U*x-p%hop%cmu*elec-p%hop%M*k
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
				x=-p%hop%ct(k)*l
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
		type(hop),pointer::cl
		type(sparse_real_type),pointer::H
		integer::i,j,orbit,dim
		orbit=basis_get_n_orbit(b)
		dim=basis_get_n_basis(b)
		H=>sparse_real_init(dim)
		call buildH(cl,dim,orbit,H,b)
	end function

	subroutine buildH(cl,dim,site,H,b)
		implicit none
		type(hop),pointer::cl
		integer,intent(in)::dim,site
		type(sparse_real_type),pointer,intent(inout)::H
		integer,pointer,dimension(:,:)::t
		type(basis_type),pointer,intent(in)::b
		real(kind=8),dimension(2)::d
		integer::i,j,k1,k,l,r,s,alpha,beta,sp,init,final,ii
		integer::n_block
		real(kind=8)::x
		logical::ok
		
		n_block=basis_get_n_block(b)
		t=>cl%cluster
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
							do ii=1,cl%nct
								if((t(1,ii)==alpha.and.t(2,ii)==beta).or.(t(1,ii)==beta.and.t(2,ii)==alpha)) exit
							end do
							if(ii==cl%nct+1) cycle
							if(t(3,ii)==0) cycle
							s=basis_cplus(b,beta,sp,i)
							if(s==0) cycle
							k=k1*s/abs(s)
							s=abs(s)
							x=k*t(3,ii)
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


end module hubbard_cluster
