module lehmann
	use hubbard_cluster
	implicit none
	
	complex(kind=8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	!real(kind=8),private,parameter::Pi=3.1415926
	integer,private,parameter::nk=10
	
	type Q_type
		real(kind=8)::w
		complex(kind=8),pointer,dimension(:)::q
		type(Q_type),pointer::next
	end type Q_type
	
		real(kind=8),pointer,dimension(:),private::Z,e
		real(kind=8),pointer,dimension(:,:),private::X
	public::lehmann_init,lehmann_Q_clean,lehmann_potthoff_functional,lehmann_green_function,hop_cluster,hop_lattice

	private::Q_matrix,newnode,findblock,partition,Trln,get_M,get_w,c_cplus

	contains

!--------------------public--------------------

	function lehmann_init(cluster,omega,spin) result(Q)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(kind=8),intent(out)::omega
		type(Q_type),pointer::Q
		integer::i,j,k,n,spin
		real(kind=8)::Tep
		Tep=cluster%Tep
		n=cluster_get_dim(cluster)
		call cluster_solver(cluster,e,X)
		Z=>partition(n,e,Tep,omega)
		Q=>Q_Matrix(cluster,Tep,spin)
		deallocate(e)
		deallocate(X)
		deallocate(Z)
	end function

	subroutine lehmann_Q_clean(Q)
		implicit none
		type(Q_type),pointer,intent(inout)::Q
		type(Q_type),pointer::p,d
		if(.not.associated(Q)) return
		p=>Q
		do
			if(.not.associated(p)) exit
			d=>p
			p=>p%next
			deallocate(d%q)
			deallocate(d)
		end do
		nullify(Q)
	end subroutine

	function lehmann_restore(cluster,G,kx,ky) result(x)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		complex(kind=8),pointer,dimension(:,:),intent(in)::G
		real(kind=8),intent(in)::kx,ky
		complex(kind=8)::x,q,p,e
		integer::i,j,k,n,dim,orbit
		real(kind=8),allocatable,dimension(:,:)::coor
		dim=2
		orbit=cluster_get_n_orbit(cluster)
		allocate(coor(dim,orbit))
		coor=cluster%hop%coordinate
!		call square_coor(dim,orbit,coor)
!		call honey_coor(dim,orbit,coor)
		x=Zero
		if(.not.associated(cluster).or..not.associated(G)) return
		n=cluster_get_n_orbit(cluster)
		do i=1,orbit
			do j=1,orbit
				q=Xi*((coor(1,j)-coor(1,i))*kx+(coor(2,j)-coor(2,i))*ky)
				e=exp(q)
				p=e*G(i,j)
				x=x+p
			end do
		end do
		x=x/orbit
	end function

	function lehmann_particle_density(cluster,Q,miu) result(x)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(in)::Q
		real(kind=8),intent(in)::miu
		real(kind=8)::kx,ky,sum,Tep,x,w,Pi
		type(Q_type),pointer::p,a
		integer::i,j,k
		x=0.0
		if(.not.associated(Q).or..not.associated(cluster)) return
		Tep=cluster%Tep
		Pi=asin(1.d0)*2
		do i=0,nk-1
			do j=0,nk-1
				kx=i*Pi/nk
				ky=j*Pi/nk
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				a=>p
				do
					if(.not.associated(a)) exit
					w=a%w
					sum=dot_product(a%q,a%q)
					if(abs(Tep)<1.d-8) then
						if(w<0) x=x+2*sum
					else
						if(w<0) x=x+2*sum/(exp(w/Tep)+1)
						if(w>0) x=x+2*sum*exp(-w/Tep)/(1+exp(-w/Tep))
					end if
					a=>a%next
				end do
				call lehmann_Q_clean(p)
			end do
		end do
		x=x/nk/nk
	end function

	function lehmann_real_green(cluster,Q,miu,alpha,xt,yt,time) result(g)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(in)::Q
		real(kind=8),intent(in)::miu,time
		integer::alpha,xt,yt
		real(kind=8)::kx,ky,sum,Tep,w,Pi
		complex(8)::s,g,x
		type(Q_type),pointer::p,a
		integer::i,j,k
		x=0.0
		if(.not.associated(Q).or..not.associated(cluster)) return
		Tep=cluster%Tep
		Pi=asin(1.d0)*2
		s=Zero
		do i=0,nk-1
			do j=0,nk-1
				kx=i*Pi/nk
				ky=j*Pi/nk
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				a=>p
				x=Zero
				do
					if(.not.associated(a)) exit
					w=a%w
					sum=conjg(a%q(1))*a%q(alpha)
					if(abs(Tep)<1.d-6) then
						if(w<0) x=x+2*sum*exp(Xi*w*time)
					else
						if(w<0) x=x+2*sum*exp(Xi*w*time)/(exp(w/Tep)+1)
						if(w>0) x=x+2*sum*exp(-w/Tep)*exp(Xi*w*time)/(1+exp(-w/Tep))
					end if
					a=>a%next
				end do
				s=s+x/nk/nk*exp(Xi*(kx*xt+ky*yt))
				call lehmann_Q_clean(p)
			end do
		end do
		g=s
	end function


	function lehmann_potthoff_functional(cluster,Q,omega,miu) result(func)
		implicit none
		real(kind=8),intent(in)::omega
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(in)::Q
		real(kind=8),intent(in)::miu
		real(kind=8)::func,kx,ky,sum,o1,o2,Tep,Pi
		integer::i,j,k,nw,orbit
		real(kind=8),pointer,dimension(:)::wprime,w
		real(kind=8),pointer,dimension(:,:)::tprime
		complex(kind=8),pointer,dimension(:,:)::t,V,M
		if(.not.associated(Q).or..not.associated(cluster)) return
		Tep=cluster%Tep
		orbit=cluster_get_n_orbit(cluster)
		Pi=asin(1.d0)*2
		allocate(t(orbit,orbit))
		allocate(tprime(orbit,orbit))
		allocate(V(orbit,orbit))
		call hop_cluster(cluster,tprime,1)
		call get_w(Q,nw,wprime)
		o1=Trln(nw,wprime,Tep)
		allocate(w(nw))
		o2=0.0
		do i=0,nk-1
			do j=0,nk-1
				kx=i*Pi/nk
				ky=j*Pi/nk
				call hop_lattice(cluster,miu,t,kx,ky,1)
				V=t-tprime
				M=>get_M(Q,V,nw,orbit,wprime)
				call diagnal_exact_complex(nw,M,w,M)
				o2=o2+Trln(nw,w,Tep)
				deallocate(M)
			end do
		end do
		func=omega-2*o1+2*o2/nk/nk
		deallocate(w)
		deallocate(wprime)
		deallocate(V)
		deallocate(t)
		deallocate(tprime)
	end function

	subroutine lehmann_transform_Q(cluster,Q,miu,kx,ky)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(inout)::Q
		real(kind=8),intent(in)::miu,kx,ky
		real(kind=8),pointer,dimension(:)::w,wprime
		real(kind=8),pointer,dimension(:,:)::tprime
		complex(kind=8),pointer,dimension(:,:)::t,V,M
		complex(kind=8),pointer,dimension(:,:)::tran
		complex(kind=8)::z
		type(Q_type),pointer::a
		integer::i,j,k,nw,orbit
		orbit=cluster_get_n_orbit(cluster)
		allocate(tprime(orbit,orbit))
		allocate(t(orbit,orbit))
		allocate(V(orbit,orbit))
		call hop_cluster(cluster,tprime,1)
		call hop_lattice(cluster,miu,t,kx,ky,1)
		V=t-tprime
		call get_w(Q,nw,wprime)
		allocate(w(nw))
		M=>get_M(Q,V,nw,orbit,wprime)
		call diagnal_exact_complex(nw,M,w,M)
		allocate(tran(nw,orbit))
		tran(1:nw,1:orbit)=0.0
		a=>Q
		do k=1,nw
			do j=1,orbit
				do i=1,nw
					tran(i,j)=tran(i,j)+(M(k,i)*a%q(j))!conjg
				end do
			end do
			a=>a%next
		end do
		a=>Q
		do i=1,nw
			a%w=w(i)
			a%q=tran(i,1:orbit)
			a=>a%next
		end do
		deallocate(tran)
		deallocate(M)
		deallocate(w)
		deallocate(wprime)
		deallocate(V)
		deallocate(t)
		deallocate(tprime)
	end subroutine

	function lehmann_clone_Q(Q) result(p)
		implicit none
		type(Q_type),pointer,intent(in)::Q
		type(Q_type),pointer::p,a,b,node
		complex(kind=8),pointer,dimension(:)::y
		integer::n
		nullify(p)
		if(.not.associated(Q))return
		a=>Q
		n=sizeof(a%q)/sizeof(a%q(1))
		allocate(y(n))
		y(1:n)=a%q(1:n)
		node=>newnode(a%w,y)
		p=>node
		b=>p
		do
			a=>a%next
			if(.not.associated(a)) exit
			allocate(y(n))
			y(1:n)=a%q(1:n)
			node=>newnode(a%w,y)
			b%next=>node
			b=>node
		end do 
	end function

	subroutine lehmann_green_function(w,orbit,Q,G)
		implicit none
		complex(kind=8),intent(in)::w
		integer,intent(in)::orbit
		type(Q_type),pointer,intent(in)::Q
		complex(kind=8),dimension(:,:)::G
		type(Q_type),pointer::a,b
		complex(kind=8),pointer,dimension(:)::x,y
		complex(kind=8)::sum,z
		integer::i,j,k,alpha,beta,n
		i=0
		j=0
		n=0
		a=>Q
		G(1:orbit,1:orbit)=Zero
		do
			if(.not.associated(a)) exit
			i=i+1
			x=>a%q
			do alpha=1,orbit
				do beta=1,orbit
					G(alpha,beta)=G(alpha,beta)+conjg(x(alpha))*x(beta)/(w-a%w)
				end do
			end do
		a=>a%next
		end do
	end subroutine

!--------------------------private----------------------------

	function Q_Matrix(cluster,Tep,spin) result(Q)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(kind=8),intent(in)::Tep
		type(Q_type),pointer::Q
		type(basis_type),pointer::bas
		type(Q_type),pointer::node,d
		real(kind=8),pointer,dimension(:)::a,b
		complex(kind=8),pointer,dimension(:)::test,test2
		real(kind=8),pointer,dimension(:,:)::c
		integer,pointer,dimension(:)::block
		integer::i,j,k,n,m,nb,orbit,plus,alpha,ai,af,xi,xf,spin
		real(kind=8)::sum,w,apt
		bas=>cluster_show_basis(cluster)
		n=basis_get_n_basis(bas)
		orbit=basis_get_n_orbit(bas)
		nb=basis_get_n_block(bas)
		allocate(block(nb+1))
		do i=1,nb
			block(i)=basis_get_block(i,bas)
		end do
		block(nb+1)=n+1
		m=0
		do plus=0,1
			do i=1,n
				if(Z(i)<1.d-4) cycle
				allocate(c(n,orbit))
				c(1:n,1:orbit)=0.0
				ai=n
				af=1
				xi=n
				xf=1
				do alpha=1,orbit
					a=>X(1:n,i)
					b=>c(1:n,alpha)
					call c_cplus(alpha,spin,a,b,bas,plus)
					call findblock(b,xi,xf,nb,block)
					if(xi<ai) ai=xi
					if(xf>af) af=xf
				end do
				do j=ai,af
					w=(e(j)-e(i))*(-1)**(plus+1)
					allocate(test(orbit))
					sum=0.0
					do alpha=1,orbit
						test(alpha)=dot_product(c(ai:af,alpha),X(ai:af,j))*sqrt(Z(i))
						sum=sum+real(test(alpha)*test(alpha))
					end do
					if(sum<1.d-8) then
						deallocate(test)
						cycle
					end if
					node=>newnode(w,test)
					m=m+1
					if(m==1) then
						Q=>node
						d=>node
					else
						d%next=>node
						d=>node
					end if
				end do
				deallocate(c)
			end do
		end do
		deallocate(block)
	end function

	function newnode(w,q) result(node)
		implicit none
		real(kind=8)::w
		complex(kind=8),pointer,dimension(:)::q
		type(Q_type),pointer::node
		allocate(node)
		node%w=w
		node%q=>q
		nullify(node%next)
	end function

	subroutine findblock(y,ai,af,nb,block)
		implicit none
		integer,intent(in)::nb
		integer,intent(out)::ai,af
		real(kind=8),pointer,dimension(:),intent(in)::y
		integer,pointer,dimension(:),intent(in)::block
		integer::i,j
		c1: do i=1,nb
			do j=block(i),block(i+1)-1
				if(abs(y(j))>1.d-8) then
					ai=block(i)
					exit c1
				end if
			end do
		end do c1
		c2: do i=nb,1,-1
			do j=block(i),block(i+1)-1
				if(abs(y(j))>1.d-8) then
					af=block(i+1)-1
					exit c2
				end if
			end do
		end do c2
	end subroutine


	function partition(n,e,T,omega) result(Z)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::T
		real(kind=8),intent(out)::omega
		real(kind=8),pointer,dimension(:),intent(in)::e
		real(kind=8),pointer,dimension(:)::Z
		integer::i,j,k
		real*8::sum,w0
		allocate(Z(n))
		w0=999
		do i=1,n
			if(e(i)<w0) w0=e(i)
		end do
		sum=0
		if(abs(T)<1.d-6) then
			omega=w0
			k=0
			do i=1,n
				Z(i)=0.0
				if(abs(e(i)-w0).lt.1.d-8) then
					Z(i)=1.0
					k=k+1
				else
				end if
			end do
			do i=1,n
				Z(i)=Z(i)/k
			end do
		else
			do i=1,n
				sum=sum+exp(-(e(i)-w0)/T)
			end do
			omega=w0-T*log(sum)
			do i=1,n
				Z(i)=exp(-(e(i)-w0)/T)/sum
			end do
		end if
	end function

	function Trln(n,w,Tep) result(sum)
		implicit none
		integer,intent(in)::n
		real(kind=8),pointer,dimension(:),intent(in)::w
		real(kind=8),intent(in)::Tep
		real(kind=8)::sum
		integer::i
		sum=0.0
		do i=1,n
			if(abs(Tep)<1.d-6) then
				if(w(i)<0) sum=sum+w(i)
			else
				if(w(i)>0) sum=sum-Tep*log(1+exp(-w(i)/Tep))
				if(w(i)<0) sum=sum-Tep*log(1+exp(w(i)/Tep))+w(i)
			end if
		end do
	end function

	function get_nw(Q) result(n)
		implicit none
		type(Q_type),pointer,intent(in)::Q
		type(Q_type),pointer::a
		integer::n
		n=0
		if(.not.associated(Q)) return
		a=>Q
		do
			if(.not.associated(a)) exit
			n=n+1
			a=>a%next
		end do
	end function

	function get_M(Q,V,n,orbit,w) result(M)
		implicit none
		type(Q_type),pointer,intent(in)::Q
		complex(kind=8),pointer,dimension(:,:),intent(in)::V
		integer,intent(in)::n,orbit
		real(kind=8),pointer,dimension(:),intent(in)::w
		complex(kind=8),pointer,dimension(:,:)::M
		type(Q_type),pointer::a,b
		integer::i,j,alpha,beta
		complex(kind=8)::sum
		allocate(M(n,n))
		a=>Q
		do i=1,n
			b=>Q
			do j=1,i
					sum=Zero
					do alpha=1,orbit
						do beta=1,orbit
							sum=sum+a%q(alpha)*V(alpha,beta)*b%q(beta)
						end do
					end do
					M(i,j)=sum
					M(j,i)=conjg(sum)
				b=>b%next
			end do
			M(i,i)=M(i,i)+w(i)
			a=>a%next
		end do
	end function

	subroutine get_w(Q,n,w)
		implicit none
		type(Q_type),pointer,intent(in)::Q
		integer,intent(out)::n
		real(kind=8),pointer,dimension(:),intent(out)::w
		type(Q_type),pointer::a
		integer::i
		a=>Q
		n=0
		do
			if(.not.associated(a)) exit
			n=n+1
			a=>a%next
		end do
		allocate(w(n))
		i=0
		a=>Q
		do
			if(.not.associated(a)) exit
			i=i+1
			w(i)=a%w
			a=>a%next
		end do
	end subroutine

	subroutine c_cplus(alpha,spin,x,y,bas,plus)
		implicit none
		integer,intent(in)::alpha,spin,plus
		real(kind=8),pointer,dimension(:),intent(in)::x
		real(kind=8),pointer,dimension(:),intent(out)::y
		type(basis_type),pointer::bas
		integer::n,i,j,k
		n=basis_get_n_basis(bas)
		do i=1,n
			if(abs(x(i))<1.d-10) cycle
			if(plus==1) then
				j=basis_cplus(bas,alpha,spin,i)
			else
				j=basis_c(bas,alpha,spin,i)
			end if
			if(j==0) cycle
			k=j/abs(j)
			j=abs(j)
			y(j)=k*x(i)
		end do
	end subroutine

	subroutine hop_cluster(cluster,tprime,spin)
		implicit none
		integer,intent(in)::spin
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(kind=8),pointer,dimension(:,:),intent(inout)::tprime
		integer,allocatable,dimension(:,:)::hop
		real(kind=8),dimension(2)::d
		integer::i,j,k,n,site
		site=cluster_get_n_orbit(cluster)
		allocate(hop(site,site))
!		call honey_cluster(hop,site)
!		call square_cluster(hop,site)
		hop(1:site,1:site)=cluster%hop%hop(1:site,1:site,0)
		do i=1,site
			do j=1,site
				k=hop(i,j)
				if(k==0) tprime(i,j)=0.0
				if(k==1) then
					d=cluster%hop%coordinate(1:2,i)-cluster%hop%coordinate(1:2,j)
					if(abs(d(1))>abs(d(2))) then
						tprime(i,j)=-cluster%tx
					else
						tprime(i,j)=-cluster%ty
					end if
				end if
				if(k==2) tprime(i,j)=0.0
			end do
			tprime(i,i)=-cluster%miu-cluster%M*(-1)**(spin+i)
		end do
	end subroutine

	subroutine hop_lattice(cluster,miu,lattice,kx,ky,spin)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		integer,intent(in)::spin
		complex(kind=8),pointer,dimension(:,:),intent(inout)::lattice
		real(kind=8)::kx,ky,miu
		integer,allocatable,dimension(:)::st
		real(kind=8),allocatable,dimension(:)::p,f
		real(kind=8)::q
		integer::i,j,k,n,l,site,dim
		dim=cluster%hop%dim
		site=cluster_get_n_orbit(cluster)
		l=cluster%hop%translation
		allocate(st(dim))
		allocate(p(dim))
		allocate(f(dim))
!		call honey_lattice(lattice,site,kx,ky,miu)
!		call square_lattice(lattice,site,kx,ky,miu)
		lattice=cluster%hop%hop(1:site,1:site,0)*(-One)
		p(1)=kx
		p(2)=ky
		do i=1,8
			q=dot_product(cluster%hop%vector(1:dim,i),p)
			lattice=lattice-cluster%hop%hop(1:site,1:site,i)*exp(Xi*q)
		end do
		do i=1,site
			lattice(i,i)=lattice(i,i)-miu
		end do
	end subroutine


!===========================================================================================

	subroutine honey_coor(dim,orbit,coor)
		integer::dim,orbit
		real(kind=8),dimension(dim,orbit)::coor
		coor(1,1)=0.00
		coor(2,1)=0.00
		coor(1,2)=-0.50
		coor(2,2)=1.0*sqrt(3.0)/2
		coor(1,3)=0.00
		coor(2,3)=1.0*sqrt(3.0)
		coor(1,4)=1.00
		coor(2,4)=1.0*sqrt(3.0)
		coor(1,5)=1.50
		coor(2,5)=1.0*sqrt(3.0)/2
		coor(1,6)=1.00
		coor(2,6)=0.00
	end subroutine

	subroutine square_coor(dim,orbit,coor)
		integer::dim,orbit
		real(kind=8),dimension(dim,orbit)::coor
		coor(1,1)=0.00
		coor(2,1)=0.00
		coor(1,2)=0.00
		coor(2,2)=1.00
		coor(1,3)=1.00
		coor(2,3)=1.00
		coor(1,4)=1.00
		coor(2,4)=0.00
	end subroutine

	subroutine honey_cluster(t,n)
		implicit none
		integer,intent(in)::n
		integer,dimension(n,n),intent(inout)::t
		t(1:n,1:n)=0
		t(1,2)=1
		t(1,6)=1
		t(2,1)=1
		t(2,3)=1
		t(3,2)=1
		t(3,4)=1
		t(4,3)=1
		t(4,5)=1
		t(5,4)=1
		t(5,6)=1
		t(6,5)=1
		t(6,1)=1
	end subroutine
    
	subroutine honey_lattice(t,n,kx,ky,miu)
		implicit none
		complex(kind=8),pointer,dimension(:,:),intent(out)::t
		real(kind=8),intent(in)::miu,kx,ky
		real(kind=8),dimension(2,2)::trans
		integer::i,j,n
		trans(1,1)=3.00
		trans(2,1)=0.00
		trans(1,2)=1.50
		trans(2,2)=3*sqrt(3.0)/2
		do i=1,n
			do j=1,n
				t(i,j)=Zero
			end do
			t(i,i)=-One*miu
		end do

		t(1,2)=-One
		t(1,6)=-One
		t(1,4)=-exp(Xi*(kx*trans(1,2)+ky*trans(2,2)))
		t(2,1)=-One
		t(2,3)=-One
		t(2,5)=-exp(Xi*(kx*trans(1,1)+ky*trans(2,1)))
		t(3,2)=-One
		t(3,4)=-One
		t(3,6)=-exp(Xi*(kx*(trans(1,1)-trans(1,2))+ky*(trans(2,1)-trans(2,2))))
		t(4,3)=-One
		t(4,5)=-One
		t(4,1)=-exp(-Xi*(kx*trans(1,2)+ky*trans(2,2)))
		t(5,4)=-One
		t(5,6)=-One
		t(5,2)=-exp(-Xi*(kx*trans(1,1)+ky*trans(2,1)))
		t(6,5)=-One
		t(6,1)=-One
		t(6,3)=-exp(-Xi*(kx*(trans(1,1)-trans(1,2))+ky*(trans(2,1)-trans(2,2))))
	end subroutine

	subroutine square_lattice(t,site,kx,ky,miu)
		implicit none
		complex(kind=8),pointer,dimension(:,:),intent(out)::t
		real(kind=8),intent(in)::miu,kx,ky
		integer::i,j,site
		do i=1,site
			do j=1,site
				t(i,j)=Zero
			end do
			t(i,i)=-One*miu
		end do

		t(1,2)=-(One+exp(2*Xi*ky))
		t(1,4)=-(One+exp(2*Xi*kx))
		t(2,1)=-(One+exp(-2*Xi*ky))
		t(2,3)=-(One+exp(2*Xi*kx))
		t(3,2)=-(One+exp(-2*Xi*kx))
		t(3,4)=-(One+exp(-2*Xi*ky))
		t(4,1)=-(One+exp(-2*Xi*kx))
		t(4,3)=-(One+exp(2*Xi*ky))
	end subroutine


	subroutine square_cluster(t,n)
		implicit none
		integer,intent(in)::n
		integer,dimension(n,n),intent(inout)::t
		t(1:n,1:n)=0
		t(1,2)=1
		t(1,4)=1
		t(2,1)=1
		t(2,3)=1
		t(3,2)=1
		t(3,4)=1
		t(4,3)=1
		t(4,1)=1
	end subroutine




end	module lehmann
