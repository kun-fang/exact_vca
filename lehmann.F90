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
		Tep=cluster%hop%Tep
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

	function lehmann_particle_density(cluster,Q) result(x)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(in)::Q
		real(kind=8)::kx,ky,sum,Tep,x,w,Pi
		type(Q_type),pointer::p,a
		integer::i,j,k
		x=0.0
		if(.not.associated(Q).or..not.associated(cluster)) return
		Tep=cluster%hop%Tep
		Pi=asin(1.d0)*2
		do i=0,nk-1
			do j=0,nk-1
				kx=i*Pi/nk
				ky=j*Pi/nk
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,kx,ky)
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

	function lehmann_real_green(cluster,Q,alpha,xt,yt,time) result(g)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(in)::Q
		real(kind=8),intent(in)::time
		integer::alpha,xt,yt
		real(kind=8)::kx,ky,sum,Tep,w,Pi
		complex(8)::s,g,x
		type(Q_type),pointer::p,a
		integer::i,j,k
		x=0.0
		if(.not.associated(Q).or..not.associated(cluster)) return
		Tep=cluster%hop%Tep
		Pi=asin(1.d0)*2
		s=Zero
		do i=0,nk-1
			do j=0,nk-1
				kx=i*Pi/nk
				ky=j*Pi/nk
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,kx,ky)
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


	function lehmann_potthoff_functional(cluster,Q,omega) result(func)
		implicit none
		real(kind=8),intent(in)::omega
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(in)::Q
		real(kind=8)::func,kx,ky,sum,o1,o2,Tep,Pi
		integer::i,j,k,nw,orbit
		real(kind=8),pointer,dimension(:)::wprime,w
		real(kind=8),pointer,dimension(:,:)::tprime
		complex(kind=8),pointer,dimension(:,:)::t,V,M
		if(.not.associated(Q).or..not.associated(cluster)) return
		Tep=cluster%hop%Tep
		orbit=cluster_get_n_orbit(cluster)
		Pi=asin(1.d0)*2
		allocate(t(orbit,orbit))
		allocate(tprime(orbit,orbit))
		allocate(V(orbit,orbit))
		call hop_cluster(cluster%hop,tprime,1)
		call get_w(Q,nw,wprime)
		o1=Trln(nw,wprime,Tep)
		allocate(w(nw))
		o2=0.0
		do i=0,nk-1
			do j=0,nk-1
				kx=i*Pi/nk
				ky=j*Pi/nk
				call hop_lattice(cluster%hop,t,kx,ky,1)
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

	subroutine lehmann_transform_Q(cluster,Q,kx,ky)
		implicit none
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		type(Q_type),pointer,intent(inout)::Q
		real(kind=8),intent(in)::kx,ky
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
		call hop_cluster(cluster%hop,tprime,1)
		call hop_lattice(cluster%hop,t,kx,ky,1)
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


	subroutine hop_cluster(cl,tprime,spin)
		implicit none
		type(hop),pointer::cl
		real(8),pointer,dimension(:,:)::tprime
		integer::site,i,j,k,l,spin
		site=cl%site
		allocate(tprime(site,site))
		tprime(1:site,1:site)=0.0
		do i=1,cl%nct
			j=cl%cluster(1,i)
			k=cl%cluster(2,i)
			l=cl%cluster(3,i)
			tprime(j,k)=tprime(j,k)-cl%ct(l)
			tprime(k,j)=tprime(j,k)
		end do
		do i=1,site
			tprime(i,i)=tprime(i,i)-cl%cmu-cl%M*(-1)**(i+spin)
		end do
	end subroutine

	subroutine hop_lattice(cl,lattice,kx,ky,spin)
		implicit none
		type(hop),pointer::cl
		complex(8),pointer,dimension(:,:)::lattice
		real(8),dimension(2)::q,kk
		real(8)::kx,ky
		integer::site,i,j,k,l,spin
		site=cl%site
		kk(1)=kx
		kk(2)=ky
		allocate(lattice(site,site))
		lattice(1:site,1:site)=Zero
		do i=1,cl%nlt
			j=cl%lattice(1,i)
			k=cl%lattice(2,i)
			l=cl%lattice(3,i)
			q=cl%lattice(5,i)*cl%vector(1:2,1)+cl%lattice(6,i)*cl%vector(1:2,2)
			lattice(j,k)=lattice(j,k)-cl%lt(l)*exp(-Xi*dot_product(q,kk))
			lattice(k,j)=conjg(lattice(j,k))
		end do
		do i=1,site
			lattice(i,i)=lattice(i,i)-cl%lmu
		end do
	end subroutine





end	module lehmann
