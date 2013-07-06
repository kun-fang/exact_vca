module Efunc
	use hubbard_cluster
	use lehmann
	use cmatrix
	use MCDOS

	type(Q_type),private,pointer::commq
	type(hubbard_cluster_type),private,pointer::commcl
	real(kind=8),private::commiu
	
	complex(kind=8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	real(kind=8),private,parameter::eps=0.015

contains

	subroutine Efunc_init(n,Emin,Emax,qmat,clust,miu)
		integer::n
		real(8)::Emin,Emax,dE
		type(Q_type),pointer::qmat
		type(hubbard_cluster_type),pointer::clust
		real(kind=8)::miu
		commq=>qmat
		commcl=>clust
		commiu=miu
		call initDOS(n,Emin,Emax)
	end subroutine
	
	function Efunc_position(d,x,hit) result(k)
		implicit none
		type(Q_type),pointer::p
		integer::orbit,i,j,d,k,n,hit(nE)
		real(kind=8)::omega,x(d),kx,ky,Pi,gl,a(3)
		complex(kind=8),pointer,dimension(:,:)::G
		if(.not.associated(commcl)) return
		Pi=asin(1.d0)*2
		n=nE
		orbit=cluster_get_n_orbit(commcl)
		allocate(G(orbit,orbit))
		p=>lehmann_clone_Q(commq)
		kx=x(1)
		ky=x(2)
		a(1:3)=0
		call lehmann_transform_Q(commcl,p,commiu,kx,ky)
		hit(1:n)=0
		k=0
		do i=1,n
			call lehmann_green_function(Ek(i)+eps*Xi,orbit,p,G)
			gl=-imag(lehmann_restore(commcl,G,kx,ky))/Pi
			a(1)=a(2)
			a(2)=a(3)
			a(3)=gl
			if(i<3) cycle
			if(a(2)>a(1).and.a(2)>a(3)) then
				hit(i-1)=1
				k=k+1
				!print *,k
			end if
		end do
		deallocate(G)
		call lehmann_Q_clean(p)
	end function

	subroutine Efunc_clean()
		deallocate(Ek)
	end subroutine

end module Efunc


module VCA
	use hubbard_cluster
	use lehmann
	use optimal
	use random
	use cmatrix
	use MCDOS
	use Efunc
	!use cmatrix
	complex(kind=8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	real(kind=8),private,parameter::eps=0.001
	integer,private,parameter::nk=20


	contains

	subroutine VCA_optimal(cluster,miu,omega,density,err)
		implicit none
		real(kind=8),intent(in)::miu
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		real(kind=8),intent(out)::omega,density
		real(kind=8),pointer,dimension(:)::c,d,var
		real(kind=8),pointer,dimension(:,:)::q
		real(kind=8)::f,r,fin,fout
		type(optimal_type),pointer::p
		type(random_type),pointer::rand
		integer::i,j,n,orbit,ie
		logical::OK,err
		err=.false.
		OK=.false.
		orbit=cluster_get_n_orbit(cluster)
		call varinit(cluster,n,var)
		fin=sqrt(dot_product(var,var))
		allocate(c(n))
		allocate(d(n))
		allocate(q(n,n))
		open(unit=9,file='optim.data',status='old',iostat=ie)
		if(ie==0) then
			do while(ie==0)
				read(9,*,iostat=ie) cluster%tx,cluster%miu,omega
			end do
			close(9)
		end if
		p=>optimal_init(n,var)
		rand=>random_init()
		open(unit=9,file='optim.data',access='append')
		do
			if(OK) exit
			call update(cluster,p%v)
			f=VCA_potthoff_functional(cluster,miu)
			fout=sqrt(dot_product(p%v,p%v))
			if(abs(fout-fin)>150) then
				err=.true.
				print *,"error!"
				exit
			end if
			write(9,*) p%v,f
			do i=1,n
				c(1:n)=p%v(1:n)
				c(i)=c(i)+p%dx
				call update(cluster,c)
				d(i)=VCA_potthoff_functional(cluster,miu)
				do j=1,n
					c(1:n)=p%v(1:n)
					c(i)=c(i)+p%dx
					c(j)=c(j)+p%dx
					call update(cluster,c)
					q(i,j)=VCA_potthoff_functional(cluster,miu)
				end do
			end do
			r=random_new(rand)
			OK=optimal_stationary(p,f,d,q,r)
		end do
		if(.not.err) then
			call update(cluster,p%v)
			omega=VCA_potthoff_functional(cluster,miu)
			density=VCA_particle_density(cluster,miu)
		end if
		call optimal_clean(p)
		call random_clean(rand)
		deallocate(var)
		deallocate(c)
		deallocate(q)
		deallocate(d)
		close(9)
	end subroutine

	function VCA_potthoff_functional(cluster,miu) result(x)
		implicit none
		real(kind=8),intent(in)::miu
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q
		real(kind=8)::omega,x
		integer::orbit
		if(.not.associated(cluster)) return
		orbit=cluster_get_n_orbit(cluster)
		Q=>lehmann_init(cluster,omega,1)
		if(.not.associated(Q)) return
		x=lehmann_potthoff_functional(cluster,Q,omega,miu)/orbit
		!write(7,'(4F20.15)') cluster%tx,cluster%ty,cluster%miu,x
		call lehmann_Q_clean(Q)
	end function

	function VCA_green(cluster,miu,w,kx,ky) result(g)
		implicit none
		type(hubbard_cluster_type),pointer::cluster
		real(kind=8)::miu,w,kx,ky,omega,pkx,pky
		integer::orbit,i,j,k,dim
		complex(kind=8)::g
		complex(kind=8),pointer,dimension(:,:)::Gp
		type(Q_type),pointer::Q
		logical::OK
		if(.not.associated(cluster)) return
		dim=2
		orbit=cluster_get_n_orbit(cluster)
		allocate(Gp(orbit,orbit))
		Q=>lehmann_init(cluster,omega,1)
		call lehmann_transform_Q(cluster,Q,miu,kx,ky)
		call lehmann_green_function(w+Xi*eps,orbit,Q,Gp)
		g=lehmann_restore(cluster,Gp,kx,ky)
		call lehmann_Q_clean(Q)
		deallocate(Gp)
	end function
	
	function VCA_real_green(cluster,miu,x,y,time,spin) result(g)
		implicit none
		type(hubbard_cluster_type),pointer::cluster
		real(kind=8)::miu,time,omega
		integer::alpha,orbit,i,j,k,x,y,xt,yt,xr,yr,s,spin
		complex(8)::g
		type(Q_type),pointer::Q
		logical::OK
		if(.not.associated(cluster)) return
		i=0
		if(x>=0) then
			s=1
		else
			s=-1
		end if
		do
			if(i*2>abs(x)) exit
			i=i+1
		end do
		xt=2*(i-(1+s)/2)*s
		xr=x-xt
		i=0
		if(y>=0) then
			s=1
		else
			s=-1
		end if
		do
			if(i*2>abs(y)) exit
			i=i+1
		end do
		yt=2*(i-(1+s)/2)*s
		yr=y-yt
		if(xr==0.and.yr==0) alpha=1
		if(xr==0.and.yr==1) alpha=2
		if(xr==1.and.yr==1) alpha=3
		if(xr==1.and.yr==0) alpha=4
		orbit=cluster_get_n_orbit(cluster)
		Q=>lehmann_init(cluster,omega,spin)
		g=lehmann_real_green(cluster,Q,miu,alpha,xt,yt,time)
		call lehmann_Q_clean(Q)
	end function

	function VCA_particle_density(cluster,miu) result(density)
		implicit none
		real(kind=8),intent(in)::miu
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q
		real(kind=8)::density,omega
		integer::orbit
		if(.not.associated(cluster)) return
		orbit=cluster_get_n_orbit(cluster)
		Q=>lehmann_init(cluster,omega,1)
		if(.not.associated(Q)) return
		density=lehmann_particle_density(cluster,Q,miu)/orbit
		call lehmann_Q_clean(Q)
	end function

	subroutine VCA_spectral_function(cluster,miu)
		implicit none
		real(kind=8),intent(in)::miu
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q,p
		real(kind=8)::omega,x,kx,ky,w,Pi
		complex(kind=8)::gl
		integer::orbit,i,j,k
		complex(kind=8),pointer,dimension(:,:)::G
		if(.not.associated(cluster)) return
		orbit=cluster_get_n_orbit(cluster)
		allocate(G(orbit,orbit))
		Q=>lehmann_init(cluster,omega,1)
		if(.not.associated(Q)) return
		open(unit=7,file='spectr.data',status='replace')
		k=0
		kx=0.d0
		ky=0.d0
		Pi=3.1415926d0

		do
			if(ky-Pi>1.d-6) then
				ky=Pi
				exit
			end if
			w=-5
			k=k+1
			write(*,*) kx,ky
			do
				if(w>5) exit
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				call lehmann_green_function(w+Xi*eps,orbit,p,G)
				gl=lehmann_restore(cluster,G,kx,ky)
				x=-imag(gl)/Pi
				write(7,*) k,w,x
				call lehmann_Q_clean(p)
				w=w+0.1
			end do
			ky=ky+Pi/nk/2.0
			!kx=0.d0
			kx=ky
		end do

		do
			if(ky<-1.d-6) then
				kx=0
				exit
			end if
			w=-5
			k=k+1
			write(*,*) kx,ky
			do
				if(w>5) exit
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				call lehmann_green_function(w+Xi*eps,orbit,p,G)
				gl=lehmann_restore(cluster,G,kx,ky)
				x=-imag(gl)/Pi
				write(7,*) k,w,x
				call lehmann_Q_clean(p)
				w=w+0.1
			end do
			!kx=kx+sqrt(3.d0)*Pi/nk/2.d0
			!ky=-kx/sqrt(3.d0)+4.d0*Pi/3/sqrt(3.d0)
			ky=ky-Pi/nk/2.0
		end do

		do
			if(kx<-1.d-6) then
				kx=0
				exit
			end if
			w=-5
			k=k+1
			write(*,*) kx,ky
			do
				if(w>5) exit
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				call lehmann_green_function(w+Xi*eps,orbit,p,G)
				gl=lehmann_restore(cluster,G,kx,ky)
				x=-imag(gl)/Pi
				write(7,*) k,w,x
				call lehmann_Q_clean(p)
				w=w+0.1
			end do
			kx=kx-Pi/nk/2.d0
			!ky=sqrt(3.d0)*kx
		end do

		call lehmann_Q_clean(Q)
		deallocate(G)
		close(7)
	end subroutine

	subroutine VCA_trans(cluster,miu,v)
		implicit none
		real(kind=8),intent(in)::miu,v
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q,p,p1,p2,a,a1,a2
		real(kind=8)::omega,x,y,kx,ky,qx,qy,coor(2,4),Pi,broad,w,wi,wo,ws
		real(kind=8),allocatable::wn(:),wn1(:),wn2(:)
		complex(kind=8)::gl,qq,tt,gg
		integer::orbit,i,j,k,ii,jj,s,r,l,nw,nw1,nw2,alpha,beta,nn
		logical::OK
		complex(kind=8),allocatable,dimension(:,:)::Gq,qn,qn1,qn2
		complex(kind=8),pointer::G(:,:)
		if(.not.associated(cluster)) return
		Pi=asin(1.d0)*2
		orbit=cluster_get_n_orbit(cluster)
		allocate(G(orbit,orbit))
		Q=>lehmann_init(cluster,omega,1)
		if(.not.associated(Q)) return
		coor=cluster%hop%coordinate
		nn=200
		gl=Zero
		do i=0,nn
			kx=i*2*Pi/nn+Pi-0.1*Pi
			ky=0.d0
			p=>lehmann_clone_Q(Q)
			call lehmann_transform_Q(cluster,p,miu,kx,ky)
			call lehmann_green_function(v+eps*Xi,orbit,p,G)
			gl=lehmann_restore(cluster,G,kx,ky)
			call lehmann_Q_clean(p)
			x=-imag(gl)/Pi
			p=>lehmann_clone_Q(Q)
			call lehmann_transform_Q(cluster,p,miu,-kx,-ky)
			call lehmann_green_function(-v+eps*Xi,orbit,p,G)
			gl=lehmann_restore(cluster,G,-kx,-ky)
			call lehmann_Q_clean(p)
			x=x-imag(gl)/Pi
			kx=Pi
			ky=i*2*Pi/nn-0.1*Pi
			p=>lehmann_clone_Q(Q)
			call lehmann_transform_Q(cluster,p,miu,kx,ky)
			call lehmann_green_function(v+eps*Xi,orbit,p,G)
			gl=lehmann_restore(cluster,G,kx,ky)
			call lehmann_Q_clean(p)
			y=-imag(gl)/Pi
			p=>lehmann_clone_Q(Q)
			call lehmann_transform_Q(cluster,p,miu,-kx,-ky)
			call lehmann_green_function(-v+eps*Xi,orbit,p,G)
			gl=lehmann_restore(cluster,G,-kx,-ky)
			call lehmann_Q_clean(p)
			y=y-imag(gl)/Pi
			write(*,*) ky,x,y
            if(i==50) exit
		end do
		deallocate(G)
	end subroutine


	function VCA_op(cluster,miu,v) result(x)
		implicit none
		real(kind=8),intent(in)::miu,v
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q1,Q2,p,p1,p2,a,a1,a2
		real(kind=8)::omega,x,kx,ky,qx,qy,coor(2,4),Pi,broad,w,wi,wo,ws
		real(kind=8),allocatable::wn(:),wn1(:),wn2(:)
		complex(kind=8)::glx,gly,qq,tt,gg
		integer::orbit,i,j,k,ii,jj,s,r,l,nw,nw1,nw2,alpha,beta,nn
		logical::OK
		complex(kind=8),allocatable,dimension(:,:)::G,Gq,qn,qn1,qn2
		if(.not.associated(cluster)) return
		Pi=asin(1.d0)*2
		broad=1.d-4
		orbit=cluster_get_n_orbit(cluster)
		allocate(G(orbit,orbit))
		allocate(Gq(orbit,orbit))
		Q1=>lehmann_init(cluster,omega,1)
		Q2=>lehmann_init(cluster,omega,2)
		if(.not.associated(Q1)) return
		if(.not.associated(Q2)) return
		coor=cluster%hop%coordinate
		wi=Pi*0.0001
		wo=5
		ws=2*Pi*0.0001
		nn=10
		glx=Zero
		do ii=0,nn
			do jj=0,nn
				kx=ii*Pi/nn
				ky=jj*Pi/nn
				p=>lehmann_clone_Q(Q1)
				p1=>lehmann_clone_Q(Q2)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				call lehmann_transform_Q(cluster,p1,miu,kx,ky)
				call lehmann_green_function(w+Xi*eps,orbit,p,G)
				glx=glx+(G(1,4)+G(4,1)*exp(-Xi*2*kx))*(G(4,1)+G(1,4)*exp(Xi*2*kx))
				deallocate(G)
				call lehmann_green_function(w+Xi*eps,orbit,p1,G)
				glx=glx+(G(1,4)+G(4,1)*exp(-Xi*2*kx))*(G(4,1)+G(1,4)*exp(Xi*2*kx))
				deallocate(G)
		call lehmann_Q_clean(p)
		call lehmann_Q_clean(p1)
			end do
		end do
		gly=Zero
		do ii=0,nn
			do jj=0,nn
				kx=ii*Pi/nn
				ky=jj*Pi/nn
				p=>lehmann_clone_Q(Q1)
				p1=>lehmann_clone_Q(Q2)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				call lehmann_transform_Q(cluster,p1,miu,kx,ky)
				call lehmann_green_function(w+Xi*eps,orbit,p,G)
				gly=gly+(G(1,2)+G(2,1)*exp(-Xi*2*kx))*(G(2,1)+G(1,2)*exp(Xi*2*kx))
				deallocate(G)
				call lehmann_green_function(w+Xi*eps,orbit,p1,G)
				gly=gly+(G(1,2)+G(2,1)*exp(-Xi*2*kx))*(G(2,1)+G(1,2)*exp(Xi*2*kx))
				deallocate(G)
			end do
		end do
		call lehmann_Q_clean(p)
		call lehmann_Q_clean(p1)
		x=glx-gly
		print *,glx,gly
	end function


	subroutine VCA_spin_pair(cluster,miu,v)
		implicit none
		real(kind=8),intent(in)::miu,v
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q1,Q2,p,p1,p2,a,a1,a2
		real(kind=8)::omega,x,kx,ky,qx,qy,coor(2,4),Pi,broad,w,wi,wo,ws
		real(kind=8),allocatable::wn(:),wn1(:),wn2(:)
		complex(kind=8)::glx,gly,qq,tt,gg
		integer::orbit,i,j,k,ii,jj,s,r,l,nw,nw1,nw2,alpha,beta,nn
		logical::OK
		complex(kind=8),allocatable,dimension(:,:)::G,Gq,qn,qn1,qn2
		if(.not.associated(cluster)) return
		Pi=asin(1.d0)*2
		broad=1.d-4
		orbit=cluster_get_n_orbit(cluster)
		allocate(G(orbit,orbit))
		allocate(Gq(orbit,orbit))
		Q1=>lehmann_init(cluster,omega,1)
		Q2=>lehmann_init(cluster,omega,2)
		if(.not.associated(Q1)) return
		if(.not.associated(Q2)) return
		coor=cluster%hop%coordinate
		wi=Pi*0.0001
		wo=5
		ws=2*Pi*0.0001
		nn=40
		do ii=0,nn
			do jj=0,nn
				qx=ii*0.4*Pi/nn !+0.8*Pi
				qy=jj*0.4*Pi/nn+0.8*Pi
				glx=Zero
				do i=0,nk-1
					do j=0,nk-1
						kx=i*Pi/nk
						ky=j*Pi/nk
						!spin up
						p=>lehmann_clone_Q(Q1)
						call lehmann_transform_Q(cluster,p,miu,kx,ky)
						p1=>lehmann_clone_Q(Q1)
						call lehmann_transform_Q(cluster,p1,miu,kx+qx,ky+qy)
!                        a=>p
!                        glx=0.0
!                        do
!                            if(.not.associated(a)) exit
!                            a1=>p1
!                            do
!                                if(.not.associated(a1)) exit
!                                if(abs(a1%w+v-a%w)<broad) then
!                                    do alpha=1,orbit
!                                        do beta=1,orbit
!                                            qq=-xi*((coor(1,beta)-coor(1,alpha))*qx+(coor(2,beta)-coor(2,alpha))*qy)
!                                            glx=glx+exp(qq)*a%q(alpha)*conjg(a%q(beta))*a1%q(beta)*conjg(a1%q(alpha))
!                                        end do
!                                    end do
!                                end if
!                                a1=>a1%next
!                            end do
!                            a=>a%next
!                        end do
						w=wi
						do
							if(w>wo) exit
							call lehmann_green_function(Xi*w,orbit,p,G)
							call lehmann_green_function(Xi*w+v,orbit,p1,Gq)
							do alpha=1,orbit
								do beta=1,orbit
									qq=-xi*((coor(1,beta)-coor(1,alpha))*qx+(coor(2,beta)-coor(2,alpha))*qy)
									glx=glx+exp(qq)*G(beta,alpha)*Gq(alpha,beta)*exp(Xi*w*0.001)
								end do
							end do
							w=w+ws
						end do
						call lehmann_Q_clean(p)
						call lehmann_Q_clean(p1)
						!spin down
						p=>lehmann_clone_Q(Q2)
						call lehmann_transform_Q(cluster,p,miu,kx,ky)
						p1=>lehmann_clone_Q(Q2)
						call lehmann_transform_Q(cluster,p1,miu,kx+qx,ky+qy)
!                        a=>p
!                        glx=0.0
!                        do
!                            if(.not.associated(a)) exit
!                            a1=>p1
!                            do
!                                if(.not.associated(a1)) exit
!                                if(abs(a1%w+v-a%w)<broad) then
!                                    do alpha=1,orbit
!                                        do beta=1,orbit
!                                            qq=-xi*((coor(1,beta)-coor(1,alpha))*qx+(coor(2,beta)-coor(2,alpha))*qy)
!                                            glx=glx+exp(qq)*a%q(alpha)*conjg(a%q(beta))*a1%q(beta)*conjg(a1%q(alpha))
!                                        end do
!                                    end do
!                                end if
!                                a1=>a1%next
!                            end do
!                            a=>a%next
!                        end do
						w=wi
						do
							if(w>wo) exit
							call lehmann_green_function(Xi*w,orbit,p,G)
							call lehmann_green_function(Xi*w+v,orbit,p1,Gq)
							do alpha=1,orbit
								do beta=1,orbit
									qq=-Xi*((coor(1,beta)-coor(1,alpha))*qx+(coor(2,beta)-coor(2,alpha))*qy)
									glx=glx+exp(qq)*G(beta,alpha)*Gq(alpha,beta)*exp(Xi*w*0.001)
								end do
							end do
							w=w+ws
						end do
						call lehmann_Q_clean(p)
						call lehmann_Q_clean(p1)
					end do
				end do
				write(*,*) qy/2/Pi,real(-2*glx/nk/nk/orbit/2/Pi)
			end do
                    exit
		end do
		call lehmann_Q_clean(Q1)
		call lehmann_Q_clean(Q2)
	end subroutine

	subroutine VCA_fermi_surface(cluster,miu)
		implicit none
		real(kind=8),intent(in)::miu
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q,p,a
		real(kind=8)::omega,x,kx,ky,Pi
		complex(kind=8)::gl
		integer::orbit,i,j,n
		logical::OK
		complex(kind=8),pointer,dimension(:,:)::G
		if(.not.associated(cluster)) return
		orbit=cluster_get_n_orbit(cluster)
		allocate(G(orbit,orbit))
		Pi=asin(1.d0)*2
		Q=>lehmann_init(cluster,omega,1)
		if(.not.associated(Q)) return
		open(unit=8,file='fermi.data',status='replace')
		do i=0,nk-1
			do j=0,nk-1
				kx=i*2*Pi/nk-Pi
				ky=j*2*Pi/nk-Pi
				p=>lehmann_clone_Q(Q)
				call lehmann_transform_Q(cluster,p,miu,kx,ky)
				call lehmann_green_function(eps*Xi,orbit,p,G)
				!call matrix_inverse(orbit,G,OK)
				gl=lehmann_restore(cluster,G,kx,ky)
				x=-imag(gl)/Pi
				write(8,*) kx,ky,x
				call lehmann_Q_clean(p)
			end do
		end do
		deallocate(G)
		call lehmann_Q_clean(Q)
		close(8)
	end subroutine

	subroutine VCA_DOS(cluster,miu)
		implicit none
		real(kind=8),intent(in)::miu
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		type(Q_type),pointer::Q,p,a
		real(kind=8)::omega,kx,ky,Pi,w
		integer,allocatable::S(:)
		complex(kind=8)::gl
		integer::orbit,i,j,n
		logical::OK
		if(.not.associated(cluster)) return
		orbit=cluster_get_n_orbit(cluster)
		Q=>lehmann_init(cluster,omega,1)
		if(.not.associated(Q)) return
		n=100
		allocate(S(n))
		call Efunc_init(n,-4.0d0,4.0d0,Q,cluster,miu)
		call calDOS(n,2,Ekxy,S)
		open(unit=10,file='DOS.data',status='replace')
		do i=1,n
			write(10,*) Ek(i),S(i)
		end do
		close(10)
		call lehmann_Q_clean(Q)
		call Efunc_clean()
	end subroutine

	function Ekxy(d,x,hit) result(k)
		integer::d,hit(nE)
		real(8)::x(d)
		k=Efunc_position(d,x,hit)
	end function

	subroutine update(cluster,v)
		implicit none
		type(hubbard_cluster_type),pointer,intent(inout)::cluster
		real(kind=8),pointer,dimension(:),intent(in)::v
		cluster%miu=v(2)
		cluster%tx=v(1)
		cluster%ty=v(1)
	end subroutine

	subroutine varinit(cluster,n,v)
		implicit none
		integer,intent(out)::n
		type(hubbard_cluster_type),pointer,intent(in)::cluster
		real(kind=8),pointer,dimension(:),intent(out)::v
		n=2
		allocate(v(n))
		v(2)=cluster%miu
		v(1)=cluster%tx
		!v(2)=cluster%ty
		!write(7,*) "n=",n
end subroutine

end module VCA
