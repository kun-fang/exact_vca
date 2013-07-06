module MCDOS
	use random
	type(random_type),pointer,private::rp
	integer::nE
	real(8),allocatable::Ek(:)
	real(8),parameter,private::epsilon=1.d-8
	integer,parameter,private::MAXSTEP=40000

contains

	subroutine initDOS(n,Emin,Emax)
		implicit none
		integer::n,i
		real(8)::Emin,Emax,dE
		nE=n
		rp=>random_init()
		allocate(Ek(n))
		dE=(Emax-Emin)/n
		do i=1,n
			Ek(i)=Emin+dE*(i-1)
		end do
	end subroutine
	
	subroutine cleanDOS()
		implicit none
		call random_clean(rp)
		deallocate(Ek)
	end subroutine

	subroutine calDOS(n,d,calE,S)
		implicit none
		integer::n,d,hit(n),S(n),nh,i,k
		real(8)::x(d),f
		integer,external::calE
		logical::OK
		S(1:n)=0
		x(1:d)=0.0
		f=1.0
		do i=1,MAXSTEP
			call evolve(d,x)
			nh=calE(d,x,hit)
			S=S+hit*f
		end do
		k=0
		do i=1,n
			k=k+S(n)
		end do
		do i=1,n
			S(n)=S(n)*2/n
		end do
	end subroutine

	subroutine evolve(d,x)
		integer::d,i
		real(8)::x(d),pi,dx
		pi=2*asin(1.d0)
		do i=1,d
			dx=random_new(rp,0.d0,2*pi)
			x(i)=x(i)+dx
			if(x(i)>2*pi) x(i)=x(i)-2*pi
		end do
	end subroutine

end module MCDOS
