module random
	implicit none
	
	type random_type
		integer,private::seed
		real(kind=8),private::rand
	end type
	
	public::random_init,random_new,random_vector
	
	private::init_seed,new
	
	contains
	
!-----------------------------public-----------------------------

	function random_init() result(x)
		implicit none
		type(random_type),pointer::x
		allocate(x)
		call init_seed(x)
	end function

	subroutine random_clean(p)
		type(random_type),pointer,intent(inout)::p
		if(.not.associated(p)) return
		deallocate(p)
		nullify(p)
	end subroutine

	function random_new(p,a,b) result(x)
		type(random_type),pointer,intent(inout)::p
		real(kind=8),intent(in),optional::a,b
		real(kind=8)::x
		if(.not.associated(p)) return
		call new(p)
		if(present(a).and.present(b)) then
			x=a+p%rand*(b-a)
		else
			x=p%rand
		end if
	end function

	function random_vector(p,n,a) result(x)
		type(random_type),pointer,intent(inout)::p
		real(kind=8),intent(in)::a
		integer,intent(in)::n
		real(kind=8),pointer,dimension(:)::x
		integer::i
		real(kind=8)::sum
		if(.not.associated(p)) return
		if(a<0.or.n<=0) return
		allocate(x(n))
		sum=0.0
		do i=1,n
			call new(p)
			x(i)=p%rand
			sum=sum+x(i)*x(i)
		end do
		sum=sqrt(sum)
		do i=1,n
			x(i)=x(i)*a/sum
		end do
	end function

!----------------------------private-----------------------------

	!Get initial seed for random routine
	subroutine init_seed(p)
		type(random_type),pointer,intent(inout)::p
		integer::seed
		integer*4,dimension(3)::today,now
		if(.not.associated(p)) return
		call itime(now)
		call idate(today)
		p%seed=today(3)+70*(today(2)+12*(today(1)+31*(now(1)+23*(now(2)+59*now(3)))))
	end subroutine
	
	!Get a new seed
	subroutine new(p)
		type(random_type),pointer,intent(inout)::p
		integer::a,b,m,q,r,Ixx
		if(.not.associated(p)) return
		a = 16807
		b = 0
		m = 2147483647
		q=int(m/a)
		r=mod(m,a)
		Ixx = a*(mod(p%seed,q))-r*int(p%seed/q)
		if (Ixx.lt.0) Ixx = Ixx + m
		p%seed=Ixx
		p%rand=1.0*p%seed/m
	end subroutine

end module random