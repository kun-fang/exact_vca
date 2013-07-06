module timer
	implicit none
	type timer_type
		real,private::init
		real,private::now
	end type timer_type
	
	public::timer_init,timer_runtime
	
	private::time_now
	
	contains

!----------------------public-----------------------
	function timer_init() result(p)
		implicit none
		type(timer_type),pointer::p
		real,dimension(2)::t
		real::a
		allocate(p)
		call etime(t,a)
		p%init=t(1)
		!print *,p%init
		!print *,t(1)
		!print *,t(2)
	end function
	
	subroutine timer_clean(p)
		implicit none
		type(timer_type),pointer::p
		deallocate(p)
		nullify(p)
	end subroutine
	
	function timer_runtime(p) result(x)
		implicit none
		type(timer_type),pointer::p
		real::x
		call time_now(p)
		x=p%now-p%init
	end function
	
!----------------------private----------------------
	subroutine time_now(p)
		implicit none
		real,dimension(2)::t
		type(timer_type),pointer::p
		real::a
		call etime(t,a)
		p%now=t(1)
		!print *,"------------------------------------------"
		!print *,p%now
		!print *,t(1)
		!print *,t(2)
	end subroutine

end module timer
