module optimal
	use rmatrix
	implicit none

	real(kind=8),parameter,private::Tol=5.d-5,Step=5.d-2
	
	type optimal_type
		real(kind=8)::dx
		integer::n
		real(kind=8),pointer,dimension(:)::v
	end type optimal_type

	public::optimal_init,optimal_clean,optimal_stationary

	contains
	
	function optimal_init(n,v,dx) result(p)
		implicit none
		integer,intent(in)::n
		real(kind=8),pointer,dimension(:),intent(in)::v
		real(kind=8),intent(in),optional::dx
		type(optimal_type),pointer::p
		allocate(p)
		p%n=n
		allocate(p%v(n))
		p%v(1:n)=v(1:n)
		p%dx=Step/10
		if(present(dx)) p%dx=dx
	end function

	subroutine optimal_clean(p)
		implicit none
		type(optimal_type),pointer::p
		if(.not.associated(p)) return
		deallocate(p%v)
		deallocate(p)
	end subroutine

	function optimal_stationary(p,f,d,q,scaling) result(conv)
		implicit none
		type(optimal_type),pointer,intent(in)::p
		real(kind=8),intent(in)::f,scaling
		real(kind=8),pointer,dimension(:),intent(in)::d
		real(kind=8),pointer,dimension(:,:),intent(in)::q
		real(kind=8),pointer,dimension(:)::df,move
		real(kind=8),pointer,dimension(:,:)::d2f
		real(kind=8)::dx,Z
		logical::OK,conv
		integer::n,i,j
		conv=.false.
		if(.not.associated(p)) return
		n=p%n
		dx=p%dx
		allocate(move(n))
		allocate(df(n))
		allocate(d2f(n,n))
		do i=1,n
			df(i)=(d(i)-f)/dx
		end do
		Z=sqrt(dot_product(df,df))
		if(Z<Tol) then
			conv=.true.
			deallocate(df)
			deallocate(d2f)
			deallocate(move)
			return
		end if
		do i=1,n
			do j=1,n
				d2f(i,j)=(q(i,j)+f-d(i)-d(j))/(dx*dx)
			end do
		end do
		call rmatrix_inverse(n,d2f,OK)
		move=matmul(d2f,df)
		Z=sqrt(dot_product(move,move))
		if(Z>Step) then
			move=scaling*Step*move
			p%dx=Step*0.618
		else
			move=scaling*move
			p%dx=Z*0.618
		end if
		p%v=p%v-move
		deallocate(df)
		deallocate(d2f)
		deallocate(move)
	end function


end module optimal
