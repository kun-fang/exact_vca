module rmatrix
	implicit none
	
	private::LU,determinant,inverse
	
	public::rmatrix_determinant,rmatrix_inverse
	
	contains
	
!-----------------public------------------

	function rmatrix_determinant(n,a) result(x)
		implicit none
		integer,intent(in)::n
		real(kind=8),pointer,dimension(:,:),intent(in)::a
		real(kind=8)::x
		real(kind=8),pointer,dimension(:,:)::b,c
		real(kind=8),pointer,dimension(:)::y,z
		integer,pointer,dimension(:)::pivot
		integer::i,j,sgn
		logical::OK
		allocate(b(n,n))
		allocate(pivot(n))
		sgn=1
		do i=1,n
			pivot(i)=i
			do j=1,n
				b(j,i)=a(j,i)
			end do
		end do
		call LU(n,b,pivot,OK,sgn)
		
		if(OK) then
			x=determinant(n,b,sgn)
		else
			x=0.0
		end if
		allocate(c(n,n))
		do i=1,n
			do j=1,n
				c(i,j)=0.0
				if(j>i) then
					c(i,j)=b(i,j)
					b(i,j)=0.0
				end if
			end do
			c(i,i)=1.0
		end do
		deallocate(pivot)
		deallocate(b)
	end function

	subroutine rmatrix_inverse(n,a,OK)
		implicit none
		integer,intent(in)::n
		real(kind=8),pointer,dimension(:,:),intent(inout)::a
		integer,pointer,dimension(:)::pivot
		real(kind=8),pointer,dimension(:,:)::b
		integer::i,j,sgn
		logical::OK
		allocate(pivot(n))
		sgn=1
		do i=1,n
			pivot(i)=i
		end do
		call LU(n,a,pivot,OK,sgn)
		if(OK) then
			call Inverse(n,a,pivot)
		else
			print *,"fatal error! Determinent of the matrix is zero."
			return
		end if
		deallocate(pivot)
	end subroutine
	
	
!----------------private------------------

	subroutine LU(n,a,pivot,OK,sgn)
		implicit none
		integer,intent(in)::n
		integer,intent(out)::sgn
		real(kind=8),pointer,dimension(:,:),intent(inout)::a
		integer,pointer,dimension(:),intent(inout)::pivot
		real(kind=8),pointer,dimension(:,:)::dcmp
		logical::OK
		integer::i,j,k,pi,pj,m
		real*8::max_diag
		OK=.true.
		allocate(dcmp(n,n))
		do i=1,n
			pivot(i)=i
		end do
		m=0
		do k=1,n
			pi=k
			pj=k
			max_diag=0.0
			do i=k,n
				dcmp(i,k)=a(i,k)-dot_product(dcmp(i,1:k-1),dcmp(1:k-1,k))
				if(abs(dcmp(i,k))>abs(max_diag)) then
					max_diag=dcmp(i,k)
					pi=i
					pj=k
				end if
			end do
			if(pi/=pj) then
				call pivot_row(n,a,dcmp,pivot,pi,pj)
				m=m+1
			end if
			if(abs(dcmp(k,k))<1.d-10) then
				OK=.false.
				return
			end if
			do j=k+1,n
				dcmp(k,j)=(a(k,j)-dot_product(dcmp(k,1:k-1),dcmp(1:k-1,j)))/dcmp(k,k)
			end do
		end do
		deallocate(a)
		a=>dcmp
		sgn=sgn*((-1)**m)
	end subroutine
	
	subroutine pivot_row(n,a,b,pivot,i,j)
		implicit none
		integer,intent(in)::i,j,n
		real(kind=8),pointer,dimension(:,:),intent(inout)::a,b
		integer,pointer,dimension(:),intent(inout)::pivot
		real(kind=8),pointer,dimension(:)::swap
		integer::k
		allocate(swap(n))
		swap=a(i,1:n)
		a(i,1:n)=a(j,1:n)
		a(j,1:n)=swap
		swap=b(i,1:n)
		b(i,1:n)=b(j,1:n)
		b(j,1:n)=swap
		k=pivot(i)
		pivot(i)=pivot(j)
		pivot(j)=k
		deallocate(swap)
	end subroutine

	function determinant(n,a,sgn) result(det)
		implicit none
		integer,intent(in)::n,sgn
		real(kind=8),pointer,dimension(:,:),intent(in)::a
		real(kind=8)::det
		integer::i
		det=1.0*sgn
		do i=1,n
			det=det*a(i,i)
		end do
	end function

	subroutine inverse(n,a,pivot)
		implicit none
		integer,intent(in)::n
		real(kind=8),pointer,dimension(:,:),intent(inout)::a
		integer,pointer,dimension(:),intent(in)::pivot
		real(kind=8),pointer,dimension(:,:)::b
		real(kind=8),pointer,dimension(:)::y
		integer::i,j,k,l
		allocate(y(n))
		allocate(b(n,n))
		do k=1,n
			do i=1,n
				if(pivot(i)==k) then
					l=i
					exit
				end if
			end do
			do i=1,n
				if(i==l) then
					y(i)=(1.0-dot_product(a(i,1:i-1),y(1:i-1)))/a(i,i)
				else
					y(i)=-dot_product(a(i,1:i-1),y(1:i-1))/a(i,i)
				end if
			end do
			do i=n,1,-1
				b(i,k)=y(i)-dot_product(a(i,i+1:n),b(i+1:n,k))
			end do
		end do
		deallocate(y)
		deallocate(a)
		a=>b
	end subroutine
	

end module rmatrix

module cmatrix
	implicit none
	
	complex(kind=8),parameter,private::Zero=(0.0,0.0),Xi=(0.0,1.0),One=(1.0,0.0)

	private::LU,determinant,inverse,vector_product
	
	public::matrix_determinant,matrix_inverse
	
	contains
	
!-----------------public------------------

	function matrix_determinant(n,a) result(x)
		implicit none
		integer,intent(in)::n
		complex(kind=8),pointer,dimension(:,:),intent(in)::a
		complex(kind=8)::x
		complex(kind=8),pointer,dimension(:,:)::b,c
		complex(kind=8),pointer,dimension(:)::y,z
		integer,pointer,dimension(:)::pivot
		integer::i,j,sgn
		logical::OK
		allocate(b(n,n))
		allocate(pivot(n))
		sgn=1
		do i=1,n
			pivot(i)=i
			do j=1,n
				b(j,i)=a(j,i)
			end do
		end do
		call LU(n,b,pivot,OK,sgn)
		
		if(OK) then
			x=determinant(n,b,sgn)
		else
			x=0.0
		end if
		allocate(c(n,n))
		do i=1,n
			do j=1,n
				c(i,j)=Zero
				if(j>i) then
					c(i,j)=b(i,j)
					b(i,j)=Zero
				end if
			end do
			c(i,i)=One
		end do
		deallocate(pivot)
		deallocate(b)
	end function

	subroutine matrix_inverse(n,a,OK)
		implicit none
		integer,intent(in)::n
		complex(kind=8),pointer,dimension(:,:),intent(inout)::a
		integer,pointer,dimension(:)::pivot
		integer::i,j,sgn
		logical::OK
		allocate(pivot(n))
		sgn=1
		do i=1,n
			pivot(i)=i
		end do
		call LU(n,a,pivot,OK,sgn)
		if(OK) then
			call Inverse(n,a,pivot)
		else
			print *,"fatal error! Determinent of the matrix is zero."
			return
		end if
		deallocate(pivot)
	end subroutine
	
	
!----------------private------------------

	subroutine LU(n,a,pivot,OK,sgn)
		implicit none
		integer,intent(in)::n
		integer,intent(out)::sgn
		complex(kind=8),pointer,dimension(:,:),intent(inout)::a
		integer,pointer,dimension(:),intent(inout)::pivot
		complex(kind=8),pointer,dimension(:,:)::dcmp
		logical::OK
		integer::i,j,k,pi,pj,m
		complex*8::max_diag
		OK=.true.
		allocate(dcmp(n,n))
		do i=1,n
			pivot(i)=i
		end do
		m=0
		do k=1,n
			pi=k
			pj=k
			max_diag=Zero
			do i=k,n
				dcmp(i,k)=a(i,k)-vector_product(k-1,dcmp(i,1:k-1),dcmp(1:k-1,k))
				if(abs(dcmp(i,k))>abs(max_diag)) then
					max_diag=dcmp(i,k)
					pi=i
					pj=k
				end if
			end do
			if(pi/=pj) then
				call pivot_row(n,a,dcmp,pivot,pi,pj)
				m=m+1
			end if
			if(abs(dcmp(k,k))<1.d-10) then
				OK=.false.
				return
			end if
			do j=k+1,n
				dcmp(k,j)=(a(k,j)-vector_product(k-1,dcmp(k,1:k-1),dcmp(1:k-1,j)))/dcmp(k,k)
			end do
		end do
		deallocate(a)
		a=>dcmp
		sgn=sgn*((-1)**m)
	end subroutine
	
	subroutine pivot_row(n,a,b,pivot,i,j)
		implicit none
		integer,intent(in)::i,j,n
		complex(kind=8),pointer,dimension(:,:),intent(inout)::a,b
		integer,pointer,dimension(:),intent(inout)::pivot
		complex(kind=8),pointer,dimension(:)::swap
		integer::k
		allocate(swap(n))
		swap=a(i,1:n)
		a(i,1:n)=a(j,1:n)
		a(j,1:n)=swap
		swap=b(i,1:n)
		b(i,1:n)=b(j,1:n)
		b(j,1:n)=swap
		k=pivot(i)
		pivot(i)=pivot(j)
		pivot(j)=k
		deallocate(swap)
	end subroutine

	function determinant(n,a,sgn) result(det)
		implicit none
		integer,intent(in)::n,sgn
		complex(kind=8),pointer,dimension(:,:),intent(in)::a
		complex(kind=8)::det
		integer::i
		det=One*sgn
		do i=1,n
			det=det*a(i,i)
		end do
	end function

	subroutine inverse(n,a,pivot)
		implicit none
		integer,intent(in)::n
		complex(kind=8),pointer,dimension(:,:),intent(inout)::a
		integer,pointer,dimension(:),intent(in)::pivot
		complex(kind=8),pointer,dimension(:,:)::b
		complex(kind=8),pointer,dimension(:)::y
		integer::i,j,k,l
		allocate(y(n))
		allocate(b(n,n))
		do k=1,n
			do i=1,n
				if(pivot(i)==k) then
					l=i
					exit
				end if
			end do
			do i=1,n
				if(i==l) then
					y(i)=(One-vector_product(i-1,a(i,1:i-1),y(1:i-1)))/a(i,i)
				else
					y(i)=-vector_product(i-1,a(i,1:i-1),y(1:i-1))/a(i,i)
				end if
			end do
			do i=n,1,-1
				b(i,k)=y(i)-vector_product(n-i,a(i,i+1:n),b(i+1:n,k))
			end do
		end do
		deallocate(y)
		deallocate(a)
		a=>b
	end subroutine
	
	function vector_product(n,a,b) result(x)
		implicit none
		integer,intent(in)::n
		complex(kind=8),dimension(n),intent(in)::a,b
		complex(kind=8)::x
		integer::i
		x=Zero
		do i=1,n
			x=x+a(i)*b(i)
		end do
	end function



end module cmatrix
