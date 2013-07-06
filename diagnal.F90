module diagnal

	private
	
	public::diagnal_exact_real,diagnal_exact_complex,diagnal_tridiag_real
	
	private::diagm,rebakc,reducc,tql2,trbakc,tredc,ranff
	
	contains
!--------------------------public-------------------------------

	subroutine diagnal_exact_real(n,H,e,X)
		implicit none
		integer,intent(in)::n
		real(kind=8),dimension(n,n),intent(inout)::H
		real(kind=8),dimension(n,n),intent(inout),optional::X
		real(kind=8),dimension(n),intent(inout)::e
		real(kind=8),dimension(n,n)::A,Zi,Zr
		integer::i,j,k
		do i=1,n
			A(i,i)=H(i,i)
			do j=i+1,n
				A(j,i)=real(H(j,i))
				A(i,j)=0.0
			end do
		end do
		call diagm(n,n,A,A,e,Zr,Zi)
		if(present(X)) then
			do i=1,n
				do j=1,n
					X(j,i)=Zr(j,i)
				end do
			end do
		end if
	end subroutine
	
	
	subroutine diagnal_exact_complex(n,H,e,X)
		implicit none
		integer,intent(in)::n
		complex(kind=8),dimension(n,n),intent(inout)::H
		complex(kind=8),dimension(n,n),intent(inout),optional::X
		real(kind=8),dimension(n),intent(inout)::e
		real(kind=8),dimension(n,n)::A,Zi,Zr
		integer::i,j,k
		do i=1,n
			A(i,i)=real(H(i,i))
			do j=i+1,n
				A(j,i)=real(H(j,i))
				A(i,j)=imag(H(j,i))
			end do
		end do
		call diagm(n,n,A,A,e,Zr,Zi)
		if(present(X)) then
			do i=1,n
				do j=1,n
					X(j,i)=Zr(j,i)+(0.0,1.0)*Zi(j,i)
				end do
			end do
		end if
	end subroutine
	
	subroutine diagnal_tridiag_real(n,d,e,Z)
		implicit none
		integer,intent(in)::n
		real(kind=8),dimension(n),intent(inout)::d,e
		real(kind=8),dimension(n,n),intent(inout),optional::Z
		real(kind=8),pointer,dimension(:,:)::X
		integer::ierr
		if(present(Z)) then
			call tql2(n,n,d,e,Z,ierr)
		else
			allocate(X(n,n))
			call tql2(n,n,d,e,X,ierr)
			deallocate(X)
		end if
		if(ierr/=0) stop 'tql2'
	end subroutine
	
	
!--------------------------private------------------------------
	
      subroutine diagm(nm,n,h,s,e,zr,zi)
!**********************************************************************
!
!     solves the complex hermitian eigenvalue problem (h-e*s)z=0
!     for matrices h and s of order n.  only the lower triangle
!     of the hermitian matrices need be supplied.  if a(i,j) is a
!     real matrix and b is a hermitian matrix, then b is stored
!     in a in the following way (i.ge.j):
!       a(i,j) = real ( b(i,j) )
!       a(j,i) = imag ( b(i,j) )
!     since the diagonal elements of b are real, there is no need
!     to store the zero imaginary parts.
!
!                  m. weinert   july 1983
!
!     nm         first dimension of arrays
!     n          order of matrices
!     h          hamiltonian matrix (overwritten on output)
!     s          overlap matrix (overwritten on output)
!     e          eigenvalues
!     zr,zi      eigenvectors (real, imaginary parts)
!     e1,e2,tau  work arrays
!     time       time required for diagonalization
!
!     modified for use in BEST  Fall 1993 - Spring 1994
!**********************************************************************
!
      implicit none
      integer  i,ierr,j,n,nm
      real*8   t1,t2,time

!--->    hamiltonian and overlap matrices
      real*8 s(nm,n),h(nm,n)
!--->   eigenvalues and eigenvectors
      real*8 e(n),zr(nm,n),zi(nm,n)
!--->   work arrays
      real*8 e1(n),e2(n),tau(2,n)
!
!--->    reduce the general problem to the standard problem
!      call reducc(nm,n,h,s)
!--->    reduce the standard problem to real tridiagonal form
      call tredc(nm,n,h,e,e1,e2,tau)
!--->    find eigenvalues and eigenvectors of real triadiagonal matrix
      do i=1,n
         do j=1,n
            zr(j,i)=0.0d0
         enddo
         zr(i,i)=1.0d0
      enddo
      call tql2(nm,n,e,e1,zr,ierr)
      if(ierr.ne.0) stop 'tql2'
!--->    back-transform the eigenvectors to the standard problem
      call trbakc(nm,n,h,tau,n,zr,zi)
!--->    back-transform the eigenvectors to the original problem
!      call rebakc(nm,n,n,s,zr,zi)
      return
      end subroutine
      
      
      subroutine rebakc(nm,n,m,b,zr,zi)
!**********************************************************************
!     complex version of the algol procedure rebaka, linear
!     algebra, vol. ii, 1971 by wilkinson and reinsch.
!     forms the eigenvectors of the generalized hermitian
!     eigensystem by back-transforming those of the derived
!     standard matrix determined by reducc.
!     input:
!      nm     row dimension of the 2-d arrays
!      n      order of the matrix system
!      m      number of eigenvectors to back-transform
!      b      contains the cholesky decomposition obtained
!             in reducc in compact storage mode.
!      zr,zi  contain the real and imaginary parts of the
!             eigenvectors to be back-transformed in the
!             first m columns.
!     output:
!      zr,zi  contain the back-transformed eigenvectors
!                 m. weinert   july 1983
!**********************************************************************
      implicit none
      integer  i,j,k,m,n,nm
      real*8   xi,xr
      real*8   b(nm,n),zr(nm,n),zi(nm,n)

      do j=1,m
         do i=n,1,-1
            xr=zr(i,j)
            xi=zi(i,j)
            do k=i+1,n
               xr=xr - b(k,i)*zr(k,j) - b(i,k)*zi(k,j)
               xi=xi - b(k,i)*zi(k,j) + b(i,k)*zr(k,j)
            enddo
            zr(i,j)=xr/b(i,i)
            zi(i,j)=xi/b(i,i)
         enddo
      enddo
      return
      end subroutine
      
      
      subroutine reducc(nm,n,a,b)
!**********************************************************************
!     complex version of the algol procedure reduc1, linear
!     algebra, vol. ii, wilkinson and reinsch, 1971.
!
!     reduction of the general hermitian eigenvalue problem
!     a*x=lamda*b*x to the equivalent problem p*z=lamda*z
!     using cholesky decomposition.
!
!     the procedure will fail if b, perhaps due to rounding
!     errors, is not positive definite.
!
!     input:
!      nm     row dimension of 2-d arrays a and b as declared in
!             call routine
!      n      order of the eigensystem
!      a,b    the lower triangle of the hermitian matrices stored
!             in compact mode as:  (i.ge.j)
!                a(i,j) = real ( (h(i,j) )
!                a(j,i) = imag ( (h(i,j) )
!
!     output:
!      a      contains the lower triangle of the reduced problem
!             stored in compact mode
!      b      contains the lower triangular cholesky decomposition
!             stored in compact mode
!
!                     m. weinert     july 1983
!
!     note that a is in the form required in tredc and b is in the
!     form required in rebakc.
!**********************************************************************
      implicit none

      integer  i,j,k,n,nm
      real*8   a(nm,n),b(nm,n)
      real*8   xi,xr,y,zero

      parameter(zero=0.0d0)

!--->    form l in lower triangle of b
      do j=1,n
         xr=b(j,j)
         do k=1,j-1
            xr=xr - b(j,k)*b(j,k) - b(k,j)*b(k,j)
         enddo
         if(xr.le.zero) stop 'reducc'
         y=sqrt(xr)
         b(j,j)=y
         do i=j+1,n
            xr=b(i,j)
            xi=b(j,i)
            do k=1,j-1
               xr=xr - b(i,k)*b(j,k) - b(k,i)*b(k,j)
               xi=xi - b(k,i)*b(j,k) + b(i,k)*b(k,j)
            enddo
            b(j,i)=xi/y
            b(i,j)=xr/y
         enddo
      enddo
!--->    form hermitian conjugate of inv(l)*a
      do j=1,n
         y=b(j,j)
         xr=a(j,j)
         do k=1,j-1
            xr=xr - b(j,k)*a(j,k) - b(k,j)*a(k,j)
         enddo
         a(j,j)=xr/y
         do i=j+1,n
            xr=a(i,j)
            xi=a(j,i)
            do k=1,j-1
               xr=xr - b(j,k)*a(i,k) - b(k,j)*a(k,i)
               xi=xi + b(k,j)*a(i,k) - b(j,k)*a(k,i)
            enddo
            a(j,i)=xi/y
            a(i,j)=xr/y
         enddo
      enddo
!--->    premultiply by inv(l)
      do i=1,n
         y=b(i,i)
         do j=1,i-1
            xr=a(i,j) - b(i,j)*a(j,j)
            xi=a(j,i) - b(j,i)*a(j,j)
            do k=j+1,i-1
               xr=xr - b(i,k)*a(k,j) + b(k,i)*a(j,k)
               xi=xi - b(k,i)*a(k,j) - b(i,k)*a(j,k)
            enddo
            do k=1,j-1
               xr=xr - b(i,k)*a(j,k) - b(k,i)*a(k,j)
               xi=xi - b(k,i)*a(j,k) + b(i,k)*a(k,j)
            enddo
            a(j,i)=xi/y
            a(i,j)=xr/y
         enddo
         xr=a(i,i)
         do k=1,i-1
            xr=xr - b(i,k)*a(i,k) - b(k,i)*a(k,i)
         enddo
         a(i,i)=xr/y
      enddo

      return
      end subroutine
      
      
      
      

      subroutine tql2(nm,n,d,e,z,ierr)
!*********************************************************************
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!        n is the order of the matrix.
!        d contains the diagonal elements of the input matrix.
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!        e has been destroyed.
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!     questions and comments should be directed to b. s. garbow,
!     applied mathematics division, argonne national laboratory
!     ------------------------------------------------------------------
      implicit none
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real*8 d(n),e(n),z(nm,n)
      real*8 b,c,c2,c3,dl1,el1,f,g,h,p,r,s,s2
!     .......... first executable statement  tql2 .........
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
  100 e(i-1) = e(i)
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (b .lt. h) b = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            if (b + abs(e(m)) .eq. b) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = sqrt( 1.0d0 + p**2 )
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         do 140 i = l2, n
  140    d(i) = d(i) - h
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p) .lt. abs(e(i))) go to 150
            c = e(i) / p
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * p * r
            s = c / r
            c = 1.0d0 / r
            go to 160
  150       c = p / e(i)
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * e(i) * r
            s = 1.0d0 / r
            c = c * s
  160       p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
  200    continue
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         if (b + abs(e(l)) .gt. b) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
  300 continue
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
      
      
      
      subroutine trbakc(nm,n,a,tau,m,zr,zi)
!*********************************************************************
!     complex version of the algol procedure trbak3 in
!     linear algebra, vol. ii, 1971 by wilkiinson and reinsch.
!     forms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by tredc.
!     input:
!      nm     row dimension of 2-d arrays
!      n      order of the matrix system
!      a      contains information about the unitary transformations
!             used in the reduction by tredc
!      tau    contains further information about the transformations
!      m      number of eigenvectors to be back transformed
!      zr     contains the eigenvectors to be back transformed
!             in its first m columns
!     output:
!      zr, zi contain the real and imaginary parts, respectively,
!             of the transformed eigenvectors in the first m columns
!
!     the last component of each returned vector is real and the
!     vector euclidean norms are preserved.
!                m. weinert
!*********************************************************************
      implicit none
      integer  i,j,k,m,n,nm
      real*8   a(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      real*8   h,s,si
!
      if (m .eq. 0) return
!--->   transform eigenvectors of the real symmetric tridiagonal matrix
!--->   to those of the hermitian tridiagonal matrix
      do k=1,n
         do j=1,m
            zi(k,j)= -zr(k,j)*tau(2,k)
            zr(k,j)=  zr(k,j)*tau(1,k)
         enddo
      enddo
!--->   apply the householder transformations
      do i=2,n
         h=a(i,i)
         if(h.ne.0.0d0) then
            do j=1,m
               s =0.0d0
               si=0.0d0
               do k=1,i-1
                  s =s +a(i,k)*zr(k,j)
                  si=si+a(i,k)*zi(k,j)
               enddo
               do k=1,i-1
                  s =s -a(k,i)*zi(k,j)
                  si=si+a(k,i)*zr(k,j)
               enddo
               s =(s /h)/h
               si=(si/h)/h
               do k=1,i-1
                  zr(k,j)=zr(k,j)-s *a(i,k)-si*a(k,i)
                  zi(k,j)=zi(k,j)-si*a(i,k)+s *a(k,i)
               enddo
            enddo
         endif
      enddo
      return
      end subroutine
      
      
      subroutine tredc(nm,n,a,d,e,e2,tau)
!*********************************************************************
!
!     complex version of the algol procedure tred3, linear
!     algebra, vol. ii, wilkinson and reinsch, 1971.
!     reduces a complex hermitian matrix, stored as a single
!     square array, to a real symmetric tridiagonal matrix
!     using unitary similarity transformations.
!
!     on input
!        nm must be set to the row dimension of two-dimensional
!           array parameters as declared in the calling program
!           dimension statement.
!        n  is the order of the matrix.
!        a  contains the lower triangle of the complex hermitian input
!           matrix.  the real parts of the matrix elements are stored
!           in the full lower triangle of a, and the imaginary parts
!           are stored in the transposed positions of the strict upper
!           triangle of a.  no storage is required for the zero
!           imaginary parts of the diagonal elements.
!
!     on output
!        a   contains information about the unitary transformations
!            used in the reduction.
!        d   contains the diagonal elements of the the tridiagonal
!            matrix.
!        e   contains the subdiagonal elements of the tridiagonal
!            matrix in its last n-1 positions.  e(1) is set to zero.
!        e2  contains the squares of the corresponding elements of e.
!            e2 and e must be separate locations.
!        tau contains further information about the transformations.
!
!                    m. weinert  june 1990
!*********************************************************************
      implicit none
      integer  i,j,k,l,n,ii,nm,jm1,jp1
      real*8   a(nm,n),d(n),e(n),e2(n),tau(2,n)
      real*8   f,g,h,fi,gi,hh,si,scale
!
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
!
      do 300 i=n,2,-1
!--->    use d and e2 has temporary storage
        do k=1,i-1
           d(k) =a(i,k)
           e2(k)=a(k,i)
        enddo
!--->    scale rows
        scale=0.0d0
        do k=1,i-1
           scale = scale + abs(d(k)) + abs(e2(k))
        enddo
!--->    if scale is too small to guarantee orthogonality,
!--->    transformation is skipped
        h=0.0d0
        if (scale .eq. 0.0d0) then
           tau(1,i-1)=1.0d0
           tau(2,i-1)=0.0d0
           e(i) =0.0d0
           e2(i)=0.0d0
           go to 200
        endif
        do k=1,i-1
           d(k) = d(k)/scale
           e2(k)=e2(k)/scale
        enddo
        h=0.0d0
        do k=1,i-1
           h=h + d(k)*d(k) + e2(k)*e2(k)
        enddo
        e2(i)=scale*scale*h
        g=sqrt(h)
        e(i)=scale*g
        f=sqrt( d(i-1)**2 + e2(i-1)**2 )
!--->     form next diagonal element
        if (f.eq.0.0d0) then
           tau(1,i-1) = -tau(1,i)
           si=tau(2,i)
           d(i-1)=g
           a(i,i-1)=scale*d(i-1)
        else
           tau(1,i-1)=(e2(i-1) * tau(2,i) -  d(i-1) * tau(1,i))/f
           si        =( d(i-1) * tau(2,i) + e2(i-1) * tau(1,i))/f
           h = h + f * g
           g = 1.0d0 + g / f
           d(i-1) = g *  d(i-1)
           e2(i-1)= g * e2(i-1)
           a(i,i-1) = scale* d(i-1)
           a(i-1,i) = scale*e2(i-1)
        endif
!--->     form element of a*u
        f = 0.0d0
        do j=1,i-1
           g =0.0d0
           gi=0.0d0
           do k=1,j-1
              g = g  + a(j,k) *  d(k)
              gi= gi - a(j,k) * e2(k)
           enddo
           do k=1,j-1
              g = g  + a(k,j) * e2(k)
              gi= gi + a(k,j) *  d(k)
           enddo
           g = g  + a(j,j) *  d(j)
           gi= gi - a(j,j) * e2(j)
           do k = j+1,i-1
              g = g  + a(k,j) *  d(k)
              gi= gi - a(k,j) * e2(k)
           enddo
           do k = j+1,i-1
              g = g  - a(j,k) * e2(k)
              gi= gi - a(j,k) *  d(k)
           enddo
!--->    form element of p
           e(j) = g / h
           tau(2,j) = gi / h
           f = f + e(j) * d(j) - tau(2,j) * e2(j)
        enddo

        hh = f / (h + h)
!--->    form reduced a
        do j=1,i-1
           f = d(j)
           g = e(j) - hh * f
           e(j) = g
           fi = -e2(j)
           gi = tau(2,j) - hh * fi
           tau(2,j) = -gi
           a(j,j) = a(j,j) - 2.0d0 * (f * g + fi * gi)
           do k=1,j-1
              a(j,k) = a(j,k) - f * e(k) - g * d(k) + fi * tau(2,k) + gi * e2(k)
           enddo
           do k=1,j-1
              a(k,j) = a(k,j) - f * tau(2,k) - g * e2(k) - fi * e(k) - gi * d(k)
           enddo
        enddo

        tau(2,i-1) = -si
  200   d(i) = a(i,i)
        a(i,i) = scale * sqrt(h)
  300 continue
      e(1)=0.0d0
      e2(1)=0.0d0
      d(1)=a(1,1)
      a(1,1)=0.0d0
      return
      end subroutine
      
      
      function ranff(idum)
!
!     park and miller 'minimal' random number generator with 
!     bays-durham shuffle and safeguards. returns a uniform deviate
!     in (0.0,1.0). initialize with idum a negative integer, but
!     do not change with calls. based on the routine ran1 in 
!     press, et al. (2nd edition), p. 271.
!
      implicit none
      integer  ia,im,iq,ir,ntab,ndiv,j,k,idum
      real*8   ranff,am,eps,rnmx

      parameter(ia=16807,im=2147483647,am=1.d0/im,iq=127773,ir=2836)
      parameter(ntab=32,ndiv=1+(im-1)/ntab)
      parameter(eps=2.d0**(-46),rnmx=1.d0-eps)
!
      integer iv(ntab),iy
      save iv,iy
      data iy/0/
!--->    initialize and make sure idum.ne.0
      if(idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if(idum.lt.0) idum=idum+im
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
!--->    get random number and refill shuffle table
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if(idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ranff=min(am*iy,rnmx)
      return
      end function

end module diagnal
	
	
