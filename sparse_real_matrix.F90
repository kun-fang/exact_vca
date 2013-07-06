module sparse_real
	implicit none
		
	complex(kind=8),private,parameter::Zero=(0.0,0.0),One=(1.0,0.0),Xi=(0.0,1.0)
	
	type real_element
		integer::i
		real(kind=8),private::x
		type(real_element),pointer,private::next
	end type

	type sparse_real_type
		integer,private::dim
		type(real_element),pointer,dimension(:),private::matrix
		type(real_element),pointer,private::cursor
		integer,private::cur_row,cur_col
	end type

	public::sparse_real_init,sparse_real_clean,sparse_real_show_dim,sparse_real_show_element
	public::sparse_real_cursor_reset,sparse_real_cursor_next,sparse_real_cursor_read,sparse_real_cursor_assign
	public::sparse_real_add_element,sparse_real_multi_element,sparse_real_multi_vector,sparse_real_connect
	public::sparse_real_multi_matrix,sparse_real_transpose,sparse_real_clone
	
	private::newnode,insertnode,delnode

	contains

!-----------------------public----------------------

	function sparse_real_init(dim) result(p)
		implicit none
		integer,intent(in)::dim
		type(sparse_real_type),pointer::p
		integer::i
		if(dim<=0) then
			nullify(p)
			return
		end if
		allocate(p)
		p%dim=dim
		allocate(p%matrix(dim))
		do i=1,dim
			p%matrix(i)%i=i
			p%matrix(i)%x=0.0
			nullify(p%matrix(i)%next)
		end do
		p%cursor=>p%matrix(1)
		p%cur_row=1
		p%cur_col=1
	end function

	subroutine sparse_real_clean(p)
		implicit none
		type(sparse_real_type),pointer,intent(inout)::p
		type(real_element),pointer::a,d
		integer::i
		if(.not.associated(p)) return
		do i=1,p%dim
			d=>p%matrix(i)
			d=>d%next
			do
				if(.not.associated(d)) exit
				a=>d
				d=>d%next
				deallocate(a)
			end do
		end do
		deallocate(p%matrix)
		deallocate(p)
		nullify(p)
	end subroutine
	
	function sparse_real_show_dim(p) result(x)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer::x
		x=0
		if(.not.associated(p)) return
		x=p%dim
	end function
	
	function sparse_real_show_element(p,i,j) result(x)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer,intent(in)::i,j
		real(kind=8)::x
		type(real_element),pointer::d
		x=0.0
		if(.not.associated(p)) return
		d=>p%matrix(i)
		if(i==j) then
			x=d%x
			return
		end if
		d=>d%next
		do
			if(.not.associated(d)) exit
			if(d%i==j) then
				x=d%x
				return
			end if
			d=>d%next
		end do
	end function
	
	subroutine sparse_real_cursor_reset(p)
		implicit none
		type(sparse_real_type),pointer,intent(inout)::p
		if(.not.associated(p)) return
		p%cursor=>p%matrix(1)
		p%cur_row=1
		p%cur_col=1
	end subroutine
	
	subroutine sparse_real_cursor_next(p,OK)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		logical,intent(inout)::OK
		if(.not.associated(p)) then
			OK=.false.
			return
		end if
		OK=.true.
		if(associated(p%cursor%next)) then
			p%cursor=>p%cursor%next
			p%cur_col=p%cursor%i
		else
			p%cur_row=p%cur_row+1
			if(p%cur_row>p%dim) then
				OK=.false.
				return
			end if
			p%cur_col=p%cur_row
			p%cursor=>p%matrix(p%cur_col)
		end if
	end subroutine
	
	subroutine sparse_real_cursor_read(p,i,j,x)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer,intent(out)::i,j
		real(kind=8),intent(out)::x
		i=0
		j=0
		x=0.0
		if(.not.associated(p)) return
		i=p%cur_row
		j=p%cur_col
		x=p%cursor%x
	end subroutine
	
	subroutine sparse_real_cursor_assign(p,x)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		real(kind=8),intent(in)::x
		type(real_element),pointer::a,b
		if(.not.associated(p)) return
		p%cursor%x=x
		if(abs(x).lt.1.d-10.and.p%cur_col/=p%cur_row) then
			b=>p%matrix(p%cur_row)
			a=>b
			b=>b%next
			do
				if(b%i==p%cur_col) then
					p%cursor=>a
					p%cur_col=a%i
					call delnode(a,b)
					return
				end if
				a=>b
				b=>b%next
			end do
		end if
	end subroutine
	
	subroutine sparse_real_add_element(p,i,j,x)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer,intent(in)::i,j
		real(kind=8)::x
		type(real_element),pointer::a,d
		if(.not.associated(p)) return
		if(abs(x)<1.d-10) return
		d=>p%matrix(i)
		if(i==j) then
			d%x=d%x+x
			return
		end if
		a=>d
		d=>d%next
		do
			if(.not.associated(d)) exit
			if(d%i==j) then
				d%x=d%x+x
				if(abs(d%x)<1.d-10) call delnode(a,d)
				return
			end if
			if(d%i>j) exit
			a=>d
			d=>d%next
		end do
		d=>newnode()
		d%i=j
		d%x=x
		call insertnode(a,d)
	end subroutine
	
	subroutine sparse_real_multi_element(p,i,j,x)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer,intent(in)::i,j
		real(kind=8)::x
		type(real_element),pointer::a,d
		if(.not.associated(p)) return
		d=>p%matrix(i)
		if(i==j) then
			d%x=d%x*x
			return
		end if
		a=>d
		d=>d%next
		do
			if(.not.associated(d)) exit
			if(d%i==j) then
				d%x=d%x*x
				if(abs(d%x)<1.d-10) call delnode(a,d)
				return
			end if
			if(d%i>j) exit
			a=>d
			d=>d%next
		end do
	end subroutine	

	function sparse_real_multi_vector(p,x) result(y)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		real(kind=8),pointer,dimension(:),intent(in)::x
		real(kind=8),pointer,dimension(:)::y
		type(real_element),pointer::d
		integer::i
		real(kind=8)::sum
		if(.not.associated(p)) return
		allocate(y(p%dim))
		do i=1,p%dim
			d=>p%matrix(i)
			sum=0.0
			do
				if(.not.associated(d)) exit
				sum=sum+d%x*x(d%i)
				d=>d%next
			end do
			y(i)=sum
		end do
	end function
	
	function sparse_real_multi_matrix(p1,p2) result(p)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p1,p2
		type(sparse_real_type),pointer::d,p
		type(real_element),pointer::d1,d2
		integer::i,j
		real(kind=8)::sum
		if(.not.associated(p1).or..not.associated(p2)) return
		if(p1%dim/=p2%dim) return
		p=>sparse_real_init(p1%dim)
		d=>sparse_real_transpose(p2)
		do i=1,p1%dim
			do j=1,p2%dim
				sum=0.0
				d1=>p1%matrix(i)
				d2=>d%matrix(j)%next
	c1:		do
					if(.not.associated(d1)) exit c1
					if(d1%i==j) then
						sum=sum+d1%x*d%matrix(j)%x
					else
						do
							if(.not.associated(d2)) exit c1
							if(d1%i==d2%i) then
								sum=sum+d1%x*d2%x
								exit
							end if
							if(d2%i>d1%i) exit
							d2=>d2%next
						end do
					end if
					d1=>d1%next
				end do c1
				call sparse_real_add_element(p,i,j,sum)
			end do
		end do
		call sparse_real_clean(d)
	end function
	
	function sparse_real_connect(a,b) result(p)
		implicit none
		type(sparse_real_type),pointer,intent(in)::a,b
		type(sparse_real_type),pointer::p
		type(real_element),pointer::d,q
		integer::i,j,k,n
		if(.not.associated(a).or..not.associated(b)) return
		n=a%dim+b%dim
		p=>sparse_real_init(n)
		do i=1,a%dim
			p%matrix(i)%x=q%x
			do
				if(.not.associated(q%next)) exit
				q=>q%next
				j=q%i
				call sparse_real_add_element(p,i,j,q%x)
			end do
		end do
		do i=1,b%dim
			p%matrix(i)%x=q%x
			do
				if(.not.associated(q%next)) exit
				q=>q%next
				j=q%i
				call sparse_real_add_element(p,i+a%dim,j+a%dim,q%x)
			end do
		end do
	end function

	function sparse_real_transpose(p) result(pt)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		type(sparse_real_type),pointer::pt
		type(real_element),pointer::a,d,node
		integer::i,j,k
		if(.not.associated(p)) return
		pt=>sparse_real_init(p%dim)
		do i=1,p%dim
			a=>p%matrix(i)
			d=>pt%matrix(i)
			d%i=a%i
			d%x=a%x
		end do
		do i=1,p%dim
			a=>p%matrix(i)
			a=>a%next
			do
				if(.not.associated(a)) exit
				j=a%i
				call sparse_real_add_element(pt,j,i,a%x)
				a=>a%next
			end do
		end do
	end function

	function sparse_real_clone_matrix(p,ai,af) result(pt)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer,intent(in),optional::ai,af
		real(kind=8),pointer,dimension(:,:)::pt
		type(real_element),pointer::a,d,node
		integer::i,j,k,l,n,xi,xj
		if(.not.associated(p)) return
		if(present(ai).and.present(af)) then
			j=ai
			k=af
		else
			j=1
			k=p%dim
		end if
		n=k-j+1
		allocate(pt(n,n))
		pt(1:n,1:n)=0.0
		l=0
		do i=1,p%dim
			if(i<j) cycle
			if(i>k) exit
			l=l+1
			a=>p%matrix(i)
			xi=l
			pt(xi,xi)=a%x
			a=>a%next
			do
				if(.not.associated(a)) exit
				if(a%i>k) exit
				if(a%i>=j) then
					xj=a%i-j+1
					pt(xi,xj)=a%x
				end if
				a=>a%next
			end do
		end do
	end function


	function sparse_real_clone(p,ai,af) result(pt)
		implicit none
		type(sparse_real_type),pointer,intent(in)::p
		integer,intent(in),optional::ai,af
		type(sparse_real_type),pointer::pt
		type(real_element),pointer::a,d,node
		integer::i,j,k,l
		if(.not.associated(p)) return
		if(present(ai).and.present(af)) then
			j=ai
			k=af
		else
			j=1
			k=p%dim
		end if
		pt=>sparse_real_init(k-j+1)
		l=0
		do i=1,p%dim
			if(i<j) cycle
			if(i>k) exit
			l=l+1
			a=>p%matrix(i)
			d=>pt%matrix(l)
			d%i=l
			d%x=a%x
			a=>a%next
			do
				if(.not.associated(a)) exit
				if(a%i>k) exit
				if(a%i>=j) then
					node=>newnode()
					node%i=a%i-j+1
					node%x=a%x
					call insertnode(d,node)
					d=>d%next
				end if
				a=>a%next
			end do
		end do
	end function



!---------------------private-----------------------

	function newnode() result(node)
		implicit none
		type(real_element),pointer::node
		allocate(node)
		nullify(node%next)
		node%x=0.0
		node%i=0
	end function

	subroutine insertnode(prev,node)
		implicit none
		type(real_element),pointer,intent(inout)::prev,node
		node%next=>prev%next
		prev%next=>node
	end subroutine
	
	subroutine delnode(prev,node)
		implicit none
		type(real_element),pointer,intent(inout)::prev,node
		prev%next=>node%next
		deallocate(node)
	end subroutine

end module sparse_real
