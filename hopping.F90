module hop_mod
	implicit none

	type,public::cluster
		integer::site,dim,translation
		real,pointer,dimension(:,:)::coordinate
		real,pointer,dimension(:,:)::vector
		integer,pointer,dimension(:,:,:)::hop
		logical::IsInit=.false.
	end type cluster

	contains

	subroutine cluster_site(cl,n,a)
		class(cluster)::cl
		integer::n,dim
		real,dimension(*)::a
		dim=cl%dim
		a(1:dim)=cl%coordinate(1:dim,n)
	end subroutine

	function cluster_site_distance(dim,a,b) result(r)
		integer::dim,i
		real,dimension(dim)::a,b
		real::r
		r=0
		do i=1,dim
			r=r+(a(i)-b(i))*(a(i)-b(i))
		end do
		r=sqrt(r)
	end function

	function cluster_distance(cl,i,j) result(r)
		class(cluster)::cl
		integer::dim,i,j
		real(kind=8)::r
		r=-1.0
		if(.not.cl%IsInit) return
		if(i>cl%site.or.j>cl%site) return
		dim=cl%dim
		r=cluster_site_distance(dim,cl%coordinate(1:dim,i),cl%coordinate(1:dim,j))
	end function

	subroutine cluster_delete(cl)
		class(cluster)::cl
		if(.not.cl%IsInit) return
		deallocate(cl%coordinate)
		deallocate(cl%vector)
		deallocate(cl%hop)
		cl%IsInit=.false.
	end subroutine

	subroutine cluster_hopping(cl)
		class(cluster)::cl
		integer::site,dim,h,i,j,l,s,r,q
		real,allocatable,dimension(:)::a,b,c
		real::p,k
		integer,pointer,dimension(:,:,:)::hop
		site=cl%site
		dim=cl%dim
		allocate(a(dim))
		allocate(b(dim))
		allocate(c(dim))
		allocate(hop(site,site,0:8))
		hop(1:site,1:site,0:8)=0
		do i=1,site
			a=cl%coordinate(1:dim,i)
			do j=1,site
				b=cl%coordinate(1:dim,j)
				do l=0,8
					c=b-a+cl%vector(1:dim,l)
					k=dot_product(c,c)
					select case(round(k))
						case(1)
							q=1
						case(2)
							q=0
						case default
							q=0
							cycle
					end select
					hop(j,i,l)=q
				end do
			end do
		end do
		cl%hop=>hop
		cl%translation=h
	end subroutine
	
	function round(r) result(x)
		real::r
		integer::x
		x=int(r);
		if(r-x>0.5) x=x+1
	end function

	subroutine cluster_init(cl,filename)
		class(cluster)::cl
		character(len=*)::filename
		real(8),allocatable::tr(:,:)
		integer::site,dim,i,j,k
		open(unit=7,file=filename,status='old',action='read')
		read(7,*)
		read(7,*) dim
		read(7,*)
		read(7,*) site
		cl%dim=dim
		cl%site=site
		allocate(cl%coordinate(dim,site))
		allocate(tr(dim,3))
		allocate(cl%vector(dim,0:8))
		read(7,*)
		do i=1,site
			read(7,*) cl%coordinate(1:dim,i)
		end do
		read(7,*)
		do i=1,2
			read(7,*) tr(1:dim,i)
		end do
		k=0
		cl%vector(1,0)=0.d0
		cl%vector(2,0)=0.d0
		do i=-1,1
			do j=-1,1
				if(i==0.and.j==0) cycle
				k=k+1
				tr(1:dim,3)=tr(1:dim,1)*i+tr(1:dim,2)*j
				cl%vector(1:dim,k)=tr(1:dim,3)
			end do
		end do
		close(7)
		call cluster_hopping(cl)
		cl%IsInit=.true.
	end subroutine

end module hop_mod

