module hop_mod
	implicit none

	type,public::hop
		integer::site,dim,translation,nlt,nct
		real,pointer,dimension(:,:)::coordinate
		real,pointer,dimension(:,:)::vector
		integer,pointer,dimension(:,:)::lattice
		integer,pointer,dimension(:,:)::cluster
		real(8),pointer,dimension(:)::lt,ct
		real(8)::lmu,cmu,U,Tep,M
	end type hop

	contains

	subroutine cluster_site(cl,n,a)
		type(hop)::cl
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
		type(hop)::cl
		integer::dim,i,j
		real(kind=8)::r
		r=-1.0
		if(i>cl%site.or.j>cl%site) return
		dim=cl%dim
		r=cluster_site_distance(dim,cl%coordinate(1:dim,i),cl%coordinate(1:dim,j))
	end function

	subroutine cluster_delete(cl)
		type(hop),pointer::cl
		deallocate(cl%coordinate)
		deallocate(cl%vector)
		deallocate(cl%lattice)
		deallocate(cl%cluster)
		deallocate(cl%lt)
		deallocate(cl%ct)
		deallocate(cl)
	end subroutine

	function cluster_init() result(cl)
		type(hop),pointer::cl
		character(len=20)::filename,rc
		real(8),allocatable::tr(:,:)
		integer::site,dim,i,j,k,ierr
		allocate(cl)
		open(unit=8,file='cluster.input',status='old',action='read')
		read(8,*) rc
		filename=trim(adjustl(rc))//'.conf'
		open(unit=7,file=trim(adjustl(filename)),status='old',action='read')
		read(7,*)
		read(7,*) dim
		read(7,*)
		read(7,*) site
		cl%dim=dim
		cl%site=site
		allocate(cl%coordinate(dim,site))
		allocate(cl%vector(dim,dim))
		read(7,*)
		do i=1,site
			read(7,*) cl%coordinate(1:dim,i)
		end do
		read(7,*)
		do i=1,dim
			read(7,*) cl%vector(1:dim,i)
		end do
		close(7)
		read(8,*)
		read(8,*)
		read(8,*)
		i=0
		do
			i=i+1
			read(8,*,iostat=ierr) j
			if(ierr/=0) exit
		end do
		do j=1,i+1
			backspace(8)
		end do
		allocate(cl%lattice(6,i-1))
		k=0
		cl%nlt=i-1
		do j=1,cl%nlt
			read(8,*) cl%lattice(1:6,j)
			if(cl%lattice(3,j)>k) k=cl%lattice(3,j)
		end do
		read(8,*)
		allocate(cl%lt(k))
		do i=1,k
			read(8,*) rc,cl%lt(i)
		end do
		read(8,*)
		read(8,*) rc,cl%lmu
		read(8,*) rc,cl%U
		read(8,*) rc,cl%Tep
		read(8,*)
		read(8,*)
		read(8,*)
		i=0
		do
			i=i+1
			read(8,*,iostat=ierr) j
			if(ierr/=0) exit
		end do
		do j=1,i+1
			backspace(8)
		end do
		allocate(cl%cluster(4,i-1))
		k=0
		cl%nct=i-1
		do j=1,cl%nct
			read(8,*) cl%cluster(1:4,j)
			if(cl%cluster(3,j)>k) k=cl%cluster(3,j)
		end do
		read(8,*)
		allocate(cl%ct(k))
		do i=1,k
			read(8,*) rc,cl%ct(i)
		end do
		read(8,*)
		read(8,*) rc,cl%cmu
		close(8)
		cl%M=0
	end function

end module hop_mod

