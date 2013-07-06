program greenfunc
	use timer
	use hubbard_cluster
	use basis
 	use VCA
	implicit none
	type(hubbard_cluster_type),pointer::clust
	type(timer_type),pointer::t
	real(kind=8)::miu,elec,kx,ky,w,time,hx,hy,x,y
	real(kind=8),pointer,dimension(:)::v
	complex(kind=8)::g,a,b
	integer::i,j,k,ie
	logical::err
	character(len=40)::bc
	complex(kind=8)::Xi

	t=>timer_init()
	clust=>cluster_build(w,miu)
	open(unit=7,file='input.data',status='old',iostat=ie)
	if(ie==0) then
		do while(ie==0)
			read(7,*,iostat=ie)
		end do
		ie=1
		backspace(7)
		backspace(7)
		backspace(7)
		do while(ie/=0)
			!read(7,*,iostat=ie) bc
			!print *,bc
			read(7,*,iostat=ie) clust%tx,clust%ty,clust%miu
			backspace(7)
			backspace(7)
		end do
		close(7)
	end if
	write(*,*)"----------------------------------------------------------------"
	write(*,*)"VCA calculation: initial parameters are"
	write(*,*)"U=",clust%U,"mu=",miu,"t=",1.00
	write(*,*)"Cluster t=",clust%tx,"Cluster mu=",clust%miu
	write(*,*)"----------------------------------------------------------------"
	call VCA_optimal(clust,miu,x,elec,err)
	write(*,*) miu,x,elec
	write(*,*) clust%tx,clust%ty,clust%miu
	call VCA_fermi_surface(clust,miu)
	call VCA_spectral_function(clust,miu)
	call VCA_DOS(clust,miu)
	call VCA_trans(clust,miu,1.d-3)
39	write(*,*) "total time =",timer_runtime(t),"seconds"
	call cluster_clean(clust)
	deallocate(t)
end program greenfunc
