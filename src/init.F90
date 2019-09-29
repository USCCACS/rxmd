!------------------------------------------------------------------------------------------
module init
!------------------------------------------------------------------------------------------

  use base
  use atoms
  use mpi_mod
  use qeq_mod, only : QEq
  use pqeq_mod, only : PQEq, initialize_eField, initialize_pqeq
  use force_mod, only : force_reaxff
  use memory_allocator_mod
  use fileio, only : ReadBIN, ReadXYZ

  use fnn, only : fnn_param, fnn_param_obj, get_cutoff_fnn, num_networks_per_atom, &
                  network_ctor, mean_stddev_loader, num_forcecomps, num_pairs, num_types, & 
                  mddriver_fnn, get_max_cutoff

  use lists_mod, only: nbrcut, get_mesh_for_nonbonding_list 

  use fnnin_parser, only : fnn_param_ctor, set_feature_tables_fnn

  use reaxff_param_mod, only : chi, eta, mddriver_reaxff, &
      get_cutoff_bondorder, set_potentialtables_reaxff, get_forcefield_params_reaxff

  use msd_mod, only : msd_data, msd_type_ctor
  use stats_mod, only : stats, stats_type_ctor


contains

!------------------------------------------------------------------------------------------
subroutine mdcontext_base(mdbase, atype, pos, v, f, q)
!------------------------------------------------------------------------------------------
implicit none

type(mdbase_class),intent(in out) :: mdbase 

real(8),intent(in out),allocatable,dimension(:) :: atype, q
real(8),intent(in out),allocatable,dimension(:,:) :: pos, v, f

real(8) :: dns
integer :: i,j,k, ity, l(3)

!--- an error trap
if(vprocs(1)*vprocs(2)*vprocs(3) /= nprocs ) then
  if(myid==0) write(6,'(a60,3i3,i5)')  &
  "ERROR: requested/allocated # of procs are not consistent: ", vprocs(1:3), nprocs
  call MPI_FINALIZE(ierr)
  stop
endif

!--- TODO make this stat a class
!--- allocate & initialize Array size Stat variables
call allocator(maxas, 1, (ntime_step/pstep)+1, 1, nmaxas)

!--- setup the vector ID and parity for processes, in x, y and z order.
vID(1)=mod(myid,vprocs(1))
vID(2)=mod(myid/vprocs(1),vprocs(2))
vID(3)=myid/(vprocs(1)*vprocs(2))
myparity(1)=mod(vID(1),2)
myparity(2)=mod(vID(2),2)
myparity(3)=mod(vID(3),2)

k=0
do i=1,3
!--- Add (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1) to vector ID applying PBC. 
   do j=1,-1,-2

!--- save vector ID.
      l(1:3) = vID(1:3)

!--- apply PBC.
      l(i)=mod( (vID(i) + j + vprocs(i)), vprocs(i) )
     
!--- get direction index.
      k=k+1 ![123456] 

!--- convert vector ID to sequential ID.
      target_node(k) = l(1)+l(2)*vprocs(1)+l(3)*vprocs(1)*vprocs(2) 
   enddo         
                   
enddo    

if(is_fnn) then
  fnn_param_obj = mdcontext_fnn()
  mdbase%ff => fnn_param_obj
  if(myid==0) print*,'get_mdcontext_func : mdcontext_fnn'
else
  call mdcontext_reaxff()
  call set_potentialtables_reaxff()
  if(myid==0) print*,'get_mdcontext_func : mdcontext_reaxff'
endif

if(isRunFromXYZ) then
  call ReadXYZ(atype, pos, v, q, f, RunFromXYZPath)
else
  call ReadBIN(atype, pos, v, q, f, trim(DataDir)//"/rxff.bin")
endif

!--- get global number of atoms by summing up each type. 
natoms_per_type = 0
do i=1, NATOMS
   ity=nint(atype(i))
   natoms_per_type(ity)=natoms_per_type(ity)+1
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, natoms_per_type, size(natoms_per_type), &
                   MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

GNATOMS = sum(natoms_per_type)

!--- determine cutoff distances
if(is_fnn) then
!--- set cutoff distance
  call get_cutoff_fnn(rc, rc2, maxrc, get_max_cutoff(fnn_param_obj))
else
!--- get cutoff distance based on the bond-order
  call get_cutoff_bondorder(rc, rc2, maxrc, natoms_per_type)
!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
  call get_mesh_for_nonbonding_list(rctap)
endif

!--- index array for returning reaction force
call allocator(frcindx,1,NBUFFER)

!--- TODO where to put the spring force extention? 
!--- for spring force
if (isSpring) then

  write(6,'(a)') repeat('-',60)
  write(6,'(a,f8.2, $)') 'springConst [kcal/mol] ', springConst 
  do ity=1, size(hasSpringForce) 
     if(hasSpringForce(ity)) print'(i3, $)',ity
  enddo
  print*
  write(6,'(a)') repeat('-',60)
endif

!--- print out parameters and open data file
if(myid==0) then
   write(6,'(a)') repeat('-',60)
   write(6,'(a30,i9,a3,i9)') "req/alloc # of procs:", vprocs(1)*vprocs(2)*vprocs(3), "  /",nprocs
   write(6,'(a30,3i9)')      "req proc arrangement:", vprocs(1),vprocs(2),vprocs(3)
   write(6,'(a30,es12.2)')   "time step[fs]:",dt*UTIME
   write(6,'(a30,i3, i10, i10)') "MDMODE CURRENTSTEP NTIMESTEP:", &
                                  mdmode, current_step, ntime_step
   write(6,'(a30,f12.3,f8.3,i9)') 'treq,vsfact,sstep:',treq*UTEMP0, vsfact, sstep
   write(6,'(a30,2i6)') 'fstep,pstep:', fstep,pstep
   write(6,'(a30,i24,i24)') "NATOMS GNATOMS:", NATOMS, GNATOMS
   write(6,'(a30,3f12.3)') "LBOX:",LBOX(1:3)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,1,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,2,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,3,0)
   print'(a30,3f12.3)', 'lata,latb,latc:', lata,latb,latc
   print'(a30,3f12.3)', 'lalpha,lbeta,lgamma:', lalpha,lbeta,lgamma
   write(6,'(a30,2i9)') "NBUFFER, MAXNEIGHBS:", NBUFFER, MAXNEIGHBS
   write(6,'(a)') repeat('-',60)
   write(6,'(a30,a12)') "DataDir :", trim(DataDir)
   write(6,'(a30,2(a12,1x))') &
         "FFPath, ParmPath:", trim(FFPath),trim(ParmPath)
   write(6,'(a)') repeat('-',60)

endif

!--- To get density, dhtm, hmas, need mass from forcefield
dns = sum(mass*natoms_per_type(1:size(mass)))/mdbox*UDENS

!--- dt/2*mass, mass/2
call allocator(dthm, 1, size(mass))
call allocator(hmas, 1, size(mass))
do ity=1, size(mass)
   if(mass(ity) > 0.d0) then
      dthm(ity) = dt*0.5d0/mass(ity)
      hmas(ity) = 0.5d0*mass(ity)
   endif
enddo

!--- maxrc from forcefiled is necessary to update box-related variables. 
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

!--- Linked List & Near Neighb Parameters. cc is from update_box_params()
call allocator(nbrlist,1,NBUFFER,0,MAXNEIGHBS)
call allocator(llist,1,NBUFFER)
call allocator(header,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)
call allocator(nacell,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)

if(myid==0) then
   write(6,'(a)') repeat('-',60)
   write(6,'(a30,3f10.4)') "density [g/cc]:", dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", &
        maxrc,lata/cc(1)/vprocs(1),latb/cc(2)/vprocs(2),latc/cc(3)/vprocs(3)

   print'(a30, $)','# of atoms per type:'
   do ity=1, size(natoms_per_type)
      if(natoms_per_type(ity)>0) print'(i12,a2,i2, $)',natoms_per_type(ity),' -',ity
   enddo
   print*
   write(6,'(a)') repeat('-',60)
endif

!--- MSD constractor
msd_data = msd_type_ctor(atmname, msd_size=(ntime_step/pstep), unit_time=(pstep*dt*UTIME))
if(msd_data%is_msd) then
  if(myid==0) then
     write(6,fmt='(a)') 'start MSD measurement.'
     write(6,fmt='(a,i6,f8.3)') 'msd_size, unit_time: ', &
                         size(msd_data%dat,dim=2), msd_data%unit_time
  endif
endif

!--- Structure analyzer constractor
!NOTE: cutoff in stats, e.g. g(r), can be the min value of each lattice vectors,
! but use rctap for ReaxFF since the nonboded list is used in the force/energy calcs.
if(is_fnn) then
  nbrcut = minval([lata,latb,latc])
else
  nbrcut = rctap
endif
stats = stats_type_ctor(nbrcut,atmname)

!--- keep initial position
if(isSpring .or. msd_data%is_msd) then

  has_initial_pos = .true.

  call allocator(ipos,1,NBUFFER,1,3)
  ipos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)

  if(myid==0) write(6,fmt='(a,l3)') 'has_initial_pos:', has_initial_pos

endif

end subroutine

!------------------------------------------------------------------------------------------
function mdcontext_fnn() result(fp)
!------------------------------------------------------------------------------------------
implicit none

integer :: i,j, ity, num_models
integer,allocatable :: dims(:)
character(len=:),allocatable :: path
character(len=1) :: xyz_suffix(3)=['x','y','z']

logical :: verbose=.false.

type(fnn_param) :: fp

if((.not.isRunFromXYZ) .and. (myid==0)) verbose=.true.

!FIXME path needs to given from cmdline
fp = fnn_param_ctor(str_gen('fnn.in'))

!--- setup feature calculation tables
call set_feature_tables_fnn(fp)

!--- FNN specific output 
if(myid==0) call fp%print()

num_models = size(fp%models)

num_types = num_models
num_pairs = num_types*(num_types+1)/2

allocate( mass(num_models), atmname(num_models) )
do i=1, num_models
   atmname(i) = fp%models(i)%element
   mass(i) = fp%models(i)%mass
enddo

!print'(a,a3,f8.3,2i6)','atmname, mass: ', &
!       atmname, mass, size(atmname), size(mass)

do i=1, num_models
   associate ( m => fp%models(i) )
      allocate(m%networks(num_networks_per_atom))
      allocate(m%fstat(num_networks_per_atom))
   
      path = str_gen('FNN/'//m%element//'/')
      dims = m%layersize

      do j=1, num_networks_per_atom
         m%networks(j)  = network_ctor(dims, path, xyz_suffix(j), verbose)
         call mean_stddev_loader(m%fstat(j)%mean, m%fstat(j)%stddev, &
                                 dims(1), path, xyz_suffix(j), verbose)
      enddo

!--- features(3, num_features, natoms)
     call allocator(m%features, 1,3, 1,dims(1), 1,NBUFFER) 

   end associate
   !print*,i, ' : ', atmname(i), mass(i), size(fp%models(i)%networks)
enddo

!--- set md dirver function 
mddriver_func => mddriver_fnn


end function

!------------------------------------------------------------------------------------------
subroutine mdcontext_reaxff()
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
implicit none

integer :: i,j,k, ity

wt0 = MPI_WTIME()

!--- allocate variales used during reaxff energy & force calculations
call mdvariables_allocator_reaxff()

!--- set charge model
if(isPQEq) then
  charge_model_func => PQEq
  rctap = rctap0_pqeq

  call allocator(spos,1,NBUFFER,1,3)

  call initialize_pqeq(chi,eta)
  if(isEfield) call initialize_eField(myid)

else
  charge_model_func => QEq
  rctap = rctap0
endif

!--- set taper function for the vdw and coulomb terms
rctap2 = rctap**2
CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(rctap)**4, &
          84.d0/(rctap)**5, -70.d0/(rctap)**6, &
          20.d0/(rctap)**7 /)

!--- set force model
force_model_func => force_reaxff

!--- set md dirver function 
mddriver_func => mddriver_reaxff

!--- set force field parameters
call get_forcefield_params_reaxff(FFPath)

!--- ReaxFF specific output 
if(myid==0) then

   write(6,'(a)') repeat('-',60)
   write(6,'(a30,2i6)') "NMINCELL, MAXNEIGHBS10:", NMINCELL, MAXNEIGHBS10
   write(6,'(a30,i6,es10.1,i6,i6)') "isQEq,QEq_tol,NMAXQEq,qstep:", &
                                     isQEq,QEq_tol,NMAXQEq,qstep
   write(6,'(a30,f8.3,f8.3)') 'Lex_fqs,Lex_k:',Lex_fqs,Lex_k
   write(6,'(a)') repeat('-',60)

   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

end subroutine

end module
