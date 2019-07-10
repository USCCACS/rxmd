!-------------------------------------------------------------------------------------------
module cmdline_args
use atoms
implicit none
!-------------------------------------------------------------------------------------------

interface get_token_and_set_value
   module procedure set_r8, set_i4, set_l, set_str
end interface

contains

subroutine set_r8(linein, val)
   character(len=:),allocatable,intent(in out) :: linein
   real(8),intent(in out) :: val
   character(len=:),allocatable :: token

   if(getstr(linein, token) < 0) &
     stop 'Error:  token not found in set_r8'
    
   read(token,*) val
end subroutine

subroutine set_i4(linein, val)
   character(len=:),allocatable,intent(in out) :: linein
   integer,intent(in out) :: val
   character(len=:),allocatable :: token

   if(getstr(linein, token) < 0) &
     stop 'Error:  token not found in set_i4'

   read(token,*) val
end subroutine

subroutine set_l(linein, val)
   character(len=:),allocatable,intent(in out) :: linein
   logical,intent(in out) :: val
   character(len=:),allocatable :: token

   if(getstr(linein, token) < 0) &
     stop 'Error:  token not found in set_l'

   read(token,*) val
end subroutine

subroutine set_str(linein, val)
   character(len=:),allocatable,intent(in out) :: linein
   character(MAXSTRLENGTH),intent(in out) :: val
   character(len=:),allocatable :: token

   if(getstr(linein, token) < 0) &
     stop 'Error:  token not found in set_str'

   read(token,'(a)') val
end subroutine

!-------------------------------------------------------------------------------------------
subroutine get_cmdline_args(myrank, eFieldDir, eFieldStrength)
!-------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: myrank
integer :: eFieldDir
real(8) :: eFieldStrength 

integer :: ity,idx
character(MAXSTRLENGTH) :: argv

if(command_argument_count()>0) then
  if(myid==0) then
    write(6,'(a)') repeat('-',60)
    write(6,'(a $)') 'commandline args : '
    do idx=1, command_argument_count()
       call get_command_argument(idx,argv)
       write(6,'(a,a1,$)') trim(adjustl(argv)), ' '
    enddo 
    write(6,*)
    write(6,'(a)') repeat('-',60)
  endif
endif

!--- read FF file, output dir, MD parameter file paths from command line

if(find_cmdline_argc('--help',idx).or.find_cmdline_argc('-h',idx)) then
    if(myrank==0) print'(a)', "usage : ./rxmd --ffield ffield --outDir DAT --rxmdin rxmd.in"
    call MPI_FINALIZE(ierr)
    stop
endif

if(find_cmdline_argc('--rxmdin',idx).or.find_cmdline_argc('-in',idx)) then
    call get_command_argument(idx+1,argv)
    ParmPath=trim(adjustl(argv))
endif
call get_rxmd_parms(ParmPath)

if(find_cmdline_argc('--ffield',idx).or.find_cmdline_argc('-ff',idx)) then
    call get_command_argument(idx+1,argv)
    FFPath=trim(adjustl(argv))
endif

if(find_cmdline_argc('--outDir',idx).or.find_cmdline_argc('-o',idx)) then
    call get_command_argument(idx+1,argv)
    DataDir=trim(adjustl(argv))
endif

if(find_cmdline_argc('--run_from_xyz',idx).or.find_cmdline_argc('-runxyz',idx)) then
    call get_command_argument(idx+1,argv)
    isRunFromXYZ=.true.
    RunFromXYZPath=trim(adjustl(argv))
endif

!--- PQEq parameters
if(find_cmdline_argc('--pqeq',idx).or.find_cmdline_argc('-pqeq',idx)) then
    isPQEq = .true.
    call get_command_argument(idx+1,argv)
    PQEqParmPath=trim(adjustl(argv))

endif

if(isPQEq) then
    if(myrank==0) then
       print'(a60)',repeat('-',60)
       print'(2a30)','Enabling PQEq: PQEqParmPath: ', trim(adjustl(PQEqParmPath))
       print'(a60)',repeat('-',60)
    endif

    call get_pqeq_parms(PQEqParmPath)
endif


!--- electric field 
if(find_cmdline_argc('--efield',idx).or.find_cmdline_argc('-e',idx)) then
    isEfield=.true.
    call get_command_argument(idx+1,argv)
    read(argv,*) eFieldDir
    call get_command_argument(idx+2,argv)
    read(argv,*) eFieldStrength
endif

if(isEfield .and. myrank==0) then
   print'(a60)',repeat('-',60)
   print'(a30)','Enabling electric field : '
   print'(a30,i3,f10.5)', 'eField [V/A], eFiled direction : ', eFieldDir, eFieldStrength
   print'(a60)',repeat('-',60)
endif



if(find_cmdline_argc('--lg',idx).or.find_cmdline_argc('-lg',idx)) then
    if(myrank==0) print'(a30)','Enabling LG term'
    isLG=.true.
endif

if(find_cmdline_argc('--spring',idx).or.find_cmdline_argc('-s',idx)) then
    isSpring=.true.
    call get_command_argument(idx+1,argv)
    read(argv,*) springConst
    if(myrank==0) then
       print'(a60)',repeat('-',60)
       print'(a30)','Enabling spring force: springConst: '
       print'(f10.5)', springConst 
       print'(a60)',repeat('-',60)
    endif
endif

end subroutine

!-------------------------------------------------------------------------------------------
subroutine get_pqeq_parms(PQEqParmPath)
use MemoryAllocator
implicit none
!-------------------------------------------------------------------------------------------
character(MAXSTRLENGTH),intent(in) :: PQEqParmPath
character(MAXSTRLENGTH) :: linein0
character(len=:),allocatable :: linein, token 
integer :: n,nparms, istat

open(1,file=PQEqParmPath,status='old')

n=0; nparms=1

do while (n<nparms)
  read(1,'(a)',end=10) linein0
  linein = trim(adjustl(linein0))

  if(linein(1:1)=='#') cycle

  if(index(linein, 'NPARMS') /= 0) then
    istat=getstr(linein, token)

    istat=getstr(linein, token)
    read(token,*) nparms

    ntype_pqeq=nparms
    ntype_pqeq2=ntype_pqeq**2

    allocate(character(2) :: Elempqeq(ntype_pqeq))
    allocate(isPolarizable(ntype_pqeq))
    allocate(X0pqeq(ntype_pqeq), J0pqeq(ntype_pqeq))
    allocate(Rcpqeq(ntype_pqeq), Rspqeq(ntype_pqeq))
    allocate(Zpqeq(ntype_pqeq), Kspqeq(ntype_pqeq))

    allocate(alphacc(ntype_pqeq,ntype_pqeq))
    allocate(alphasc(ntype_pqeq,ntype_pqeq))
    allocate(alphass(ntype_pqeq,ntype_pqeq))

    if(myid==0) then
       print'(a40)','n,name,polarizable,X0,J0,Z,Rs,Rc,Ks: '
       print'(a60)',repeat('-',60)
    endif

    cycle
  endif

  n=n+1

  istat=getstr(linein, token); Elempqeq(n) = trim(adjustl(token))
  istat=getstr(linein, token); isPolarizable(n) = .true.
  istat=getstr(linein, token); read(token,*) X0pqeq(n)
  istat=getstr(linein, token); read(token,*) J0pqeq(n)
  istat=getstr(linein, token); read(token,*) Zpqeq(n)
  istat=getstr(linein, token); read(token,*) Rcpqeq(n)
  istat=getstr(linein, token); read(token,*) Rspqeq(n)
  istat=getstr(linein, token); read(token,*) Kspqeq(n) 

  if(myid==0) & 
  print'(i3,a3,l3,6f12.5)',n,Elempqeq(n),isPolarizable(n), &
      X0pqeq(n), J0pqeq(n), Zpqeq(n), Rcpqeq(n), Rspqeq(n), Kspqeq(n)

end do

10 close(1)

if(myid==0) print'(a60)',repeat('-',60)

return
end

!-------------------------------------------------------------------------------------------
subroutine get_rxmd_parms(rxmdParmPath)
implicit none
!-------------------------------------------------------------------------------------------
character(MAXSTRLENGTH),intent(in) :: rxmdParmPath
character(MAXSTRLENGTH) :: argv, linein0
character(len=:),allocatable :: linein, token
integer :: idx

open(1, file=trim(ParmPath), status="old")

do while (.true.)
  read(1,'(a)',end=10) linein0
  linein = trim(adjustl(linein0))

  if(getstr(linein, token) > 0) then

    select case (token)
      case ('mdmode')
         call get_token_and_set_value(linein, mdmode)
      case ('time')
         call get_token_and_set_value(linein, dt)
         call get_token_and_set_value(linein, ntime_step)
      case ('temperature')
         call get_token_and_set_value(linein, treq)
         call get_token_and_set_value(linein, vsfact)
         call get_token_and_set_value(linein, sstep)
      case ('io_step')
         call get_token_and_set_value(linein, fstep)
         call get_token_and_set_value(linein, pstep)
      case ('io_type')
         call get_token_and_set_value(linein, isBinary)
         call get_token_and_set_value(linein, isBondFile)
         call get_token_and_set_value(linein, isPDB)
         call get_token_and_set_value(linein, isXYZ)
      case ('processors')
         call get_token_and_set_value(linein, vprocs(1))
         call get_token_and_set_value(linein, vprocs(2))
         call get_token_and_set_value(linein, vprocs(3))
      case ('QEq')
         call get_token_and_set_value(linein, isQEq)
         call get_token_and_set_value(linein, NMAXQEq)
         call get_token_and_set_value(linein, QEq_tol)
         call get_token_and_set_value(linein, qstep)
      case ('exL')
         call get_token_and_set_value(linein, Lex_fqs)
         call get_token_and_set_value(linein, Lex_k)
      case ('CG_tol')
         call get_token_and_set_value(linein, ftol)
      case ('efield')
         isEfield=.true.
         call get_token_and_set_value(linein, eFieldDir)
         call get_token_and_set_value(linein, eFieldStrength)
      case ('PQEqParm')
         isPQEq=.true.
         call get_token_and_set_value(linein, PQEqParmPath)
      case default
         if(myid==0) print*,'ERROR: '//trim(token)//' is not found'
         stop
    end select

    !--- goto next line
    cycle
  endif

end do

10 close(1)

if(find_cmdline_argc('--mdmode',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) mdmode
endif

if(find_cmdline_argc('--dt',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) dt
endif

if(find_cmdline_argc('--ntime_step',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) ntime_step
endif

if(find_cmdline_argc('--treq',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) treq
endif

if(find_cmdline_argc('--vsfact',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) vsfact
endif

if(find_cmdline_argc('--sstep',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) sstep
endif

if(find_cmdline_argc('--fstep',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) fstep
endif

if(find_cmdline_argc('--pstep',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) pstep
endif

if(find_cmdline_argc('--vprocs',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) vprocs(1)
    call get_command_argument(idx+2,argv)
    read(argv,*) vprocs(2)
    call get_command_argument(idx+3,argv)
    read(argv,*) vprocs(3)
endif

if(find_cmdline_argc('--isQEq',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*)isQEq 
endif

if(find_cmdline_argc('--NMAXQEq',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) NMAXQEq
endif

if(find_cmdline_argc('--QEq_tol',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) QEq_tol
endif

if(find_cmdline_argc('--qstep',idx)) then
    call get_command_argument(idx+1,argv)
    read(argv,*) qstep
endif

if(find_cmdline_argc('--isBinary',idx)) isBinary=.true. 
if(find_cmdline_argc('--isBondFile',idx)) isBondFile=.true. 
if(find_cmdline_argc('--isPDB',idx)) isPDB=.true. 
if(find_cmdline_argc('--isXYZ',idx)) isXYZ=.true. 

end subroutine

!-------------------------------------------------------------------------------------------
integer function getstr(linein,lineout)
implicit none
!-------------------------------------------------------------------------------------------

character(len=:),allocatable,intent(in out) :: linein,lineout
integer :: pos1

!--- remove whitespace 
linein = adjustl(linein)

!--- return if black line
if(len(linein)==0) then
  getstr=-2
  return
endif

!--- return if it's a comment line or entirely whitespace
if(linein(1:1)=='#' .or. linein == repeat(' ', len(linein)) ) then
   getstr=-1
   return
endif

! find position in linein to get a token. if whitespace is not found, take entire line
pos1=index(linein,' ')-1
if(pos1==-1) pos1=len(linein)

lineout=linein(1:pos1)
linein=linein(pos1+1:)
getstr=len(lineout)

return
end

!-------------------------------------------------------------------------------------------
logical function find_cmdline_argc(key,idx)
implicit none
!-------------------------------------------------------------------------------------------
integer,intent(inout) :: idx
character(*) :: key

integer :: i
character(MAXSTRLENGTH) :: argv

do i=1, command_argument_count()
   call get_command_argument(i,argv)
   if(index(argv,trim(adjustl(key))//' ')/=0) then ! trailing zero to distinguish '-foo ' and '-fooo'
      idx=i
      find_cmdline_argc=.true.
      return
   endif
enddo

idx=-1
find_cmdline_argc=.false.

return
end function

end module
