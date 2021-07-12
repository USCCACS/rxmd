!------------------------------------------------------------------------------
#ifdef LIBRXMD
subroutine rxmd()
#else
program rxmd
#endif
!------------------------------------------------------------------------------
use base, only : mdbase_class, atype, f, pos, q, v 
use init
use mpi_mod
use cmdline_args, only : get_cmdline_args, get_rxmd_parms
!------------------------------------------------------------------------------
implicit none
integer :: i,ity,it1,it2,irt

type(mdbase_class) :: mdbase

integer :: resultlen
character(len=8) :: hostname

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!call mpi_get_processor_name(hostname, resultlen, ierr) 
!print'(a,1x,i9,1x,a)','INFO: myid,hostname ', myid, hostname

if(myid==0)  print'(a30)', 'rxmd has started'

!--- process command line arguments
call get_cmdline_args(myid, eFieldDir, eFieldStrength)

!--- get forcefield-independent parameters 
call get_rxmd_parms(ParmPath)

!--- initialize the MD system
call mdcontext_base(mdbase, atype, pos, v, f, q)

!--- Enter Main MD loop 
call system_clock(it1,irt)
call mddriver_func(mdbase, ntime_step)
call system_clock(it2,irt)
it_timer(MAXTIMERS)=(it2-it1)

if(myid==0)  print'(a30)', 'rxmd has finished successfully'

call MPI_FINALIZE(ierr)
end 
