!--------------------------------------------------------------------------------------------------------------
subroutine COPYATOMS_QEq(imode, dr, atype, pos, v, f, q)
use atoms
!use atoms, only ne, target_node;
!--------------------------------------------------------------------------------------------------------------
implicit none

integer,intent(IN) :: imode
real(8),intent(IN) :: dr(3)
real(8),target :: atype(NBUFFER), q(NBUFFER)
real(8),target :: pos(NBUFFER,3),v(NBUFFER,3),f(NBUFFER,3)

integer :: commflag(0:NBUFFER)

integer :: ixyz,tn1,tn2, dflag
integer :: ni, ity

integer :: ti,tj,tk,tti,ttj

! index of the 6 way communication
integer, parameter :: dinv(6)=(/2,1,4,3,6,5/)
! index to the copyptr. Look up the number of resident atoms before i-th send operation.
integer, parameter :: cptridx(6)=(/0,0,2,2,4,4/)
! tells the direction in which the information needs to be send 
integer, parameter :: is_xyz(6)=(/1,1,2,2,3,3/)  !<- [xxyyzz]

!--- clear total # of copied atoms, sent atoms, recieved atoms
na=0; ns=0; nr=0
copyptr(0)=NATOMS


!--- set Nr of elements during this communication. 
if(imode==MODE_QCOPY1) ne=2
if(imode==MODE_QCOPY2) ne=3

!--- Get the normalized local coordinate will be used through this function
call xu2xs_inplace(max(copyptr(6),NATOMS),pos)

do dflag=1, 6
   ! the node number send to
   tn1 = target_node(dflag)
   ! the node number receive from
   tn2 = target_node(dinv(dflag))
   ixyz = is_xyz(dflag)

   call step_preparation_qeq(dflag, dr, commflag)
   call store_atoms_qeq(tn1, dflag, imode)
   call i_send_recv_qeq(tn1, tn2, myparity(ixyz),dflag)
   call append_atoms_qeq(dflag, imode)
enddo

!--- by here, we got new atom positions in the normalized coordinate, need to update real coordinates.
call xs2xu_inplace(copyptr(6),pos)
return

CONTAINS

!--------------------------------------------------------------------------------------------------------------
! tag the atoms to be sent
subroutine step_preparation_qeq(dflag, dr, commflag)
implicit none
!--------------------------------------------------------------------------------------------------------------
! direction flag
integer, intent(in) :: dflag
! the thickness of the skin from the surface
real(8), intent(in) :: dr(3)
! commflag tag array
integer, intent(inout) :: commflag(0:NBUFFER)
integer :: i
commflag(0)=0
do i=1, copyptr(cptridx(dflag))
   ! check if an atom is within the skin
   ! Ye: FIXME commflag should not be just flags, it should be the index of the filterred atoms
   if (inBuffer_qeq(dflag,dr,pos(i,is_xyz(dflag)))) then
      commflag(0) = commflag(0) + 1
      commflag(commflag(0)) = i
   end if
enddo
return
end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine store_atoms_qeq(tn, dflag, imode)
use atoms;
!--------------------------------------------------------------------------------------------------------------
implicit none
! Ye: why not used
integer,intent(IN) :: tn
integer,intent(IN) :: dflag
integer,intent(IN) :: imode
integer :: n,ni,is,a,b,ioffset,atom_id
real(8) :: sft

!--- reset the number of atoms to be sent
ns=0
!--- # of elements to be sent. should be more than enough. 
! Ye: FIXME it seems computing the reserved space. probably more than required.
ni = commflag(0)*ne
!--- <sbuffer> will deallocated in store_atoms
call CheckSizeThenReallocate_qeq(sbuffer,ni)

do atom_id=1,commflag(0)
   n = commflag(atom_id)
   if (imode==MODE_QCOPY1) then
       sbuffer(atom_id) = qs(n)
       sbuffer(commflag(0)+atom_id) = qt(n)
   end if
   if(imode==MODE_QCOPY2) then
      sbuffer(atom_id) =  hs(n)
      sbuffer(commflag(0)+atom_id) =  ht(n)
      sbuffer(2*commflag(0)+atom_id) = q(n)
   end if
      ns = ns + ne
end do

end subroutine store_atoms_qeq

!--------------------------------------------------------------------------------------------------------------
subroutine  append_atoms_qeq(dflag, imode)
use atoms;
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, imode
integer :: m, i, ine,st,en

if( (na+nr)/ne > NBUFFER) then
    print'(a,i4,5i8)', "ERROR: over capacity in append_atoms; myid,na,nr,ne,(na+nr)/ne,NBUFFER: ", &
          myid,na,nr,ne,(na+nr)/ne, NBUFFER
    call MPI_FINALIZE(ierr)
    stop
endif

copyptr(dflag) = copyptr(dflag-1) + nr/ne

st=copyptr(dflag-1) + 1 
en=copyptr(dflag-1) + 1 + nr/ne
m=nr/ne

if (imode==MODE_QCOPY1) then
    qs(st:en) = rbuffer(1:m)
    qt(st:en) = rbuffer(m+1:2*m)
end if

if (imode==MODE_QCOPY2) then
   hs(st:en) = rbuffer(1:m)
   ht(st:en) = rbuffer(m+1:2*m)
   q(st:en)  = rbuffer(2*m+1:3*m)
end if

!--- update the total # of transfered elements.
na=na+nr

end subroutine append_atoms_qeq

!--------------------------------------------------------------------------------------------------------------
subroutine  i_send_recv_qeq(tn1, tn2, mypar,dflag)
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) ::tn1, tn2, mypar, dflag
integer :: recv_stat(MPI_STATUS_SIZE)
real(8) :: recv_size
integer :: send_request,recv_request

!--- if myid is the same of target-node ID, don't use MPI call.
!--- Just copy <sbuffer> to <rbuffer>. Because <send_recv()> will not be used,
!--- <nr> has to be updated here for <append_atoms()>.

ns = ns_atoms(dflag)*ne

if(myid==tn1) then
   if(ns>0) then
      nr=ns
      call CheckSizeThenReallocate_qeq(rbuffer,nr)
      rbuffer(1:ns) = sbuffer(1:ns)
   else
      nr=0
   endif

   return
endif

call system_clock(ti,tk)

if (ns >0) then
    call MPI_Isend(sbuffer,ns,MPI_DOUBLE_PRECISION,tn1,10,MPI_COMM_WORLD,send_request,ierr)
else
    call MPI_Isend(1,1,MPI_DOUBLE_PRECISION,tn1,10,MPI_COMM_WORLD,send_request,ierr)
end if

nr = nr_atoms(dflag)*ne
call CheckSizeThenReallocate_qeq(rbuffer,nr)
call MPI_Irecv(rbuffer,nr,MPI_DOUBLE_PRECISION,tn2,10,MPI_COMM_WORLD,recv_request,ierr)
if(nr==1) nr=0

call MPI_Wait(recv_request,recv_stat,ierr)

call system_clock(tj,tk)
it_timer(25)=it_timer(25)+(tj-ti)

end subroutine

!--------------------------------------------------------------------------------------------------------------
! check if an atom is with the box to be communicated
function inBuffer_qeq(dflag, dr, rr) result(isInside)
use atoms, only: lbox
!--------------------------------------------------------------------------------------------------------------
implicit none
integer, intent(IN) :: dflag
real(8), intent(IN) :: dr(3), rr
logical :: isInside

select case(dflag)
   case(1) 
      isInside = lbox(1) - dr(1) < rr
   case(2) 
      isInside = rr <= dr(1)
   case(3) 
      isInside = lbox(2) - dr(2) < rr
   case(4) 
      isInside = rr <= dr(2)
   case(5) 
      isInside = lbox(3) - dr(3) < rr
   case(6) 
      isInside = rr <= dr(3)
   case default
      write(6,*) "ERROR: no matching directional flag in isInside: ", dflag
end select

end function

!--------------------------------------------------------------------------------------------------------------
subroutine CheckSizeThenReallocate_qeq(buffer,nsize)
implicit none
!--------------------------------------------------------------------------------------------------------------
real(8),allocatable :: buffer(:)
integer,intent(IN) :: nsize
integer :: ast

if(allocated(buffer)) then
   if(nsize > size(buffer)) then
       deallocate(buffer)
       allocate(buffer(2*nsize), stat=ast)
   endif
else
   allocate(buffer(nsize), stat=ast)
endif

end subroutine

end subroutine COPYATOMS_QEq

