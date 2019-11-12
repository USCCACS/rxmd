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

logical :: commflag(NBUFFER)

integer :: ixyz,tn1,tn2, dflag
integer :: ni, ity

integer :: ti,tj,tk,tti,ttj

! index of the 6 way communication
integer, parameter :: dinv(6)=(/2,1,4,3,6,5/)
! number of resident atoms after i-th receive operation. 0 refers the count before any communication
integer :: natom_resident_receive(0:6)
! index to the natom_resident_receive. Look up the number of resident atoms before i-th send operation.
integer, parameter :: cptridx(6)=(/0,0,2,2,4,4/)
! tells the direction in which the information needs to be send 
integer, parameter :: is_xyz(6)=(/1,1,2,2,3,3/)  !<- [xxyyzz]

!--- clear total # of copied atoms, sent atoms, recieved atoms
na=0; ns=0; nr=0
natom_resident_receive(0)=NATOMS


!--- set Nr of elements during this communication. 
if(imode==MODE_QCOPY1) ne=2
if(imode==MODE_QCOPY2) ne=3

do dflag=1, 6
   ! the node number send to
   tn1 = target_node(dflag)
   ! the node number receive from
   tn2 = target_node(dinv(dflag))
   ixyz = is_xyz(dflag)

   call step_preparation_qeq(dflag, dr, commflag)
   call store_atoms_qeq(tn1, dflag, imode)
   call send_recv_qeq(tn1, tn2, myparity(ixyz))
   call append_atoms_qeq(dflag, imode)
enddo

return

CONTAINS

!--------------------------------------------------------------------------------------------------------------
! tag the atoms to be sent
! @param dflag direction flag
! @param dr ?
! @param commflag tag array
subroutine step_preparation_qeq(dflag, dr, commflag)
implicit none
!--------------------------------------------------------------------------------------------------------------
integer,intent(in) :: dflag
real(8),intent(in) :: dr(3)
logical,intent(inout) :: commflag(NBUFFER)
integer :: i
do i=1, natom_resident_receive(cptridx(dflag))
   commflag(i) = inBuffer_qeq(dflag,dr,pos(i,is_xyz(dflag)))
enddo
return
end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine store_atoms_qeq(tn, dflag, imode)
use atoms;
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: tn, dflag, imode
integer :: n,ni,is,a,b,ioffset
real(8) :: sft

!--- reset the number of atoms to be sent
ns=0
!--- # of elements to be sent. should be more than enough. 
ni = natom_resident_receive(cptridx(dflag))*ne
!--- <sbuffer> will deallocated in store_atoms
call CheckSizeThenReallocate_qeq(sbuffer,ni)

do n=1,natom_resident_receive(cptridx(dflag))
   if(commflag(n)) then
       if (imode==MODE_QCOPY1) then
           sbuffer(ns+1) = qs(n)
           sbuffer(ns+2) = qt(n)
       end if
       if(imode==MODE_QCOPY2) then
           sbuffer(ns+1) =  hs(n)
           sbuffer(ns+2) =  ht(n)
           sbuffer(ns+3) = q(n)
        end if
        ns = ns + ne
     end if
end do

end subroutine store_atoms_qeq

!--------------------------------------------------------------------------------------------------------------
subroutine  append_atoms_qeq(dflag, imode)
use atoms;
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, imode
integer :: m, i, ine

if( (na+nr)/ne > NBUFFER) then
    print'(a,i4,5i8)', "ERROR: over capacity in append_atoms; myid,na,nr,ne,(na+nr)/ne,NBUFFER: ", &
          myid,na,nr,ne,(na+nr)/ne, NBUFFER
    call MPI_FINALIZE(ierr)
    stop
endif

if(nr==0) then
  natom_resident_receive(dflag) = natom_resident_receive(dflag-1)
else
  natom_resident_receive(dflag) = natom_resident_receive(dflag-1) + nr/ne
end if

if (imode==MODE_QCOPY1) then
    do i=0, nr/ne-1
       ine=i*ne
       m = natom_resident_receive(dflag-1) + 1 + i
       qs(m) = rbuffer(ine+1)
       qt(m) = rbuffer(ine+2)
    end do
end if

if (imode==MODE_QCOPY2) then
   do i=0, nr/ne-1
      ine=i*ne
      m = natom_resident_receive(dflag-1) + 1 + i
      hs(m) = rbuffer(ine+1)
      ht(m) = rbuffer(ine+2)
      q(m)  = rbuffer(ine+3)
   end do
end if

!--- update the total # of transfered elements.
na=na+nr

end subroutine append_atoms_qeq

!--------------------------------------------------------------------------------------------------------------
subroutine  send_recv_qeq(tn1, tn2, mypar)
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) ::tn1, tn2, mypar
integer :: recv_stat(MPI_STATUS_SIZE)
real(8) :: recv_size

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

if (mypar == 0) then

     ! the number of elements per data packet has to be greater than 1, for
     ! example NE_COPY = 10.
     ! if ns == 0, send one double to tell remote rank that there will be no
     ! atom data to be sent. 
     if (ns > 0) then
       call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn1, 10, MPI_COMM_WORLD,ierr)
     else
       call MPI_SEND(1, 1, MPI_DOUBLE_PRECISION, tn1, 10, MPI_COMM_WORLD, ierr)
     endif

     call MPI_Probe(tn2, 11, MPI_COMM_WORLD, recv_stat, ierr)
     call MPI_Get_count(recv_stat, MPI_DOUBLE_PRECISION, nr, ierr)

     call CheckSizeThenReallocate_qeq(rbuffer,nr)

     call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn2, 11, MPI_COMM_WORLD,recv_stat, ierr)

     ! the number of elements per data packet has to be greater than 1, for
     ! example NE_COPY = 10.
     ! nr == 1 means no atom data to be received. 
     if(nr==1) nr=0

elseif (mypar == 1) then

     call MPI_Probe(tn2, 10, MPI_COMM_WORLD, recv_stat, ierr)
     call MPI_Get_count(recv_stat, MPI_DOUBLE_PRECISION, nr, ierr)

     call CheckSizeThenReallocate_qeq(rbuffer,nr)

     call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn2, 10, MPI_COMM_WORLD, recv_stat, ierr)

     ! the number of elements per data packet has to be greater than 1, for
     ! example NE_COPY = 10.
     ! nr == 1 means no atom data to be received. 
     if(nr==1) nr=0

     if (ns > 0) then
        call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn1, 11, MPI_COMM_WORLD, ierr)
     else
        call MPI_SEND(1, 1, MPI_DOUBLE_PRECISION, tn1, 11, MPI_COMM_WORLD, ierr)
     endif

endif
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

