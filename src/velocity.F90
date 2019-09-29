module velocity_modifiers_mod

  use mpi_mod
  use utils, only : UTEMP0, UTEMP, matinv
  use base, only : NATOMS, GNATOMS, mass, treq, KE, GKE, myid, ierr, dthm, hmas, atmname

contains

!------------------------------------------------------------------------------------------
subroutine simple_scaling(atype, v, scaling_factor)
!------------------------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)
real(8),intent(in) :: scaling_factor

   v(1:NATOMS,1:3)=scaling_factor*v(1:NATOMS,1:3)

end subroutine

!------------------------------------------------------------------------------------------
subroutine scale_to_target_temperature(atype, v, t_target)
!------------------------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)
real(8),intent(in) :: t_target
real(8) :: ctmp

   ctmp = (t_target*UTEMP0)/(GKE*UTEMP)
   v(1:NATOMS,1:3)=sqrt(ctmp)*v(1:NATOMS,1:3)

end subroutine

!------------------------------------------------------------------------------------------
subroutine gaussian_dist_velocity(atype, v)
! Generate gaussian distributed velocity as an initial value using Box-Muller algorithm
!------------------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)

integer :: i, k, ity
real(8) :: vv(2), vsqr, vsl, rndm(2)
real(8) :: vCM(3), GvCM(3), mm, Gmm
real(8) :: vfactor

!--- assign velocity to two atoms together with BM algoritm. 
!--- If <NATOMS> is odd, the <NATOMS> + 1 element will be the ignored in later calculations.

do i=1, NATOMS, 2

  do k=1,3 ! three directions
     !--- generate gaussian distributed velocity
     vsqr=0.d0
     do while ( (vsqr >= 1.d0) .or. (vsqr==0.d0) ) 
        call random_number(rndm)
        vv(1) = 2.d0 * rndm(1) - 1.d0
        vv(2) = 2.d0 * rndm(2) - 1.d0
        vsqr = vv(1)**2 + vv(2)**2
     enddo

     vsl = sqrt(-2.d0 * log(vsqr)/vsqr)
     v(i,k)   = vv(1)*vsl
     v(i+1,k) = vv(2)*vsl
  enddo
  
enddo

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   vCM(1:3)=vCM(1:3) + mass(ity)*v(i,1:3)
   mm = mm + mass(ity)
enddo
 
call MPI_ALLREDUCE(MPI_IN_PLACE, vCM, size(vCM), MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, mm,1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
GvCM = vCM
Gmm = mm

!--- get the global momentum
GvCM(:)=GvCM(:)/Gmm

!--- set the total momentum to be zero and get the current kinetic energy. 
KE = 0.d0
do i=1, NATOMS
   v(i,1:3) = v(i,1:3) - GvCM(1:3)

   ity = nint(atype(i))
   KE = KE + hmas(ity)*sum( v(i,1:3)*v(i,1:3) )
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, KE, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
GKE = KE/GNATOMS

!--- scale the obtained velocity to get the initial temperature.
vfactor = sqrt(1.5d0*treq/GKE)
v(:,:) = vfactor * v(:,:)

end subroutine

!-----------------------------------------------------------------------
subroutine adjust_temperature(atype, v)
!-----------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)

integer :: i,ity
real(8) :: Ekinetic, ctmp

Ekinetic=0.d0

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic=Ekinetic+0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, Ekinetic, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
Ekinetic=Ekinetic/GNATOMS

! adjust if current temperature deviates from treq by 5%
ctmp = sqrt( (treq*UTEMP0)/(Ekinetic*UTEMP) )
if( abs(ctmp-1.d0) > 0.05d0) then

   do i=1, NATOMS
      ity=nint(atype(i))
      v(i,1:3)=ctmp*v(i,1:3)
   enddo

   call linear_momentum(atype, v)

endif

return
end

!-----------------------------------------------------------------------
subroutine scale_temperature(atype, v)
!-----------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)

integer :: i,ity
integer,parameter :: MAX_ELEMENT=20
real(8) :: Ekinetic(2,MAX_ELEMENT), ctmp(MAX_ELEMENT)

Ekinetic(:,:)=0.d0

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic(1,ity)=Ekinetic(1,ity)+1.d0 
   Ekinetic(2,ity)=Ekinetic(2,ity)+0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, Ekinetic, size(Ekinetic), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

do ity=1, MAX_ELEMENT
   if(Ekinetic(1,ity)>1.d0) then
      ctmp(ity) = Ekinetic(2,ity)/Ekinetic(1,ity)
      ctmp(ity) = sqrt( (treq*UTEMP0)/(ctmp(ity)*UTEMP) )
      if(myid==0) print'(a,i4,f8.3,i9,es15.5)', &
          'atom_type, scaling_factor, num_samples, previous_Ekin: ', &
          ity, ctmp(ity), nint(Ekinetic(1,ity)), Ekinetic(2,ity)
   else
      ctmp(ity) = 0.d0
   endif
enddo

do i=1, NATOMS
   ity=nint(atype(i))
   v(i,1:3)=ctmp(ity)*v(i,1:3)
enddo

call linear_momentum(atype, v)

return
end

!-----------------------------------------------------------------------
subroutine linear_momentum(atype, v)
!-----------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)

integer :: i,ity
real(8) :: vcm(0:3)

!--- get the local momentum and mass.
vcm(:)=0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   vcm(1:3)=vcm(1:3) + mass(ity)*v(i,1:3)
   vcm(0)= vcm(0) + mass(ity)
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, vcm, size(vcm), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

!--- get the global momentum
vcm(1:3)=vcm(1:3)/vcm(0)

!--- set the total momentum to be zero 
do i=1, NATOMS
   v(i,1:3) = v(i,1:3) - vcm(1:3)
enddo

return
end

!----------------------------------------------------------------------
subroutine angular_momentum(atype, pos, v)
!----------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), pos(:,:)
real(8),allocatable,intent(in out):: v(:,:)

integer :: i,ity
real(8) :: com(3), Gcom(3), intsr(3,3), Gintsr(3,3), intsr_i(3,3), angm(3), Gangm(3), angv(3), mm, Gmm
real(8) :: dr(3), dv(3)

!--- get center of mass
com(:)=0.d0;     Gcom(:)=0.d0
mm=0.d0; Gmm=0.d0

do i=1, NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
   com(1:3) = mass(ity)*pos(i,1:3)
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, mm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gmm = mm
call MPI_ALLREDUCE(MPI_IN_PLACE, com, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gcom = com
Gcom(1:3) = Gcom(1:3)/Gmm

!--- get the angular momentum and inertia tensor from the com
angm(:)=0.d0;    Gangm(:)=0.d0
intsr(:,:)=0.d0; Gintsr(:,:)=0.d0

do i=1, NATOMS
   dr(1:3) = pos(i,1:3) - Gcom(1:3)
   
   angm(1) = mass(ity)*( dr(2)*v(i,3)-dr(3)*v(i,2) )
   angm(2) = mass(ity)*( dr(3)*v(i,1)-dr(1)*v(i,3) )
   angm(3) = mass(ity)*( dr(1)*v(i,2)-dr(2)*v(i,1) )

   intsr(1,1) = mass(ity)*( dr(2)**2+dr(3)**2 )
   intsr(2,2) = mass(ity)*( dr(3)**2+dr(1)**2 )
   intsr(3,3) = mass(ity)*( dr(1)**2+dr(2)**2 )

   intsr(1,2) =-mass(ity)*( dr(1)*dr(2) )
   intsr(1,3) =-mass(ity)*( dr(1)*dr(3) )
   intsr(2,3) =-mass(ity)*( dr(2)*dr(3) )

   intsr(2,1) = intsr(1,2)
   intsr(3,1) = intsr(1,3)
   intsr(3,2) = intsr(2,3)
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, angm, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gangm = angm
call MPI_ALLREDUCE(MPI_IN_PLACE, intsr, 9, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gintsr = intsr

!--- get angular velocity
call matinv(Gintsr, intsr_i)

angv(1) = sum(intsr_i(1,1:3)*angm(1:3))
angv(2) = sum(intsr_i(2,1:3)*angm(1:3))
angv(3) = sum(intsr_i(3,1:3)*angm(1:3))


!--- correct rotational motion wrt CoM.
do i=1,NATOMS
   dr(1:3) = pos(i,1:3) - Gcom(1:3)
   dv(1) = angv(2)*dr(3) - angv(3)*dr(2)
   dv(2) = angv(3)*dr(1) - angv(1)*dr(3)
   dv(3) = angv(1)*dr(2) - angv(2)*dr(1)

   v(i,1:3) = v(i,1:3) - dv(1:3)
enddo

end subroutine

!-----------------------------------------------------------------------
subroutine maximally_preserving_BD(atype, v, scaling_factor)
!-----------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:)
real(8),intent(in) :: scaling_factor

integer :: i,ity
real(8) :: kin 
integer :: Nsamples(size(atmname))
real(8) :: Ekinetic(0:size(atmname)), ctmp(0:size(atmname))

!--- get element-wise temperature first. 
Ekinetic = 0.d0; Nsamples=0
do i=1, NATOMS
   ity = nint(atype(i))
   kin = hmas(ity)*sum(v(i,1:3)*v(i,1:3))
   Ekinetic(ity) = Ekinetic(ity) + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
   Nsamples(ity) = Nsamples(ity) + 1
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, Ekinetic, size(Ekinetic), &
                   MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, Nsamples, size(Nsamples), &
                   MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

!--- get element-wise scaling factor. 
do ity = 1, size(atmname)
  if(Nsamples(ity)>0) then
     Ekinetic(ity) = Ekinetic(ity)/Nsamples(ity)
     ctmp(ity) = sqrt( (treq*UTEMP0)/(Ekinetic(ity)*UTEMP) )
  else
     ctmp(ity) = -1.d0
  endif
enddo

!--- adjust each atomic velocity & get system temperature after the adjustment.
do i=1, NATOMS
   ity = nint(atype(i))
   Ekinetic(0) = hmas(ity)*sum(v(i,1:3)*v(i,1:3))
   ctmp(0) = sqrt( (Ekinetic(0)*UTEMP)/(treq*UTEMP0) )

   if(ctmp(0)>scaling_factor) &
      v(i,1:3) = v(i,1:3)*sqrt(0.5d0*(scaling_factor+1.d0))/ctmp(0)

enddo

call linear_momentum(atype, v)

Ekinetic = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   Ekinetic(ity) = Ekinetic(ity) + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, Ekinetic, size(Ekinetic), &
                   MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

!--- get system temperature and its difference from treq
Ekinetic(0)=sum(Ekinetic(1:size(atmname)))/GNATOMS
ctmp(0) = sqrt((Ekinetic(0)*UTEMP)/(treq*UTEMP0))
v(1:natoms,1:3) = v(1:natoms,1:3)/ctmp(0)

if(myid==0) then
  write(6,fmt='(a,2f10.3,3x)',advance='no') &
     'T,ctmp : ', Ekinetic(0)*UTEMP, ctmp(0)

  do ity=1, size(atmname)
     if(Nsamples(ity)>0) &
        write(6,fmt='(i1,a2,2f10.3,2x)',advance='no') & 
        ity,'-',Ekinetic(ity)*UTEMP/Nsamples(ity),ctmp(ity)
  enddo
  write(6,fmt=*)
endif

return
end

!------------------------------------------------------------------------------
subroutine vkick(dtf, atype, v, f)
!------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:)
real(8),allocatable,intent(in out):: v(:,:), f(:,:)
real(8),intent(in) :: dtf

integer :: i, ity

do i=1,NATOMS
   ity = nint(atype(i))
   v(i,1:3) = v(i,1:3) + dtf*dthm(ity)*f(i,1:3)
enddo

end subroutine

end module
