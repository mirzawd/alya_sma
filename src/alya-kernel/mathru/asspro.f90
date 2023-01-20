!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine asspro(&
     itask,pnode,ndofn,pgaus,lnods,lelch,gpden,gpvis,gpvol,&
     gpsha,elrhs,prope)
  !-----------------------------------------------------------------------
  !****f* mathru/asspro
  ! NAME 
  !    asspro
  ! DESCRIPTION
  !    Assembly of properties
  ! USES
  ! USED BY
  !    *_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only        :  ip,rp
  use def_elmtyp, only        :  ELEXT
  implicit none
  integer(ip),  intent(in)    :: itask,pnode,ndofn,pgaus
  integer(ip),  intent(in)    :: lnods(pnode)
  integer(ip),  intent(in)    :: lelch
  real(rp),     intent(in)    :: gpden(pgaus)
  real(rp),     intent(in)    :: gpvis(pgaus)
  real(rp),     intent(in)    :: gpvol(pgaus)
  real(rp),     intent(in)    :: gpsha(pnode,pgaus)
  real(rp),     intent(out)   :: elrhs(ndofn,*)
  real(rp),     intent(inout) :: prope(ndofn,*)
  integer(ip)                 :: inode,ipoin,idofn,igaus
  real(rp)                    :: fact1,fact2

  do inode = 1,pnode
     do idofn = 1,ndofn
        elrhs(idofn,inode) = 0.0_rp
     end do
  end do
   
  if( itask == 10 ) then
     !
     ! Density
     !
     do igaus = 1,pgaus
        fact1 = gpden(igaus) * gpvol(igaus)
        do inode = 1,pnode
           elrhs(1,inode) = elrhs(1,inode) + fact1 * gpsha(inode,igaus) 
        end do
     end do

  else if( itask == 11 ) then
     !
     ! Viscosity
     !
     do igaus = 1,pgaus
        fact1 = gpvis(igaus) * gpvol(igaus)
        do inode = 1,pnode
           elrhs(1,inode) = elrhs(1,inode) + fact1 * gpsha(inode,igaus) 
        end do
     end do
     
  else if( itask == 12 ) then
     !
     ! Density and viscosity
     ! 
     do igaus = 1,pgaus
        fact1 = gpden(igaus) * gpvol(igaus)
        fact2 = gpvis(igaus) * gpvol(igaus)
        do inode = 1,pnode
           elrhs(1,inode) = elrhs(1,inode) + fact1 * gpsha(inode,igaus) 
           elrhs(2,inode) = elrhs(2,inode) + fact2 * gpsha(inode,igaus) 
        end do
     end do
     
  end if

  !----------------------------------------------------------------------
  !
  ! Assembly
  !
  !----------------------------------------------------------------------

  if( lelch == ELEXT ) then

     inode = 1

     if( itask == 10 .or. itask == 11 ) then
        ipoin          = lnods(inode)
        !$OMP ATOMIC
        prope(1,ipoin) = prope(1,ipoin) + elrhs(1,inode)
        
     else if( itask == 12 ) then
        
        ipoin          = lnods(inode)
        !$OMP ATOMIC
        prope(1,ipoin) = prope(1,ipoin) + elrhs(1,inode)
        !$OMP ATOMIC
        prope(2,ipoin) = prope(2,ipoin) + elrhs(2,inode)

     end if

  else

     if( itask == 10 .or. itask == 11 ) then

        do inode = 1,pnode
           ipoin          = lnods(inode)
           !$OMP ATOMIC
           prope(1,ipoin) = prope(1,ipoin) + elrhs(1,inode)
        end do
        
     else if( itask == 12 ) then

        do inode = 1,pnode
           ipoin          = lnods(inode)
           !$OMP ATOMIC
           prope(1,ipoin) = prope(1,ipoin) + elrhs(1,inode)
           !$OMP ATOMIC
           prope(2,ipoin) = prope(2,ipoin) + elrhs(2,inode)
        end do
        
     end if
  end if

end subroutine asspro
