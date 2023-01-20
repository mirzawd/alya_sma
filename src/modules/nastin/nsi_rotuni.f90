!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_rotuni(itask,unrot_nsi,ipoin)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_rotuni
  ! NAME 
  !    nsi_rotuni
  ! DESCRIPTION
  !    This routine rotates the nodal velocities using the appropiate 
  !    rotation matrix:
  !    Contrary to nsi_rotunk that rotates all ipoin, this only rotates 1 ipoin.
  !    further it only works with ndofn=ndime therfore it has only 2 tasks 1 & 2
  !    ITASK=1 ... From global to local
  !    ITASK=2 ... From local to global
  !    Modifications need to be done only if there exist image nodes or
  !    boundary conditions in skew systems.
  ! USES
  !    mbvab0
  ! USED BY
  !    nsi_bouope
  !***
  !-----------------------------------------------------------------------
  use def_domain
  use def_master
  use def_nastin
  implicit none
  integer(ip), intent(in)    :: itask,ipoin
  real(rp),    intent(inout) :: unrot_nsi(ndime) 
  integer(ip)                :: ibopo,iroty,ndofn
  real(rp)                   :: worma(3),worve(3)

  if( kfl_local_nsi /= 0 .and. INOTMASTER ) then 

     ndofn = ndime

     if( itask == 1 ) then
        !
        ! Global to local
        !
        ibopo=lpoty(ipoin)
        if(ibopo>0) then
           iroty=kfl_fixrs_nsi(ipoin)

           if( iroty == -1 .or. iroty == -4 ) then
              !
              ! Boundary conditions in the tangent skew system
              !           
              worve(1:ndime) = unrot_nsi
              call mbvatb(worma,exnor(1,1,ibopo),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           else if( iroty == -2 ) then
              !
              ! Boundary conditions in the NSI tangent skew system
              !           
              worve(1:ndime) = unrot_nsi
              call mbvatb(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           else if( iroty == -3 ) then
              !
              ! Boundary conditions in geometrical system
              !           
              worve(1:ndime) = unrot_nsi 
              call mbvatb(worma,skcos(1,1,ibopo),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           else if( iroty >= 1 ) then
              !
              ! Boundary conditions in a given skew system
              ! 
              worve(1:ndime) = unrot_nsi
              call mbvatb(worma,skcos(1,1,iroty),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           end if
        end if

     else if( itask == 2 ) then
        !
        ! Local to global
        !
        ibopo=lpoty(ipoin)
        if(ibopo>0) then
           iroty=kfl_fixrs_nsi(ipoin)
           if(iroty==-1 .or. iroty == -4 ) then
              !
              ! Boundary conditions in the tangent skew system
              !           
              worve(1:ndime) = unrot_nsi 
              call mbvab0(worma,exnor(1,1,ibopo),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           else if(iroty==-2) then
              !
              ! Boundary conditions in the NSI skew system
              !           
              worve(1:ndime) = unrot_nsi 
              call mbvab0(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           else if(iroty==-3) then
              !
              ! Boundary conditions in geometrical system
              !           
              worve(1:ndime) = unrot_nsi 
              call mbvab0(worma,skcos(1,1,ibopo),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)

           else if(iroty>=1) then
              !
              ! Boundary conditions in a given skew system
              ! 
              worve(1:ndime) = unrot_nsi 
              call mbvab0(worma,skcos(1,1,iroty),worve,ndime,ndime)
              unrot_nsi = worma(1:ndime)
           end if
        end if

     end if

  end if

end subroutine nsi_rotuni

