!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_rotunk(itask,unrot_nsi)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_rotunk
  ! NAME 
  !    nsi_rotunk
  ! DESCRIPTION
  !    This routine rotates the nodal velocities using the appropiate 
  !    rotation matrix:
  !    ITASK=1 ... From global to local
  !    ITASK=2 ... From local to global
  !    Modifications need to be done only if there exist image nodes or
  !    boundary conditions in skew systems.
  !                          Monolithic  Schur
  !      global to local :     1          3
  !      local to global:      2          4
  !
  ! USES
  !    mbvab0
  ! USED BY
  !    nsi_solite
  !***
  !-----------------------------------------------------------------------
  use def_domain
  use def_master
  use def_nastin
  use mod_local_basis, only : local_basis_global_to_local
  use mod_local_basis, only : local_basis_local_to_global
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(inout) :: unrot_nsi(*) 
  integer(ip)                :: ndofn
!  integer(ip)                :: ipoin,ibopo,iroty,itotv
!  real(rp)                   :: worma(3),worve(3)

  if( kfl_local_nsi /= 0 .and. INOTEMPTY ) then 
 
     if( itask == 1 .or. itask == 2 ) then
        ndofn = solve(1) % ndofn
     else
        ndofn = ndime

     end if

     if( itask == 1 .or. itask == 3 ) then
        !
        ! Global to local
        !
        call local_basis_global_to_local(ndime,kfl_fixrs_nsi,unrot_nsi)
        
!!$        do ipoin=1,npoin
!!$           ibopo=lpoty(ipoin)
!!$           if(ibopo>0) then
!!$              iroty=kfl_fixrs_nsi(ipoin)
!!$
!!$              if( iroty == -1 .or. iroty == -4 ) then
!!$                 !
!!$                 ! Boundary conditions in the tangent skew system
!!$                 !           
!!$                 itotv = (ipoin-1)*ndofn
!!$                 worve(1:ndime) = unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvatb(worma,exnor(1,1,ibopo),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime) = worma(1:ndime)
!!$
!!$              else if( iroty == -2 ) then
!!$                 !
!!$                 ! Boundary conditions in the NSI tangent skew system
!!$                 !           
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvatb(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$                 
!!$              else if( iroty == -3 ) then
!!$                 !
!!$                 ! Boundary conditions in geometrical system
!!$                 !           
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvatb(worma,skcos(1,1,ibopo),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$                 
!!$              else if( iroty >= 1 ) then
!!$                 !
!!$                 ! Boundary conditions in a skew system given by a field
!!$                 ! 
!!$                 !itotv=(ipoin-1)*ndofn
!!$                 !worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)                 
!!$                 !call mbvatb(worma,xfiel(iroty) % a(1,ipoin,1),worve,ndime,ndime)
!!$                 !unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$
!!$              end if
!!$           end if
!!$        end do

     else if( itask == 2 .or. itask == 4 ) then
        !
        ! Local to global
        !
        call local_basis_local_to_global(ndime,kfl_fixrs_nsi,unrot_nsi)
        
!!$        do ipoin=1,npoin
!!$           ibopo=lpoty(ipoin)
!!$           if(ibopo>0) then
!!$              iroty=kfl_fixrs_nsi(ipoin)
!!$              if(iroty==-1 .or. iroty == -4 ) then
!!$                 !
!!$                 ! Boundary conditions in the tangent skew system
!!$                 !           
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvab0(worma,exnor(1,1,ibopo),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$
!!$              else if(iroty==-2) then
!!$                 !
!!$                 ! Boundary conditions in the NSI skew system
!!$                 !           
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvab0(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$
!!$              else if(iroty==-3) then
!!$                 !
!!$                 ! Boundary conditions in geometrical system
!!$                 !           
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvab0(worma,skcos(1,1,ibopo),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$
!!$              else if(iroty>=1) then
!!$                 !
!!$                 ! Boundary conditions in a given skew system
!!$                 ! 
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvab0(worma,skcos(1,1,iroty),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$                 
!!$              else if(iroty>=1) then
!!$                 !
!!$                 ! Boundary conditions in a given skew system
!!$                 ! 
!!$                 itotv=(ipoin-1)*ndofn
!!$                 worve(1:ndime)=unrot_nsi(itotv+1:itotv+ndime)
!!$                 call mbvab0(worma,skcos(1,1,iroty),worve,ndime,ndime)
!!$                 unrot_nsi(itotv+1:itotv+ndime)=worma(1:ndime)
!!$              end if
!!$           end if
!!$        end do

     end if

  end if

end subroutine nsi_rotunk

