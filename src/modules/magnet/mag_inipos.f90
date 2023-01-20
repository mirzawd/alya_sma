!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inipos()

  use def_master, only: postp

  implicit none
  !
  ! Postprocess Variable
  ! name:     '.....'
  ! type:     'SCALA/VECTO'
  ! entity:   'NPOIN/NELEM'
  !
  ! Magnetic Intensity H @ nodes
  !
  postp(1) % wopos(1:3, 1) = [ 'MAGNE', 'VECTO', 'NPOIN' ]
  !
  ! Current Density z Jz @ elements (ONLY 2D)
  !
  postp(1) % wopos(1:3, 2) = [ 'CURCZ', 'SCALA', 'NELEM' ]
  !
  ! Current Density z Jz @ nodes (ONLY 2D)
  !
  postp(1) % wopos(1:3, 3) = [ 'CURNZ', 'SCALA', 'NPOIN' ]
  !
  ! Magnetic Intensity H @ elements
  !
  postp(1) % wopos(1:3, 4) = [ 'MAGCE', 'VECTO', 'NELEM' ]
  !
  ! Magnetic Field Density B @ nodes
  !
  postp(1) % wopos(1:3, 5) = [ 'FLUXN', 'VECTO', 'NPOIN' ]
  !
  ! Magnetic Field Density B @ elements
  !
  postp(1) % wopos(1:3, 6) = [ 'FLUXC', 'VECTO', 'NELEM' ]
  !
  ! Current Density J @ elements (ONLY 3D)
  !
  postp(1) % wopos(1:3, 7) = [ 'CURRC', 'VECTO', 'NELEM' ]
  ! 
  ! Current Density J @ nodes (ONLY 3D)
  !
  postp(1) % wopos(1:3, 8) = [ 'CURRN', 'VECTO', 'NPOIN' ]
  !
  ! EM Force F @ nodes
  !
  postp(1) % wopos(1:3, 9) = [ 'FORCN', 'VECTO', 'NPOIN' ]
  !
  ! EM Force F @ elements
  !
  postp(1) % wopos(1:3, 10) = [ 'FORCE', 'VECTO', 'NELEM' ]
  !
  ! Heat Dissipation G @ nodes
  !
  postp(1) % wopos(1:3, 11) = [ 'JOULN', 'SCALA', 'NPOIN' ]
  !
  ! Heat Dissipation G @ elements
  !
  postp(1) % wopos(1:3, 12) = [ 'JOULC', 'SCALA', 'NELEM' ]
  !
  ! Element set variables
  !
  postp(1) % woese(1) = 'VARI1'       ! Variable 1
  postp(1) % woese(2) = 'VARI2'       ! Variable 2
  postp(1) % woese(3) = 'VARI3'       ! Variable 3
  !
  ! Boundary set variables
  !
  postp(1) % wobse(1) = 'VARI1'       ! Variable 1
  postp(1) % wobse(2) = 'VARI2'       ! Variable 2
  postp(1) % wobse(3) = 'VARI3'       ! Variable 3
  !
  ! Node set variables
  !
  postp(1) % wonse(1) = 'VARI1'       ! Variable 1
  postp(1) % wonse(2) = 'VARI2'       ! Variable 2
  postp(1) % wonse(3) = 'VARI3'       ! Variable 3
  !
  ! Witness variables
  !
  postp(1) % wowit(1) = 'VARI1'       ! Variable 1
  postp(1) % wowit(2) = 'VARI2'       ! Variable 2
  postp(1) % wowit(3) = 'VARI3'       ! Variable 3

end subroutine mag_inipos
