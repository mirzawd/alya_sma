!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_reaphy()
  !------------------------------------------------------------------------
  !****f* Wavequ/lev_reaphy
  ! NAME 
  !    lev_rwaphy
  ! DESCRIPTION
  !    Read physical problem
  ! USES
  ! USED BY
  !    lev_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_levels
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_advec_lev = 0_ip
     kfl_reave_lev = 0_ip
     nmate_lev     = 1
     thicl         = -1.0_rp
     !
     ! Reach the section
     !
     call ecoute('lev_reaphy')
     do while( words(1) /= 'INITI' .and. words(1) /= 'PHYSI' )
        call ecoute('lev_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while( words(1) /= 'ENDIN' .and. words(1) /= 'ENDPH' )
        call ecoute('lev_reaphy')
        !
        ! Interface thickness
        !
        if(exists('THICK')) thicl = getrea('THICK',0.0_rp,'#Interface thickness')
        !
        ! Function definition data
        !
        if(words(1)=='VELOC') then               ! Velocity field Initialisation
           call ecoute('lev_reaphy')
           do while(words(1)/='ENDVE')
              if(words(1)=='DEFVE') then               ! Source term
                 if(exists('FUNCT')) then                 
                    kfl_advec_lev = getint('FUNCT',4_ip,  '#Source term function')
                 end if
                 if(exists('INITI')) then                 
                    kfl_reave_lev = getint('INITI',1_ip,  '# Velocity initialisation from file')
                 end if
              end if
              call ecoute('lev_reaphy')
           end do

        else if(words(1)=='CONVE') then               ! Convective term
           if(exists('ON   ')) then
              kfl_advec_lev = 1
              if(exists('FUNCT')) &
                   kfl_advec_lev = getint('FUNCT',0_ip,'#Velocity function')
              if(words(3)=='VELOC') then
                 if(words(4)=='FUNCT') then
                    kfl_advec_lev = getint('FUNCT',0_ip,'#Velocity function')
                 else
                    kfl_advec_lev = 1                       
                 end if
              end if
           end if

        end if
     end do
   
     if ( (thicl == -1.0_rp) .and. (kfl_prope == 0) ) call runend ('LEV_REAPHY: thickness must be read')
     
  end if

end subroutine lev_reaphy
