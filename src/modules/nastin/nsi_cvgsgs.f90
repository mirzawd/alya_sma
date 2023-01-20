!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_cvgsgs()
  !-----------------------------------------------------------------------
  !****f* nastin/nsi_cvgsgs
  ! NAME 
  !    nsi_cvgsgs
  ! DESCRIPTION
  !    This routine outputs SGS statistics
  ! USES
  ! USED BY
  !    nsi_endite (itask=1,2)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_iofile,         only : iofile_flush_unit
  implicit none
  integer(ip), save :: kpass=0
  integer(ip)       :: kk,ii,mcomp
  !
  ! RESGS_NSI: Compute SGS residual
  !
  call PAR_SUM(2_ip,resgs_nsi,'IN MY CODE')
  if( resgs_nsi(2) > 0.0_rp ) resgs_nsi(1) = sqrt(resgs_nsi(1)/resgs_nsi(2))
  !
  ! Maximum subgrid scale residual
  !
  if( kfl_sgsco_nsi /= 0 .and. IPARALL ) then
     call PAR_MAX(rmsgs_nsi,'IN MY CODE')
  end if
  !
  ! Subgrid scale statistics for convergence 
  !
  if( kfl_sgsco_nsi /= 0 ) then

     if( kpass == 0 .and. kfl_rstar /= 2 .and. INOTSLAVE ) then              
        write(lun_stasg_nsi,200) tosgs_nsi,misgs_nsi
        write(lun_cvgsg_nsi,300)
     end if

     call PAR_SUM(mmsgs_nsi,itsta_nsi,'IN MY CODE')

     ii = 0
     do kk = 1,mmsgs_nsi
        ii = ii + itsta_nsi(kk)
     end do
     if( ii == 0 ) ii = 1

     if( INOTSLAVE )&
          write(lun_stasg_nsi,201) &
          (100.0_rp*real(itsta_nsi(kk),rp)/real(ii,rp),kk=1,mmsgs_nsi)

     domcomp: do mcomp = mmsgs_nsi,1,-1
        if( itsta_nsi(mcomp) /= 0 ) exit domcomp
     end do domcomp

     call PAR_SUM(2_ip,mcomp,resis_nsi,'IN MY CODE')

     if( INOTSLAVE ) then
        do kk = 1,mcomp
           if( resis_nsi(2,kk) < 1.0e-10_rp ) resis_nsi(2,kk)=1.0_rp
           resis_nsi(1,kk) = sqrt(resis_nsi(1,kk)/resis_nsi(2,kk))
           write(lun_cvgsg_nsi,301) &
                ittim,itcou,itinn(modul),kk,cutim,resis_nsi(1,kk)
        end do
        call iofile_flush_unit(lun_stasg_nsi)
        call iofile_flush_unit(lun_cvgsg_nsi)
     end if

  end if
  kpass = 1 

  !----------------------------------------------------------------------
  !
  ! Formats
  !
  !----------------------------------------------------------------------

200 format(&
       & '# Percentage of the number of iterations to achieve the L2 tolerance= ',e12.6,&
       & ' with a maxiumum number of iterations= ',i3,/,&
       &      '#            1','             2','             3','             4',&
       &      '             5','             6','             7','             8',&
       &      '             9','            10','            11','            12',&
       &      '            13','            14','            15','           >16')
201 format(500(2x,e12.6))
300 format('# --| ALYA convergence '  ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time Step         2. Global Iteration  3. Inner Iteration   ',/,&
       & '# --| 4. Inner SGS         5. Current time      6. Velocity SGS      ',/,&
       & '# ','          1','          2','          3','          4',&
       &      '             5','             6') 
301 format(2x,4(2x,i9),3(2x,e12.6))

end subroutine nsi_cvgsgs
