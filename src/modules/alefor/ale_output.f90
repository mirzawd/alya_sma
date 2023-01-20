!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_output()
  !------------------------------------------------------------------------
  !****f* Alefor/ale_output
  ! NAME 
  !    ale_output
  ! DESCRIPTION
  !    Output and postprocess of solution
  ! USES
  ! USED BY
  !    ale_endite
  !    ale_endste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use mod_iofile
  use mod_outfor,              only : outfor
  use mod_output_postprocess,  only : output_postprocess_variables
  implicit none
  external             :: ale_outvar
  integer(ip)          :: iimbo,iiset,inaux
  integer(ip), save    :: ipass = 0
  real(rp),    pointer :: force(:,:),accel(:,:),velol(:,:),posil(:,:)   ! Linear motion
  real(rp),    pointer :: torqu(:,:),accea(:,:),veloa(:,:),posia(:,:)   ! Angular motion
  real(rp),    pointer :: pp_pf(:,:),pp_vf(:,:),pp_pt(:,:),pp_vt(:,:)

  real(rp),    pointer :: Mass,Rho,Vol,Momi(:),COG(:)
  !
  ! Initial solution, end of a time step and and of run
  !
  call output_postprocess_variables(ale_outvar)

  !---------------------------------------------------------------
  !
  ! Rigid Body properties
  !
  !---------------------------------------------------------------

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then

     if( INOTSLAVE ) then
        !
        ! Output results
        !
        if (kfl_crist_ale == 1_ip) then 
           inaux = 3_ip
        else
           inaux = 1_ip
        end if
        if ( kfl_rigid_ale == 1 ) then
           if ( ( ipass == 0 ) .and. (kfl_rstar /= 2) ) write(lun_outpu_ale,99) nrbod
           ipass    = ipass + 1
           ioutp(1) = ittim
           routp(1) = cutim
           call outfor(46_ip,lun_outpu_ale,' ')
           if( ipass == 1 ) then
              do iimbo = 1,nrbod
                 Mass  => rbbou(iimbo) % massa
                 Rho   => rbbou(iimbo) % densi
                 Vol   => rbbou(iimbo) % volum
                 Momi  => rbbou(iimbo) % momin
                 COG   => rbbou(iimbo) % posgr
                 write(momod(modul) % lun_outpu,102) iimbo,Vol,Rho,Mass
                 write(momod(modul) % lun_outpu,103) iimbo,COG
                 write(momod(modul) % lun_outpu,104) iimbo,Momi
              end do
              flush(momod(modul) % lun_outpu)
           end if

           do iimbo = 1,nrbod
              force  => rbbou(iimbo) % force
              accel  => rbbou(iimbo) % accel
              velol  => rbbou(iimbo) % velol
              posil  => rbbou(iimbo) % posil
              torqu  => rbbou(iimbo) % torqu
              accea  => rbbou(iimbo) % accea
              veloa  => rbbou(iimbo) % veloa
              posia  => rbbou(iimbo) % posia

              pp_pf  => rbbou(iimbo) % pp_pf   ! to postprocess forces separatelly by sets
              pp_vf  => rbbou(iimbo) % pp_vf
              pp_pt  => rbbou(iimbo) % pp_pt
              pp_vt  => rbbou(iimbo) % pp_vt

              do iiset = 1,rbbou(iimbo) % nrbse
                 write(lun_outpu_ale,101) &
                      rbbou(iimbo) % lrbse(iiset),                       &
                      posil(1,1)     , posil(2,1)     , posil(3,1)     , &
                      velol(1,1)     , velol(2,1)     , velol(3,1)     , &
                      accel(1,1)     , accel(2,1)     , accel(3,1)     , &
                      posia(1,1)     , posia(2,1)     , posia(3,1)     , &
                      veloa(1,1)     , veloa(2,1)     , veloa(3,1)     , &
                      accea(1,1)     , accea(2,1)     , accea(3,1)     , &
                      force(1,inaux) , force(2,inaux) , force(3,inaux) , &
                      pp_vf(1,iiset) , pp_vf(2,iiset) , pp_vf(3,iiset) , &
                      pp_pf(1,iiset) , pp_pf(2,iiset) , pp_pf(3,iiset) , &
                      torqu(1,inaux) , torqu(2,inaux) , torqu(3,inaux) , &
                      pp_vt(1,iiset) , pp_vt(2,iiset) , pp_vt(3,iiset) , &
                      pp_pt(1,iiset) , pp_pt(2,iiset) , pp_pt(3,iiset)
              end do
           end do
           flush(lun_outpu_ale)
        end if
     end if

  end if
  !
  ! Formats
  !
99 format('# ALYA Boundary set results',/,&
       & '#',/,&
       & '# HEADER',/,&
       & '# ISET, Column :   1',/,&
       & '# X,    Column :   2',/,&
       & '# Y,    Column :   3',/,&
       & '# Z,    Column :   4',/,&
       & '# VX,   Column :   5',/,&
       & '# VY,   Column :   6',/,&
       & '# VZ,   Column :   7',/,&
       & '# AX,   Column :   8',/,&
       & '# AY,   Column :   9',/,&
       & '# AZ,   Column :  10',/,&
       & '# S1,   Column :  11',/,&
       & '# S2,   Column :  12',/,&
       & '# S3,   Column :  13',/,&
       & '# W1,   Column :  14',/,&
       & '# W2,   Column :  15',/,&
       & '# W3,   Column :  16',/,&
       & '# Z1,   Column :  17',/,&
       & '# Z2,   Column :  18',/,&
       & '# Z3,   Column :  19',/,&
       & '# F1,   Column :  20',/,&
       & '# F2,   Column :  21',/,&
       & '# F3,   Column :  22',/,&
       & '# Fv1,  Column :  23',/,&
       & '# Fv2,  Column :  24',/,&
       & '# Fv3,  Column :  25',/,&
       & '# Fp1,  Column :  26',/,&
       & '# Fp2,  Column :  27',/,&
       & '# Fp3,  Column :  28',/,&
       & '# T1,   Column :  29',/,&
       & '# T2,   Column :  30',/,&
       & '# T3,   Column :  31',/,&
       & '# Tv1,  Column :  32',/,&
       & '# Tv2,  Column :  33',/,&
       & '# Tv3,  Column :  34',/,&
       & '# Tp1,  Column :  35',/,&
       & '# Tp2,  Column :  36',/,&
       & '# Tp3,  Column :  37',/,&
       & '# NUMVARIABLES :  37',/,&
       & '# NUMSETS      :  ',i4,/,&
       & '# START')

101 format(4x,i9,50(2x,e12.6))

102 format('# ISET, volum,densi,massa: ',i6,3(1x,e12.5))
103 format('# ISET, posgr: ',i6,3(1x,e12.5))
104 format('# ISET, : MOM. INERC. ',i6,6(1x,e12.5))

end subroutine ale_output
