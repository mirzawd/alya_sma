!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_openfi
  ! NAME 
  !    nsi_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names and open them to be used by 
  !    the module in two possible ways:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !       encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as
  !       argument when the binary file Alya is launched "naked".
  ! USES
  !    iofile
  ! USED BY
  !    nsi_turnon
  !***
  !-----------------------------------------------------------------------
  use def_nastin
  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask  
  character(150)          :: fil_psmat_nsi
  character(150)          :: fil_stasg_nsi
  character(150)          :: fil_cvgsg_nsi
  character(150)          :: fil_refer_nsi
  character(150)          :: fil_recvg_nsi
  character(150)          :: fil_lmach_nsi
  character(150)          :: fil_dynlo_nsi
  character(150)          :: fil_dynre_nsi
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then
     !
     ! Define unit opening option if this is a restart run
     !
     if(kfl_rstar==2) then 
        statu='old'
        forma='formatted'
        posit='append'
     else
        statu='unknown'
        forma='formatted'
        posit='asis'
     end if

     select case (itask)
 
     case (2_ip)
        !
        ! Open files needed occasionally
        !
        !ALYA-ADAN
        call nsi_cadan(1_ip)
        !ALYA-ADAN
        if(kfl_naked==0) then 
           call GET_ENVIRONMENT_VARIABLE('FOR112',fil_stasg_nsi)
           call GET_ENVIRONMENT_VARIABLE('FOR113',fil_cvgsg_nsi)
           call GET_ENVIRONMENT_VARIABLE('FOR114',fil_psmat_nsi)
           call GET_ENVIRONMENT_VARIABLE('FOR117',fil_refer_nsi)     
           call GET_ENVIRONMENT_VARIABLE('FOR118',fil_recvg_nsi)     
           call GET_ENVIRONMENT_VARIABLE('FOR120',fil_lmach_nsi)     
           call GET_ENVIRONMENT_VARIABLE('FOR122',fil_dynin_nsi)      
           call GET_ENVIRONMENT_VARIABLE('FOR123',fil_dynou_nsi)      
           call GET_ENVIRONMENT_VARIABLE('FOR124',fil_dynlo_nsi)      
           call GET_ENVIRONMENT_VARIABLE('FOR125',fil_dynre_nsi)     
        else
           fil_stasg_nsi = adjustl(trim(namda))//'-sgs.'                //exmod(modul)//'.log'
           fil_cvgsg_nsi = adjustl(trim(namda))//'-sgs.'                //exmod(modul)//'.cvg'
           fil_psmat_nsi = adjustl(trim(namda))//'-matrix.'             //exmod(modul)//'.ps'
           fil_refer_nsi = adjustl(trim(namda))//'-reference.'          //exmod(modul)//'.res'
           fil_recvg_nsi = adjustl(trim(namda))//'-reference.'          //exmod(modul)//'.cvg'   
           fil_lmach_nsi = adjustl(trim(namda))//'-low-mach.'           //exmod(modul)//'.cvg'    
           fil_dynin_nsi = adjustl(trim(namda))//'-dynamic-coupling.'   //exmod(modul)//'.in'     
           fil_dynou_nsi = adjustl(trim(namda))//'-dynamic-coupling.'   //exmod(modul)//'.out'    
           fil_dynlo_nsi = adjustl(trim(namda))//'-dynamic-coupling.'   //exmod(modul)//'.log'    
           fil_dynre_nsi = adjustl(trim(namda))//'-dynamic-coupling.'   //exmod(modul)//'.res' 
        end if
        ! 
        ! Subgrid scales
        !
        if(kfl_sgsco_nsi==1) then
           call iofile(zero,lun_stasg_nsi,fil_stasg_nsi,'NASTIN SUBGRID SCALE STAT.',statu,forma,posit)
           call iofile(zero,lun_cvgsg_nsi,fil_cvgsg_nsi,'NASTIN SUBGRID SCALE CONV.',statu,forma,posit)
        end if
        !
        ! Matrix profile
        !
        if(kfl_psmat_nsi>0) &
             call iofile(zero,lun_psmat_nsi,fil_psmat_nsi,'NASTIN MATRIX',statu,forma,posit)
        !
        ! Reference solution
        !
        if(kfl_refer_nsi==1) then
           call iofile(zero,lun_refer_nsi,fil_refer_nsi,'NASTIN REFERENCE SOLUTION',statu,forma,posit)
           call iofile(zero,lun_recvg_nsi,fil_recvg_nsi,'NASTIN REFERENCE CONVERGENCE',statu,forma,posit)
        end if
        !
        ! Low-Mach model
        !
        if(kfl_regim_nsi==3) &
             call iofile(zero,lun_lmach_nsi,fil_lmach_nsi,'NASTIN LOW-MACH',statu,forma,posit) 

     case (4_ip)
        !
        ! Close output file
        !
        if(kfl_refer_nsi/=0) call iofile(two,lun_recvg_nsi,' ','NASTIN REFERENCE CONVERGENCE')

     case (5_ip)

     case (6_ip)

     case (7_ip)

     case (8_ip)
        !
        ! Close matrix file
        !
        call iofile(two,lun_psmat_nsi,' ','NASTIN MATRIX PROFILE')

     case (9_ip)
        !
        ! Close reference solution file
        !
        if(kfl_refer_nsi/=0) &
             call iofile(two,lun_refer_nsi,' ','NASTIN REFERENCE SOLUTION')

     end select

  end if

end subroutine nsi_openfi

