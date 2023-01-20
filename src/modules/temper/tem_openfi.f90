!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_openfi(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_openfi
  ! NAME 
  !    tem_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names and open them to be used by 
  !    the module in two possible ways:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked".  
  ! USES
  ! USED BY
  !    tem_turnon
  !------------------------------------------------------------------------
  use def_temper
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(150)          :: fil_splot_tem,fil_psmat_tem
  character(150)          :: fil_bound_tem,fil_funck_tem
  character(150)          :: fil_funcc_tem,fil_intbc_tem
  character(150)          :: fil_dynlo_tem,fil_dynre_tem,fil_lmach_tem
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit
  character(20)           :: wmate

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

     if( INOTSLAVE ) then
        !
        ! Open files needed occasionally
        !
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR210',fil_bound_tem)   
           call GET_ENVIRONMENT_VARIABLE('FOR221',fil_intbc_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR222',fil_dynin_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR223',fil_dynou_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR224',fil_dynlo_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR225',fil_dynre_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR232',fil_splot_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR233',fil_ramsh_tem)      
           call GET_ENVIRONMENT_VARIABLE('FOR234',fil_rares_tem)       
           call GET_ENVIRONMENT_VARIABLE('FOR235',fil_lmach_tem)     
        else
           fil_bound_tem = adjustl(trim(namda))//'.'                 //exmod(modul)//'.bcs'
           fil_intbc_tem = adjustl(trim(namda))//'-bcinterpolation.' //exmod(modul)//'.fix'   
           fil_dynin_tem = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.in'     
           fil_dynou_tem = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.out'    
           fil_dynlo_tem = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.log'    
           fil_dynre_tem = adjustl(trim(namda))//'-dynamic-coupling.'//exmod(modul)//'.res'    
           fil_splot_tem = adjustl(trim(namda))//'-specificheat.'    //exmod(modul)//'.fun'     
           fil_ramsh_tem = adjustl(trim(namda))//'-radiation.'       //exmod(modul)//'.post.msh'     
           fil_rares_tem = adjustl(trim(namda))//'-radiation.'       //exmod(modul)//'.post.res'  
           fil_lmach_tem = adjustl(trim(namda))//'-low-mach.'        //exmod(modul)//'.cvg'
        end if
        !
        ! Radiation files
        !
        if(kfl_radia_tem==1.and.kfl_viewf_tem(1)/=0.and.kfl_rstar/=2) then
           call iofile(zero,lun_ramsh_tem,fil_ramsh_tem,'TEMPER RADIATION MESH')
           call iofile(zero,lun_rares_tem,fil_rares_tem,'TEMPER RADIATION RESULT')     
        end if
        !
        ! Surface plot file
        !
        if(kfl_splot_tem==1) &
             call iofile(zero,lun_splot_tem,fil_splot_tem,'TEMPER SURFACE PLOT ',statu,forma,posit)
        !
        ! Boundary conditions
        !
        if(npp_bound_tem>0) &
             call iofile(zero,lun_bound_tem,fil_bound_tem,'TEMPER BOUND. COND. ',statu,forma,posit)
        !
        ! Bc interpolation file
        !
        if(kfl_intbc_tem/=0) &
             call iofile(zero,lun_intbc_tem,fil_intbc_tem,'TEMPER BC INTERPOLATION','old')  
        !
        ! Dynamic coupling
        !
        if(kfl_dynco_tem/=0) then
           call iofile(zero,lun_dynlo_tem,fil_dynlo_tem,'TEMPER DYNAMIC COUPLING LOG')  
           call iofile(zero,lun_dynre_tem,fil_dynre_tem,'TEMPER DYNAMIC COUPLING RESULT')  
        end if
        !
        ! Low-Mach model
        !
        if(kfl_regim_tem>=3) then
           call iofile(zero,lun_lmach_tem,fil_lmach_tem,'TEMPER LOW-MACH',statu,forma,posit)
        end if
     end if
     !
     ! Matrix profile
     !
     if( INOTMASTER .and. kfl_psmat_tem /=0 ) then
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR214',fil_psmat_tem) 
        else
           if( kfl_psmat_tem < 0 ) then
              fil_psmat_tem = adjustl(trim(namda))//'-matrix.'//exmod(modul)//'.post.msh'                  
           else
              if( ISEQUEN ) then
                 fil_psmat_tem = adjustl(trim(namda))//'-matrix.'//exmod(modul)//'.ps'    
              else
                 fil_psmat_tem = adjustl(trim(namda))//'-matrix-'//trim(intost(kfl_paral))//'.'//exmod(modul)//'.ps'    
              end if
           end if
        end if
        call iofile(zero,lun_psmat_tem,fil_psmat_tem,'TEMPER MATRIX',statu,forma,posit)
     end if

  case(4_ip)
     !
     ! Close matrix file
     !
     if( INOTMASTER ) call iofile(two,lun_psmat_tem,' ','TEMPER MATRIX PROFILE')

  case(5_ip)
     !
     ! Interpolation file for k
     !
     if( INOTMASTER ) then
        if (kfl_naked==0) then 
           call GET_ENVIRONMENT_VARIABLE('FOR215',fil_funck_tem)       
        else
           fil_funck_tem = adjustl(trim(namda))//'-conductivity.'//exmod(modul)//'.fun'   
        end if
        wmate=intost(imate_tem)
        fil_funck_tem=trim(fil_funck_tem)//trim(wmate)
        call iofile(zero,lun_funck_tem,trim(fil_funck_tem),'TEMPER K INTERPOLATION','old')
     end if

  case(6_ip)
     !
     ! Interpolation file for Cp
     !
     if( INOTMASTER ) then 
        if (kfl_naked==0) then 
           call GET_ENVIRONMENT_VARIABLE('FOR216',fil_funcc_tem)       
        else
           fil_funcc_tem = adjustl(trim(namda))//'-specificheat.'//exmod(modul)//'.fun'   
        end if
        wmate=intost(imate_tem)
        fil_funcc_tem=trim(fil_funcc_tem)//trim(wmate)
        call iofile(zero,lun_funcc_tem,trim(fil_funcc_tem),'TEMPER CP INTERPOLATION','old')
     end if

  end select


end subroutine tem_openfi

