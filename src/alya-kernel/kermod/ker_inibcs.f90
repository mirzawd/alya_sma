!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_inibcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Impose boundary conditions
!> @details Impose boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine ker_inibcs()
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_opebcs
  use mod_tubes
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  integer(ip)  :: ipoin,iboun,inodb,dummi
!  integer(ip)  :: ielem

  !-------------------------------------------------------------
  !
  ! WALLD: wall distance
  !
  !-------------------------------------------------------------

  if( kfl_walld /= 0 .and. solve(2) % kfl_algso /= -999 ) then

     dummi = 0
     
     if( INOTEMPTY ) then
        !
        ! Allocate memory
        !
        call ker_memory(1_ip)       ! KFL_FIXNO_WALLD_KER
        call ker_memory(2_ip)       ! KFL_FIXBO_WALLD_KER
        !
        ! Boundary codes
        !
        if( kfl_icodb > 0 .and. nboun > 0 ) then
           iffun     =  0
           kfl_fixbo => kfl_fixbo_walld_ker
           tbcod     => tbcod_ker(2:) 
           call memgen(0_ip,1_ip,nboun)
           bvnat     => gevec
           call reacod(20_ip)
           call memgen(2_ip,1_ip,nboun)
           do iboun = 1,nboun
              if( kfl_fixbo_walld_ker(iboun) == 1 ) then
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    kfl_fixno_walld_ker(1,ipoin) = 1
                 end do
              end if
           end do
        end if
        call ker_memory(-2_ip)      ! KFL_FIXBO_WALLD_KER
        !
        ! Node codes
        !
        if( kfl_icodn > 0 ) then 
           iffun     =  0
           ifbop     =  0
           ifloc     =  0
           ifbes     =  0
           kfl_fixno => kfl_fixno_walld_ker
           tncod     => tncod_ker(2:)
           call reacod(IMPOSE_NODE_CODES)
        end if

        call PAR_INTERFACE_NODE_EXCHANGE(kfl_fixno_walld_ker,'SUM')
        do ipoin = 1,npoin
           if( kfl_fixno_walld_ker(1,ipoin) /= 0 ) dummi = dummi + 1
           kfl_fixno_walld_ker(1,ipoin) = min(kfl_fixno_walld_ker(1,ipoin),1_ip)
        end do        
        !
        ! Wall Node numbers
        !
        kfl_nwall = 0
        do ipoin = 1,npoi1
          if (kfl_fixno_walld_ker(1,ipoin) == 1) then
            kfl_nwall = kfl_nwall + 1
          endif  
        enddo
        do ipoin = npoi2,npoi3
          if (kfl_fixno_walld_ker(1,ipoin) == 1) then
            kfl_nwall = kfl_nwall + 1
          endif  
        enddo
        
     end if
     
     call PAR_SUM(kfl_nwall,'IN MY CODE' )
     call PAR_SUM(dummi,    'IN MY CODE' )

     if( dummi == 0 ) then
        !do ipoin = 1,npoin
        !   ibopo = lpoty(ipoin)
        !   if( ibopo > 0 ) then
        !      kfl_fixno_walld_ker(1,ipoin) = 1
        !   end if
        !end do
        call runend('KER_INIBCS: NO BOUNDARY CONDITIONS PRESCRIBED FOR WALL DISTANCE')
     end if

  end if

  !-------------------------------------------------------------
  !
  ! WALLN: wall normal
  !
  !-------------------------------------------------------------

  if( kfl_walln /= 0 .and. solve(5) % kfl_algso /= -999 ) then

     dummi = 0

     if( INOTEMPTY ) then
        !
        ! Allocate memory
        !
        call ker_memory(10_ip)       ! KFL_FIXNO_WALLN_KER
        call ker_memory(11_ip)       ! KFL_FIXBO_WALLN_KER
        !
        ! Boundary codes
        !
        if( kfl_icodb > 0 .and. nboun > 0 ) then
           iffun     =  0
           kfl_fixbo => kfl_fixbo_walln_ker
           tbcod     => tbcod_ker(5:) 
           call memgen(0_ip,1_ip,nboun)
           bvnat     => gevec
           call reacod(20_ip)
           call memgen(2_ip,1_ip,nboun)
           do iboun = 1,nboun
              if( kfl_fixbo_walln_ker(iboun) == 1 ) then
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    kfl_fixno_walln_ker(1,ipoin) = 1
                 end do
              end if
           end do
        end if
        call ker_memory(-11_ip)      ! KFL_FIXBO_WALLN_KER
        !
        ! Node codes
        !
        if( kfl_icodn > 0 ) then 
           iffun     =  0
           ifbop     =  0
           ifloc     =  0
           ifbes     =  0
           kfl_fixno => kfl_fixno_walln_ker
           tncod     => tncod_ker(5:)
           call reacod(IMPOSE_NODE_CODES)
        end if

        call PAR_INTERFACE_NODE_EXCHANGE(kfl_fixno_walln_ker,'SUM')
        
        do ipoin = 1,npoin
           if( kfl_fixno_walln_ker(1,ipoin) /= 0 ) dummi = dummi + 1
           kfl_fixno_walln_ker(1,ipoin) = min(kfl_fixno_walln_ker(1,ipoin),1_ip)
        end do

     end if

     call PAR_SUM(dummi)

     if( dummi == 0 ) then
        call runend('KER_INIBCS: NO BOUNDARY CONDITIONS PRESCRIBED FOR WALL NORMAL')
     end if

  end if

  !-------------------------------------------------------------
  !
  ! ROUGH: roughness extension
  !
  !-------------------------------------------------------------

  if( kfl_rough /= -1 .and. kfl_extro > 0 ) then

     dummi = 0

     if( INOTEMPTY ) then
        !
        ! Allocate memory
        !
        call ker_memory(3_ip)       ! KFL_FIXNO_ROUGH_KER
        call ker_memory(4_ip)       ! KFL_FIXBO_ROUGH_KER
        !
        ! Boundary codes
        !
        if( kfl_icodb > 0 .and. nboun > 0 ) then
           iffun     =  0
           kfl_fixbo => kfl_fixbo_rough_ker
           tbcod     => tbcod_ker(1:) 
           call memgen(0_ip,1_ip,nboun)
           bvnat     => gevec
           call reacod(20_ip)
           call memgen(2_ip,1_ip,nboun)
           do iboun = 1,nboun
              if( kfl_fixbo_rough_ker(iboun) == 1 ) then
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    kfl_fixno_rough_ker(1,ipoin) = 1
                 end do
              end if
           end do
        end if
        call ker_memory(-4_ip)      ! KFL_FIXBO_ROUGH_KER
        !
        ! Node codes
        !
        if( kfl_icodn > 0 ) then 
           iffun     =  0
           ifbop     =  0
           ifloc     =  0
           ifbes     =  0
           kfl_fixno => kfl_fixno_rough_ker
           tncod     => tncod_ker(2:)
           call reacod(IMPOSE_NODE_CODES)
        end if

        call PAR_INTERFACE_NODE_EXCHANGE(kfl_fixno_rough_ker,'SUM')
        
        do ipoin = 1,npoin
           if( kfl_fixno_rough_ker(1,ipoin) /= 0 ) dummi = dummi + 1
           kfl_fixno_rough_ker(1,ipoin) = min(kfl_fixno_rough_ker(1,ipoin),1_ip)
        end do

     end if

     call PAR_SUM(dummi)

     if( dummi == 0 ) call runend('KER_WALGEN: NO BOUNDARY CONDITIONS PRESCRIBED FOR ROUGHNESS EXTENSION')

  end if
  
  !-------------------------------------------------------------
  !
  ! NO_SL_WALL_LAW: boundary code for no slip wall law
  !
  !-------------------------------------------------------------

  if( kfl_noslw_ker > 0 ) then

     call ker_memory(12_ip)       ! KFL_FIXBO_NSW_KER 
     dummi = 0

     if( INOTEMPTY ) then
        !
        ! Allocate memory
        !
        !        call ker_memory(12_ip)       ! KFL_FIXBO_NSW_KER
        !        call memory_alloca(mem_modul(1:2,modul),'LNSW_EXCH','ker_memory' , lnsw_exch , max(1_ip,nelem) )

        !allocate(lnsw_exch(max(1_ip,nelem) ))
        !do ielem = 1,nelem
        !   lnsw_exch(ielem) % velav = 0.0_rp
        !end do
        !
        ! Boundary codes
        !
        if( kfl_icodb > 0 .and. nboun > 0 ) then
           iffun     =  0
           kfl_fixbo => kfl_fixbo_nsw_ker
           tbcod     => tbcod_ker(7:) 
           call memgen(0_ip,1_ip,nboun)
           bvnat     => gevec
           call reacod(20_ip)
           call memgen(2_ip,1_ip,nboun)
           dummi = dummi + 1
        end if
!        call ker_memory(-12_ip)      ! KFL_FIXBO_NSW_KER  - no lo dealoco lo necesito!!
        
     end if
     
     call PAR_SUM(dummi)

     if( dummi == 0 ) call runend('KER_WALGEN: NO BOUNDARY CONDITIONS PRESCRIBED FOR NO SLIP WALL LAW')

  end if

  !-------------------------------------------------------------
  !
  ! Mesh deformation
  !
  !-------------------------------------------------------------

  if( kfl_defor /= 0 .and. INOTEMPTY ) then
     !
     ! Allocate memory: KFL_FIXNO_DEFOR_KER and BVESS_DEFOR_KER
     !
     call ker_memory(9_ip) 
     !
     ! Node codes
     !
     if( kfl_icodn > 0 ) then 
        iffun     =  0
        ifbop     =  0
        ifloc     =  0
        kfl_fixno => kfl_fixno_defor_ker
        bvess     => bvess_defor_ker
        tncod     => tncod_ker(3:)
        call reacod(IMPOSE_NODE_CODES)
     end if

  end if

  !-------------------------------------------------------------
  !
  ! SUPPO: Support geometry
  !
  !-------------------------------------------------------------

  if( kfl_suppo /= 0 .and. INOTEMPTY ) then
     
     call ker_memory(6_ip)
     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then
           kfl_fixno_suppo_ker(1:ndime,ipoin) = 1
        end if
     end do
     
  end if

end subroutine ker_inibcs
