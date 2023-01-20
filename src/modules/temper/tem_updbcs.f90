!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_updbcs
  ! NAME 
  !    tem_updbcs
  ! DESCRIPTION
  !    This routine updates the temperature boundary conditions:
  !    1. Before a time step begins
  !    2. Before a global iteration begins
  !    3. Before an inner iteration begins
  ! USED BY
  !    tem_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame,        only : pi
  use def_domain
  use def_kermod
  use def_temper
  use mod_physics,       only : physics_T_2_HCp
  use mod_ker_functions, only : ker_functions
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ibopo,iboun,ipnat,ifunc
  integer(ip)             :: itype
  real(rp)                :: tenew,venor,teold
  real(rp)                :: tbnew(npnat_tem)
  real(rp), external      :: funcre
  real(rp)                :: dummy
  real(rp), pointer       :: bvess_loc(:,:,:)

  select case ( itask )

  case ( ITASK_INIUNK , ITASK_BEGSTE )
     !  
     ! Before a time step.
     ! Update THERM(:,1), BVESS_TEM(1,:,1)
     !    
     if( kfl_intbc_tem /= 0 ) then
        call tem_intbcs()
     end if

     if( kfl_conbc_tem == 0 ) then
 
        if( kfl_regim_tem == 4 ) then
           bvess_loc => bvtem_tem
        else
           bvess_loc => bvess_tem
        endif
        do ipoin = 1,npoin

           ifunc = kfl_funno_tem(ipoin)
           itype = kfl_funtn_tem(ipoin)

           if( kfl_fixno_tem(1,ipoin) == 1 ) then

              call ker_functions(ipoin,ifunc,itype,bvess_loc(1,ipoin,2),tenew)

              if( itype /= 0 ) then
                if( kfl_timei_tem /= 0 .and. kfl_tiacc_tem == 2 .and. kfl_tisch_tem == 1 ) then
                    teold = therm(ipoin,ncomp_tem)
                    bvess_loc(1,ipoin,1) = 0.50_rp*(tenew+teold)
                 else
                    bvess_loc(1,ipoin,1) = tenew
                 end if                                  
              end if
              
           end if
        end do


        do iboun = 1,nboun
           
          ifunc = kfl_funbo_tem(iboun)
          itype = kfl_funtb_tem(iboun)
          
          if( kfl_fixbo_tem(iboun) /= 0 ) then

             call ker_functions(iboun,ifunc,itype,bvnat_tem(:,iboun,2),tbnew,WHEREIN='ON BOUNDARIES')
             
             !if( itype == FUNCTION_TIME ) then
             !   !
             !   ! Time function
             !   !
             !   do ipnat = 1,npnat_tem
             !      tenew = bvnat_tem(ipnat,iboun,2)&
             !           *funcre(funpa_tem(1,kfl_funbo_tem(iboun)),6,&
             !           kfl_funty_tem(kfl_funbo_tem(iboun)),cutim)
             !   end do                
             ! end if

              if( itype /= 0 ) then
                 if( kfl_timei_tem /= 0 .and. kfl_tiacc_tem == 2 .and. kfl_tisch_tem == 1 ) then
                    bvnat_tem(ipnat,iboun,1) = 0.50_rp*(tbnew(ipnat)+bvnat_tem(ipnat,iboun,1))
                 else
                    bvnat_tem(ipnat,iboun,1) = tbnew(ipnat)                       
                 end if
              end if
              
           end if
        end do

        if( kfl_exist_fixi7_tem == 1_ip .and. itask == ITASK_INIUNK ) then
           !
           ! Open or closes point depending on angle between normal and VELOC
           !
           call memgen(1_ip,npoin,0_ip)   !allocate gisca
           if( INOTMASTER )  call open_close(2_ip,kfl_fixno_tem,dummy,1_ip) ! calculates gisca(1:npoin)  
           do ipoin=1,npoin
              if(abs(kfl_fixno_tem(1,ipoin)) == 7_ip ) then
                 if( gisca(ipoin) > 0_ip ) then  
                    kfl_fixno_tem(1,ipoin) = 7_ip
                 else
                    kfl_fixno_tem(1,ipoin) = -7_ip
                 end if
            !     kfl_fixno_tem(1,ipoin) = 7_ip
              end if
           end do
           call memgen(3_ip,npoin,0_ip)   !deallocate gisca
        end if

     end if
     !
     ! Impose Dirichlet bc
     !
     if( kfl_regim_tem == 4_ip ) then
        call tem_calcEnthalpyBC()
        do ipoin = 1,npoin
           if( kfl_fixno_tem(1,ipoin) > 0 ) then
              tempe(ipoin,1) = bvtem_tem(1,ipoin,1)
              therm(ipoin,1) = bvess_tem(1,ipoin,1)
           end if
        end do
     else
        do ipoin = 1,npoin
           if( kfl_fixno_tem(1,ipoin) > 0 ) then
              therm(ipoin,1) = bvess_tem(1,ipoin,1)
           end if
        end do
     endif
     
  case ( ITASK_BEGITE )
     !
     ! Adaptive boundary condition
     !  
     do ipoin=1,npoin
        if(kfl_fixno_tem(1,ipoin)==4.or.kfl_fixno_tem(1,ipoin)==-4) then
           ibopo=lpoty(ipoin)
           if(ibopo>0) then
              !venor=dot_product(veloc(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo))
              venor=dot_product(advec(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo))
              if(venor<=0.0_rp) then
                 kfl_fixno_tem(1,ipoin)= 4
              else
                 kfl_fixno_tem(1,ipoin)=-4
              end if
           else
              kfl_fixno_tem(1,ipoin)=-4
           end if
        end if
     end do

  end select

end subroutine tem_updbcs
