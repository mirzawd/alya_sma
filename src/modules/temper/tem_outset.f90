!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outset()

  !-----------------------------------------------------------------------
  !****f* Temper/tem_outset
  ! NAME 
  !    tem_outset
  ! DESCRIPTION
  !    Compute and write results on sets
  ! USES
  ! USED BY
  !    tem_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_memory, only : memory_alloca, memory_deallo
  use mod_output_postprocess, only : output_postprocess_node_sets_parall
  use mod_output_postprocess, only : output_postprocess_boundary_sets_parall
  use mod_output_postprocess, only : output_postprocess_element_sets_parall
  implicit none
  integer(ip) :: ieset,ibset,inset

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setse) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepelset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
         if( INOTMASTER ) then
            do ieset=1,neset        
               call tem_elmset( lesec(ieset), ieset)
            end do
         end if
         call output_postprocess_element_sets_parall()
         
         !
         ! Normalize results if necessary
         !
         if( INOTSLAVE ) then
         endif
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb)>0 ) then
      if(( mod(ittim, postp(1) % npp_stepboset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
         if( INOTMASTER ) then
            if( postp(1) % npp_setsb(3) /= 0 ) then
               call memory_alloca(mem_modul(1:2,modul),'HEATF_TEM','tem_outset',heatf_tem,nboun)
               call memory_alloca(mem_modul(1:2,modul),'MASSB_TEM','tem_outset',massb_tem,nboun)
               if (ittim.gt.0) call tem_bouope(2_ip) ! return total heat flux heatf_tem and surface massb_Tem
            end if
            
            do ibset=1,nbset
               call tem_bouset(lbsec(ibset),ibset)
            end do

            if( postp(1) % npp_setsb(3) /= 0 ) then
               call memory_deallo(mem_modul(1:2,modul),'HEATF_TEM','tem_outset',heatf_tem)
               call memory_deallo(mem_modul(1:2,modul),'MASSB_TEM','tem_outset',massb_tem)
            end if
         end if
         call output_postprocess_boundary_sets_parall()
         if(postp(1) % npp_setsb(1)/=0) then
            do ibset=1,nbset
               if (postp(1) % vbset(postp(1) % nvabs+1,ibset)/=0.0_rp) then 
                     postp(1) % vbset(1,ibset)=postp(1) % vbset(1,ibset)&
                     /postp(1) % vbset(postp(1) % nvabs+1,ibset)
                     !
                     ! Density:
                     !
                     postp(1) % vbset(4,ibset)=postp(1) % vbset(4,ibset)&
                     /postp(1) % vbset(postp(1) % nvabs+1,ibset)

               endif
            end do   
         end if
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn)>0 ) then
      if(( mod(ittim, postp(1) % npp_stepnoset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
         if( INOTMASTER ) then
            do inset=1,nnset
               if(lnsec(inset)/=0) then
                  if(postp(1) % npp_setsn(1)/=0) postp(1) % vnset(1,inset)=therm(lnsec(inset),1)
               end if
            end do
         end if
         call output_postprocess_node_sets_parall()
      end if
  end if

end subroutine tem_outset
