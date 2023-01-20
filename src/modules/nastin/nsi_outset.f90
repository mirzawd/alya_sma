!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outset()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_outset
  ! NAME 
  !    nsi_outset
  ! DESCRIPTION
  !    Compute and write results on sets:
  !    - Element sets:
  !      1.    VELOC: mean velocity =  sqrt{ (int_W u^2 dw)/(2*meas(W)) }
  !      2.    VORTI: mean vorticity = sqrt{ (int_W w^2 dw)/(2*meas(W)) }
  !                                    with w=dv/dx-du/dy
  !      3.    VORTI: mean vorticity = int_W rho*u^2/2 dw
  !
  !    - Boundary sets:
  !      1.    MEANP: mean pressure  =  int_W p/meas(W)
  !      2.    MASS:  mass           =  int_W rho*u.n 
  !      3-8.  FORCE: force          =  int_W [-pI+2 mu E(u)].n 
  !      9-14. TORQU: torque         =  int_W r x ( [-pI+2 mu E(u)].n )
  !      15.   MEANY: mean y+        =  int_W y+/meas(W)
  !      Force and torque are exerted by the solid on the fluid.
  !
  !   - Node sets:
  !      1-3.  VELOC: velocity components
  !      4.    PRESS: pressure
  !      5.    YPLUS: y+
  ! USES
  ! USED BY
  !    nsi_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_iofile
  use mod_communications, only : PAR_BROADCAST
  use mod_ker_nsw_visc2,  only : ker_nod_nsw_visc_0
  use mod_frivel,         only : frivel
  use mod_wall_exchange,  only : ker_waexlo_getval
  use mod_output_postprocess, only : output_postprocess_node_sets_parall
  use mod_output_postprocess, only : output_postprocess_boundary_sets_parall
  use mod_output_postprocess, only : output_postprocess_element_sets_parall
  implicit none
  integer(ip) :: keset,kbset,knset,idime,ipoin,ibopo
!  integer(ip) :: iiset
  real(rp)    :: tveno,ustar,yplus,nu,ubulk

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setse) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepelset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef) ) then
         if( INOTMASTER ) then
            do keset = 1,neset
               call nsi_elmset(lesec(keset),keset)
            end do
         end if
         call output_postprocess_element_sets_parall()
         !
         ! Normalize results of necessary
         !
         do keset = 1,neset
            if(postp(1) % npp_setse(1)/=0.and.veset(postp(1) % nvaes+1,keset)>0.0_rp) then
               veset(1,keset)=sqrt(veset(1,keset)/(2.0_rp*veset(postp(1) % nvaes+1,keset)))
            end if
            if(postp(1) % npp_setse(2)/=0.and.veset(postp(1) % nvaes+1,keset)>0.0_rp) then
               veset(2,keset)=sqrt(veset(2,keset)/(2.0_rp*veset(postp(1) % nvaes+1,keset)))
            end if
            if(postp(1) % npp_setse(8)/=0.and.veset(postp(1) % nvaes+1,keset)>0.0_rp) then
               veset(8,keset)=veset(8,keset)/(veset(postp(1) % nvaes+1,keset))
            end if
         end do
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepboset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef) ) then
         if ( kfl_waexl_ker == 1_ip ) then  ! Wall with exchange location - this was inside bouope - moved outside so that master enters
            !
            ! Interpolate
            !
            call ker_waexlo_getval(veloc,velel_ker)            
         end if
         if( kfl_noslw_ker /= 0 ) call ker_nod_nsw_visc_0   ! obtains el_nsw_visc(1:nelem)
         do kbset = 1,nbset 
            call nsi_bouset(lbsec(kbset),kbset)
         end do
         call output_postprocess_boundary_sets_parall() ! PAR_SUM of postp(1)%vbset
         !
         ! Normalize results of necessary
         !
         do kbset = 1,nbset
            if(postp(1) % npp_setsb(1)/=0.and.vbset(postp(1) % nvabs+1,kbset)>0.0_rp) then
               vbset( 1,kbset)=vbset( 1,kbset)/vbset(postp(1) % nvabs+1,kbset)
            end if
            if(postp(1) % npp_setsb(15)/=0.and.vbset(postp(1) % nvabs+1,kbset)>0.0_rp) then
               vbset(15,kbset)=vbset(15,kbset)/vbset(postp(1) % nvabs+1,kbset)
            end if
            if(postp(1) % npp_setsb(16)/=0.and.vbset(postp(1) % nvabs+1,kbset)>0.0_rp) then
               vbset(16,kbset)=vbset(16,kbset)/vbset(postp(1) % nvabs+1,kbset)
            end if
            if(postp(1) % npp_setsb(35)/=0.and.vbset(postp(1) % nvabs+1,kbset)>0.0_rp) then
               vbset(35,kbset)=vbset(35,kbset)/vbset(postp(1) % nvabs+1,kbset)
            end if
         end do
         !
         ! Mass flow rate control
         !
         if ( kfl_mfrco_nsi > 0 ) then

            ubulk = vbset(16,mfrse_nsi)
            if (IPARALL) call PAR_BROADCAST(ubulk,'IN MY ZONE')

            if ( INOTMASTER ) then

               if ( ittim .gt. 1) then
                     grnor_nsi = grnor_nsi + mfccf_nsi*((mfrub_nsi) + 2.0_rp*(ubulk) -(ubpre_nsi))/dtime
                     ubpre_nsi = ubulk
               endif

            end if

         end if
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then
      if(( mod(ittim, postp(1) % npp_stepnoset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef) ) then
         if( INOTMASTER ) then
            do knset = 1,nnset

               if(lnsec(knset)/=0) then
                  do idime=1,ndime
                     if(postp(1) % npp_setsn(idime)/=0)&
                           vnset(idime,knset)=veloc(idime,lnsec(knset),1)
                  end do
                  if(ndime==2) vnset(3,knset)=0.0_rp
                  if(postp(1) % npp_setsn(4)/=0) vnset(4,knset)=press(lnsec(knset),1)
                  if(postp(1) % npp_setsn(5)/=0) then
                     ipoin=lnsec(knset)
                     ibopo=lpoty(ipoin)
                     if(ibopo==0) then
                        yplus=0.0_rp
                     else
                        call vecnor(veloc(1,ipoin,1),ndime,tveno,2_ip)
                        nu=visco(ipoin,1)/densi(ipoin,1)
                        if( kfl_rough > 0 ) rough_dom = rough(ipoin)  ! non uniform roughness
                        if (kfl_delta == 1) then
                           call frivel(kfl_ustar,ywalp(ibopo),rough_dom,tveno,nu,ustar)
                           yplus = ywalp(ibopo)*ustar/nu
                        else
                           call frivel(kfl_ustar,delta_nsi,rough_dom,tveno,nu,ustar)
                           yplus = delta_nsi*ustar/nu
                        end if
                     end if
                     vnset(5,knset)=yplus
                  end if
                  if(postp(1) % npp_setsn(6)/=0) then
                     call nsi_outpmv(-1_ip,lnsec(knset),lnsec(knset),vnset(6,knset))
                  end if
                  if(postp(1) % npp_setsn(7)/=0) then
                     call nsi_outpmv(-2_ip,lnsec(knset),lnsec(knset),vnset(7,knset))
                  end if
               end if
            end do
         end if
         call output_postprocess_node_sets_parall()
      end if
  end if

  !----------------------------------------------------------------------
  !
  ! Immersed boundary sets
  !
  !----------------------------------------------------------------------

  !if(maxval(npp_setsb_nsi)>0.and.niset>0) then
  !   !
  !   ! Header
  !   !
  !   call outset(&
  !        4_ip,lun_setsi_nsi,nvars_nsi,npp_setsb_nsi,wobse_nsi,viset_nsi)
  !   !
  !   ! Variables on sets
  !   !
  !   if( INOTMASTER ) then
  !      do iiset=1,nimbo
  !         call nsi_bkbset(&
  !              2_ip,&
  !              iiset,viset_nsi(nvars_nsi+1,iiset),&
  !              viset_nsi(1,iiset),viset_nsi(2,iiset),viset_nsi(3,iiset),&
  !              viset_nsi(9,iiset),viset_nsi(15,iiset),viset_nsi(16,iiset))
  !      end do
  !   end if
  !   if( IPARALL ) then
  !      nparr =  niset*(nvars_nsi+1)
  !      call par_operat(3_ip)
  !   end if
  !   !
  !   ! Normalize results of necessary
  !   !
  !   if( INOTSLAVE ) then
  !      do iiset=1,niset
  !         if(npp_setsb_nsi(1)/=0.and.viset_nsi(nvars_nsi+1,iiset)>0.0_rp) then
  !            viset_nsi( 1,iiset)=viset_nsi( 1,iiset)/viset_nsi(nvars_nsi+1,iiset)
  !         end if
  !         if(npp_setsb_nsi(15)/=0.and.viset_nsi(nvars_nsi+1,iiset)>0.0_rp) then
  !            viset_nsi(15,iiset)=viset_nsi(15,iiset)/viset_nsi(nvars_nsi+1,iiset)
  !         end if
  !         if(npp_setsb_nsi(16)/=0.and.viset_nsi(nvars_nsi+1,iiset)>0.0_rp) then
  !            viset_nsi(16,iiset)=viset_nsi(16,iiset)/viset_nsi(nvars_nsi+1,iiset)
  !         end if
  !      end do
  !   end if
  !   !
  !   ! Write results
  !   !
  !   call outset(&
  !        40_ip,lun_setsi_nsi,nvars_nsi,npp_setsb_nsi,wobse_nsi,viset_nsi)!
  !
  !end if

end subroutine nsi_outset
