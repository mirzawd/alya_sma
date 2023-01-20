!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_outwit()
  !------------------------------------------------------------------------
  !****f* chemic/chm_outwit
  ! NAME
  !    chm_outwit
  ! DESCRIPTION
  !    Output chemic species on witness points
  ! USES
  ! USED BY
  !    chm_output
  !***
  !------------------------------------------------------------------------
  use def_master,      only     : INOTMASTER, postp, witne, conce, therm
  use def_kermod,      only     : nwitn, lewit, shwit
  use def_domain,      only     : ltype, nnode, lnods
  use def_kintyp,      only     : ip, rp
  use def_chemic,      only     : React_ind, table_fw, nclas_chm, xZr_chm, xZs_chm, xYr_chm, xYs_chm, d32_chm, Sigma_chm,&
                                  Sigm0_chm, hrr_chm
  use mod_interp_tab,  only     : fw_scale_cont_var
  use mod_interp_tab,  only     : max_lookup_dim
  use mod_ker_proper,  only     : ker_proper

  implicit none
  integer(ip)       :: iwitn,ielem,inode,pnode,pelty,ipoin,iclas,idimt,ivawi,dummi

  real(rp)          :: control(max_lookup_dim)   ! input of table lookup function
  real(rp)          :: scale_control(max_lookup_dim)
  real(rp)          :: lim_control(max_lookup_dim,2_ip)
  integer(ip)       :: ind(max_lookup_dim)

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
     !
     ! Results on witness points
     !
     witne => postp(1) % witne

     if( INOTMASTER ) then

       do iwitn = 1,nwitn
          ielem = lewit(iwitn)
          if( ielem > 0 ) then
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             do ivawi = 1,postp(1) % nvawi
                witne(ivawi,iwitn) = 0.0_rp
             end do

             if( postp(1) % npp_witne(1) == 1 ) then
                do iclas = 1,min(8_ip,nclas_chm)
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      witne(iclas,iwitn) = witne(iclas,iwitn) + shwit(inode,iwitn) * conce(ipoin,iclas,1) * postp(1) % witne_dt(1)
                   end do
                end do
             end if

             if( postp(1) % npp_witne(8+1) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+1,iwitn) = witne(8+1,iwitn) + shwit(inode,iwitn) * xZr_chm(ipoin) * postp(1) % witne_dt(9)
                end do
             end if

             if( postp(1) % npp_witne(8+2) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+2,iwitn) = witne(8+2,iwitn) + shwit(inode,iwitn) * xZs_chm(ipoin) * postp(1) % witne_dt(10)
                end do
             end if

             if( postp(1) % npp_witne(8+3) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+3,iwitn) = witne(8+3,iwitn) + shwit(inode,iwitn) * xYr_chm(ipoin) * postp(1) % witne_dt(11)
                end do
             end if

             if( postp(1) % npp_witne(8+4) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+4,iwitn) = witne(8+4,iwitn) + shwit(inode,iwitn) * xYs_chm(ipoin) * postp(1) % witne_dt(12)
                end do
             end if

             if( postp(1) % npp_witne(8+5) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+5,iwitn) = witne(8+5,iwitn) + shwit(inode,iwitn) * d32_chm(ipoin) * postp(1) % witne_dt(13)
                end do
             end if

             if( postp(1) % npp_witne(8+6) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+6,iwitn) = witne(8+6,iwitn) + shwit(inode,iwitn) * Sigma_chm(ipoin)  * postp(1) % witne_dt(14)
                end do
             end if

             if( postp(1) % npp_witne(8+7) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+7,iwitn) = witne(8+7,iwitn) + shwit(inode,iwitn) * Sigm0_chm(ipoin) * postp(1) % witne_dt(15)
                end do
             end if

             if( postp(1) % npp_witne(8+17) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+17,iwitn) = witne(8+17,iwitn) + shwit(inode,iwitn) * hrr_chm(ipoin) * postp(1) % witne_dt(17)
                end do
             end if

             if( postp(1) % npp_witne(8+18) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+18,iwitn) = witne(8+18,iwitn) + shwit(inode,iwitn) * real(sum(React_ind(ipoin,:)),rp)&
                       * postp(1) % witne_dt(18)
                end do
             end if

             if( postp(1) % npp_witne(8+19) == 1 ) then
                !
                ! Density
                !
                call ker_proper('DENSI','IGAUS',dummi,ielem,witne(8+19,iwitn),pnode,1_ip,shwit(:,iwitn))   ! Fluid density
             end if


             if( postp(1) % npp_witne(8+8) == 1 ) then
                !
                ! Scaled control variables
                ! Assume CONCE is output as well, so the control variables are
                ! interpolated already except therm
                !
                control = 0.0_rp
                ind     = 1_ip
                do idimt = 1, table_fw % main_table % ndim
                   if (table_fw % kfl_chm_control(idimt) > 0) then
                      !
                      ! >0: one of the conces
                      !
                      control(idimt) = witne(table_fw % kfl_chm_control(idimt),iwitn)
                   else
                      if (table_fw % kfl_chm_control(idimt) == -1) then
                         !
                         ! -1: enthalpy
                         !
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            control(idimt) = control(idimt) + shwit(inode,iwitn) * therm(ipoin,1) * postp(1) % witne_dt(16)
                         enddo
                      elseif (table_fw % kfl_chm_control(idimt) == -2) then
                         !
                         ! -2: scalar dissipation rate
                         !
                         control(idimt) = witne(8+1,iwitn) +  witne(8+2,iwitn)
                      endif
                   endif
                enddo

                call fw_scale_cont_var( control, scale_control, lim_control, table_fw, ind)

                do idimt = 1, table_fw % main_table % ndim
                    witne(8+7+idimt,iwitn) = scale_control(idimt)
                end do
             end if


          end if

       end do


     end if

  end if

end subroutine chm_outwit

