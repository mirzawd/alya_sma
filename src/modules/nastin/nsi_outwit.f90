!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outwit()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_outwit
  ! NAME 
  !    nsi_outwit
  ! DESCRIPTION
  !    Output results on witness points
  ! USES
  ! USED BY
  !    nsi_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper,     only : ker_proper
  use mod_communications, only : PAR_MAX,PAR_MIN
  use mod_witness,        only : WITNESS_INSTANTANEAOUS
  implicit none
  integer(ip) :: iwitn,ielem,inode,pnode,pelty,ipoin
  real(rp)    :: raux

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0.and. &
       ( postp(1) % witne_kfl_time(1)/= WITNESS_INSTANTANEAOUS.or.mod(ittim, postp(1) % npp_stepw) == 0)  ) then

     if( INOTMASTER ) then 

        do iwitn = 1, nwitn
           ielem = lewit(iwitn)
           if( ielem > 0 ) then

              pelty = ltype(ielem)
              pnode = nnode(pelty)

              if( postp(1) % npp_witne(1) == 1 ) then
                 !
                 ! u
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) * veloc(1,ipoin,1) * postp(1) % witne_dt(1)
                 end do
              end if

              if( postp(1) % npp_witne(2) == 1 ) then
                 !
                 ! v
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(2,iwitn) = witne(2,iwitn) + shwit(inode,iwitn) * veloc(2,ipoin,1) * postp(1) % witne_dt(2)
                 end do
              end if

              if( postp(1) % npp_witne(3) == 1 ) then
                 !
                 ! w
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(3,iwitn) = witne(3,iwitn) + shwit(inode,iwitn) * veloc(ndime,ipoin,1) * postp(1) % witne_dt(3)
                 end do
              end if

              if( postp(1) % npp_witne(4) == 1 ) then
                 !
                 ! P
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(4,iwitn) = witne(4,iwitn) + shwit(inode,iwitn) * press(ipoin,1) * postp(1) % witne_dt(4)
                 end do
              end if

              if( postp(1) % npp_witne(5) == 1 ) then  
                 !
                 ! S11
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(5,iwitn) = witne(5,iwitn) + dewit(1,inode,iwitn) * veloc(1,ipoin,1) * postp(1) % witne_dt(5)
                 end do
              end if

              if( postp(1) % npp_witne(6) == 1 ) then  
                 !
                 ! S22
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(6,iwitn) = witne(6,iwitn) + dewit(2,inode,iwitn) * veloc(2,ipoin,1) * postp(1) % witne_dt(6)
                 end do
              end if

              if( postp(1) % npp_witne(7) == 1 ) then  
                 !
                 ! S12
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(7,iwitn) = witne(7,iwitn) &
                         + 0.5_rp*(   dewit(2,inode,iwitn) * veloc(1,ipoin,1) &
                         &          + dewit(1,inode,iwitn) * veloc(2,ipoin,1) ) * postp(1) % witne_dt(7)
                 end do
              end if

              if( postp(1) % npp_witne(8) == 1 ) then  
                 !
                 ! S33
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(8,iwitn) = witne(8,iwitn) + dewit(ndime,inode,iwitn) * veloc(ndime,ipoin,1) * postp(1) % witne_dt(8)
                 end do
              end if

              if( postp(1) % npp_witne(9) == 1 ) then  
                 !
                 ! S13
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(9,iwitn) = witne(9,iwitn) &
                         + 0.5_rp*(   dewit(ndime,inode,iwitn) * veloc(1,ipoin,1) &
                         &          + dewit(1,inode,iwitn)     * veloc(ndime,ipoin,1) ) * postp(1) % witne_dt(9)
                 end do
              end if

              if( postp(1) % npp_witne(10) == 1 ) then  
                 !
                 ! S23
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(10,iwitn) = witne(10,iwitn) &
                         + 0.5_rp*(   dewit(ndime,inode,iwitn) * veloc(2,ipoin,1) &
                         &          + dewit(2,inode,iwitn)     * veloc(ndime,ipoin,1) ) * postp(1) % witne_dt(10)
                 end do
              end if

              if( postp(1) % npp_witne(11) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(11,iwitn) = witne(11,iwitn) + shwit(inode,iwitn) * coord(1,ipoin)
                 end do
              end if

              if( postp(1) % npp_witne(12) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(12,iwitn) = witne(12,iwitn) + shwit(inode,iwitn) * coord(2,ipoin)
                 end do
              end if

              if( postp(1) % npp_witne(13) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(13,iwitn) = witne(13,iwitn) + shwit(inode,iwitn) * coord(ndime,ipoin)
                 end do
              end if
              !
              ! TURBU ! postproces nu_t!!!
              !              
              if( postp(1) % npp_witne(14) == 1 ) then
                 call ker_proper('TURBU','IGAUS',1_ip,ielem,witne(14,iwitn),pnode,1_ip,shwit(:,iwitn))
                 witne(14,iwitn) = witne(14,iwitn) * postp(1) % witne_dt(14)
              end if
 
              if( postp(1) % npp_witne(15) == 1 ) then
                 !
                 ! VELMX
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(15,iwitn) = witne(15,iwitn) + shwit(inode,iwitn) * velom(1,ipoin) * postp(1) % witne_dt(15)
                 end do
              end if
              
              if( postp(1) % npp_witne(16) == 1 ) then
                 !
                 ! VELMY
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(16,iwitn) = witne(16,iwitn) + shwit(inode,iwitn) * velom(2,ipoin) * postp(1) % witne_dt(16)
                 end do
              end if

              if( postp(1) % npp_witne(17) == 1 .and. ndime == 3_ip ) then
                 !
                 ! VELMZ
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(17,iwitn) = witne(17,iwitn) + shwit(inode,iwitn) * velom(3,ipoin) * postp(1) % witne_dt(17)
                 end do
              end if
              
           end if
        end do

     end if
     !
     ! Convergence with witness
     !
     if ( kfl_stop_by_wit_nsi == 1 ) call nsi_convergence_witness()
  end if
  !
  ! Sends to slaves so that everybody stops correctly
  !
  if( kfl_stop_by_wit_nsi == 1 ) then
     !
     ! PAR_MIN & MAX do not have the option INCLUDE_MASTER for integers
     ! It would be difficult to add it (see my file Bad_programming_style).
     ! Therefore I will solve it with a chapuza (convert to real and back) 
     !
     raux = real(kfl_stead_nsi,rp)
     call PAR_MAX(raux,'IN MY CODE',INCLUDE_ROOT=.true.)
     kfl_stead_nsi = nint(raux,ip)

     raux = real(kfl_gotim,rp)
     call PAR_MIN(raux,'IN MY CODE',INCLUDE_ROOT=.true.)
     kfl_gotim = nint(raux,ip)
  end if

end subroutine nsi_outwit

