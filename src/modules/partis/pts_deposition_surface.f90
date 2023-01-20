!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_deposition_surface.f90
!> @author  houzeaux
!> @date    2018-06-11
!> @brief   Surface of deposition
!> @details Surface of deposition
!> @}
!-----------------------------------------------------------------------

subroutine pts_deposition_surface(depos_surface)

  use def_master
  use def_domain
  use def_partis
  use mod_communications, only : PAR_SUM
  use mod_bouder
 
  implicit none
  real(rp), intent(out) :: depos_surface(*)        !< Surface deposition
  integer(ip)           :: iboun,pblty,pnodb,igaub
  integer(ip)           :: inodb,itype,ipoin
  real(rp)              :: gbsur,eucta,baloc(3,3)
  real(rp)              :: bocod(ndime,mnodb)
  
  depos_surface(1:ntyla_pts+1) = 0.0_rp

  if( kfl_depos_surface_pts == 1 ) then
    
     do iboun = 1,nboun

        if( maxval(depob_pts(:,iboun)) > 0.0_rp ) then
           pblty = ltypb(iboun)
           pnodb = nnode(pblty)
           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
           end do

           gauss_points: do igaub=1,ngaus(pblty)
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),&
                   bocod,baloc,eucta)
              gbsur = elmar(pblty) % weigp(igaub) * eucta
              do itype = 1,ntyla_pts
                 if( parttyp(itype) % kfl_exist /= 0 ) then                
                    if( depob_pts(itype,iboun) > 0.0_rp ) then
                       depos_surface(itype)       = depos_surface(itype)       + gbsur
                       depos_surface(ntyla_pts+1) = depos_surface(ntyla_pts+1) + gbsur
                    end if
                 end if
              end do
           end do gauss_points

        end if

     end do

     call PAR_SUM(ntyla_pts+1_ip,depos_surface)
    
  end if

end subroutine pts_deposition_surface
