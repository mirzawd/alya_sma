!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_outflow_mass.f90
!> @author  Guillaume Houzeaux
!> @brief   Outflew mass
!> @details Compute mass ath outflow (boundaries with code=20)
!>
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_outflow_mass()

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_communications, only : PAR_SUM
  use mod_bouder
  implicit none

  real(rp)    :: baloc(ndime,ndime)    
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: bovel(ndime,mnodb)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gbsur(mgaub),eucta  
  real(rp)    :: gbvel(3)
  integer(ip) :: ielem,ipoin,inode
  integer(ip) :: pnode,pgaus,iboun,igaub,inodb
  integer(ip) :: pelty,pblty,pnodb,pgaub,iflow

  outflow_mass = 0.0_rp

  if( INOTMASTER ) then

     boundaries: do iboun = 1,nboun

        if( kfl_fixbo_nsi(iboun) == 20 ) then               
           !
           ! Element properties and dimensions
           !
           pblty = ltypb(iboun) 
           pnodb = nnode(pblty)
           ielem = lelbo(iboun)
           pelty = ltype(ielem)

           if( pelty > 0 ) then

              pnode = nnode(pelty)
              pgaub = ngaus(pblty) 
              pgaus = ngaus(pelty)
              iflow = int(bvnat_nsi(4,iboun,1),ip)
              !
              ! Gather operations: ELVEL, ELCOD, BOVEL
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)    
              end do
              do inodb = 1,pnodb
                 ipoin                = lnodb(inodb,iboun)
                 bovel(1:ndime,inodb) = veloc(1:ndime,ipoin,1)
                 bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
              end do

              gauss_points: do igaub = 1,pgaub

                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                      bocod,baloc,eucta)                                   ! and Jacobian
                 gbsur(igaub) = elmar(pblty)%weigp(igaub)*eucta 
                 call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
                 !
                 !  +-
                 !  | u.n ds 
                 ! -+S
                 !
                 gbvel = 0.0_rp
                 do inodb = 1,pnodb
                    gbvel(1:ndime) = gbvel(1:ndime) + bovel(1:ndime,inodb) * elmar(pblty) % shape(inodb,igaub)
                 end do
                 outflow_mass(iflow) = outflow_mass(iflow) + dot_product(gbvel(1:ndime),baloc(1:ndime,ndime)) * gbsur(igaub)

              end do gauss_points

           end if

        end if

     end do boundaries

  end if

  call PAR_SUM(outflow_mass)

end subroutine nsi_outflow_mass
