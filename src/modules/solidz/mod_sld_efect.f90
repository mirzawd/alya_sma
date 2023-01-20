!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_efect
!> @author  David Oks 
!> @date    2022-01-12
!> @brief   Embedded Finite Element Coupling Technique (EFECT) for FSI
!> @details Embedded Finite Element Coupling Technique (EFECT) for FSI
!> @}
!------------------------------------------------------------------------

module mod_sld_efect

  use def_kintyp, only : ip, rp, lg
  use def_master, only : displ, INOTMASTER, TIME_N, ITER_K
  use def_domain

  implicit none

  public                            :: &
       sld_compute_total_force

contains
  
  
  subroutine sld_compute_total_force(icoup, force, kfl_how_forces, totfo)
  
    use def_coupli,         only : coupling_type
    use mod_communications, only : PAR_SUM
  
    implicit none
  
    integer(ip), intent(in)  :: icoup
    integer(ip), intent(in)  :: kfl_how_forces 
    real(rp),    intent(in)  :: force(ndime,*)
    real(rp),    intent(out) :: totfo(ndime)
    integer(ip)              :: ipoin, iboun, kboun, idime, inodb, pnodb, bcode_target, pblty
  
    totfo = 0.0_rp
  
    if (kfl_how_forces == 1_ip) then
       ! Compute total force on nodes
       if (INOTMASTER) then
          ! Identify target coe
          bcode_target = coupling_type(icoup) % where_number
          ! Loop over wet boundaries
          do kboun = 1,coupling_type(icoup) % wet % nboun_wet
             iboun = coupling_type(icoup) % wet % lboun_wet(kboun)
             ! If wet boundary is the target boundary => add to total force
             if( kfl_codbo(iboun) == bcode_target ) then
                pblty = ltypb(iboun)
                pnodb = nnode(pblty)
                ! Loop over nodes
                do inodb = 1,pnodb
                   ipoin = lnodb(inodb,iboun)
                   ! Only consider my own subdomain nodes
                   if (ipoin <= npoi1 .or. (ipoin >= npoi2 .and. ipoin <= npoi3)) then
                      ! Loop over dimensions
                      do idime = 1,ndime
                         ! Add contribution to total force
                         totfo(idime) = totfo(idime) + force(idime,ipoin)
                      end do
                   end if
                end do 
             end if
          end do
       end if
       ! MPI Sum
       call PAR_SUM(ndime,totfo,'IN MY CODE')
       !
    else if (kfl_how_forces == 2_ip) then
       ! Compute total force on boundary elements
       if (INOTMASTER) then
          do iboun=1,nboun
             do idime = 1,ndime
                totfo(idime) = totfo(idime) + force(idime,iboun)
             end do
          end do
       end if
       ! MPI Sum
       call PAR_SUM(ndime,totfo,'IN MY CODE')
       !
    end if
  
  end subroutine sld_compute_total_force


  subroutine sld_integrate_tractions(tract, stres, kfl_how_forces)

  use def_solidz,        only :  kfl_fixno_sld
  use def_solidz,        only :  kfl_fixbo_sld
  use mod_bouder

  implicit none 
  
  real(rp),    intent(inout) :: tract(ndime,*)
  real(rp),    intent(in)    :: stres(ntens,*)
  integer(ip), intent(in)    :: kfl_how_forces 
  real(rp)                   :: baloc(ndime,ndime)
  real(rp)                   :: bocod(ndime,mnodb)
  real(rp)                   :: elcod(ndime,mnode)
  real(rp)                   :: streb(ndime,ndime)
  real(rp)                   :: tractb(ndime)
  real(rp)                   :: eucta
  integer(ip)                :: ltens(ndime,ndime) 
  integer(ip)                :: pblty, pnodb, pgaub, ielem, pelty, pnode 
  integer(ip)                :: iboun, igaub, inodb, itens, inode 
  integer(ip)                :: idime, jdime
  integer(ip)                :: ipoin
    !
    ! Translate tensor to matrix indices
    !
    if ( ndime == 2 ) then
       ltens(1,1) = 1 
       ltens(2,2) = 2 
       ltens(1,2) = 3
       ltens(2,1) = 3
    else if ( ndime == 3 ) then
       ltens(1,1) = 1 
       ltens(2,2) = 2 
       ltens(1,2) = 3
       ltens(2,1) = 3
       ltens(3,3) = 4
       ltens(1,3) = 5 
       ltens(3,1) = 5 
       ltens(2,3) = 6 
       ltens(3,2) = 6
    end if
    !
    ! Compute fluid tractions at solid boundary 
    !
    boundaries: do iboun = 1,nboun

       pblty = ltypb(iboun)
       pnodb = nnode(pblty)
       pgaub = ngaus(pblty)
       ielem = lelbo(iboun)
       pelty = ltype(ielem)
       pnode = nnode(pelty)
       streb = 0.0_rp

       gauss_points: do igaub = 1,pgaub
          !
          ! Interpolate stresses to solid boundary gauss points
          !
          do inodb = 1,pnodb
             ipoin = lnodb(inodb,iboun)
             do idime = 1,ndime
                bocod(idime,inodb) = coord(idime,ipoin)
                do jdime = 1,ndime
                   itens = ltens(idime,jdime)
                   streb(idime,jdime) = streb(idime,jdime) &
                                      + stres(itens,ipoin) &
                                      * elmar(pblty) % shape(inodb,igaub)
                end do
             end do 
          end do
          !
          ! Compute normal to surface at gauss point
          !
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                elcod(idime,inode) = coord(idime,ipoin)
             end do
          end do
          call bouder(pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),bocod,baloc,eucta)
          call chenor(pnode,baloc,bocod,elcod)
          !
          ! Calculate tractions ( t = sigma . n ) at gauss points
          !
          tractb = 0.0_rp
          do idime = 1,ndime
             do jdime = 1,ndime
                tractb(idime) = tractb(idime) + streb(idime,jdime) * baloc(jdime,ndime)
             end do
          end do
          if (kfl_how_forces == 1_ip) then
             !
             ! Spread tractions to nodes 
             !
             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                do idime = 1,ndime
                   if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
                      tract(idime,ipoin) = tract(idime,ipoin) &
                                         + tractb(idime) * elmar(pblty) % shape(inodb,igaub)
                   end if
                end do
             end do
          else if (kfl_how_forces == 2_ip) then
             !
             ! Average tractions on boundary element
             !
             do idime = 1,ndime
                if( kfl_fixbo_sld(iboun) == 6 ) then
                   tract(idime,iboun) = tract(idime,iboun) + tractb(idime) !/pgaub
                end if
             end do
          end if

       end do gauss_points

    end do boundaries

  end subroutine sld_integrate_tractions


end module mod_sld_efect
