!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_tubes.f90
!> @author  houzeaux
!> @date    2020-09-16
!> @brief   Coupling with rubes
!> @details Coupling of alya with tubes
!-----------------------------------------------------------------------

module mod_nsi_tubes

  use def_master
  use def_domain
  use def_nastin
  use mod_tubes,                 only : tubes_number
  use mod_tubes,                 only : tubes_total_number
  use mod_tubes,                 only : tubes_recv_value_from_module
  use mod_tubes,                 only : tubes_give_value_to_module
  use mod_communications_global, only : PAR_SUM
  use mod_bouder
  implicit none
  private

  public :: nsi_nastin_to_tubes
  public :: nsi_tubes_to_nastin

contains
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Tubes to Nastin
  !> @details Get pressure from tubes
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_tubes_to_nastin()

    integer(ip)             :: num_tubes,itube
    integer(ip)             :: iboun,icode
    real(rp),   allocatable :: tubes_pressure(:)

    num_tubes = tubes_total_number()

    if( num_tubes > 0 ) then
       allocate(tubes_pressure(num_tubes))
       do iboun = 1,nboun
          icode = kfl_codbo(iboun)
          itube = tubes_number(icode)
          if( itube > 0 ) then
             bvnat_nsi(1,iboun,1) = tubes_give_value_to_module(icode)
          end if
       end do
    end if
    
  end subroutine nsi_tubes_to_nastin
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Nastin to tubes
  !> @details Compute input parameters for tubes (mass)
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_nastin_to_tubes()

    integer(ip)             :: num_tubes,itube
    integer(ip)             :: iboun,pblty,pnodb
    integer(ip)             :: pgaub,ipoin,inodb
    integer(ip)             :: igaub,idime
    real(rp)                :: bocod(ndime,mnodb)
    real(rp)                :: bovel(ndime,mnodb)
    real(rp)                :: gbvel(ndime)
    real(rp)                :: baloc(ndime,ndime)
    real(rp)                :: gbsur,norma
    real(rp),   allocatable :: tubes_mass(:)

    num_tubes = tubes_total_number()
    
    if( num_tubes > 0 ) then
       allocate(tubes_mass(num_tubes))
       tubes_mass(:) = 0.0_rp

       boundaries: do iboun = 1,nboun

          itube = tubes_number(kfl_codbo(iboun))
          if( itube > 0 ) then
             pblty = ltypb(iboun) 
             pnodb = nnode(pblty)
             pgaub = ngaus(pblty)
             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
                bovel(1:ndime,inodb) = veloc(1:ndime,ipoin,1)
             end do
             do igaub = 1,pgaub
                gbvel = 0.0_rp
                norma = 0.0_rp
                do inodb = 1,pnodb
                   gbvel(1:ndime) = gbvel(1:ndime) + elmar(pblty)%shape(inodb,igaub) * bovel(1:ndime,inodb)
                end do
                call bouder(&
                     pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),& 
                     bocod,baloc,gbsur)
                do idime =1,ndime 
                  norma = norma + sqrt(baloc(idime,ndime)*baloc(idime,ndime))
                end do
                gbsur = elmar(pblty)%weigp(igaub)*gbsur
                tubes_mass(itube) = tubes_mass(itube) + dot_product(gbvel(1:ndime),baloc(1:ndime,ndime)/norma)*gbsur
             end do
          end if

       end do boundaries
 
       call PAR_SUM(num_tubes,tubes_mass)
       !
       ! Pass values to tubes
       !
       do itube = 1,num_tubes
          call tubes_recv_value_from_module(tubes_mass)
       end do
       deallocate(tubes_mass)

    end if

    
  end subroutine nsi_nastin_to_tubes

end module mod_nsi_tubes
!> @}
