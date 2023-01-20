!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_nsi_immersed_boundary

  !------------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    mod_nsi_immersed_boundary.f90
  !> @date    24/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Immersed boundary
  !> @details Immersed boundary
  !>
  !> @{
  !------------------------------------------------------------------------

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : INOTMASTER
  use def_master,         only : veloc,mem_modul
  use def_master,         only : modul,zeror
  use def_master,         only : veloc
  use def_domain,         only : npoin,ndime
  use def_domain,         only : ntens
  use def_nastin,         only : kfl_immer_nsi
  use def_nastin,         only : lagra_nsi
  use def_nastin,         only : tauib_nsi
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_ker_proper,     only : ker_proper
  use mod_communications, only : PAR_SUM
  use mod_gradie,         only : graten
  implicit none

  public :: nsi_ib_lagrange_multiplier_update
  
contains 

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/04/2017
  !> @brief   Update Lagrange multiplier
  !> @details Immersed boundary Lagrange multiplier update
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_ib_lagrange_multiplier_update()

    integer(ip)         :: ipoin,dummi,idime,jdime
    real(rp),   pointer :: mixin(:)
    real(rp),   pointer :: defor(:,:)
    real(rp)            :: resid(2,ndime),new_value(ndime),denom,suma(ndime)
    real(rp)            :: tauib_new(ndime,ndime),rr=1.0_rp

    nullify(mixin)
    nullify(defor)
          
    if( kfl_immer_nsi == 1 ) then
       !
       ! Velocity based
       !
       if( INOTMASTER ) then
          
          call memory_alloca(mem_modul(1:2,modul),'MIXIN','nsi_ib_lagrange_multiplier_update',mixin,npoin)
          call ker_proper('MIXIN','NPOIN',dummi,dummi,mixin)

          resid = 0.0_rp
          suma  = 0.0_rp

          do ipoin = 1,npoin
             !new_value            = lagra_nsi(:,ipoin,1) + 30.0_rp * rr * mixin(ipoin) * veloc(:,ipoin,1)
             new_value            = lagra_nsi(:,ipoin,1) + 30.0_rp * rr * veloc(:,ipoin,1)
             do idime = 1,ndime
                resid(1,idime) = resid(1,idime) + (new_value(idime)-lagra_nsi(idime,ipoin,1))**2
                resid(2,idime) = resid(2,idime) + new_value(idime)**2
             end do
             lagra_nsi(:,ipoin,1) = new_value
             suma                 = suma + new_value
          end do
          suma = suma / real(npoin,rp)
          !do idime = 1,ndime
          !   suma(idime) = maxval(lagra_nsi(idime,:,1))
          !end do

          !suma =  lagra_nsi(:,1,1)
          do ipoin = 1,npoin
             !lagra_nsi(:,ipoin,1) = lagra_nsi(:,ipoin,1) - suma !- lagra_nsi(:,1,1)
             !lagra_nsi(:,ipoin,1) = lagra_nsi(:,ipoin,1) - 20.0_rp
          end do
          
          call memory_deallo(mem_modul(1:2,modul),'MIXIN','nsi_ib_lagrange_multiplier_update',mixin)
       end if
       
       !call PAR_SUM(resid)

       idime = 1
       denom = sqrt(resid(2,idime))
       resid(1,idime) = sqrt( resid(1,idime) / max(zeror,resid(2,idime)) )
       write(90,*) resid(1,idime),denom,maxval(abs(lagra_nsi(1,:,1))),maxval(abs(lagra_nsi(2,:,1))) ; flush(90)

    else if( kfl_immer_nsi == 2 ) then
       !
       ! Deformation based
       !
       if( INOTMASTER ) then

          call memory_alloca(mem_modul(1:2,modul),'DEFOR','nsi_ib_lagrange_multiplier_update',defor,ntens,npoin)
          call graten(veloc,defor)

          resid = 0.0_rp
          
          if( ndime == 2 ) then
             do ipoin = 1,npoin                
                tauib_new(1,1) = 0.5_rp * defor(1,ipoin)
                tauib_new(2,2) = 0.5_rp * defor(2,ipoin)
                tauib_new(1,2) = 0.5_rp * defor(3,ipoin)
                tauib_new(2,1) = 0.5_rp * defor(3,ipoin)
                do jdime = 1,ndime
                   do idime = 1,ndime
                      resid(1,1)                     = resid(1,1) + (tauib_new(idime,jdime)-tauib_nsi(idime,jdime,ipoin))**2
                      resid(2,1)                     = resid(2,1) + (tauib_new(idime,jdime))**2
                      tauib_nsi(idime,jdime,ipoin) = tauib_nsi(idime,jdime,ipoin) + 2.0_rp * 10.0_rp * rr * tauib_new(idime,jdime)
                   end do
                end do
             end do
          else
             do ipoin = 1,npoin               
                tauib_new(1,1) = 0.5_rp * defor(1,ipoin)
                tauib_new(2,2) = 0.5_rp * defor(2,ipoin)
                tauib_new(1,2) = 0.5_rp * defor(3,ipoin)
                tauib_new(2,1) = 0.5_rp * defor(3,ipoin)             
                tauib_new(3,3) = 0.5_rp * defor(4,ipoin)             
                tauib_new(1,3) = 0.5_rp * defor(5,ipoin)             
                tauib_new(3,1) = 0.5_rp * defor(5,ipoin)             
                tauib_new(2,3) = 0.5_rp * defor(6,ipoin)             
                tauib_new(3,2) = 0.5_rp * defor(6,ipoin)
                do jdime = 1,ndime
                   do idime = 1,ndime
                      resid(1,1)                     = resid(1,1) + (tauib_new(idime,jdime)-tauib_nsi(idime,jdime,ipoin))**2
                      resid(2,1)                     = resid(2,1) + (tauib_new(idime,jdime))**2
                      tauib_nsi(idime,jdime,ipoin) = tauib_nsi(idime,jdime,ipoin) + 2.0_rp * 10.0_rp * rr * tauib_new(idime,jdime)
                   end do
                end do
              end do
           end if
           
          call memory_deallo(mem_modul(1:2,modul),'DEFOR','nsi_ib_lagrange_multiplier_update',defor)

       end if
       
       !call PAR_SUM(resid)

       denom = sqrt(resid(2,1))
       resid(1,1) = sqrt( resid(1,1) / max(zeror,resid(2,1)) )
       write(90,*) resid(1,1),denom ; flush(90)
       
    end if

  end subroutine nsi_ib_lagrange_multiplier_update

end module mod_nsi_immersed_boundary
!> @}
