!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup nastin
!> @{
!> @file    mod_nsi_algebraic_forces.f90
!> @author  houzeaux
!> @date    2020-03-19
!> @brief   Algebraic forces
!> @details Toolbox for algebaric forces
!-----------------------------------------------------------------------

module mod_nsi_algebraic_forces

  use def_kintyp
  use def_domain
  use def_elmtyp
  use def_master
  use def_nastin
  use mod_memory
  use mod_nsi_multi_step_fs, only : Grad_nsi => Grad 
  implicit none

  private

  public :: nsi_algebraic_forces_allocate
  public :: nsi_algebraic_forces

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-19
  !> @brief   Allocate memory
  !> @details ???
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_algebraic_forces_allocate()

    integer(ip) :: ipoin

    if( kfl_intfo_nsi > 0 ) then
       if( associated(intfo_nsi) ) then
          do ipoin = 1,npoin
             call memory_deallo(mem_modul(1:2,modul),'INTFO_NSI % AUU','nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % Auu)
             call memory_deallo(mem_modul(1:2,modul),'INTFO_NSI % AUP','nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % Aup)
             call memory_deallo(mem_modul(1:2,modul),'INTFO_NSI % BU' ,'nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % bu )
             intfo_nsi(ipoin) % kfl_exist = 0
          end do
          deallocate(intfo_nsi)
       end if
       allocate(intfo_nsi(npoin))
       do ipoin = 1,npoin
          nullify(intfo_nsi(ipoin) % Auu)
          nullify(intfo_nsi(ipoin) % Aup)
          nullify(intfo_nsi(ipoin) % bu)
          intfo_nsi(ipoin) % kfl_exist = 0
       end do
    end if

  end subroutine nsi_algebraic_forces_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-19
  !> @brief   Compute internal forces
  !> @details Compute algebaric forces
  !>          F = bu - Auu u - Aup p 
  !>          F is the force of the fluid on the solid
  !>          F has unit [RHO*U*^2L]
  !>          sig.n = F / mass
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_algebraic_forces(Auu,Aup,bu)

    real(rp),    intent(in) :: Auu(ndime,ndime,*) !> Auu
    real(rp),    intent(in) :: Aup(ndime,*)       !> Aup
    real(rp),    intent(in) :: bu(ndime,*)        !> bu
    integer(ip)             :: ipoin,jzdom,idime,jdime,izdom

    if( kfl_intfo_nsi >= 1 ) then

       !----------------------------------------------------------------
       !
       ! Determine where force should be computed
       !
       !----------------------------------------------------------------

       do ipoin = 1,npoin
          if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then
             if( kfl_intfo_nsi == 2 ) then
                call memory_deallo(mem_modul(1:2,modul),'INTFO_NSI % AUU','nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % Auu)
                call memory_deallo(mem_modul(1:2,modul),'INTFO_NSI % AUP','nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % Aup)
                call memory_deallo(mem_modul(1:2,modul),'INTFO_NSI % BU' ,'nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % bu )
             end if
          end if
       end do
       !
       ! On boundary nodes
       !
       do ipoin = 1,npoin
          if( ( lpoty(ipoin) > 0 ) .or. any( kfl_fixno_nsi(1:ndime,ipoin) /= 0 ) ) then
             intfo_nsi(ipoin) % kfl_exist = 1
          else
             intfo_nsi(ipoin) % kfl_exist = 0
          end if
       end do
       !
       ! Allocate
       !
       do ipoin = 1,npoin           
          if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then             
             jzdom = r_dom(ipoin+1)-r_dom(ipoin)
             call memory_alloca(mem_modul(1:2,modul),'INTFO_NSI % AUU','nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % Auu,ndime,ndime,jzdom)
             call memory_alloca(mem_modul(1:2,modul),'INTFO_NSI % AUP','nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % Aup,ndime,jzdom)
             call memory_alloca(mem_modul(1:2,modul),'INTFO_NSI % BU' ,'nsi_algebraic_forces_allocate',intfo_nsi(ipoin) % bu ,ndime )
          end if
       end do
       
       !----------------------------------------------------------------
       !
       ! Save Auu, Aup and bu in INTFO_NSI(:) % 
       !        
       !----------------------------------------------------------------

       if(      NSI_FRACTIONAL_STEP .and. kfl_grad_div_nsi == 0 ) then
          !
          ! Only Aup has been assembled
          !
          do ipoin = 1,npoin                      
             if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then
                do idime = 1,ndime
                   intfo_nsi(ipoin) % bu(idime) = bu(idime,ipoin)
                end do
                jzdom = 0
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jzdom = jzdom + 1
                   do idime = 1,ndime
                      intfo_nsi(ipoin) % Aup(idime,jzdom) = Aup(idime,izdom)
                   end do
                end do
             end if
          end do

       else if( NSI_FRACTIONAL_STEP .and. kfl_grad_div_nsi /= 0 ) then
          !
          ! Momentum matrices are not assembled
          !
          do ipoin = 1,npoin                      
             if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then
                do idime = 1,ndime
                   intfo_nsi(ipoin) % bu(idime) = bu(idime,ipoin)
                end do
                jzdom = 0
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jzdom = jzdom + 1
                   do idime = 1,ndime
                      intfo_nsi(ipoin) % Aup(idime,jzdom) = Grad_nsi(idime,izdom)
                   end do
                end do
             end if
          end do
       else
          !
          ! Matrices have been assembled
          !
          do ipoin = 1,npoin           

             if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then 
                jzdom = 0
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jzdom = jzdom + 1
                   do idime = 1,ndime
                      do jdime = 1,ndime
                         intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) = Auu(jdime,idime,izdom)
                      end do
                   end do
                   do idime = 1,ndime
                      intfo_nsi(ipoin) % Aup(idime,jzdom) = Aup(idime,izdom)
                   end do
                end do
                do idime = 1,ndime
                   intfo_nsi(ipoin) % bu(idime) = bu(idime,ipoin)
                end do
             end if
          end do
       end if

       kfl_intfo_nsi = 2

    end if

  end subroutine nsi_algebraic_forces

end module mod_nsi_algebraic_forces
!> @}
