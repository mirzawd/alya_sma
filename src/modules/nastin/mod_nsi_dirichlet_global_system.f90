!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_dirichlet_global_system
!> @author  Guillaume Houzeaux
!> @date    10/10/2016
!> @brief   Dirichlet conditions on Navier-Stokes
!> @details Impose Dirichlet conditions in global system
!>          For fractional, Auu and App are not assembled.
!>
!------------------------------------------------------------------------

module mod_nsi_dirichlet_global_system

  use def_kintyp
  use def_domain
  use def_master
  use def_nastin
  use mod_matrix,      only : matrix_rotate_system
  use mod_local_basis, only : local_basis_matrix
  use mod_nsi_schur_operations

  implicit none
  
  private

  public :: nsi_dirichlet_global_system  
  public :: nsi_dirichlet_matrix_rhs
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Impose Dirichlet on global system
  !> @details Impose Dirichlet boundary conditions on all the matrices
  !>          of the global system just after assembly
  !> 
  !-----------------------------------------------------------------------

   subroutine nsi_dirichlet_global_system()

     if( INOTEMPTY .and. kfl_matdi_nsi == NSI_DIRICHLET_MATRIX ) then
       if(      NSI_FRACTIONAL_STEP .and. kfl_grad_div_nsi == 0 ) then
          !
          ! Prescibe Laplacian and rotate matrices Aup and Apu
          !
          call nsi_dirichlet_matrix_rhs(&
               Q         = lapla_nsi,&
               ROTATION  = .false.,&
               DIRICHLET = .true.)
          call nsi_dirichlet_matrix_rhs(&
               Aup       = amatr(poaup_nsi:),&
               Apu       = amatr(poapu_nsi:),&
               bu        = rhsid,&
               uu        = unkno,&
               ROTATION  = .true.,&
               DIRICHLET = .false.)

       else if( NSI_FRACTIONAL_STEP .and. kfl_grad_div_nsi /= 0 ) then
          !
          ! Only bu
          !
          call nsi_dirichlet_matrix_rhs(&
               bu        = rhsid,&
               ROTATION  = .true.,&
               DIRICHLET = .false.)

          if(kfl_fsgrb_nsi/=0) call nsi_dirichlet_matrix_rhs(&
               bu        = rhsid_gravb,&
               ROTATION  = .true.,&
               DIRICHLET = .true.)

       else if( NSI_SCHUR_COMPLEMENT ) then
          !
          ! All matrices
          !
          call nsi_dirichlet_matrix_rhs(&
               Aup       = amatr(poaup_nsi:),&
               Apu       = amatr(poapu_nsi:),&
               Q         = lapla_nsi,&
               bu        = rhsid,&
               bp        = rhsid(ndbgs_nsi+1:),&
               uu        = unkno,&
               pp        = unkno(ndbgs_nsi+1:),&
               Auu       = amatr(poauu_nsi:),&
               App       = amatr(poapp_nsi:),&
               ROTATION  = .true.,&
               DIRICHLET = .true.)

       end if
    end if

  end subroutine nsi_dirichlet_global_system

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Impose Dirichlet on system
  !> @details Impose Dirichlet boundary conditions matrices and rhs
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_dirichlet_matrix_rhs(Aup,Apu,Q,bu,bp,uu,pp,Auu,App,ROTATION,DIRICHLET)

    real(rp),    intent(inout), optional :: Aup(ndime,nzdom)
    real(rp),    intent(inout), optional :: Apu(ndime,nzdom)
    real(rp),    intent(inout), optional :: Q(nzdom)
    real(rp),    intent(inout), optional :: bu(ndime,npoin)
    real(rp),    intent(inout), optional :: bp(npoin)
    real(rp),    intent(inout), optional :: uu(ndime,npoin)
    real(rp),    intent(inout), optional :: pp(npoin)
    real(rp),    intent(inout), optional :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(inout), optional :: App(nzdom)
    logical(lg), intent(in),    optional :: ROTATION
    logical(lg), intent(in),    optional :: DIRICHLET
    real(rp)                             :: Auud(3),Qd,Appd
    integer(ip)                          :: ipoin,jzdom,idime,jdime,izdom,jpoin
    integer(ip)                          :: izdod,ibopo,jbopo,kpoin,iroty,kdime
    integer(ip)                          :: idofn
    real(rp)                             :: worma(ndime,ndime)
    logical(lg)                          :: if_rotation
    logical(lg)                          :: if_dirichlet
    real(rp)                             :: rotma(ndime,ndime)

    if_rotation  = optional_argument(.true.,ROTATION)
    if_dirichlet = optional_argument(.true.,DIRICHLET)

    !----------------------------------------------------------------------
    !
    ! Rotate matrices of NS system Auu, Aup, Apu and bu
    !
    !----------------------------------------------------------------------

    if( kfl_local_nsi == 1 .and. if_rotation ) then

       do ipoin = 1,npoin
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             iroty =  kfl_fixrs_nsi(ipoin)
             
             if( iroty /= 0 ) then
 
                call local_basis_matrix(ipoin,ibopo,iroty,rotma)
                !
                ! Modifies column number IPOIN of AMATR ( A_j,imodi <-- A_j,imodi R )
                !
                do jzdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jpoin = c_dom(jzdom)
                   do izdom = r_dom(jpoin),r_dom(jpoin+1)-1
                      kpoin = c_dom(izdom)
                      if( kpoin == ipoin ) then
                         if( present(Auu) ) then
                            !
                            ! Auu
                            !
                            do idime = 1,ndime
                               do jdime = 1,ndime
                                  worma(idime,jdime) = 0.0_rp
                                  do kdime = 1,ndime
                                     worma(idime,jdime) = worma(idime,jdime) &
                                          + Auu(kdime,idime,izdom) * rotma(kdime,jdime)
                                  end do
                               end do
                            end do
                            do idime = 1,ndime
                               do jdime = 1,ndime
                                  Auu(jdime,idime,izdom) = worma(idime,jdime)
                               end do
                            end do
                         end if
                         if( present(Apu) ) then
                            !
                            ! Apu
                            !
                            do jdime = 1,ndime
                               worma(1,jdime) = 0.0_rp
                               do kdime = 1,ndime
                                  worma(1,jdime) = worma(1,jdime)&
                                       + Apu(kdime,izdom) * rotma(kdime,jdime)
                               end do
                            end do
                            do jdime = 1,ndime
                               Apu(jdime,izdom) = worma(1,jdime)
                            end do
                         end if
                      end if
                   end do
                end do

                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   !
                   ! Modifies row number IPOIN of AMATR ( A_imodi,j <-- R^t A_imodi,j )
                   !
                   if( present(Auu) ) then
                      !
                      ! Auu
                      !
                      jpoin = c_dom(izdom)
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            worma(idime,jdime) = 0.0_rp
                            do kdime = 1,ndime
                               worma(idime,jdime) = worma(idime,jdime) &
                                    + Auu(jdime,kdime,izdom) * rotma(kdime,idime)
                            end do
                         end do
                      end do
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            Auu(jdime,idime,izdom) = worma(idime,jdime)
                         end do
                      end do
                   end if
                   if( present(Aup) ) then
                      !
                      ! Aup
                      !
                      do idime = 1,ndime
                         worma(idime,1) = 0.0_rp
                         do kdime = 1,ndime
                            worma(idime,1) = worma(idime,1) &
                                 + rotma(kdime,idime) * Aup(kdime,izdom)
                         end do
                      end do
                      do idime = 1,ndime
                         Aup(idime,izdom) = worma(idime,1)
                      end do
                   end if

                end do
                if( present(bu) ) then
                   !
                   ! bu
                   !
                   do idime = 1,ndime
                      worma(idime,1) = 0.0_rp
                      do kdime = 1,ndime
                         worma(idime,1) = worma(idime,1) &
                              + rotma(kdime,idime) * bu(kdime,ipoin)
                      end do
                   end do
                   do idime = 1,ndime
                      bu(idime,ipoin) = worma(idime,1)
                   end do
                end if
             end if
          end if
       end do
    end if

    if( if_dirichlet ) then

       !----------------------------------------------------------------------
       !
       ! Impose velocity
       !
       !----------------------------------------------------------------------

       do ipoin = 1,npoin

          do idime = 1,ndime

             if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
                !
                ! Eliminate dof of IPOIN from other equations (JPOIN)
                !
                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                   jpoin = c_dom(izdom)
                   if( ipoin /= jpoin ) then

                      do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                         kpoin = c_dom(jzdom)
                         if( kpoin == ipoin ) then
                            if( present(Auu) .and. present(bu) ) then
                               !
                               ! bu <= bu - Auu * u
                               !
                               do jdime = 1,ndime
                                  bu(jdime,jpoin) = bu(jdime,jpoin) - Auu(idime,jdime,jzdom) * bvess_nsi(idime,ipoin,1)
                                  Auu(idime,jdime,jzdom) = 0.0_rp
                               end do
                            end if
                            if( present(Apu) .and. present(bp) ) then
                               !
                               ! bp <= bp - Apu * u
                               !
                               if( kfl_grad_div_nsi == 0 ) then
                                  bp(jpoin) = bp(jpoin) - Apu(idime,jzdom) * bvess_nsi(idime,ipoin,1)
                                  Apu(idime,jzdom) = 0.0_rp
                               end if
                            end if
                         end if
                      end do

                   end if
                end do
                !
                ! IZDOD: Diagonal
                !
                izdod = r_dom(ipoin) - 1
                jpoin = 0
                do while( jpoin /= ipoin )
                   izdod = izdod + 1
                   jpoin = c_dom(izdod)
                end do
                if( present(Auu) ) then
                   Auud(idime) = Auu(idime,idime,izdod)
                else
                   Auud(idime) = 1.0_rp
                end if
                if( abs(Auud(idime)) < zeror ) Auud(idime) = 1.0_rp
                !
                ! Set line to zero
                !
                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                   if( present(Auu) ) then
                      do jdime = 1,ndime
                         Auu(jdime,idime,izdom) = 0.0_rp
                      end do
                   end if
                   if( present(Aup) ) Aup(idime,izdom) = 0.0_rp
                end do
                !
                ! Prescribe value
                !
                idofn = (ipoin-1)*ndime + idime
                if( present(Auu) ) then
                   if( present(Auu) ) Auu(idime,idime,izdod) = Auud(idime)
                   if( present(bu)  ) bu(idime,ipoin)        = bvess_nsi(idime,ipoin,1) * Auud(idime)
                   if( present(uu)  ) uu(idime,ipoin)        = bvess_nsi(idime,ipoin,1)
                else
                   !
                   ! In the case when Auu is eliminated (bu <= bu - Auu*u), bu
                   ! should be zero on boundary condition
                   !
                   if( present(bu)  ) bu(idime,ipoin)        = 0.0_rp
                   if( present(uu)  ) uu(idime,ipoin)        = bvess_nsi(idime,ipoin,1)
                end if

             end if

          end do

       end do

       !----------------------------------------------------------------------
       !
       ! Impose pressure Schur complement preconditioner
       !
       !----------------------------------------------------------------------


       if( present(Q) ) then

          if( kfl_confi_nsi /= 2 ) then
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                   jpoin = c_dom(izdom)
                   if( ipoin /= jpoin ) then
                      jbopo = lpoty(jpoin)
                      if( jbopo /= 0 ) then
                         if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                            Q(izdom) = 0.0_rp
                         end if
                      end if
                   end if
                end do
             end do
          end if

          do ipoin = 1,npoin

             ibopo = lpoty(ipoin)
             if( ibopo > 0 ) then

                if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                   !
                   ! IZDOD: Diagonal
                   !
                   izdod = r_dom(ipoin) - 1
                   jpoin = 0
                   do while( jpoin /= ipoin )
                      izdod = izdod + 1
                      jpoin = c_dom(izdod)
                   end do
                   Qd = Q(izdod)
                   if( abs(Qd) < zeror ) Qd = 1.0_rp

                   if( kfl_confi_nsi == 2 ) then
                      !
                      ! Increase diagonal for weak imposition of pressure
                      !
                      Q(izdod) = mulpr_nsi*Qd
                      
                   else
                      !
                      ! Set line to zero and maintain diagonal
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         Q(izdom) = 0.0_rp
                      end do
                      Q(izdod) = Qd
                      
                   end if

                end if

             end if
          end do

       end if

       !----------------------------------------------------------------------
       !
       ! Impose pressure in Apu, App, Aup, bp, bu
       !
       !----------------------------------------------------------------------

       if( present(App) .or. present(bp) .or. present(pp) .or. ( present(Aup) .and. present(bu) ) ) then

          do ipoin = 1,npoin

             if( solve(1) % block_array(2) % kfl_fixno(1,ipoin) > 0 ) then

                if( solve(2) % kfl_symme == 1 ) then
                   call runend('NOT CODED: CHECK IT')
                else
                   !
                   ! Eliminate pressure at IPOIN on all lines
                   !
                   do jpoin = 1,npoin
                      if( ipoin /= jpoin ) then
                         izdom = r_dom(jpoin)
                         do while( izdom < r_dom(jpoin+1) )
                            kpoin = c_dom(izdom)
                            if( kpoin == ipoin ) then
                               if( present(App) .and. present(bp) ) then
                                  bp(jpoin)  = bp(jpoin) - App(izdom) * solve(1) % block_array(2) % bvess(1,ipoin)
                                  App(izdom) = 0.0_rp
                               end if
                               izdom      = r_dom(jpoin+1)
                            end if
                            izdom = izdom + 1
                         end do
                      end if
                   end do
                   !
                   ! Diagonal
                   !
                   izdod = r_dom(ipoin) - 1
                   jpoin = 0
                   do while( jpoin /= ipoin )
                      izdod = izdod + 1
                      jpoin = c_dom(izdod)
                   end do
                   if( present(App) ) then
                      Appd = App(izdod)
                   else
                      Appd = 0.0_rp
                   end if
                   if( abs(Appd) < zeror ) Appd = 1.0_rp
                   !
                   ! Set line to zero
                   !
                   do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                      if( present(App) ) App(izdom) = 0.0_rp
                      if( present(Apu) ) then
                         do idime = 1,ndime
                            Apu(idime,izdom) = 0.0_rp
                         end do
                      end if
                   end do
                   !
                   ! Presrcibe value
                   !
                   if( present(App) ) App(izdod) = Appd
                   if( present(bp)  ) bp(ipoin)  = Appd * solve(1) % block_array(2) % bvess(1,ipoin)
                   if( present(pp)  ) pp(ipoin)  = solve(1) % block_array(2) % bvess(1,ipoin)

                end if
                !
                ! Eliminate prescribed pressure from IPOIN momentum equation
                !
                if( present(Aup) .and. present(bu) ) then
                   do jpoin = 1,npoin
                      izdom = r_dom(jpoin)
                      do while( izdom < r_dom(jpoin+1) )
                         kpoin = c_dom(izdom)
                         if( kpoin == ipoin ) then
                            do idime = 1,ndime
                               bu(idime,jpoin)  = bu(idime,jpoin) - Aup(idime,izdom) * solve(1) % block_array(2) % bvess(1,ipoin)
                               Aup(idime,izdom) = 0.0_rp
                            end do
                            izdom = r_dom(jpoin+1)
                         end if
                         izdom = izdom + 1
                      end do
                   end do
                end if

             end if

          end do

       end if

    end if

  end subroutine nsi_dirichlet_matrix_rhs

end module mod_nsi_dirichlet_global_system
!> @}
