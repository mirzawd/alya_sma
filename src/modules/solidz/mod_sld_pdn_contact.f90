!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_pdn_contact.f90
!> @author  Matias Rivero
!> @date    October, 2018
!>          - Subroutine written
!> @brief   PDN contact unilateral and bilateral
!> @details
!> @author  Miguel Zavala and Gerard Guillamet
!> @date    February, 2019
!>          - RBO subroutines written
!> @brief   PDN contact for solving RBO + deformable problem.
!> @details
!>
!>          Generalized version of PDN contact algorithm:
!>          It is compatible to solve the UNILATERAL, BILATERAL AND RBO_DEFORMABLE problems
!>          Refactoring of the following subroutines:
!>            - plugin
!>            - projection send
!>            - projection recieve
!>
!> @}
!-----------------------------------------------------------------------

module mod_sld_pdn_contact
  
  use def_kintyp,            only : ip, rp, lg
  use def_kintyp_solvers,    only : soltyp
  use def_master,            only : IMASTER, INOTMASTER, npoi3, zeror
  use def_master,            only : momod, modul, intost
  use def_domain,            only : npoin, nboun, ndime, coord
  use mod_communications,    only : PAR_SUM, PAR_MAX
  use mod_messages,          only : messages_live
#if defined COMMDOM && COMMDOM == 2
  use mod_plepp_pdn_contact, only : CNT_CPLNG
  use mod_plepp_pdn_contact, only : PLEPP_CPLNG
  use mod_plepp_pdn_contact, only : CNT_SENDRECV
  use mod_plepp_pdn_contact, only : plepp_pdn_driver_exchange02
  use mod_plepp_pdn_contact, only : plepp_pdn_dynamic_reduce_sum
#endif
  use def_solidz,            only : kfl_fixno_sld
  use def_solidz,            only : release_nodes, new_to_release
  use def_solidz,            only : frxid_sld, bvess_sld, bvnat_sld
  use def_solidz,            only : kfl_rigid_sld
  use def_solidz,            only : kfl_conta_sld,SLD_PDN_BILATERAL
  use mod_sld_rbo,           only : sld_rbo_updunk
  use mod_sld_csys,          only : sld_csys_rotuni

  implicit none

  type(soltyp), pointer :: solve(:)

#if defined COMMDOM && COMMDOM == 2
 
  public  :: commdom_sld_plugin_pdn        ! Original PDN contact
  public  :: sld_rbo_contact_plugin_rbo    ! New Generalized version of PDN contact
  public  :: commdom_sld_nodes_release

  private :: commdom_proje_nodes_2d
  private :: commdom_proje_nodes_3d
  private :: sld_pdn_conta_projection_send ! New projection send
  private :: sld_pdn_conta_projection_recv ! New projection recv
  private :: commdom_sld_calc_reaction
  private :: commdom_sld_release_nodes
  private :: instr2
  private :: commdom_sld_bilateral_neumann
  private :: sld_rbo_bilateral_send        ! New RBO bilateral send
  private :: sld_rbo_bilateral_recv        ! New RBO bilateral recv
!!!  private :: commdom_sld_impose_friction   ! Not used

#endif
  
contains

#if defined COMMDOM && COMMDOM == 2
  
  !-----------------------------------------------------------------------
  !>
  !> @author  Miguel Zavala and Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Main subroutine for PDN contact coupling (Generalized)
  !> @details Main subroutine for PDN contact coupling (Generalized)
  !>
  !-----------------------------------------------------------------------

  subroutine sld_rbo_contact_plugin_rbo()

    implicit none

    integer(ip)       :: nsend              !< Localized nodes sent
    integer(ip)       :: nrecv              !< Localized nodes recieved
    integer(ip)       :: ndata              !< Dimensions data
    real(rp), pointer :: data_send(:,:) => null()
    real(rp), pointer :: data_recv(:,:) => null()

    if ( any(CNT_SENDRECV) ) then
       !
       ! Data from PLEPP
       !
       nsend = PLEPP_CPLNG%n_ij
       nrecv = PLEPP_CPLNG%n_ji
       !
       if ( ndime == 2_ip ) then
          ndata = 8_ip
       else
          ndata = 11_ip
       end if
       !
       ! Allocate memory
       !
       allocate(data_send(ndata,nsend))
       allocate(data_recv(ndata,nrecv))
       !
       ! Initialize
       !
       data_send(:,:) = 0.0_rp
       data_recv(:,:) = 0.0_rp

       code_i: if(CNT_CPLNG%current_code == CNT_CPLNG%code_i) then

          !------------------------------------------------------------------
          !
          ! Rigid body (code_i==1)
          !
          !------------------------------------------------------| U--> |---!

          if( CNT_SENDRECV(7) ) call sld_pdn_conta_projection_send(ndata,nsend,nrecv,data_send(:,:),data_recv(:,:))
          !
          ! Bilateral
          !---------------------------------------------------| dUdn<-- |---!
          if( CNT_SENDRECV(12) ) call sld_rbo_bilateral_recv()

       end if code_i

       code_j: if(CNT_CPLNG%current_code == CNT_CPLNG%code_j) then

          !------------------------------------------------------------------
          !
          ! Deformable body (code_j==2)
          !
          !------------------------------------------------------| U<-- |---!

          if( CNT_SENDRECV(7) ) call sld_pdn_conta_projection_recv(ndata,nsend,nrecv,data_send(:,:),data_recv(:,:))
          !
          ! Bilateral
          !---------------------------------------------------| dUdn--> |---!
          if( CNT_SENDRECV(12) ) call sld_rbo_bilateral_send()

       end if code_j
       !
       ! Deallocate memory arrays
       !
       deallocate(data_send)
       deallocate(data_recv)

    endif

  end subroutine sld_rbo_contact_plugin_rbo

  !-----------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details Subroutine plugin PDN
  !-----------------------------------------------------------------------

  subroutine commdom_sld_plugin_pdn()

    use def_master,          only : displ, TIME_N, iblok, ITER_K
    use def_domain,          only : ndime, kfl_codno, lpoty
    use def_coupli,          only : coupling_driver_iteration, coupling_driver_max_iteration, kfl_gozon
    use mod_communications,  only : PAR_SUM
    use def_solidz,          only : jacrot_du_dq_sld, jacrot_dq_du_sld, contactbou_sld
    use def_solidz,          only : release_nodes, new_to_release, SLD_PDN_BILATERAL
    use def_solidz,          only : kfl_conta_sld, vect_proje_sld, petol_sld
    use def_solidz,          only : kfl_timet_sld, SLD_IMPLICIT_SCHEME
    use def_solidz,          only : kfl_conta_sld, SLD_PDN_UNILATERAL
    
    implicit none

    integer(ip)       :: ipoin, idime
    real(rp), pointer :: norm_proje(:,:) => null()
    real(rp), pointer :: norm_proje2(:,:) => null()
    real(rp), pointer :: tangent_proje(:,:) => null()
    real(rp), pointer :: norm_dist(:) => null()
    real(rp), pointer :: in_contact(:) => null()
    real(rp)          :: e_recta(ndime)
    integer(ip)       :: nsend, nrecv
    real(rp), pointer :: data_send(:,:) => null()
    real(rp), pointer :: data_recv(:,:) => null()
    real(rp)          :: tmp_tgx, tmp_tgy, tmp_tgz
    !
    integer(ip)       :: ibopo
    real(rp)          :: resid_tmp(ndime), dummy_matrix(ndime,ndime)
    real(rp)          :: prop_in, prop_out

    solve => momod(modul) % solve(1:)

    nsend = PLEPP_CPLNG%n_ij
    nrecv = PLEPP_CPLNG%n_ji

    !
    ! Allocate memory
    !
    if ( ndime == 2_ip ) then
       allocate(data_send(8,nsend))
       allocate(data_recv(8,nrecv))
    else if ( ndime == 3_ip ) then
       allocate(data_send(11,nsend))
       allocate(data_recv(11,nrecv))
    end if
    allocate(norm_proje(ndime,nsend))
    allocate(norm_proje2(ndime,nsend))
    allocate(tangent_proje(ndime,nsend))
    allocate(norm_dist(nsend))
    allocate(in_contact(nsend))
    !
    ! Initialize
    !
    norm_proje    = 0.0_rp
    norm_proje2   = 0.0_rp
    tangent_proje = 0.0_rp
    norm_dist     = 0.0_rp
    in_contact    = 666.0_rp
    data_send     = 0.0_rp
    data_recv     = 0.0_rp

    if ( any(CNT_SENDRECV) ) then

       if ( CNT_CPLNG%current_code == CNT_CPLNG%code_i ) then

          !------------------------------------------------------------------
          !
          !  Rigid body (Neumann)
          !
          !------------------------------------------------------------------

          if ( CNT_SENDRECV(7) ) then

             if ( ndime == 2_ip ) then  ! 2-d problems
                !
                ! Projection vector --> CHANGE THIS: MUST BE DYNAMICALLY DEFINED
                !
                e_recta(1) = vect_proje_sld(1)
                e_recta(2) = vect_proje_sld(2)
                !
                ! Projection of nodes detected inside indenter body to the boundary
                !
                call commdom_proje_nodes_2d(tangent_proje,norm_proje,norm_proje2,norm_dist,in_contact,e_recta)
                do ipoin = 1,nsend
                   data_send(1,ipoin) = tangent_proje(1,ipoin) ! tangent vector
                   data_send(2,ipoin) = tangent_proje(2,ipoin)
                   data_send(3,ipoin) = norm_proje(1,ipoin)    ! normal
                   data_send(4,ipoin) = norm_proje(2,ipoin)
                   data_send(5,ipoin) = norm_dist(ipoin)       ! distance
                   data_send(6,ipoin) = norm_proje2(1,ipoin)
                   data_send(7,ipoin) = norm_proje2(2,ipoin)
                   data_send(8,ipoin) = in_contact(ipoin)
                end do
                call commdom_locator_exchange_double_stride(data_send(1:8,1:nsend),data_recv(1:8,1:nrecv),8)

             else if ( ndime == 3_ip ) then  ! 3-d problems
                !
                ! Projection vector --> CHANGE THIS: MUST BE DYNAMICALLY DEFINED
                !
                e_recta(1) = vect_proje_sld(1)
                e_recta(2) = vect_proje_sld(2)
                e_recta(3) = vect_proje_sld(3)
                !
                ! Projection of nodes detected inside indenter body to the boundary
                !
                call commdom_proje_nodes_3d(tangent_proje,norm_proje,norm_proje2,norm_dist,in_contact,e_recta)
                do ipoin = 1,nsend
                   data_send(1, ipoin) = tangent_proje(1,ipoin)
                   data_send(2, ipoin) = tangent_proje(2,ipoin)
                   data_send(3, ipoin) = tangent_proje(3,ipoin)
                   data_send(4, ipoin) = norm_proje(1,ipoin)
                   data_send(5, ipoin) = norm_proje(2,ipoin)
                   data_send(6, ipoin) = norm_proje(3,ipoin)
                   data_send(7, ipoin) = norm_dist(ipoin)
                   data_send(8, ipoin) = norm_proje2(1,ipoin)
                   data_send(9, ipoin) = norm_proje2(2,ipoin)
                   data_send(10,ipoin) = norm_proje2(3,ipoin)
                   data_send(11,ipoin) = in_contact(ipoin)
                end do
                call commdom_locator_exchange_double_stride(data_send(1:11,1:nsend),data_recv(1:11,1:nrecv),11)

             end if

          end if

          if ( CNT_SENDRECV(12) .and. kfl_conta_sld == SLD_PDN_BILATERAL ) then
             !
             ! BILATERAL
             !
             call commdom_sld_bilateral_neumann()

          end if

       endif

       if ( CNT_CPLNG%current_code == CNT_CPLNG%code_j ) then

          !------------------------------------------------------------------
          !
          !  Deformable body (Dirichlet)
          !
          !------------------------------------------------------------------

          if ( CNT_SENDRECV(7) ) then

             if ( INOTMASTER ) then

                do ipoin = 1,npoin
                   do idime = 1,ndime
                      if ( kfl_fixno_sld(idime,ipoin) == 3_ip ) kfl_fixno_sld(idime,ipoin) = 0_ip
                   end do
                end do

                if ( ndime == 2_ip ) then ! 2-d problem

                   call commdom_locator_exchange_double_stride(data_send(1:8,1:nsend),data_recv(1:8,1:nrecv),8)
                   !TODO: debo eliminar los nodos que no son de la boundary de contacto del deformable
                   !para eso debo ver la interior list
                   do ipoin = 1,nrecv
                      if( any(kfl_codno(:,PLEPP_CPLNG%interior_list_j(ipoin)) == contactbou_sld) ) then
                       

                         jacrot_du_dq_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)) =  data_recv(3,ipoin) !nx
                         jacrot_du_dq_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin)) =  data_recv(4,ipoin) !ny
                         jacrot_du_dq_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin)) = -data_recv(4,ipoin) !tx
                         jacrot_du_dq_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin)) =  data_recv(3,ipoin) !ty

                         jacrot_dq_du_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin))

                         !
                         ! Contact nodes with fixity in the first dof
                         !
                         if( kfl_fixno_sld(1,ipoin) == 1_ip ) then
                            do idime = 1,ndime-1
                               kfl_fixno_sld(idime+1,ipoin) = 1_ip
                            end do
                         end if
                         
                         kfl_fixno_sld(1,PLEPP_CPLNG%interior_list_j(ipoin)) = 3_ip
                         
                         if( kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then
                            displ(1,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N)     = data_recv(5,ipoin)*data_recv(6,ipoin) + &
                                 displ(1,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                            displ(2,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N)     = data_recv(5,ipoin)*data_recv(7,ipoin) + &
                                 displ(2,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                         else
                            bvess_sld(1,PLEPP_CPLNG%interior_list_j(ipoin),ITER_K) = data_recv(5,ipoin)*data_recv(6,ipoin) + &
                                 displ(1,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                            bvess_sld(2,PLEPP_CPLNG%interior_list_j(ipoin),ITER_K) = data_recv(5,ipoin)*data_recv(7,ipoin) + &
                                 displ(2,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                         end if
                         !kfl_fixno_sld(2,1) = 3_ip   !<-bilateral-cylinder
                         !displ(1,1,TIME_N)  = 0.0_rp

                      end if
                   end do

                else if ( ndime == 3_ip ) then !---------> NDIME = 3

                   call commdom_locator_exchange_double_stride(data_send(1:11,1:nsend),data_recv(1:11,1:nrecv),11)
                   !TODO: debo eliminar los nodos que no son de la boundary de contacto del deformable
                   !para eso debo ver la interior list
                   do ipoin = 1,nrecv
                      if( any(kfl_codno(:,PLEPP_CPLNG%interior_list_j(ipoin)) == contactbou_sld) ) then

                         jacrot_du_dq_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)) = data_recv(4,ipoin) !nx
                         jacrot_du_dq_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin)) = data_recv(5,ipoin) !ny
                         jacrot_du_dq_sld(3,1,PLEPP_CPLNG%interior_list_j(ipoin)) = data_recv(6,ipoin) !nz
                         jacrot_du_dq_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin)) = data_recv(1,ipoin) !t1x
                         jacrot_du_dq_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin)) = data_recv(2,ipoin) !t1y
                         jacrot_du_dq_sld(3,2,PLEPP_CPLNG%interior_list_j(ipoin)) = data_recv(3,ipoin) !t1z

                         tmp_tgx =   data_recv(5,ipoin)*data_recv(3,ipoin) - data_recv(6,ipoin)*data_recv(2,ipoin)
                         tmp_tgy = -(data_recv(4,ipoin)*data_recv(3,ipoin) - data_recv(6,ipoin)*data_recv(1,ipoin))
                         tmp_tgz =   data_recv(4,ipoin)*data_recv(2,ipoin) - data_recv(5,ipoin)*data_recv(1,ipoin)

                         jacrot_du_dq_sld(1,3,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              tmp_tgx/sqrt(tmp_tgx*tmp_tgx + tmp_tgy*tmp_tgy + tmp_tgz*tmp_tgz) !t2x
                         jacrot_du_dq_sld(2,3,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              tmp_tgy/sqrt(tmp_tgx*tmp_tgx + tmp_tgy*tmp_tgy + tmp_tgz*tmp_tgz) !t2y
                         jacrot_du_dq_sld(3,3,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              tmp_tgz/sqrt(tmp_tgx*tmp_tgx + tmp_tgy*tmp_tgy + tmp_tgz*tmp_tgz) !t2z

                         jacrot_dq_du_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(1,3,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(3,1,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(2,3,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(3,2,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(3,1,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(1,3,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(3,2,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(2,3,PLEPP_CPLNG%interior_list_j(ipoin))
                         jacrot_dq_du_sld(3,3,PLEPP_CPLNG%interior_list_j(ipoin)) = &
                              jacrot_du_dq_sld(3,3,PLEPP_CPLNG%interior_list_j(ipoin))

                         !
                         ! Contact nodes with fixity in the first dof
                         !
                         if( kfl_fixno_sld(1,ipoin) == 1_ip ) then
                            do idime = 1,ndime-1
                               kfl_fixno_sld(idime+1,ipoin) = 1_ip
                            end do
                         end if

                         kfl_fixno_sld(1,PLEPP_CPLNG%interior_list_j(ipoin)) = 3_ip

                         if( kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then
                            displ(1,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N)     = data_recv(7,ipoin)*data_recv(8,ipoin)  + &
                                 displ(1,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                            displ(2,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N)     = data_recv(7,ipoin)*data_recv(9,ipoin)  + &
                                 displ(2,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                            displ(3,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N)     = data_recv(7,ipoin)*data_recv(10,ipoin) + &
                                 displ(3,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                         else
                            bvess_sld(1,PLEPP_CPLNG%interior_list_j(ipoin),ITER_K) = data_recv(7,ipoin)*data_recv(8,ipoin)  + &
                                 displ(1,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                            bvess_sld(2,PLEPP_CPLNG%interior_list_j(ipoin),ITER_K) = data_recv(7,ipoin)*data_recv(9,ipoin)  + &
                                 displ(2,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                            bvess_sld(3,PLEPP_CPLNG%interior_list_j(ipoin),ITER_K) = data_recv(7,ipoin)*data_recv(10,ipoin) + &
                                 displ(3,PLEPP_CPLNG%interior_list_j(ipoin),TIME_N) + petol_sld
                         end if

                      end if

                   end do

                end if

                if(.not. associated(release_nodes)) then
                   allocate(release_nodes(nrecv))
                   release_nodes = 0_ip
                   new_to_release = 0_ip
                end if
                
             end if

          end if

          if ( CNT_SENDRECV(12) ) then
             !
             ! Bilateral
             !
             if( kfl_conta_sld == SLD_PDN_BILATERAL ) then
                if (inotmaster) then
                   !
                   CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0_rp
                   !
                   do ipoin = 1,nrecv
                      !
                      resid_tmp(1:ndime) = solve(1) % reaction(1:ndime,PLEPP_CPLNG%interior_list_j(ipoin))
                      !
                      ibopo = lpoty(PLEPP_CPLNG%interior_list_j(ipoin))
                      if( ibopo > 0_ip ) then
                         !
                         ! Local --> Global
                         !
                         call sld_rotsys(2_ip, 1_ip, 1_ip, ndime, ndime, dummy_matrix, resid_tmp(1:ndime), &
                              jacrot_du_dq_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)), &
                              jacrot_dq_du_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)))
                      end if
                      ! Save it
                      CNT_CPLNG%var_ij(1:ndime,PLEPP_CPLNG%interior_list_j(ipoin)) = resid_tmp(1:ndime)
                      !
                   end do
                   !
                end if
                !
                call plepp_pdn_driver_exchange02( CNT_CPLNG )
                !
                prop_in = 0.0_rp
                prop_out = 0.0_rp
                call plepp_pdn_dynamic_reduce_sum(prop_in, prop_out)
                if(prop_out > 0_rp) then
                   kfl_gozon = 0_ip
                   coupling_driver_iteration(iblok) = coupling_driver_max_iteration(iblok) - 1 !< !!???
                end if

             else
                
                if(associated(release_nodes)) deallocate(release_nodes)
                
             end if
             
          end if

       endif
       !
    endif

    !
    ! Deallocate arrays
    !
    deallocate(data_send)
    deallocate(data_recv)
    deallocate(norm_proje)
    deallocate(norm_proje2)
    deallocate(tangent_proje)
    deallocate(norm_dist)
    deallocate(in_contact)

  end subroutine commdom_sld_plugin_pdn

  !----------------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief   Projection nodes for 2-d problems
  !> @details Algorithm 5 PhD thesis Matias Rivero
  !----------------------------------------------------------------------------

  subroutine commdom_proje_nodes_2d(tangent_proje,normal_proje,normal_proje2,norm_dist,in_contact,e_recta)

    use def_domain,         only: coord, lnodb, kfl_codbo, nboun, ndime, mnodb, ltypb, nnode
    use def_master,         only: displ, IPARALL
    use def_kermod,         only: kfl_conta
    use mod_parall,         only: PAR_COMM_MY_CODE
    use mod_communications, only: PAR_ALLGATHER, PAR_ALLGATHERV, PAR_SUM, PAR_COMM_RANK_AND_SIZE, PAR_MAX
    use def_solidz,         only: contactbou_sld
    use def_solidz,         only: SLD_PDN_UNILATERAL

    implicit none

    real(rp), pointer, intent(inout)           :: tangent_proje(:,:)
    real(rp), pointer, intent(inout)           :: normal_proje(:,:)
    real(rp), pointer, intent(inout)           :: normal_proje2(:,:)
    real(rp), pointer, intent(inout)           :: norm_dist(:)
    real(rp), pointer, intent(inout)           :: in_contact(:)              !< Nodes in contact
    real(rp), intent(in)                       :: e_recta(ndime)             !< Projection direction
    real(rp)                                   :: coord_proje_tmp(ndime)
    integer(ip)                                :: nsend,nrecv
    real(rp), pointer                          :: dist_coords_j(:) => null() !< Coor. nodes body_j (A2) which are inside body_i (A1)
    real(rp), pointer                          :: dist_props_j(:) => null()
    real(rp), pointer                          :: coord_j(:,:) => null()
    real(rp), pointer                          :: orig_coord_j(:,:) => null()
    integer(ip)                                :: ipoin, idime, iboun, stride
    integer(ip)                                :: boidx(mnodb), pnodb, pblty, tmp, tmp2
    real(rp)                                   :: vec1(ndime), vec2(ndime), vec3(ndime), lambda
    real(rp)                                   :: dxc, dyc, cross, d1, d2, d3, d2x, d2y
    real(rp)                                   :: ax, ay, bx, by, x1, x2, y1, y2
    integer(ip)                                :: sizearray, contactbou_total, rankempi, sizempi
    integer(ip)                                :: contactbou_counter, dummy_count, mnodb_max
    integer(ip), pointer                       :: recvcounter(:) => null()
    real(rp), pointer                          :: bocod_total(:,:,:) => null()
    real(rp), pointer                          :: bocod_local(:,:,:) => null()

    nsend  = PLEPP_CPLNG%n_ij                     ! Nodos localizados del deformable respeto al rigido (indentador)
    nrecv  = PLEPP_CPLNG%n_ji                     ! Nodos localizados del rigido respeto al deformable
    stride = PLEPP_CPLNG%stride
    dist_coords_j => PLEPP_CPLNG%dist_coords_j    !
    dist_props_j  => PLEPP_CPLNG%dist_props_j     !

    if (stride /= ndime ) call runend("runend in commdom_proje_nodes_2d -- stride /= ndime")

    allocate(coord_j(ndime,nsend))           ! coordenadas deformable actuales
    allocate(orig_coord_j(ndime,nsend))      ! coordenadas deformable malla referencia.
    !
    ! Coordinates of ALL nodes of the other instance which have penetrated my current instance
    !
    if (inotmaster .and. (nsend > 0_ip)) then
       if ( kfl_conta == SLD_PDN_UNILATERAL ) then
          ! Old way (implicit/explicit unilateral)
          do ipoin = 1,nsend
             do idime = 1,ndime
                coord_j(idime,ipoin) = dist_coords_j( ndime*(ipoin-1)+idime ) ! coordenades del deformable
             end do
          end do
       else
          ! New way (bilateral)
          do ipoin = 1,nsend
             do idime = 1,ndime
                coord_j(idime,ipoin)      = dist_coords_j( ndime*(ipoin-1)+idime ) + dist_props_j( ndime*(ipoin-1)+idime )
                orig_coord_j(idime,ipoin) = dist_coords_j( ndime*(ipoin-1)+idime )
             end do
          end do
       end if
    end if

    !at this point I have the coordinates of ALL the nodes of the other code which have penetrated my code domain
    !now is time to project those points (restricted to the boundary) and all the planes defined by my boundary elements
    !last check is to see if the projection is inside the triangle or quadrange defined by the boundary element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--PARALLEL TWEAK - SEND CONTACT BOUNDARY TO ALL SLAVES--!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Get total number of boundary elements
    !
    contactbou_counter = 0_ip
    if (INOTMASTER) then
       do iboun = 1,nboun
          if (kfl_codbo(iboun) == contactbou_sld) then
             contactbou_counter = contactbou_counter + 1_ip
          end if
       end do
       mnodb_max = mnodb
    else
       contactbou_counter = 0_ip
       mnodb_max = 0_ip
    end if

    call PAR_MAX(mnodb_max,'IN MY CODE')
    sizearray = ndime*mnodb_max*contactbou_counter

    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE,rankempi,sizempi)
    allocate(recvcounter(sizempi))
    recvcounter = -666_ip
    if (IPARALL) then
       call PAR_ALLGATHER(sizearray,recvcounter,1_4,'IN MY CODE')
    else
       recvcounter = sizearray
    end if

    contactbou_total = contactbou_counter
    call PAR_SUM(contactbou_total, 'IN MY CODE') !this par_sum does not includes the master
    allocate(bocod_total(ndime,mnodb_max,contactbou_total))
    bocod_total = 0_rp

    if (INOTMASTER) then
       allocate(bocod_local(ndime,mnodb_max,contactbou_counter))
       bocod_local = 0_rp
       dummy_count = 0_ip
       do iboun = 1,nboun
          if (kfl_codbo(iboun) == contactbou_sld) then
             dummy_count = dummy_count + 1_ip
             pblty = ltypb(iboun) !boundary element type
             pnodb = nnode(pblty) !boundary number of nodes
             boidx(1:pnodb) = lnodb(1:pnodb, iboun) !boundary list of nodes
             bocod_local(1:ndime,1:pnodb,dummy_count) = coord(1:ndime,boidx(1:pnodb)) + displ(1:ndime,boidx(1:pnodb),1)
          end if
       end do
    else
       allocate(bocod_local(ndime,mnodb_max,0_ip))
    end if
    if (IPARALL) then
       call PAR_ALLGATHERV(bocod_local, bocod_total, recvcounter, 'IN MY CODE')
    else
       bocod_total = bocod_local
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( INOTMASTER ) then
       !
       ! Projections
       !
       tmp2 = 0_ip
       do ipoin = 1, nsend ! for each of the nodes that belong to the other code (those which plepp has located)
          tmp = 0_ip

          do iboun = 1, contactbou_total ! for each of all the boundaries defined as "contact boundary"
             !as each slave has all the information, I can do that
             !yb = vec1(1)x + vec1(2) --> line created using the boundary element
             vec1(1) = (bocod_total(2,1,iboun) - bocod_total(2,2,iboun)) / (bocod_total(1,1,iboun) - bocod_total(1,2,iboun)) !slope
             vec1(2) = bocod_total(2,1,iboun) - vec1(1)*bocod_total(1,1,iboun) !origin ordinate

             vec2(1) = bocod_total(1,1,iboun) - bocod_total(1,2,iboun) !x component - boundary tangent
             vec2(2) = bocod_total(2,1,iboun) - bocod_total(2,2,iboun) !y component - boundary tangent

             vec3(1) = -vec2(2) !x component - boundary normal
             vec3(2) = vec2(1) !y component - boundary normal

             !write the projection line in parametric way: (x,y) = (x0,y0) + lambda*(ex,ey)
             !then replace the parametric expression in the line descripted by vec1 to find the lambda
             lambda = (vec1(1)*coord_j(1,ipoin) + vec1(2) - coord_j(2,ipoin)) / (e_recta(2) - vec1(1)*e_recta(1))

             !then replace the lambda in the parametric expression to find the projection point
             coord_proje_tmp(1) = coord_j(1,ipoin) + lambda*e_recta(1) !x component of projected point
             coord_proje_tmp(2) = coord_j(2,ipoin) + lambda*e_recta(2) !y component of projected point
             !
             ! now I must check if the projected point is inside the boundary element
             !
             dxc = coord_proje_tmp(1) - bocod_total(1,2,iboun) !x component - vector projected point - boundary node
             dyc = coord_proje_tmp(2) - bocod_total(2,2,iboun) !y component - vector projected point - boundary node
             cross = dxc*vec2(2) - dyc*vec2(1) ! cross product to see if the vectors are colinear
             if (cross < 1e-6_rp) then !if they are colinear, check if the projected point is inside the boundary
                d1 = sqrt(dxc*dxc + dyc*dyc)
                d2x = coord_proje_tmp(1) - bocod_total(1,1,iboun)
                d2y = coord_proje_tmp(2) - bocod_total(2,1,iboun)
                d2 = sqrt(d2x*d2x + d2y*d2y)                                            !            _________ projected point
                d3 = sqrt(vec2(1)*vec2(1) + vec2(2)*vec2(2))                            !           /
                if ((d1 <= d3) .and. (d2 <= d3)) then  !the node is inside the boundary      x------o------x
                   tmp = tmp + 1_ip                                                      !   |--d2--|--d1--|
                   tangent_proje(1,ipoin) = vec1(1)                                      !   |-----d3------|
                   tangent_proje(2,ipoin) = vec1(2)
                   normal_proje(1,ipoin) = vec3(1)/sqrt(vec3(1)*vec3(1) + vec3(2)*vec3(2))
                   normal_proje(2,ipoin) = vec3(2)/sqrt(vec3(1)*vec3(1) + vec3(2)*vec3(2))
                   !now I have to calculate the norm_dist, which is the normal distance from the point to the line
                   !is the distance from the bounday point that I have to project and the tangent line that passes by its projection
                   norm_dist(ipoin) = abs(vec1(1)*coord_j(1,ipoin) - coord_j(2,ipoin) + vec1(2))/sqrt(vec1(1)*vec1(1) + 1)
                   tmp2 = tmp2 + 1_ip
                   !
                   if ( kfl_conta == SLD_PDN_UNILATERAL ) then
                      ! Old way: unilateral implicit/explicit
                      normal_proje2(1,ipoin) = normal_proje(1,ipoin)
                      normal_proje2(2,ipoin) = normal_proje(2,ipoin)
                      in_contact(ipoin) = 1.0_rp
                   else
                      ! New way (bilateral)
                      ax = orig_coord_j(1,ipoin)
                      ay = orig_coord_j(2,ipoin)
                      bx = coord_j(1,ipoin)
                      by = coord_j(2,ipoin)
                      x1 = bocod_total(1,1,iboun)
                      y1 = bocod_total(2,1,iboun)
                      x2 = bocod_total(1,2,iboun)
                      y2 = bocod_total(2,2,iboun)
                      if ( ((y1-y2)*(ax-x1)+(x2-x1)*(ay-y1))*((y1-y2)*(bx-x1)+(x2-x1)*(by-y1)) < 0_rp ) then
                         ! LADO OPUESTO, INVIERTO SIGNOS
                         normal_proje2(1,ipoin) = -1_ip*normal_proje(1,ipoin) !Gg
                         normal_proje2(2,ipoin) = -1_ip*normal_proje(2,ipoin) !GG
                         in_contact(ipoin) = 0.0_rp
                      else
                         ! DEL MISMO LADO, DEJO TODO COMO ESTA
                         normal_proje2(1,ipoin) = normal_proje(1,ipoin)
                         normal_proje2(2,ipoin) = normal_proje(2,ipoin)
                         in_contact(ipoin) = 1.0_rp
                      end if
                   end if

                end if

             end if

             if(tmp > 1_ip) call runend("ERROR IN COMMDOM_PROJE_NODES_2D: a projection is inside of more than one boundary element")

          end do

       end do

       if (tmp2 /= nsend) then
          call runend("ERROR IN COMMDOM_PROJE_NODES_2D: ME SOBRAN O ME FALTAN PROYECCIONES DE NODOS... BYE!!!!")
       end if

    end if

    deallocate(coord_j)
    deallocate(orig_coord_j)
    deallocate(recvcounter)
    deallocate(bocod_total)
    deallocate(bocod_local)

  end subroutine commdom_proje_nodes_2d

  !----------------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief   Projection nodes for 3-d problems
  !> @details Algorithm 5 PhD thesis Matias Rivero
  !----------------------------------------------------------------------------

  subroutine commdom_proje_nodes_3d(tangent_proje,normal_proje,normal_proje2,norm_dist,in_contact,e_recta)

    use def_domain,         only: coord, lnodb, kfl_codbo, nboun, ndime, mnodb, ltypb, nnode
    use def_master,         only: displ, IPARALL
    use def_kermod,         only: kfl_conta
    use def_elmtyp,         only: TRI03, QUA04
    use mod_parall,         only: PAR_COMM_MY_CODE
    use mod_communications, only: PAR_ALLGATHER, PAR_ALLGATHERV, PAR_SUM, PAR_COMM_RANK_AND_SIZE, PAR_MAX
    use def_solidz,         only: contactbou_sld
    use def_solidz,         only: SLD_PDN_UNILATERAL

    implicit none

    real(rp), pointer, intent(inout)           :: tangent_proje(:,:)
    real(rp), pointer, intent(inout)           :: normal_proje(:,:)
    real(rp), pointer, intent(inout)           :: normal_proje2(:,:)
    real(rp), pointer, intent(inout)           :: norm_dist(:)
    real(rp), pointer, intent(inout)           :: in_contact(:)
    real(rp), intent(in)                       :: e_recta(ndime)
    real(rp)                                   :: coord_proje_tmp(ndime)
    integer(ip)                                :: nsend,nrecv
    real(rp), pointer                          :: dist_coords_j(:) => null() !< Coor. nodes body_j (A2) which are inside body_i (A1)
    real(rp), pointer                          :: dist_props_j(:) => null()
    real(rp), pointer                          :: coord_j(:,:) => null()
    real(rp), pointer                          :: orig_coord_j(:,:) => null()
    integer(ip)                                :: ipoin, idime, iboun, ifoun, stride
    real(rp)                                   :: bari1, bari2
    integer(ip)                                :: boidx(mnodb), pnodb, pblty, ntria, tmp, tmp2
    real(rp)                                   :: vec1(ndime), vec2(ndime), vec3(ndime+1), t1, t2
    integer(ip)                                :: sizearray, sizearray2, contactbou_total, rankempi, sizempi
    integer(ip)                                :: contactbou_counter, dummy_count, mnodb_max
    integer(ip), pointer                       :: recvcounter(:)  => null()
    integer(ip), pointer                       :: recvcounter2(:) => null()
    real(rp), pointer                          :: bocod_total(:,:,:) => null()
    real(rp), pointer                          :: bocod_local(:,:,:) => null()
    integer(ip), pointer                       :: pblty_total(:) => null()
    integer(ip), pointer                       :: pblty_local(:) => null()
    integer(ip), pointer                       :: pnodb_total(:) => null()
    integer(ip), pointer                       :: pnodb_local(:) => null()
    real(rp)                                   :: pi_actual_coord, pi_origin_coord

    nsend  = PLEPP_CPLNG%n_ij
    nrecv  = PLEPP_CPLNG%n_ji
    stride = PLEPP_CPLNG%stride
    dist_coords_j => PLEPP_CPLNG%dist_coords_j
    dist_props_j  => PLEPP_CPLNG%dist_props_j

    if ( stride /= ndime ) call runend("runend in commdom_proje_nodes_3d -- stride /= ndime")

    allocate(coord_j(ndime,nsend))
    allocate(orig_coord_j(ndime,nsend))

    if ( INOTMASTER .and. nsend > 0_ip ) then
       if (kfl_conta == SLD_PDN_UNILATERAL ) then ! Old way (implicit/explicit unilateral)
          do ipoin = 1,nsend
             do idime = 1,ndime
                coord_j(idime,ipoin) = dist_coords_j( ndime*(ipoin-1)+idime )
             end do
          end do
       else                                       ! New way (bilateral)
          do ipoin = 1,nsend
             do idime = 1,ndime
                coord_j(idime,ipoin) = dist_coords_j( ndime*(ipoin-1)+idime ) + dist_props_j( ndime*(ipoin-1)+idime )
                orig_coord_j(idime,ipoin) = dist_coords_j( ndime*(ipoin-1)+idime )
             end do
          end do
       end if
    end if

    !at this point I have the coordinates of ALL the nodes of the other code which have penetrated my code domain
    !now is time to project those points (restricted to the boundary) and all the planes defined by my boundary elements
    !last check is to see if the projection is inside the triangle or quadrange defined by the boundary element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--PARALLEL TWEAK - SEND CONTACT BOUNDARY TO ALL SLAVES--!!!!!!!!!!!!!!!!!!!!!!!!!!
    contactbou_counter = 0_ip
    if ( INOTMASTER ) then
       do iboun = 1,nboun
          if (kfl_codbo(iboun) == contactbou_sld) then
             contactbou_counter = contactbou_counter + 1_ip
          end if
       end do
       mnodb_max = mnodb
    else
       contactbou_counter = 0_ip
       mnodb_max = 0_ip
    end if

    call PAR_MAX(mnodb_max,'IN MY CODE')
    sizearray  = ndime*mnodb_max*contactbou_counter
    sizearray2 = contactbou_counter

    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE,rankempi,sizempi)
    allocate(recvcounter( sizempi))
    allocate(recvcounter2(sizempi))
    recvcounter  = -666_ip
    recvcounter2 = -666_ip
    if (IPARALL) then
       call PAR_ALLGATHER(sizearray, recvcounter, 1_4,'IN MY CODE')
       call PAR_ALLGATHER(sizearray2,recvcounter2,1_4,'IN MY CODE')
    else
       recvcounter  = sizearray
       recvcounter2 = sizearray2
    end if

    contactbou_total = contactbou_counter
    call PAR_SUM(contactbou_total, 'IN MY CODE') !this par_sum does not includes the master
    allocate(pblty_total(contactbou_total))
    allocate(pnodb_total(contactbou_total))
    allocate(bocod_total(ndime,mnodb_max,contactbou_total))
    pblty_total = 0
    pnodb_total = 0
    bocod_total = 0.0_rp
    if( INOTMASTER ) then
       allocate(pblty_local(contactbou_counter))
       allocate(pnodb_local(contactbou_counter))
       allocate(bocod_local(ndime,mnodb_max,contactbou_counter))
       pblty_local = 0_ip
       pnodb_local = 0_ip
       bocod_local = 0.0_rp
       dummy_count = 0_ip
       do iboun = 1,nboun
          if( kfl_codbo(iboun) == contactbou_sld ) then
             dummy_count = dummy_count + 1_ip
             pblty = ltypb(iboun)                   ! boundary element type
             pnodb = nnode(pblty)                   ! boundary number of nodes
             boidx(1:pnodb) = lnodb(1:pnodb, iboun) ! boundary list of nodes
             pblty_local(dummy_count) = pblty
             pnodb_local(dummy_count) = pnodb
             bocod_local(1:ndime,1:pnodb,dummy_count) = coord(1:ndime,boidx(1:pnodb)) + displ(1:ndime,boidx(1:pnodb),1)
          end if
       end do
    else
       allocate(pblty_local(0_ip))
       allocate(pnodb_local(0_ip))
       allocate(bocod_local(ndime,mnodb_max,0_ip))
       pblty_local = 0_ip
       pnodb_local = 0_ip
       bocod_local = 0.0_rp
    end if
    
    if( IPARALL ) then
       call PAR_ALLGATHERV(pblty_local, pblty_total, recvcounter2, 'IN MY CODE')
       call PAR_ALLGATHERV(pnodb_local, pnodb_total, recvcounter2, 'IN MY CODE')
       call PAR_ALLGATHERV(bocod_local, bocod_total, recvcounter,  'IN MY CODE')
    else
       pblty_total = pblty_local
       pnodb_total = pnodb_local
       bocod_total = bocod_local
    end if
    
    if( INOTMASTER ) then

       tmp2 = 0_ip
       do ipoin = 1, nsend !for each of the nodes that belong to the other code (those which plepp has located)

          tmp = 0_ip
          do iboun = 1, contactbou_total ! For each of all the boundaries defined as "contact boundary"
             pblty = pblty_total(iboun)  ! boundary element type
             pnodb = pnodb_total(iboun)  ! boundary number of nodes
                             
             ! As each slave has all the information, I can do that
             vec1(1:ndime) = bocod_total(1:ndime,2,iboun) - bocod_total(1:ndime,1,iboun) !tangent vector 1
             vec2(1:ndime) = bocod_total(1:ndime,3,iboun) - bocod_total(1:ndime,1,iboun) !tangent vector 2

             ! cross product between vec1 and vec2 (vec3 = vec1 x vec2)
             ! vec3 is the normal vector of the plane defined by the boundary element
             ! PI: ax + by + cz + d = 0
             vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2) !x component (a)
             vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3) !y component (b)
             vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1) !z component (c)
             !independent term (d)
             vec3(4) = -(vec3(1)*bocod_total(1,1,iboun) + vec3(2)*bocod_total(2,1,iboun) + vec3(3)*bocod_total(3,1,iboun)) 

             t1 = -(vec3(1)*coord_j(1,ipoin) + vec3(2)*coord_j(2,ipoin) + vec3(3)*coord_j(3,ipoin) + vec3(4)) !-PI(coord_j)
             t2 = vec3(1)*e_recta(1) + vec3(2)*e_recta(2) + vec3(3)*e_recta(3) !n.u

             !if ((t2 .eq. 0_rp) .or. (t1 .eq. 0_rp)) then
             !  call runend("runend in commdom_proje_nodes_3d - point to project is already on the plane")
             !end if
             !L(PI) = L(coord_j) + (t1/t2)*u
             !coord_proje_tmp(1) = coord_j(1,ipoin) + (t1/t2)*e_recta(1) !x component of projected point
             !coord_proje_tmp(2) = coord_j(2,ipoin) + (t1/t2)*e_recta(2) !y component of projected point
             !coord_proje_tmp(3) = coord_j(3,ipoin) + (t1/t2)*e_recta(3) !z component of projected point

             if ((t2 .eq. 0_rp) .or. (t1 .eq. 0_rp)) then
                !L(PI) = L(coord_j) + (t1/t2)*u
                coord_proje_tmp(1) = coord_j(1,ipoin)  !x component of projected point
                coord_proje_tmp(2) = coord_j(2,ipoin)  !y component of projected point
                coord_proje_tmp(3) = coord_j(3,ipoin)  !z component of projected point
             else
                !L(PI) = L(coord_j) + (t1/t2)*u
                coord_proje_tmp(1) = coord_j(1,ipoin) + (t1/t2)*e_recta(1) !x component of projected point
                coord_proje_tmp(2) = coord_j(2,ipoin) + (t1/t2)*e_recta(2) !y component of projected point
                coord_proje_tmp(3) = coord_j(3,ipoin) + (t1/t2)*e_recta(3) !z component of projected point
             end if

             ! Now I must check if the projected point is inside the boundary element
             ! usage of instr2 (kernel/domain/instri.f90), which is valid for linear traingles and quadrangles
             ifoun = 0_ip
             if ( (pblty .eq. TRI03) .or. (pblty .eq. QUA04) ) then
                call instr2(pnodb,coord_proje_tmp(:),bocod_total(:,:,iboun),ifoun,bari1,bari2,ntria)
             else
                call runend("runend in commdom_proje_nodes_3d -- projections other than to TRI03 or QUA04 are not implemented")
             end if
             if ( ifoun == 1_ip ) then !the projected point is inside the triangle or quad
                tmp = tmp + 1_ip
                tangent_proje(1,ipoin) = vec1(1)/sqrt(vec1(1)*vec1(1) + vec1(2)*vec1(2) + vec1(3)*vec1(3))
                tangent_proje(2,ipoin) = vec1(2)/sqrt(vec1(1)*vec1(1) + vec1(2)*vec1(2) + vec1(3)*vec1(3))
                tangent_proje(3,ipoin) = vec1(3)/sqrt(vec1(1)*vec1(1) + vec1(2)*vec1(2) + vec1(3)*vec1(3))
                normal_proje(1,ipoin) = vec3(1)/sqrt(vec3(1)*vec3(1) + vec3(2)*vec3(2) + vec3(3)*vec3(3))
                normal_proje(2,ipoin) = vec3(2)/sqrt(vec3(1)*vec3(1) + vec3(2)*vec3(2) + vec3(3)*vec3(3))
                normal_proje(3,ipoin) = vec3(3)/sqrt(vec3(1)*vec3(1) + vec3(2)*vec3(2) + vec3(3)*vec3(3))
                !
                if ( kfl_conta == SLD_PDN_UNILATERAL ) then
                   ! Old way: unilateral implicit/explicit
                   normal_proje2(1,ipoin) = normal_proje(1,ipoin)
                   normal_proje2(2,ipoin) = normal_proje(2,ipoin)
                   normal_proje2(3,ipoin) = normal_proje(3,ipoin)
                   in_contact(ipoin) = 1.0_rp
                else
                   ! New way (bilateral)
                   pi_actual_coord = vec3(1)*coord_j(1,ipoin) + vec3(2)*coord_j(2,ipoin) + vec3(3)*coord_j(3,ipoin) + vec3(4)
                   pi_origin_coord = vec3(1)*orig_coord_j(1,ipoin) + &
                        vec3(2)*orig_coord_j(2,ipoin) + vec3(3)*orig_coord_j(3,ipoin) + vec3(4)
                   if ( pi_origin_coord*pi_actual_coord < 0_rp ) then           !opposite side
                      normal_proje2(1,ipoin) = -1_ip*normal_proje(1,ipoin)
                      normal_proje2(2,ipoin) = -1_ip*normal_proje(2,ipoin)
                      normal_proje2(3,ipoin) = -1_ip*normal_proje(3,ipoin)
                      in_contact(ipoin) = 0.0_rp
                   else       !same side
                      normal_proje2(1,ipoin) = normal_proje(1,ipoin)
                      normal_proje2(2,ipoin) = normal_proje(2,ipoin)
                      normal_proje2(3,ipoin) = normal_proje(3,ipoin)
                      in_contact(ipoin) = 1.0_rp
                   end if
                end if

                !now I have to calculate the norm_dist, which is the normal distance from the point to the plane
                !is the distance from the bounday point that I have to project and the tangent plane that passes by its projection
                !the point to project is given by coord_j and the plane has a normal given by vec3(:) and passes through
                !coord_proje(:)
                !recalculate the independent term of tangent plane
                vec3(4) = -(vec3(1)*coord_proje_tmp(1) + vec3(2)*coord_proje_tmp(2) + vec3(3)*coord_proje_tmp(3)) 
                norm_dist(ipoin) = abs(vec3(1)*coord_j(1,ipoin) + vec3(2)*coord_j(2,ipoin) + vec3(3)*coord_j(3,ipoin) + &
                     vec3(4))/sqrt(vec3(1)*vec3(1) + vec3(2)*vec3(2) + vec3(3)*vec3(3))
                tmp2 = tmp2 + 1_ip
             end if
             if (tmp > 1_ip) then
                call runend("ERROR IN COMMDOM_PROJE_NODES_3D -- a projection is inside of more than one boundary element")
             end if
          end do

       end do

       if (tmp2 /= nsend) then
          call runend("ERROR IN COMMDOM_PROJE_NODES_3D: ME SOBRAN O ME FALTAN PROYECCIONES DE NODOS... BYE!!!!")
       end if

    end if

    deallocate(coord_j)
    deallocate(orig_coord_j)
    deallocate(recvcounter)
    deallocate(pblty_total)
    deallocate(pblty_local)
    deallocate(pnodb_total)
    deallocate(pnodb_local)
    deallocate(bocod_total)
    deallocate(bocod_local)

  end subroutine commdom_proje_nodes_3d

  !-----------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details Bilateral (Neumann part)
  !-----------------------------------------------------------------------

  subroutine commdom_sld_bilateral_neumann()

    use def_master,          only : iblok, IPARALL
    use def_domain,          only : ndime
    use def_coupli,          only : coupling_driver_iteration, coupling_driver_max_iteration, kfl_gozon
    use mod_communications,  only : PAR_SUM
    use def_solidz,          only : release_nodes, reaction_ant
    use def_solidz,          only : neumann_relax, coupling_contact_tol
    use def_solidz,          only : relax_bilateral

    implicit none

    integer(ip)       :: ipoin, idime
    integer(ip)       :: nrecv
    real(rp)          :: prop_in, prop_out
    real(rp)          :: res_diff(ndime), res_reac(ndime)
    real(rp)          :: numer, denom, rip1, rip2

    solve => momod(modul) % solve(1:)
    nrecv = PLEPP_CPLNG%n_ji
    !
    ! Initializations
    !
    do ipoin = 1,npoin
       CNT_CPLNG%var_ji(1:ndime,ipoin) = 0.0_rp
    end do
    !
    nullify(reaction_ant)
    if(.not. associated(reaction_ant)) then
       allocate(reaction_ant(ndime,npoin,4_ip))
       reaction_ant = 0.0_rp
    end if
    !
    call plepp_pdn_driver_exchange02( CNT_CPLNG )
    !
    ! Relaxation
    !
    do ipoin = 1,npoin
       reaction_ant(1:ndime,ipoin,3_ip) = reaction_ant(1:ndime,ipoin,2_ip)
       reaction_ant(1:ndime,ipoin,2_ip) = reaction_ant(1:ndime,ipoin,1_ip)
    end do
    !
    if ( coupling_driver_iteration(iblok) < 3_ip*10000_ip ) then
       relax_bilateral = neumann_relax
    else
       numer = 0.0_rp
       denom = 0.0_rp
       do ipoin = 1,nrecv
          do idime = 1,ndime
             rip1 = reaction_ant(idime,PLEPP_CPLNG%interior_list_j(ipoin),4_ip) - &
                  reaction_ant(idime,PLEPP_CPLNG%interior_list_j(ipoin),3_ip)
             rip2 = CNT_CPLNG%var_ji(idime,PLEPP_CPLNG%interior_list_j(ipoin))  - &
                  reaction_ant(idime,PLEPP_CPLNG%interior_list_j(ipoin),2_ip)
             numer = numer + rip1 * (rip2-rip1)
             denom = denom + (rip2-rip1)*(rip2-rip1)
          end do
       end do
       relax_bilateral = relax_bilateral * (-1_rp*numer / (denom+zeror))
       !
    end if
    !
    ! Initialize force received from other instance
    !
    do ipoin = 1,npoin
       solve(1) % bvnat(1:ndime,ipoin) = 0.0_rp
    end do
    if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
       if( INOTMASTER ) then
          solve(1) % bvnat(1:ndime,PLEPP_CPLNG%interior_list_j) = &
               relax_bilateral * CNT_CPLNG%var_ji(1:ndime,PLEPP_CPLNG%interior_list_j) + &
               (1.0_rp - relax_bilateral) * reaction_ant(1:ndime,PLEPP_CPLNG%interior_list_j,2_ip)
       end if
    end if
    !
    if( IPARALL .and. INOTMASTER ) then
       solve(1) % bvnat(:,(npoi3+1):npoin) = 0.0_rp
    end if
    !
    ! Calculate residual
    !
    if ( INOTMASTER ) then
       do idime = 1,ndime
          res_diff(idime) = dot_product( (solve(1) % bvnat(idime,PLEPP_CPLNG%interior_list_j) - &
               reaction_ant(idime,PLEPP_CPLNG%interior_list_j,1_ip)), &
               (solve(1) % bvnat(idime,PLEPP_CPLNG%interior_list_j) - &
               reaction_ant(idime,PLEPP_CPLNG%interior_list_j,1_ip)) )
          res_reac(idime) = &
               dot_product(solve(1)%bvnat(idime,PLEPP_CPLNG%interior_list_j),solve(1)%bvnat(idime,PLEPP_CPLNG%interior_list_j))
       end do
    end if
    !
    ! Assign reaction to the other instance
    !
    do ipoin = 1,npoin
       reaction_ant(1:ndime,ipoin,1_ip) = solve(1) % bvnat(1:ndime,ipoin)
       reaction_ant(1:ndime,ipoin,4_ip) = CNT_CPLNG%var_ji(1:ndime,ipoin)
    end do
    !
    ! Parall
    !
    call PAR_SUM(ndime, res_diff(1:ndime), 'IN MY CODE')
    call PAR_SUM(ndime, res_reac(1:ndime), 'IN MY CODE')
    call PAR_SUM(nrecv, 'IN MY CODE')
    !
    ! Residual coupling contact
    !
    res_diff(1:ndime) = res_diff(1:ndime)/(res_reac(1:ndime) + zeror)
    prop_in = 0.0_rp
    prop_out = 0.0_rp
    if(nrecv == 0_ip) then
       prop_in = 1.0_rp
    else
       if( all(res_diff(1:ndime) <= coupling_contact_tol) ) prop_in = 1.0_rp
    end if
    !
    ! Live info
    !
    print *,"ITERACION DE ACOPLE -- RESIDUO ACOPLE CONTACTO --> ", coupling_driver_iteration(iblok),&
         " ------ ", res_diff(1:ndime), "---", relax_bilateral
    !
    ! Parall
    !
    call plepp_pdn_dynamic_reduce_sum(prop_in, prop_out)
    if(prop_out > 0_rp) then
       kfl_gozon = 0_ip
       coupling_driver_iteration(iblok) = coupling_driver_max_iteration(iblok) - 1 !< !!???
    end if
    !
    if(associated(release_nodes)) deallocate(release_nodes)
   
  end subroutine commdom_sld_bilateral_neumann

  !----------------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details
  !>
  !    Determine if a point is inside a triangle using the same side technique
  !    point: point
  !    facoo(dimension,vertices): triangle coordinates
  !    ifoun: If is equal to 1, the point is inside the triangle, 0 otherwise
  ! USED BY
  !    pofadi,faliib
  !***
  !----------------------------------------------------------------------------

  subroutine instr2(pnodb,point,facoo,ifoun,bari1,bari2,ntria)

    use def_kintyp, only           :  ip,rp
    use def_domain, only           :  ndime
    use def_master, only           :  zeror

    implicit   none
    integer(ip),   intent(in)     :: pnodb
    integer(ip),   intent(out)    :: ifoun
    integer(ip),   intent(out)    :: ntria
    real(rp),      intent(in)     :: point(ndime)
    real(rp),      intent(in)     :: facoo(ndime,*)
    real(rp),      intent(out)    :: bari1,bari2
    real(rp)                      :: v0(3),v1(3),v2(3)
    real(rp)                      :: dot00,dot01,dot02,dot11,dot12
    real(rp)                      :: invDenom
    !
    ! 3D
    !
    v0(1) = facoo(1,3) - facoo(1,1)
    v0(2) = facoo(2,3) - facoo(2,1)
    v0(3) = facoo(3,3) - facoo(3,1)

    v1(1) = facoo(1,2) - facoo(1,1)
    v1(2) = facoo(2,2) - facoo(2,1)
    v1(3) = facoo(3,2) - facoo(3,1)

    v2(1) = point(1)   - facoo(1,1)
    v2(2) = point(2)   - facoo(2,1)
    v2(3) = point(3)   - facoo(3,1)

    dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
    dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
    dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3)
    dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
    dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    !
    ! Compute barycentric coordinates
    !
    ntria = 0
    if (abs(dot00 * dot11 - dot01 * dot01) > zeror) then
       invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
       bari1 = (dot11 * dot02 - dot01 * dot12) * invDenom
       bari2 = (dot00 * dot12 - dot01 * dot02) * invDenom
       !
       ! Check
       !
       if( bari1 >= 0.0_rp .and. bari2 >= 0.0_rp .and. bari1 + bari2 <= 1.0_rp ) then
          ntria = 1
          ifoun = 1
       else
          ifoun = 0
       end if
    else
       ifoun = 0
    end if
    if ( pnodb == 4 .and. ifoun == 0 ) then
       v0(1) = facoo(1,4) - facoo(1,1)
       v0(2) = facoo(2,4) - facoo(2,1)
       v0(3) = facoo(3,4) - facoo(3,1)

       v1(1) = facoo(1,3) - facoo(1,1)
       v1(2) = facoo(2,3) - facoo(2,1)
       v1(3) = facoo(3,3) - facoo(3,1)

       v2(1) = point(1)   - facoo(1,1)
       v2(2) = point(2)   - facoo(2,1)
       v2(3) = point(3)   - facoo(3,1)

       dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
       dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
       dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3)
       dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
       dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

       if (abs(dot00 * dot11 - dot01 * dot01) > zeror) then
          invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
          bari1 = (dot11 * dot02 - dot01 * dot12) * invDenom
          bari2 = (dot00 * dot12 - dot01 * dot02) * invDenom
          !
          ! Check
          !
          if( bari1 >= 0.0_rp .and. bari2 >= 0.0_rp .and. bari1 + bari2 <= 1.0_rp ) then
             ifoun = 1
             ntria = 2
          else
             ifoun = 0
          end if
       else
          ifoun = 0
       end if
    end if

  end subroutine instr2

  !-----------------------------------------------------------------------
  !> @author  Matias Rivero and Gerard Guillamet
  !> @date    2019, March
  !> @brief   Calculation of reaction and nodes release algorithm
  !> @details Calculation of reaction and nodes release algorithm
  !-----------------------------------------------------------------------

  subroutine commdom_sld_nodes_release()

    use mod_communications,   only : PAR_SUM
    use def_solidz,           only : kfl_goite_sld
    use def_solidz,           only : new_to_release
    use def_solidz,           only : kfl_mrele_sld
    use def_solidz,           only : jacrot_du_dq_sld,jacrot_dq_du_sld
    use def_master,           only : ITER_K
    use mod_sld_pdn_contact_plepp, only : narco_sld,nacco_sld,ntoco_sld,itrel_sld
    
    implicit none
    
    integer(ip)                   :: ipoin,nrecv,irecv
    !
    ! Check existance of sticky nodes
    !
    if( kfl_mrele_sld == 1_ip ) then

       if( kfl_goite_sld /= 1_ip ) then     ! if converged
 
          ! Compute contact reaction and mark nodes in traction
          call commdom_sld_calc_reaction()
          if( INOTMASTER ) then
             if( associated(release_nodes) ) narco_sld = sum(release_nodes)
          end if
          ! Release nodes (if any) and repeat iteration
          if( new_to_release == 1_ip ) then
             ! Ini flag
             new_to_release = 0_ip
             if( INOTMASTER ) then
                nrecv = PLEPP_CPLNG%n_ji
                do irecv = 1,nrecv
                   ipoin = PLEPP_CPLNG%interior_list_j(irecv)
                   if( kfl_fixno_sld(1,ipoin) == 3_ip .and. release_nodes(irecv) == 1_ip  ) then
                      kfl_fixno_sld(1:ndime,ipoin)            = 0_ip
                      bvess_sld(1:ndime,ipoin,ITER_K)         = 0.0_rp
                      jacrot_du_dq_sld(1:ndime,1:ndime,ipoin) = 0.0_rp
                      jacrot_dq_du_sld(1:ndime,1:ndime,ipoin) = 0.0_rp
                      nacco_sld                               = nacco_sld - 1_ip
                   end if
                end do
             end if
             ! Repeat iteration (not converged)
             kfl_goite_sld = 1_ip
             itrel_sld = itrel_sld + 1_ip
          end if
          
       end if
       
    else
       
       if( kfl_goite_sld /= 1_ip ) then     ! if converged

          ! Compute contact reaction and mark nodes in traction
          call commdom_sld_calc_reaction()
          if( INOTMASTER ) then
             if( associated(release_nodes) ) narco_sld = sum(release_nodes)
          end if
          ! Release nodes (if any) and repeat iteration
          if( new_to_release == 1_ip ) then
             ! Ini flag
             new_to_release = 0_ip
             call commdom_sld_release_nodes()
             ! Repeat iteration (not converged)
             kfl_goite_sld = 1_ip
             itrel_sld = itrel_sld + 1_ip
          end if

       end if
       
    end if
    !
    ! Check if all subdomains are converged to solution
    ! If the sum is greater than 1 (not converged) repeat iteration
    call PAR_SUM(kfl_goite_sld, 'IN MY CODE')
    if( kfl_goite_sld >= 1_ip ) kfl_goite_sld = 1_ip
    !
    ! Parall of some convergence parameters
    !
    call PAR_SUM(nacco_sld, 'IN MY CODE')
    call PAR_SUM(ntoco_sld, 'IN MY CODE')
    call PAR_SUM(narco_sld, 'IN MY CODE')
    call PAR_MAX(itrel_sld, 'IN MY CODE')
       
  end subroutine commdom_sld_nodes_release

  !-----------------------------------------------------------------------
  !> @author  Matias Rivero and Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Calculation of contact force
  !> @details Calculation of contact force and others
  !-----------------------------------------------------------------------

  subroutine commdom_sld_calc_reaction()

    use def_domain,          only : lpoty, ndime
    use def_solidz,          only : frxid_sld, fcont_sld
    use def_solidz,          only : SLD_IMPLICIT_SCHEME,kfl_timet_sld
    use def_solidz,          only : release_nodes, new_to_release
    
    implicit none

    integer(ip) :: ipoin, irecv, nrecv, itott, ibopo
    real(rp)    :: fcont(ndime)
    logical(lg) :: algebraic
    
    solve => momod(modul) % solve(1:)

    !------------------------------------------------------------------
    !
    !  Deformable body (Dirichlet)
    !
    !------------------------------------------------------------------
    !
    ! Initialize
    !
    fcont = 0.0_rp

    if( INOTMASTER ) then

       nrecv = PLEPP_CPLNG%n_ji
       !
       ! Get contact force and mark nodes
       !
       algebraic = .false.
       do irecv = 1, nrecv
          ipoin = PLEPP_CPLNG%interior_list_j(irecv)
          ibopo = lpoty(ipoin)
          itott = (ipoin - 1) * ndime
          !
          ! Nodes in contact
          !
          fcont = 0.0_rp
          if( kfl_fixno_sld(1,ipoin) == 3_ip ) then 
             !
             ! Reaction method
             !
             if( kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then
                if( algebraic ) then
                   ! Reaction from solver
                   fcont(1:ndime) = solve(1) % reaction(1:ndime,ipoin)  ! already exchanged by solver
                else
                   ! Reaction from solidz
                   fcont(1:ndime) = frxid_sld(itott+1)            ! already exchanged
                end if
             else
                fcont(1:ndime) = frxid_sld(itott+1)
             end if
             !
             ! Mark nodes in traction
             !
             if( fcont(1) > 0.0_rp .and. release_nodes(irecv) == 0_ip ) then ! not marked
                release_nodes(irecv) = 1_ip
                new_to_release    = 1_ip
             end if
             !
             ! Rotate back contact force (Local --> Global)
             !
             call sld_csys_rotuni(2_ip,ndime,ipoin,fcont_sld(1:ndime,ipoin)) ! Not exchanged yet!

          end if

       end do

    end if

  end subroutine commdom_sld_calc_reaction

  !-----------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details Release nodes
  !-----------------------------------------------------------------------

  subroutine commdom_sld_release_nodes()

    use def_solidz,                only : release_nodes, kfl_fixno_sld
    use mod_sld_pdn_contact_plepp, only : nacco_sld
    
    implicit none

    integer(ip) :: ipoin, nrecv, irecv
    
    !------------------------------------------------------------------
    !
    !  Deformable body (Dirichlet)
    !
    !------------------------------------------------------------------
    
    if( INOTMASTER ) then
       nrecv = PLEPP_CPLNG%n_ji

       do irecv = 1,nrecv
          ipoin = PLEPP_CPLNG%interior_list_j(irecv)

          if( kfl_fixno_sld(1,ipoin) == 3_ip .and. release_nodes(irecv) == 1_ip ) then    
             !
             ! Set free node
             !
             kfl_fixno_sld(1:ndime,ipoin) = 0_ip
             nacco_sld                    = nacco_sld - 1_ip

          end if

       end do
    end if

  end subroutine commdom_sld_release_nodes

  !-----------------------------------------------------------------------
  !>
  !> @author  Gerard Guillamet and Miguel Zavala
  !> @date    March, 2019
  !> @brief   Projections from RBO to the DBO
  !> @details Nodes detected inside RBO are projected to the boundary of DBO
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_conta_projection_send(ndata,nsend,nrecv,data_send,data_recv)

    use def_solidz,          only : vect_proje_sld

    implicit none
    !
    integer(ip), intent(in)    :: ndata
    integer(ip), intent(in)    :: nsend
    integer(ip), intent(in)    :: nrecv
    real(rp),    intent(inout) :: data_send(ndata,nsend)
    real(rp),    intent(inout) :: data_recv(ndata,nrecv)
    !
    integer(ip)       :: isend
    real(rp)          :: e_recta(ndime)
    real(rp), pointer :: tangent_proje(:,:) => null()
    real(rp), pointer :: norm_proje(:,:) => null()
    real(rp), pointer :: norm_proje2(:,:) => null()
    real(rp), pointer :: norm_dist(:) => null()
    real(rp), pointer :: in_contact(:) => null()

    !------------------------------------------------------------------
    !
    ! Rigid body
    !
    !------------------------------------------------------------------
    !
    ! Allocate memory for arrays
    !
    allocate(tangent_proje(ndime,nsend))
    allocate(norm_proje(ndime,nsend))
    allocate(norm_proje2(ndime,nsend))
    allocate(norm_dist(nsend))
    allocate(in_contact(nsend))
    !
    ! Initialize
    !
    tangent_proje(:,:) = 0.0_rp
    norm_proje(:,:)    = 0.0_rp
    norm_dist(:)       = 0.0_rp
    norm_proje2(:,:)   = 0.0_rp
    in_contact(:)      = 666.0_rp
    !
    ! Projection vector --> CHANGE THIS: MUST BE DYNAMICALLY DEFINED
    !
    e_recta(1:ndime) = vect_proje_sld(1:ndime)
    !
    ! Projections: Nodes detected inside RBO to the boundary
    !
    if ( ndime == 2_ip ) then
       call commdom_proje_nodes_2d(tangent_proje,norm_proje,norm_proje2,norm_dist,in_contact,e_recta)
    else
       call commdom_proje_nodes_3d(tangent_proje,norm_proje,norm_proje2,norm_dist,in_contact,e_recta)
    end if
    !
    do isend = 1,nsend
       data_send(            1:ndime,isend) = tangent_proje(1:ndime,isend)
       data_send(    ndime+1:ndime*2,isend) = norm_proje(1:ndime,isend)
       data_send(          ndime*2+1,isend) = norm_dist(isend)
       data_send(ndime*2+2:ndime*3+1,isend) = norm_proje2(1:ndime,isend)
       data_send(          ndime*3+2,isend) = in_contact(isend)
    end do
    !
    call commdom_locator_exchange_double_stride(data_send(1:ndata,1:nsend),data_recv(1:ndata,1:nrecv),ndata)
    !
    ! Deallocate memory for arrays
    !
    deallocate(norm_proje)
    deallocate(norm_proje2)
    deallocate(tangent_proje)
    deallocate(norm_dist)
    deallocate(in_contact)
    !
  end subroutine sld_pdn_conta_projection_send

  !-----------------------------------------------------------------------
  !>
  !> @author  Gerard Guillamet and Miguel Zavala
  !> @date    March, 2019
  !> @brief   Recieve projections at DBO, mark nodes and apply Dirichlet
  !> @details Recieve projections at DBO, mark nodes and apply Dirichlet
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_conta_projection_recv(ndata,nsend,nrecv,data_send,data_recv)

    use def_master,           only : displ, TIME_N, ITER_K
    use def_domain,           only : ndime, kfl_codno, lpoty
    use mod_maths,            only : maths_vectorial_product
    use mod_maths,            only : maths_normalize_vector
    use def_solidz,           only : jacrot_du_dq_sld, jacrot_dq_du_sld, contactbou_sld, petol_sld
    use def_solidz,           only : kfl_timet_sld, SLD_IMPLICIT_SCHEME
    use mod_sld_pdn_contact_plepp, only : nacco_sld,ntoco_sld
    
    implicit none

    integer(ip),intent(in)    :: ndata
    integer(ip),intent(in)    :: nsend
    integer(ip),intent(in)    :: nrecv
    real(rp),   intent(inout) :: data_send(ndata,nsend)
    real(rp),   intent(inout) :: data_recv(ndata,nrecv)

    integer(ip)               :: irecv, ipoin, idime, ibopo
    real(rp)                  :: bvess_aux(ndime)
    real(rp)                  :: N(ndime),T1(ndime),T2(ndime),rho

    !------------------------------------------------------------------
    !
    ! Deformable body
    !
    !------------------------------------------------------------------
    
    if( INOTMASTER ) then
       !
       ! Initialize all nodes
       !
       do ipoin = 1,npoin
          do idime = 1,ndime
             if ( kfl_fixno_sld(idime,ipoin) == 3_ip ) kfl_fixno_sld(idime,ipoin) = 0_ip
          end do
       end do
       !
       call commdom_locator_exchange_double_stride(data_send(1:ndata,1:nsend),data_recv(1:ndata,1:nrecv),ndata)
       !
       ! Loop nodes recieved from PLEPP
       !
       do irecv = 1, nrecv

          ipoin = PLEPP_CPLNG%interior_list_j(irecv)
          ibopo = lpoty(ipoin)

          if( any(kfl_codno(:,ipoin) == contactbou_sld) ) then       ! nodes that belong to the boundary
             !
             ! Get normal and tangent vectors
             !
             N(:) = data_recv(ndime+1:ndime*2,irecv)
             if ( ndime == 2 ) then
                T1(1) = -data_recv(ndime*2,irecv)
                T1(2) =  data_recv(ndime+1,irecv)
             else
                T1(:) = data_recv(1:ndime,irecv)
                ! T2 =  N x T1
                call maths_vectorial_product(N(:), T1(:), T2(:), ndime)
                call maths_normalize_vector(ndime,T2(:),rho)
             end if
             !
             ! Build rotation matrix
             !
             jacrot_du_dq_sld(1:ndime,1,ipoin) = N(1:ndime)
             jacrot_du_dq_sld(1:ndime,2,ipoin) = T1(1:ndime)
             if( ndime == 3_ip ) jacrot_du_dq_sld(1:ndime,3,ipoin) = T2(1:ndime)
             ! Transpose of R
             jacrot_dq_du_sld(:,:,ipoin) = transpose(jacrot_du_dq_sld(:,:,ipoin))
             !
             ! Contact nodes with fixity in the first dof
             !
             if( kfl_fixno_sld(1,ipoin) == 1_ip ) then
                do idime = 1,ndime-1
                   kfl_fixno_sld(idime+1,ipoin) = 1_ip
                end do
             end if
             !
             ! Mark contact node
             !
             kfl_fixno_sld(1,ipoin) = 3_ip
             ntoco_sld              = ntoco_sld + 1
             nacco_sld              = nacco_sld + 1
             !
             ! Prescribe Dirichlet
             !
             if ( kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then
                bvess_aux(1:ndime) = data_recv(ndime*2+1,irecv)*data_recv(ndime*2+2:ndime*3+1,irecv) + petol_sld
                displ(1:ndime,ipoin,TIME_N) = bvess_aux(1:ndime) + displ(1:ndime,ipoin,TIME_N)
             else
                bvess_sld(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,TIME_N) + &
                     data_recv(ndime*2+1,irecv)*data_recv(ndime*2+2:ndime*3+1,irecv) - &
                     data_recv(ndime*2+2:ndime*3+1,irecv)*abs(petol_sld)
             end if
                          
          end if

       end do
       !
       ! Allocate memory for release nodes array
       !
       if (.not. associated(release_nodes)) then
          allocate(release_nodes(nrecv))
          release_nodes(:)  = 0_ip
          new_to_release = 0_ip
       end if
       !
    end if
    !
  end subroutine sld_pdn_conta_projection_recv

  !-----------------------------------------------------------------------
  !>
  !> @author  Miguel Zavala and Gerard Guillamet
  !> @date    2019-03-08
  !> @brief   Send reaction force from DBO to RBO
  !> @details Send reaction force from DBO to RBO
  !>
  !-----------------------------------------------------------------------

  subroutine sld_rbo_bilateral_send()

    use def_master,   only : mem_modul
    use mod_memory,   only : memory_alloca
    use mod_memory,   only : memory_deallo
    use def_solidz,   only : fcont_sld

    implicit none

    integer(ip)       :: ipoin
    real(rp)          :: force_dbo(ndime)

    !------------------------------------------------------------------
    !
    ! Deformable body
    !
    !------------------------------------------------------------------
    !
    ! Compatibility checks
    !
    solve(1_ip:) => momod(modul) % solve(1_ip:)
    if(solve(1)%kfl_bvnat==1) call runend("solve(1)%kfl_bvnat==1")  ! si es Neumaan sale!!
    !
    ! Calculate total contact force
    !
    force_dbo = 0.0_rp
    if( INOTMASTER ) then
       force_dbo(1:ndime) = sum(fcont_sld(1:ndime,:), dim=2)
       !
       ! Parall exchange
       !
       call rhsmod(ndime,fcont_sld)
    end if
    !
    ! Parall service
    !
    call PAR_SUM(ndime,force_dbo(1:ndime),'IN MY CODE')
    !
    ! We fill this matrix with the obtained value
    !
    CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0_rp
    if( INOTMASTER ) then
       if( PLEPP_CPLNG%n_ij > 0_ip ) then
          do ipoin = 1,npoin
             CNT_CPLNG%var_ij(1:ndime,ipoin) = -force_dbo(1:ndime)
          end do
       end if
    end if
    !
    ! Exchange data between instancess
    !
    call plepp_pdn_driver_exchange02( CNT_CPLNG ) !-> SEND
    !
    ! Live info
    !
    call messages_live('-> SEND '//trim('FORCE')//' AS A SOURCE FOR COUPLING '//trim(intost(1_ip)))
    !
    ! Deallocate release node vector
    !
    if(associated(release_nodes)) deallocate(release_nodes)

  end subroutine sld_rbo_bilateral_send

  !-----------------------------------------------------------------------
  !>
  !> @author  Miguel Zavala and Gerard Guillamet
  !> @date    2019-03-08
  !> @brief   Receive force from DBO to RBO
  !> @details Receive force from DBO to RBO
  !>
  !-----------------------------------------------------------------------

  subroutine sld_rbo_bilateral_recv( )

    use mod_communications, only : PAR_MIN, PAR_MAX
    use def_solidz,         only : rbfor_sld

    implicit none

    real(rp)    :: dummin, dummax

    !------------------------------------------------------------------
    !
    ! Rigid body
    !
    !------------------------------------------------------------------
    !
    ! Initialize
    !
    CNT_CPLNG%var_ji(1:ndime,1:npoin) = 0.0_rp
    !
    ! Recieve force from other instance
    !
    call plepp_pdn_driver_exchange02( CNT_CPLNG ) !<- RECV
    !
    ! Compatibility checks
    !
    solve(1:) => momod(modul) % solve(1:)
    if(solve(1)%kfl_react==1) call runend("solve(1)%kfl_react==1")
    !
    ! Get force
    !
    rbfor_sld(:) = 0.0_rp
    if( INOTMASTER ) then
       !!if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
       if ( PLEPP_CPLNG%n_ji > 0_ip ) then
          rbfor_sld(1:ndime) = sum(CNT_CPLNG%var_ji(1:ndime,PLEPP_CPLNG%interior_list_j), dim=2) / real(PLEPP_CPLNG%n_ji,rp)
       end if
    end if
    !
    ! Parall service to all MPIs
    !
    dummax = rbfor_sld(1)
    dummin = rbfor_sld(1)
    call PAR_MIN(dummin, 'IN MY CODE')
    call PAR_MAX(dummax, 'IN MY CODE')
    if( dummin == 0.0_rp ) rbfor_sld(1) = dummax
    if( dummax == 0.0_rp ) rbfor_sld(1) = dummin
    dummax = rbfor_sld(2)
    dummin = rbfor_sld(2)
    call PAR_MIN(dummin, 'IN MY CODE')
    call PAR_MAX(dummax, 'IN MY CODE')
    if( dummin == 0.0_rp ) rbfor_sld(2) = dummax
    if( dummax == 0.0_rp ) rbfor_sld(2) = dummin
    dummax = rbfor_sld(3)
    dummin = rbfor_sld(3)
    call PAR_MIN(dummin, 'IN MY CODE')
    call PAR_MAX(dummax, 'IN MY CODE')
    if( dummin == 0.0_rp ) rbfor_sld(3) = dummax
    if( dummax == 0.0_rp ) rbfor_sld(3) = dummin
    !
    ! Live info
    !
    call messages_live('<- RECV '//trim('FORCE')//' AS A TARGET FOR COUPLING '//trim(intost(1_ip)))

  end subroutine sld_rbo_bilateral_recv

  !-----------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details Impose friction
  !-----------------------------------------------------------------------
  
!!!  subroutine commdom_sld_impose_friction()
!!!!    use def_domain,          only : kfl_codno, lpoty
!!!!    use def_solidz,          only : contactbou_sld, jacrot_du_dq_sld, jacrot_dq_du_sld, ndofn_sld, rhsid1st_it
!!!!    use def_master,          only : rhsid, itinn
!!!!
!!!!    implicit none
!!!!
!!!!    integer(ip)                :: ipoin, nrecv, idx
!!!!    real(rp)                   :: mu, fric_force
!!!!    character(5)               :: legend
!!!!    integer(ip)                :: ievat,jevat,ibopo
!!!!    real(rp)                   :: rhsid_TMP(npoin*ndofn_sld),dummy_matrix(ndofn_sld,ndofn_sld)
!!!!
!!!!    rhsid_TMP = rhsid
!!!!    dummy_matrix = 0.0_rp
!!!!    if ( INOTMASTER ) then
!!!!       do ipoin=1,npoin
!!!!          ibopo = lpoty(ipoin)
!!!!          if ( ibopo > 0 ) then
!!!!             if( kfl_fixno_sld(1,ipoin) == 2 .or. kfl_fixno_sld(1,ipoin) == 3 ) then
!!!!                ievat = (ipoin-1)*ndofn_sld + 1
!!!!                jevat = (ipoin-1)*ndofn_sld + ndofn_sld
!!!!                call sld_rotsys(2_ip, &
!!!!                     1_ip,1_ip,ndofn_sld,ndofn_sld,&
!!!!                     dummy_matrix,rhsid_TMP(ievat:jevat),&
!!!!                     jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
!!!!             end if
!!!!          end if
!!!!       end do
!!!!    end if
!!!!
!!!!    mu = 0.1_rp
!!!!
!!!!    if ( CNT_CPLNG%current_code == CNT_CPLNG%code_j ) then
!!!!       if ( INOTMASTER ) then
!!!!
!!!!          if (.not. associated(rhsid1st_it)) then
!!!!             allocate(rhsid1st_it(npoin*ndofn_sld))
!!!!          end if
!!!!          if(itinn(modul) == 1_ip) then
!!!!             rhsid1st_it = rhsid
!!!!          end if
!!!!
!!!!          nrecv = PLEPP_CPLNG%n_ji
!!!!          do ipoin = 1,nrecv
!!!!             idx = PLEPP_CPLNG%interior_list_j(ipoin)
!!!!             if( any(kfl_codno(:,PLEPP_CPLNG%interior_list_j(ipoin)) == contactbou_sld) ) then
!!!!                if(ndime == 2_ip) then
!!!!                   if(rhsid(2*idx) /= 0_rp) then
!!!!                      !jacrot_du_dq_sld(1,1,PLEPP_CPLNG%interior_list_j(ipoin)) --> nx
!!!!                      !jacrot_du_dq_sld(2,1,PLEPP_CPLNG%interior_list_j(ipoin)) --> ny
!!!!                      !jacrot_du_dq_sld(1,2,PLEPP_CPLNG%interior_list_j(ipoin)) --> tx
!!!!                      !jacrot_du_dq_sld(2,2,PLEPP_CPLNG%interior_list_j(ipoin)) --> ty
!!!!                      if(abs(rhsid(2*idx)) > mu*abs(rhsid(2*idx-1))) then !slip
!!!!                         rhsid(2*idx  ) = rhsid(2*idx ) - mu*abs(rhsid(2*idx-1))*(rhsid(2*idx)/abs(rhsid(2*idx)))
!!!!                         legend = "slip_"
!!!!                         fric_force = - mu*abs(rhsid(2*idx-1))*(rhsid(2*idx)/abs(rhsid(2*idx)))
!!!!                      else !stick
!!!!                         kfl_fixno_sld(2,PLEPP_CPLNG%interior_list_j(ipoin)) = 3_ip
!!!!                         legend = "stick"
!!!!                         fric_force = 0.0_rp
!!!!                      end if
!!!!                      !if( (0.8199999_rp < coord(1,idx)) .and. (coord(1,idx) < 0.82010000_rp) ) then
!!!!                      !write(789,*) coord(1,idx), rhsid1st_it(2*idx-1), rhsid1st_it(2*idx), cutim, ittim
!!!!                      !end if
!!!!                   end if
!!!!                else if(ndime == 3_ip) then
!!!!                   !rhsid(3*idx-2)
!!!!                   !rhsid(3*idx-1)
!!!!                   !rhsid(3*idx  )
!!!!                end if
!!!!             end if
!!!!          end do
!!!!
!!!!       end if
!!!!    end if
!!!  end subroutine commdom_sld_impose_friction

#endif

end module mod_sld_pdn_contact

