!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_vect_assembly_elfem.f90
!> @author  Adria Quintanas-Corominas
!> @date    November 2021
!> @brief   Vectorised assembly for linear elements 
!------------------------------------------------------------------

module mod_sld_vect_assembly_elfem

#include "def_solidz_vect_dim.inc" 
   use def_kintyp, only: ip, rp, lg
   use def_domain, only: ndime
   use def_domain, only: mnode
   use def_solidz, only: ndofn_sld 
   use mod_sld_vect_csm 
   use mod_sld_vect_maths

   implicit none

   integer(ip), protected           :: lelem_loc(DVS)
   integer(ip), protected           :: lmate_loc(DVS)
   integer(ip), protected           :: lelch_loc(DVS)
   logical(lg), protected           :: lmask_loc(DVS)
   integer(ip), protected,  pointer :: lnods_loc(:,:)
   real(rp),    protected,  pointer :: vdsde_sld(:,:,:,:,:)

   public                           :: sld_vect_element_operations_ELFEM_impl
   public                           :: sld_vect_element_operations_ELFEM_expl
   
 contains
   
     !-----------------------------------------------------------------------
     !> 
     !> @author  aquintanas 
     !> @date    November 2021
     !> @brief   Vectorized element assembly for Implicit schemes
     !> @details Vectorized element assembly for Implicit schemes 
     !> 
     !-----------------------------------------------------------------------
     
      subroutine sld_vect_element_operations_ELFEM_impl(itask,pnode,pgaus,list_elements)
         use def_elmtyp,                 only : BAR3D, SHELL, ELINT, ELFEM
         use def_domain,                 only : lnods, ltype, llapl, lorde, lelch, lmate
         use def_domain,                 only : ntens
         use def_domain,                 only : hnatu
         use def_domain,                 only : elmar
         use def_master,                 only : solve
         use mod_element_integration,    only : element_shape_function_derivatives_jacobian
         use def_solidz,                 only : lawst_sld   
         use def_solidz,                 only : kfl_conta_sld, kfl_local_sld
         use def_solidz,                 only : kfl_timei_sld, SLD_DYNAMIC_PROBLEM 
         !
         ! Input and output variables
         !
         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pnode
         integer(ip), intent(in)    :: pgaus
         integer(ip), intent(in)    :: list_elements(DVS)
         !
         ! Element matrices and vectors
         !
         real(rp)    :: elrhs(DVS,ZDOFU)
         real(rp)    :: elfin(DVS,ZDOFU)
         real(rp)    :: elfex(DVS,ZDOFU)
         real(rp)    :: elfrx(DVS,ZDOFU)
         real(rp)    :: elfco(DVS,ZDOFU)
         real(rp)    :: elmas(DVS,ZNODE,ZNODE)
         real(rp)    :: elmat(DVS,ZDOFU,ZDOFU)
         real(rp)    :: elmuu(DVS,ZNODE)
         !
         ! Gather
         !
         real(rp)    :: elcod(DVS,ZDIME,ZNODE)
         real(rp)    :: eldis(DVS,ZDIME,ZNODE)
         real(rp)    :: elvel(DVS,ZDIME,ZNODE)
         real(rp)    :: elacc(DVS,ZDIME,ZNODE)
         real(rp)    :: elddi(DVS,ZDIME,ZNODE)
         real(rp)    :: elsys(DVS,ZDIME,ZDIME)
         real(rp)    :: elcel(DVS)
         !
         ! Variables at Gauss points
         !
         real(rp)    :: gpvol(DVS,ZGAUS)
         real(rp)    :: gpdet(DVS,ZGAUS)                           
         real(rp)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp)    :: gpidg(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp)    :: gpsha(DVS,ZNODE,ZGAUS)
         real(rp)    :: gpder(DVS,ZDIME,ZNODE,ZGAUS)
         real(rp)    :: gpcar(DVS,ZDIME,ZNODE,ZGAUS)
         real(rp)    :: gphes(DVS,ntens,ZNODE,ZGAUS)
         real(rp)    :: gppio(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp)    :: gptmo(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
         !
         ! Indices and dimensions
         !
         integer(ip) :: ielem,imate,iv
         integer(ip) :: pelty, pmate, plapl, porde
         !
         ! Element characteristics
         !
         real(rp)    :: hleng(DVS,ZDIME)
         real(rp)    :: rorde
   
         ielem = list_elements(1)
         imate = lmate(ielem)
         pelty = ltype(ielem)
         plapl = llapl(ZELTY)
         porde = lorde(ZELTY)
         pmate = lmate(ielem)

         !--------------------------------------------------------------------
         !
         ! Gather: global to local
         !
         !--------------------------------------------------------------------

         allocate(lnods_loc(DVS,ZNODE))

         do iv = 1,DVS
            ielem = abs(list_elements(iv))
            if( ielem /= 0 ) then
               lelem_loc(iv)         = list_elements(iv)
               lnods_loc(iv,1:ZNODE) = lnods(1:ZNODE,ielem)
               lelch_loc(iv)         = lelch(ielem)
               lmate_loc(iv)         = lmate(ielem)
               lmask_loc(iv)         = .true.
            else
               lelem_loc(iv)         = list_elements(1)
               lnods_loc(iv,1:ZNODE) = lnods(1:ZNODE,1)
               lelch_loc(iv)         = lelch(1)
               lmate_loc(iv)         = lmate(1)
               lmask_loc(iv)         = .false.
            end if
         end do

         !--------------------------------------------------------------------
         !
         ! Get elemental residual (forces) and matrix
         !
         !--------------------------------------------------------------------

         elrhs(:,:)     = 0.0_rp
         elfin(:,:)     = 0.0_rp
         elfex(:,:)     = 0.0_rp
         elfrx(:,:)     = 0.0_rp
         elfco(:,:)     = 0.0_rp
         elmuu(:,:)     = 0.0_rp
         elmas(:,:,:)   = 0.0_rp
         elmat(:,:,:)   = 0.0_rp
         
         ! TODO: AQC : Assuming that if elements are packed according to material, they are also packed
         !             acording to lelch
         elem_selector: if( lelch_loc(1) == ELFEM .and. lawst_sld(imate) /= 200_ip )then

            ! Gather some variables
            call sld_vect_gather_nodal(pnode,lnods_loc,eldis,elddi,elvel,elacc,elcod)
  
            ! Get interpolation scheme at GP
            call element_shape_function_derivatives_jacobian( &
               ZNODE, ZGAUS, plapl, &
               elmar(ZELTY) % weigp, &
               elmar(ZELTY) % shape, &
               elmar(ZELTY) % deriv, &
               elmar(ZELTY) % heslo, &
               elcod, gpvol, gpsha, gpder, gpcar, gphes, &
               list_elements=lelem_loc )
            
            if( plapl == 0_ip ) gphes = 0.0_rp

            ! Compute element length
            call sld_vect_element_length( &
               ZNODE, hnatu(ZELTY), elcod, &
               elmar(ZELTY) % dercg, &
               hleng )

            ! Compute Characteristic element length 
            rorde =  1.0_rp/real(porde,rp)
            elcel(1:DVS) = hleng(1:DVS,ZDIME) * rorde

            ! Interpolate variables at GP
            call sld_vect_interpolate( &
               ZNODE, ZGAUS, &
               eldis, &
               gpcar, gpgdi, gpidg, gpdet )

            ! TODO: Improve this
            call sld_vect_gather_elemental_csys(lelem_loc,elsys)

            ! Evaluate material law (stress model)
            call sld_builtin_materials_vect( &
            itask, pgaus, pmate, &
            elsys, &
            gpgdi, gpidg, gpdet, &
            gppio, gptmo )

            ! Compute elemental RHS and AMAT
            if( itask == 2_ip ) then
               call sld_vect_elmmat_RHS( &
                    pgaus, pmate, pnode, &
                    gpvol, gpsha, gpcar, gppio, &
                    elfin, elfex, elrhs, &
                    elmuu)
               call sld_vect_elmmat_AMAT(&
                    pgaus,pnode,&
                    gpvol,gpcar,gptmo,&
                    elmat)
            endif
            
            ! Add inertial contribution
            if( itask == 2_ip .and. kfl_timei_sld == SLD_DYNAMIC_PROBLEM )then
               call sld_vect_add_inertial_contribution( &
                  pgaus, pmate, pnode, gpvol, gpsha, elacc, elrhs, elmat ) 
            endif

            ! From global to local coordinate system 
            if( kfl_conta_sld /= 0_ip .or. kfl_local_sld /= 0_ip  )then 
               if( itask == 2_ip )then
                  call sld_vect_change_csys_from_global_to_local_elrhs_elmat( &
                     pnode, lmask_loc, lnods_loc, &
                     elrhs, elmat )
               else
                  call sld_vect_change_csys_from_global_to_local( &
                     pnode, lmask_loc, lnods_loc, &
                     elrhs )
               endif  
            endif

            ! Get contact forces
            if( kfl_conta_sld /= 0_ip )then
               call sld_vect_get_contact_force_at_ELEM( &
                  pnode, lmask_loc, lnods_loc, &
                  elrhs, elfco )
            endif

            ! Get reaction forces
            call sld_vect_get_reaction_force_at_ELEM( &
               pnode, lmask_loc, lnods_loc, &
               elrhs, elfrx )

            ! Apply elemental boundary conditions
            if( solve(1) % kfl_iffix == 0 )then
               call sld_vect_apply_essential_BC_at_ELEM( &
                  pnode, lmask_loc, lnods_loc, &
                  elrhs, elmat )
            endif
            
            ! Scatter some variables
            call sld_vect_scatter( lelem_loc, elcel )

            ! Assemble to global RHS and AMAT
            call sld_vect_assemble_global_RHS_AMAT( &
               itask, pnode, lmask_loc, lelem_loc, lnods_loc, &
               elfin, elfex, elfrx, elfco, elrhs, elmuu, elmat )

         end if elem_selector

         deallocate(lnods_Loc)

      end subroutine sld_vect_element_operations_ELFEM_impl

      !-----------------------------------------------------------------------
      !> 
      !> @author  aquintanas and gguillam
      !> @date    2022-11-04
      !> @brief   Vectorized element assembly for Explicit schemes
      !> @details Vectorized element assembly for Explicit schemes
      !> 
      !-----------------------------------------------------------------------
      
      subroutine sld_vect_element_operations_ELFEM_expl(itask,pnode,pgaus,list_elements)

        use def_elmtyp,                 only : BAR3D, SHELL, ELINT, ELFEM
        use def_domain,                 only : lnods, ltype, llapl, lelch, lmate
        use def_domain,                 only : ntens
        use def_domain,                 only : elmar
        use mod_element_integration,    only : element_shape_function_derivatives_jacobian 
        use def_solidz,                 only : kfl_conta_sld, kfl_local_sld
        !
        ! Input and output variables
        !
        integer(ip), intent(in)    :: itask
        integer(ip), intent(in)    :: pnode
        integer(ip), intent(in)    :: pgaus
        integer(ip), intent(in)    :: list_elements(DVS)
        !
        ! Element matrices and vectors
        !
        real(rp)    :: elrhs(DVS,ZDOFU)
        real(rp)    :: elfin(DVS,ZDOFU)
        real(rp)    :: elfex(DVS,ZDOFU)
        real(rp)    :: elfrx(DVS,ZDOFU)
        real(rp)    :: elfco(DVS,ZDOFU)
        real(rp)    :: elmuu(DVS,ZNODE)
        !
        ! Gather
        !
        real(rp)    :: elcod(DVS,ZDIME,ZNODE)
        real(rp)    :: eldis(DVS,ZDIME,ZNODE)
        real(rp)    :: elvel(DVS,ZDIME,ZNODE)
        real(rp)    :: elacc(DVS,ZDIME,ZNODE)
        real(rp)    :: elddi(DVS,ZDIME,ZNODE)
        real(rp)    :: elsys(DVS,ZDIME,ZDIME)
        !
        ! Variables at Gauss points
        !
        real(rp)    :: gpvol(DVS,ZGAUS)
        real(rp)    :: gpdet(DVS,ZGAUS)                           
        real(rp)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
        real(rp)    :: gpidg(DVS,ZDIME,ZDIME,ZGAUS)
        real(rp)    :: gpsha(DVS,ZNODE,ZGAUS)
        real(rp)    :: gpder(DVS,ZDIME,ZNODE,ZGAUS)
        real(rp)    :: gpcar(DVS,ZDIME,ZNODE,ZGAUS)
        real(rp)    :: gphes(DVS,ntens,ZNODE,ZGAUS)
        real(rp)    :: gppio(DVS,ZDIME,ZDIME,ZGAUS)
        real(rp)    :: gptmo(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
        !
        ! Indices and dimensions
        !
        integer(ip) :: ielem, iv
        integer(ip) :: pelty, pmate, plapl

        ielem = list_elements(1)
        pelty = ltype(ielem)
        plapl = llapl(ZELTY)
        pmate = lmate(ielem)

        !--------------------------------------------------------------------
        !
        ! Gather: global to local
        !
        !--------------------------------------------------------------------

        allocate(lnods_loc(DVS,ZNODE))

        do iv = 1,DVS
           ielem = abs(list_elements(iv))
           if( ielem /= 0 ) then
              lelem_loc(iv)         = list_elements(iv)
              lnods_loc(iv,1:ZNODE) = lnods(1:ZNODE,ielem)
              lelch_loc(iv)         = lelch(ielem)
              lmate_loc(iv)         = lmate(ielem)
              lmask_loc(iv)         = .true.
           else
              lelem_loc(iv)         = list_elements(1)
              lnods_loc(iv,1:ZNODE) = lnods(1:ZNODE,1)
              lelch_loc(iv)         = lelch(1)
              lmate_loc(iv)         = lmate(1)
              lmask_loc(iv)         = .false.
           end if
        end do

        !--------------------------------------------------------------------
        !
        ! Get elemental residual (forces)
        !
        !--------------------------------------------------------------------

        elrhs(:,:)     = 0.0_rp
        elfin(:,:)     = 0.0_rp
        elfex(:,:)     = 0.0_rp
        elfrx(:,:)     = 0.0_rp
        elfco(:,:)     = 0.0_rp
        elmuu(:,:)     = 0.0_rp

        ! Gather some variables
        call sld_vect_gather_nodal(pnode,lnods_loc,eldis,elddi,elvel,elacc,elcod)

        ! Get interpolation scheme at GP
        call element_shape_function_derivatives_jacobian( &
             ZNODE, ZGAUS, plapl, &
             elmar(ZELTY) % weigp, &
             elmar(ZELTY) % shape, &
             elmar(ZELTY) % deriv, &
             elmar(ZELTY) % heslo, &
             elcod, gpvol, gpsha, gpder, gpcar, gphes, &
             list_elements=lelem_loc )

        if( plapl == 0_ip ) gphes = 0.0_rp

        ! Interpolate variables at GP
        call sld_vect_interpolate( &
             ZNODE, ZGAUS, &
             eldis, &
             gpcar, gpgdi, gpidg, gpdet )

        ! TODO: Improve this
        call sld_vect_gather_elemental_csys(lelem_loc,elsys)

        ! Evaluate material law (stress model)
        call sld_builtin_materials_vect( &
             itask, pgaus, pmate, &
             elsys, &
             gpgdi, gpidg, gpdet, &
             gppio, gptmo )

        ! Compute elemental RHS
        call sld_vect_elmmat_RHS( &
             pgaus, pmate, pnode, &
             gpvol, gpsha, gpcar, gppio, &
             elfin, elfex, elrhs, &
             elmuu)

        ! From global to local coordinate system 
        if( kfl_conta_sld /= 0_ip .or. kfl_local_sld /= 0_ip  )then 
           call sld_vect_change_csys_from_global_to_local( &
                pnode, lmask_loc, lnods_loc, &
                elrhs )
        endif

        ! Get contact forces
        if( kfl_conta_sld /= 0_ip )then
           call sld_vect_get_contact_force_at_ELEM( &
                pnode, lmask_loc, lnods_loc, &
                elrhs, elfco )
        endif

        ! Get reaction forces
        call sld_vect_get_reaction_force_at_ELEM( &
             pnode, lmask_loc, lnods_loc, &
             elrhs, elfrx )

        ! Assemble to global RHS
        call sld_vect_assemble_global_RHS( &
             itask, pnode, lmask_loc, lnods_loc, &
             elfin, elfex, elfrx, elfco, elrhs, elmuu )

        deallocate(lnods_Loc)

      end subroutine sld_vect_element_operations_ELFEM_expl

      
      pure subroutine sld_vect_gather_nodal(pnode,lnods_loc,eldis,elddi,elvel,elacc,elcod)
         use def_master, only : displ, ITER_K
         use def_domain, only : coord
         use def_solidz, only : veloc_sld, accel_sld, ddisp_sld
     
         integer(ip), intent(in)  :: pnode
         integer(ip), intent(in)  :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(out) :: elcod(DVS,ZDIME,ZNODE)
         real(rp),    intent(out) :: eldis(DVS,ZDIME,ZNODE)
         real(rp),    intent(out) :: elddi(DVS,ZDIME,ZNODE)
         real(rp),    intent(out) :: elvel(DVS,ZDIME,ZNODE)
         real(rp),    intent(out) :: elacc(DVS,ZDIME,ZNODE)
         integer(ip)              :: inode, ipoin, iv
     
         do iv = 1,DVS
            do inode = 1,ZNODE
               ipoin = lnods_loc(iv,inode)
               elcod(iv,1:ZDIME,inode) = coord(1:ZDIME,ipoin)
               eldis(iv,1:ZDIME,inode) = displ(1:ZDIME,ipoin,ITER_K)
               elddi(iv,1:ZDIME,inode) = ddisp_sld(1:ZDIME,ipoin,ITER_K)
               elvel(iv,1:ZDIME,inode) = veloc_sld(1:ZDIME,ipoin,ITER_K)
               elacc(iv,1:ZDIME,inode) = accel_sld(1:ZDIME,ipoin,ITER_K)
            end do
         end do
     
      end subroutine sld_vect_gather_nodal


      pure subroutine sld_vect_gather_elemental_csys(lelem_loc,elsys)
         use def_solidz, only : kfl_fiber_sld, axis1_sld, axis2_sld, axis3_sld 

         integer(ip), intent(in)  :: lelem_loc(DVS)
         real(rp),    intent(out) :: elsys(DVS,ZDIME,ZDIME)
         integer(ip)              :: iv, ielem 

         ! More information in mod_sld_csys
         if( kfl_fiber_sld >= 4_ip .and. kfl_fiber_sld <= 7_ip )then
            do iv = 1,DVS
               ielem = lelem_loc(iv)
               elsys(iv,1,1:ZDIME) = real(axis1_sld(1:ZDIME,ielem),rp)
               elsys(iv,2,1:ZDIME) = real(axis2_sld(1:ZDIME,ielem),rp)
               elsys(iv,3,1:ZDIME) = real(axis3_sld(1:ZDIME,ielem),rp)
            end do
         endif

      end subroutine sld_vect_gather_elemental_csys


      subroutine sld_vect_scatter(lelem_loc,elcel)
         use def_solidz, only : celen_sld

         integer(ip), intent(in)  :: lelem_loc(DVS)
         real(rp),    intent(in)  :: elcel(DVS)
         integer(ip)              :: iv
     
         do iv = 1,DVS
            celen_sld(lelem_loc(iv)) = elcel(iv)
         end do
     
      end subroutine sld_vect_scatter


      pure subroutine sld_vect_interpolate_pure(pnode,pgaus,eldis,&
         gpcar,gpgdi,gpidg,gpdet)

         integer(ip), intent(in)    :: pnode
         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: eldis(DVS,ZDIME,ZNODE)
         real(rp),    intent(in)    :: gpcar(DVS,ZDIME,ZNODE,ZGAUS)
         real(rp),    intent(out)   :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpidg(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpdet(DVS,ZGAUS)
         integer(ip)                :: ig,ii

         ! Deformation Gradient tensor related variables
         do ig = 1, ZGAUS
            ! F = I + dU/dX
            gpgdi(:,:,:,ig) = &
               vmath_MxT(eldis(:,:,:),gpcar(:,:,:,ig))
            do ii = 1, ZDIME
               gpgdi(1:DVS,ii,ii,ig) = gpgdi(1:DVS,ii,ii,ig) + 1.0_rp
            enddo
            ! F^-1
            gpidg(:,:,:,ig) = vmath_INV(ZDIME,gpgdi(:,:,:,ig))
            ! J = det(F)
            gpdet(:,ig) = vmath_DET(ZDIME,gpgdi(:,:,:,ig))
         enddo
         
      end subroutine sld_vect_interpolate_pure


      pure subroutine sld_vect_interpolate(pnode,pgaus,eldis,&
         gpcar,gpgdi,gpidg,gpdet)

         integer(ip), intent(in)    :: pnode
         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: eldis(DVS,ZDIME,ZNODE)
         real(rp),    intent(in)    :: gpcar(DVS,ZDIME,ZNODE,ZGAUS)
         real(rp),    intent(out)   :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpidg(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpdet(DVS,ZGAUS)
         integer(ip)                :: ig,ii,jj,kk

         ! Deformation Gradient tensor related variables
         do ig = 1, ZGAUS
            ! F = I + dU/dX
            gpgdi(:,:,:,ig) = 0.0_rp
            do jj = 1, ZDIME
               do ii = 1, ZDIME
                  do kk = 1, ZNODE
                     gpgdi(1:DVS,ii,jj,ig) = gpgdi(1:DVS,ii,jj,ig) + eldis(1:DVS,ii,kk) * gpcar(1:DVS,jj,kk,ig)
                  enddo
               enddo
            enddo
            do ii = 1, ZDIME
               gpgdi(1:DVS,ii,ii,ig) = gpgdi(1:DVS,ii,ii,ig) + 1.0_rp
            enddo
            ! F^-1
            gpidg(:,:,:,ig) = vmath_INV(ZDIME,gpgdi(:,:,:,ig))
            ! J = det(F)
            gpdet(:,ig) = vmath_DET(ZDIME,gpgdi(:,:,:,ig))
         enddo

      end subroutine sld_vect_interpolate

      
      pure subroutine sld_vect_element_length(pnode,hnatu,elcod,gpder,hleng)

         integer(ip), intent(in)    :: pnode
         real(rp),    intent(in)    :: hnatu
         real(rp),    intent(in)    :: elcod(DVS,ZDIME,ZNODE)
         real(rp),    intent(in)    :: gpder(ZDIME,ZNODE)
         real(rp),    intent(out)   :: hleng(DVS,ZDIME)
         real(rp)                   :: enor0(DVS)
         real(rp)                   :: xjacm(DVS,ZDIME,ZDIME)
         real(rp)                   :: xjacd(DVS)
         real(rp)                   :: xjaci(DVS,ZDIME,ZDIME)
         integer(ip)                :: ii, jj, kk, iv
         real(rp)                   :: h_tmp

         ! Jacobian matrix (J)
         do jj = 1, ZDIME
            do ii = 1, ZDIME
               xjacm(1:DVS,ii,jj) = 0.0_rp
            enddo
         enddo
         do jj = 1, ZDIME
            do ii = 1, ZDIME
               do kk = 1, ZNODE
                  xjacm(1:DVS,ii,jj) = xjacm(1:DVS,ii,jj) + elcod(1:DVS,ii,kk) * gpder(jj,kk)
               enddo
            enddo
         enddo

         ! det(J) and J^-1
         xjacd = vmath_DET(ZDIME,xjacm)
         xjaci = vmath_INV(ZDIME,xjacm)

         ! Length according to dimension
         dimension_selector: if( ZDIME == 2_ip )then
      
            enor0(1:DVS)   =  xjaci(1:DVS,1,1) * xjaci(1:DVS,1,1) &
              + xjaci(1:DVS,1,2) * xjaci(1:DVS,1,2)
            hleng(1:DVS,1) =  hnatu/sqrt(enor0(1:DVS))

            enor0(1:DVS)   =  xjaci(1:DVS,2,1) * xjaci(1:DVS,2,1) &
               + xjaci(1:DVS,2,2) * xjaci(1:DVS,2,2)
            hleng(1:DVS,2) =  hnatu/sqrt(enor0(1:DVS))
         
            do iv = 1,DVS
               if( hleng(iv,2) > hleng(iv,1) )then
                  h_tmp    = hleng(iv,2)
                  hleng(iv,2) = hleng(iv,1)
                  hleng(iv,1) = h_tmp
               end if
            end do
            
         else if( ZDIME == 3_ip )then
        
            enor0(1:DVS)     = xjaci(1:DVS,1,1) * xjaci(1:DVS,1,1) &
               + xjaci(1:DVS,1,2) * xjaci(1:DVS,1,2) &
               + xjaci(1:DVS,1,3) * xjaci(1:DVS,1,3)
            hleng(1:DVS,1)   = hnatu/sqrt(enor0(1:DVS))
            !hleng(1:DVS,1)   = hnatu/max(sqrt(enor0(1:DVS)),epsilon(1.0_rp))
            enor0(1:DVS)     = xjaci(1:DVS,2,1) * xjaci(1:DVS,2,1) &
               + xjaci(1:DVS,2,2) * xjaci(1:DVS,2,2) &
               + xjaci(1:DVS,2,3) * xjaci(1:DVS,2,3)
            hleng(1:DVS,2)   = hnatu/sqrt(enor0(1:DVS))
            !hleng(1:DVS,2)   = hnatu/max(sqrt(enor0(1:DVS)),epsilon(1.0_rp))
            enor0(1:DVS)     = xjaci(1:DVS,3,1) * xjaci(1:DVS,3,1) &
               + xjaci(1:DVS,3,2) * xjaci(1:DVS,3,2) &
               + xjaci(1:DVS,3,3) * xjaci(1:DVS,3,3)
            hleng(1:DVS,3)  = hnatu/sqrt(enor0(1:DVS))
            !hleng(1:DVS,3)   = hnatu/max(sqrt(enor0(1:DVS)),epsilon(1.0_rp))
            do iv = 1, DVS
        
               if( hleng(iv,2) > hleng(iv,1) )then
                  h_tmp       = hleng(iv,2)
                  hleng(iv,2) = hleng(iv,1)
                  hleng(iv,1) = h_tmp
               end if
               if( hleng(iv,3) > hleng(iv,1) )then
                  h_tmp    = hleng(iv,3)
                  hleng(iv,3) = hleng(iv,1)
                  hleng(iv,1) = h_tmp
               end if
               if( hleng(iv,3) > hleng(iv,2) )then
                  h_tmp       = hleng(iv,3)
                  hleng(iv,3) = hleng(iv,2)
                  hleng(iv,2) = h_tmp
               end if
            end do 
   
         end if dimension_selector

      end subroutine sld_vect_element_length

      !-----------------------------------------------------------------------
      !> 
      !> @author  aquintanas
      !> @date    2022-11-15
      !> @brief   Vectorised built-in materials 
      !> @details Vectorised built-in materials
      !>
      !>          PKN : Nominal stress tensor  PKN = P^T. In Bely's book is P
      !<          PK1 : 1st Piola Kirchoff stress tensor. 
      !<          PK2 : 2nd Piola Kirchoff stress tensor. In Bely's book is S
      !<          SIG : Cauchy stress tensor. In Bely's book is sigma
      !>
      !-----------------------------------------------------------------------
      
      subroutine sld_builtin_materials_vect(itask,pgaus,pmate,&
         mcsys,gpgdi,gpidg,gpdet,gppio,gptmo)
         use def_solidz,                    only : lawst_sld
         use def_solidz,                    only : kfl_strai_sld
         use def_solidz,                    only : SLD_GREEN,SLD_INFINITESIMAL 
         use mod_sld_vect_stress_models
         use mod_sld_vect_stress_model_154, only : nsdv154
         use mod_sld_vect_stress_model_154, only : sld_vect_stress_model_154
         use mod_sld_vect_stress_model_154, only : sld_vect_stress_model_154_getSDV
         use mod_sld_vect_stress_model_154, only : sld_vect_stress_model_154_setSDV

         external                   :: runend
         
         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         integer(ip), intent(in)    :: pmate
         real(rp),    intent(in)    :: mcsys(DVS,ZDIME,ZDIME,*)     ! Material coordinate system
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(in)    :: gpidg(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(in)    :: gpdet(DVS,ZGAUS) 
         real(rp),    intent(out)   :: gppio(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gptmo(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
         real(rp)                   :: gpcau(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp)                   :: gpgre(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp)                   :: gpstr(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp)                   :: gpdds(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
         real(rp),    pointer       :: gpsdv(:,:,:)
         logical(lg)                :: from_PK2_to_PK1
         logical(lg)                :: from_CAU_to_PK1
         !
         ! Initializations
         !
         from_PK2_to_PK1 = .false.    
         from_CAU_to_PK1 = .false.    
         gppio(1:DVS,:,:,:)     = 0.0_rp
         gptmo(1:DVS,:,:,:,:,:) = 0.0_rp
         !
         ! Strain 
         if(      kfl_strai_sld == SLD_GREEN         ) then
            from_PK2_to_PK1 = .true.
            call sld_vcsm_get_lagrange_strain_tensors_at_GP(ZGAUS, gpgdi, gpcau, gpgre)
         else if( kfl_strai_sld == SLD_INFINITESIMAL ) then
            from_CAU_to_PK1 = .true.
            call sld_vcsm_get_infinitesimal_strain_tensors_at_GP(ZGAUS, gpgdi, gpgre)
         end if
               
         select case( lawst_sld(pmate) )

         case( SM_ISOLIN )
            ! Isotropic linear Saint-Venant model
            call sld_vect_stress_model_100(itask,pmate,ZGAUS,gpgre,gpstr,gpdds)

         case( SM_ORTLIN )
            ! Othotropic linear Saint-Venant model
            call sld_vect_stress_model_151b(itask,pmate,ZGAUS,mcsys(:,:,:,1),gpgre,gpstr,gpdds)

         case( SM_BESSA )
            ! Transversally isotropic model with damage for composites (Maimi et al. 2007)
            allocate(gpsdv(DVS,nsdv154,ZGAUS))
            call sld_vect_stress_model_154_getSDV(lelem_loc, ZGAUS, gpsdv)
            call sld_vect_stress_model_154(lelem_loc, pmate, ZGAUS, mcsys(:,:,:,1), gpgre, gpstr, gpsdv)
            call sld_vect_stress_model_154_setSDV(lelem_loc, ZGAUS, gpsdv)
            deallocate(gpsdv)

         case default
            ! Material not implemented
            call runend('SLD_BUILTIN_MATERIALS_VECT: MATERIAL MODEL NOT IMPLEMENTED')

         end select

         ! If PK2 and dSdE are provided compute PK1 and dPdF
         if( from_PK2_to_PK1 )then
            ! Transform PK2 -> PK1
            call sld_vcsm_transform_PK2_to_PK1_at_GP(ZGAUS,gpgdi,gpstr,gppio)
            if( itask == 2_ip )then
               ! Transform dSdE -> dPdF
               call sld_vcsm_transform_dPK2dE_to_dPK1dF_at_GP(ZGAUS,gpgdi,gpstr,gpdds,gptmo)
            endif
         endif
         
         ! If SIG is provided compute PK1
         if( from_CAU_to_PK1 ) then
            ! Transform SIG (Cauchy) -> PK1  
            !call sld_vcsm_transform_sigma_to_PK1_at_GP(ZGAUS,gpidg,gpdet,gpstr,gppio)
            gppio(1:DVS,:,:,:) = gpstr(1:DVS,:,:,:)
         endif
         
      end subroutine sld_builtin_materials_vect


      pure subroutine sld_vect_elmmat_AMAT(pgaus,pnode,&
            gpvol,gpcar,gptmo,elmat)
  
         integer(ip), intent(in)  :: pgaus
         integer(ip), intent(in)  :: pnode
         real(rp),    intent(in)  :: gpvol(DVS,ZGAUS)                         !> Weight of gauss point
         real(rp),    intent(in)  :: gpcar(DVS,ZDIME,ZNODE,ZGAUS)             !> Cartesian derivatives of the shape functions
         real(rp),    intent(in)  :: gptmo(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS) !> First elasticity tensor at the gauss point
         real(rp),    intent(out) :: elmat(DVS,ZDOFU,ZDOFU)                   !> Elemental AMAT
         integer(ip)              :: pn,jn,ii,jj,kk,ll,iv,jv,ig
         !
         ! Initialization
         !
         elmat = 0.0_rp
         !
         ! Adding contributions to ELMAT
         !         
         ! Compute elemental matrix
         do ig = 1,ZGAUS
            do pn = 1,ZNODE
               do ii = 1,ZDIME
                  iv = ( pn - 1_ip ) * ZDOFN + ii
                  do jn = 1,ZNODE
                     do kk = 1,ZDIME
                        jv = ( jn - 1_ip ) * ZDOFN + kk
                        do jj = 1,ZDIME
                           do ll = 1,ZDIME
                              ! Internal force linearitzation
                              !    K_iakb = int_\Omega (dP_iJ / dF_kL * dNa/dxJ * dNb/dxL * (1 + beta_rayleigh / dt) d\Omega
                              elmat(1:DVS,iv,jv) = elmat(1:DVS,iv,jv) + gpvol(1:DVS,ig) * &
                                   ( gptmo(1:DVS,ii,jj,kk,ll,ig) * gpcar(1:DVS,jj,pn,ig) * gpcar(1:DVS,ll,jn,ig) )
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
     
      end subroutine sld_vect_elmmat_AMAT

      pure subroutine sld_vect_elmmat_RHS(pgaus,pmate,pnode,&
            gpvol,gpsha,gpcar,gppio,&
            elfin,elfex,elrhs,elmuu)
         use def_solidz, only : densi_sld,grnor_sld
         use def_solidz, only : gravi_sld,bodyf_sld
   
         integer(ip), intent(in)  :: pgaus
         integer(ip), intent(in)  :: pmate
         integer(ip), intent(in)  :: pnode
         real(rp),    intent(in)  :: gpvol(DVS,ZGAUS)             !> Weight of gauss point
         real(rp),    intent(in)  :: gpsha(DVS,ZNODE,ZGAUS)       !> Shape functions at gauss points
         real(rp),    intent(in)  :: gpcar(DVS,ZDIME,ZNODE,ZGAUS) !> Cartesian derivatives of the shape functions at gauss points
         real(rp),    intent(in)  :: gppio(DVS,ZDIME,ZDIME,ZGAUS) !> First Piola Kirchoff stress tensor at the gauss point
         real(rp),    intent(out) :: elfin(DVS,ZDOFU)             !> Internal force vector
         real(rp),    intent(out) :: elfex(DVS,ZDOFU)             !> External force vector
         real(rp),    intent(out) :: elrhs(DVS,ZDOFU)             !> Elemental RHSID
         real(rp),    intent(out) :: elmuu(DVS,ZNODE)             !> Lumped mass matrix
         integer(ip)              :: pn,jn,ii,jj,iv,ig
         !
         ! Initialization
         !
         elfin = 0.0_rp
         elfex = 0.0_rp
         elrhs = 0.0_rp
         elmuu = 0.0_rp
         !
         ! Adding contributions to ELRHSID
         !
         do ig = 1, ZGAUS
            do pn = 1, ZNODE
               do ii = 1, ZDIME
                  iv = ( pn - 1_ip ) * ZDOFN + ii
                  do jj = 1,ZDIME
                     ! Internal force contribution
                     !    fint_kb = int_\Omega P_kL * dNb/dxL d\Omega
                     elfin(1:DVS,iv) = elfin(1:DVS,iv) + &
                        gppio(1:DVS,ii,jj,ig) * gpcar(1:DVS,jj,pn,ig) * gpvol(1:DVS,ig)
                  enddo
                  ! External force contribution
                  !    fext_kb = int_\Omega b_0*N d\Omega + int_\dOmega tract*N d\S
                  !        where b_0 are bodyforces and  gravity
                  !              tract are traction forces (Neumann)
                  elfex(1:DVS,iv) = elfex(1:DVS,iv) + &
                     gpsha(1:DVS,pn,ig) * densi_sld(1,pmate) * gpvol(1:DVS,ig) * &
                     ( grnor_sld * gravi_sld(ii) + bodyf_sld(ii) )
               enddo
               ! Mass matrix
               ! M_ab = int_\Omega /rho * Na * Nb d\Omega
               do jn = 1,ZNODE
                  elmuu(1:DVS,pn) = elmuu(1:DVS,pn) + &
                     densi_sld(1,pmate) * gpvol(1:DVS,ig) * gpsha(1:DVS,pn,ig) * gpsha(1:DVS,jn,ig) 
               enddo
            enddo
         enddo
         
         elrhs = elfex - elfin
         
      end subroutine sld_vect_elmmat_RHS
      
      pure subroutine sld_vect_add_inertial_contribution( &
            pgaus,pmate,pnode,gpvol,gpsha,elacc,elrhs,elmat)
         use def_master, only : dtime
         use def_solidz, only : densi_sld, tifac_sld

         integer(ip), intent(in)    :: pgaus
         integer(ip), intent(in)    :: pmate
         integer(ip), intent(in)    :: pnode
         real(rp),    intent(in)    :: gpvol(DVS,ZGAUS)        !> Weight of gauss point
         real(rp),    intent(in)    :: gpsha(DVS,ZNODE,ZGAUS)  !> Shape functions at gauss points
         real(rp),    intent(in)    :: elacc(DVS,ZDIME,ZNODE)  !> Elemental accelerations 
         real(rp),    intent(inout) :: elrhs(DVS,ZDOFU)        !> Elemental RHSID
         real(rp),    intent(inout) :: elmat(DVS,ZDOFU,ZDOFU)  !> Elemental AMAT
         real(rp)                   :: elfine(DVS,ZDOFU)
         real(rp)                   :: elmas(DVS,ZDOFU,ZDOFU)
         real(rp)                   :: rho, betanewmark
         integer(ip)                :: ig, aa, bb, iv, jv, ii 

         ! Rho
          rho = densi_sld(1,pmate)
          betanewmark= 1.0_rp / (tifac_sld(1) * dtime * dtime)

         ! Consistent mass matrix
         !   M_ab = int_\Omega /rho * Na * Nb d\Omega
         elmas = 0.0_rp
         do ig = 1, ZGAUS 
            do aa = 1, ZNODE
               do bb = 1, ZNODE
                  elmas(1:DVS,aa,bb) = elmas(1:DVS,aa,bb) + rho * gpvol(1:DVS,ig) * gpsha(1:DVS,aa,ig) * gpsha(1:DVS,bb,ig)
               end do
            end do
         end do

         ! Inertial forces
         !    f_kb += M_ba * accel_ak
         elfine = 0.0_rp
         do ii = 1, ZDIME
            do aa = 1, ZNODE
               iv = ( aa - 1_ip ) * ZDOFN + ii
               do bb = 1, ZNODE
                  elfine(1:DVS,iv) = elfine(1:DVS,iv) + elmas(1:DVS,aa,bb) * elacc(1:DVS,ii,bb)
               end do
            end do
         end do
         
         ! ELRHS
         elrhs = elrhs - elfine

         ! ELMAT: Inertial contribution
         !    A = 1/(beta*dtime^2)*M + dFint/dU
         do ii = 1, ZDIME
            do aa = 1, ZNODE
               iv = ( aa - 1 ) * ZDOFN + ii
               do bb = 1, ZNODE
                  jv = ( bb - 1 ) * ZDOFN + ii
                  elmat(1:DVS,iv,jv) = elmat(1:DVS,iv,jv) + elmas(1:DVS,aa,bb) * betanewmark 
               end do
            end do
         end do

      end subroutine sld_vect_add_inertial_contribution


      subroutine sld_vect_assemble_global_RHS_AMAT(itask,pnode, &
         lmask_loc,lelem_loc,lnods_loc,elfin,elfex,elfrx,elfco, &
         elrhs,elmuu,elmat) 
         use def_solidz, only : kfl_timei_sld, SLD_DYNAMIC_PROBLEM, kfl_conta_sld
         use def_solidz, only : frxid_sld, finte_sld, fexte_sld, fcont_sld, vmass_sld
         use mod_matrix, only : matrix_assemble_element_matrix_to_CSR
         use mod_matrix, only : matrix_assemble_element_RHS
         use mod_solver, only : solver_assemble_element_matrix_vector
         use def_kermod, only : kfl_element_to_csr
         use def_master, only : rhsid, amatr, solve 
         use def_domain, only : r_dom, c_dom, lezdo

         integer(ip), intent(in)  :: itask
         integer(ip), intent(in)  :: pnode
         logical(lg), intent(in)  :: lmask_loc(DVS)
         integer(ip), intent(in)  :: lelem_loc(DVS)
         integer(ip), intent(in)  :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(in)  :: elfin(DVS,ZDOFU)
         real(rp),    intent(in)  :: elfex(DVS,ZDOFU)
         real(rp),    intent(in)  :: elfrx(DVS,ZDOFU)
         real(rp),    intent(in)  :: elfco(DVS,ZDOFU)
         real(rp),    intent(in)  :: elrhs(DVS,ZDOFU)
         real(rp),    intent(in)  :: elmuu(DVS,ZNODE)
         real(rp),    intent(in)  :: elmat(DVS,ZDOFU,ZDOFU)
         integer(ip)              :: iv

         not_postprocess: if( itask == 2_ip )then

            ! Assembly RHS
            do iv = 1, DVS
               if( lmask_loc(iv) ) then
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elrhs(iv,:), rhsid )
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elfrx(iv,:), frxid_sld)
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elfin(iv,:), finte_sld )
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elfex(iv,:), fexte_sld )
               endif
            enddo

            ! Lumped contact force
            if( kfl_conta_sld /= 0_ip )then
               do iv = 1, DVS
                  if( lmask_loc(iv) ) then
                     call matrix_assemble_element_RHS( &
                        ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                        elfco(iv,:), fcont_sld )
                  endif
               enddo
            endif

            ! Assemble LUMPED MASS MATRIX
            if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
               do iv = 1, DVS
                  if( lmask_loc(iv) ) then
                     call matrix_assemble_element_RHS( &
                           1_ip, 1_ip, ZNODE, lnods_loc(iv,:), &
                           elmuu(iv,:), vmass_sld )
                  endif
               enddo
            end if
            
            ! Assemble matrix
            do iv = 1, DVS
               if( lmask_loc(iv) ) then
                  call matrix_assemble_element_matrix_to_CSR( & 
                       kfl_element_to_csr, solve(1) % ndofn, ZNODE, ZDOFU, &
                       lelem_loc(iv), lnods_loc(iv,:), &
                       elmat(iv,:,:), r_dom, c_dom, amatr, lezdo)
               endif
            enddo

         endif not_postprocess

      end subroutine sld_vect_assemble_global_RHS_AMAT

      subroutine sld_vect_assemble_global_RHS(itask,pnode, &
         lmask_loc,lnods_loc,elfin,elfex,elfrx,elfco, &
         elrhs,elmuu) 
         use def_solidz, only : kfl_conta_sld
         use def_solidz, only : frxid_sld, finte_sld, fexte_sld, fcont_sld, vmass_sld
         use mod_matrix, only : matrix_assemble_element_RHS
         use def_master, only : rhsid

         integer(ip), intent(in)  :: itask
         integer(ip), intent(in)  :: pnode
         logical(lg), intent(in)  :: lmask_loc(DVS)
         integer(ip), intent(in)  :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(in)  :: elfin(DVS,ZDOFU)
         real(rp),    intent(in)  :: elfex(DVS,ZDOFU)
         real(rp),    intent(in)  :: elfrx(DVS,ZDOFU)
         real(rp),    intent(in)  :: elfco(DVS,ZDOFU)
         real(rp),    intent(in)  :: elrhs(DVS,ZDOFU)
         real(rp),    intent(in)  :: elmuu(DVS,ZNODE)
         integer(ip)              :: iv

         not_postprocess: if( itask == 1_ip ) then

            ! Assembly RHS
            do iv = 1, DVS
               if( lmask_loc(iv) ) then
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elrhs(iv,:), rhsid )
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elfrx(iv,:), frxid_sld)
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elfin(iv,:), finte_sld )
                  call matrix_assemble_element_RHS( &
                     ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                     elfex(iv,:), fexte_sld )
               endif
            enddo
            
            ! Lumped contact force
            if( kfl_conta_sld /= 0_ip )then
               do iv = 1, DVS
                  if( lmask_loc(iv) ) then
                     call matrix_assemble_element_RHS( &
                        ZDOFN, ZDOFN, ZNODE, lnods_loc(iv,:), &
                        elfco(iv,:), fcont_sld )
                  endif
               enddo
            endif

            ! Assemble LUMPED MASS MATRIX
            do iv = 1, DVS
               if( lmask_loc(iv) ) then
                  call matrix_assemble_element_RHS( &
                       1_ip, 1_ip, ZNODE, lnods_loc(iv,:), &
                       elmuu(iv,:), vmass_sld )
               endif
            enddo

         endif not_postprocess

      end subroutine sld_vect_assemble_global_RHS

      pure subroutine sld_vect_apply_essential_BC_at_ELEM(pnode,lmask_loc,lnods_loc,elrhs,elmat)
         use def_solidz, only : kfl_fixno_sld
   
         integer(ip), intent(in)    :: pnode
         logical(lg), intent(in)    :: lmask_loc(DVS)
         integer(ip), intent(in)    :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(inout) :: elrhs(DVS,ZDOFU)
         real(rp),    intent(inout) :: elmat(DVS,ZDOFU,ZDOFU)
         integer(ip)                :: iv, in, jn, id, jd, ii, jj, ipo
         real(rp)                   :: adiag

         do iv = 1, DVS
            if( lmask_loc(iv) ) then
               do in = 1, ZNODE
                  ipo = lnods_loc(iv,in)
                  do ii = 1, ZDIME
                     if( kfl_fixno_sld(ii,ipo) > 0_ip ) then
                        id = ( in - 1_ip) * ZDOFN + ii
                        adiag = elmat(iv,id,id)
                        do jn = 1, ZNODE
                           do jj = 1, ZDIME
                              jd = (jn - 1_ip) * ZDOFN + jj
                              elmat(iv,id,jd) = 0.0_rp
                              elmat(iv,jd,id) = 0.0_rp
                           end do
                        end do
                        elrhs(iv,id)    = 0.0_rp
                        elmat(iv,id,id) = adiag
                     end if
                  end do
               end do
            endif
         enddo

      end subroutine sld_vect_apply_essential_BC_at_ELEM


      pure subroutine sld_vect_get_reaction_force_at_ELEM(pnode,lmask_loc,lnods_loc,elrhs,elfrx)
         use def_solidz, only : kfl_fixno_sld
   
         integer(ip), intent(in)    :: pnode
         logical(lg), intent(in)    :: lmask_loc(DVS)
         integer(ip), intent(in)    :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(inout) :: elrhs(DVS,ZDOFU)
         real(rp),    intent(out)   :: elfrx(DVS,ZDOFU)
         integer(ip)                :: iv, ino, id, ii, ipo

         do iv = 1, DVS
            elfrx(iv,:) = 0.0_rp
           if( lmask_loc(iv) ) then
               do ino = 1, ZNODE
                  ipo = lnods_loc(iv,ino)
                  do ii = 1, ZDIME
                     if( kfl_fixno_sld(ii,ipo) > 0_ip ) then
                        id = ( ino - 1_ip) * ZDOFN + ii
                        elfrx(iv,id) = elrhs(iv,id) 
                     end if
                  end do
               end do
           endif
         enddo

      end subroutine sld_vect_get_reaction_force_at_ELEM


      pure subroutine sld_vect_get_contact_force_at_ELEM(pnode,lmask_loc,lnods_loc,elrhs,elfco)
         use def_solidz, only : kfl_fixno_sld
         use mod_sld_bclocal, only : sld_bclocal_elmat_and_elrhs

         integer(ip), intent(in)    :: pnode
         logical(lg), intent(in)    :: lmask_loc(DVS)
         integer(ip), intent(in)    :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(inout) :: elrhs(DVS,ZDOFU)
         real(rp),    intent(out)   :: elfco(DVS,ZDOFU)
         integer(ip)                :: iv, ino, ipo, ii, jj

         do iv = 1, DVS
            elfco(iv,:) = 0.0_rp
            if( lmask_loc(iv) )then
               do ino = 1, ZNODE
                  ipo = lnods_loc(iv,ino)
                  ii = ( ino - 1_ip ) * ZDOFN + 1_ip
                  jj = ( ino - 1_ip ) * ZDOFN + ZDIME
                  if( kfl_fixno_sld(1,ipo) == 3_ip ) then
                     elfco(iv,ii:jj) = elrhs(iv,ii:jj)
                  end if
               end do
            endif
         enddo

      end subroutine sld_vect_get_contact_force_at_ELEM


      subroutine sld_vect_change_csys_from_global_to_local(pnode,lmask_loc,lnods_loc,elrhs)
         use mod_sld_bclocal, only: sld_bclocal_elmat_and_elrhs

         integer(ip), intent(in)    :: pnode
         logical(lg), intent(in)    :: lmask_loc(DVS)
         integer(ip), intent(in)    :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(inout) :: elrhs(DVS,ZDOFU)
         integer(ip)                :: iv

         do iv = 1, DVS
            if( lmask_loc(iv) )then
               call sld_bclocal_elmat_and_elrhs( &
                  ZDIME,ZNODE,lnods_loc(iv,:),elrhs(iv,:))
            endif
         enddo

      end subroutine sld_vect_change_csys_from_global_to_local


      subroutine sld_vect_change_csys_from_global_to_local_elrhs_elmat(pnode,lmask_loc,lnods_loc,elrhs,elmat)
         use def_domain, only: lpoty
         use def_solidz, only: kfl_fixno_sld, kfl_fixrs_sld, jacrot_du_dq_sld, jacrot_dq_du_sld
         
         external                   :: sld_rotsys
         
         integer(ip), intent(in)    :: pnode
         logical(lg), intent(in)    :: lmask_loc(DVS)
         integer(ip), intent(in)    :: lnods_loc(DVS,ZNODE)
         real(rp),    intent(inout) :: elrhs(DVS,ZDOFU)
         real(rp),    intent(inout) :: elmat(DVS,ZDOFU,ZDOFU)
         integer(ip)                :: iv, ipoin, ibopo, ievat, inode

         do iv = 1, DVS
            if( lmask_loc(iv) )then
               do inode = 1,ZNODE
                  ipoin = lnods_loc(iv,inode)
                  ibopo = lpoty(ipoin)
                  ievat = (inode - 1_ip) * ZDOFN
                  if (ibopo > 0) then
                     if( kfl_fixno_sld(1,ipoin) == 2_ip .or. kfl_fixno_sld(1,ipoin) == 3_ip .or. & 
                         kfl_fixrs_sld(ipoin) /= 0_ip ) then
                         call sld_rotsys(1_ip, inode, ZNODE, ZDOFN, ZDOFU, elmat(iv,:,:), elrhs(iv,:), &
                            jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
                     end if
                  end if
                end do
            endif
         enddo
      end subroutine sld_vect_change_csys_from_global_to_local_elrhs_elmat
      
end module mod_sld_vect_assembly_elfem
!> @}
