!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Mass_matrix_Toolbox
!> Toolbox for mass matrix calculations
!> @{
!> @file    mod_mass_matrix.f90
!> @author  Guillaume Houzeaux
!> @date    15/04/2019
!> @brief   Compute different mass matrices
!>
!>          Lumped mass open rule
!>          Lumped mass close rule
!>          Consistent mass
!>          Weighted Consistent mass
!>          Weighted lumped mass open rule by diagonal scale
!>          Inverse weighted mass aprroximation
!>
!------------------------------------------------------------------------

module mod_mass_matrix

  use def_kintyp,              only : ip,rp,lg
  use def_parame,              only : twopi
  use def_elmtyp,              only : ELEXT
  use def_elmtyp,              only : NOHOL
  use def_elmtyp,              only : PYR05  
  use def_elmtyp,              only : BAR3D  
  use def_master,              only : IMASTER
  use def_master,              only : INOTMASTER
  use def_master,              only : INOTEMPTY
  use def_master,              only : zeror,intost
  use def_master,              only : lninv_loc
  use def_kermod,              only : ndivi
  use def_kermod,              only : kfl_grpro
  use def_kermod,              only : kfl_conma
  use def_kermod,              only : kfl_conma_weighted
  use def_kermod,              only : kfl_approx_inv_mass
  use def_kermod,              only : kfl_element_to_csr
  use def_kermod,              only : kfl_dmass
  use mod_ker_proper,          only : ker_proper
  use def_domain,              only : memor_dom
  use def_domain,              only : vmass,meshe
  use def_domain,              only : vmasc,ndimb
  use def_domain,              only : cmass,cmass_weighted
  use def_domain,              only : nelem,lnods
  use def_domain,              only : ltype,ldime
  use def_domain,              only : lelch,lnoch
  use def_domain,              only : nnode,ngaus
  use def_domain,              only : nhang,npoin
  use def_domain,              only : coord,elmar
  use def_domain,              only : ndime,mnode
  use def_domain,              only : lorde,lexis
  use def_domain,              only : kfl_naxis
  use def_domain,              only : kfl_spher
  use def_domain,              only : kfl_horde
  use def_domain,              only : mgaus
  use def_domain,              only : c_dom,r_dom,lezdo
  use def_domain,              only : nzdom,dmass
  use def_domain,              only : mass_right_fact
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use mod_matrix,              only : matrix_assemble_element_matrix_to_CSR
  use mod_matrix,              only : matrix_add_diagonal_CSR
  use mod_matrix,              only : matrix_scaling_CSR
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_memory,              only : memory_copy
  use mod_messages,            only : messages_live
  use mod_messages,            only : livinf
  use mod_outfor,              only : outfor
  use mod_domain,              only : domain_memory_reallocate
  use mod_domain,              only : domain_memory_deallocate
  implicit none

  integer(ip)             :: ipoin,inode,igaus,ielem,pdime,knode
  integer(ip)             :: pgaus,pnode,pelty,ichkm
  real(rp)                :: gpdet,xjacm(9)
  real(rp)                :: aux,time1,time2,T,d
#if defined(_OPENMP)
  real(rp)                :: xfact
#endif
  character(20)           :: chnod

  private

  public :: mass_matrix_open_lumped
  public :: mass_matrix_close
  public :: mass_matrix_consistent
  public :: mass_matrix_approx_inv_mass
  public :: mass_matrix_diagonal_lumped
  public :: mass_matrix_destructor

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-15
  !> @brief   Destroy
  !> @details Destroy all mass matrices
  !> 
  !-----------------------------------------------------------------------

  subroutine mass_matrix_destructor()

    call domain_memory_deallocate('VMASS')
    call domain_memory_deallocate('VMASC')
    call domain_memory_deallocate('CMASS')
    call domain_memory_deallocate('CMASS_WEIGHTED')
    call domain_memory_deallocate('MASS_RIGHT_FACT')
    call domain_memory_deallocate('DMASS')

  end subroutine mass_matrix_destructor
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-15
  !> @brief   Lumped mass matrix
  !> @details Compute lumped diagonal mass matrix VMASS(NPOIN)
  !> 
  !-----------------------------------------------------------------------

  subroutine mass_matrix_open_lumped()

    real(rp) :: elcod(ndime,mnode)
    real(rp) :: numer(mnode),gpvol

    call cputim(time1)

    call messages_live('COMPUTE LUMPED MASS MATRIX')

    if( INOTEMPTY ) then
       !
       ! Allocate memory
       !
       call domain_memory_reallocate('VMASS')       
       do ipoin = 1,npoin
          vmass(ipoin) = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       if( kfl_naxis == 0 .and. kfl_spher == 0 ) then

          if( kfl_horde == 0 ) then
             !
             ! Only linear elements
             !
             !$OMP  PARALLEL DO SCHEDULE (STATIC)                                     &
             !$OMP  DEFAULT (NONE)                                                    &
             !$OMP  PRIVATE ( aux, elcod, gpdet, gpvol,                               &
             !$OMP            ielem, igaus, inode, ipoin, numer, pelty, pgaus, pnode, &
             !$OMP            xjacm, xfact, pdime )                                   &
             !$OMP  SHARED  ( coord, elmar, lnods, ltype, nelem, ngaus, nnode,        &
#ifndef NDIMEPAR
             !$OMP            ndime,                                                  &
#endif
             !$OMP            vmass, kfl_naxis, kfl_spher , lelch, ldime )
             do ielem = 1,nelem
                !
                ! Element properties and dimensions
                !
                pelty = ltype(ielem)
                if( pelty > 0 ) then

                   pdime = ldime(pelty)
                   pnode = nnode(pelty)
                   pgaus = ngaus(pelty)
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                    !  if(ipoin==0) print*,'a=',pelty,lnods(:,ielem)
                      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                   end do

                   do inode = 1,pnode
                      numer(inode) = 0.0_rp
                   end do

                   if( pelty == BAR3D ) then
                      call mass_matrix_BARD3D(elcod,numer)
                   else
                      do igaus = 1,pgaus
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty) % deriv(1,1,igaus),&
                              xjacm,gpdet)
                         gpvol = elmar(pelty) % weigp(igaus) * gpdet
                         do inode = 1,pnode
                            numer(inode) = numer(inode) &
                                 + gpvol * elmar(pelty) % shape(inode,igaus)
                         end do
                      end do
                   end if

                   if( lelch(ielem) == ELEXT ) then
                      inode        = 1
                      ipoin        = lnods(inode,ielem)
                      aux          = numer(inode)
                      !$OMP         ATOMIC
                      vmass(ipoin) = vmass(ipoin) + aux
                   else
                      do inode = 1,pnode
                         ipoin        = lnods(inode,ielem)
                         aux          = numer(inode)
                         !$OMP         ATOMIC
                         vmass(ipoin) = vmass(ipoin) + aux
                      end do
                   end if

                end if

             end do

          else
             !
             ! There are quadratic elements
             !
             !$OMP  PARALLEL DO SCHEDULE (STATIC)                                     &
             !$OMP  DEFAULT (NONE)                                                    &
             !$OMP  PRIVATE ( aux,   elcod, gpdet, gpvol,                             &
             !$OMP            ielem, igaus, inode, ipoin, numer, pelty, pgaus, pnode, &
             !$OMP            xjacm, xfact, t, D, pdime )                             &
             !$OMP  SHARED  ( coord, elmar, lnods, ltype, nelem, ngaus, nnode,        &
#ifndef NDIMEPAR
             !$OMP            ndime,                                                  &
#endif
             !$OMP            vmass, kfl_naxis, kfl_spher ,lelch, lorde, ldime )
             do ielem = 1,nelem

                pelty = ltype(ielem)

                if( pelty > 0 ) then

                   pdime = ldime(pelty)
                   pnode = nnode(pelty)
                   pgaus = ngaus(pelty)
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                      numer(inode) = 0.0_rp
                   end do

                   if( lorde(pelty) == 1 ) then

                      if( pelty == BAR3D ) then
                         call mass_matrix_BARD3D(elcod,numer)
                      else
                         do igaus = 1,pgaus
                            call jacdet(&
                                 pdime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                                 xjacm,gpdet)
                            gpvol = elmar(pelty)%weigp(igaus) * gpdet
                            do inode = 1,pnode
                               numer(inode) = numer(inode) &
                                    + gpvol * elmar(pelty)%shape(inode,igaus)
                            end do
                         end do
                         if( lelch(ielem) == ELEXT ) then
                            inode        = 1
                            ipoin        = lnods(inode,ielem)
                            aux          = numer(inode)
                            !$OMP         ATOMIC
                            vmass(ipoin) = vmass(ipoin) + aux
                         else
                            do inode = 1,pnode
                               ipoin        = lnods(inode,ielem)
                               aux          = numer(inode)
                               !$OMP         ATOMIC
                               vmass(ipoin) = vmass(ipoin) + aux
                            end do
                         end if
                      end if

                   else
                      !
                      ! This is a more than linear element: Use D-matrix
                      !
                      T = 0.0_rp
                      d = 0.0_rp
                      do igaus = 1,pgaus
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                              xjacm,gpdet)
                         gpvol = elmar(pelty)%weigp(igaus) * gpdet
                         T     = T + gpvol
                         do inode = 1,pnode
                            numer(inode) = numer(inode) &
                                 + gpvol &
                                 * elmar(pelty)%shape(inode,igaus) &
                                 * elmar(pelty)%shape(inode,igaus)
                         end do
                      end do
                      do inode = 1,pnode
                         d = d + numer(inode)
                      end do
                      d = 1.0_rp / d
                      do inode = 1,pnode
                         numer(inode) = numer(inode) * T * d
                      end do
                      if( lelch(ielem) == ELEXT ) then
                         inode        = 1
                         ipoin        = lnods(inode,ielem)
                         aux          = numer(inode)
                         !$OMP         ATOMIC
                         vmass(ipoin) = vmass(ipoin) + aux
                      else
                         do inode = 1,pnode
                            ipoin        = lnods(inode,ielem)
                            aux          = numer(inode)
                            !$OMP         ATOMIC
                            vmass(ipoin) = vmass(ipoin) + aux
                         end do
                      end if

                   end if

                end if

             end do

          end if

       else

          write(6,*) 'SPHERICAL COORDINATES'
          call runend('LUMPED MASS MATRIX NOT CODED IN NON-CARTESIAN COORDINATE SYSTEM')

       end if
       !
       ! Modify RHS due Parall service
       !
       call rhsmod(1_ip,vmass)
       !
       ! Loop over nodes to control zero-volume points
       !
       ichkm = 0
       nodes: do ipoin = 1,npoin
          if( lnoch(ipoin) == NOHOL ) then
             vmass(ipoin) = 1.0_rp
          else if( vmass(ipoin) < zeror ) then
             ichkm = ichkm + 1
             chnod = intost(ipoin)
             !write(6,*) 'pipi',vmass(ipoin),zeror
             write(6,'(a,i10,1(1x,e12.6))') 'MASSMA: lninv_loc(ipoin),vmass(ipoin)= ',lninv_loc(ipoin),vmass(ipoin)
             call livinf(10000_ip, &
                  'NODE NUMBER ( '//adjustl(trim(chnod)) &
                  // ' HAS ZERO MASS)',0_ip)
             vmass(ipoin) = 1.0_rp
          end if
       end do nodes

       if( ichkm > 0 ) then
          !   call runend('MASSMA: ZERO MASS NODES FOUND')
       end if

    end if

    !-------------------------------------------------------------------
    !
    ! Fill in mesh type
    !
    !-------------------------------------------------------------------

    meshe(ndivi) % vmass => vmass

    call cputim(time2)

  end subroutine mass_matrix_open_lumped

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-15
  !> @brief   Mass matrix with close rule
  !> @details Compute diagonal mass matrix with close rule VMASC(NPOIN)
  !> 
  !-----------------------------------------------------------------------

  subroutine mass_matrix_close()

    real(rp) :: elcod(ndime,mnode)
    real(rp) :: numer(mnode),gpvol
    real(rp) :: elmas(mnode)

    call cputim(time1)

    call messages_live('COMPUTE CLOSED MASS MATRIX')

    if( INOTEMPTY ) then
       !
       ! Allocate memory
       !
       call domain_memory_reallocate('VMASC')              
       do ipoin = 1,npoin
          vmasc(ipoin) = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       if( ( ( ndime == 3 .and. lexis(PYR05) == 0 ) .or. ndime == 2 ) .and. kfl_naxis == 0 .and. kfl_spher == 0 ) then
          !
          ! 3D: Avoid if and go fast
          !
          !$OMP  PARALLEL DO SCHEDULE (STATIC)                                     &
          !$OMP  DEFAULT (NONE)                                                    &
          !$OMP  PRIVATE ( aux, elcod, elmas, gpdet, gpvol, xjacm, ielem, igaus,   &
          !$OMP            inode, ipoin, numer, pelty, pgaus, pnode, pdime )       &
          !$OMP  SHARED  ( coord, elmar, lnods, ltype, nelem, ngaus, nnode,        &
#ifndef NDIMEPAR
          !$OMP            ndime,                                                  &
#endif
          !$OMP            vmasc, kfl_naxis, kfl_spher ,lelch, ldime )
          !
          do ielem = 1,nelem
             pelty = ltype(ielem)

             if( pelty > 0 ) then

                pdime = ldime(pelty)
                pnode = nnode(pelty)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do
                if( lelch(ielem) == ELEXT ) then
                   inode = 1
                   call jacdet(&
                        pdime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                        xjacm,gpdet)
                   elmas(inode) = elmar(pelty)%weigc(inode)*gpdet
                   ipoin = lnods(inode,ielem)
                   !$OMP ATOMIC
                   vmasc(ipoin) = vmasc(ipoin) + elmas(inode)
                else
                   if( pelty == BAR3D ) then
                      call mass_matrix_BARD3D(elcod,elmas)
                   else
                      do inode = 1,pnode
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                              xjacm,gpdet)
                         elmas(inode) = elmar(pelty)%weigc(inode)*gpdet
                      end do
                   end if
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      !$OMP ATOMIC
                      vmasc(ipoin) = vmasc(ipoin) + elmas(inode)
                   end do
                end if
             end if

          end do

       else
          !
          !$OMP  PARALLEL DO SCHEDULE (STATIC)                                     &
          !$OMP  DEFAULT (NONE)                                                    &
          !$OMP  PRIVATE ( aux,   elcod, elmas, gpdet, gpvol,                      &
          !$OMP            ielem, igaus, inode, ipoin, numer, pelty, pgaus, pnode, &
          !$OMP            xjacm, knode, pdime )                                   &
          !$OMP  SHARED  ( coord, elmar, lnods, ltype, nelem, ngaus, nnode,        &
#ifndef NDIMEPAR
          !$OMP            ndime,                                                  &
#endif
          !$OMP            vmasc, kfl_naxis, kfl_spher, lelch, ldime )
          !
          do ielem = 1,nelem
             pelty = ltype(ielem)
             if( pelty > 0 ) then
                pdime = ldime(pelty)
                pnode = nnode(pelty)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do
                !
                ! Treat extension elements
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if

                if( pelty == PYR05 ) then
                   !
                   ! PYR05: Pyramid element
                   !
                   pgaus = ngaus(pelty)
                   numer(1:pnode) = 0.0_rp

                   if( pelty == BAR3D ) then
                      call mass_matrix_BARD3D(elcod,elmas)
                   else
                      do igaus = 1,pgaus
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                              xjacm,gpdet)
                         gpvol = elmar(pelty) % weigp(igaus) * gpdet
                         do inode = 1,knode
                            numer(inode) = numer(inode) + gpvol * elmar(pelty) % shape(inode,igaus)
                         end do
                      end do
                   end if
                   
                   do inode = 1,knode
                      ipoin        = lnods(inode,ielem)
                      aux          = numer(inode)
                      !$OMP         ATOMIC
                      vmasc(ipoin) = vmasc(ipoin) + aux
                   end do

                else
                   !
                   ! Other elements: loop over Gauss points (which are nodes)
                   !
                   if( kfl_naxis == 1 ) then

                      do inode=1,pnode
                         ipoin=lnods(inode,ielem)
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                              xjacm,gpdet)
                         gpvol=elmar(pelty)%weigc(inode)*gpdet
                         if( elcod(1,inode) == 0.0_rp ) then
                            gpvol = gpvol * twopi * 1.0e-12_rp
                         else
                            gpvol = gpvol * twopi * elcod(1,inode)
                         end if
                         !$OMP ATOMIC
                         vmasc(ipoin) = vmasc(ipoin) + gpvol
                      end do

                   else if( kfl_spher == 1 ) then

                      do inode=1,pnode
                         ipoin=lnods(inode,ielem)
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                              xjacm,gpdet)
                         gpvol=elmar(pelty)%weigc(inode)*gpdet
                         if( elcod(1,inode) == 0.0_rp ) then
                            gpvol = gpvol * twopi * 1.0e-12_rp
                         else
                            gpvol = gpvol * 2.0_rp * twopi * elcod(1,inode) * elcod(1,inode)
                         end if
                         !$OMP ATOMIC
                         vmasc(ipoin) = vmasc(ipoin) + gpvol
                      end do

                   else

                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         call jacdet(&
                              pdime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                              xjacm,gpdet)
                         elmas(inode) = elmar(pelty)%weigc(inode)*gpdet
                      end do
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         !$OMP ATOMIC
                         vmasc(ipoin) = vmasc(ipoin) + elmas(inode)
                      end do
                   end if

                end if
             end if
          end do
       end if
       !
       ! Modify RHS due to periodicity and Parall service
       !
       call rhsmod(1_ip,vmasc)
       !
       ! Loop over nodes to control zero-volume points
       !
       if( kfl_horde == 0 ) then
          ichkm = 0
          aux   = epsilon(1.0_rp)
          nodes: do ipoin = 1,npoin
             if( lnoch(ipoin) == NOHOL ) then
                vmasc(ipoin) = 1.0_rp
             else if( vmasc(ipoin) < aux ) then
                ichkm = 1
                chnod = intost(ipoin)
                write(*,*)'massmc:lninv_loc(ipoin),vmasc(ipoin),aux',lninv_loc(ipoin),vmasc(ipoin),aux
                call livinf(10000_ip, &
                     'NODE NUMBER ( '//adjustl(trim(chnod)) &
                     // ' HAS ZERO MASS)',0_ip)
                vmasc(ipoin) = 1.0_rp
             end if
          end do nodes

          !        if( ichkm == 1 ) call runend('MASSMC: ZERO MASS NODES FOUND')
       else

       end if

    end if
    !
    ! If there are quadratic elements, use always lumped matrix
    !
    if( kfl_horde /= 0 ) then
       if( kfl_grpro == 1 ) then
          call outfor(2_ip,0_ip,'THERE ARE QUADRATIC ELEMENT: USE ALWAYS LUMPED MATRIX')
          kfl_grpro = 0
       end if
    end if

    !-------------------------------------------------------------------
    !
    ! Fill in mesh type
    !
    !-------------------------------------------------------------------

    meshe(ndivi) % vmasc => vmasc

    call cputim(time2)

  end subroutine mass_matrix_close

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-15
  !> @brief   Consistent mass matrix 
  !> @details Compute consistent mass matrix in CSR format
  !> 
  !-----------------------------------------------------------------------

  subroutine mass_matrix_consistent(CONSISTENT_MASS,CONSISTENT_WEIGHTED_MASS)

    logical(lg), optional, intent(in) :: CONSISTENT_MASS
    logical(lg), optional, intent(in) :: CONSISTENT_WEIGHTED_MASS
    integer(ip)                       :: ipoin,idime,inode,ielem
    integer(ip)                       :: pgaus,pnode,pelty,dummi
    real(rp)                          :: gpsha(mnode,mgaus)
    real(rp)                          :: gpder(ndime,mnode,mgaus)
    real(rp)                          :: gpcar(ndime,mnode,mgaus)
    real(rp)                          :: gpvol(mgaus)
    real(rp)                          :: gpden(mgaus)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: elmat(mnode,mnode)
    real(rp)                          :: elmat_weighted(mnode,mnode)
    logical(lg)                       :: if_consistent_mass
    logical(lg)                       :: if_consistent_weighted_mass
    
    if_consistent_mass          = .false.
    if_consistent_weighted_mass = .false.
    if( present(CONSISTENT_MASS) )          if_consistent_mass          = CONSISTENT_MASS 
    if( present(CONSISTENT_WEIGHTED_MASS) ) if_consistent_weighted_mass = CONSISTENT_WEIGHTED_MASS
    if_consistent_mass          = if_consistent_mass          .and. ( kfl_conma /= 0 )
    if_consistent_weighted_mass = if_consistent_weighted_mass .and. ( kfl_conma_weighted /= 0 )    

    if( if_consistent_mass ) then
       call domain_memory_reallocate('CMASS')       
       call messages_live('COMPUTE CONSISTENT MASS MATRIX')       
    end if
    if( if_consistent_weighted_mass ) then
       call domain_memory_reallocate('CMASS_WEIGHTED')  
       call messages_live('COMPUTE CONSISTENT WEIGHTED MASS MATRIX')       
    end if    
 
    gpden = 1.0_rp
    
    if( ( if_consistent_mass .or. if_consistent_weighted_mass ) .and. INOTEMPTY ) then
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem

          pelty = ltype(ielem) 
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)

          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                elcod(idime,inode) = coord(idime,ipoin)
             end do
          end do

          call element_shape_function_derivatives_jacobian(&
               pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
               elmar(pelty) % deriv,elmar(pelty) % heslo,&
               elcod,gpvol,gpsha,gpder,gpcar,IELEM=ielem)

          if( if_consistent_weighted_mass ) then
             call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)  
             call mass_matrix_elemental(pnode,pgaus,gpvol,gpsha,gpden,elmat,elmat_weighted)
          else if( if_consistent_mass ) then
             call mass_matrix_elemental(pnode,pgaus,gpvol,gpsha,gpden,elmat)
          end if

          if( if_consistent_mass ) then
             call matrix_assemble_element_matrix_to_CSR(&
                  kfl_element_to_csr,1_ip,pnode,pnode,&
                  ielem,lnods(:,ielem),elmat,r_dom,c_dom,cmass,lezdo)
          end if
          if( if_consistent_weighted_mass ) then
             call matrix_assemble_element_matrix_to_CSR(&
                  kfl_element_to_csr,1_ip,pnode,pnode,&
                  ielem,lnods(:,ielem),elmat,r_dom,c_dom,cmass_weighted,lezdo)             
          end if

       end do elements

    end if

  end subroutine mass_matrix_consistent

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-16
  !> @brief   Element mass matrix
  !> @details Compute the element mass matrix
  !> 
  !-----------------------------------------------------------------------

  subroutine mass_matrix_elemental(pnode,pgaus,gpvol,gpsha,gpden,elmat,elmat_weighted)

    integer(ip),           intent(in)  :: pnode,pgaus
    real(rp),              intent(in)  :: gpvol(pgaus)
    real(rp),              intent(in)  :: gpsha(pnode,pgaus)
    real(rp),              intent(in)  :: gpden(pgaus)
    real(rp),              intent(out) :: elmat(pnode,pnode)
    real(rp),    optional, intent(out) :: elmat_weighted(pnode,pnode)
    integer(ip)                        :: igaus,inode,jnode
    real(rp)                           :: fact1,fact2

    elmat = 0.0_rp
    if( present(elmat_weighted) ) elmat_weighted = 0.0_rp

    do igaus = 1,pgaus
       do inode = 1,pnode
          fact1 = gpvol(igaus) * gpsha(inode,igaus)
          if( present(elmat_weighted) ) fact2 = fact1 * gpden(igaus)
          do jnode = 1,pnode
             elmat(inode,jnode) = elmat(inode,jnode) + fact1 * gpsha(jnode,igaus)
             if( present(elmat_weighted) ) elmat_weighted(inode,jnode) = elmat_weighted(inode,jnode) &
                  + fact2 * gpsha(jnode,igaus)
          end do
       end do
    end do

  end subroutine mass_matrix_elemental

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-14
  !> @brief   Close mass matrix of BAR3D element
  !> @details Lumped mass matrix of BAR3D element of length l
  !>
  !>          1                  2
  !>          o------------------o
  !>                  +-   -+
  !>                  | 1 0 |
  !>          M = 1/2 | 0 1 |
  !>                  +-   -+
  !> 
  !-----------------------------------------------------------------------

  subroutine mass_matrix_BARD3D(bocod,numer)

    real(rp),    intent(in)  :: bocod(ndime,2)
    real(rp),    intent(out) :: numer(2)
    real(rp)                 :: ll

    ll       = sqrt(dot_product(bocod(:,2)-bocod(:,1),bocod(:,2)-bocod(:,1)))
    numer(1) = 0.5_rp * ll
    numer(2) = 0.5_rp * ll
    
  end subroutine mass_matrix_BARD3D

  !-----------------------------------------------------------------------
  !> 
  !> @author  Lucas Gasparino
  !> @date    2019-12-02
  !> @brief   Approximate inverse mass matrix
  !> @details Computes an approximation to inv(Mc) according to Guermond and Pasquetti reference
  !>
  !>                         inv(Mc) = (I+A)*inv(Ml)
  !>                         A = inv(Ml)*(Ml+Mc)
  !> 
  !>          A is stored in case A**2 is to be added to the approximation
  !>          Uses weighted matrices (Ml_weigh and Mc_weigh)
  !>
  !-----------------------------------------------------------------------

  subroutine mass_matrix_approx_inv_mass()

    ! Have to change vmass and cmass to corresp. weighted forms

    implicit none

    integer(ip)             :: ipoin,iz
    real(rp),    pointer    :: diag_mass(:)

    nullify(diag_mass)

    if( kfl_approx_inv_mass /= 0 ) then

       call messages_live('COMPUTE APPROXIMATE INVERSE MASS MATRIX')
       !
       ! Compute weighted lumped mass
       !
       call mass_matrix_diagonal_lumped(EXCHANGE=.false.,DIAG_MASS=diag_mass)
       !
       ! Compute A in CSR format
       !
       call domain_memory_reallocate('MASS_RIGHT_FACT')       
       do iz = 1,nzdom
          mass_right_fact(iz) = -1.0_rp*cmass(iz)
       end do
       
       if( INOTEMPTY ) then
          !
          ! Ml* - Mc* 
          !
          call matrix_add_diagonal_CSR(npoin,1_ip,r_dom,c_dom,diag_mass,mass_right_fact,1.0_rp)
          !
          ! Exhange lumped mass Ml*
          !
          call rhsmod(1_ip,diag_mass)
          !
          ! Scaling: A = Ml^-1 ( Ml* - Mc* )
          !
          do ipoin = 1,npoin
             if( diag_mass(ipoin) /= 0.0_rp ) diag_mass(ipoin) = 1.0_rp / diag_mass(ipoin) 
          end do
          call matrix_scaling_CSR(1_ip,npoin,1_ip,r_dom,c_dom,mass_right_fact,diag_mass,LEFT_SCALING=.true.)

       end if
       
       call memory_deallo(memor_dom,'DIAG_MASS','mod_mass_matrix',diag_mass)                 

    end if

  end subroutine mass_matrix_approx_inv_mass

  !-----------------------------------------------------------------------
  !> 
  !> @author  Lucas Gasparino
  !> @date    2019-12-02
  !> @brief   Weighted lumped mass matrix
  !> @details Computes  the weighted lumped mass by diagonal scaling
  !>
  !>                         alpha = 1/((gpsha(inode,igaus)**2)
  !>                         Ma)_l = w_k*gpden_k*|Je_k|*((Na_k)**2)
  !> 
  !>          See T. J. R. Hughes for a better description.
  !>          This one is useful for any element type/order.
  !-----------------------------------------------------------------------

  subroutine mass_matrix_diagonal_lumped(EXCHANGE,diag_mass)

    logical(lg),          optional, intent(in)    :: EXCHANGE
    real(rp),    pointer, optional, intent(inout) :: diag_mass(:)
    integer(ip)                                   :: ipoin,inode,ielem
    integer(ip)                                   :: pgaus,pnode,pelty
    real(rp)                                      :: gpsha(mnode,mgaus)
    real(rp)                                      :: gpder(ndime,mnode,mgaus)
    real(rp)                                      :: gpcar(ndime,mnode,mgaus)
    real(rp)                                      :: gpvol(mgaus), gpvol_weighted(mgaus)
    real(rp)                                      :: elcod(ndime,mnode)
    real(rp)                                      :: alpha, aux
    logical(lg)                                   :: if_exchange
    real(rp),   pointer                           :: dmass_loc(:)
    
    if( present(EXCHANGE) ) then
       if_exchange = EXCHANGE
    else       
       if_exchange = .true.
    end if
    
    if( present(diag_mass) .or. kfl_dmass /= 0 ) then
       
       if( present(diag_mass) ) then
          if( .not. associated(diag_mass) ) then
             call memory_alloca(memor_dom,'DIAG_MASS','mod_mass_matrix',diag_mass,npoin)                 
          end if
          dmass_loc => diag_mass
       else
          call domain_memory_reallocate('DMASS')       
          dmass_loc => dmass
       end if

       if( .not. present(diag_mass) ) call messages_live('COMPUTE DIAGONALLY SCALED LUMPED MASS MATRIX')
       !
       ! Loop over elements
       !
       do ielem = 1,nelem

          !
          ! Element data
          !
          pelty = ltype(ielem)
          pdime = ldime(pelty)
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(1:ndime,inode) = coord(1:ndime,ipoin)
          end do

          !
          ! Compute gpvol_weighted, gpvol, gpdet and gpden at ALL GauB points
          !
          call element_shape_function_derivatives_jacobian(&
               pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
               elmar(pelty) % deriv,elmar(pelty) % heslo,&
               elcod,gpvol,gpsha,gpder,gpcar,IELEM=ielem)
          do igaus = 1,pgaus
             gpvol_weighted(igaus) = gpvol(igaus)
          end do
          !
          ! Compute alpha = 1/(Na_k*Na_k)
          !
          alpha = 0.0_rp
          do inode = 1,pnode
             do igaus = 1,pgaus
                alpha = alpha+(gpsha(inode,igaus)**2)*gpvol_weighted(igaus)
             end do
          end do
          alpha = sum(gpvol_weighted)/alpha

          !
          ! Compute global entries of dmass
          !
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             aux = 0.0_rp
             do igaus = 1,pgaus
                aux = aux+gpvol_weighted(igaus)*(gpsha(inode,igaus)**2) ! W_k*gpden_k*|Je_k|*((Na_k)**2)
             end do
             dmass_loc(ipoin) = dmass_loc(ipoin)+alpha*aux
          end do
       end do
       !
       ! Modify RHS due Parall service
       !
       if( if_exchange .and. INOTEMPTY ) call rhsmod(1_ip,dmass_loc)

    end if

  end subroutine mass_matrix_diagonal_lumped
                                      
end module mod_mass_matrix
!> @}


