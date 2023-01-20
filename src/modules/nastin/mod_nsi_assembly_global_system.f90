!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_nsi_assembly_global_system

  !------------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    nsi_elmmat.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Navier-Stokes system element assembly and other element
  !>          calculations
  !> @details Compute element matrices
  !>
  !>          \verbatim
  !>
  !>          Without enrichement
  !>          -------------------
  !>            +-        +  +- -+     +-  -+
  !>            | Auu Aup |  | u |     | bu |
  !>            |         |  |   |  =  |    |
  !>            | Apu App |  | p |     | bp |
  !>            +-       -+  +- -+     +-  -+
  !>
  !>          With enrichement
  !>          ----------------
  !>            +-            +  +- -+     +-  -+
  !>            | Auu Aup Auq |  | u |     | bu |
  !>            |             |  |   |     |    |
  !>            | Apu App Apq |  | p |  =  | bp |
  !>            |             |  |   |     |    |
  !>            | Aqu Aqp Aqq |  | q |     | bq |
  !>            +-           -+  +- -+     +-  -+
  !>
  !>            q = Aqq^-1 ( bq - Aqu u - Aqp p )
  !>
  !>            Auu <= Auu - Auq Aqq^-1 Aqu
  !>            Aup <= Aup - Auq Aqq^-1 Aqp
  !>            bu  <= bu  - Auq Aqq^-1 bq
  !>
  !>            Apu <= Apu - Apq Aqq^-1 Aqu
  !>            App <= App - Apq Aqq^-1 Aqp
  !>            bp  <= bp  - Apq Aqq^-1 bq
  !>
  !>          \endverbatim
  !>
  !------------------------------------------------------------------------

#include "def_vector_size.inc"
  use def_kintyp, only :  ip,rp,lg
  use def_master, only : solve
  use def_domain, only : ndime
  use def_domain, only : nzdom,r_dom,c_dom
  use def_domain, only : r_sol,c_sol,lezdo
  use def_domain, only : r_sym,c_sym
  use def_elmtyp, only : ELFEM
  use def_elmtyp, only : ELEXT
  use def_solver, only : SOL_MATRIX_HAS_CHANGED
  use def_solver, only : SOL_YES
  use def_kermod, only : kfl_element_to_csr
  use def_nastin, only : pabdf_nsi
  use def_nastin, only : kfl_predi_nsi
  use def_nastin, only : NSI_FRACTIONAL_STEP
  use def_nastin, only : NSI_DIRICHLET_ELEMENT
  use def_nastin, only : NSI_SOLVER_VISCOUS_TERM

  implicit none
  private

  interface nsi_assembly_schur_method
     module procedure nsi_assembly_schur_method_scalar,&
          &           nsi_assembly_schur_method_vector
  end interface nsi_assembly_schur_method

  public :: nsi_assembly_monolithic
  public :: nsi_assembly_schur_method
  public :: nsi_assembly_fractional_step
  public :: nsi_assembly_fractional_step_rhs
  public :: nsi_assembly_fractional_step_scalar
  public :: nsi_assembly_fractional_step_boundary_scalar
  public :: nsi_assembly_algebraic_split_oss
  public :: nsi_assembly_dt_rho_tau_nu
  public :: nsi_assembly_scaling_pressure_equation
  public :: nsi_assembly_semi_implicit_method
  
contains

  !-----------------------------------------------------------------------
  !
  !> @date    21/02/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Stabilize pressure equation
  !> @details Algebraic split method. Pressure is stabilized so that it
  !>          gives the same pressure stabilization if the stabilization
  !>          parameter is the same as the one of the FS method.
  !>
  !>          Substitute the pressure equation
  !>          Apu u = bp
  !>
  !>          by
  !>
  !>          Apu u + App p = bp + S p^{i-1}, where
  !>
  !>          App  = Stab L
  !>          S    = - ( Apu Stab M^{-1} Aup ) p^{i-1}
  !>          Stab = Tim or Tau
  !>          L    =        \int_V ( grad(Ni).grad(Nj) ) dV
  !>          M    =        \int_V ( Ni.Nj ) dV
  !>          Tim  = M^{-1} \int_V ( dt/rho * Ni.Nj ) dV
  !>          Tau  = M^{-1} \int_V ( tau * Ni.Nj ) dV
  !>
  !>          The preconditioner will be computed further on
  !>          and should approximate
  !>
  !>          Q    = App + Uzawa
  !>               = Stab L + ( 1 / Tim + 1 / Tau )^{-1} L
  !>               = Scal L
  !>          Scal = Stab + ( 1 / Tim + 1 / Tau )^{-1}
  !>
  !>          Thus, we eventually solve:
  !>
  !>          Scal^{-1} ( Aup u + App p ) = Scal^{-1} ( bp + S p )
  !>          with Q = L
  !>
  !>          Note: Laplacian scaling must be cast to the RHS so
  !>          that matrix Q is symmetric. In fact: Scal L is not
  !>          symmetric in general
  !
  !-----------------------------------------------------------------------

  subroutine nsi_assembly_algebraic_split_oss(Aup,Apu,App,L,bp,pp)

    use def_kintyp,               only : ip,rp
    use def_master,               only : INOTMASTER
    use def_master,               only : modul,mem_modul
    use def_domain,               only : ndime,npoin
    use def_domain,               only : r_dom
    use def_domain,               only : vmass,nzdom
    use mod_memory,               only : memory_alloca
    use mod_memory,               only : memory_deallo
    use def_nastin,               only : kfl_stabi_nsi
    use def_nastin,               only : NSI_ALGEBRAIC_SPLIT_OSS
    use def_nastin,               only : dt_rho_nsi
    use def_nastin,               only : tau_nsi
    use mod_nsi_schur_operations, only : nsi_apuvec
    use mod_nsi_schur_operations, only : nsi_aupvec

    use mod_solver

    implicit none

    real(rp),    intent(inout)         :: Aup(ndime,nzdom)
    real(rp),    intent(inout)         :: Apu(ndime,nzdom)
    real(rp),    intent(inout)         :: App(nzdom)
    real(rp),    intent(inout)         :: L(nzdom)
    real(rp),    intent(inout)         :: bp(npoin)
    real(rp),    intent(out)           :: pp(npoin)

    integer(ip)                        :: ipoin,izdom
    real(rp),    pointer               :: Aupp(:,:)
    real(rp),    pointer               :: ww(:)
    real(rp),    pointer               :: Mass(:)
    real(rp),    pointer               :: Stab(:)
    real(rp),    pointer               :: Tim(:)
    real(rp),    pointer               :: Tau(:)

    if( kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS .and. INOTMASTER ) then

       nullify(Aupp)
       nullify(ww)
       nullify(Mass)
       nullify(Stab)
       call memory_alloca(mem_modul(1:2,modul),'AUPP','nsi_assembly_algebraic_split_oss',Aupp,ndime,npoin)
       call memory_alloca(mem_modul(1:2,modul),'WW'  ,'nsi_assembly_algebraic_split_oss',ww,  npoin)

       Mass => vmass
       Tim  => dt_rho_nsi
       Tau  => tau_nsi

       if( 1 == 1 ) then
          !
          ! Pressure is stabilized with Tau L + Apu Tau M^{-1} Aup
          !
          Stab => Tau

       else
          !
          ! Pressure is stabilized with Tim L + Apu Tim M^{-1} Aup
          !
          ! This stabilization is useful to check the result
          ! is the same as the FS result using a fixed time step
          !
          Stab => Tim

       end if
       !
       ! bp = bp - ( Apu Stab M^{-1} Aup ) p
       !
       call nsi_aupvec(1_ip,Aup,pp,Aupp)
       call rhsmod(ndime,Aupp)
       do ipoin = 1,npoin
          Aupp(:,ipoin) = Stab(ipoin) / Mass(ipoin) * Aupp(:,ipoin)
       end do
       call nsi_apuvec(1_ip,Apu,Aupp,ww)
       !
       ! bp  <= ( bp - ( Apu Stab M^{-1} Aup ) p )
       ! App <= Stab L
       ! Apu <= Apu
       !
       do ipoin = 1,npoin
          bp(ipoin) = ( bp(ipoin) - ww(ipoin) )
          do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
             App(izdom)   = Stab(ipoin) * L(izdom)
             Apu(:,izdom) = Apu(:,izdom)
          end do
       end do
       !
       ! Deallocate
       !
       call memory_deallo(mem_modul(1:2,modul),'AUPP','nsi_assembly_algebraic_split_oss',Aupp)
       call memory_deallo(mem_modul(1:2,modul),'WW'  ,'nsi_assembly_algebraic_split_oss',ww)

    end if

  end subroutine nsi_assembly_algebraic_split_oss

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    6/10/2016
  !> @brief   Fractional step assembly
  !> @details Send some terms to the right and assemble matrices
  !>          for the fractional step method
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_assembly_fractional_step(&
       pnode,pevat,list_elements,lnods,elvel,elpre,&
       elauu,elaup,elapp,elapu,elmap,elrbu,elrbp,Aup,Apu,&
       bu,bp,Q,Auu,App)

    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  nzdom,ndime,r_sol
    use def_domain, only       :  c_sol,lezdo
    use def_domain, only       :  lpoty,exnor,skcos
    use def_master, only       :  solve
    use def_kermod, only       :  kfl_element_to_csr
    use def_nastin, only       :  kfl_fixrs_nsi
    use def_nastin, only       :  kfl_matdi_nsi
    use def_nastin, only       :  kfl_local_nsi
    use def_nastin, only       :  skcos_nsi


    integer(ip), intent(in)              :: pnode
    integer(ip), intent(in)              :: pevat
    integer(ip), intent(in)              :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)              :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)              :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)              :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)              :: elauu(VECTOR_SIZE,pevat,pevat)
    real(rp),    intent(in)              :: elaup(VECTOR_SIZE,pevat,pnode)
    real(rp),    intent(in)              :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)              :: elapu(VECTOR_SIZE,pnode,pevat)
    real(rp),    intent(in)              :: elmap(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(inout)           :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(inout)           :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout)           :: Aup(ndime,nzdom)
    real(rp),    intent(inout)           :: Apu(ndime,nzdom)
    real(rp),    intent(inout)           :: bu(ndime,*)
    real(rp),    intent(inout)           :: bp(*)
    real(rp),    intent(inout), optional :: Q(solve(2)%nzmat)
    real(rp),    intent(inout), optional :: Auu(ndime,ndime,*)
    real(rp),    intent(inout), optional :: App(*)
    integer(ip)                          :: ndofn,inode,jnode,iposi,jposi,ielem
    integer(ip)                          :: idime,jdime,izsol,jpoin,ipoin,jcolu,ivect
    integer(ip)                          :: iposi0,jposi0,ievat,jevat,iroty,ibopo
    real(rp)                             :: elvel_tmp(1:VECTOR_SIZE,ndime,pnode)
    real(rp)                             :: worve(3),worma(3)

    !----------------------------------------------------------------
    !
    ! Compute bu^(e) <= bu^(e) - Auu^(e) u^(e)
    ! Compute bp^(e) <= bp^(e) - App^(e) p^{e}
    !
    !----------------------------------------------------------------

    if( .not. present(Auu) ) then

       elvel_tmp(1:VECTOR_SIZE,1:ndime,1:pnode) = elvel(1:VECTOR_SIZE,1:ndime,1:pnode,1)

       if( kfl_local_nsi /= 0 .and. kfl_matdi_nsi == NSI_DIRICHLET_ELEMENT ) then !.and. ( .not. NSI_FRACTIONAL_STEP ) ) then
          !
          ! Rotate the velocity to local system, as Auu has been rotated previously
          !
          do ivect = 1,VECTOR_SIZE
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                if( list_elements(ivect) /= 0 ) then
                   ibopo = lpoty(ipoin)
                else
                   ibopo = 0
                end if
                if( ibopo > 0 ) then
                   iroty = kfl_fixrs_nsi(ipoin)

                   if( iroty == -1 ) then
                      !
                      ! Boundary conditions in the tangent skew system
                      !
                      worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                      call mbvatb(worma,exnor(1,1,ibopo),worve,ndime,ndime)
                      elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                   else if( iroty == -2 ) then
                      !
                      ! Boundary conditions in the NSI tangent skew system
                      !
                      worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                      call mbvatb(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
                      elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                   else if( iroty == -3 ) then
                      !
                      ! Boundary conditions in geometrical system
                      !
                      worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                      call mbvatb(worma,skcos(1,1,ibopo),worve,ndime,ndime)
                      elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                   else if( iroty >= 1 ) then
                      !
                      ! Boundary conditions in a given skew system
                      !
                      worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                      call mbvatb(worma,skcos(1,1,iroty),worve,ndime,ndime)
                      elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                   end if
                end if
             end do
          end do
       end if

       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndime+idime
                   elrbu(1:VECTOR_SIZE,idime,inode) = elrbu(1:VECTOR_SIZE,idime,inode) &
                        - elauu(1:VECTOR_SIZE,ievat,jevat) * elvel_tmp(1:VECTOR_SIZE,jdime,jnode)
                end do
             end do
          end do
       end do
       do jnode = 1,pnode
          do inode = 1,pnode
             elrbp(1:VECTOR_SIZE,inode) = elrbp(1:VECTOR_SIZE,inode) &
                  - elapp(1:VECTOR_SIZE,inode,jnode) * elpre(1:VECTOR_SIZE,jnode,1)
          end do
       end do
    end if

    !if( present(Q) ) solve(2) % kfl_assem = SOL_MATRIX_HAS_CHANGED

    ndofn = ndime + 1

    !----------------------------------------------------------------
    !
    ! Assemble the following matrices:
    !
    ! Aup <= Aup^(e)
    ! Apu <= Apu^(e)
    ! bu  <= bu^(e)
    ! bp  <= bp^(e)
    !
    !----------------------------------------------------------------

    if( solve(2) % kfl_symme == 1 ) then
       !
       ! App requires symmetric assembly
       !
       call runend('SYMMETRIC ASSEMBLY NOT CODED')
    else
       !
       ! App requires asymmetric assembly
       !
       if( kfl_element_to_csr == 1 ) then

          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then

                do inode = 1,pnode
                   iposi0 = (inode-1) * ndime
                   do jnode = 1,pnode
                      jposi0 = (jnode-1) * ndime
                      izsol = lezdo(inode,jnode,ielem)
                      do idime = 1,ndime
                         iposi = iposi0+idime
                         if( present(Auu) ) then
                            do jdime = 1,ndime
                               jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                               !$OMP ATOMIC
#endif
                               Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                                    &                   + elauu(ivect,iposi,jposi)          ! Auu
                            end do
                         end if
                      end do
                      do idime = 1,ndime
                         iposi = iposi0 + idime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Aup(idime,izsol) = Aup(idime,izsol) + elaup(ivect,iposi,jnode) ! Aup
                      end do
                      do jdime = 1,ndime
                         jposi = jposi0 + jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(ivect,inode,jposi)  ! Apu
                      end do
                      if( present(Q) ) then
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Q(izsol) = Q(izsol) + elmap(ivect,inode,jnode)                     ! Q
                      end if
                      if( present(App) ) then
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         App(izsol) = App(izsol) + elapp(ivect,inode,jnode)              ! App
                      end if
                   end do
                end do

             end if
          end do

       else

          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then

                do inode = 1,pnode
                   ipoin = lnods(ivect,inode)
                   do jnode = 1,pnode
                      jpoin = lnods(ivect,jnode)
                      izsol = r_sol(ipoin)
                      jcolu = c_sol(izsol)
                      do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                         izsol = izsol + 1
                         jcolu = c_sol(izsol)
                      end do
                      if( jcolu == jpoin ) then
                         iposi0 = (inode-1) * ndime
                         jposi0 = (jnode-1) * ndime
                         do idime = 1,ndime
                            iposi = iposi0 + idime
                            if( present(Auu) ) then
                               do jdime=1,ndime
                                  jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                                  !$OMP ATOMIC
#endif
                                  Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                                       &                   + elauu(ivect,iposi,jposi)          ! Auu
                               end do
                            end if
                         end do
                         do idime = 1,ndime
                            iposi = iposi0 + idime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Aup(idime,izsol) = Aup(idime,izsol) + elaup(ivect,iposi,jnode)  ! Aup
                            
                         end do
                         do jdime = 1,ndime
                            jposi = jposi0 + jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(ivect,inode,jposi)  ! Apu
                         end do
                         if( present(Q) ) then
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Q(izsol) = Q(izsol) + elmap(ivect,inode,jnode)                     ! Q
                         end if
                         if( present(App) ) then
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            App(izsol) = App(izsol) + elapp(ivect,inode,jnode)                 ! App
                         end if
                      end if
                   end do
                end do
             end if

          end do

       end if

    end if
    !
    ! RHS
    !
    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
             do idime = 1,ndime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                bu(idime,ipoin) = bu(idime,ipoin) + elrbu(ivect,idime,inode)
             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             bp(ipoin) = bp(ipoin) + elrbp(ivect,inode)
          end do
       end if
    end do

  end subroutine nsi_assembly_fractional_step

  subroutine nsi_assembly_fractional_step_scalar(&
       pnode,pevat,ielem,lnods,elvel,elpre,&
       elauu,elaup,elapp,elapu,elmap,elrbu,elrbp,Aup,Apu,&
       bu,bp,Q,Auu,App)

    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  nzdom,ndime,r_sol
    use def_domain, only       :  c_sol,lezdo
    use def_master, only       :  solve
    use def_kermod, only       :  kfl_element_to_csr

    integer(ip), intent(in)              :: pnode
    integer(ip), intent(in)              :: pevat
    integer(ip), intent(in)              :: ielem
    integer(ip), intent(in)              :: lnods(pnode)
    real(rp),    intent(in)              :: elvel(ndime,pnode,*)
    real(rp),    intent(in)              :: elpre(pnode,*)
    real(rp),    intent(in)              :: elauu(pevat,pevat)
    real(rp),    intent(in)              :: elaup(pevat,pnode)
    real(rp),    intent(in)              :: elapp(pnode,pnode)
    real(rp),    intent(in)              :: elapu(pnode,pevat)
    real(rp),    intent(in)              :: elmap(pnode,pnode)
    real(rp),    intent(inout)           :: elrbu(ndime,pnode)
    real(rp),    intent(inout)           :: elrbp(pnode)
    real(rp),    intent(inout)           :: Aup(ndime,nzdom)
    real(rp),    intent(inout)           :: Apu(ndime,nzdom)
    real(rp),    intent(inout)           :: bu(ndime,*)
    real(rp),    intent(inout)           :: bp(*)
    real(rp),    intent(inout), optional :: Q(solve(2)%nzmat)
    real(rp),    intent(inout), optional :: Auu(ndime,ndime,*)
    real(rp),    intent(inout), optional :: App(*)
    integer(ip)                          :: ndofn,inode,jnode,iposi,jposi
    integer(ip)                          :: idime,jdime,izsol,jpoin,ipoin,jcolu
    integer(ip)                          :: iposi0,jposi0,ievat,jevat


    !----------------------------------------------------------------
    !
    ! Compute bu^(e) <= bu^(e) - Auu^(e) u^(e)
    ! Compute bp^(e) <= bp^(e) - App^(e) p^{e}
    !
    !----------------------------------------------------------------

    if( .not. present(Auu) ) then
       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndime+idime
                   elrbu(idime,inode) = elrbu(idime,inode) &
                        - elauu(ievat,jevat) * elvel(jdime,jnode,1)
                end do
             end do
          end do
       end do
       do jnode = 1,pnode
          do inode = 1,pnode
             elrbp(inode) = elrbp(inode) &
                  - elapp(inode,jnode) * elpre(jnode,1)
          end do
       end do
    end if

    ndofn = ndime + 1

    !----------------------------------------------------------------
    !
    ! Assemble the following matrices:
    !
    ! Aup <= Aup^(e)
    ! Apu <= Apu^(e)
    ! bu  <= bu^(e)
    ! bp  <= bp^(e)
    !
    !----------------------------------------------------------------

    if( solve(2) % kfl_symme == 1 ) then
       !
       ! App requires symmetric assembly
       !
       call runend('SYMMETRIC ASSEMBLY NOT CODED')
    else
       !
       ! App requires asymmetric assembly
       !
       if( kfl_element_to_csr == 1 ) then

          do inode = 1,pnode
             iposi0 = (inode-1) * ndime
             do jnode = 1,pnode
                jposi0 = (jnode-1) * ndime
                izsol = lezdo(inode,jnode,ielem)
                do idime = 1,ndime
                   iposi = iposi0+idime
                   if( present(Auu) ) then
                      do jdime = 1,ndime
                         jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                              &                   + elauu(iposi,jposi)          ! Auu
                      end do
                   end if
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode) ! Aup
                end do
                do jdime=1,ndime
                   jposi = jposi0+jdime
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu
                end do
                if( present(Q) ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   Q(izsol) = Q(izsol) + elmap(inode,jnode)                 ! App
                end if
                if( present(App) ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   App(izsol) = App(izsol) + elapp(inode,jnode)              ! App
                end if
             end do
          end do

       else

          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = r_sol(ipoin)
                jcolu = c_sol(izsol)
                do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = c_sol(izsol)
                end do
                if( jcolu == jpoin ) then
                   iposi0 = (inode-1) * ndime
                   jposi0 = (jnode-1) * ndime
                   do idime = 1,ndime
                      iposi = iposi0 + idime
                      if( present(Auu) ) then
                         do jdime=1,ndime
                            jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                                 &                   + elauu(iposi,jposi)          ! Auu
                         end do
                      end if
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode) ! Aup
                   end do

                   do jdime=1,ndime
                      jposi = jposi0+jdime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu
                   end do
                   if( present(Q) ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Q(izsol) = Q(izsol) + elmap(inode,jnode)                 ! App
                   end if
                   if( present(App) ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      App(izsol) = App(izsol) + elapp(inode,jnode)              ! App
                   end if
                end if
             end do
          end do

       end if

    end if
    !
    ! RHS
    !
    do inode = 1,pnode
       ipoin = lnods(inode)
       do idime = 1,ndime
#ifdef NO_COLORING
          !$OMP ATOMIC
#endif
          bu(idime,ipoin) = bu(idime,ipoin) + elrbu(idime,inode)
       end do
#ifdef NO_COLORING
       !$OMP ATOMIC
#endif
       bp(ipoin) = bp(ipoin) + elrbp(inode)
    end do

  end subroutine nsi_assembly_fractional_step_scalar

  subroutine nsi_assembly_fractional_step_boundary_scalar(&
       pnode,pevat,ielem,lnods,elvel,elpre,elmat,elrhs,&
       Aup,Apu,bu,bp,Auu,App)

    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  nzdom,ndime,r_sol
    use def_domain, only       :  c_sol,lezdo
    use def_master, only       :  solve
    use def_kermod, only       :  kfl_element_to_csr

    integer(ip), intent(in)              :: pnode
    integer(ip), intent(in)              :: pevat
    integer(ip), intent(in)              :: ielem
    integer(ip), intent(in)              :: lnods(pnode)
    real(rp),    intent(in)              :: elvel(ndime,pnode,*)
    real(rp),    intent(in)              :: elpre(pnode,*)
    real(rp),    intent(in)              :: elmat(pevat,pevat)
    real(rp),    intent(inout)           :: elrhs(pevat)
    real(rp),    intent(inout), optional :: Aup(ndime,nzdom)
    real(rp),    intent(inout), optional :: Apu(ndime,nzdom)
    real(rp),    intent(inout), optional :: bu(ndime,*)
    real(rp),    intent(inout), optional :: bp(*)
    real(rp),    intent(inout), optional :: Auu(ndime,ndime,*)
    real(rp),    intent(inout), optional :: App(*)
    integer(ip)                          :: ndofn,inode,jnode,iposi,jposi
    integer(ip)                          :: idime,jdime,izsol,jpoin,ipoin,jcolu
    integer(ip)                          :: iposi0,jposi0,ievat,jevat

    !----------------------------------------------------------------
    !
    ! Compute bu^(e) <= bu^(e) - Auu^(e) u^(e)
    ! Compute bp^(e) <= bp^(e) - App^(e) p^{e}
    !
    !----------------------------------------------------------------

    ndofn = ndime + 1

    if( .not. present(Auu) ) then
       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndofn+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndofn+idime
                   elrhs(ievat) = elrhs(ievat) &
                        - elmat(ievat,jevat) * elvel(jdime,jnode,1)
                end do
             end do
          end do
       end do
       do jnode = 1,pnode
          jevat = jnode * ndofn
          do inode = 1,pnode
             ievat = inode * ndofn
             elrhs(ievat) = elrhs(ievat) &
                  - elmat(ievat,jevat) * elpre(jnode,1)
          end do
       end do
    end if

    !----------------------------------------------------------------
    !
    ! Assemble the following matrices:
    !
    ! Aup <= Aup^(e)
    ! Apu <= Apu^(e)
    ! bu  <= bu^(e)
    ! bp  <= bp^(e)
    !
    !----------------------------------------------------------------

    if( solve(2) % kfl_symme == 1 ) then
       !
       ! App requires symmetric assembly
       !
       call runend('SYMMETRIC ASSEMBLY NOT CODED')
    else
       !
       ! App requires asymmetric assembly
       !
       if( kfl_element_to_csr == 1 ) then

          do inode = 1,pnode
             iposi0 = (inode-1) * ndofn
             do jnode = 1,pnode
                jposi0 = (jnode-1) * ndofn
                izsol = lezdo(inode,jnode,ielem)
                do idime = 1,ndime
                   iposi = iposi0+idime
                   if( present(Auu) ) then
                      do jdime = 1,ndime
                         jposi = (jnode-1)*ndofn+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                              &                   + elmat(iposi,jposi)          ! Auu
                      end do
                   end if
                   if( present(Aup) ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Aup(idime,izsol) = Aup(idime,izsol) + elmat(iposi,jnode*ndofn) ! Aup
                   end if
                end do
                if( present(Apu) ) then
                   do jdime=1,ndime
                      jposi = jposi0+jdime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Apu(jdime,izsol) = Apu(jdime,izsol) + elmat(inode*ndofn,jposi)  ! Apu
                   end do
                end if
                if( present(App) ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   App(izsol) = App(izsol) + elmat(inode*ndofn,jnode*ndofn)  ! App
                end if
             end do
          end do

       else

          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = r_sol(ipoin)
                jcolu = c_sol(izsol)
                do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = c_sol(izsol)
                end do
                if( jcolu == jpoin ) then
                   iposi0 = (inode-1) * ndofn
                   jposi0 = (jnode-1) * ndofn
                   do idime = 1,ndime
                      iposi = iposi0 + idime
                      if( present(Auu) ) then
                         do jdime=1,ndime
                            jposi                  = (jnode-1)*ndofn+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                                 &                   + elmat(iposi,jposi)          ! Auu
                         end do
                      end if
                      if( present(Aup) ) then
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Aup(idime,izsol) = Aup(idime,izsol) + elmat(iposi,jnode*ndofn) ! Aup
                      end if
                   end do

                   if( present(Apu) ) then
                      do jdime=1,ndime
                         jposi = jposi0+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Apu(jdime,izsol) = Apu(jdime,izsol) + elmat(inode*ndofn,jposi)  ! Apu
                      end do
                   end if
                   if( present(App) ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      App(izsol) = App(izsol) + elmat(inode*ndofn,jnode*ndofn)              ! App
                   end if
                end if
             end do
          end do

       end if

    end if
    !
    ! RHS
    !
    do inode = 1,pnode
       ipoin = lnods(inode)
       do idime = 1,ndime
          iposi = (inode-1)*ndofn+idime
          if( present(bu) ) then
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             bu(idime,ipoin) = bu(idime,ipoin) + elrhs(iposi)
          end if
       end do
#ifdef NO_COLORING
       !$OMP ATOMIC
#endif
       if( present(bp) ) then
          bp(ipoin) = bp(ipoin) + elrhs(inode*ndofn)
       end if
    end do

  end subroutine nsi_assembly_fractional_step_boundary_scalar

  subroutine nsi_assembly_monolithic(&
       pnode,pevat,list_elements,lnods,elauu,elaup,elapp,elapu,&
       elrbu,elrbp,A,b)
    !-----------------------------------------------------------------------
    !****f* mathru/assma3
    ! NAME
    !    assma3
    ! DESCRIPTION
    !    Assembly an elemental matrix ELMAT in global matrix AMATR
    ! USES
    ! USED BY
    !    ***_elmope
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pevat
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elauu(VECTOR_SIZE,pevat,pevat)
    real(rp),    intent(in)    :: elaup(VECTOR_SIZE,pevat,pnode)
    real(rp),    intent(in)    :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)    :: elapu(VECTOR_SIZE,pnode,pevat)
    real(rp),    intent(in)    :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout) :: A(ndime+1,ndime+1,nzdom)
    real(rp),    intent(inout) :: b(ndime+1,*)
    integer(ip)                :: ndofn,inode,jnode
    integer(ip)                :: iposi,jposi,ielem
    integer(ip)                :: idime,jdime,izsol
    integer(ip)                :: jpoin,ipoin,jcolu
    integer(ip)                :: ivect

    ndofn = ndime + 1
    solve(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then

          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
             !
             ! Matrix
             !
             do jnode = 1,pnode
                jpoin = lnods(ivect,jnode)
                izsol = r_sol(ipoin)
                jcolu = c_sol(izsol)
                do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = c_sol(izsol)
                end do
                if( jcolu == jpoin ) then
                   do idime = 1,ndime
                      iposi = (inode-1) * ndime + idime
                      do jdime = 1,ndime
                         jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         A(jdime,idime,izsol) = A(jdime,idime,izsol) + elauu(ivect,iposi,jposi)   ! Auu
                      end do
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      A(ndofn,idime,izsol) = A(ndofn,idime,izsol) + elaup(ivect,iposi,jnode)      ! Aup
                   end do
                   do jdime=1,ndime
                      jposi = (jnode-1) * ndime + jdime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      A(jdime,ndofn,izsol) = A(jdime,ndofn,izsol) + elapu(ivect,inode,jposi)      ! Apu
                   end do
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   A(ndofn,ndofn,izsol) = A(ndofn,ndofn,izsol) + elapp(ivect,inode,jnode)         ! App
                end if
             end do
             !
             ! RHS
             !
             do idime = 1,ndime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                b(idime,ipoin) = b(idime,ipoin) + elrbu(ivect,idime,inode)
             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             b(ndofn,ipoin)   = b(ndofn,ipoin)   + elrbp(ivect,inode)
          end do
       end if
    end do

  end subroutine nsi_assembly_monolithic

  !------------------------------------------------------------------------
  !> @addtogroup NastinMatrixAssembly
  !> @{
  !> @file    nsi_assemble_schur.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Matrix and RHS assembly
  !> @details Assembly of matrix and RHS:
  !>          1. Auu, Aup, Apu, App, bu, bp
  !>          2. Only App
  !>          3. CMM consitent mass matrix
  !> @}
  !------------------------------------------------------------------------

  subroutine nsi_assembly_schur_method_scalar(&
       itask,pnode,pevat,ielem,lnods,elauu,elaup,elapp,elapu,&
       elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

    integer(ip), intent(in)    :: itask,pnode,pevat,ielem
    integer(ip), intent(in)    :: lnods(pnode)
    real(rp),    intent(in)    :: elauu(pevat,pevat)
    real(rp),    intent(in)    :: elaup(pevat,pnode)
    real(rp),    intent(in)    :: elapp(pnode,pnode)
    real(rp),    intent(in)    :: elapu(pnode,pevat)
    real(rp),    intent(in)    :: elrbu(ndime,pnode)
    real(rp),    intent(in)    :: elrbp(pnode)
    real(rp),    intent(inout) :: Auu(*)
    real(rp),    intent(inout) :: Aup(*)
    real(rp),    intent(inout) :: Apu(*)
    real(rp),    intent(inout) :: App(*)
    real(rp),    intent(inout) :: bu(*),bp(*)

    call nsi_assembly_schur_method_scalar_go(&
         itask,pnode,pevat,ielem,lnods,elauu,elaup,elapp,elapu,&
         elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

  end subroutine nsi_assembly_schur_method_scalar

  subroutine nsi_assembly_schur_method_vector(&
       itask,pnode,pevat,list_elements,lnods,elauu,elaup,elapp,elapu,&
       elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

    integer(ip), intent(in)    :: itask,pnode,pevat
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elauu(VECTOR_SIZE,pevat,pevat)
    real(rp),    intent(in)    :: elaup(VECTOR_SIZE,pevat,pnode)
    real(rp),    intent(in)    :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)    :: elapu(VECTOR_SIZE,pnode,pevat)
    real(rp),    intent(in)    :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout) :: Auu(*)
    real(rp),    intent(inout) :: Aup(*)
    real(rp),    intent(inout) :: Apu(*)
    real(rp),    intent(inout) :: App(*)
    real(rp),    intent(inout) :: bu(*),bp(*)

!    integer(ip)                :: lnods_tmp(pnode),ivect,ielem
!    real(rp)                   :: elauu_tmp(pevat,pevat)
!    real(rp)                   :: elaup_tmp(pevat,pnode)
!    real(rp)                   :: elapp_tmp(pnode,pnode)
!    real(rp)                   :: elapu_tmp(pnode,pevat)
!    real(rp)                   :: elrbu_tmp(ndime,pnode)
!    real(rp)                   :: elrbp_tmp(pnode)

    if( itask == 1 ) then
       !do ivect = 1,VECTOR_SIZE
       !   ielem = list_elements(ivect)
       !   if( ielem > 0 ) then
       !      lnods_tmp(:)   = lnods(ivect,:)
       !      elauu_tmp(:,:) = elauu(ivect,:,:)
       !      elaup_tmp(:,:) = elaup(ivect,:,:)
       !      elapp_tmp(:,:) = elapp(ivect,:,:)
       !      elapu_tmp(:,:) = elapu(ivect,:,:)
       !      elrbu_tmp(:,:) = elrbu(ivect,:,:)
       !      elrbp_tmp(:)   = elrbp(ivect,:)
       !      call nsi_assembly_schur_method_scalar_go(&
       !           itask,pnode,pevat,ielem,lnods_tmp,elauu_tmp,elaup_tmp,elapp_tmp,elapu_tmp,&
       !           elrbu_tmp,elrbp_tmp,Auu,Aup,App,Apu,bu,bp)
       !   end if
       !end do
       call nsi_assembly_schur_method_vector_go(&
            itask,pnode,pevat,list_elements,lnods,elauu,elaup,elapp,elapu,&
            elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

    else

       call nsi_assembly_schur_method_vector_go(&
            itask,pnode,pevat,list_elements,lnods,elauu,elaup,elapp,elapu,&
            elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

       !do ivect = 1,VECTOR_SIZE
       !   ielem             = list_elements(ivect)
       !   if( ielem > 0 ) then
       !      lnods_tmp(:)   = lnods(ivect,:)
       !      elapp_tmp(:,:) = elapp(ivect,:,:)
       !      call nsi_assembly_schur_method_scalar_go(&
       !           itask,pnode,pevat,ielem,lnods_tmp,elauu_tmp,elaup_tmp,elapp_tmp,elapu_tmp,&
       !           elrbu_tmp,elrbp_tmp,Auu,Aup,App,Apu,bu,bp)
       !   end if
       !end do
    end if

    !call nsi_assembly_schur_method_vector_go(&
    !   itask,pnode,pevat,list_elements,lnods,elauu,elaup,elapp,elapu,&
    !   elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

  end subroutine nsi_assembly_schur_method_vector

  !
  ! In Finite volume
  ! PNODE          =  2
  ! LNODS(1:PNODE) =  nodes of the edge
  ! LEZDO          <= fv_face_graph
  ! IELEM          <= IFACG
  ! r_sol,c_sol    <= r_elm_2, c_elm_2
  !
  subroutine nsi_assembly_schur_method_scalar_go(&
       itask,pnode,pevat,ielem,lnods,elauu,elaup,elapp,elapu,&
       elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)

    implicit none
    integer(ip), intent(in)    :: itask,pnode,pevat,ielem
    integer(ip), intent(in)    :: lnods(pnode)
    real(rp),    intent(in)    :: elauu(pevat,pevat)
    real(rp),    intent(in)    :: elaup(pevat,pnode)
    real(rp),    intent(in)    :: elapp(pnode,pnode)
    real(rp),    intent(in)    :: elapu(pnode,pevat)
    real(rp),    intent(in)    :: elrbu(ndime,pnode)
    real(rp),    intent(in)    :: elrbp(pnode)
    real(rp),    intent(inout) :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(inout) :: Aup(ndime,nzdom)
    real(rp),    intent(inout) :: Apu(ndime,nzdom)
    real(rp),    intent(inout) :: App(solve(2)%nzmat)
    real(rp),    intent(inout) :: bu(ndime,*),bp(*)
    integer(ip)                :: ndofn,inode,jnode,iposi,jposi,izsym
    integer(ip)                :: idime,jdime,izsol,jpoin,ipoin,jcolu

    ndofn = ndime + 1

    select case(itask)

    case (1_ip)

       !----------------------------------------------------------------
       !
       ! Assemble the following matrices:
       !
       ! Auu <= Auu^(e)
       ! Aup <= Aup^(e)
       ! Apu <= Apu^(e)
       ! App <= App^(e)
       ! bu  <= bu^(e)
       ! bp  <= bp^(e)
       !
       !----------------------------------------------------------------

       solve(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED

       if( solve(2) % kfl_symme == 1 ) then
          !
          ! App requires symmetric assembly
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = r_sol(ipoin)
                jcolu = c_sol(izsol)
                do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = c_sol(izsol)
                end do
                if( jcolu == jpoin ) then
                   do idime=1,ndime
                      iposi=(inode-1)*ndime+idime
                      do jdime=1,ndime
                         jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                              &                   + elauu(iposi,jposi)          ! Auu
                      end do
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode)  ! Aup
                   end do
                   iposi = iposi+1
                   do jdime=1,ndime
                      jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu
                   end do
                end if
                if( jpoin <= ipoin ) then
                   izsym = r_sym(ipoin)
                   jcolu = c_sym(izsym)
                   do while( jcolu /= jpoin .and. izsym < r_sym(ipoin+1)-1 )
                      izsym = izsym + 1
                      jcolu = c_sym(izsym)
                   end do
                   if( jcolu == jpoin ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      App(izsym) = App(izsym) + elapp(inode,jnode)              ! App
                   end if
                end if
             end do
          end do
       else
          !
          ! App requires assymmetric assembly
          !
          if( kfl_element_to_csr == 1 ) then
             do inode = 1,pnode
                do jnode = 1,pnode
                   izsol = lezdo(inode,jnode,ielem)
                   do idime = 1,ndime
                      iposi = (inode-1)*ndime+idime
                      do jdime = 1,ndime
                         jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol)&
                              +                   elauu(iposi,jposi)           ! Auu
                      end do
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode) ! Aup
                   end do
                   do jdime=1,ndime
                      jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu
                   end do
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   App(izsol) = App(izsol) + elapp(inode,jnode)                 ! App
                end do
             end do
          else
             do inode = 1,pnode
                ipoin = lnods(inode)
                do jnode = 1,pnode
                   jpoin = lnods(jnode)
                   izsol = r_sol(ipoin)
                   jcolu = c_sol(izsol)
                   do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                      izsol = izsol + 1
                      jcolu = c_sol(izsol)
                   end do
                   if( jcolu == jpoin ) then
                      do idime = 1,ndime
                         iposi = (inode-1)*ndime+idime
                         do jdime = 1,ndime
                            jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol)&
                                 +                   elauu(iposi,jposi)           ! Auu
                         end do
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode) ! Aup
                      end do

                      do jdime=1,ndime
                         jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu
                      end do
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      App(izsol) = App(izsol) + elapp(inode,jnode)                 ! App
                   end if
                end do
             end do
          end if
       end if
       !
       ! RHS
       !
       do inode = 1,pnode
          ipoin = lnods(inode)
          do idime = 1,ndime
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             bu(idime,ipoin) = bu(idime,ipoin) + elrbu(idime,inode)
          end do
#ifdef NO_COLORING
          !$OMP ATOMIC
#endif
          bp(ipoin) = bp(ipoin) + elrbp(inode)
       end do

    case (2_ip)

       !----------------------------------------------------------------
       !
       ! Assemble only (App is pressure Schur preconditioner Q):
       ! App <= App^(e)
       !
       !----------------------------------------------------------------

       if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
          continue
       else
          solve(2) % kfl_assem = SOL_MATRIX_HAS_CHANGED
       end if
       
       if( solve(2) % kfl_symme == 1 ) then
          !
          ! App requires symmetric assembly
          !
          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                izsol = r_sol(ipoin)
                jcolu = c_sol(izsol)
                do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                   izsol = izsol + 1
                   jcolu = c_sol(izsol)
                end do
                if( jpoin <= ipoin ) then
                   izsym = r_sym(ipoin)
                   jcolu = c_sym(izsym)
                   do while( jcolu /= jpoin .and. izsym < r_sym(ipoin+1)-1 )
                      izsym = izsym + 1
                      jcolu = c_sym(izsym)
                   end do
                   if( jcolu == jpoin ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      App(izsym) = App(izsym) + elapp(inode,jnode) ! App
                   end if
                end if
             end do
          end do
       else
          !
          ! App requires asymmetric assembly
          !
          if( kfl_element_to_csr == 1 ) then
             do inode = 1,pnode
                do jnode = 1,pnode
                   izsol = lezdo(inode,jnode,ielem)
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   App(izsol) = App(izsol) + elapp(inode,jnode) ! App
                end do
             end do
          else
             do inode = 1,pnode
                ipoin = lnods(inode)
                do jnode = 1,pnode
                   jpoin = lnods(jnode)
                   izsol = r_sol(ipoin)
                   jcolu = c_sol(izsol)
                   do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                      izsol = izsol + 1
                      jcolu = c_sol(izsol)
                   end do
                   if( jcolu == jpoin ) then
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      App(izsol) = App(izsol) + elapp(inode,jnode) ! App
                   end if
                end do
             end do
          end if

       end if

    case (3_ip)

       !----------------------------------------------------------------
       !
       ! Assemble only Auu <= Auu^(e)  - Actually I will use it to assemble the consistent mass matrix that has the sam shape
       !
       !----------------------------------------------------------------

       solve(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED

       do inode = 1,pnode
          ipoin = lnods(inode)
          do jnode = 1,pnode
             jpoin = lnods(jnode)
             izsol = r_sol(ipoin)
             jcolu = c_sol(izsol)
             do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                izsol = izsol + 1
                jcolu = c_sol(izsol)
             end do
             if( jcolu == jpoin ) then
                do idime=1,ndime
                   iposi=(inode-1)*ndime+idime
                   do jdime=1,ndime
                      jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                           &                   + elauu(iposi,jposi)          ! Auu
                   end do
                end do

             end if
          end do
       end do

    end select

  end subroutine nsi_assembly_schur_method_scalar_go

  subroutine nsi_assembly_schur_method_vector_go(&
       itask,pnode,pevat,list_elements,lnods,elauu,elaup,elapp,elapu,&
       elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)
    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  nzdom,ndime,r_sol
    use def_domain, only       :  r_sym,c_sym,c_sol,lezdo
    use def_master, only       :  solve
    use def_kermod, only       :  kfl_element_to_csr

    integer(ip), intent(in)    :: itask,pnode,pevat,list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elauu(VECTOR_SIZE,pevat,pevat)
    real(rp),    intent(in)    :: elaup(VECTOR_SIZE,pevat,pnode)
    real(rp),    intent(in)    :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)    :: elapu(VECTOR_SIZE,pnode,pevat)
    real(rp),    intent(in)    :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout) :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(inout) :: Aup(ndime,nzdom)
    real(rp),    intent(inout) :: Apu(ndime,nzdom)
    real(rp),    intent(inout) :: App(solve(2)%nzmat)
    real(rp),    intent(inout) :: bu(ndime,*),bp(*)
    integer(ip)                :: ndofn,inode,jnode,iposi,jposi,izsym,ielem
    integer(ip)                :: idime,jdime,izsol,jpoin,ipoin,jcolu,ivect
    integer(ip)                :: iposi0,jposi0

    ndofn = ndime + 1

    if( itask == 1 ) then

       !----------------------------------------------------------------
       !
       ! Assemble the following matrices:
       !
       ! Auu <= Auu^(e)
       ! Aup <= Aup^(e)
       ! Apu <= Apu^(e)
       ! App <= App^(e)
       ! bu  <= bu^(e)
       ! bp  <= bp^(e)
       !
       !----------------------------------------------------------------

       solve(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED
       
       if( solve(2) % kfl_symme == 1 ) then
          !
          ! App requires symmetric assembly
          !
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then

                do inode = 1,pnode
                   ipoin = lnods(ivect,inode)
                   do jnode = 1,pnode
                      jpoin = lnods(ivect,jnode)
                      izsol = r_sol(ipoin)
                      jcolu = c_sol(izsol)
                      do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                         izsol = izsol + 1
                         jcolu = c_sol(izsol)
                      end do
                      if( jcolu == jpoin ) then
                         do idime=1,ndime
                            iposi=(inode-1)*ndime+idime
                            do jdime=1,ndime
                               jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                               !$OMP ATOMIC
#endif
                               Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                                    &                   + elauu(ivect,iposi,jposi)          ! Auu
                            end do
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Aup(idime,izsol) = Aup(idime,izsol) + elaup(ivect,iposi,jnode)  ! Aup
                         end do
                         iposi = iposi+1
                         do jdime=1,ndime
                            jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(ivect,inode,jposi)  ! Apu
                         end do
                      end if
                      if( jpoin <= ipoin ) then
                         izsym = r_sym(ipoin)
                         jcolu = c_sym(izsym)
                         do while( jcolu /= jpoin .and. izsym < r_sym(ipoin+1)-1 )
                            izsym = izsym + 1
                            jcolu = c_sym(izsym)
                         end do
                         if( jcolu == jpoin ) then
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            App(izsym) = App(izsym) + elapp(ivect,inode,jnode)              ! App
                         end if
                      end if
                   end do
                end do
             end if
          end do

       else
          !
          ! App requires asymmetric assembly
          !
          if( kfl_element_to_csr == 1 ) then

             do ivect = 1,VECTOR_SIZE
                ielem = list_elements(ivect)
                if( ielem > 0 ) then

                   do inode = 1,pnode
                      iposi0 = (inode-1) * ndime
                      do jnode = 1,pnode
                         jposi0 = (jnode-1) * ndime
                         izsol = lezdo(inode,jnode,ielem)
                         do idime = 1,ndime
                            iposi = iposi0+idime
                            do jdime = 1,ndime
                               jposi                  = jposi0+jdime
#ifdef NO_COLORING
                               !$OMP ATOMIC
#endif
                               Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol)&
                                    +                   elauu(ivect,iposi,jposi)           ! Auu
                            end do
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Aup(idime,izsol) = Aup(idime,izsol) + elaup(ivect,iposi,jnode) ! Aup
                         end do
                         do jdime=1,ndime
                            jposi = jposi0+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(ivect,inode,jposi)  ! Apu
                         end do
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         App(izsol) = App(izsol) + elapp(ivect,inode,jnode)                 ! App
                      end do
                   end do

                end if
             end do

          else

             do ivect = 1,VECTOR_SIZE
                ielem = list_elements(ivect)
                if( ielem > 0 ) then

                   do inode = 1,pnode
                      ipoin = lnods(ivect,inode)
                      do jnode = 1,pnode
                         jpoin = lnods(ivect,jnode)
                         izsol = r_sol(ipoin)
                         jcolu = c_sol(izsol)
                         do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                            izsol = izsol + 1
                            jcolu = c_sol(izsol)
                         end do
                         if( jcolu == jpoin ) then
                            iposi0 = (inode-1) * ndime
                            jposi0 = (jnode-1) * ndime
                            do idime = 1,ndime
                               iposi = iposi0 + idime
                               do jdime = 1,ndime
                                  jposi                  = jposi0+jdime
#ifdef NO_COLORING
                                  !$OMP ATOMIC
#endif
                                  Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol)&
                                       +                   elauu(ivect,iposi,jposi)           ! Auu
                               end do
#ifdef NO_COLORING
                               !$OMP ATOMIC
#endif
                               Aup(idime,izsol) = Aup(idime,izsol) + elaup(ivect,iposi,jnode) ! Aup
                            end do

                            do jdime=1,ndime
                               jposi = jposi0+jdime
#ifdef NO_COLORING
                               !$OMP ATOMIC
#endif
                               Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(ivect,inode,jposi)  ! Apu
                            end do
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            App(izsol) = App(izsol) + elapp(ivect,inode,jnode)                 ! App
                         end if
                      end do
                   end do
                end if

             end do

          end if

       end if
       !
       ! RHS
       !
       do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                do idime = 1,ndime
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   bu(idime,ipoin) = bu(idime,ipoin) + elrbu(ivect,idime,inode)
                end do
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                bp(ipoin) = bp(ipoin) + elrbp(ivect,inode)
             end do
          end if
       end do

    else if( itask == 2 ) then

       !----------------------------------------------------------------
       !
       ! Assemble only (App is pressure Schur preconditioner):
       ! App <= App^(e)
       !
       !----------------------------------------------------------------

       if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
          continue
       else
          solve(2) % kfl_assem = SOL_MATRIX_HAS_CHANGED
       end if
       
       if( solve(2) % kfl_symme == 1 ) then
          !
          ! App requires symmetric assembly
          !
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                do inode = 1,pnode
                   ipoin = lnods(ivect,inode)
                   do jnode = 1,pnode
                      jpoin = lnods(ivect,jnode)
                      izsol = r_sol(ipoin)
                      jcolu = c_sol(izsol)
                      do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                         izsol = izsol + 1
                         jcolu = c_sol(izsol)
                      end do
                      if( jpoin <= ipoin ) then
                         izsym = r_sym(ipoin)
                         jcolu = c_sym(izsym)
                         do while( jcolu /= jpoin .and. izsym < r_sym(ipoin+1)-1 )
                            izsym = izsym + 1
                            jcolu = c_sym(izsym)
                         end do
                         if( jcolu == jpoin ) then
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            App(izsym) = App(izsym) + elapp(ivect,inode,jnode) ! App
                         end if
                      end if
                   end do
                end do
             end if
          end do

       else
          !
          ! App requires asymmetric assembly
          !
          if( kfl_element_to_csr == 1 ) then

             do ivect = 1,VECTOR_SIZE
                ielem = list_elements(ivect)
                if( ielem > 0 ) then

                   do inode = 1,pnode
                      do jnode = 1,pnode
                         izsol = lezdo(inode,jnode,ielem)
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         App(izsol) = App(izsol) + elapp(ivect,inode,jnode) ! App
                      end do
                   end do
                end if
             end do

          else

             do ivect = 1,VECTOR_SIZE
                ielem = list_elements(ivect)
                if( ielem > 0 ) then

                   do inode = 1,pnode
                      ipoin = lnods(ivect,inode)
                      do jnode = 1,pnode
                         jpoin = lnods(ivect,jnode)
                         izsol = r_sol(ipoin)
                         jcolu = c_sol(izsol)
                         do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                            izsol = izsol + 1
                            jcolu = c_sol(izsol)
                         end do
                         if( jcolu == jpoin ) then
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            App(izsol) = App(izsol) + elapp(ivect,inode,jnode) ! App
                         end if
                      end do
                   end do

                end if
             end do

          end if

       end if

    else if( itask == 3 ) then

       !----------------------------------------------------------------
       !
       ! Assemble only Auu <= Auu^(e)  - Actually I will use it to assemble the consistent mass matrix that has the sam shape
       !
       !----------------------------------------------------------------

       do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then

             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                do jnode = 1,pnode
                   jpoin = lnods(ivect,jnode)
                   izsol = r_sol(ipoin)
                   jcolu = c_sol(izsol)
                   do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                      izsol = izsol + 1
                      jcolu = c_sol(izsol)
                   end do
                   if( jcolu == jpoin ) then
                      do idime=1,ndime
                         iposi=(inode-1)*ndime+idime
                         do jdime=1,ndime
                            jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                                 &                   + elauu(ivect,iposi,jposi)          ! Auu
                         end do
                      end do

                   end if
                end do
             end do

          end if
       end do

    end if

  end subroutine nsi_assembly_schur_method_vector_go

  !------------------------------------------------------------------------
  !> @addtogroup NastinMatrixAssembly
  !> @{
  !> @file    nsi_assembly_dt_rho_tau.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Assemble projections
  !> @details Assembly of following projections:
  !>          DT_RHO_NSI = Projection of dt / rho
  !>          MASS_RHO_NSI = Projection rho
  !>          TAU_NSI    = Projection of tau
  !> @}
  !------------------------------------------------------------------------

  subroutine nsi_assembly_dt_rho_tau_nu(&
       pnode,list_elements,lnods,eldtrho,eltau,elmurho,&
       dt_rho_nsi,tau_nsi,mass_rho_nsi)

    integer(ip), intent(in)           :: pnode
    integer(ip), intent(in)           :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)           :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)           :: eldtrho(VECTOR_SIZE,pnode)
    real(rp),    intent(in)           :: elmurho(VECTOR_SIZE,pnode)
    real(rp),    intent(in)           :: eltau(VECTOR_SIZE,pnode)
    real(rp),    intent(inout), pointer :: dt_rho_nsi(:)
    real(rp),    intent(inout), pointer :: mass_rho_nsi(:,:)
    real(rp),    intent(inout), pointer :: tau_nsi(:)
    integer(ip)                       :: ivect,ielem,inode,ipoin

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             mass_rho_nsi(ipoin,1) = mass_rho_nsi(ipoin,1) + elmurho(ivect,inode)

             tau_nsi(ipoin) = tau_nsi(ipoin) + eltau(ivect,inode)
          end do
       end if
    end do

  end subroutine nsi_assembly_dt_rho_tau_nu

  !-----------------------------------------------------------------------
  !
  !> @date    21/02/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Scaling of pressure equation
  !> @details Scaling of pressure equation. This is useful
  !>          in order to solve the pressure Schur complement
  !>          solver without changing the Schur complement
  !>          preconditioner:
  !>
  !>          Q = Scal L Dp = rc
  !>          where (Scal L) approximates the App + Uzawa operator.
  !>
  !>          Last system is then
  !>
  !>          L Dp = Scal^{-1} rc
  !>
  !>          with
  !>
  !>          KFL_PREDI_NSI == 8 ... Scal = 1.0_rp / ( 1.0_rp / Tau + pabdf_nsi(1) / Tim ) + Tau
  !>          KFL_PREDI_NSI == 9 ... Scal = Tim / pabdf_nsi(1) + Tim
  !>
  !
  !-----------------------------------------------------------------------

  subroutine nsi_assembly_scaling_pressure_equation(Apu,App,bp)

    use def_kintyp,               only : ip,rp
    use def_master,               only : INOTMASTER
    use def_master,               only : modul,mem_modul
    use mod_memory,               only : memory_alloca
    use mod_memory,               only : memory_deallo
    use def_domain,               only : npoin
    use def_domain,               only : r_dom,nzdom
    use def_nastin,               only : kfl_predi_nsi
    use def_nastin,               only : NSI_SCHUR_COMPLEMENT
    use def_nastin,               only : dt_rho_nsi
    use def_nastin,               only : tau_nsi

    use mod_solver

    implicit none

    real(rp),    intent(inout)         :: Apu(ndime,nzdom)
    real(rp),    intent(inout)         :: App(nzdom)
    real(rp),    intent(inout)         :: bp(npoin)

    integer(ip)                        :: ipoin,izdom
    real(rp),    pointer               :: Tim(:)
    real(rp),    pointer               :: Tau(:)
    real(rp),    pointer               :: Scal(:)

    if( NSI_SCHUR_COMPLEMENT .and. INOTMASTER .and. &
         ( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) ) then

       nullify(Scal)
       call memory_alloca(mem_modul(1:2,modul),'SCAL','nsi_assembly_scaling_pressure_equation',Scal,npoin)
       Tim  => dt_rho_nsi
       Tau  => tau_nsi
       !
       ! Schur complement approximation
       !
       if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 ) then

          Scal =  1.0_rp / ( 1.0_rp / Tau + pabdf_nsi(1) / Tim ) + Tau

       else if( kfl_predi_nsi == 9 ) then

          Scal =  Tim / pabdf_nsi(1) + Tim

       end if
       !
       ! bp  <= Scal^{-1} bp
       ! App <= Scal^{-1} App
       ! Apu <= Scal^{-1} Apu
       !
       do ipoin = 1,npoin
          bp(ipoin) = bp(ipoin)          / Scal(ipoin)
          do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
             App(izdom)   = App(izdom)   / Scal(ipoin)
             Apu(:,izdom) = Apu(:,izdom) / Scal(ipoin)
          end do
       end do
       
       call memory_deallo(mem_modul(1:2,modul),'SCAL','nsi_assembly_scaling_pressure_equation',Scal)

    end if

  end subroutine nsi_assembly_scaling_pressure_equation

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    6/10/2016
  !> @brief   Fractional step assembly
  !> @details Send some terms to the right and assemble matrices
  !>          for the fractional step method
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_assembly_fractional_step_rhs(&
       pnode,pevat,list_elements,lnods,elvel,elpre,&
       elauu,elaup,elapp,elapu,elmap,elrbu,elrbp,bu,bp,Q)

    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime,r_sol
    use def_domain, only       :  c_sol,lezdo
    use def_domain, only       :  lpoty,exnor,skcos
    use def_master, only       :  solve
    use def_kermod, only       :  kfl_element_to_csr
    use def_nastin, only       :  kfl_grad_div_nsi
    use def_nastin, only       :  kfl_fixrs_nsi
    use def_nastin, only       :  kfl_matdi_nsi
    use def_nastin, only       :  kfl_local_nsi
    use def_nastin, only       :  skcos_nsi

    integer(ip), intent(in)              :: pnode
    integer(ip), intent(in)              :: pevat
    integer(ip), intent(in)              :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)              :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)              :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)              :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)              :: elauu(VECTOR_SIZE,pevat,pevat)
    real(rp),    intent(in)              :: elaup(VECTOR_SIZE,pevat,pnode)
    real(rp),    intent(in)              :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)              :: elapu(VECTOR_SIZE,pnode,pevat)
    real(rp),    intent(in)              :: elmap(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(inout)           :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(inout)           :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout)           :: bu(ndime,*)
    real(rp),    intent(inout)           :: bp(*)
    real(rp),    intent(inout), optional :: Q(solve(2)%nzmat)
    integer(ip)                          :: inode,jnode,ielem
    integer(ip)                          :: idime,jdime,izsol,jpoin,ipoin,jcolu,ivect
    integer(ip)                          :: ievat,jevat,iroty,ibopo
    real(rp)                             :: elvel_tmp(1:VECTOR_SIZE,ndime,pnode)
    real(rp)                             :: worve(3),worma(3)

    !----------------------------------------------------------------
    !
    ! Compute bu^(e) <= bu^(e) - Auu^(e) u^(e)
    ! Compute bp^(e) <= bp^(e) - App^(e) p^{e}
    !
    !----------------------------------------------------------------

    elvel_tmp(1:VECTOR_SIZE,1:ndime,1:pnode) = elvel(1:VECTOR_SIZE,1:ndime,1:pnode,1)

    if( kfl_local_nsi /= 0 .and. kfl_matdi_nsi == NSI_DIRICHLET_ELEMENT ) then !.and. ( .not. NSI_FRACTIONAL_STEP ) ) then
       !
       ! Rotate the velocity to local system, as Auu has been rotated previously
       !
       do ivect = 1,VECTOR_SIZE
          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
             ibopo = lpoty(ipoin)
             if( ibopo > 0 ) then
                iroty = kfl_fixrs_nsi(ipoin)

                if( iroty == -1 ) then
                   !
                   ! Boundary conditions in the tangent skew system
                   !
                   worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                   call mbvatb(worma,exnor(1,1,ibopo),worve,ndime,ndime)
                   elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                else if( iroty == -2 ) then
                   !
                   ! Boundary conditions in the NSI tangent skew system
                   !
                   worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                   call mbvatb(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
                   elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                else if( iroty == -3 ) then
                   !
                   ! Boundary conditions in geometrical system
                   !
                   worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                   call mbvatb(worma,skcos(1,1,ibopo),worve,ndime,ndime)
                   elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                else if( iroty >= 1 ) then
                   !
                   ! Boundary conditions in a given skew system
                   !
                   worve(1:ndime)=elvel_tmp(ivect,1:ndime,inode)
                   call mbvatb(worma,skcos(1,1,iroty),worve,ndime,ndime)
                   elvel_tmp(ivect,1:ndime,inode)=worma(1:ndime)

                end if
             end if
          end do
       end do
    end if

    do jnode = 1,pnode
       do jdime = 1,ndime
          jevat = (jnode-1)*ndime+jdime
          do inode = 1,pnode
             do idime = 1,ndime
                ievat = (inode-1)*ndime+idime
                elrbu(1:VECTOR_SIZE,idime,inode) = elrbu(1:VECTOR_SIZE,idime,inode) &
                     - elauu(1:VECTOR_SIZE,ievat,jevat) * elvel_tmp(1:VECTOR_SIZE,jdime,jnode)
             end do
          end do
       end do
    end do
    do jnode = 1,pnode
       do inode = 1,pnode
          elrbp(1:VECTOR_SIZE,inode) = elrbp(1:VECTOR_SIZE,inode) &
               - elapp(1:VECTOR_SIZE,inode,jnode) * elpre(1:VECTOR_SIZE,jnode,1)
       end do
    end do

    !if( present(Q) ) solve(2) % kfl_assem = SOL_MATRIX_HAS_CHANGED

    if( present(Q) ) then

       !----------------------------------------------------------------
       !
       ! Assemble the following matrices:
       !
       ! Aup <= Aup^(e)
       ! Apu <= Apu^(e)
       ! bu  <= bu^(e)
       ! bp  <= bp^(e)
       !
       !----------------------------------------------------------------

       if( kfl_grad_div_nsi /= 2 ) then

          if( solve(2) % kfl_symme == 1 ) then
             !
             ! App requires symmetric assembly
             !
             call runend('SYMMETRIC ASSEMBLY NOT CODED')
          else
             !
             ! App requires asymmetric assembly
             !
             if( kfl_element_to_csr == 1 ) then

                do ivect = 1,VECTOR_SIZE
                   ielem = list_elements(ivect)
                   if( ielem > 0 ) then

                      do inode = 1,pnode
                         do jnode = 1,pnode
                            izsol = lezdo(inode,jnode,ielem)
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Q(izsol) = Q(izsol) + elmap(ivect,inode,jnode)
                         end do
                      end do

                   end if
                end do

             else

                do ivect = 1,VECTOR_SIZE
                   ielem = list_elements(ivect)

                   if( ielem > 0 ) then

                      do inode = 1,pnode
                         ipoin = lnods(ivect,inode)
                         do jnode = 1,pnode
                            jpoin = lnods(ivect,jnode)
                            izsol = r_sol(ipoin)
                            jcolu = c_sol(izsol)
                            do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                               izsol = izsol + 1
                               jcolu = c_sol(izsol)
                            end do
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Q(izsol) = Q(izsol) + elmap(ivect,inode,jnode) ! Q
                         end do
                      end do
                   end if

                end do

             end if

          end if
       end if
    end if
    !
    ! RHS
    !
    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
             do idime = 1,ndime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                bu(idime,ipoin) = bu(idime,ipoin) + elrbu(ivect,idime,inode)
             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             bp(ipoin) = bp(ipoin) + elrbp(ivect,inode)
          end do
       end if
    end do

  end subroutine nsi_assembly_fractional_step_rhs

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-22
  !> @brief   Semi-implicit
  !> @details Assembly of semi-implicit method
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_assembly_semi_implicit_method(&
       pnode,list_elements,lnods,elauu,Auu)

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(inout) :: Auu(ndime,ndime,*)
    integer(ip)                :: ielem,ivect,idime,izsol,jposi,jdime,iposi
    integer(ip)                :: inode,jnode,jcolu,ipoin,jpoin

    if( solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly == SOL_YES ) then

       if( kfl_element_to_csr == 1 ) then
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                do inode = 1,pnode
                   do jnode = 1,pnode
                      izsol = lezdo(inode,jnode,ielem)
                      do idime = 1,ndime
                         iposi = (inode-1)*ndime+idime
                         do jdime = 1,ndime
                            jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                            !$OMP ATOMIC
#endif
                            Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) + elauu(ivect,iposi,jposi)
                         end do
                      end do
                   end do
                end do
             end if
          end do
       else
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                do inode = 1,pnode
                   ipoin = lnods(ivect,inode)
                   do jnode = 1,pnode
                      jpoin = lnods(ivect,jnode)
                      izsol = r_dom(ipoin)
                      jcolu = c_dom(izsol)
                      do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                         izsol = izsol + 1
                         jcolu = c_sol(izsol)
                      end do
                      if( jcolu == jpoin ) then
                         do idime = 1,ndime
                            iposi = (inode-1)*ndime+idime
                            do jdime = 1,ndime
                               jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                               !$OMP ATOMIC
#endif
                               Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) + elauu(ivect,iposi,jposi)
                               !if(idime/=jdime) then
                               !   if( abs(elauu(ivect,iposi,jposi)) > 1.0e-16_rp ) print*,'CACA=',Auu(jdime,idime,izsol),elauu(ivect,iposi,jposi)
                               !   if( abs(Auu(jdime,idime,izsol)) > 1.0e-16_rp )   print*,'CACA=',Auu(jdime,idime,izsol),elauu(ivect,iposi,jposi)
                               !end if
                            end do
                         end do
                      end if
                   end do
                end do
             end if
          end do
       end if
    end if

  end subroutine nsi_assembly_semi_implicit_method
  
end module mod_nsi_assembly_global_system
!> @}
