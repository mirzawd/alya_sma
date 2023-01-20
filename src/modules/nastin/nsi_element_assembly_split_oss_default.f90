!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_element_assembly_split_oss_default(&
       vector_size_2, kfl_lumped, ndime,&
       mnode, ntens, kfl_stabi_nsi,&
       fvins_nsi, kfl_regim_nsi,&
       kfl_press_nsi,&
       kfl_linea_nsi, pabdf_nsi, nbdfp_nsi,&
       kfl_sgsti_nsi, kfl_nota1_nsi, kfl_limit_nsi,&
       penal_nsi, kfl_convection_type_nsi, &
       NSI_GALERKIN, NSI_ALGEBRAIC_SPLIT_OSS, &
       NSI_FRACTIONAL_STEP_int, &
       NSI_CONVECTION_CONSERVATIVE, NSI_CONVECTION_SKEW, &
       NSI_CONVECTION_EMAC, &
       kfl_noslw_ker, &
       pnode,pgaus,gpden_nsi,gpvis_nsi,gpvis_nsw,gppor,gpsp1,gpsp2,gpvol,   &
       gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,gprhc,gpvel,   &
       gpgve,gpsgs,elvel,elpre,elbub,elauu,elaup,elapp,   &
       elapu,elrbu,elrbp,dtinv_loc,dtsgs,pbubl,           &
       gpsha_bub,gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq, &
       elrbq,densi,elavv,elporrbu)
   

#include "def_vector_size.inc"
    use def_kintyp,         only :  ip,rp
    use def_nastin,         only :  kfl_surte_nsi
    use mod_memory,         only :  memory_alloca
    use mod_memory,         only :  memory_deallo
    use def_master,         only :  mem_modul, modul

#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    implicit none

    !----sauvageons
    integer(ip), intent(in)    ::   VECTOR_SIZE_2
    integer(ip), intent(in)    ::   kfl_lumped
    integer(ip), intent(in)    ::   ndime,mnode,ntens
    integer(ip), intent(in)    ::   kfl_stabi_nsi
    real(rp),    intent(in)    ::   fvins_nsi
    integer(ip), intent(in)    ::   kfl_regim_nsi
    integer(ip), intent(in)    ::   kfl_press_nsi
    integer(ip), intent(in)    ::   kfl_linea_nsi
    real(rp),    intent(in)    ::   pabdf_nsi(10)
    integer(ip), intent(in)    ::   nbdfp_nsi,kfl_sgsti_nsi
    integer(ip), intent(in)    ::   kfl_nota1_nsi,kfl_limit_nsi
    real(rp),    intent(in)    ::   penal_nsi
    integer(ip), intent(in)    ::   kfl_convection_type_nsi
    integer(ip), intent(in)    ::   NSI_GALERKIN
    integer(ip), intent(in)    ::   NSI_ALGEBRAIC_SPLIT_OSS
    integer(ip), intent(in)    ::   NSI_FRACTIONAL_STEP_int
    integer(ip), intent(in)    ::   NSI_CONVECTION_CONSERVATIVE
    integer(ip), intent(in)    ::   NSI_CONVECTION_SKEW
    integer(ip), intent(in)    ::   NSI_CONVECTION_EMAC
    integer(ip), intent(in)    ::   kfl_noslw_ker
    !----premiers nés
    integer(ip), intent(in)    :: pnode,pgaus
    real(rp),    intent(in)    :: gpden_nsi(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis_nsi(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis_nsw(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpgrp(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvel(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)     !< Velocity gradient
    real(rp),    intent(in)    :: gpsgs(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)    :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    ! Matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elporrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    ! Others
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    real(rp),    intent(in)    :: dtsgs(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: densi(VECTOR_SIZE,pgaus,nbdfp_nsi)
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)
    ! Local arrays
    real(rp)                   :: gpden(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpvis(VECTOR_SIZE,pgaus)
    real(rp)                   :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)                   :: agrau(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: gpsp1_p(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpsp1_v(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpsp2_v(VECTOR_SIZE,pgaus)
    real(rp)                   :: c1(VECTOR_SIZE)
    real(rp)                   :: c2(VECTOR_SIZE)
    real(rp)                   :: c3(VECTOR_SIZE)
    real(rp)                   :: c4(VECTOR_SIZE)
    real(rp)                   :: alpha(VECTOR_SIZE)
    real(rp)                   :: beta(VECTOR_SIZE)
    real(rp)                   :: fact0(VECTOR_SIZE)
    real(rp)                   :: fact1(VECTOR_SIZE)
    real(rp)                   :: fact2(VECTOR_SIZE)
    real(rp)                   :: fact3(VECTOR_SIZE)
    real(rp)                   :: fact4(VECTOR_SIZE)
    real(rp)                   :: fact5(VECTOR_SIZE)
    real(rp)                   :: fact6(VECTOR_SIZE)
    real(rp)                   :: fact7(VECTOR_SIZE)
    real(rp)                   :: fact8(VECTOR_SIZE)
    real(rp)                   :: gpveo(VECTOR_SIZE,3)
    real(rp)                   :: fact1_p(VECTOR_SIZE)
    real(rp)                   :: dtinv_mod(VECTOR_SIZE)
    real(rp),pointer           :: elauu_nsw(:,:,:)
    integer(ip)                :: inode,jnode,jdime
    integer(ip)                :: idofv,jdof2,jdof3
    integer(ip)                :: idof1,idof3,idof2,igaus
    integer(ip)                :: idime,jdof1,jdofv,itime
    integer(ip)                :: ievat,jevat

#ifdef OPENACC
    integer(ip)                :: ivect
#endif

    real(rp), parameter        :: zeror = epsilon(1.0_rp)

    if ( kfl_noslw_ker /= 0_ip ) then   ! falta hacer todo con uprom y mu_eff_prom
       call memory_alloca(mem_modul(1:2,modul),'ELAUU_NSW','nsi_element_assembly_split_oss_default',elauu_nsw,int(VECTOR_SIZE,ip),pnode*ndime,pnode*ndime)
    end if
    
    !----------------------------------------------------------------------
    !
    ! Fractional step: time derivatives is computed using mass matrix
    ! later on, when updating the velocity using AB or RK schemes
    !
    !----------------------------------------------------------------------

    if( NSI_FRACTIONAL_STEP_int == 1_ip) then
       dtinv_mod(:) = 0.0_rp
    else
       dtinv_mod(:) = dtinv_loc(:)
    end if

    if(kfl_surte_nsi == 2_ip) then
       do igaus = 1,pgaus
          gpvis(DEF_VECT, igaus) = gpvis_nsi(DEF_VECT, igaus)/gpden_nsi(DEF_VECT, igaus)
          gpden(DEF_VECT, igaus) = 1.0_rp
       end do
    else
       gpvis(DEF_VECT, 1:pgaus) = gpvis_nsi(DEF_VECT, 1:pgaus)
       gpden(DEF_VECT, 1:pgaus) = gpden_nsi(DEF_VECT, 1:pgaus)
    end if

    !----------------------------------------------------------------------
    !
    ! possibility of using only pressure stabilization - not ready with limiter - nor with shock capturing
    !
    !----------------------------------------------------------------------

    gpsp1_p = gpsp1
    gpsp1_v = gpsp1
    gpsp2_v = gpsp2
    !
    ! Momentum is not stabilized
    !
    if( kfl_nota1_nsi == 1 ) gpsp1_v = 0.0_rp
    !
    ! Algebraic split OSS: stabilization is introduced at the algebraic level
    ! on pressure equation only
    !
    if( kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then
       gpsp1_p = 0.0_rp
       gpsp1_v = 0.0_rp
       gpsp2_v = 0.0_rp
    end if

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    elrbp = 0.0_rp
    elrbu = 0.0_rp
    elapp = 0.0_rp
    elauu = 0.0_rp
    elaup = 0.0_rp
    elapu = 0.0_rp

    !----------------------------------------------------------------------
    !
    ! Test functions
    !
    !----------------------------------------------------------------------

    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !

    if (kfl_regim_nsi==3 .and. kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then
       agrau(DEF_VECT,:,:)   = 0.0_rp
       wgrgr(DEF_VECT,:,:,:) = 0.0_rp
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                   &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
             agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus) / sqrt(gpden(DEF_VECT,igaus))
             do jnode = 1,pnode
                do idime = 1,ndime
                   wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                      &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
       ! grad of the modified velocity

    else 
       agrau(DEF_VECT,:,:)   = 0.0_rp
       wgrgr(DEF_VECT,:,:,:) = 0.0_rp
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                   &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
             agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus) 
             do jnode = 1,pnode
                do idime = 1,ndime
                   wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                      &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Auu
    !
    !----------------------------------------------------------------------
    !
    ! Galerkin + ( tau2 * div(u) , div(v) ) + ( tau1 * rho*a.grad(u), rho*a.grad(v) )
    !
    if (kfl_regim_nsi==3 .and. kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then

       do igaus = 1,pgaus

          fact0(DEF_VECT) = gpsp2_v(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact6(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
          fact7(DEF_VECT) = gpsp1_v(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact8(DEF_VECT) = pabdf_nsi(1) * gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT) + gppor(DEF_VECT,igaus)

          do inode = 1,pnode
             do idime = 1,ndime

                idofv           = (inode-1)*ndime+idime
                fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)

                do jnode = 1,pnode
                   !
                   ! div(u) * tau2' * div(v)
                   !
                   do jdime = 1,ndime
                      jdofv                       = (jnode-1)*ndime+jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,jnode,igaus)
                   end do
                   !
                   ! ( rho/dt N_j + s Nj + rho*(a.grad)Nj ) Ni
                   ! + mu * grad(Ni) . grad(Nj)
                   ! + t1 * rho*(a.grad)Nj * rho*(a.grad)Ni
                   !
                   jdofv           = (jnode-1)*ndime+idime
                   fact4(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                   fact5(DEF_VECT) = fact4(DEF_VECT) * (  fact8(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) ) &
                      &           + fact6(DEF_VECT) *   wgrgr(DEF_VECT,inode,jnode,igaus) &
                      &           + fact7(DEF_VECT) *   agrau(DEF_VECT,jnode,igaus) * agrau(DEF_VECT,inode,igaus)
                   fact5(DEF_VECT) = fact5(DEF_VECT) + fact4(DEF_VECT) * gpden(DEF_VECT,igaus) * agrau(DEF_VECT,jnode,igaus) / sqrt(gpden(DEF_VECT,igaus))  
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact5(DEF_VECT)

                end do
             end do
          end do
       end do
    else
       do igaus = 1,pgaus

          fact0(DEF_VECT) = gpsp2_v(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact6(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
          fact7(DEF_VECT) = gpsp1_v(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact8(DEF_VECT) = pabdf_nsi(1) * gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT) + gppor(DEF_VECT,igaus)

          do inode = 1,pnode
             do idime = 1,ndime

                idofv           = (inode-1)*ndime+idime
                fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)

                do jnode = 1,pnode
                   !
                   ! div(u) * tau2' * div(v)
                   !
                   do jdime = 1,ndime
                      jdofv                       = (jnode-1)*ndime+jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,jnode,igaus)
                   end do
                   !
                   ! ( rho/dt N_j + s Nj + rho*(a.grad)Nj ) Ni
                   ! + mu * grad(Ni) . grad(Nj)
                   ! + t1 * rho*(a.grad)Nj * rho*(a.grad)Ni
                   !
                   jdofv           = (jnode-1)*ndime+idime
                   fact4(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                   fact5(DEF_VECT) = fact4(DEF_VECT) * ( agrau(DEF_VECT,jnode,igaus) + fact8(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) ) &
                      &           + fact6(DEF_VECT) *   wgrgr(DEF_VECT,inode,jnode,igaus) &
                      &           + fact7(DEF_VECT) *   agrau(DEF_VECT,jnode,igaus) * agrau(DEF_VECT,inode,igaus)
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact5(DEF_VECT)

                end do
             end do
            ! elporrbu(DEF_VECT,idime,inode) = elporrbu(DEF_VECT,idime,inode) - fact4(DEF_VECT) * gpvel(DEF_VECT,idime,igaus,1) * gppor(DEF_VECT,igaus)
          end do  ! idime
       end do   ! inode
    end if

    if( kfl_convection_type_nsi == NSI_CONVECTION_CONSERVATIVE ) then
       !
       ! Conservative: non-conservative + rho * (div u)u
       !
       call runend('CONSERVATIVE FORM NOT CODED')
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact2(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             fact2(DEF_VECT) = fact2(DEF_VECT) + gpgve(DEF_VECT,idime,idime,igaus)
          end do
          fact2(DEF_VECT) = fact2(DEF_VECT) * fact0(DEF_VECT)
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime+idime
                   jdofv = (jnode-1)*ndime+idime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
       end do

    else if( kfl_regim_nsi==3 .and.  kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then
       !
       ! Skew-symmetric convection: non-conservative + 1/2 (div rho * u)u
       ! gpgve already contains div(rho u)
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = 0.5_rp * gpvol(DEF_VECT,igaus)
          fact2(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             fact2(DEF_VECT) = fact2(DEF_VECT) + gpgve(DEF_VECT,idime,idime,igaus)
          end do
          fact2(DEF_VECT) = fact2(DEF_VECT) * fact0(DEF_VECT)
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime+idime
                   jdofv = (jnode-1)*ndime+idime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
       end do

       !do igaus = 1,pgaus
       !   fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       !   do inode = 1,pnode
       !      do jnode = 1,pnode
       !         fact1(DEF_VECT) = 0.0_rp
       !         fact2(DEF_VECT) = 0.0_rp
       !         do idime = 1,ndime
       !            fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
       !           fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
       !        end do
       !        fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
       !        fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)
       !        do idime = 1,ndime
       !           idof1 = (inode-1)*ndime+idime
       !           jdof1 = (jnode-1)*ndime+idime
       !           elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact1(DEF_VECT) - fact2(DEF_VECT)
       !        end do
       !     end do
       !  end do
       !end do
       !
       ! Low-Mach with skew symmetric (adding temporal term: d(rho*u)/dt)
       !
       if (kfl_regim_nsi==3 .and. kfl_stabi_nsi /= NSI_GALERKIN ) then
          do itime =2, nbdfp_nsi         ! only rhs with Galerkin needs to be modified
             do igaus =1, pgaus
                gpveo(DEF_VECT,1:ndime) = 0.0_rp  ! velocity at itime
                do inode =1, pnode
                   do idime=1,ndime
                      gpveo(DEF_VECT,idime) = gpveo(DEF_VECT,idime) + elvel(DEF_VECT,idime,inode,itime) * gpsha(DEF_VECT,inode,igaus)
                   end do
                end do
                fact0(DEF_VECT) = 0.5_rp * (gpden(DEF_VECT,igaus) - densi(DEF_VECT,igaus, itime))*pabdf_nsi(itime) &
                     * dtinv_loc(DEF_VECT) * gpvol(DEF_VECT,igaus)
                do inode =1, pnode
                   do idime=1,ndime
                      elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime, inode) &
                           + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpveo(DEF_VECT,idime)
                   end do
                end do
             end do
          end do
       end if


    else if( kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then
       !
       ! Skew-symmetric convection: non-conservative + 1/2 rho * (div u)u
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = 0.5_rp * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact2(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             fact2(DEF_VECT) = fact2(DEF_VECT) + gpgve(DEF_VECT,idime,idime,igaus)
          end do
          fact2(DEF_VECT) = fact2(DEF_VECT) * fact0(DEF_VECT)
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime+idime
                   jdofv = (jnode-1)*ndime+idime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
       end do

    else if( kfl_regim_nsi/=3 .and. kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then
       !
       ! EMAC convection: 2*eps(u) + (div u) u = non-conservative + rho * u.grad(u)^t + rho * (div u) u - 0.5*grad(u·u) [this term
       ! to recover p]
       !
       ! convective_i = rho * uj * ( dui/dxj + duj/dxi ) + rho * (div u) ui
       !              = rho * uj * dui/dxj + rho * uj * duj/dxi + rho * (div u) ui
       !                      +-                            -+         +-                            -+         +-         -+
       !                      | u du/dx + v du/dy + w du/dz  |         | u du/dx + v dv/dx + w dw/dx  |         | (div u) u |
       !              = rho * | u dv/dx + v dv/dy + w dw/dz  | + rho * | u du/dy + v dv/dy + w dw/dy  | + rho * | (div u) v |
       !                      | u dw/dx + v dw/dy + w dw/dz  |         | u du/dz + v dv/dz + w dw/dz  |         | (div u) w |
       !                      +-                            -+         +-                            -+         +-         -+
       !
       do igaus = 1,pgaus

          fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact2(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             fact2(DEF_VECT) = fact2(DEF_VECT) + gpgve(DEF_VECT,idime,idime,igaus)
          end do
          fact2(DEF_VECT) = fact2(DEF_VECT) * fact0(DEF_VECT)
          !
          ! rho * (div u) u
          !
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime+idime
                   jdofv = (jnode-1)*ndime+idime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
          !
          ! rho * u.grad(u)^t
          !
          do inode = 1,pnode
             do idime = 1,ndime
                idofv = (inode-1)*ndime+idime
                do jnode = 1,pnode
                   do jdime = 1,ndime
                      jdofv = (jnode-1)*ndime+jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                           + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus,1)
                   end do
                end do
             end do
          end do

           !
           ! 0.5*grad(u**2)
           !
           do inode = 1,pnode
              do idime = 1,ndime
                 do jdime = 1,ndime
                    elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                       +  fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus,1)*gpgve(DEF_VECT,idime,jdime,igaus)
                 end do
              end do
           end do

       end do

    else if( kfl_regim_nsi==3 .and. kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then
       !
       ! EMAC convection: 2*eps(u) + (div u) u = non-conservative + rho * u.grad(u)^t + rho * (div u) u - 0.5*grad(u·u) [this term
       ! to recover p]
       !
       ! convective_i = rho * uj * ( dui/dxj + duj/dxi ) + rho * (div u) ui
       !              = rho * uj * dui/dxj + rho * uj * duj/dxi + rho * (div u) ui
       !                      +-                            -+         +-                            -+         +-         -+
       !                      | u du/dx + v du/dy + w du/dz  |         | u du/dx + v dv/dx + w dw/dx  |         | (div u) u |
       !              = rho * | u dv/dx + v dv/dy + w dw/dz  | + rho * | u du/dy + v dv/dy + w dw/dy  | + rho * | (div u) v |
       !                      | u dw/dx + v dw/dy + w dw/dz  |         | u du/dz + v dv/dz + w dw/dz  |         | (div u) w |
       !                      +-                            -+         +-                            -+         +-         -+
       !
       do igaus = 1,pgaus

          fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) / sqrt(gpden(DEF_VECT,igaus))
          fact2(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             fact2(DEF_VECT) = fact2(DEF_VECT) + gpgve(DEF_VECT,idime,idime,igaus)
          end do
          fact2(DEF_VECT) = fact2(DEF_VECT) * fact0(DEF_VECT)  
          !
          ! rho * (div u) u
          !
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime+idime
                   jdofv = (jnode-1)*ndime+idime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                      + fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
          !
          ! rho * u.grad(u)^t
          !
          do inode = 1,pnode
             do idime = 1,ndime
                idofv = (inode-1)*ndime+idime
                do jnode = 1,pnode
                   do jdime = 1,ndime
                      jdofv = (jnode-1)*ndime+jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                         + fact0(DEF_VECT) *  (gpsha(DEF_VECT,inode,igaus) * gpden(DEF_VECT,igaus) / sqrt(gpden(DEF_VECT,igaus)))  &
                         & * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus,1)
                   end do
                end do
             end do
          end do

          !
          ! 0.5*grad(u**2)
          !
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                      +  fact0(DEF_VECT) * ( gpsha(DEF_VECT,inode,igaus) * gpden(DEF_VECT,igaus) / sqrt(gpden(DEF_VECT,igaus))) & 
                      &* gpvel(DEF_VECT,jdime,igaus,1)*gpgve(DEF_VECT,idime,jdime,igaus)
                end do
             end do
          end do
       end do

    end if
    !
    ! ( mu*duj/dxi , dv/dxj ) (only div form)
    !
    if( fvins_nsi > 0.9_rp ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idofv = (inode-1)*ndime + idime
                do jnode = 1,pnode
                   fact1(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                   do jdime = 1,ndime
                      jdofv                       = (jnode-1)*ndime + jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,inode,igaus)
                   end do
                end do
                if( fvins_nsi == 2.0_rp ) then
                   fact1(DEF_VECT) = -2.0_rp / 3.0_rp * gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                   do jnode = 1,pnode
                      do jdime = 1,ndime
                         jdofv                       = (jnode-1)*ndime + jdime
                         elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,jnode,igaus)
                      end do
                   end do
                end if
             end do
          end do
       end do
    end if
    if ( kfl_noslw_ker /= 0_ip ) then   ! falta hacer todo con uprom y mu_eff_prom       ! OJO lo de vector size hay que ver como manejarlo supongo por ahora para probar corre con vector size=1
       ! SUPONGO QUE NO JOIDERÁ SI LA VISCO QUEDA A 0 PENSAR BIEN
       !       if ( kfl_nswel_ker(ielem) /= 0_ip ) then   ! this I will avoid it because it complicates for vector case - it will do it for all elemenst but in most gpvis_nsw = 0
       !
       ! mu * grad(Ni) . grad(Nj)
       !
       do igaus = 1,pgaus

          fact6(DEF_VECT) = gpvis_nsw(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)                        

          do inode = 1,pnode
             do idime = 1,ndime
                idofv           = (inode-1)*ndime+idime
                do jnode = 1,pnode
                   !
                   !  mu * grad(Ni) . grad(Nj)
                   !
                   jdofv           = (jnode-1)*ndime+idime
                   fact5(DEF_VECT) =  fact6(DEF_VECT) *   wgrgr(DEF_VECT,inode,jnode,igaus)
                   elauu_nsw(DEF_VECT,idofv,jdofv) = elauu_nsw(DEF_VECT,idofv,jdofv) + fact5(DEF_VECT)
                end do
             end do
          end do
       end do
       !
       ! ( mu*duj/dxi , dv/dxj ) (only div form)
       !
       if( fvins_nsi > 0.9_rp ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime + idime
                   do jnode = 1,pnode
                      fact1(DEF_VECT) = gpvis_nsw(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                      do jdime = 1,ndime
                         jdofv                       = (jnode-1)*ndime + jdime
                         elauu_nsw(DEF_VECT,idofv,jdofv) = elauu_nsw(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,inode,igaus)
                      end do
                   end do
                   if( fvins_nsi == 2.0_rp ) then
                      fact1(DEF_VECT) = -2.0_rp / 3.0_rp * gpvis_nsw(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                      do jnode = 1,pnode
                         do jdime = 1,ndime
                            jdofv                       = (jnode-1)*ndime + jdime
                            elauu_nsw(DEF_VECT,idofv,jdofv) = elauu_nsw(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,jnode,igaus)
                         end do
                      end do
                   end if
                end do
             end do
          end do
       end if
       !
       !send to RHS - I do it here to avoid getting elauu_nsw outside of this subroutine
       !
       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndime+idime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        - elauu_nsw(DEF_VECT,ievat,jevat) * elavv(DEF_VECT,jdime,jnode)
                end do
             end do
          end do
       end do
    end if
 
    !call nsi_element_system_output(&
    !     pnode,elauu(1,:,:),elaup(1,:,:),elapp(1,:,:),elapu(1,:,:),elrbu(1,:,:),elrbp(1,:),&
    !     elauq(1,:,:),elapq(1,:,:),elaqu(1,:,:),elaqp(1,:,:),elaqq(1,:,:),elrbq(1,:))
    !stop
    !
    ! Lumped evolution matrix (only backward euler)
    !
    if( kfl_lumped == 1 ) then
       !
       ! Remove Galerkin term and add lumped term
       !
       if( ndime == 2 ) then
          call runend('PREGUNTAR A MATIAS QUE LO PROGRAME')
       else
          do igaus = 1,pgaus
             gpveo(DEF_VECT,1:3) = 0.0_rp
             do inode = 1,pnode
                do idime = 1,ndime
                   gpveo(DEF_VECT,idime) = gpveo(DEF_VECT,idime) + elvel(DEF_VECT,idime,inode,2) * gpsha(DEF_VECT,inode,igaus)
                end do
             end do
             do inode = 1,pnode
                idof1                       = 3*inode-2
                idof2                       = 3*inode-1
                idof3                       = 3*inode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_mod(DEF_VECT)
                elauu(DEF_VECT,idof1,idof1) = elauu(DEF_VECT,idof1,idof1) + fact0(DEF_VECT)
                elauu(DEF_VECT,idof2,idof2) = elauu(DEF_VECT,idof2,idof2) + fact0(DEF_VECT)
                elauu(DEF_VECT,idof3,idof3) = elauu(DEF_VECT,idof3,idof3) + fact0(DEF_VECT)
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - fact0(DEF_VECT) * gpveo(DEF_VECT,idime)
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact0(DEF_VECT) * elvel(DEF_VECT,idime,inode,2)
                end do
                do jnode = 1,pnode
                   jdof1                       = 3*jnode-2
                   jdof2                       = 3*jnode-1
                   jdof3                       = 3*jnode
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
       end if

    else if( kfl_lumped == 2 ) then
       !
       ! No time term have been added up to now: add Galerkin term
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT)
          do inode = 1, pnode
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             do idime = 1,ndime
                idof1                       = (inode-1) * ndime + idime
                elauu(DEF_VECT,idof1,idof1) = elauu(DEF_VECT,idof1,idof1) + fact1(DEF_VECT)
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT) * elvel(DEF_VECT,idime,inode,2)
             end do
          end do
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Apu and Aup
    !
    !----------------------------------------------------------------------
    !
    ! ( div(u) , q ) and - ( p , div(v) )
    !
    if( ndime == 2 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             idof1 = 2*inode-1
             idof2 = 2*inode
             do jnode = 1,pnode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus)       * gpsha(DEF_VECT,jnode,igaus)
                fact1(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,1,inode,igaus)
                fact2(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,2,inode,igaus)
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + fact1(DEF_VECT)
                elapu(DEF_VECT,jnode,idof2) = elapu(DEF_VECT,jnode,idof2) + fact2(DEF_VECT)
                elaup(DEF_VECT,idof1,jnode) = elaup(DEF_VECT,idof1,jnode) - fact1(DEF_VECT)
                elaup(DEF_VECT,idof2,jnode) = elaup(DEF_VECT,idof2,jnode) - fact2(DEF_VECT)
             end do
          end do
       end do
    else
       do igaus = 1,pgaus
          do inode = 1,pnode
             idof1 = 3*inode-2
             idof2 = 3*inode-1
             idof3 = 3*inode
             do jnode = 1,pnode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus)       * gpsha(DEF_VECT,jnode,igaus)
                fact1(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,1,inode,igaus)
                fact2(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,2,inode,igaus)
                fact3(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,3,inode,igaus)
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + fact1(DEF_VECT)
                elapu(DEF_VECT,jnode,idof2) = elapu(DEF_VECT,jnode,idof2) + fact2(DEF_VECT)
                elapu(DEF_VECT,jnode,idof3) = elapu(DEF_VECT,jnode,idof3) + fact3(DEF_VECT)
                elaup(DEF_VECT,idof1,jnode) = elaup(DEF_VECT,idof1,jnode) - fact1(DEF_VECT)
                elaup(DEF_VECT,idof2,jnode) = elaup(DEF_VECT,idof2,jnode) - fact2(DEF_VECT)
                elaup(DEF_VECT,idof3,jnode) = elaup(DEF_VECT,idof3,jnode) - fact3(DEF_VECT)
             end do
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! App
    !
    !----------------------------------------------------------------------
    !
    ! Pressure: ( tau1' * grad(p) , grad(q) )
    !
    ! In the case of Fractional step and algebaric split methods, the pressure
    ! stabilization is carried out further on
    !
    if( kfl_stabi_nsi /= NSI_GALERKIN .and. kfl_stabi_nsi /= NSI_ALGEBRAIC_SPLIT_OSS ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do jnode = inode+1,pnode
                fact1(DEF_VECT)             = gpsp1_p(DEF_VECT,igaus) * wgrgr(DEF_VECT,jnode,inode,igaus) * gpvol(DEF_VECT,igaus)
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + fact1(DEF_VECT)
                elapp(DEF_VECT,inode,jnode) = elapp(DEF_VECT,inode,jnode) + fact1(DEF_VECT)
             end do
             fact1(DEF_VECT)             = gpsp1_p(DEF_VECT,igaus) * wgrgr(DEF_VECT,inode,inode,igaus) * gpvol(DEF_VECT,igaus)
             elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT)
          end do
       end do
    end if
    !
    ! Penalization
    !
    do igaus = 1,pgaus
       fact1(DEF_VECT) = penal_nsi * gpvol(DEF_VECT,igaus)
       do inode = 1,pnode
          elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
          elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elpre(DEF_VECT,inode,1)
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! bu and bp
    !
    ! P1  = P [ tau1' * rho * a.grad(u) ]
    ! P1' = P1 + tau1' * rho * u'n / dt
    !
    ! P2  = P [ tau1' * ( grad(p) - f ) ]
    ! P2' = P2 + tau1' * rho * u'n / dt + tau1' * f
    !
    !----------------------------------------------------------------------
    !
    ! Limiter
    !
    if( kfl_limit_nsi == -1 .or. kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then

       gpvep(DEF_VECT,:,:) = 0.0_rp

    else if( kfl_limit_nsi > 0 ) then

       do igaus = 1,pgaus
          c1(DEF_VECT) = 0.0_rp
          c2(DEF_VECT) = 0.0_rp
          c3(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             c4(DEF_VECT) = 0.0_rp
             do inode = 1,pnode
                c4(DEF_VECT) = c4(DEF_VECT) + agrau(DEF_VECT,inode,igaus) * elvel(DEF_VECT,idime,inode,1)
             end do
             c4(DEF_VECT) = gpsp1(DEF_VECT,igaus) * c4(DEF_VECT)
             c1(DEF_VECT) = c1(DEF_VECT) + ( gpvep(DEF_VECT,idime,igaus) - c4(DEF_VECT) )**2
             c3(DEF_VECT) = c3(DEF_VECT) + gpvep(DEF_VECT,idime,igaus) * gpvep(DEF_VECT,idime,igaus)
             c2(DEF_VECT) = c2(DEF_VECT) + c4(DEF_VECT) * c4(DEF_VECT)
          end do
          c3(DEF_VECT)   = sqrt( c2(DEF_VECT) ) + sqrt( c3(DEF_VECT) )
          c1(DEF_VECT)   = sqrt( c1(DEF_VECT) )
          beta(DEF_VECT) = c1(DEF_VECT) / ( c3(DEF_VECT) + epsilon(1.0_rp) )
          if( kfl_limit_nsi == 1 ) then
             alpha(DEF_VECT) = min(1.0_rp,2.0_rp*(1.0_rp-beta(DEF_VECT)))
          else if( kfl_limit_nsi == 2 ) then
             alpha(DEF_VECT) = 0.5_rp*(tanh(20.0_rp*(beta(DEF_VECT)-0.8_rp))+1.0_rp)
          end if
          do idime = 1,ndime
             gpvep(DEF_VECT,idime,igaus) = alpha(DEF_VECT) * gpvep(DEF_VECT,idime,igaus)
          end do
       end do

    end if
    !
    ! P2 <= P2 + tau1' * f
    !
    if( kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then
       gpgrp(DEF_VECT,:,:) = 0.0_rp
    else
       do igaus = 1,pgaus
          do idime = 1,ndime
             gpgrp(DEF_VECT,idime,igaus) = gpgrp(DEF_VECT,idime,igaus) + gpsp1_p(DEF_VECT,igaus) * gprhs(DEF_VECT,idime,igaus)
          end do
       end do
       !
       ! P1 <= P1 + tau1' * rho * u'n / dt
       ! P2 <= P2 + tau1' * rho * u'n / dt
       !
       if( kfl_sgsti_nsi == 1 ) then
          do igaus = 1,pgaus
             fact1(DEF_VECT)    = gpden(DEF_VECT,igaus) * dtsgs(DEF_VECT) * gpsp1_v(DEF_VECT,igaus)
             fact1_p (DEF_VECT) = gpden(DEF_VECT,igaus) * dtsgs(DEF_VECT) * gpsp1_p(DEF_VECT,igaus)
             do idime = 1,ndime
                gpvep(DEF_VECT,idime,igaus) = gpvep(DEF_VECT,idime,igaus) + fact1(DEF_VECT)   * gpsgs(DEF_VECT,idime,igaus,2)
                gpgrp(DEF_VECT,idime,igaus) = gpgrp(DEF_VECT,idime,igaus) + fact1_p(DEF_VECT) * gpsgs(DEF_VECT,idime,igaus,2)
             end do
          end do
       end if
    end if
    !
    ! bu = ( f + rho*u^n/dt , v ) + ( rho * a.grad(v) , tau1' * rho u'^n/dt + P1 )
    !    = ( f + rho*u^n/dt , v ) + ( rho * a.grad(v) , P1' )
    !
    ! bp = ( f + rho*u'^n/dt , tau1' grad(q) ) + ( P2 , grad(q) )
    !    = ( P2' , grad(q) )
    !

    do igaus = 1,pgaus
       fact4(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT)
       do itime = 2,nbdfp_nsi
          do idime = 1,ndime
             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - pabdf_nsi(itime) * fact4(DEF_VECT) * gpvel(DEF_VECT,idime,igaus,itime)
          end do
       end do
       do inode = 1,pnode
          fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f + rho*u^n/dt , v )
          fact3(DEF_VECT) = gpvol(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)  ! ( rho * a.grad(v) , P1' )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT) * gprhs(DEF_VECT,idime,igaus) &
                &                                                    + fact3(DEF_VECT) * gpvep(DEF_VECT,idime,igaus)
          end do
          elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * gprhc(DEF_VECT,igaus)                ! ( rhs, q )
          do idime = 1,ndime
             elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus) * gpgrp(DEF_VECT,idime,igaus) ! ( P2' , grad(q) )
          end do
       end do 
    end do

    !--------------------------------------------------------------------
    !
    ! Pressure bubble
    !
    !--------------------------------------------------------------------

    if( maxval(pbubl) == 1 ) then
       if( kfl_stabi_nsi /= NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) call runend('BUBBLE NOT CODED FOR SPLIT OSS')
       !
       ! Initialization
       !
       elauq = 0.0_rp
       elapq = 0.0_rp
       elaqu = 0.0_rp
       elaqp = 0.0_rp
       elaqq = 0.0_rp
       elrbq = 0.0_rp
       !
       ! Auq and Aqu
       !
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) - fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
                   elaqu(DEF_VECT,1,idofv) = elaqu(DEF_VECT,1,idofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                   elaqu(DEF_VECT,1,idofv) = elaqu(DEF_VECT,1,idofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
                end do
             end do
          end do
       end if
       !
       ! Penalization and others
       !
       do igaus = 1,pgaus
          elaqq(DEF_VECT,1,1) = elaqq(DEF_VECT,1,1) + gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * penal_nsi
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * penal_nsi * elbub(DEF_VECT)
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * gprhc(DEF_VECT,igaus)
       end do

    end if


    if ( kfl_noslw_ker /= 0_ip ) then   ! falta hacer todo con uprom y mu_eff_prom
       call memory_deallo(mem_modul(1:2,modul),'ELAUU_NSW','nsi_element_assembly_split_oss_default',elauu_nsw)
    end if


  end subroutine

