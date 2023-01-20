!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!





!- parte de mod_nsi_element_operations_hh80  -- esto es solo paar P1

module mod_nsi_element_operations_hh80

#include "def_vector_size.inc"
   use def_kintyp_basic, only: ip, rp         !in_const
   use def_domain, only: ndime         !in_const_scalar  ! actually paarmeter if ndimepar

#ifdef _OPENACC
   use openacc
#endif

   use def_master, only: kfl_paral

   implicit none

   private
   public :: nsi_element_operations_hh80

contains

   !-----------------------------------------------------------------------
   !>
   !> @author  houzeaux
   !> @date    2020-05-06
   !> @brief   Fast assembly of NS
   !> @details Assembly of Navier-Stokes for CPU and GPU, with
   !>          restricted options
   !>
   !-----------------------------------------------------------------------

#if defined PNODE_VALUE && defined PGAUS_VALUE
   subroutine nsi_element_operations_hh80( &
      VECTOR_DIM, qnode, qgaus, list_elements, time1)
#else
      subroutine nsi_element_operations_hh80( &
         VECTOR_DIM, pnode, pgaus, list_elements, time1)
#endif
!@@@    !$acc routine(vecnor,frivel) seq
         use def_kintyp, only: ip, rp                                     ! in_const_scalar
         use def_master, only: rhsid                                     ! out       ! real(rp), pointer     :: rhsid(:)
         use def_master, only: veloc                                     ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
         use def_domain, only: coord                                     ! in_const  ! real(rp), pointer     :: coord(:,:)
         use def_domain, only: ltype                                     ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
         use def_domain, only: lnods                                     ! in_const  ! integer(ip), pointer  :: lnods(:,:)
         use def_nastin, only: grnor_nsi                                 ! in_const_scalar
         use def_nastin, only: gravi_nsi                                 ! in_const  ! real(rp)         ::gravi_nsi(3)
         use def_nastin, only: densi_aux, visco_aux
         !
         implicit none
         !
         ! Input and output variables
         !
         integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
#if defined PNODE_VALUE && defined PGAUS_VALUE
         integer(ip), intent(in)          :: qnode                        !< Number of nodes
         integer(ip), intent(in)          :: qgaus                        !< Number of Gauss points
         integer(ip), parameter           :: pnode = PNODE_VALUE
         integer(ip), parameter           :: pgaus = PGAUS_VALUE
#else
         integer(ip), intent(in)          :: pnode                        !< Number of nodes
         integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
#endif
         integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

         real(rp), intent(inout)       :: time1(10)                                        ! Timings
         !
         ! Element matrices and vectors (stiffness and preconditioner)
         !
         real(rp)                         :: elrbu(VECTOR_DIM, ndime, pnode)                    ! bu
         !
         ! Gather
         !
         real(rp)                         :: elvel(VECTOR_DIM, ndime, pnode)                    ! u
         real(rp)                         :: elcod(VECTOR_DIM, ndime, pnode)                    ! x
         !
         ! Indices and dimensions
         !
         integer(ip)                      :: ielem, inode, ivect
         integer(ip)                      :: pelty
         integer(ip)                      :: ipoin, igaus
         integer(ip)                      :: lnods_loc(VECTOR_DIM, pnode)
         !
         ! Gauss point values
         !
         real(rp)                         :: gpsha(VECTOR_DIM, pnode, pgaus)                    ! N
         real(rp)                         :: gpcar(VECTOR_DIM, ndime, pnode)                    ! dN/dxi
         real(rp)                         :: gpvol(VECTOR_DIM)                          ! w*|J|, |J|
         real(rp)                         :: gpadv(VECTOR_DIM, ndime, pgaus)                    ! u+u'
         real(rp)                         :: gprhs(VECTOR_DIM, ndime, pgaus)                    ! RHS
         real(rp)                         :: gpvel(VECTOR_DIM, ndime, pgaus)                    ! u
         real(rp)                         :: gpgve(VECTOR_DIM, ndime, ndime)              ! grad(u)
         real(rp)                         :: xjaci(VECTOR_DIM, ndime, ndime)
         real(rp)                         :: gpdet(VECTOR_DIM)

         real(rp)                         :: timea, timeb
         integer(ip)                      :: idime
         integer(ip)                      :: jdime
         integer(ip)                      :: iauxi
         !
         ! For Vreman
         !
         real(rp), parameter        :: hnatu = 1.0   ! would need to correct for elem /=1!!!!
         integer(ip)               :: kdime
         real(rp)                  :: xnutu(VECTOR_DIM), xmile(VECTOR_DIM)
         real(rp)                  :: G__ij(VECTOR_DIM, 3, 3)
         real(rp)                  :: Aalph(VECTOR_DIM), Bbeta(VECTOR_DIM)
         real(rp), parameter        :: const = 0.1_rp   ! seria mejor poner el valor leido
         real(rp), parameter        :: zeror = epsilon(1.0_rp) ! podria usarlo de def_master

! P1

         real(rp), parameter                   :: weigp = 0.041666666666666664D+00
         real(rp), parameter                   :: epsilgeo_div = epsilon(1.0_rp) !epsilgeo usad to avoid divisions by zero
         real(rp), parameter                   :: sha_aux(4, 4) = reshape((/0.585410196624969D+00, 0.138196601125011D+00, &
                                                              0.138196601125011D+00, 0.138196601125011D+00, 0.138196601125010D+00, &
                                                              0.585410196624969D+00, 0.138196601125011D+00, 0.138196601125011D+00, &
                                                              0.138196601125010D+00, 0.138196601125011D+00, 0.585410196624969D+00, &
                                                                            0.138196601125011D+00, 0.138196601125010D+00, &
                                                    0.138196601125011D+00, 0.138196601125011D+00, 0.585410196624969D+00/), (/4, 4/))

         real(rp)                  :: visco_tot(VECTOR_DIM)
         !
         ! Internal
         !
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

#ifdef OPENACCHHH

#define FACT1X     fact1
#define FACT2X     fact2
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT1X     fact1(1:VECTOR_SIZE)
#define FACT2X     fact2(1:VECTOR_SIZE)
#define T1X        t1(1:VECTOR_SIZE)
#define T2X        t2(1:VECTOR_SIZE)
#define T3X        t3(1:VECTOR_SIZE)
#define DENOMX     denom(1:VECTOR_SIZE)

#endif

         real(rp)    :: FACT2X
         real(rp)    :: T1X
         real(rp)    :: T2X
         real(rp)    :: T3X
         real(rp)    :: DENOMX

         !--------------------------------------------------------------------
         !
         ! Gather: global to local
         !
         !--------------------------------------------------------------------

         call cputim(timea)

         ielem = list_elements(1)
         pelty = abs(ltype(ielem))


         call cputim(timeb)
         time1(3) = time1(3) + timeb - timea
         call cputim(timea)

         !--------------------------------------------------------------------
         !
         ! Shape function and derivatives   !ver si gpvis_nsw tine que entrar aca
         ! added avelavv  becaise it was giving me error at run time not sure if it is correct
         ! tambien avupo_ker,elavv
         ! Aunque estoy corriendo un caso sin no slip walllaw y esas solo aparcen dentro de if las quiere igual
         !
         !--------------------------------------------------------------------

         !$acc enter data create(  lnods_loc , gpvol     ,                     &
         !$acc              elvel     ,  gpadv , gpvel ,                       &
         !$acc              gpgve   , elcod     ,  gpcar ,                     &
         !$acc              gpdet   ,                                          &
         !$acc              xjaci   ,                                          &
         !$acc              xmile,xnutu, G__ij    , Bbeta,              &
         !$acc              gprhs   , elrbu, visco_tot, aalph,gpsha                         )          &

         !$acc copyin(                                                         &
         !$acc              list_elements,gravi_nsi                            )
         !
         ! This do starts here both in the openacc version and in the nonopenacc
         ! for the non opencc it ends 30  lines later,
         ! in the openacc case it covers all the subroutine.
         ! Similarly the scatter in the nonopenacc case needs a do ivect
         !
         !$acc parallel loop gang vector default(present)
         !
         do ivect = 1, VECTOR_DIM
            ielem = abs(list_elements(ivect))
            if (ielem /= 0) then
               lnods_loc(ivect, 1:pnode) = lnods(1:pnode, ielem)
               ielem = list_elements(ivect)
            else
               lnods_loc(ivect, 1:pnode) = lnods(1:pnode, list_elements(1))
               ielem = list_elements(1)
            end if
            !
            ! Transient
            !
            do inode = 1, pnode
               ipoin = lnods_loc(ivect, inode)
               elvel(ivect, 1:ndime, inode) = veloc(1:ndime, ipoin, 1)
               elcod(ivect, 1:ndime, inode) = coord(1:ndime, ipoin)
            end do
#ifndef OPENACCHHH
         end do

         call cputim(timeb)
         time1(1) = time1(1) + timeb - timea
         call cputim(timea)
#endif

         !--------------------------------------------------------------------
         !
         ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
         !
         !--------------------------------------------------------------------

         do igaus = 1, pgaus
            do inode = 1, pnode
               gpsha(DEF_VECT, inode, igaus) = sha_aux(inode, igaus)
            end do
         end do
         !
         ! GPCAR (from elmgeo_cartesian_derivatives_jacobian_vector), and GPVOL
         !
         !
         ! 3D P1 element
         !
         gpcar(DEF_VECT, 1, 1) = elcod(DEF_VECT, 1, 2) - elcod(DEF_VECT, 1, 1)
         gpcar(DEF_VECT, 1, 2) = elcod(DEF_VECT, 1, 3) - elcod(DEF_VECT, 1, 1)
         gpcar(DEF_VECT, 1, 3) = elcod(DEF_VECT, 1, 4) - elcod(DEF_VECT, 1, 1)
         gpcar(DEF_VECT, 2, 1) = elcod(DEF_VECT, 2, 2) - elcod(DEF_VECT, 2, 1)
         gpcar(DEF_VECT, 2, 2) = elcod(DEF_VECT, 2, 3) - elcod(DEF_VECT, 2, 1)
         gpcar(DEF_VECT, 2, 3) = elcod(DEF_VECT, 2, 4) - elcod(DEF_VECT, 2, 1)
         gpcar(DEF_VECT, 3, 1) = elcod(DEF_VECT, 3, 2) - elcod(DEF_VECT, 3, 1)
         gpcar(DEF_VECT, 3, 2) = elcod(DEF_VECT, 3, 3) - elcod(DEF_VECT, 3, 1)
         gpcar(DEF_VECT, 3, 3) = elcod(DEF_VECT, 3, 4) - elcod(DEF_VECT, 3, 1)
         T1X = gpcar(DEF_VECT, 2, 2)*gpcar(DEF_VECT, 3, 3) - gpcar(DEF_VECT, 3, 2)*gpcar(DEF_VECT, 2, 3)
         T2X = -gpcar(DEF_VECT, 2, 1)*gpcar(DEF_VECT, 3, 3) + gpcar(DEF_VECT, 3, 1)*gpcar(DEF_VECT, 2, 3)
         T3X = gpcar(DEF_VECT, 2, 1)*gpcar(DEF_VECT, 3, 2) - gpcar(DEF_VECT, 3, 1)*gpcar(DEF_VECT, 2, 2)
         gpdet(DEF_VECT) = gpcar(DEF_VECT, 1, 1)*T1X + gpcar(DEF_VECT, 1, 2)*T2X + gpcar(DEF_VECT, 1, 3)*T3X

         DENOMX = 1.0_rp/(sign(1.0_rp, gpdet(DEF_VECT))*max(abs(gpdet(DEF_VECT)), epsilgeo_div))

         xjaci(DEF_VECT, 1, 1) = T1X*DENOMX
         xjaci(DEF_VECT, 2, 1) = T2X*DENOMX
         xjaci(DEF_VECT, 3, 1) = T3X*DENOMX
         xjaci(DEF_VECT, 2, 2) = (gpcar(DEF_VECT, 1, 1)*gpcar(DEF_VECT, 3, 3) - gpcar(DEF_VECT, 3, 1)*gpcar(DEF_VECT, 1, 3))*DENOMX
         xjaci(DEF_VECT, 3, 2) = (-gpcar(DEF_VECT, 1, 1)*gpcar(DEF_VECT, 3, 2) + gpcar(DEF_VECT, 1, 2)*gpcar(DEF_VECT, 3, 1))*DENOMX
         xjaci(DEF_VECT, 3, 3) = (gpcar(DEF_VECT, 1, 1)*gpcar(DEF_VECT, 2, 2) - gpcar(DEF_VECT, 2, 1)*gpcar(DEF_VECT, 1, 2))*DENOMX
         xjaci(DEF_VECT, 1, 2) = (-gpcar(DEF_VECT, 1, 2)*gpcar(DEF_VECT, 3, 3) + gpcar(DEF_VECT, 3, 2)*gpcar(DEF_VECT, 1, 3))*DENOMX
         xjaci(DEF_VECT, 1, 3) = (gpcar(DEF_VECT, 1, 2)*gpcar(DEF_VECT, 2, 3) - gpcar(DEF_VECT, 2, 2)*gpcar(DEF_VECT, 1, 3))*DENOMX
         xjaci(DEF_VECT, 2, 3) = (-gpcar(DEF_VECT, 1, 1)*gpcar(DEF_VECT, 2, 3) + gpcar(DEF_VECT, 2, 1)*gpcar(DEF_VECT, 1, 3))*DENOMX

         gpcar(DEF_VECT, 1, 1) = -xjaci(DEF_VECT, 1, 1) - xjaci(DEF_VECT, 2, 1) - xjaci(DEF_VECT, 3, 1)
         gpcar(DEF_VECT, 1, 2) = xjaci(DEF_VECT, 1, 1)
         gpcar(DEF_VECT, 1, 3) = xjaci(DEF_VECT, 2, 1)
         gpcar(DEF_VECT, 1, 4) = xjaci(DEF_VECT, 3, 1)
         gpcar(DEF_VECT, 2, 1) = -xjaci(DEF_VECT, 1, 2) - xjaci(DEF_VECT, 2, 2) - xjaci(DEF_VECT, 3, 2)
         gpcar(DEF_VECT, 2, 2) = xjaci(DEF_VECT, 1, 2)
         gpcar(DEF_VECT, 2, 3) = xjaci(DEF_VECT, 2, 2)
         gpcar(DEF_VECT, 2, 4) = xjaci(DEF_VECT, 3, 2)
         gpcar(DEF_VECT, 3, 1) = -xjaci(DEF_VECT, 1, 3) - xjaci(DEF_VECT, 2, 3) - xjaci(DEF_VECT, 3, 3)
         gpcar(DEF_VECT, 3, 2) = xjaci(DEF_VECT, 1, 3)
         gpcar(DEF_VECT, 3, 3) = xjaci(DEF_VECT, 2, 3)
         gpcar(DEF_VECT, 3, 4) = xjaci(DEF_VECT, 3, 3)

         gpvol(DEF_VECT) = weigp*gpdet(DEF_VECT)
!!!!!!!!!!!!!!!!!!!
         xmile(DEF_VECT) = ((hnatu/sqrt(xjaci(DEF_VECT, 1, 1)*xjaci(DEF_VECT, 1, 1) &
               &     + xjaci(DEF_VECT, 1, 2)*xjaci(DEF_VECT, 1, 2)           &
               &     + xjaci(DEF_VECT, 1, 3)*xjaci(DEF_VECT, 1, 3)))*      &
               &     (hnatu/sqrt(xjaci(DEF_VECT, 2, 1)*xjaci(DEF_VECT, 2, 1) &
               &     + xjaci(DEF_VECT, 2, 2)*xjaci(DEF_VECT, 2, 2)           &
               &     + xjaci(DEF_VECT, 2, 3)*xjaci(DEF_VECT, 2, 3)))*       &
               &     (hnatu/sqrt(xjaci(DEF_VECT, 3, 1)*xjaci(DEF_VECT, 3, 1) &
               &     + xjaci(DEF_VECT, 3, 2)*xjaci(DEF_VECT, 3, 2)           &
               &     + xjaci(DEF_VECT, 3, 3)*xjaci(DEF_VECT, 3, 3))))**0.3333333_rp

!!!!!!!!!!!!!!!!!!!

         !--------------------------------------------------------------------
         !
         ! Properties
         !
         !--------------------------------------------------------------------

         !----------------------------------------------------------------------
         !
         ! Gauss point values
         !
         !----------------------------------------------------------------------

         gpgve(DEF_VECT, :, :) = 0.0_rp
         gprhs(DEF_VECT, :, :) = 0.0_rp
         gpvel(DEF_VECT, :, :) = 0.0_rp

         do idime = 1, ndime
            do inode = 1, pnode
               do jdime = 1, ndime
                  gpgve(DEF_VECT, jdime, idime) = gpgve(DEF_VECT, jdime, idime) &
                                                  + gpcar(DEF_VECT, jdime, inode)*elvel(DEF_VECT, idime, inode)
               end do
            end do
         end do
         do igaus = 1, pgaus
            FACT2X = densi_aux*grnor_nsi
            do idime = 1, ndime
               do inode = 1, pnode
         gpvel(DEF_VECT, idime, igaus) = gpvel(DEF_VECT, idime, igaus) + elvel(DEF_VECT, idime, inode)*gpsha(DEF_VECT, inode, igaus)
               end do
               gpadv(DEF_VECT, idime, igaus) = gpvel(DEF_VECT, idime, igaus)
               gprhs(DEF_VECT, idime, igaus) = gprhs(DEF_VECT, idime, igaus) + FACT2X*gravi_nsi(idime)
            end do
         end do

!!!!!!!!!!!!!!!!!!!   VREMAN  - I need to add igaus & DEF_VEC   - beware I also need to obtain xmile a priori
! en ker_turbul.f90   pense que se usba DEF_VECT pero no supongo que caca va a ser todo mucho mas rapido
         G__ij(DEF_VECT, :, :) = 0.0_rp ! G = g^T*g
         do idime = 1_ip, 3_ip
            do jdime = 1_ip, 3_ip
               do kdime = 1_ip, 3_ip
                  G__ij(DEF_VECT, idime, jdime) = G__ij(DEF_VECT, idime, jdime) + (gpgve(DEF_VECT, idime, kdime) &
                                                                     *gpgve(DEF_VECT, jdime, kdime)*xmile(DEF_VECT)*xmile(DEF_VECT))
               end do
            end do
         end do

!    I guess I do not need to add DEF_VECT  and igauss everywhere it shoul be don automaticaly if Aalph,Bbeta & xmile are of size(DEF_VECT,pgaus)

         Aalph(DEF_VECT) = (G__ij(DEF_VECT, 1_ip, 1_ip) + G__ij(DEF_VECT, 2_ip, 2_ip) + G__ij(DEF_VECT, 3_ip, 3_ip))/(xmile(DEF_VECT)*xmile(DEF_VECT))
    Bbeta(DEF_VECT) =  G__ij(DEF_VECT,1_ip,1_ip)*G__ij(DEF_VECT,2_ip,2_ip) + G__ij(DEF_VECT,2_ip,2_ip)*G__ij(DEF_VECT,3_ip,3_ip) + G__ij(DEF_VECT,3_ip,3_ip)*G__ij(DEF_VECT,1_ip,1_ip) &
             - G__ij(DEF_VECT,1_ip,2_ip)*G__ij(DEF_VECT,1_ip,2_ip) - G__ij(DEF_VECT,2_ip,3_ip)*G__ij(DEF_VECT,2_ip,3_ip) - G__ij(DEF_VECT,1_ip,3_ip)*G__ij(DEF_VECT,1_ip,3_ip)

#ifndef OPENACCHHH


    xnutu = const
    where ( Aalph > 10.0_rp*zeror )                        ! Avoid divide by zero
       xnutu = xnutu  * sqrt (max(( Bbeta ) / ( Aalph ), zeror))
    elsewhere ( Aalph <= 10.0_rp*zeror )
       xnutu = 0.0_rp
    endwhere
#else
         xnutu(DEF_VECT) = const
         if (Aalph(DEF_VECT) > 10.0_rp*zeror) then                        ! Avoid divide by zero
            xnutu(DEF_VECT) = xnutu(DEF_VECT)*sqrt(max((Bbeta(DEF_VECT))/(Aalph(DEF_VECT)), zeror))
         else if (Aalph(DEF_VECT) <= 10.0_rp*zeror) then
            xnutu(DEF_VECT) = 0.0_rp
         end if

#endif

         !    do igaus=1,pgaus
!       gpvis(:,igaus) = gpvis(:,igaus) + densi_aux * xnutu   ! creo  correcto  : gpvis(VECTOR_DIM,pgaus), xnutu(VECTOR_DIM,pgaus)
!    end do
         ! note taht now viscotot is not of size pgaus , just vectorsize
         visco_tot(DEF_VECT) = visco_aux + densi_aux*xnutu(DEF_VECT)

!!!!!!!

         !----------------------------------------------------------------------
         !
         ! Element matrices
         !
         !----------------------------------------------------------------------

         elrbu(DEF_VECT, :, :) = 0.0_rp

         !--------------------------------------------------------------------
         !
         ! Convective,viscous term and RHS  -- todo tidy up
         !
         !--------------------------------------------------------------------

         do inode = 1, pnode
            do idime = 1, ndime
               elrbu(DEF_VECT, idime, inode) = 0.0
               do igaus = 1, pgaus
                elrbu(DEF_VECT, idime, inode) = elrbu(DEF_VECT, idime, inode) + gpvol(DEF_VECT)*gpsha(DEF_VECT,inode, igaus)*gprhs(DEF_VECT, idime, igaus)
                  do jdime = 1, ndime
                     elrbu(DEF_VECT, idime, inode) = elrbu(DEF_VECT, idime, inode) &
- densi_aux*gpvol(DEF_VECT)*gpsha(DEF_VECT, inode, igaus)*gpvel(DEF_VECT, jdime, igaus)*gpgve(DEF_VECT, jdime, idime) &     ! non-conservative - (u.grad) u
                        - visco_tot(DEF_VECT) * gpvol(DEF_VECT) * gpcar(DEF_VECT,jdime,inode) * (gpgve(DEF_VECT,jdime,idime) + gpgve(DEF_VECT,idime,jdime))         ! viscous fvins_nsi==1 leads to bulkv=1.0
                  end do
               end do
            end do
         end do
!
!    rho * (div u) u
!
         do inode = 1, pnode
            do idime = 1, ndime
               do igaus = 1, pgaus
                  elrbu(DEF_VECT, idime, inode) = elrbu(DEF_VECT, idime, inode) &
                                          - densi_aux*gpvol(DEF_VECT)*gpsha(DEF_VECT, inode, igaus)*gpvel(DEF_VECT, idime, igaus)* &
                                                  (gpgve(DEF_VECT, 1, 1) + gpgve(DEF_VECT, 2, 2) + gpgve(DEF_VECT, 3, 3))  ! Addit oriol
               end do
            end do
         end do

         !--------------------------------------------------------------------
         !
         ! Assembly in global syste,
         !
         !--------------------------------------------------------------------
         !
         ! Scatter element matrix to global one
         !
#ifndef OPENACCHHH
         call cputim(timeb)
         time1(7) = time1(7) + timeb - timea
         call cputim(timea)
         do ivect = 1, VECTOR_DIM
#endif
            ielem = list_elements(ivect)
            if (ielem > 0) then
               do inode = 1, pnode
                  ipoin = lnods_loc(ivect, inode)
                  do idime = 1, ndime
                     iauxi = idime + (ipoin - 1)*ndime
                     !$acc atomic update
#ifdef NO_COLORING
                     !$OMP ATOMIC
#endif
                     rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect, idime, inode)

                  end do

               end do
            end if
         end do
         !$acc end parallel loop

         !$acc exit data delete( lnods_loc , gpvol,                    &
         !$acc              elvel     ,  gpadv , gpvel ,                       &
         !$acc              gpgve   , elcod     ,  gpcar ,                     &
         !$acc              gpdet   ,                                          &
         !$acc              xjaci     ,                                        &
         !$acc              xmile,xnutu, G__ij    , Bbeta,              &
         !$acc              gprhs   , elrbu, visco_tot, aalph,gpsha  )

#ifndef OPENACCHHH
         call cputim(timeb)
         time1(8) = time1(8) + timeb - timea
#endif

      end subroutine nsi_element_operations_hh80

      end module mod_nsi_element_operations_hh80
