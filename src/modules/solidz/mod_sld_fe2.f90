!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_fe2.f90
!> @author  Guido Giuntoli
!> @date    July, 2018
!>
!> @brief   Interface to call MicroPP for FE2 multiscale calculations
!> @details
!>
!> @}
!------------------------------------------------------------------------

module mod_sld_fe2

#ifdef MICROPP
  use libmicropp
#endif
  use def_kintyp
  use def_domain
  use def_master
  use def_solidz
  use mod_sld_stress_model_comput
  use mod_parall,         only: PAR_COMM_MY_CODE
  use mod_communications, only: PAR_COMM_RANK_AND_SIZE
    
  implicit none
  !
  ! Parameters
  !
  integer(ip), parameter              ::  &
       MAT_MICRO_NO_COUPLING       = 666, &
       MAT_MICRO_ONE_WAY           = 667, &
       MAT_MICRO_FULL              = 668
    
  public                ::          &
    fe2_micropp_create            , &
    fe2_set_strains               , &
    fe2_needs_micropp             , &
    fe2_homogenize                , &
    fe2_update_vars               , &
    fe2_strain_tensor             , &
    fe2_get_stress_and_ctan       , &
    fe2_get_cost                  , &
    fe2_has_converged             , &
    fe2_is_non_linear             , &
    fe2_write_micropp             , &
    fe2_write_restart             , &
    fe2_read_restart              , &
    fe2_write_profiling           , &
    fe2_get_dF_gp

  public               ::           &
       MAT_MICRO_NO_COUPLING,       &
       MAT_MICRO_ONE_WAY,           &
       MAT_MICRO_FULL

#ifdef MICROPP
  type(micropp3)      :: micro
  type(material_base) :: mat_params(3)
  integer(ip)         :: ngp_total
  integer(ip)         :: ngp_per_elem
#endif

  private

  contains

!-----------------------------------------------------------------------
!>
!> @brief Inits the micropp class to start performing micropp
!>        calculations.
!>
!-----------------------------------------------------------------------

subroutine fe2_micropp_create( )

  implicit none

#ifdef MICROPP
  integer(ip) :: sizes(3), micro_type
  integer(ip) :: imate, needs_micropp
  integer(ip) :: pmate
  integer(4)  :: M1, M2, M3
  integer(ip) :: i, j
  real(rp)    :: micro_params(4), params(20)
  real(rp)    :: E1, nu1, Ka1, Sy1, Xt1
  real(rp)    :: E2, nu2, Ka2, Sy2, Xt2
  real(rp)    :: E3, nu3, Ka3, Sy3, Xt3
  integer(ip) :: mpi_rank, sizempi
  integer(ip), allocatable :: coupling(:)
  integer(ip) :: num_no_coupling
  integer(ip) :: num_one_way
  integer(ip) :: num_full

  call fe2_needs_micropp (imate, needs_micropp)
  if (needs_micropp == 0 .or. IMASTER) then
     return
  end if

  sizes = (/ parco_sld(1, imate), parco_sld(2, imate), parco_sld(3, imate) /)
  micro_type = parco_sld(4, imate)
  do i = 1, 4
     micro_params(i) = parco_sld(4 + i, imate)
  end do

  M1 = parco_sld(9, imate)
  M2 = parco_sld(10, imate)
  M3 = parco_sld(11, imate)

  E1  = parco_sld(12, imate)
  nu1 = parco_sld(13, imate)
  Ka1 = parco_sld(14, imate)
  Sy1 = parco_sld(15, imate)
  Xt1 = parco_sld(16, imate)

  E2  = parco_sld(17, imate)
  nu2 = parco_sld(18, imate)
  Ka2 = parco_sld(19, imate)
  Sy2 = parco_sld(20, imate)
  Xt2 = parco_sld(21, imate)

  E3  = parco_sld(22, imate)
  nu3 = parco_sld(23, imate)
  Ka3 = parco_sld(24, imate)
  Sy3 = parco_sld(25, imate)
  Xt3 = parco_sld(26, imate)

  call material_set(mat_params(1), M1, E1, nu1, Ka1, Sy1, Xt1)
  call material_set(mat_params(2), M2, E2, nu2, Ka2, Sy2, Xt2)
  call material_set(mat_params(3), M3, E3, nu3, Ka3, Sy3, Xt3)

  ngp_per_elem = lgaus(1) ! Number of Gauss points
  ngp_total = nelem * ngp_per_elem
  !allocate coupling vector
  allocate(coupling(ngp_total))
  call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE, mpi_rank, sizempi)
  ! call set_coupling_vector(coupling,ngp_total,)

  num_no_coupling = 0
  num_one_way = 0
  num_full = 0
  do i = 1, nelem

      pmate = lmate_sld(i)
      if (lawst_sld(pmate)==MAT_MICRO_NO_COUPLING) then

        do j =1,ngp_per_elem
           coupling((i - 1) * ngp_per_elem + j) = 0
           num_no_coupling = num_no_coupling + 1
        end do
           
     else if (lawst_sld(pmate)==MAT_MICRO_ONE_WAY) then

        do j =1,ngp_per_elem
           coupling((i - 1) * ngp_per_elem + j) = 1
           num_one_way = num_one_way + 1
        end do

     else if (lawst_sld(pmate)==MAT_MICRO_FULL) then
 
        do j =1,ngp_per_elem
           coupling((i - 1) * ngp_per_elem + j) = 2
           num_full = num_full + 1
        end do

     end if
   
  end do

  print *, "RANK = ", mpi_rank, " nelem = ", nelem, " ngp_total = ", ngp_total
  write(6,*) "RANK = ", mpi_rank, "NO_COUPLING = ", num_no_coupling, &
          "ONE_WAY = ", num_one_way, "FULL = ", num_full
  
  call micropp3_new(micro, ngp_total, sizes, micro_type, micro_params, &
          mat_params, coupling, 1, mpi_rank)

  !if (INOTSLAVE) then
     call micropp3_print_info(micro)
  !end if
  deallocate(coupling)
#endif
        
end subroutine fe2_micropp_create

!-----------------------------------------------------------------------
!>
!> @brief   Returns the average stress
!>          from micropp.
!>
!> @param[in] elem_id
!> @param[in] num_gp
!> @param[out] Sig
!> @param[out] Ct
!>
!-----------------------------------------------------------------------

subroutine fe2_get_stress_and_ctan (e, num_gp, Sig, Ct)

  implicit none

  integer(ip), intent(in)    :: num_gp
  integer(ip), intent(in)    :: e
  real(rp),    intent(out)   :: Sig (ndime, ndime, num_gp, 2)
  real(rp),    intent(out)   :: Ct (ndime, ndime, ndime, ndime, num_gp)
#ifdef MICROPP
  real(rp)                   :: aux, Sig_voi (6)
  real(rp)                   :: Ct_gp (ndime, ndime, ndime, ndime)
  real(rp)                   :: Ct_voi (36), Ct_voi_mat (6,6)
  integer(ip)                :: gp, i, j, gp_id

  do gp = 1, num_gp

     gp_id = (e - 1) * ngp_per_elem + (gp - 1)

     call micropp3_get_stress(micro, gp_id, Sig_voi)

     Sig(1,1,gp,1) = Sig_voi(1)
     Sig(2,1,gp,1) = Sig_voi(6)
     Sig(3,1,gp,1) = Sig_voi(5)

     Sig(1,2,gp,1) = Sig_voi(6)
     Sig(2,2,gp,1) = Sig_voi(2)
     Sig(3,2,gp,1) = Sig_voi(4)

     Sig(1,3,gp,1) = Sig_voi(5)
     Sig(2,3,gp,1) = Sig_voi(4)
     Sig(3,3,gp,1) = Sig_voi(3)

     call micropp3_get_ctan(micro, gp_id, Ct_voi)

     do i = 1, nvoig_sld
      do j = 1, nvoig_sld
       Ct_voi_mat (i,j) = Ct_voi((i-1)*nvoig_sld + j)
      end do
     end do

     aux = Ct_voi_mat(4,4)
     Ct_voi_mat(4,4) = Ct_voi_mat(6,6)
     Ct_voi_mat(6,6) = aux

     call SM_tensor_to_voigt_fourth(ndime, 1_ip, Ct_gp, Ct_voi_mat) ! Ct_voi -> Ct

     Ct(:,:,:,:,gp) = Ct_gp

  end do
#endif

end subroutine fe2_get_stress_and_ctan


!-----------------------------------------------------------------------
!>
!> @brief   Returns the total cost per element MicroPP
!>          integration point
!>
!> @param[in] e
!> @param[out] cost
!>
!-----------------------------------------------------------------------

function fe2_get_cost(e)

  implicit none
  integer(ip)                :: e
  integer(ip)                :: fe2_get_cost
#ifdef MICROPP
  integer(ip)                :: cost_a
  integer(ip)                :: gp
  integer(ip)                :: gp_id
#endif

  fe2_get_cost = 0_ip

#ifdef MICROPP
  do gp = 1, ngp_per_elem
     gp_id = (e - 1) * ngp_per_elem + (gp - 1)
     cost_a = micropp3_get_cost(micro, gp_id)
     fe2_get_cost = fe2_get_cost + cost_a
  end do
#endif

end function fe2_get_cost

!-----------------------------------------------------------------------
!>
!> @brief   Returns if an element has at least one integration 
!>          point which has not converged
!>
!> @param[in] e
!> @param[out] convergence
!>
!-----------------------------------------------------------------------

function fe2_has_converged(e)

  implicit none
  integer(ip)  :: e
  logical      :: fe2_has_converged
#ifdef MICROPP
  logical      :: converged
  integer(ip)  :: gp
  integer(ip)  :: gp_id
#endif

  fe2_has_converged = .true.

#ifdef MICROPP
  do gp = 1, ngp_per_elem
     gp_id = (e - 1) * ngp_per_elem + (gp - 1)
     converged = micropp3_has_converged(micro, gp_id)
     fe2_has_converged = (fe2_has_converged .and. converged)
  end do
#endif

end function fe2_has_converged


!-----------------------------------------------------------------------
!>
!> @brief   Returns if the element has at least one non-linear
!>          integration point
!>
!> @param[in] e
!> @param[out] nl_flag
!>
!-----------------------------------------------------------------------

function fe2_is_non_linear(e)

  implicit none
  integer(ip)    :: e
  logical        :: fe2_is_non_linear
#ifdef MICROPP
  logical        :: non_linear
  integer(ip)    :: gp
  integer(ip)    :: gp_id
#endif

  fe2_is_non_linear = .false.

#ifdef MICROPP
  do gp = 1, ngp_per_elem
     gp_id = (e - 1) * ngp_per_elem + (gp - 1)
     non_linear = micropp3_is_non_linear(micro, gp_id)
     fe2_is_non_linear = (fe2_is_non_linear .or. non_linear)
  end do
#endif

end function fe2_is_non_linear


!-----------------------------------------------------------------------
!>
!> @brief Returns a flag that says if Alya needs micropp, at least one
!>        element should be MICRO type.
!>
!> @param[in] imate
!> @param[out] needs_micropp
!>
!-----------------------------------------------------------------------

subroutine fe2_needs_micropp (imate, needs_micropp)

  implicit none

  integer(ip), intent(out) :: needs_micropp, imate
  integer(ip)              :: i

  needs_micropp = 0
  do i = 1,nmate_sld
    if (lawst_sld(i) == MAT_MICRO_NO_COUPLING .or. &
        lawst_sld(i) == MAT_MICRO_ONE_WAY .or. &
        lawst_sld(i) == MAT_MICRO_FULL) then
      needs_micropp = 1
      imate = i
      exit
    end if
  end do

end subroutine fe2_needs_micropp

!-----------------------------------------------------------------------
!>
!> @brief Goes through all elements and Gauss points and calculates the
!>        strains passing it to micropp.
!>
!-----------------------------------------------------------------------

subroutine fe2_set_strains()

  implicit none

#ifdef MICROPP
  integer(ip) :: e, gp, gp_id
  integer(ip) :: imate, needs_micropp
  real(rp)    :: strain_gp(3,3), strain_gp_voi(6)

  call fe2_needs_micropp (imate, needs_micropp)

  if (needs_micropp == 1) then

    do e = 1, nelem
      do gp = 1, ngp_per_elem
        gp_id = (e - 1) * ngp_per_elem + (gp - 1)
        call fe2_strain_tensor (e, gp, strain_gp)
        strain_gp_voi(1) = strain_gp(1,1)
        strain_gp_voi(2) = strain_gp(2,2)
        strain_gp_voi(3) = strain_gp(3,3)
        strain_gp_voi(4) = strain_gp(2,3) * 2.0_rp
        strain_gp_voi(5) = strain_gp(1,3) * 2.0_rp
        strain_gp_voi(6) = strain_gp(1,2) * 2.0_rp
        call micropp3_set_strain(micro, gp_id, strain_gp_voi)
      end do
    end do

  end if
#endif

end subroutine fe2_set_strains

!-----------------------------------------------------------------------
!>
!> @brief Orders micropp to start converting the strains at each
!>        Gauss point into an average stress on the microstructure
!>
!-----------------------------------------------------------------------

subroutine fe2_homogenize()

#ifdef MICROPP
  integer(ip) :: imate, needs_micropp

  call fe2_needs_micropp (imate, needs_micropp)

  if (needs_micropp == 0 .or. IMASTER) then
          return
  end if

  call micropp3_homogenize(micro)
#endif

end subroutine fe2_homogenize

!-----------------------------------------------------------------------
!>
!> @brief Ask micropp to update internal variables.
!>
!-----------------------------------------------------------------------

subroutine fe2_update_vars()

#ifdef MICROPP

  integer(ip) :: imate, needs_micropp

  call fe2_needs_micropp (imate, needs_micropp)

  if (needs_micropp == 0 .or. IMASTER) then
          return
  end if

  call micropp3_update_vars(micro)

#endif

end subroutine fe2_update_vars

!-----------------------------------------------------------------------
!>
!> @brief Returns a the strain in the Gauss point of an element
!>
!> @param[in] e
!> @param[in] gp
!> @param[out] strain_gp
!>
!-----------------------------------------------------------------------

subroutine fe2_strain_tensor (e, gp, strain_gp)

  implicit none

  integer(ip),    intent(in)  :: e, gp
  real(rp),       intent(out) :: strain_gp(3,3)
  real(rp)                    :: dF(3,3)

  call fe2_get_dF_gp (e, gp, dF)
  call SM_strain_tensor (0_ip, dF, strain_gp)

end subroutine fe2_strain_tensor

!-----------------------------------------------------------------------
!>
!> @brief Returns the deformation gradient at a gauss point of
!>        an element.
!>
!> @param[in] e
!> @param[in] gp
!> @param[out] dF_gp
!>
!-----------------------------------------------------------------------

subroutine fe2_get_dF_gp(e, gp, dF_gp)

  implicit none

  integer(ip), intent(in)    :: e, gp
  real(rp)   , intent(out)   :: dF_gp(3,3)
  integer(ip)                :: i, j, n, pelty
  integer(ip)                :: pgaus, pnode, ipoin
  real(rp)                   :: eldis(3,8)
  real(rp)                   :: dsh_gp(3,8)

  pelty = ltype(e)
  pnode = nnode(pelty)
  pgaus = ngaus(pelty)

  do n = 1, pnode
    ipoin = lnods(n, e)
    eldis(1:ndime,n) = displ(1:ndime,ipoin,ITER_K)
  end do

  call fe2_get_dsh_gp (e, gp, dsh_gp)

  dF_gp = 0.0_rp
  dF_gp(1,1) = 1.0_rp
  dF_gp(2,2) = 1.0_rp
  dF_gp(3,3) = 1.0_rp
  do i=1,ndime
   do j=1,ndime
    do n=1,pnode
      dF_gp(i,j) = dF_gp(i,j) + eldis(i,n) * dsh_gp(j,n)
    end do
   end do
  end do

end subroutine fe2_get_dF_gp

!-----------------------------------------------------------------------
!>
!> @brief Returns the derivatives of the shape functions at a Gauss
!>        point of an element
!>
!> @param[in] e
!> @param[in] gp
!> @param[out] dsh_gp
!>
!-----------------------------------------------------------------------

subroutine fe2_get_dsh_gp (e, gp, dsh_gp)

  implicit none

  integer(ip), intent(in)    :: e, gp
  real(rp)   , intent(out)   :: dsh_gp(3,8)
  integer(ip)                :: i, j, n, pelty, pnode, ipoin
  real(rp)                   :: jac(3,3), ijac(3,3), dsh_0(3,8,8)
  real(rp)                   :: c11, c21, c31, c12, c22, c32, c13, c23, c33, det
  real(rp)                   :: elem_coor(3,8)

  pelty = ltype(e)
  pnode = nnode(pelty)

  dsh_0 = elmar(pelty)%deriv

  do n = 1, pnode
    ipoin = lnods(n, e)
    elem_coor(1:ndime,n) = coord(1:ndime,ipoin)
  end do

  jac = 0.0_rp
  do i=1,ndime
     do j=1,ndime
        do n=1,pnode
           jac(i,j) = jac(i,j) + elem_coor(j,n) * dsh_0(i,n,gp)
        end do
     end do
  end do

  c11 = +jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3); c12 = -jac(2,1)*jac(3,3)+jac(3,1)*jac(2,3); c13 = +jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2);
  c21 = -jac(1,2)*jac(3,3)+jac(3,2)*jac(1,3); c22 = +jac(1,1)*jac(3,3)-jac(3,1)*jac(1,3); c23 = -jac(1,1)*jac(3,2)+jac(3,1)*jac(1,2);
  c31 = +jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3); c32 = -jac(1,1)*jac(2,3)+jac(2,1)*jac(1,3); c33 = +jac(1,1)*jac(2,2)-jac(2,1)*jac(1,2);

  det = jac(1,1)*c11 + jac(1,2)*c12 + jac(1,3)*c13;

  ijac(1,1) = c11/det; ijac(1,2) = c21/det; ijac(1,3) = c31/det;
  ijac(2,1) = c12/det; ijac(2,2) = c22/det; ijac(2,3) = c32/det;
  ijac(3,1) = c13/det; ijac(3,2) = c23/det; ijac(3,3) = c33/det;

  do n=1,pnode
     do i=1,ndime
        dsh_gp(i,n) = 0.0_rp
        do j=1,ndime
           dsh_gp(i,n) = dsh_gp(i,n) + ijac(i,j) * dsh_0(j,n,gp)
        end do
     end do
  end do

end subroutine fe2_get_dsh_gp

!-----------------------------------------------------------------------
!>
!> @brief This sub create and write the output files for the micro-scale
!>
!-----------------------------------------------------------------------

subroutine fe2_write_micropp(ielem, ielem_global, time_step)

        integer(ip)  , intent(in) :: ielem
        integer(ip)  , intent(in) :: ielem_global
        integer(ip)  , intent(in) :: time_step
#ifdef MICROPP
        integer(ip)               :: gp_id


        gp_id = (ielem - 1) * ngp_per_elem + 0
        
        call micropp3_output2(micro, gp_id, ielem_global, time_step)
#endif

end subroutine fe2_write_micropp

!-----------------------------------------------------------------------
!>
!> @brief Writes the restart files for the micro-scale
!>
!-----------------------------------------------------------------------

subroutine fe2_write_restart()

   integer(ip) :: imate, needs_micropp

   call fe2_needs_micropp (imate, needs_micropp)
   if (needs_micropp == 0 .or. IMASTER) then
           return
   end if

#ifdef MICROPP
   call micropp3_write_restart(micro, 16)
#endif

end subroutine fe2_write_restart


!-----------------------------------------------------------------------
!>
!> @brief Reads the restart files for the micro-scale
!>
!-----------------------------------------------------------------------

subroutine fe2_read_restart()

   integer(ip) :: imate, needs_micropp

   call fe2_needs_micropp (imate, needs_micropp)
   if (needs_micropp == 0 .or. IMASTER) then
           return
   end if

#ifdef MICROPP
      call micropp3_read_restart(micro, 16)
#endif

end subroutine fe2_read_restart


!-----------------------------------------------------------------------
!>
!> @brief Writes the micropp profiling
!>
!-----------------------------------------------------------------------

subroutine fe2_write_profiling(time_step)

   integer(ip), intent(in) :: time_step
   integer(ip)             :: imate, needs_micropp

   call fe2_needs_micropp (imate, needs_micropp)
   if (needs_micropp == 0 .or. IMASTER) then
           return
   end if

#ifdef MICROPP
      call micropp3_write_profiling(micro, time_step)
#endif

end subroutine fe2_write_profiling

end module mod_sld_fe2
