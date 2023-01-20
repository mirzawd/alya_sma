!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_updtcc_flamLet_dist(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtcc_flamLet_dist
  ! NAME
  !    chm_updtcc_flamLet_dist
  ! DESCRIPTION
  !    This routine call a subroutine that calcualates the critical time
  !      step size in flamelet models using a vectorized implementation
  ! USED BY
  !    chm_updtss
  !***
  !-----------------------------------------------------------------------
  use def_master,                 only : INOTMASTER
  use def_domain,                 only : lnnod, lgaus
  use def_kintyp,                 only : ip, rp
  use mod_parall,                 only : list_elements_par     !@
  use mod_parall,                 only : typ_list_elements_par !@
  use mod_parall,                 only : num_subd_par
  use mod_parall,                 only : num_pack_par
  use mod_communications, only : PAR_MIN


  implicit none
  real(rp),   intent(inout) :: dtmin


  integer(ip)                          :: pnode,pgaus
  integer(ip)                          :: isubd,ipack,ielem
  integer(ip)                          :: num_subd
  integer(ip)                          :: VECTOR_DIM
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)

  external                             :: chm_updtcc_flamLet_fast

  if( INOTMASTER ) then
    num_subd      =  num_subd_par
    num_pack      => num_pack_par
    list_elements => list_elements_par


    do isubd = 1,num_subd
      do ipack = 1,num_pack(isubd)

          ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
          pnode = lnnod(ielem)                                ! Number of nodes
          pgaus = lgaus(ielem)                                ! Number of Gauss points

          VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)

          call chm_updtcc_flamLet_fast(VECTOR_DIM,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,dtmin)

      end do
    end do

  end if

  !
  ! Look for minimum over subdomains
  !
  call PAR_MIN(dtmin,'IN MY CODE')

end subroutine chm_updtcc_flamLet_dist
