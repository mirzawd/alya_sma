!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    mod_gus_transmssion_conditions.f90
!> @author  houzeaux
!> @date    2020-10-22
!> @brief   Compute and impose transmission conditions
!> @details ???
!>          ???
!-----------------------------------------------------------------------

  module mod_gus_transmssion_conditions

    use def_kintyp_basic, only : ip,rp,i1p
    use def_elmtyp
    use def_master
    use def_domain
    use mod_elmgeo
    use mod_communications_point_to_point
    use mod_memory_basic
    use def_gusano
    implicit none

    private

    public :: gus_transmssion_conditions
    public :: gus_transmssion_conditions_initialization

  contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-22
  !> @brief   Transmission conditions initialization
  !> @details Transmission conditions initializatio
  !> 
  !-----------------------------------------------------------------------
  
  subroutine gus_transmssion_conditions_initialization()

    integer(ip)               :: ipoin,ielem,inode
    integer(ip)               :: ipoin_dir,pelty
    !
    ! Dirichlet
    !
    do ielem = 1,nelem
       pelty = ltype(ielem)
       if( pelty == DDDNE ) then
          ipoin_dir = lnods(1,ielem)
          call memory_alloca(mem_modul(1:2,modul),'DIRICH_GUS % L','gus_memall',dirich_gus(ipoin_dir) % l,elmgeo_number_nodes(pelty,lnods(:,ielem))-1_ip)
          do inode = 2,elmgeo_number_nodes(pelty,lnods(:,ielem))
             ipoin = lnods(inode,ielem)
             dirich_gus(ipoin_dir) % l(inode-1) = ipoin
          end do
       end if
    end do
    !
    ! Neumann
    !
    do ielem = 1,nelem
       pelty = ltype(ielem)
       if( pelty == DDDNE ) then
          ipoin_dir = lnods(1,ielem)
          do inode = 2,elmgeo_number_nodes(pelty,lnods(:,ielem))
             ipoin = lnods(inode,ielem)
             neuman_gus(ipoin) = ipoin_dir
          end do
       end if
    end do
    
  end subroutine gus_transmssion_conditions_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-22
  !> @brief   Assemble Dirichlet
  !> @details Assemble Dirichlet condition
  !> 
  !-----------------------------------------------------------------------
  
  subroutine gus_transmssion_conditions(aa,bb)

    real(rp),   pointer, intent(inout) :: aa(:)
    real(rp),   pointer, intent(inout) :: bb(:)

    if( associated(aa) .and. associated(bb) ) call gus_transmssion_conditions_go(aa,bb)
    
  end subroutine gus_transmssion_conditions
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-22
  !> @brief   Assemble Dirichlet
  !> @details Assemble Dirichlet condition
  !> 
  !-----------------------------------------------------------------------
  
  subroutine gus_transmssion_conditions_go(aa,bb)

    real(rp),   intent(inout) :: aa(2,2,*)
    real(rp),   intent(inout) :: bb(2,*)
    integer(ip)               :: ipoin,iz,jpoin,inode
    !
    ! Dirichlet velocity Q_dir + Q2 + Q3  = 0 => A_dir u_dir + A_2 u_2 + A_3 u_3 = 0
    !    
    do ipoin = 1,npoin
       if( kfl_fixno_gus(1,ipoin) == -1 ) then
          bb(1,ipoin) = 0.0_rp
          do iz = r_dom(ipoin),r_dom(ipoin+1)-1
             jpoin = c_dom(iz)
             aa(1,1,iz) = 0.0_rp
             aa(2,1,iz) = 0.0_rp
          end do
       end if
    end do
    
    do ipoin = 1,npoin
       if( associated(dirich_gus(ipoin) % l) ) then
          if( kfl_fixno_gus(1,ipoin) == -1 ) then
             do iz = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(iz)
                if( jpoin == ipoin ) then
                   !aa(1,1,iz) = areas_gus(ipoin) * exn1d_gus(jpoin)
                   aa(1,1,iz) = exn1d_gus(jpoin)
                else
                   do inode = 1,memory_size(dirich_gus(ipoin) % l)
                      if( dirich_gus(ipoin) % l(inode) == jpoin ) then
                         !aa(1,1,iz) = areas_gus(jpoin) * exn1d_gus(jpoin)
                         aa(1,1,iz) = exn1d_gus(jpoin)
                      end if
                   end do
                end if
             end do
          end if
       end if
    end do
    !
    ! Neumann pressure. (L(u),v) + p_dir n = 0
    !
    do ipoin = 1,npoin
       if( kfl_fixno_gus(1,ipoin) == -2 .and. neuman_gus(ipoin) /= 0 ) then
          do iz = r_dom(ipoin),r_dom(ipoin+1)-1
             jpoin = c_dom(iz)
             if(      neuman_gus(ipoin)==  jpoin ) then
                aa(2,1,iz) = exn1d_gus(ipoin)
             end if
          end do
       end if
    end do
        
  end subroutine gus_transmssion_conditions_go

end module mod_gus_transmssion_conditions
!> @}
