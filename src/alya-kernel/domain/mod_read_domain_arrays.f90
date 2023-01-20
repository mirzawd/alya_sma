!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_readom.f90
!> @author  houzeaux
!> @date    2018-11-16
!> @brief   Module for domain reading
!> @details Bridge to different readings
!-----------------------------------------------------------------------

module mod_read_domain_arrays

  use def_kintyp,               only : ip,rp
  use def_master,               only : cpu_start
  use def_master,               only : CPU_READ_GEO
  use def_master,               only : CPU_READ_SETS
  use def_master,               only : CPU_READ_BCS
  use def_master,               only : CPU_READ_FIELDS
  use def_master,               only : intost
  use def_domain,               only : nelty
  use def_inpout,               only : words
  use def_inpout,               only : param
  use def_inpout,               only : nunit
  use def_inpout,               only : getcha

  use mod_mpio_config,          only : mpio_config

  use mod_reaset,               only : reaset_seq
  use mod_reabcs,               only : reabcs_seq
  use mod_reafie,               only : reafie_seq

  use mod_mpio_par_readom,      only : reageo_par_mpio
  use mod_mpio_par_readom,      only : reaset_par_mpio
  use mod_mpio_par_readom,      only : reabcs_par_mpio
  use mod_mpio_par_readom,      only : reafie_par_mpio
  use mod_mpio_par_readom,      only : par_readom_finalize

  use mod_mpio_seq_readom,      only : reageo_seq_mpio
  use mod_mpio_seq_readom,      only : reaset_seq_mpio
  use mod_mpio_seq_readom,      only : reabcs_seq_mpio
  use mod_mpio_seq_readom,      only : reafie_seq_mpio

  use mod_elmgeo,               only : elmgeo_element_name_to_type
  use mod_ecoute,               only : ecoute
 
  implicit none

  private
  
  real(rp) :: time1,time2

  public :: read_domain_arrays_reageo
  public :: read_domain_arrays_reaset
  public :: read_domain_arrays_reabcs
  public :: read_domain_arrays_reafie
  public :: read_domain_arrays_finalize
  public :: read_domain_arrays_types
  
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Read geometry
  !> @details Read geometry in one of the four possible modes
  !>
  !-----------------------------------------------------------------------

  subroutine read_domain_arrays_reageo()

    call cputim(time1)

    if (mpio_config%input%enabled) then

        if (mpio_config%input%parallel) then
            
            call reageo_par_mpio()

        else

            call reageo_seq_mpio()

        end if

    else

       call reageo_seq()

    end if


    call cputim(time2)
    cpu_start(CPU_READ_GEO) = time2 - time1

  end subroutine read_domain_arrays_reageo

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Read sets
  !> @details Read sets in one of the four possible modes
  !>
  !-----------------------------------------------------------------------

  subroutine read_domain_arrays_reaset()

    call cputim(time1)

    if (mpio_config%input%enabled) then

        if (mpio_config%input%parallel) then
            
            call reaset_par_mpio()

        else

            call reaset_seq_mpio()

        end if

    else

       call reaset_seq()

    end if


    call cputim(time2)
    cpu_start(CPU_READ_SETS) = time2 - time1

  end subroutine read_domain_arrays_reaset

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Read boundary conditions
  !> @details Read boundary conditions in one of the four possible modes
  !>
  !-----------------------------------------------------------------------

  subroutine read_domain_arrays_reabcs()

    call cputim(time1)

    if (mpio_config%input%enabled) then

        if (mpio_config%input%parallel) then
            
            call reabcs_par_mpio()

        else

            call reabcs_seq_mpio()

        end if

    else

       call reabcs_seq()

    end if

    call cputim(time2)
    cpu_start(CPU_READ_BCS) = time2 - time1

  end subroutine read_domain_arrays_reabcs

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Read fields
  !> @details Read fields in one of the four possible modes
  !>
  !-----------------------------------------------------------------------

  subroutine read_domain_arrays_reafie()

    call cputim(time1)

    if (mpio_config%input%enabled) then

        if (mpio_config%input%parallel) then
            
            call reafie_par_mpio()

        else

            call reafie_seq_mpio()

        end if

    else

       call reafie_seq()

    end if

    call cputim(time2)
    cpu_start(CPU_READ_FIELDS) = time2 - time1

  end subroutine read_domain_arrays_reafie

  subroutine read_domain_arrays_finalize()
      use def_domain,                    only : nfiel, xfiel
      use mod_memory,                    only : memory_size
      implicit none   

      if (mpio_config%input%parallel) then

          call par_readom_finalize()

      end if
  end subroutine read_domain_arrays_finalize

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-17
  !> @brief   Read types
  !> @details Read types
  !> 
  !-----------------------------------------------------------------------

  subroutine read_domain_arrays_types(kfl_icgns,neleg,ktype,ltypg,lexig)

    integer(ip), intent(in)  :: kfl_icgns,neleg
    integer(ip), intent(out) :: ktype
    integer(ip), intent(out) :: ltypg(:)
    integer(ip), intent(out) :: lexig(:)
    integer(ip)              :: ielty,ielem,dummi
    integer(ip)              :: kelem
    character(5)             :: welem

    if( words(2) == 'ALL  ' ) then
       !
       ! All elements are of the same type
       !
       ktype = 1
       ielty = int(param(2))
       if( ielty == 0 ) then
          welem = getcha('ALL  ','NULL ','#Element type')
          call elmgeo_element_name_to_type(welem,ielty)
          if( ielty >= 2 .and. ielty <= nelty ) lexig(ielty)=1
       end if
       do ielem = 1,neleg
          ltypg(ielem) = ielty
       end do
       ktype = neleg
       do while( words(1) /= 'ENDTY' )
          call ecoute('reageo')
       end do

    else if( words(2) == 'NONOR' ) then
       ktype=0
       if(kfl_icgns==1) then
          !
          ! CGNS type, old Alya type
          !
          do ielem=1,neleg
             read(nunit,*,err=1) kelem,ielty
             ielty=ltnew(ielty)
             lexig(ielty)=1
             ltypg(kelem)=ielty
          end do
          ktype=neleg
       else
          !
          ! Alya type as defined in def_elmtyp
          !
          do ielem=1,neleg
             read(nunit,*,err=1) kelem,ielty
             lexig(ielty)=1
             ltypg(kelem)=ielty
          end do
          ktype=neleg        
       end if
       call ecoute('reageo')
       if(words(1)/='ENDTY')&
            call runend('REAGEO: WRONG TYPES OF ELEMENT FIELD')

    else

       ktype=0
       if(kfl_icgns==1) then
          !
          ! CGNS type, old Alya type
          !
          do ielem=1,neleg
             read(nunit,*,err=1) dummi,ielty
             ielty=ltnew(ielty)
             lexig(ielty)=1
             ltypg(ielem)=ielty
          end do
          ktype=neleg
       else
          !
          ! Alya type as defined in def_elmtyp
          !
          do ielem = 1,neleg
             read(nunit,*,err=1) dummi,ielty
             if( ielty > 0 ) lexig(ielty)=1
             ltypg(ielem) = ielty
          end do
          ktype = neleg        
       end if
       call ecoute('reageo')
       if(words(1)/='ENDTY')&
            call runend('REAGEO: WRONG TYPES OF ELEMENT FIELD')
    end if

    return

1   call runend('REAGEO: WRONG NUMBER OF NODES FOR ELEMENT '//trim(intost(ielem)))

  end subroutine read_domain_arrays_types

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-17
  !> @brief   Converts old type to new type
  !> @details CGNS type
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function ltnew(ityol)
    
    integer(ip) :: ityol

    if (ityol ==  2) ltnew= 2
    if (ityol ==  3) ltnew= 3
    if (ityol ==  4) ltnew= 4
    if (ityol ==  5) ltnew= 10
    if (ityol ==  6) ltnew= 11
    if (ityol ==  7) ltnew= 12
    if (ityol ==  8) ltnew= 13
    if (ityol ==  9) ltnew= 14
    if (ityol == 10) ltnew= 30
    if (ityol == 11) ltnew= 31
    if (ityol == 12) ltnew= 32
    if (ityol == 13) ltnew= 33
    if (ityol == 14) ltnew= 34
    if (ityol == 15) ltnew= 35
    if (ityol == 16) ltnew= 36
    if (ityol == 17) ltnew= 37
    if (ityol == 18) ltnew= 38
    if (ityol == 19) ltnew= 39  

  end function ltnew

end module mod_read_domain_arrays
!> @}
