!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outmat.f90
!> @author  Gerard Guillamet
!> @date    November, 2022
!>          - Subroutine written
!> @brief   Initialization of some materials for the pre-computing of properties
!> @details
!> @}
!------------------------------------------------------------------------

subroutine sld_inimat

  use def_kintyp,                    only : ip, rp, lg
  use def_master,                    only : INOTEMPTY
  use def_elmtyp,                    only : ELINT, ELFEM
  use def_domain,                    only : lmate, ltype, lnnod, lnods, lorde, lelch
  use def_domain,                    only : elmar, mnode, coord, nelem, ndime
  use def_domain,                    only : hnatu
  use def_solidz,                    only : kfl_ellen_sld
  use def_solidz,                    only : lawst_sld, celen_sld
  use mod_sld_vect_stress_models,    only : sld_stress_model_151_compute_properties
  use mod_sld_vect_stress_model_154, only : sld_stress_model_154_initialisations
  use mod_sld_vect_stress_model_154, only : sld_stress_model_154_compute_properties
  use mod_sld_vect_stress_model_154, only : sld_stress_model_154_message
  
  implicit none

  external    :: elmlen
  external    :: elmchl

  integer(ip) :: ielem, inode, ipoin
  integer(ip) :: pelty, pnode, porde, pmate 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: tragl(9),hleng(3)
  real(rp)    :: chale(2),dumm1(ndime,mnode),dumm2(ndime,2)
  integer(ip) :: ireducedXT
  integer(ip) :: ireducedXTO
  integer(ip) :: ireducedXC
  integer(ip) :: ireducedXCO
  integer(ip) :: ireducedYT
  integer(ip) :: ireducedYC
  integer(ip) :: ireducedSL
  !
  ! Initializations
  !
  if( any(lawst_sld(:) == 154_ip) ) then
     call sld_stress_model_154_initialisations(&
          ireducedXT,ireducedXTO,ireducedXC,ireducedXCO,ireducedYT,ireducedYC,ireducedSL)
  end if

  if( INOTEMPTY ) then

     do ielem = 1,nelem

        pelty = ltype(ielem)
        pnode = lnnod(ielem)
        porde = lorde(pelty)
        pmate = lmate(ielem)
               
        if( pelty > 0 .and. lelch(ielem) == ELFEM ) then
           !
           ! Element coordinates
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
           !
           ! Characteristic element length 
           !
           call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
           call elmchl(tragl,hleng,elcod,dumm1,dumm2,chale,pelty,pnode,porde,hnatu(pelty),0_ip,kfl_ellen_sld)
           !
           ! Material properties computations of some models
           !
           if(      lawst_sld(pmate) == 151_ip ) then
              
              call sld_stress_model_151_compute_properties(pmate)
              
           else if( lawst_sld(pmate) == 154_ip ) then
              
              call sld_stress_model_154_compute_properties(ielem,pmate,chale(1),&
                   ireducedXT,ireducedXTO,ireducedXC,ireducedXCO,ireducedYT,ireducedYC,ireducedSL)
              
           end if

        else if( lelch(ielem) == ELINT ) then
           chale(1) = 1.0_rp

        end if
        !
        ! Save global variables
        !
        celen_sld(ielem) = chale(1)
        
     end do

  end if
  !
  ! Live info
  !
  if( any(lawst_sld(:) == 154_ip) ) then
     call sld_stress_model_154_message(&
          ireducedXT,ireducedXTO,ireducedXC,ireducedXCO,ireducedYT,ireducedYC,ireducedSL)
  end if

end subroutine sld_inimat
