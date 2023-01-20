!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_updtss.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Computes the proper time step size 
!> @details Computes the proper time step size 
!> @} 
!-----------------------------------------------------------------------
subroutine exm_updtss
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi
  use mod_eccoupling,          only : kfl_exmsld_ecc
  use mod_communications,      only : PAR_MIN
  use mod_exm_diffusivity
  implicit none 

  integer(ip) :: ielem,pmate,idime,jdime,inode
  integer(ip) :: pnode,ipoin,pelty,igaus,noion
  real(rp)    :: dtcri,dtmin,dtaux,hlmin,cndmx,xfact
  real(rp)    :: difex(ndime,ndime,mgaus), difin(ndime,ndime,mgaus)
  real(rp)    :: hleng(ndime),tragl(ndime,ndime),elcod(ndime,mnode)
  !
  ! Time step computed from the critical one. : 
  !
  ! The time step is computed from the diffusive limit... this should be
  ! IMPROVED to take into account:
  ! 
  ! 1. reaction and non-linear terms
  ! 2. tensorial character of the conduction
  ! 3. coupling
  !
  if (kfl_timco == 0 ) then                             ! Externally fixed time step
     dtmin = dtext_exm
     dtcri = dtmin
     return
  else
     dtmin = 0.0_rp
     dtcri = 1000000000.0_rp
  end if

  call times(2) % ini()

  if(INOTMASTER) then

     difex = 0.0_rp
     difin = 0.0_rp


     xfact =  0.5_rp 

!!!! ojo: corregir este xfact para elementos cuadraticos...

     xfact = xfact * safet_exm


     elements: do ielem = 1,nelem

        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pmate = 1
        if( nmate > 1 ) then
           pmate = lmate(ielem)
        end if

        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime= 1,ndime
              elcod(idime,inode)  = coord(idime,ipoin)
           end do
           if( kfl_exmsld_ecc ) then
              if (kfl_gcoup_exm == 1) then
                 do idime=1,ndime
                    elcod(idime,inode)= elcod(idime,inode) + displ(idime,ipoin,1)
                 end do
              end if
           end if

        end do

        igaus = 1

        ! 1. calcular la hleng mas corta de cada elemento...      

        call elmlen(ndime,pnode,elmar(pelty)%dercg(1,1),tragl,elcod,hnatu(pelty),hleng)
        hlmin = hleng(ndime)

        ! 2. calcular la conductividad maxima de cada elemento... 
        ! Was commented for single material/conductivity problems, except for multiple conductivities it does make sense
        call exm_diffusivity_at_gp_comcnd(-ielem,pmate,difin,noion)

        cndmx=0.0_rp
        do jdime=1,ndime
           do idime=1,ndime
              if (difin(idime,jdime,1) >= cndmx) cndmx = difin(idime,jdime,1)
              if (kfl_cemod_exm == 2) then
                 if (difex(idime,jdime,1) >= cndmx) cndmx = difex(idime,jdime,1)
              end if
           end do
        end do

        !cndmx= gdiff_exm(1,1,pmate)

        ! 3. sacar un dtcri

        if (cndmx > 1.0e-10_rp) then
           dtaux =   xfact * hlmin * hlmin  / cndmx  !conductivity is diffusion already
        else
           dtaux= 1.0e+10_rp
        end if

        if (dtaux <= dtcri) dtcri = dtaux

     end do elements

     ! compare dtcri with dtext_exm

     if (kfl_genal_exm >= 2) then
        call runend('EXM_UPDTSS: USE ONLY EXPLICIT SCHEMES')
     end if

     !     dtcri_exm = dtcri * safet_exm  !!! THIS WAS WRONG
     dtcri_exm = dtcri

     if (kfl_timco == 1) then 
        dtmin = dtcri
     end if

  end if
  !
  ! Look for minimum over whole mesh
  !
  call PAR_MIN(dtmin,'IN MY CODE')

  dtcri_exm = dtmin
  if(dtcri_exm/=0.0_rp) dtinv_exm = 1.0_rp/dtcri_exm
  if(kfl_timco==1) dtinv=max(dtinv,dtinv_exm)
  
  call times(2) % add()

end subroutine exm_updtss
