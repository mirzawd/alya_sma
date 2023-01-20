!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_post_flame_index()
   !------------------------------------------------------------------------
   ! NAME 
   !    chm_flame_index
   ! DESCRIPTION
   !    Elemental operations for post-processing Flame index
   !    in flamelet combustion model
   ! USES
   ! USED BY
   !    chm_outvar
   !***
   !------------------------------------------------------------------------
   use def_kintyp, only      : ip,rp
   use def_domain, only      : ndime,mnode,mgaus,nelem,ltype,   &
                               nnode,ngaus,llapl,lorde,ltopo,   &
                               lnods,elmar,hnatu,ntens

   use def_chemic, only      : kfl_advec_chm,kfl_ellen_chm,     &
                               ADR_chm,nclas_chm,zgrad_gp,      &
                               flame_index_gp,kfl_izmean_chm,kfl_icmean_chm

   implicit none

   integer(ip) :: ielem,igaus,inode
   integer(ip) :: pelty,pnode
   integer(ip) :: pgaus,plapl,porde,ptopo

   real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
   real(rp)    :: elcod(ndime,mnode)
   real(rp)    :: elvel(ndime,mnode)
   real(rp)    :: gpvol(mgaus)
   real(rp)    :: gpcar(ndime,mnode,mgaus)                 
   real(rp)    :: gphes(ntens,mnode,mgaus)
   real(rp)    :: gp_z(mgaus)
   real(rp)    :: gp_c(mgaus)
 

   real(rp)    :: dummr(mgaus*ndime)
   real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)

   external    :: elmlen
   external    :: elmchl
   external    :: elmcar
   external    :: chm_calc_flame_index
   external    :: chm_post_gather

   elements: do ielem = 1,nelem

      !
      ! Element dimensions
      !
      pelty = ltype(ielem)
      pnode = nnode(pelty)
      pgaus = ngaus(pelty)
      plapl = llapl(pelty) 
      porde = lorde(pelty)
      ptopo = ltopo(pelty)

      !
      ! Gather all
      !
      call chm_post_gather(&
           pnode,lnods(1:pnode,ielem),elcon(1:pnode,:,:),elcod,elvel)
      !
      ! CHALE, HLENG and TRAGL 
      !
      call elmlen(&
           ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
           hleng)
      call elmchl(&
           tragl,hleng,elcod,dummr,chave,chale,pelty,pnode,&
           porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)

      !
      ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
      !
      call elmcar(&
           pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
           elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)

      gp_z(:) = 0.0_rp
      gp_c(:) = 0.0_rp
      do igaus = 1,pgaus
         do inode = 1,pnode
            gp_z(igaus) = gp_z(igaus)+ elmar(pelty)%shape(inode,igaus) * elcon(inode,kfl_izmean_chm,1)
            gp_c(igaus) = gp_c(igaus)+ elmar(pelty)%shape(inode,igaus) * elcon(inode,kfl_icmean_chm,1)
         end do
      end do
      !
      ! Compute scalar dissipation rate at gauss points
      !
      call chm_calc_flame_index(&
           pnode,pgaus,elcon(1:pnode,kfl_izmean_chm,1),elcon(1:pnode,kfl_icmean_chm,1),&
           gpcar,flame_index_gp(ielem)%a(1:pgaus,1),&
           zgrad_gp(ielem) % a (1:pgaus,1) ,gp_z(1:pgaus), gp_c(1:pgaus))

   end do elements
end subroutine chm_post_flame_index
