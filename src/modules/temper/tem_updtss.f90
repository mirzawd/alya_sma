!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_updtss()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_updtss
  ! NAME 
  !    tem_updtss
  ! DESCRIPTION
  !    This routine computes the time step size for the temperature
  !    equation.
  ! USED BY
  !    tem_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper
  use mod_ADR,            only : ADR_critical_time_step
  use mod_ADR,            only : FROM_CRITICAL
  use mod_communications, only : PAR_MIN
  use mod_tem_turbul
  implicit none 
  integer(ip)      :: ielem,inode
  integer(ip)      :: pnode,pmate,pelty,porde
  real(rp)         :: dtmin,dtcri(2)
  real(rp)         :: chale(2),hleng(3),tragl(9),chave(6)
  real(rp)         :: gpcar(ndime,mnode),gphes(ntens,mnode)
  real(rp)         :: elcod(ndime,mnode),elvel(ndime,mnode)
  real(rp)         :: eltem(mnode,ADR_tem % ntime)
  real(rp)         :: elmsh(ndime,mnode)
  real(rp)         :: gpden(2),gpcon(2),gpdif(2)
  real(rp)         :: gpsph(2),gptur(1),gprea(2)
  real(rp)         :: gpvel(3),gpvol(2),gpgrd(3,1)

  if( ADR_tem % kfl_time_integration /= 0 ) then

#ifdef EVENT
     call mpitrace_user_function(1)
#endif

     dtmin = huge(1.0_rp)

     do ielem = 1,nelem
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           porde = lorde(pelty)
           pmate = lmate(ielem)
           !
           ! Gather operations 
           !
           call tem_elmgat(&
                ielem,pnode,lnods(:,ielem),eltem,elvel,elcod,elmsh)
           !
           ! Shape function and derivatives
           !
           call elmcar(&
                pnode,1_ip,0_ip,elmar(pelty) % weicg,elmar(pelty) % shacg,&
                elmar(pelty) % dercg,elmar(pelty) % hescg,elcod,gpvol,gpcar,&
                gphes,ielem)
           !
           ! HLENG and TRAGL at center of gravity
           !
           call elmlen(ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
                hnatu(pelty),hleng)
           !
           ! Compute the characteristic length CHALE
           !
           call elmchl(tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,& 
                porde,hnatu(pelty),kfl_advec_tem,kfl_ellen_tem)
           !
           ! Properties
           !
           call ker_proper('DENSI','COG  ',1_ip,ielem,gpden)
           call ker_proper('CONDU','COG  ',1_ip,ielem,gpcon)
           call ker_proper('SPHEA','COG  ',1_ip,ielem,gpsph)
           call ker_proper('TURBU','COG  ',1_ip,ielem,gptur)
           ! 
           ! Coupling with turbul
           !
           call tem_turbul(& 
                pnode,1_ip,1_ip,1_ip,gpcon,gpsph,gpdif,gpgrd,gpden,gptur)
           !
           ! Reaction
           !
           call tem_elmrea(1_ip,pnode,1_ip,1_ip,1_ip,elvel,gpden,gpcar,gprea)
           !
           ! Time step at element center of gravity
           !
           gpden(1) = gpden(1) * gpsph(1)
           gpvel    = 0.0_rp

           if( kfl_advec_tem /= 0 ) then
              do inode = 1,pnode
                 gpvel(1:ndime) = gpvel(1:ndime) + elmar(pelty) % shacg(inode) * elvel(1:ndime,inode)
              end do
              if( associated(velom) ) then
                 do inode = 1,pnode
                    gpvel(1:ndime) = gpvel(1:ndime) - elmar(pelty) % shacg(inode) * elmsh(1:ndime,inode)
                 end do
              end if
           end if
           call ADR_critical_time_step(ADR_tem,gpden,gpvel,gpdif,gprea,dtcri,chale(1),chale(2))
           !
           ! Take minimum
           !
           dtmin = min(dtmin,dtcri(1))
        end if
     end do
     !
     ! Look for minimum over whole mesh
     !
     call PAR_MIN(dtmin,'IN MY CODE')

     dtcri_tem = dtmin

     if( dtcri_tem /= 0.0_rp ) dtinv_tem = 1.0_rp / (dtcri_tem * safet_tem)
     if( ADR_tem % kfl_time_step_strategy == FROM_CRITICAL ) dtinv = max(dtinv,dtinv_tem)
     ADR_tem % dtinv = dtinv_tem
     
#ifdef EVENT
     call mpitrace_user_function(0)
#endif

  end if

end subroutine tem_updtss
