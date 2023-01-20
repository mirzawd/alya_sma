!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_heatfl(gradt)
  !------------------------------------------------------------------------
  !
  ! Case of a scalar
  !
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper
  implicit none
  real(rp),    intent(out), target :: gradt(ndime,npoin)
  integer(ip)                      :: ipoin,idime,inode,ielem
  integer(ip)                      :: pnode,pelty,pmate,dummi,jnode
  real(rp)                         :: gpdet,gpvol,gpcar(ndime,mnode) 
  real(rp)                         :: eltem(mnode,2),elcod(ndime,mnode)
  real(rp)                         :: elvel(ndime,mnode)
  real(rp)                         :: xjaci(9),xjacm(9),gpcon(1)
  real(rp)                         :: gptem(mnode),fact1,fact2
  real(rp)                         :: gpden(2)

  if(kfl_paral/=0) then
     !
     ! Initialization
     !
     gradt = 0.0_rp
     !
     ! Loop over elements
     !
     elements: do ielem = 1,nelem
        pelty = ltype(ielem) 
        pnode = nnode(pelty)
        pmate = 1
        if( nmate > 1 ) pmate = lmate(ielem)
        !
        ! Gather vectors
        !
        call tem_elmgah(&
             pnode,lnods(1,ielem),eltem,elvel,elcod)
        !
        ! Compute GPTEM
        !
        call gather(&
             2_ip,pnode,pnode,1_ip,dummi,&
             elmar(pelty)%shapc,eltem,gptem)
        !
        ! Loop over Gauss points (which are nodes)
        !
        gauss_points: do inode = 1,pnode
           !
           !  Properties: GPCON (k)
           ! 
           call ker_proper('DENSI','IGAUS',1_ip,ielem,gpden,pnode,1_ip,elmar(pelty)%shapc(:,inode))
           call ker_proper('CONDU','IGAUS',1_ip,ielem,gpcon,pnode,1_ip,elmar(pelty)%shapc(:,inode))

           ipoin = lnods(inode,ielem)
           call elmder(&
                pnode,ndime,elmar(pelty) % deric(1,1,inode),&
                elcod,gpcar,gpdet,xjacm,xjaci)
           gpvol = elmar(pelty) % weigc(inode) * gpdet
           if( kfl_naxis == 1 ) gpvol = gpvol * twopi * elcod(1,inode)
           !
           ! Gradient
           !
           fact1 = gpvol * gpcon(1)
           do jnode = 1,pnode
              fact2 = fact1 * eltem(jnode,1)
              do idime = 1,ndime
                 gradt(idime,ipoin) = gradt(idime,ipoin) &
                      + fact2 * gpcar(idime,jnode)
              end do
           end do
        end do gauss_points
     end do elements
     !
     ! Parall service
     !
     call pararr('SLX',NPOIN_TYPE,ndime*npoin,gradt)
     !
     ! Solve diagonal system
     !
     do ipoin = 1,npoin
        do idime = 1,ndime
           gradt(idime,ipoin) = gradt(idime,ipoin) / vmasc(ipoin)
        end do
     end do
  end if

end subroutine tem_heatfl
