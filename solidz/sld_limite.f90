!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_limite
  !------------------------------------------------------------------------
  !****f* Solidz/sld_limite
  ! NAME 
  !    sld_limite
  ! DESCRIPTION

  ! USES
  ! USED BY
  !    sld_solite
  !------------------------------------------------------------------------

  use def_master   ! general global variables
  use def_domain   ! geometry information
  use def_solidz   ! general solidz module information

  implicit none
  integer(ip) :: ielem,idime,jdime,igaus,inode,ichek,ipoin,itott   
  integer(ip) :: pelty,pnode,pgaus,plapl,pevat
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elmof(mnode)                           ! Modulating fields 
  real(rp)    :: eldis(ndime,mnode,ncomp_sld)
  real(rp)    :: elvel(ndime,mnode,ncomp_sld)
  real(rp)    :: eldix(ndime,mnode,ncomp_sld)
  real(rp)    :: elvex(ndime,mnode,ncomp_sld)
  real(rp)    :: gpvol(mgaus),gpdet(mgaus)
  real(rp)    :: gpgdi(ndime,ndime,mgaus)              ! Gradient of displacement F
  real(rp)    :: gpigd(ndime,ndime,mgaus)              ! Inverse Deformation gradient
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj



  !
  ! Loop over elements
  !
  elements: do ielem=1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)
     plapl = llapl(pelty)
     pevat = ndofn_sld*pnode

     !
     ! Gather operations
     !
     call sld_elmgat(&
          pnode,lnods(1,ielem),eldis,elvel,elcod,eldix,elvex,elmof)
     !
     ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
     !
     call elmcar(&
          pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
          elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
          gphes,ielem)

     !
     ! GPGDI: Deformation gradients : WRONG IF HYBRID MESH
     !
     do idime=1,ndime
        do igaus=1,pgaus
           do jdime=1,ndime
              gpgdi(idime,jdime,igaus)=0.0_rp
              do inode=1,pnode
                 gpgdi(idime,jdime,igaus)=gpgdi(idime,jdime,igaus)&
                      +eldis(idime,inode,1)*gpcar(jdime,inode,igaus)
              end do
           end do
        end do
     end do

     ichek=0_ip
     do igaus=1,pgaus
        gpgdi(1,1,igaus)= gpgdi(1,1,igaus) + 1.0_rp
        gpgdi(2,2,igaus)= gpgdi(2,2,igaus) + 1.0_rp
        if (ndime==3) gpgdi(3,3,igaus)= gpgdi(3,3,igaus) + 1.0_rp
        !Calcul of J
        call invmtx(gpgdi(1,1,igaus),gpigd(1,1,igaus),gpdet(igaus),ndime)  
        if (gpdet(igaus)<0.0001_rp) then 
           write(*,*) 'SLD_LIMITER:  **WARNING** LIMITER < 0.0001 Check elemennt ', ielem
           ichek=ichek+1_ip
        end if               
     end do

     if (ielem==12041) then 
     write(4111,100) cutim,gpdet(1)
     end if 
     if (ielem==7072) then 
     write(4112,100) cutim,gpdet(1)
     end if  

     !If the J is to small, keep the last displcament calculated
     if (ichek>0_ip) then
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           itott=(ipoin-1)*ndime     
           do idime=1,ndime
               itott=itott+1_ip
               rhsid(itott)= 0.0_rp
           end do
        end do
     end if 



  end do elements

100 format (2(F16.8,','))  
end subroutine sld_limite
