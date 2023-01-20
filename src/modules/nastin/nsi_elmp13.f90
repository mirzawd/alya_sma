!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmp13(&
     pnode,gpden,gpvis,gppor,gpgvi,&
     gpsp1,gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,&
     gpadv,gpvep,gprhs,rmomu,rcont,p1vec,p2vec,p2sca,&
     wgrgr,wgrvi,elauu,elaup,elapp,elapu,elrbu,elrbp,&
     rmom2,p1ve2,gpst1,gpsgs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmp13
  ! NAME 
  !    nsi_elmp13
  ! DESCRIPTION
  !    Compute element matrix and RHS
  ! USES
  ! USED BY
  !    nsi_elmop3
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_nastin, only       :  kfl_stabi_nsi,dtsgs_nsi
  use def_nastin, only       :  fvela_nsi,corio_nsi
  implicit none
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(in)    :: gpden,gpvis
  real(rp),    intent(in)    :: gppor,gpgvi(ndime)
  real(rp),    intent(inout) :: gpsp1,gptt1
  real(rp),    intent(inout) :: gpsp2,gptt2
  real(rp),    intent(in)    :: gpvol,gpsha(pnode)
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(in)    :: gpadv(ndime)
  real(rp),    intent(inout) :: gpvep(ndime)
  real(rp),    intent(in)    :: gprhs(ndime)
  real(rp),    intent(in)    :: rmomu(pnode)
  real(rp),    intent(out)   :: rcont(ndime,pnode)
  real(rp),    intent(out)   :: p1vec(pnode)
  real(rp),    intent(out)   :: p2vec(ndime,pnode)
  real(rp),    intent(out)   :: p2sca(pnode)
  real(rp),    intent(out)   :: wgrgr(pnode,pnode)
  real(rp),    intent(out)   :: wgrvi(pnode)
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elrbu(ndime,pnode)
  real(rp),    intent(out)   :: elrbp(pnode)
  real(rp),    intent(out)   :: rmom2(ndime,ndime,pnode)
  real(rp),    intent(out)   :: p1ve2(ndime,ndime,pnode)
  real(rp),    intent(in)    :: gpst1
  real(rp),    intent(in)    :: gpsgs(ndime,*)
  integer(ip)                :: inode,jdime,idime,idofv,jdofv,jnode
  real(rp)                   :: fact0,fact1,fact2,fact3,fact4,fact5,fact6
  real(rp)                   :: fact41,fact42,fact43,fact44,fact8
  real(rp)                   :: factx,facty,factz

  !----------------------------------------------------------------------
  !
  ! Test functions
  !
  !----------------------------------------------------------------------

  !
  ! P1VEC =  v * [ (tau1'/tau1) - tau1' * sig ] + tau1' * rho * (uc.grad)v
  ! P2SCA =  ( tau2^{-1} * tau2' ) * q
  ! P2VEC =  tau1' * grad(q)
  ! RCONT =  div(u)
  ! WGRVI =  grad(Ni) . grad(mu)
  ! WGRGR =  grad(Ni) . grad(Nj)
  !
  !gptt1      = 1.0_rp
  !gptt2      = 1.0_rp

  fact1      = gpsp1*gpden         * gpvol
  fact2      = (gptt1-gpsp1*gppor) * gpvol
  fact3      = gpsp1               * gpvol
  fact4      = gptt2               * gpvol

  p1vec(1)   = fact2*gpsha(1) + fact1 * ( gpadv(1)*gpcar(1,1) + gpadv(2)*gpcar(2,1) + gpadv(3)*gpcar(3,1) )
  p1vec(2)   = fact2*gpsha(2) + fact1 * ( gpadv(1)*gpcar(1,2) + gpadv(2)*gpcar(2,2) + gpadv(3)*gpcar(3,2) )
  p1vec(3)   = fact2*gpsha(3) + fact1 * ( gpadv(1)*gpcar(1,3) + gpadv(2)*gpcar(2,3) + gpadv(3)*gpcar(3,3) )
  p1vec(4)   = fact2*gpsha(4) + fact1 * ( gpadv(1)*gpcar(1,4) + gpadv(2)*gpcar(2,4) + gpadv(3)*gpcar(3,4) )

  p2sca(1)   = fact4*gpsha(1)
  p2sca(2)   = fact4*gpsha(2)
  p2sca(3)   = fact4*gpsha(3)
  p2sca(4)   = fact4*gpsha(4)

  p2vec(1,1) = fact3*gpcar(1,1)
  p2vec(2,1) = fact3*gpcar(2,1)
  p2vec(3,1) = fact3*gpcar(3,1)
  p2vec(1,2) = fact3*gpcar(1,2)
  p2vec(2,2) = fact3*gpcar(2,2)
  p2vec(3,2) = fact3*gpcar(3,2)
  p2vec(1,3) = fact3*gpcar(1,3)
  p2vec(2,3) = fact3*gpcar(2,3)
  p2vec(3,3) = fact3*gpcar(3,3)
  p2vec(1,4) = fact3*gpcar(1,4)
  p2vec(2,4) = fact3*gpcar(2,4)
  p2vec(3,4) = fact3*gpcar(3,4)

  rcont(1,1) = gpcar(1,1)
  rcont(2,1) = gpcar(2,1)
  rcont(3,1) = gpcar(3,1)
  rcont(1,2) = gpcar(1,2)
  rcont(2,2) = gpcar(2,2)
  rcont(3,2) = gpcar(3,2)
  rcont(1,3) = gpcar(1,3)
  rcont(2,3) = gpcar(2,3)
  rcont(3,3) = gpcar(3,3)
  rcont(1,4) = gpcar(1,4)
  rcont(2,4) = gpcar(2,4)
  rcont(3,4) = gpcar(3,4)

  wgrvi(1)   = gpgvi(1)*gpcar(1,1)   + gpgvi(2)*gpcar(2,1)   + gpgvi(3)*gpcar(3,1)
  wgrvi(2)   = gpgvi(1)*gpcar(1,2)   + gpgvi(2)*gpcar(2,2)   + gpgvi(3)*gpcar(3,2)
  wgrvi(3)   = gpgvi(1)*gpcar(1,3)   + gpgvi(2)*gpcar(2,3)   + gpgvi(3)*gpcar(3,3)
  wgrvi(4)   = gpgvi(1)*gpcar(1,4)   + gpgvi(2)*gpcar(2,4)   + gpgvi(3)*gpcar(3,4)

  wgrgr(1,1) = gpcar(1,1)*gpcar(1,1) + gpcar(2,1)*gpcar(2,1) + gpcar(3,1)*gpcar(3,1)   
  wgrgr(1,2) = gpcar(1,1)*gpcar(1,2) + gpcar(2,1)*gpcar(2,2) + gpcar(3,1)*gpcar(3,2)   
  wgrgr(1,3) = gpcar(1,1)*gpcar(1,3) + gpcar(2,1)*gpcar(2,3) + gpcar(3,1)*gpcar(3,3)   
  wgrgr(1,4) = gpcar(1,1)*gpcar(1,4) + gpcar(2,1)*gpcar(2,4) + gpcar(3,1)*gpcar(3,4)   

  wgrgr(2,1) = wgrgr(1,2) 
  wgrgr(2,2) = gpcar(1,2)*gpcar(1,2) + gpcar(2,2)*gpcar(2,2) + gpcar(3,2)*gpcar(3,2)   
  wgrgr(2,3) = gpcar(1,2)*gpcar(1,3) + gpcar(2,2)*gpcar(2,3) + gpcar(3,2)*gpcar(3,3)   
  wgrgr(2,4) = gpcar(1,2)*gpcar(1,4) + gpcar(2,2)*gpcar(2,4) + gpcar(3,2)*gpcar(3,4)   

  wgrgr(3,1) = wgrgr(1,3)
  wgrgr(3,2) = wgrgr(2,3) 
  wgrgr(3,3) = gpcar(1,3)*gpcar(1,3) + gpcar(2,3)*gpcar(2,3) + gpcar(3,3)*gpcar(3,3)   
  wgrgr(3,4) = gpcar(1,3)*gpcar(1,4) + gpcar(2,3)*gpcar(2,4) + gpcar(3,3)*gpcar(3,4)   

  wgrgr(4,1) = wgrgr(1,4) 
  wgrgr(4,2) = wgrgr(2,4) 
  wgrgr(4,3) = wgrgr(3,4)
  wgrgr(4,4) = gpcar(1,4)*gpcar(1,4) + gpcar(2,4)*gpcar(2,4) + gpcar(3,4)*gpcar(3,4)   

  !----------------------------------------------------------------------
  !
  ! Auu
  !
  !----------------------------------------------------------------------

  !
  ! div(u)*tau2'*div(v) + ( rmomu(u) , p1vec(v) ) 
  !       + ( mu*dui/dxk , dv/dxk ) + ( mu*d^2ui/dxk^2, vj )
  !       + ( grad(mut).grad(u), v )
  !
  fact0    = gpsp2 * gpvol
  fact6    = gpvis * gpvol

  ! inode = 1 -------------------------------------------------

  fact1  = fact0 * gpcar(1,1) 
  fact2  = fact0 * gpcar(2,1)
  fact3  = fact0 * gpcar(3,1)
  fact5  = gpsha(1) * gpvol

  fact41 = p1vec(1) * rmomu(1) + fact5 * wgrvi(1) + fact6 * wgrgr(1,1)
  fact42 = p1vec(1) * rmomu(2) + fact5 * wgrvi(2) + fact6 * wgrgr(1,2)   
  fact43 = p1vec(1) * rmomu(3) + fact5 * wgrvi(3) + fact6 * wgrgr(1,3)  
  fact44 = p1vec(1) * rmomu(4) + fact5 * wgrvi(4) + fact6 * wgrgr(1,4)
  ! jnode = 1
  elauu(1,1)  = fact1 * rcont(1,1)+fact41
  elauu(2,1)  = fact2 * rcont(1,1)
  elauu(3,1)  = fact3 * rcont(1,1)
  elauu(1,2)  = fact1 * rcont(2,1)
  elauu(2,2)  = fact2 * rcont(2,1)+fact41
  elauu(3,2)  = fact3 * rcont(2,1)
  elauu(1,3)  = fact1 * rcont(3,1)
  elauu(2,3)  = fact2 * rcont(3,1)
  elauu(3,3)  = fact3 * rcont(3,1)+fact41
  ! jnode = 2
  elauu(1,4)  = fact1 * rcont(1,2)+fact42
  elauu(2,4)  = fact2 * rcont(1,2)
  elauu(3,4)  = fact3 * rcont(1,2)
  elauu(1,5)  = fact1 * rcont(2,2)
  elauu(2,5)  = fact2 * rcont(2,2)+fact42
  elauu(3,5)  = fact3 * rcont(2,2)
  elauu(1,6)  = fact1 * rcont(3,2)
  elauu(2,6)  = fact2 * rcont(3,2)
  elauu(3,6)  = fact3 * rcont(3,2)+fact42
  ! jnode = 3  
  elauu(1,7)  = fact1 * rcont(1,3)+fact43
  elauu(2,7)  = fact2 * rcont(1,3)
  elauu(3,7)  = fact3 * rcont(1,3)
  elauu(1,8)  = fact1 * rcont(2,3)
  elauu(2,8)  = fact2 * rcont(2,3)+fact43
  elauu(3,8)  = fact3 * rcont(2,3)
  elauu(1,9)  = fact1 * rcont(3,3)
  elauu(2,9)  = fact2 * rcont(3,3)
  elauu(3,9)  = fact3 * rcont(3,3)+fact43
  ! jnode = 4  
  elauu(1,10) = fact1 * rcont(1,4)+fact44
  elauu(2,10) = fact2 * rcont(1,4)
  elauu(3,10) = fact3 * rcont(1,4)
  elauu(1,11) = fact1 * rcont(2,4)
  elauu(2,11) = fact2 * rcont(2,4)+fact44
  elauu(3,11) = fact3 * rcont(2,4)
  elauu(1,12) = fact1 * rcont(3,4)
  elauu(2,12) = fact2 * rcont(3,4)
  elauu(3,12) = fact3 * rcont(3,4)+fact44

  ! inode = 2 -------------------------------------------------

  fact1  = fact0*gpcar(1,2) 
  fact2  = fact0*gpcar(2,2)
  fact3  = fact0*gpcar(3,2)
  fact5  = gpsha(2)*gpvol

  fact41 = p1vec(2) * rmomu(1) + fact5 * wgrvi(1) + fact6 * wgrgr(2,1)
  fact42 = p1vec(2) * rmomu(2) + fact5 * wgrvi(2) + fact6 * wgrgr(2,2)   
  fact43 = p1vec(2) * rmomu(3) + fact5 * wgrvi(3) + fact6 * wgrgr(2,3)  
  fact44 = p1vec(2) * rmomu(4) + fact5 * wgrvi(4) + fact6 * wgrgr(2,4)
  ! jnode = 1
  elauu(4,1)  = fact1 * rcont(1,1)+fact41
  elauu(5,1)  = fact2 * rcont(1,1)
  elauu(6,1)  = fact3 * rcont(1,1)
  elauu(4,2)  = fact1 * rcont(2,1)
  elauu(5,2)  = fact2 * rcont(2,1)+fact41
  elauu(6,2)  = fact3 * rcont(2,1)
  elauu(4,3)  = fact1 * rcont(3,1)
  elauu(5,3)  = fact2 * rcont(3,1)
  elauu(6,3)  = fact3 * rcont(3,1)+fact41
  ! jnode = 2
  elauu(4,4)  = fact1 * rcont(1,2)+fact42
  elauu(5,4)  = fact2 * rcont(1,2)
  elauu(6,4)  = fact3 * rcont(1,2)
  elauu(4,5)  = fact1 * rcont(2,2)
  elauu(5,5)  = fact2 * rcont(2,2)+fact42
  elauu(6,5)  = fact3 * rcont(2,2)
  elauu(4,6)  = fact1 * rcont(3,2)
  elauu(5,6)  = fact2 * rcont(3,2)
  elauu(6,6)  = fact3 * rcont(3,2)+fact42
  ! jnode = 3  
  elauu(4,7)  = fact1 * rcont(1,3)+fact43
  elauu(5,7)  = fact2 * rcont(1,3)
  elauu(6,7)  = fact3 * rcont(1,3)
  elauu(4,8)  = fact1 * rcont(2,3)
  elauu(5,8)  = fact2 * rcont(2,3)+fact43
  elauu(6,8)  = fact3 * rcont(2,3)
  elauu(4,9)  = fact1 * rcont(3,3)
  elauu(5,9)  = fact2 * rcont(3,3)
  elauu(6,9)  = fact3 * rcont(3,3)+fact43
  ! jnode = 4  
  elauu(4,10) = fact1 * rcont(1,4)+fact44
  elauu(5,10) = fact2 * rcont(1,4)
  elauu(6,10) = fact3 * rcont(1,4)
  elauu(4,11) = fact1 * rcont(2,4)
  elauu(5,11) = fact2 * rcont(2,4)+fact44
  elauu(6,11) = fact3 * rcont(2,4)
  elauu(4,12) = fact1 * rcont(3,4)
  elauu(5,12) = fact2 * rcont(3,4)
  elauu(6,12) = fact3 * rcont(3,4)+fact44

  ! inode = 3 -------------------------------------------------

  fact1  = fact0*gpcar(1,3) 
  fact2  = fact0*gpcar(2,3)
  fact3  = fact0*gpcar(3,3)
  fact5  = gpsha(3)*gpvol

  fact41 = p1vec(3) * rmomu(1) + fact5*wgrvi(1) + fact6*wgrgr(3,1)
  fact42 = p1vec(3) * rmomu(2) + fact5*wgrvi(2) + fact6*wgrgr(3,2)   
  fact43 = p1vec(3) * rmomu(3) + fact5*wgrvi(3) + fact6*wgrgr(3,3)  
  fact44 = p1vec(3) * rmomu(4) + fact5*wgrvi(4) + fact6*wgrgr(3,4)
  ! jnode = 1
  elauu(7,1)  = fact1 * rcont(1,1)+fact41
  elauu(8,1)  = fact2 * rcont(1,1)
  elauu(9,1)  = fact3 * rcont(1,1)
  elauu(7,2)  = fact1 * rcont(2,1)
  elauu(8,2)  = fact2 * rcont(2,1)+fact41
  elauu(9,2)  = fact3 * rcont(2,1)
  elauu(7,3)  = fact1 * rcont(3,1)
  elauu(8,3)  = fact2 * rcont(3,1)
  elauu(9,3)  = fact3 * rcont(3,1)+fact41
  ! jnode = 2
  elauu(7,4)  = fact1 * rcont(1,2)+fact42
  elauu(8,4)  = fact2 * rcont(1,2)
  elauu(9,4)  = fact3 * rcont(1,2)
  elauu(7,5)  = fact1 * rcont(2,2)
  elauu(8,5)  = fact2 * rcont(2,2)+fact42
  elauu(9,5)  = fact3 * rcont(2,2)
  elauu(7,6)  = fact1 * rcont(3,2)
  elauu(8,6)  = fact2 * rcont(3,2)
  elauu(9,6)  = fact3 * rcont(3,2)+fact42
  ! jnode = 3  
  elauu(7,7)  = fact1 * rcont(1,3)+fact43
  elauu(8,7)  = fact2 * rcont(1,3)
  elauu(9,7)  = fact3 * rcont(1,3)
  elauu(7,8)  = fact1 * rcont(2,3)
  elauu(8,8)  = fact2 * rcont(2,3)+fact43
  elauu(9,8)  = fact3 * rcont(2,3)
  elauu(7,9)  = fact1 * rcont(3,3)
  elauu(8,9)  = fact2 * rcont(3,3)
  elauu(9,9)  = fact3 * rcont(3,3)+fact43
  ! jnode = 4  
  elauu(7,10) = fact1 * rcont(1,4)+fact44
  elauu(8,10) = fact2 * rcont(1,4)
  elauu(9,10) = fact3 * rcont(1,4)
  elauu(7,11) = fact1 * rcont(2,4)
  elauu(8,11) = fact2 * rcont(2,4)+fact44
  elauu(9,11) = fact3 * rcont(2,4)
  elauu(7,12) = fact1 * rcont(3,4)
  elauu(8,12) = fact2 * rcont(3,4)
  elauu(9,12) = fact3 * rcont(3,4)+fact44

  ! inode = 4 -------------------------------------------------

  fact1  = fact0*gpcar(1,4) 
  fact2  = fact0*gpcar(2,4)
  fact3  = fact0*gpcar(3,4)
  fact5  = gpsha(4)*gpvol

  fact41 = p1vec(4) * rmomu(1) + fact5*wgrvi(1) + fact6*wgrgr(4,1)
  fact42 = p1vec(4) * rmomu(2) + fact5*wgrvi(2) + fact6*wgrgr(4,2)   
  fact43 = p1vec(4) * rmomu(3) + fact5*wgrvi(3) + fact6*wgrgr(4,3)  
  fact44 = p1vec(4) * rmomu(4) + fact5*wgrvi(4) + fact6*wgrgr(4,4)
  ! jnode = 1
  elauu(10,1)  = fact1 * rcont(1,1)+fact41
  elauu(11,1)  = fact2 * rcont(1,1)
  elauu(12,1)  = fact3 * rcont(1,1)
  elauu(10,2)  = fact1 * rcont(2,1)
  elauu(11,2)  = fact2 * rcont(2,1)+fact41
  elauu(12,2)  = fact3 * rcont(2,1)
  elauu(10,3)  = fact1 * rcont(3,1)
  elauu(11,3)  = fact2 * rcont(3,1)
  elauu(12,3)  = fact3 * rcont(3,1)+fact41
  ! jnode = 2
  elauu(10,4)  = fact1 * rcont(1,2)+fact42
  elauu(11,4)  = fact2 * rcont(1,2)
  elauu(12,4)  = fact3 * rcont(1,2)
  elauu(10,5)  = fact1 * rcont(2,2)
  elauu(11,5)  = fact2 * rcont(2,2)+fact42
  elauu(12,5)  = fact3 * rcont(2,2)
  elauu(10,6)  = fact1 * rcont(3,2)
  elauu(11,6)  = fact2 * rcont(3,2)
  elauu(12,6)  = fact3 * rcont(3,2)+fact42
  ! jnode = 3  
  elauu(10,7)  = fact1 * rcont(1,3)+fact43
  elauu(11,7)  = fact2 * rcont(1,3)
  elauu(12,7)  = fact3 * rcont(1,3)
  elauu(10,8)  = fact1 * rcont(2,3)
  elauu(11,8)  = fact2 * rcont(2,3)+fact43
  elauu(12,8)  = fact3 * rcont(2,3)
  elauu(10,9)  = fact1 * rcont(3,3)
  elauu(11,9)  = fact2 * rcont(3,3)
  elauu(12,9)  = fact3 * rcont(3,3)+fact43
  ! jnode = 4  
  elauu(10,10) = fact1 * rcont(1,4)+fact44
  elauu(11,10) = fact2 * rcont(1,4)
  elauu(12,10) = fact3 * rcont(1,4)
  elauu(10,11) = fact1 * rcont(2,4)
  elauu(11,11) = fact2 * rcont(2,4)+fact44
  elauu(12,11) = fact3 * rcont(2,4)
  elauu(10,12) = fact1 * rcont(3,4)
  elauu(11,12) = fact2 * rcont(3,4)
  elauu(12,12) = fact3 * rcont(3,4)+fact44

  !----------------------------------------------------------------------
  !
  ! Aup
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: - ( p , div(v) ) + ( grad(p) , p1vec(v)-v )
  !  
  ! jnode = 1 -------------------------------------------------

  fact2      =   gpvol*gpsha(1)
  ! inode = 1
  fact1      = - gpvol*gpsha(1)+p1vec(1)
  elaup(1,1) = - fact2*gpcar(1,1) + fact1*gpcar(1,1)
  elaup(2,1) = - fact2*gpcar(2,1) + fact1*gpcar(2,1)
  elaup(3,1) = - fact2*gpcar(3,1) + fact1*gpcar(3,1)
  ! inode = 2  
  fact1      = - gpvol*gpsha(2)+p1vec(2)
  elaup(4,1) = - fact2*gpcar(1,2) + fact1*gpcar(1,1)
  elaup(5,1) = - fact2*gpcar(2,2) + fact1*gpcar(2,1)
  elaup(6,1) = - fact2*gpcar(3,2) + fact1*gpcar(3,1)
  ! inode = 3  
  fact1      = - gpvol*gpsha(3)+p1vec(3)
  elaup(7,1) = - fact2*gpcar(1,3) + fact1*gpcar(1,1)
  elaup(8,1) = - fact2*gpcar(2,3) + fact1*gpcar(2,1)
  elaup(9,1) = - fact2*gpcar(3,3) + fact1*gpcar(3,1)
  ! inode = 4  
  fact1       = - gpvol*gpsha(4)+p1vec(4)
  elaup(10,1) = - fact2*gpcar(1,4) + fact1*gpcar(1,1)
  elaup(11,1) = - fact2*gpcar(2,4) + fact1*gpcar(2,1)
  elaup(12,1) = - fact2*gpcar(3,4) + fact1*gpcar(3,1)

  ! jnode = 2 -------------------------------------------------

  fact2      =   gpvol*gpsha(2)
  ! inode = 1
  fact1      = - gpvol*gpsha(1)+p1vec(1)
  elaup(1,2) = - fact2*gpcar(1,1) + fact1*gpcar(1,2)
  elaup(2,2) = - fact2*gpcar(2,1) + fact1*gpcar(2,2)
  elaup(3,2) = - fact2*gpcar(3,1) + fact1*gpcar(3,2)
  ! inode = 2  
  fact1      = - gpvol*gpsha(2)+p1vec(2)
  elaup(4,2) = - fact2*gpcar(1,2) + fact1*gpcar(1,2)
  elaup(5,2) = - fact2*gpcar(2,2) + fact1*gpcar(2,2)
  elaup(6,2) = - fact2*gpcar(3,2) + fact1*gpcar(3,2)
  ! inode = 3  
  fact1      = - gpvol*gpsha(3)+p1vec(3)
  elaup(7,2) = - fact2*gpcar(1,3) + fact1*gpcar(1,2)
  elaup(8,2) = - fact2*gpcar(2,3) + fact1*gpcar(2,2)
  elaup(9,2) = - fact2*gpcar(3,3) + fact1*gpcar(3,2)
  ! inode = 4  
  fact1       = - gpvol*gpsha(4)+p1vec(4)
  elaup(10,2) = - fact2*gpcar(1,4) + fact1*gpcar(1,2)
  elaup(11,2) = - fact2*gpcar(2,4) + fact1*gpcar(2,2)
  elaup(12,2) = - fact2*gpcar(3,4) + fact1*gpcar(3,2)

  ! jnode = 3 -------------------------------------------------

  fact2      =   gpvol*gpsha(3)
  ! inode = 1
  fact1      = - gpvol*gpsha(1)+p1vec(1)
  elaup(1,3) = - fact2*gpcar(1,1) + fact1*gpcar(1,3)
  elaup(2,3) = - fact2*gpcar(2,1) + fact1*gpcar(2,3)
  elaup(3,3) = - fact2*gpcar(3,1) + fact1*gpcar(3,3)
  ! inode = 2  
  fact1          = - gpvol*gpsha(2)+p1vec(2)
  elaup(4,3) = - fact2*gpcar(1,2) + fact1*gpcar(1,3)
  elaup(5,3) = - fact2*gpcar(2,2) + fact1*gpcar(2,3)
  elaup(6,3) = - fact2*gpcar(3,2) + fact1*gpcar(3,3)
  ! inode = 3  
  fact1          = - gpvol*gpsha(3)+p1vec(3)
  elaup(7,3) = - fact2*gpcar(1,3) + fact1*gpcar(1,3)
  elaup(8,3) = - fact2*gpcar(2,3) + fact1*gpcar(2,3)
  elaup(9,3) = - fact2*gpcar(3,3) + fact1*gpcar(3,3)
  ! inode = 4  
  fact1           = - gpvol*gpsha(4)+p1vec(4)
  elaup(10,3) = - fact2*gpcar(1,4) + fact1*gpcar(1,3)
  elaup(11,3) = - fact2*gpcar(2,4) + fact1*gpcar(2,3)
  elaup(12,3) = - fact2*gpcar(3,4) + fact1*gpcar(3,3)

  ! jnode = 4 -------------------------------------------------

  fact2          =   gpvol*gpsha(4)
  ! inode = 1
  fact1          = - gpvol*gpsha(1)+p1vec(1)
  elaup(1,4) = - fact2*gpcar(1,1) + fact1*gpcar(1,4)
  elaup(2,4) = - fact2*gpcar(2,1) + fact1*gpcar(2,4)
  elaup(3,4) = - fact2*gpcar(3,1) + fact1*gpcar(3,4)
  ! inode = 2  
  fact1          = - gpvol*gpsha(2)+p1vec(2)
  elaup(4,4) = - fact2*gpcar(1,2) + fact1*gpcar(1,4)
  elaup(5,4) = - fact2*gpcar(2,2) + fact1*gpcar(2,4)
  elaup(6,4) = - fact2*gpcar(3,2) + fact1*gpcar(3,4)
  ! inode = 3  
  fact1          = - gpvol*gpsha(3)+p1vec(3)
  elaup(7,4) = - fact2*gpcar(1,3) + fact1*gpcar(1,4)
  elaup(8,4) = - fact2*gpcar(2,3) + fact1*gpcar(2,4)
  elaup(9,4) = - fact2*gpcar(3,3) + fact1*gpcar(3,4)
  ! inode = 4  
  fact1           = - gpvol*gpsha(4)+p1vec(4)
  elaup(10,4) = - fact2*gpcar(1,4) + fact1*gpcar(1,4)
  elaup(11,4) = - fact2*gpcar(2,4) + fact1*gpcar(2,4)
  elaup(12,4) = - fact2*gpcar(3,4) + fact1*gpcar(3,4)

  !----------------------------------------------------------------------
  !
  ! Apu
  !
  !----------------------------------------------------------------------

  !
  ! ( div(u) , (tau2^{-1}*tau2')*q )
  ! + ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )
  !
  ! jnode = 1
  elapu(1,1) = rcont(1,1) * p2sca(1) + rmomu(1)   * p2vec(1,1)
  elapu(1,2) = rcont(2,1) * p2sca(1) + rmomu(1)   * p2vec(2,1)
  elapu(1,3) = rcont(3,1) * p2sca(1) + rmomu(1)   * p2vec(3,1)
  ! jnode = 2
  elapu(2,1) = rcont(1,1) * p2sca(2) + rmomu(1)   * p2vec(1,2)
  elapu(2,2) = rcont(2,1) * p2sca(2) + rmomu(1)   * p2vec(2,2)
  elapu(2,3) = rcont(3,1) * p2sca(2) + rmomu(1)   * p2vec(3,2)
  ! jnode = 3
  elapu(3,1) = rcont(1,1) * p2sca(3) + rmomu(1)   * p2vec(1,3)
  elapu(3,2) = rcont(2,1) * p2sca(3) + rmomu(1)   * p2vec(2,3)
  elapu(3,3) = rcont(3,1) * p2sca(3) + rmomu(1)   * p2vec(3,3)
  ! jnode = 4
  elapu(4,1) = rcont(1,1) * p2sca(4) + rmomu(1)   * p2vec(1,4)
  elapu(4,2) = rcont(2,1) * p2sca(4) + rmomu(1)   * p2vec(2,4)
  elapu(4,3) = rcont(3,1) * p2sca(4) + rmomu(1)   * p2vec(3,4)

  ! jnode = 1
  elapu(1,4) = rcont(1,2) * p2sca(1) + rmomu(2)   * p2vec(1,1)
  elapu(1,5) = rcont(2,2) * p2sca(1) + rmomu(2)   * p2vec(2,1)
  elapu(1,6) = rcont(3,2) * p2sca(1) + rmomu(2)   * p2vec(3,1)
  ! jnode = 2
  elapu(2,4) = rcont(1,2) * p2sca(2) + rmomu(2)   * p2vec(1,2)
  elapu(2,5) = rcont(2,2) * p2sca(2) + rmomu(2)   * p2vec(2,2)
  elapu(2,6) = rcont(3,2) * p2sca(2) + rmomu(2)   * p2vec(3,2)
  ! jnode = 3
  elapu(3,4) = rcont(1,2) * p2sca(3) + rmomu(2)   * p2vec(1,3)
  elapu(3,5) = rcont(2,2) * p2sca(3) + rmomu(2)   * p2vec(2,3)
  elapu(3,6) = rcont(3,2) * p2sca(3) + rmomu(2)   * p2vec(3,3)
  ! jnode = 4
  elapu(4,4) = rcont(1,2) * p2sca(4) + rmomu(2)   * p2vec(1,4)
  elapu(4,5) = rcont(2,2) * p2sca(4) + rmomu(2)   * p2vec(2,4)
  elapu(4,6) = rcont(3,2) * p2sca(4) + rmomu(2)   * p2vec(3,4)

  ! jnode = 1
  elapu(1,7) = rcont(1,3) * p2sca(1) + rmomu(3)   * p2vec(1,1)
  elapu(1,8) = rcont(2,3) * p2sca(1) + rmomu(3)   * p2vec(2,1)
  elapu(1,9) = rcont(3,3) * p2sca(1) + rmomu(3)   * p2vec(3,1)
  ! jnode = 2
  elapu(2,7) = rcont(1,3) * p2sca(2) + rmomu(3)   * p2vec(1,2)
  elapu(2,8) = rcont(2,3) * p2sca(2) + rmomu(3)   * p2vec(2,2)
  elapu(2,9) = rcont(3,3) * p2sca(2) + rmomu(3)   * p2vec(3,2)
  ! jnode = 3
  elapu(3,7) = rcont(1,3) * p2sca(3) + rmomu(3)   * p2vec(1,3)
  elapu(3,8) = rcont(2,3) * p2sca(3) + rmomu(3)   * p2vec(2,3)
  elapu(3,9) = rcont(3,3) * p2sca(3) + rmomu(3)   * p2vec(3,3)
  ! jnode = 4
  elapu(4,7) = rcont(1,3) * p2sca(4) + rmomu(3)   * p2vec(1,4)
  elapu(4,8) = rcont(2,3) * p2sca(4) + rmomu(3)   * p2vec(2,4)
  elapu(4,9) = rcont(3,3) * p2sca(4) + rmomu(3)   * p2vec(3,4)

  ! jnode = 1
  elapu(1,10) = rcont(1,4) * p2sca(1) + rmomu(4)   * p2vec(1,1)
  elapu(1,11) = rcont(2,4) * p2sca(1) + rmomu(4)   * p2vec(2,1)
  elapu(1,12) = rcont(3,4) * p2sca(1) + rmomu(4)   * p2vec(3,1)
  ! jnode = 2
  elapu(2,10) = rcont(1,4) * p2sca(2) + rmomu(4)   * p2vec(1,2)
  elapu(2,11) = rcont(2,4) * p2sca(2) + rmomu(4)   * p2vec(2,2)
  elapu(2,12) = rcont(3,4) * p2sca(2) + rmomu(4)   * p2vec(3,2)
  ! jnode = 3
  elapu(3,10) = rcont(1,4) * p2sca(3) + rmomu(4)   * p2vec(1,3)
  elapu(3,11) = rcont(2,4) * p2sca(3) + rmomu(4)   * p2vec(2,3)
  elapu(3,12) = rcont(3,4) * p2sca(3) + rmomu(4)   * p2vec(3,3)
  ! jnode = 4
  elapu(4,10) = rcont(1,4) * p2sca(4) + rmomu(4)   * p2vec(1,4)
  elapu(4,11) = rcont(2,4) * p2sca(4) + rmomu(4)   * p2vec(2,4)
  elapu(4,12) = rcont(3,4) * p2sca(4) + rmomu(4)   * p2vec(3,4)

  !----------------------------------------------------------------------
  !
  ! App
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: ( grad(p) , tau1' grad(q) )
  ! 
  fact1      = p2vec(1,2) * gpcar(1,1) + p2vec(2,2) * gpcar(2,1) + p2vec(3,2) * gpcar(3,1)
  elapp(2,1) = fact1
  elapp(1,2) = fact1
  fact1      = p2vec(1,3) * gpcar(1,1) + p2vec(2,3) * gpcar(2,1) + p2vec(3,3) * gpcar(3,1)
  elapp(3,1) = fact1
  elapp(1,3) = fact1
  fact1      = p2vec(1,4) * gpcar(1,1) + p2vec(2,4) * gpcar(2,1)  + p2vec(3,4) * gpcar(3,1)
  elapp(4,1) = fact1
  elapp(1,4) = fact1
  fact1      = p2vec(1,1) * gpcar(1,1)  + p2vec(2,1) * gpcar(2,1) + p2vec(3,1) * gpcar(3,1)
  elapp(1,1) = fact1
  fact1      = p2vec(1,3) * gpcar(1,2) + p2vec(2,3) * gpcar(2,2) + p2vec(3,3) * gpcar(3,2)
  elapp(3,2) = fact1
  elapp(2,3) = fact1
  fact1      = p2vec(1,4) * gpcar(1,2) + p2vec(2,4) * gpcar(2,2) + p2vec(3,4) * gpcar(3,2)
  elapp(4,2) = fact1
  elapp(2,4) = fact1
  fact1      = p2vec(1,2) * gpcar(1,2) + p2vec(2,2) * gpcar(2,2) + p2vec(3,2) * gpcar(3,2)
  elapp(2,2) = fact1
  fact1      = p2vec(1,4) * gpcar(1,3) + p2vec(2,4) * gpcar(2,3) + p2vec(3,4) * gpcar(3,3)
  elapp(4,3) = fact1
  elapp(3,4) = fact1
  fact1      = p2vec(1,3) * gpcar(1,3) + p2vec(2,3) * gpcar(2,3) + p2vec(3,3) * gpcar(3,3)
  elapp(3,3) = fact1
  fact1      = p2vec(1,4) * gpcar(1,4) + p2vec(2,4) * gpcar(2,4) + p2vec(3,4) * gpcar(3,4)
  elapp(4,4) = fact1

  !----------------------------------------------------------------------
  !
  ! bu and bp
  !
  !----------------------------------------------------------------------
  !
  ! bu = ( f , p1vec(v) ) 
  ! bp = ( rho*f , tau1' grad(q) )
  !

  elrbu(1,1) = p1vec(1)   * gprhs(1)
  elrbu(2,1) = p1vec(1)   * gprhs(2)
  elrbu(3,1) = p1vec(1)   * gprhs(3)
  elrbp(1)   = p2vec(1,1) * gprhs(1) + p2vec(2,1) * gprhs(2) + p2vec(3,1) * gprhs(3)
  elrbu(1,2) = p1vec(2)   * gprhs(1)
  elrbu(2,2) = p1vec(2)   * gprhs(2)
  elrbu(3,2) = p1vec(2)   * gprhs(3)
  elrbp(2)   = p2vec(1,2) * gprhs(1) + p2vec(2,2) * gprhs(2) + p2vec(3,2) * gprhs(3)
  elrbu(1,3) = p1vec(3)   * gprhs(1)
  elrbu(2,3) = p1vec(3)   * gprhs(2)
  elrbu(3,3) = p1vec(3)   * gprhs(3)
  elrbp(3)   = p2vec(1,3) * gprhs(1) + p2vec(2,3) * gprhs(2) + p2vec(3,3) * gprhs(3)
  elrbu(1,4) = p1vec(4)   * gprhs(1)
  elrbu(2,4) = p1vec(4)   * gprhs(2)
  elrbu(3,4) = p1vec(4)   * gprhs(3)
  elrbp(4)   = p2vec(1,4) * gprhs(1) + p2vec(2,4) * gprhs(2) + p2vec(3,4) * gprhs(3)
  !
  ! Projection: ( L*(v)/tau1' , P )
  !
  if( kfl_stabi_nsi == 1 ) then

     fact1  = gpden * dtsgs_nsi
     fact2  = gpvol / gpst1
     factx  = fact1 * gpsgs(1,2) + gpvep(1) / gpst1 ! rho/dt * u'^n + tau^-1 P
     facty  = fact1 * gpsgs(2,2) + gpvep(2) / gpst1 
     factz  = fact1 * gpsgs(3,2) + gpvep(3) / gpst1 
     fact41 = fact2 * gpvep(1)                                   ! tau^-1 P
     fact42 = fact2 * gpvep(2)
     fact43 = fact2 * gpvep(3)

     elrbu(1,1) = elrbu(1,1) + p1vec(1)   * factx - gpsha(1) * fact41
     elrbu(2,1) = elrbu(2,1) + p1vec(1)   * facty - gpsha(1) * fact42
     elrbu(3,1) = elrbu(3,1) + p1vec(1)   * factz - gpsha(1) * fact43
     elrbp(1)   = elrbp(1)   + p2vec(1,1) * factx + p2vec(2,1) * facty + p2vec(3,1) * factz 
     
     elrbu(1,2) = elrbu(1,2) + p1vec(2)   * factx - gpsha(2) * fact41
     elrbu(2,2) = elrbu(2,2) + p1vec(2)   * facty - gpsha(2) * fact42
     elrbu(3,2) = elrbu(3,2) + p1vec(2)   * factz - gpsha(2) * fact43
     elrbp(2)   = elrbp(2)   + p2vec(1,2) * factx + p2vec(2,2) * facty + p2vec(3,2) * factz 
     
     elrbu(1,3) = elrbu(1,3) + p1vec(3)   * factx - gpsha(3) * fact41
     elrbu(2,3) = elrbu(2,3) + p1vec(3)   * facty - gpsha(3) * fact42
     elrbu(3,3) = elrbu(3,3) + p1vec(3)   * factz - gpsha(3) * fact43
     elrbp(3)   = elrbp(3)   + p2vec(1,3) * factx + p2vec(2,3) * facty + p2vec(3,3) * factz 
     
     elrbu(1,4) = elrbu(1,4) + p1vec(4)   * factx - gpsha(4) * fact41
     elrbu(2,4) = elrbu(2,4) + p1vec(4)   * facty - gpsha(4) * fact42
     elrbu(3,4) = elrbu(3,4) + p1vec(4)   * factz - gpsha(4) * fact43
     elrbp(4)   = elrbp(4)   + p2vec(1,4) * factx + p2vec(2,4) * facty + p2vec(3,4) * factz 
     
   
  end if
  !
  ! Coriolis
  !   
  if( corio_nsi > 1.0e-10_rp ) then
     fact0 = 2.0_rp * gpden 
     fact8 = 2.0_rp * gpden * gpsp1
     fact1 = fact0  * fvela_nsi(1)
     fact2 = fact0  * fvela_nsi(2)
     fact3 = fact0  * fvela_nsi(3)
     factx = fact8  * fvela_nsi(1)
     facty = fact8  * fvela_nsi(2)
     factz = fact8  * fvela_nsi(3)
     do inode = 1,pnode
        p1ve2(1,1,inode) =  0.0_rp
        p1ve2(2,2,inode) =  0.0_rp
        p1ve2(3,3,inode) =  0.0_rp
        p1ve2(1,2,inode) =  factz * gpsha(inode)
        p1ve2(1,3,inode) = -facty * gpsha(inode)
        p1ve2(2,1,inode) = -factz * gpsha(inode)
        p1ve2(2,3,inode) =  factx * gpsha(inode)
        p1ve2(3,1,inode) =  facty * gpsha(inode)
        p1ve2(3,2,inode) = -factx * gpsha(inode)                
        rmom2(1,1,inode) =  rmomu(inode) 
        rmom2(2,2,inode) =  rmomu(inode) 
        rmom2(3,3,inode) =  rmomu(inode) 
        rmom2(1,2,inode) = -fact3 * gpsha(inode)
        rmom2(1,3,inode) =  fact2 * gpsha(inode)
        rmom2(2,1,inode) =  fact3 * gpsha(inode)
        rmom2(2,3,inode) = -fact1 * gpsha(inode)
        rmom2(3,1,inode) = -fact2 * gpsha(inode)
        rmom2(3,2,inode) =  fact1 * gpsha(inode)                
     end do
     do inode = 1,pnode
        idofv = (inode-1) * ndime
        do idime = 1,ndime
           idofv = idofv + 1
           do jnode = 1,pnode
              jdofv = (jnode-1) * ndime
              do jdime = 1,ndime
                 jdofv = jdofv + 1
                 elauu(idofv,jdofv) = elauu(idofv,jdofv) + gpvol &
                      * p1ve2(idime,    1,inode)                 &
                      * rmom2(    1,jdime,jnode)
                 elauu(idofv,jdofv) = elauu(idofv,jdofv) + gpvol &
                      * p1ve2(idime,    2,inode)                 &
                      * rmom2(    2,jdime,jnode)
                 elauu(idofv,jdofv) = elauu(idofv,jdofv) + gpvol &
                      * p1ve2(idime,    3,inode)                 &
                      * rmom2(    3,jdime,jnode)
              end do
           end do
           elrbu(idime,inode) = elrbu(idime,inode) + gpvol     &
                &   * (   p1ve2(idime,1,inode) * gprhs(1)      &
                &       + p1ve2(idime,2,inode) * gprhs(2)      &
                &       + p1ve2(idime,3,inode) * gprhs(3) )
        end do
     end do
  end if

end subroutine nsi_elmp13

