!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_potent

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_potent
  ! NAME 
  !    nsi_optent
  ! DESCRIPTION
  !    Compute potential flow
  !    This routine computes the velocity u of an irrotational flow:
  !
  !    1. irrotational flow: grad x u = 0
  !    2. => u=grad(phi)
  !    3. Apply continuity div(u)=0
  !    4. => div(grad(phi))=0
  !
  !    * Let phi be the solution of the following Laplace equation
  !
  !      +-
  !      | div[grad(phi)] =  0                    in Omega
  !     -|           phi  =  u(x)*x+u(y)*y+u(z)*z where u prescribed
  !      |   grad(phi).n  =  0                    on walls
  !      +-
  !
  !    * Then the velocity is computed as
  !
  !      u=grad(phi) 
  !
  ! USES
  ! USED BY
  !    tem_matrix
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_nastin
  use      def_solver
  use      mod_gradie
  use      mod_bouder
  implicit none

  real(rp)    :: cartd(ndime,mnode) 
  real(rp)    :: cartb(ndime,mnode)
  real(rp)    :: xjaci(ndime,ndime) 
  real(rp)    :: xjacm(ndime,ndime) 

  real(rp)    :: elmat(mnode,mnode),wmatr(mnode,mnode) ! Element matrices
  real(rp)    :: elrhs(mnode),wrhsi(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: shapp(mnode)

  real(rp)    :: elcod(ndime,mnode)

  integer(ip) :: ielem,inode,ipoin,jnode              ! Indices and dimensions
  integer(ip) :: iboun,inodb,igaub,pnodb,jnodb
  integer(ip) :: igaus,idime,knodb
  integer(ip) :: pelty,pnode,pgaus,pblty
  integer(ip) :: kfl_poten_nsi(npoin)
  integer(ip) :: kfl_boten_nsi(nboun)

  real(rp)    :: dvolu,detjm,adiag,udotn         ! Values at Gauss points
  real(rp)    :: eucta,dsurf,acvel(3)
 
  if( IMASTER ) return
  call runend('REPROGRAM SOLVER')
  !
  ! Boundary conditions
  !
  do ipoin=1,npoin
     if(&
          kfl_fixno_nsi(    1,ipoin)==1.or.&
          kfl_fixno_nsi(    2,ipoin)==1.or.&
          kfl_fixno_nsi(ndime,ipoin)==1) then
        kfl_poten_nsi(ipoin)=1
     else
        kfl_poten_nsi(ipoin)=0
     end if
  end do
  do iboun=1,nboun
     pblty=ltypb(iboun)
     pnodb=nnode(pblty)
     knodb=0
     do inodb=1,pnodb
        ipoin=lnodb(inodb,iboun)
        knodb=knodb+kfl_poten_nsi(ipoin)
     end do
     if(knodb==pnodb) then
        kfl_boten_nsi(iboun)=0
     else
        kfl_boten_nsi(iboun)=1
     end if
  end do
  !
  ! Prepare solver
  !
  call runend('NSI_POTENT: PREPARE SOLVER')
  !call nsi_inisol()
  !
  ! Loop over elements
  !
  elements: do ielem=1,nelem
     ! Initialize
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)
     elmat=0.0_rp
     elrhs=0.0_rp
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        do idime=1,ndime
           elcod(idime,inode)=coord(idime,ipoin)
        end do
     end do
     do igaus=1,pgaus       
        ! Cartesian derivatives and Jacobian
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
             elcod,cartd,detjm,xjacm,xjaci)
        dvolu=elmar(pelty)%weigp(igaus)*detjm
        do inode=1,pnode
           do jnode=1,pnode
              do idime=1,ndime
                 elmat(inode,jnode)=elmat(inode,jnode)&
                      +dvolu*cartd(idime,inode)*cartd(idime,jnode)
              end do
           end do
        end do
     end do
     ! Impose phi(1)=0  
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        if(ipoin==1) then
           adiag=elmat(inode,inode)
           elmat(inode,1:pnode)=0.0_rp
           elmat(inode,inode)=adiag
           elrhs(inode)=0.0_rp
        end if
     end do
     ! Assembly
     call csrase(elmat,amatr,1_ip,pnode,mnode,lnods(1,ielem),2_ip)
  end do elements

  goto 10
  !      
  ! Loop over boundaries
  !
  boundaries: do iboun=1,nboun     
     pblty=ltypb(iboun)
     pnodb=nnode(pblty)
     ielem=lelbo(iboun)
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)
     ! Inititalize
     elmat=0.0_rp
     elrhs=0.0_rp        
     ! Gather operations
     bocod(1:ndime,1:pnodb)=coord(1:ndime,lnodb(1:pnodb,iboun))        
     do igaub=1,ngaus(pblty)           
        wmatr=0.0_rp
        wrhsi=0.0_rp                      
        ! Jacobian
        call bouder(pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
             bocod,baloc,eucta)
        dsurf=elmar(pblty)%weigp(igaub)*eucta
        if(kfl_boten_nsi(iboun)==1) then
           do idime=1,ndime
              acvel(idime)=0.0_rp
           end do
           do inodb=1,pnodb
              ipoin=lnodb(inodb,iboun)
              do idime=1,ndime
                 acvel(idime)=acvel(idime)+&
                      veloc(idime,ipoin,1)*elmar(pblty)%shape(inodb,igaub)
              end do
           end do
           do idime=1,ndime
              udotn=udotn+acvel(idime)*baloc(idime,ndime)
           end do
           do inodb=1,pnodb  
              inode=lboel(inodb,iboun)
              wrhsi(inode)=wrhsi(inode)&
                   +elmar(pblty)%shape(inodb,igaub)*udotn
           end do
        else
           call cartbo(&
                1_ip,lboel(1,iboun),elmar(pblty)%shape,elmar(pelty)%shaga,&
                cartd,elmar(pelty)%shape,shapp,cartb,&
                pnodb,pnode,pgaus)
           do inodb=1,pnodb
              inode=lboel(inodb,iboun)                 
              do jnodb=1,pnodb
                 jnode=lboel(jnodb,iboun)
                 do idime=1,ndime
                    wmatr(inode,jnode)=wmatr(inode,jnode)&
                         +cartb(idime,jnode)*baloc(idime,ndime)&
                         *elmar(pblty)%shape(inodb,igaub)
                 end do
              end do
           end do
        end if
        elmat=elmat+wmatr*dsurf           
        elrhs=elrhs+wrhsi*dsurf           
     end do

     ! Prescribe Dirichlet boundary conditions
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        if(ipoin==1) then
           adiag=elmat(inode,inode)
           elmat(inode,1:pnode)=0.0_rp
           elmat(inode,inode)=adiag
           elrhs(inode)=0.0_rp
        end if
     end do
     ! Assembly
     call csrase(elmat,amatr,1_ip,pnode,mnode,lnods(1,ielem),2_ip)

  end do boundaries
10 continue
  !
  ! Solve the algebraic system.
  !
  call solver(rhsid,unkno,amatr,pmatr)
  !
  ! Compute u=grad(phi)
  !
  call gradie(unkno(1:npoin),veloc(1:ndime,1:npoin,nprev_nsi))
  !
  ! Rotate unknown
  !     
  !call nsi_rotunk
  do ipoin=1,npoin
     write(*,*) (veloc(idime,ipoin,nprev_nsi),idime=1,ndime)
  end do
  stop

end subroutine nsi_potent
