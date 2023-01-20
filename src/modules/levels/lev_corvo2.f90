!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_corvo2
  !-----------------------------------------------------------------------
  !****f* Levels/lev_corvo2
  ! NAME 
  !    lev_corvo2
  !    
  ! DESCRIPTION
  !    Imposition of constant mass (cf Loehner 2006)                                
  ! USES
  !    
  ! USED BY
  !    lev_updunk
  !***
  !-----------------------------------------------------------------------

  use  def_parame
  use  def_master
  use  def_levels
  use  def_domain
  use  def_solver
  use  mod_gradie    
  use mod_messages, only : livinf
  implicit none

  integer(ip)             :: ielem,idime            ! Indices and dimensions
  integer(ip)             :: inode,ipoin,compt
  integer(ip)             :: pelty,pnode
  integer(ip)             :: itcor
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode)
  real(rp)                :: norm
  real(rp)                :: corvol, errvol, ervolr 
  real(rp)                :: eps, ervold, p, x, y, z, vx, vy, vz, t

  call livinf(59_ip,' VOLUME CORRECTION LOEHNER ',0_ip)

  call grasca(fleve,norml_lev)

  if(INOTMASTER) then

     do ipoin = 1,npoin

        norm = 0.0_rp
        do idime = 1,ndime
           norm = norm + norml_lev(idime,ipoin)*norml_lev(idime,ipoin)
        enddo
        norm= sqrt(norm)
        do idime = 1,ndime
           norml_lev(idime,ipoin) = norml_lev(idime,ipoin)/norm
        enddo

        normv_lev(ipoin) = 0.0_rp
        icupt_lev(ipoin) = 0_ip

     enddo

     do ielem=1,nelem

        !
        ! Element dimensions
        !
        pelty=ltype(ielem)
        pnode=nnode(pelty)

        !
        ! Gather operations
        !
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           ellev(inode)=fleve(ipoin,1)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do

        compt=0_ip

        !
        ! Determine type of element
        ! type= 1 empty element (phi<0) in each node
        ! type= 2 full element (phi>0) in each node
        ! type= 3 partly full element (phi<0) in some nodes and (phi>0) in some nodes
        !
        do inode=1,pnode
           if(ellev(inode)>=0.0_rp) then
              compt=compt+1
           endif
        end do

        if((compt/=0_ip).and.(compt/=pnode)) then

           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              icupt_lev(ipoin)=1_ip
           end do

        endif

     enddo

     if(kfl_advec_lev==7) then
        p=3.141592653589793_rp
        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           vx=-sin(p*x)*sin(p*x)*2*sin(p*y)*cos(p*y)
           vy=sin(p*y)*sin(p*y)*2*sin(p*x)*cos(p*x)
           normv_lev(ipoin)=normv_lev(ipoin)+vx*norml_lev(1,ipoin)
           normv_lev(ipoin)=normv_lev(ipoin)+vy*norml_lev(2,ipoin)
           normv_lev(ipoin)=abs(normv_lev(ipoin))
        end do
     else if(kfl_advec_lev==8) then
        do ipoin=1,npoin
           x = coord(1,ipoin)
           y = coord(2,ipoin)
           vx = (50_rp-y)/100_rp
           vy = (x-50_rp)/100_rp
           normv_lev(ipoin)=normv_lev(ipoin)+vx*norml_lev(1,ipoin)
           normv_lev(ipoin)=normv_lev(ipoin)+vy*norml_lev(2,ipoin)
           normv_lev(ipoin)=abs(normv_lev(ipoin))
        enddo
     else if(kfl_advec_lev==10_ip) then
        t=(ittim-1)/dtinv
        do ipoin=1,npoin
           x = coord(1,ipoin)
           y = coord(2,ipoin)
           z = coord(3,ipoin)
           vx = 2*sin(pi*x)*sin(pi*x)*sin(2*pi*y)*sin(2*pi*z)*cos(pi*t/3.0)
           vy = -sin(2*pi*x)*sin(pi*y)*sin(pi*y)*sin(2*pi*z)*cos(pi*t/3.0)
           vz = -sin(2*pi*x)*sin(2*pi*y)*sin(pi*z)*sin(pi*z)*cos(pi*t/3.0)
           normv_lev(ipoin)=normv_lev(ipoin)+vx*norml_lev(1,ipoin)
           normv_lev(ipoin)=normv_lev(ipoin)+vy*norml_lev(2,ipoin)
           normv_lev(ipoin)=normv_lev(ipoin)+vz*norml_lev(3,ipoin)
           normv_lev(ipoin)=abs(normv_lev(ipoin))
        enddo
     else 
        do ipoin = 1,npoin
           do idime=1,ndime
              normv_lev(ipoin)=normv_lev(ipoin)+veloc(idime,ipoin,1)*norml_lev(idime,ipoin)
           end do
           normv_lev(ipoin)=abs(normv_lev(ipoin))
        enddo
     endif


  endif

  call lev_calvol()

  errvol=volit_lev-volrf_lev
  ervolr=errvol/volrf_lev

  itcor=0_ip
  eps=1e-06_rp
  corvol=0.0_rp

  do while((abs(ervolr)>1.e-08).and.(itcor<10_ip))


     itcor=itcor+1_ip

     if(INOTMASTER) then
        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)+eps*normv_lev(ipoin)
        enddo
     endif
     call lev_calvol()

     ervold=(volit_lev-volrf_lev-errvol)/eps

     if(ervold/=0.0) then 

        if(INOTMASTER) then
           do ipoin=1,npoin
              if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)-(eps+corvol)*normv_lev(ipoin)
           enddo
        endif

        corvol=corvol-errvol/ervold
        if(INOTMASTER) then
           do ipoin=1,npoin
              if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)+corvol*normv_lev(ipoin)
           enddo
        endif
        call lev_calvol()
        errvol=volit_lev-volrf_lev
        ervolr=errvol/volrf_lev

     else

        itcor = 10_ip
        if(INOTMASTER) then
           do ipoin=1,npoin
              if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)-(eps+corvol)*normv_lev(ipoin)
           enddo
        endif
        call livinf(59_ip,' LEV_CORVO2: VOLUME CORRECTION FAILED ',0_ip)

     endif

  enddo

end subroutine lev_corvo2
