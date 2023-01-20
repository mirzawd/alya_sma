!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_radvuf 
!-----------------------------------------------------------------------
!****f* Temper/tem_radvuf
! NAME 
!    tem_radvuf
! DESCRIPTION
!    This routine prepares for a new time step of the temperature
!    equation      
! USED BY
!    tem_begste
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master 
  use      def_domain
  use      def_temper
  use      mod_memchk
  implicit none
  integer(ip)  :: iboun,jboun,kboun,lboun
  integer(ip)  :: ipoi1,ipoi2,ipoi3,ipok1,ipok2
  integer(8)   :: memor_rad(2)
  integer(4)   :: istat
  real(rp)     :: bcooc(ndime,nboun),barea(nboun)
  real(rp)     :: vect1(ndime,nboun),vect2(ndime,nboun)
  real(rp)     :: exwor(ndime,nboun)
  real(rp)     :: oneth,x1,y1,z1,x2,y2,z2,cosbi,cosbj,exwno
  real(rp)     :: venij,deter,vecij(ndime),vncij(ndime)
  real(rp)     :: vectt(ndime)
  real(rp)     :: a,b,c,d,e,f,xk1,yk1,xk2,yk2,xi,xj,yi,yj,zi,zj
  real(rp)     :: tij,tk,t,x,y,z,Fnorm,sumFn,Fij
  real(rp)     :: nx,ny,nz,cpu_radia,rbyte
  real(rp)     :: cpu_refe1,cpu_refe2
  character(6) :: lbyte 
!
! Allocate memory
!
  memor_rad=0
  call memchk(zero,istat,memor_rad,'VIEWF_TEM','tem_radvuf',viewf_tem)
  allocate(viewf_tem(nboun),stat=istat)
  kboun=nboun
  do iboun=1,nboun
     kboun=kboun-1
     allocate(viewf_tem(iboun)%a(kboun),stat=istat)
     call memchk(zero,istat,memor_rad,'VIEWF_TEM','dombou',viewf_tem(iboun)%a)
  end do
!
! Initialization
! 
  oneth=1.0_rp/3.0_rp 
  bcooc=0.0_rp
  if(ndime==3) then
     if((nelty/=1.or.nnode(1)/=3))&
          call runend('tem_radvuf: ELEMENT NOT IMPLEMENTED')
  end if

  call cputim(cpu_refe1)

  if(ndime==2) then
     do iboun=1,nboun
        ipoi1=lnodb(1,iboun)
        ipoi2=lnodb(2,iboun)
        barea(iboun)=sqrt(&
             &  (coord(1,ipoi1)-coord(1,ipoi2)) &
             & *(coord(1,ipoi1)-coord(1,ipoi2)) &
             & +(coord(2,ipoi1)-coord(2,ipoi2)) &
             & *(coord(2,ipoi1)-coord(2,ipoi2)))
        exwor(1,iboun)= coord(2,ipoi1)-coord(2,ipoi2)
        exwor(2,iboun)=-coord(1,ipoi1)+coord(1,ipoi2)
        call vecuni(2_ip,exwor,exwno)
        bcooc(1,iboun)=0.50_rp*(coord(1,ipoi1)+coord(1,ipoi2))
        bcooc(2,iboun)=0.50_rp*(coord(2,ipoi1)+coord(2,ipoi2))
     end do
     do iboun=1,nboun
        xi=bcooc(1,iboun)
        yi=bcooc(2,iboun)
        lboun=0
        jboun2d: do jboun=iboun+1,nboun
           lboun=lboun+1
           xj=bcooc(1,jboun)
           yj=bcooc(2,jboun)
           vecij(1)=bcooc(1,jboun)-bcooc(1,iboun)
           vecij(2)=bcooc(2,jboun)-bcooc(2,iboun)
           call vecuni(2_ip,vecij,venij)
           cosbi= exwor(1,iboun)*vecij(1)+exwor(2,iboun)*vecij(2)
           cosbj=-exwor(1,jboun)*vecij(1)-exwor(2,jboun)*vecij(2)
           if(cosbi<=zetem.or.cosbj<=zetem) cycle jboun2d
           kboun2d: do kboun=1,nboun
              if(kboun==iboun.or.kboun==jboun) cycle kboun2d
              ipok1 = lnodb(1,kboun)
              ipok2 = lnodb(2,kboun)
              xk1   = coord(1,ipok1)
              yk1   = coord(2,ipok1)
              xk2   = coord(1,ipok2)
              yk2   = coord(2,ipok2)
              a     =  xj  - xi
              b     = -xk2 + xk1
              c     =  yj  - yi
              d     = -yk2 + yk1
              deter =  a*d - b*c
              if(abs(deter)>zetem) then
                 deter = 1.0_rp/deter
                 e     = xk1-xi
                 f     = yk1-yi
                 tij   = deter*( d*e-b*f)
                 tk    = deter*(-c*e+a*f)
                 if(    tij>0.0_rp.and.tij<1.0_rp.and.&
                      &  tk>0.0_rp.and. tk<1.0_rp) cycle jboun2d
              end if
           end do kboun2d
           viewf_tem(iboun)%a(lboun)=cosbi*cosbj/(pi*venij*venij)
        end do jboun2d
     end do
  else
     !
     ! Compute centers and areas
     !  
     do iboun=1,nboun
        ipoi1=lnodb(1,iboun)
        ipoi2=lnodb(2,iboun)
        ipoi3=lnodb(3,iboun)          
        ! e1 and e2
        x1=coord(1,ipoi2)-coord(1,ipoi1)
        y1=coord(2,ipoi2)-coord(2,ipoi1)
        z1=coord(3,ipoi2)-coord(3,ipoi1)
        x2=coord(1,ipoi3)-coord(1,ipoi1)
        y2=coord(2,ipoi3)-coord(2,ipoi1)
        z2=coord(3,ipoi3)-coord(3,ipoi1)
        vect1(1,iboun)=x1
        vect1(2,iboun)=y1
        vect1(3,iboun)=z1
        vect2(1,iboun)=x2
        vect2(2,iboun)=y2
        vect2(3,iboun)=z2
        ! Normals
        call vecpro(vect1(1,iboun),vect2(1,iboun),exwor(1,iboun),3)     
        call vecuni(3_ip,exwor(1,iboun),exwno)
        ! Area
        barea(iboun)=0.50_rp*exwno
        ! Central point
        bcooc(1:3,iboun)=oneth*&
             (coord(1:3,ipoi1)+coord(1:3,ipoi2)+coord(1:3,ipoi3))
     end do
     !
     ! Compute shadowing
     !
     do iboun=1,nboun
        lboun=0
        xi=bcooc(1,iboun)
        yi=bcooc(2,iboun)
        zi=bcooc(3,iboun)
        jboun3d: do jboun=iboun+1,nboun
           lboun=lboun+1
           xj=bcooc(1,jboun)
           yj=bcooc(2,jboun)
           zj=bcooc(3,jboun)
           if(iboun==jboun) cycle jboun3d
           ! Scalar product ni.rij
           vecij(1)=xj-xi
           vecij(2)=yj-yi
           vecij(3)=zj-zi
           vncij   =vecij
           call vecuni(3_ip,vncij,venij)
           cosbi=dot_product(exwor(1:3,iboun), vncij(1:3))
           cosbj=dot_product(exwor(1:3,jboun),-vncij(1:3))
           
           if(cosbi<=zetem.or.cosbj<=zetem) cycle jboun3d

           ! Check if boundaries iboun and jboun see each others
           kboun3d: do kboun=1,nboun              
              if(kboun/=iboun.and.kboun/=jboun) then
                 nx    = exwor(1,kboun)
                 ny    = exwor(2,kboun)
                 nz    = exwor(3,kboun)
                 ipoi1 = lnodb(1,kboun)
                 f     = nx*vecij(1)+ny*vecij(2)+nz*vecij(3)
                 if(f>=zetem) then
                    t=(nx*(coord(1,ipoi1)-xi)+ny*(coord(2,ipoi1)-yi)+nz*(coord(3,ipoi1)-zi))/f
                    if(t>=0.0_rp.and.t<=1.0_rp) then
                       x=vecij(1)*t+xi
                       y=vecij(2)*t+yi
                       z=vecij(3)*t+zi
                       vectt(1)=x-coord(1,ipoi1)
                       vectt(2)=y-coord(2,ipoi1)
                       vectt(3)=z-coord(3,ipoi1)
                       if(abs(nz)>zetem) then
                          deter=  1.0_rp/( vect1(1,kboun)*vect2(2,kboun)&
                               &          -vect1(2,kboun)*vect2(1,kboun))
                          a=deter*( vect2(2,kboun)*vectt(1)-vect2(1,kboun)*vectt(2))
                          if(a<0.0_rp.or.a>1.0_rp) cycle kboun3d
                          b=deter*(-vect1(2,kboun)*vectt(1)+vect1(1,kboun)*vectt(2))
                          if(b<0.0_rp.or.b>1.0_rp.or.a+b>1) cycle kboun3d                        
                       else if(abs(ny)>zetem) then
                          deter=  1.0_rp/( vect1(1,kboun)*vect2(3,kboun)&
                               &          -vect1(3,kboun)*vect2(1,kboun))
                          a=deter*( vect2(3,kboun)*vectt(1)-vect2(1,kboun)*vectt(3))
                          if(a<0.0_rp.or.a>1.0_rp) cycle kboun3d
                          b=deter*(-vect1(3,kboun)*vectt(1)+vect1(1,kboun)*vectt(3))
                          if(b<0.0_rp.or.b>1.0_rp.or.a+b>1) cycle kboun3d                       
                       else
                          deter=  1.0_rp/( vect1(2,kboun)*vect2(3,kboun)&
                               &          -vect1(3,kboun)*vect2(2,kboun))
                          a=deter*( vect2(3,kboun)*vectt(2)-vect2(2,kboun)*vectt(3))
                          if(a<0.0_rp.or.a>1.0_rp) cycle kboun3d
                          b=deter*(-vect1(3,kboun)*vectt(2)+vect1(2,kboun)*vectt(3))
                          if(b<0.0_rp.or.b>1.0_rp.or.a+b>1) cycle kboun3d                        
                       end if
                       cycle jboun3d
                    end if
                 end if
              end if
           end do kboun3d
           ! Symmetric view factor: F'ij
           viewf_tem(iboun)%a(lboun)=cosbi*cosbj/(pi*venij*venij)
           if(viewf_tem(iboun)%a(lboun)<0) then
              write(*,*) 'WARNING: F('&
                   //adjustl(trim(intost(iboun)))//','&
                   //adjustl(trim(intost(jboun)))//') IS NEGATIVE'
              viewf_tem(iboun)%a(lboun)=0.0_rp
           end if

        end do jboun3d
     end do
  endif
!
! Modifies Fij' such that Sum_i=1^N F'ij*Aj=1
!
  sumFn=0.0_rp
  do iboun=1,nboun
     Fnorm=0.0_rp 
     do jboun=1,nboun        
        if(jboun>iboun) then
           lboun=jboun-iboun
           Fij=viewf_tem(iboun)%a(lboun)
        else if(jboun<iboun) then
           lboun=iboun-jboun
           Fij=viewf_tem(jboun)%a(lboun)
        else
           Fij=0.0_rp
        end if
        Fnorm=Fnorm+barea(jboun)*Fij
     end do
     sumFn=sumFn+Fnorm
     Fnorm=1.0_rp/Fnorm
     do jboun=1,nboun
        if(jboun>iboun) then
           lboun=jboun-iboun
           viewf_tem(iboun)%a(lboun)=viewf_tem(iboun)%a(lboun)*Fnorm
        else if(jboun<iboun) then
           lboun=iboun-jboun
           viewf_tem(jboun)%a(lboun)=viewf_tem(jboun)%a(lboun)*Fnorm
        end if        
     end do
  end do
  sumFn=sumFn/real(nboun,rp)
!
! Write results
!
  call cputim(cpu_refe2)
  cpu_radia=cpu_refe2-cpu_refe1
  mem_modul(1,modul)=mem_modul(1,modul)+memor_rad(1)
  mem_modul(2,modul)=max(mem_modul(1,modul),mem_modul(2,modul))
  if(memor_rad(1)>=1024*1024*1024) then
     rbyte=1024.0_rp*1024.0_rp*1024.0_rp
     lbyte='Gbytes'
  else if(memor_rad(1)>=1024*1024) then 
     rbyte=1024.0_rp*1024.0_rp
     lbyte='Mbytes'     
  else if(memor_rad(1)>=1024) then 
     rbyte=1024.0_rp
     lbyte='kbytes'          
  else  
     rbyte=1.0_rp
     lbyte=' bytes'     
  end if
  write(momod(modul)%lun_outpu,110)cpu_radia,real(memor_rad(1),rp)/rbyte,lbyte,sumFn
  
110 format(/,&
         & 5x, '>>>  RADIATION CALCULATION:',//,&
         & 10x,'View factors:',/,&
         & 10x,'-------------',/,&
         & 10x,'CPU time=                    ',e12.6,/,&
         & 10x,'Memory=                      ',f6.1,' ',a6,/,&
         & 10x,'1/N Sum_i^N sum_j^N F''ijAj= ',f5.2)
111 format(&
         & 10x,a,10e12.6)
  
end subroutine tem_radvuf
