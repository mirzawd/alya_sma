!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_angleh(&
     kangl,pgaus,pmate,igaus,ielem,elcod,nfibe0,norma0,nshet0,pnode,lnods,gpsha)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_134
  ! NAME
  !     sld_angleh
  ! DESCRIPTION
  !    Calculate the INITIAL angle of the fiber for a spheroid geometry (with "z" being the longitufinal axis)
  !    as describe in Hunter (1989)
  !
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
  ! **WARNING** Depend on the goemetry of the spheroid used (ex. see "a_end")
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
  !
  ! USES
  ! USED BY
  !    sld_stress_model_xxx
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,ltype
  use def_solidz, only       :  parsp_sld
  use def_master, only       :  ID_EXMEDI,ITASK_ENDRUN
  

  
  implicit none

  integer(ip), intent(in)    :: pgaus,pmate,kangl,igaus,ielem,pnode,lnods(pnode)
  real(rp),    intent(out)   :: nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus)

  real(rp),    intent(in)    :: elcod(ndime,mnode),gpsha(pnode,pgaus)  

  real(rp)                   :: xcord,ycord,zcord,vuni1(3),uuni1(3),bidon

  real(rp)                   :: lenght_pnt,omega,a_epi,b_epi,l_star,&
       rotma(3,3),nfibe0T(ndime,pgaus),a_end,b_end,rcord,&
       b_qua,d_qua,c_end,a_hyp,nu,mu_end,mu_epi,mu,helix,x_endo,y_endo,z_endo,phi   

  integer(ip)                :: idime,jdime,nnode,pelty,inode,ipoin

  real(rp),parameter         :: PI = 3.141592654_rp 


  ! Coordinate of the Gauss point 

  if (ltype(ielem)==37) nnode = 8
  if (ltype(ielem)==34) nnode = 6
  if (ltype(ielem)==30) nnode = 4

  xcord=0.0_rp
  ycord=0.0_rp
  zcord=0.0_rp  

  pelty=ltype(ielem)
  do inode=1,nnode

     ipoin=lnods(inode)

     xcord = xcord + elcod(1,inode) * gpsha(inode,igaus)
     ycord = ycord + elcod(2,inode) * gpsha(inode,igaus)
     zcord = zcord + elcod(3,inode) * gpsha(inode,igaus)   

  end do



  

  ! * * * * * * * * *    
  if(kangl==1) then   !defined  !Checker ca: on le veut TI ou ORTHO?   
  ! * * * * * * * * *    

     !f0
     nfibe0(1,igaus) = 1.0_rp
     nfibe0(2,igaus) = 0.0_rp
     nfibe0(3,igaus) = 0.0_rp 
     !n0
     norma0(1,igaus) = 0.0_rp
     norma0(2,igaus) = 1.0_rp   
     norma0(3,igaus) = 0.0_rp   
     !s0    
     nshet0(1,igaus) = 0.0_rp 
     nshet0(2,igaus) = 0.0_rp        
     nshet0(3,igaus) = 1.0_rp
  
  ! * * * * * * * * *    
  else if(kangl==2) then   !defined          
  ! * * * * * * * * *    

     !f0
     nfibe0(1,igaus) = 0.0_rp
     nfibe0(2,igaus) = 1.0_rp
     nfibe0(3,igaus) = 0.0_rp 
     !n0
     norma0(1,igaus) = 0.0_rp
     norma0(2,igaus) = 0.0_rp   
     norma0(3,igaus) = 1.0_rp   
     !s0    
     nshet0(1,igaus) = 1.0_rp 
     nshet0(2,igaus) = 0.0_rp        
     nshet0(3,igaus) = 0.0_rp

  ! * * * * * * * * *    
  else if(kangl==3) then   !defined          
  ! * * * * * * * * *    

     !f0
     nfibe0(1,igaus) = 0.0_rp
     nfibe0(2,igaus) = 0.0_rp
     nfibe0(3,igaus) = 1.0_rp 
     !n0
     norma0(1,igaus) = 0.0_rp
     norma0(2,igaus) = 1.0_rp   
     norma0(3,igaus) = 0.0_rp   
     !s0    
     nshet0(1,igaus) = 1.0_rp 
     nshet0(2,igaus) = 0.0_rp        
     nshet0(3,igaus) = 0.0_rp


  ! * * * * * * * * *    
  else if(kangl==4) then !prolate spheroid with "Z" as the longitudinal axis - 
  ! angle varies in the transmural direction
  ! see http://en.wikipedia.org/wiki/Prolate_spheroidal_coordinates
  ! * * * * * * * * *    

     !geometry of the spheroid (Hunter89)   
     a_epi=parsp_sld(1,pmate)  !2.95_rp*0.47569_rp
     b_epi=parsp_sld(2,pmate)  !4.73_rp*0.47569_rp

     a_end=parsp_sld(3,pmate)  !1.65_rp*0.47569_rp 
     b_end=parsp_sld(4,pmate) !4.13_rp*0.47569_rp
    
     c_end=sqrt(b_end**2.0_rp-a_end**2.0_rp)

     if (a_epi<1.0e-10_rp) call runend("sld_angleh: Enter the dimensions of the SPHEROID option")

     ! 1) Find e1 and f_0 as usual

     lenght_pnt = sqrt(xcord**2.0_rp + ycord**2.0_rp + zcord**2.0_rp) 

     !Norme of e1
     bidon = sqrt(xcord**2.0_rp + ycord**2.0_rp)
     if (bidon==0.0_rp) write(*,*) 'ERROR MATERIAL 134  ',ielem


     !e1
     vuni1(1)=xcord/bidon
     vuni1(2)=ycord/bidon
     vuni1(3)=0.0_rp
     !nfibe0 = k x e1
     nfibe0T(1,igaus)=-vuni1(2)
     nfibe0T(2,igaus)= vuni1(1)
     nfibe0T(3,igaus)= 0.0_rp

     !2) Rotate in the Prolate spheroidal CS

     !mu ("ellipse" direction) of the endo and epi
     mu_end=log((a_end/c_end)+sqrt((a_end/c_end)**2.0_rp +1.0_rp))
     mu_epi=log((a_epi/c_end)+sqrt((a_epi/c_end)**2.0_rp +1.0_rp))

     !transform in a 2D problem in the plane passing by the point 
     rcord = sqrt(xcord**2.0_rp + ycord**2.0_rp)

     !find the point a_hyp using the eq of the hyperbola 
     b_qua=-c_end**2.0_rp-zcord**2.0_rp-rcord**2.0_rp

     d_qua=(-b_qua + sqrt(b_qua**2.0_rp-(4.0_rp*c_end**2.0_rp*zcord**2.0_rp)))/2.0_rp
     if(d_qua<0.0_rp) write(*,*) 'ERROR MATERIAL SLD_ANGEH  '
     a_hyp=sqrt(d_qua)

     if (a_hyp>c_end) then
        d_qua=(-b_qua - sqrt(b_qua**2.0_rp-(4.0_rp*c_end**2.0_rp*zcord**2.0_rp)))/2.0_rp
        if(d_qua<0.0_rp) write(*,*) 'ERROR MATERIAL SLD_ANGEH  '
        a_hyp=sqrt(d_qua)
        if (a_hyp>c_end) write(*,*) 'ERROR MATERIAL SLD_ANGEH (2) '      
     end if

     ! nu (hyperbolic direction) for this GP location (in rad)
     if(zcord>0.0_rp) then 
        nu=acos(a_hyp/c_end)
     else
        nu=acos(-a_hyp/c_end)
     end if

     ! mu ("ellipse" direction) for this GP location
     bidon=rcord/(sin(nu)*c_end)
     mu=log(bidon+sqrt(bidon**2.0_rp +1.0_rp))


     !helix angle 
     l_star= (mu-mu_end)/(mu_epi-mu_end)  
     helix=parsp_sld(5,pmate)*2.0_rp*3.141593_rp/360.0_rp 
     omega=-(2.0_rp*helix)*l_star + helix 
     !omega=0.0


     ! Angle of P w/r to +x axis
     if (xcord>=0.0_rp .and. ycord<=0.0_rp) then
        phi=atan(ycord/xcord)
     else if (xcord<=0.0_rp .and. ycord<=0.0_rp) then    
        phi=-atan(xcord/ycord)-PI/2.0_rp
     else if (xcord<=0.0_rp .and. ycord>=0.0_rp) then    
        phi= -PI+atan(ycord/xcord)   
     else if (xcord>=0.0_rp .and. ycord>=0.0_rp) then       
        phi= -1.5_rp*PI-atan(xcord/ycord)         
     end if


     !coordinate of P_endo (intersection of endocardium and the hyperbolic plan that pass by P)  
     x_endo=c_end * (0.5_rp*(exp(mu_end)-exp(-mu_end))) * sin(nu) * cos(phi)
     y_endo=c_end * (0.5_rp*(exp(mu_end)-exp(-mu_end))) * sin(nu) * sin(phi) 
     z_endo=c_end * (0.5_rp*(exp(mu_end)+exp(-mu_end))) * cos(nu)   



     ! 3) Rotate f_0T around u (the unit vector "u" is from P_endo to P)

     !The rotation matrix
     uuni1(1)= xcord-x_endo
     uuni1(2)= ycord-y_endo
     uuni1(3)= zcord-z_endo 
     bidon=sqrt(uuni1(1)**2.0_rp + uuni1(2)**2.0_rp + uuni1(3)**2.0_rp)
     uuni1(1)=uuni1(1)/bidon
     uuni1(2)=uuni1(2)/bidon
     uuni1(3)=uuni1(3)/bidon     

     rotma(1,1) = uuni1(1)**2.0_rp + (1.0_rp-uuni1(1)**2.0_rp)*cos(omega)
     rotma(1,2) = uuni1(1)*uuni1(2)*(1.0_rp-cos(omega)) - uuni1(3)*sin(omega)
     rotma(1,3) = uuni1(3)*uuni1(1) + uuni1(2)*sin(omega)

     rotma(2,1) = uuni1(1)*uuni1(2)*(1.0_rp-cos(omega)) + uuni1(3)*sin(omega)
     rotma(2,2) = uuni1(2)**2.0_rp + (1.0_rp-uuni1(2)**2.0_rp)*cos(omega)
     rotma(2,3) = uuni1(2)*uuni1(3)*(1.0_rp-cos(omega))-uuni1(1)*sin(omega)

     rotma(3,1) = uuni1(1)*uuni1(3)*(1.0_rp-cos(omega)) - uuni1(2)*sin(omega)
     rotma(3,2) = uuni1(2)*uuni1(3)*(1.0_rp-cos(omega))+ uuni1(1)*sin(omega)
     rotma(3,3) = uuni1(3)**2.0_rp + (1.0_rp-uuni1(3)**2.0_rp)*cos(omega)   


     do idime=1,ndime
        nfibe0(idime,igaus)=0.0_rp
        do jdime=1,ndime
           nfibe0(idime,igaus) = nfibe0(idime,igaus) + rotma(idime,jdime)* nfibe0T(jdime,igaus)
        end do
     end do


     ! 4) Find s_0 anf n_0 (same equations)

     !norma0 
     norma0(1,igaus) = uuni1(1)
     norma0(2,igaus) = uuni1(2)
     norma0(3,igaus) = uuni1(3)      


     !nshet0= n x f     
     nshet0(1,igaus) = (norma0(2,igaus)*nfibe0(3,igaus)) -&
          (norma0(3,igaus)*nfibe0(2,igaus))

     nshet0(2,igaus) = -(norma0(1,igaus)*nfibe0(3,igaus)) +&
          (norma0(3,igaus)*nfibe0(1,igaus))

     nshet0(3,igaus) = (norma0(1,igaus)*nfibe0(2,igaus)) -&
          (norma0(2,igaus)*nfibe0(1,igaus))        


  end if !kangl


100 format (5(E16.8,','))
101 format (6(F16.8,','))

end subroutine sld_angleh
