!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_angle2(imate,kangl)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_134
  ! NAME
  !     sld_angleh
  ! DESCRIPTION
  !    Produce an input file of fiber to be use with exmedi
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
  use def_domain  
  use def_parame
  use def_master
  use def_solver
  use def_solidz
  
  implicit none

  integer(ip), intent(in)    :: kangl,imate
  
  real(rp)                   :: xcord,ycord,zcord,vuni1(3),uuni1(3),bidon,nfibe0(ndime)
  
  real(rp)                   :: lenght_pnt,omega,a_epi,b_epi,l_star,&
       rotma(3,3),nfibe0T(ndime),a_end,b_end,rcord,&
       b_qua,d_qua,c_end,a_hyp,nu,mu_end,mu_epi,mu,helix,x_endo,y_endo,z_endo,phi   
  
  integer(ip)                :: idime,jdime,inode,ipoin
  
!  real(rp),parameter         :: PI = 3.141592654_rp 
  



  ! * * * * * * * * *    
  if(kangl==1) then   !defined          
  ! * * * * * * * * *    

     !f0
     nfibe0(1) = 1.0_rp
     nfibe0(2) = 0.0_rp
     nfibe0(3) = 0.0_rp 
     write(7777,100) inode,nfibe0(1),nfibe0(2),nfibe0(3)
  
  ! * * * * * * * * *    
  else if(kangl==2) then   !defined          
  ! * * * * * * * * *    

     !f0
     nfibe0(1) = 0.0_rp
     nfibe0(2) = 1.0_rp
     nfibe0(3) = 0.0_rp 
     write(7777,100) inode,nfibe0(1),nfibe0(2),nfibe0(3)

  ! * * * * * * * * *    
  else if(kangl==3) then   !defined          
  ! * * * * * * * * *    

     !f0
     nfibe0(1) = 0.0_rp
     nfibe0(2) = 0.0_rp
     nfibe0(3) = 1.0_rp 
     write(7777,100) inode,nfibe0(1),nfibe0(2),nfibe0(3)  
     
     ! * * * * * * * * *    
  else if(kangl==4) then !prolate spheroid with "Z" as the longitudinal axis - 
     ! angle varies in the transmural direction
     ! see http://en.wikipedia.org/wiki/Prolate_spheroidal_coordinates
     ! * * * * * * * * *    
     

     !geometry of the spheroid (Hunter89)   


     a_epi=parsp_sld(1,imate)  
     b_epi=parsp_sld(2,imate) 

     a_end=parsp_sld(3,imate)  
     b_end=parsp_sld(4,imate) 
    
     c_end=sqrt(b_end**2.0_rp-a_end**2.0_rp)

     do ipoin = 1,npoin

        
        !coordinate of the node (if called by a "nodal loop")
        xcord = coord(1,inode)
        ycord = coord(2,inode)
        zcord = coord(3,inode)
                
        ! 1) Find e1 and f_0 as usual
        
        lenght_pnt = sqrt(xcord**2.0_rp + ycord**2.0_rp + zcord**2.0_rp) 
        
        !Norme of e1
        bidon = sqrt(xcord**2.0_rp + ycord**2.0_rp)
        if (bidon==0.0_rp) then
           nfibe0(1)= 1.0_rp
           nfibe0(2)= 0.0_rp
           nfibe0(3)= 0.0_rp

           !write in a file 
           write(7777,100) inode,nfibe0(1),nfibe0(2),nfibe0(3)
           
        else
           
           !e1
           vuni1(1)=xcord/bidon
           vuni1(2)=ycord/bidon
           vuni1(3)=0.0_rp
           !nfibe0 = k x e1
           nfibe0T(1)=-vuni1(2)
           nfibe0T(2)= vuni1(1)
           nfibe0T(3)= 0.0_rp
           
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
           helix=1.0472_rp  !+/- 60 deg (1.047 rad)
           !helix=0.1745_rp  !+/- 10 deg (1.047 rad)
           omega=-(2.0_rp*helix)*l_star + helix  !+/- 60 deg (1.047 rad)
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
              nfibe0(idime)=0.0_rp
              do jdime=1,ndime
                 nfibe0(idime) = nfibe0(idime) + rotma(idime,jdime)* nfibe0T(jdime)
              end do
           end do
           !Norm - to check 
           !bidon = sqrt(nfibe0(1)**2.0_rp + nfibe0(2)**2.0_rp +nfibe0(3)**2.0_rp)       
           !write(*,*) 'NORME= ',bidon           

           !write in a file 
           write(7777,100) inode,nfibe0(1),nfibe0(2),nfibe0(3)
                  
        end if
        
     end do !loop on all the nodes 
     
  end if !kangl
  
  
100 format (i9,3(E16.8,' '))
  
  
end subroutine sld_angle2
