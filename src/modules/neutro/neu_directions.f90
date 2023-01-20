!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_directions.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Directions
!> @details Directions
!> @} 
!-----------------------------------------------------------------------

subroutine neu_directions()

  use def_kintyp,          only : ip,rp
  use def_master,          only : momod,modul,INOTSLAVE!,INOTMASTER
  use def_domain,          only : ndime
  use mod_maths_geometry,  only : maths_angle_of_a_vector
  use def_neutro
  !, only : kfl_icosa_neu,direc_neu
  !use def_neutro, only : kfl_snord_neu,num_directions_neu
  implicit none
  integer(ip) :: idire
  real(rp)    :: angl1,angl2!,funtita,funphi
  external :: neu_ordico, neu_snordi
  real(rp), external :: funtita, funphi


  if( kfl_icosa_neu == 1 )  then 
     !
     ! ICOSAHEDRON DIRECTION
     !
     call neu_ordico()  
     
  else if( kfl_snord_neu /= 0 ) then
     !
     ! SN8 DIRECTIONS
     !
     call neu_snordi()
     
  end if
  !
  ! Write direction in output file
  !
!  if( INOTSLAVE ) then
!     allocate(phi_neu(num_directions_neu),tita_neu(num_directions_neu))
     do idire = 1,num_directions_neu 
        call maths_angle_of_a_vector(ndime,direc_neu(1:ndime,idire),angl1,angl2)
       
      !   tita_neu(idire) = funtita(direc_neu(1,idire),direc_neu(2,idire),direc_neu(3,idire))
      !   phi_neu(idire) = funphi(direc_neu(1,idire),direc_neu(2,idire),direc_neu(3,idire))
        tita_neu(idire) = funtita(direc_neu(3,idire))
        phi_neu(idire) = funphi(direc_neu(1,idire),direc_neu(2,idire))
        if( INOTSLAVE ) write(momod(modul)%lun_outpu,'(A, I3, A, 2F20.15)') '        Direction: ',idire,&
                                                         ' angles= ',phi_neu(idire),tita_neu(idire)
     end do
     if( INOTSLAVE ) flush(momod(modul)%lun_outpu)
!  end if


  
end subroutine neu_directions


! real(rp) function funtita(omex,omey,omez)
real(rp) function funtita(omez)
  use def_kintyp_basic, only : rp
implicit none
real(rp) :: omez!,omex,omey

   funtita =  acos(omez)

end function funtita

! real(rp) function funphi(omex,omey,omez)
real(rp) function funphi(omex,omey)
  use def_kintyp_basic, only : rp
  use def_parame, only : pi
  implicit none
  real(rp) :: omex,omey!,omez

!   funphi=  atan(omey/omex)
  funphi=  atan2(omey, omex)+pi

end function funphi

subroutine neu_snordi()  

  !-----------------------------------------------------------------------
  ! constructs SN ordinates
  !-----------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_domain, only : ndime
  use def_neutro, only : kfl_snord_neu
  use def_neutro, only : num_directions_neu
  use def_neutro, only : direc_neu
  use def_neutro, only : weigd_neu
  implicit none

  real(rp)    :: dirpr(12),weid(24) 
  integer(ip) :: idire,idime,N8
  external :: LDFE

  select case ( kfl_snord_neu )

  case ( 4_ip )

     !-------------------------------------------------------------------
     !
     ! SN4
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.2958759_rp
     dirpr(2) = 0.9082483_rp
     !
     !   first octant
     !
     direc_neu (1, 1) = dirpr(1)
     direc_neu (1, 2) = dirpr(1)
     direc_neu (1, 3) = dirpr(2)

     direc_neu (2, 1) = dirpr(1)
     direc_neu (2, 2) = dirpr(2)
     direc_neu (2, 3) = dirpr(1)

     direc_neu (3, 1) = dirpr(2)
     direc_neu (3, 2) = dirpr(1)
     direc_neu (3, 3) = dirpr(1)


     weigd_neu(1)=0.523598775598299_rp
     weigd_neu(2)=0.523598775598299_rp
     weigd_neu(3)=0.523598775598299_rp
     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =4, 6
        direc_neu(1,idire)= - direc_neu(1,idire-3)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-3)
        end do
        weigd_neu(idire) = weigd_neu(idire-3)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     ! ¡CHANGE! (E. Goldberg, 27/06/2022): INSTEAD OF OCTANTS 3 AND 4, 
     !            I'M USING HERE OCTANTS 7 AND 8 FOR ALL SN. 
     !            This is just a change in the order
     !            from 1,2,3,4,5,6,7,8 to 1,2,7,8,5,6,3,4, i.e. swapping 3,4 and 7,8.
     !            this way, in 2D we can use half the directions using only the first
     !            4 octants (1,2,7,8) and the results are still correct.
     !
     do idire =7, 12
        direc_neu(2,idire) = - direc_neu(2,idire-6)
        direc_neu(1,idire) = direc_neu(1,idire-6)
      !   direc_neu(3,idire) =  direc_neu(3,idire-6)    
        direc_neu(3,idire) = - direc_neu(3,idire-6)    
        weigd_neu(idire)   =  weigd_neu(idire-6)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire =13, 24
           direc_neu(2,idire) =  direc_neu(2,idire-12)
           direc_neu(1,idire) =  direc_neu(1,idire-12)
           direc_neu(3,idire) = -direc_neu(3,idire-12)    
           weigd_neu(idire)   =  weigd_neu(idire-12)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire) = 2.0_rp * weigd_neu(idire)
        end do
     end if

  case ( 8_ip ) 

     !-------------------------------------------------------------------
     !
     ! S8-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.1422555_rp
     dirpr(2) = 0.5773503_rp
     dirpr(3) = 0.8040087_rp
     dirpr(4) = 0.9795543_rp
     !
     ! first octant
     !
     do idire=1, 4
        direc_neu(1,idire) = dirpr(1)   
     end do
     do idire=5, 7
        direc_neu(1,idire) = dirpr(2)   
     end do
     do idire=8,9
        direc_neu(1,idire) = dirpr(3)   
     end do
     direc_neu(1,10) = dirpr(4)   

     do idire=1, 3
        direc_neu(2,idire)    = dirpr(idire)   
        direc_neu(2,idire +4) = dirpr(idire)   
        direc_neu(2,idire +7) = dirpr(idire)   
     end do

     direc_neu(2, 4) = dirpr(4)   
     direc_neu(2,10) = dirpr(1) 

     do idire=1, 4
        direc_neu(3,idire)    = dirpr(5-idire)
        direc_neu(3,idire +3) = dirpr(5-idire)
     end do
     direc_neu(3,8 )  = dirpr(2)
     direc_neu(3,9 )  = dirpr(1)
     direc_neu(3,10)  = dirpr(1)

     weigd_neu( 1) = 0.1712359_rp
     weigd_neu( 2) = 0.0992284_rp
     weigd_neu( 3) = 0.0992284_rp
     weigd_neu( 4) = 0.1712359_rp
     weigd_neu( 5) = 0.0992284_rp
     weigd_neu( 6) = 0.4617179_rp
     weigd_neu( 7) = 0.0992284_rp
     weigd_neu( 8) = 0.0992284_rp
     weigd_neu( 9) = 0.0992284_rp
     weigd_neu(10) = 0.1712359_rp
     !
     !octant II z>0, y>0, x<0 
     !
     do idire =11, 20
        direc_neu(1,idire)= - direc_neu(1,idire-10)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-10)
        end do
        weigd_neu(idire) = weigd_neu(idire-10)
     end do
     !
     !octants III and IV z>0, y<0, x>0 and x<0
     ! ORDER CHANGE, USING 7 AND 8 HERE INSTEAD, 
     ! SEE SN4 ABOVE FOR EXPLANATION
     !
     do idire =21, 40
        direc_neu(2,idire)= - direc_neu(2,idire-20)
        direc_neu(1,idire)= direc_neu(1,idire-20)
      !   direc_neu(3,idire)=  direc_neu(3,idire-20)    
        direc_neu(3,idire)= - direc_neu(3,idire-20)    
        weigd_neu(idire) = weigd_neu(idire-20)
     end do
     !
     !octants with z<0
     !
     if( ndime == 3 ) then
        do idire = 41,80
           direc_neu(2,idire) =  direc_neu(2,idire-40)
           direc_neu(1,idire) =  direc_neu(1,idire-40)
           direc_neu(3,idire) = -direc_neu(3,idire-40)    
           weigd_neu(idire)       =  weigd_neu(idire-40)
        end do
     else if( ndime == 2 ) then  
        do idire = 1,num_directions_neu
           weigd_neu(idire) = 2*weigd_neu(idire)
        end do
     end if

  case ( 6_ip )

     !-------------------------------------------------------------------
     !  
     ! SN6: S6-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.1838670_rp
     dirpr(2) = 0.6950514_rp
     dirpr(3) = 0.9656013_rp
     !
     ! first octant
     !
     do idire=1, 3
        direc_neu(1,idire) = dirpr(1)   
     end do
     do idire=4, 5
        direc_neu(1,idire) = dirpr(2)   
     end do

     direc_neu(1,6) = dirpr(3)   

     do idire=1, 3
        direc_neu(2,idire)    = dirpr(idire)   
        direc_neu(2,idire +3) = dirpr(idire)   

     end do

     direc_neu(2, 6) = dirpr(1)   

     do idire=1, 3
        direc_neu(3,idire)    = dirpr(4-idire)   
     end do

     direc_neu(3,4) = dirpr(2)
     direc_neu(3,5) = dirpr(1)
     direc_neu(3,6) = dirpr(1)

     weigd_neu(1) = 0.1609517_rp
     weigd_neu(2) = 0.3626469_rp
     weigd_neu(3) = 0.1609517_rp
     weigd_neu(4) = 0.3626469_rp
     weigd_neu(5) = 0.3626469_rp
     weigd_neu(6) = 0.1609517_rp
     !
     ! Octant II z>0, y>0, x<0 
     !
     do idire =7, 12
        direc_neu(1,idire)= - direc_neu(1,idire-6)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-6)
        end do
        weigd_neu(idire) = weigd_neu(idire-6)
     end do
     !
     ! Octants III and IV z>0, y<0, x>0 and x<0
     ! ORDER CHANGE, USING 7 AND 8 HERE INSTEAD, 
     ! SEE SN4 ABOVE FOR EXPLANATION
     !
     do idire =13, 24
        direc_neu(2,idire)= - direc_neu(2,idire-12)
        direc_neu(1,idire)= direc_neu(1,idire-12)
      !   direc_neu(3,idire)=  direc_neu(3,idire-12)    
        direc_neu(3,idire)= - direc_neu(3,idire-12)    
        weigd_neu(idire) = weigd_neu(idire-12)
     end do
     !
     ! Octants with z<0
     !
     if( ndime == 3 ) then
        do idire =25, 48
           direc_neu(2,idire)= direc_neu(2,idire-24)
           direc_neu(1,idire)= direc_neu(1,idire-24)
           direc_neu(3,idire)= -direc_neu(3,idire-24)    
           weigd_neu(idire) = weigd_neu(idire-24)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire)=2*weigd_neu(idire)
        end do
     end if

  case ( 10_ip )   

     !-------------------------------------------------------------------
     !
     ! SN10, S10-approximation
     !
     !-------------------------------------------------------------------
     !
     !cosen directions
     !
     dirpr(1) = 0.13727193312_rp
     dirpr(2) = 0.50468891002_rp
     dirpr(3) = 0.70041288408_rp
     dirpr(4) = 0.85231773445_rp
     dirpr(5) = 0.98097544961_rp
     !
     !   first octant
     !
     do idire = 1,3
        do idime = 1,3
           direc_neu(idime,idire) = dirpr(1)   
        end do
        direc_neu(idire, idire)= dirpr(5)
     end do
     do idire = 4,6
        do idime = 1,3
           direc_neu(idime,idire) = dirpr(2)   
        end do
        direc_neu(idire-3, idire)= dirpr(3)
     end do
     do idire = 13, 15
        do idime = 1,3
           direc_neu(idime,idire) = dirpr(3)   
        end do
        direc_neu(idire-12, idire)= dirpr(1)
     end do
     do idire = 7,9
        direc_neu(idire-6, idire)   = dirpr(1)
        direc_neu(idire-6, idire+3) = dirpr(1)
     end do
     do idire = 7,8
        direc_neu(idire-5, idire)   = dirpr(2)
        direc_neu(idire-5, idire+3) = dirpr(4)
     end do
     direc_neu(1, 9)  = dirpr(2)
     direc_neu(1, 12) = dirpr(4)
     do idire = 8,9
        direc_neu(idire-7, idire)   = dirpr(4)
        direc_neu(idire-7, idire+3) = dirpr(2)
     end do
     direc_neu(3, 7)  = dirpr(4)
     direc_neu(3, 10) = dirpr(2)

     do idire =1,3
        weigd_neu(idire) = 0.0944411600_rp
     end do
     do idire =4,6
        weigd_neu(idire) = 0.1149971656_rp
     end do
     do idire =7,12
        weigd_neu(idire) = 0.1483951164_rp
     end do
     do idire =13,15
        weigd_neu(idire) = 0.0173702170_rp
     end do
     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =16, 30
        direc_neu(1,idire)= - direc_neu(1,idire-15)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-15)
        end do
        weigd_neu(idire) = weigd_neu(idire-15)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     ! ORDER CHANGE, USING 7 AND 8 HERE INSTEAD, 
     ! SEE SN4 ABOVE FOR EXPLANATION
     !
     do idire =31, 60 
        direc_neu(2,idire) = - direc_neu(2,idire-30)
        direc_neu(1,idire) = direc_neu(1,idire-30)
      !   direc_neu(3,idire) =  direc_neu(3,idire-30)    
        direc_neu(3,idire) = - direc_neu(3,idire-30)    
        weigd_neu(idire)   =  weigd_neu(idire-30)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire =61, 120
           direc_neu(2,idire) =  direc_neu(2,idire-60)
           direc_neu(1,idire) =  direc_neu(1,idire-60)
           direc_neu(3,idire) = -direc_neu(3,idire-60)    
           weigd_neu(idire)   =  weigd_neu(idire-60)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire) = 2.0_rp * weigd_neu(idire)
        end do
     end if

  case ( 12_ip )

     !-------------------------------------------------------------------
     !
     ! SN12: S12-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.12028966644554817_rp
     dirpr(2) = 0.45363844804142206_rp
     dirpr(3) = 0.6301635337190538_rp ! 0.63016353371905381_rp
     dirpr(4) = 0.7670882067384059_rp ! 0.76708820673840594_rp
     dirpr(5) = 0.8830303248501696_rp ! 0.88303032485016963_rp
     dirpr(6) = 0.98542416871763827_rp
     !
     ! first octant
     !
     do idire=1, 3
        do idime =1,3
           direc_neu(idime,idire) = dirpr(1)   
        end do
        direc_neu(idire, idire)= dirpr(6)
     end do
     do idire=4, 6
        do idime =1,3
           direc_neu(idime,idire) = dirpr(2)   
        end do
        direc_neu(idire-3, idire)= dirpr(4)
     end do
     do idire=7, 9
        do idime =1,3
           direc_neu(idime,idire) = dirpr(3)   
        end do
        direc_neu(idire-6, idire)= dirpr(2)
     end do

     do idire=10, 12
        direc_neu(idire-9, idire)   = dirpr(1)
        direc_neu(idire-9, idire+3) = dirpr(1)
     end do
     do idire=10, 11
        direc_neu(idire-8, idire)   = dirpr(2)
        direc_neu(idire-8, idire+3) = dirpr(5)
     end do
     direc_neu(1, 12)  = dirpr(2)
     direc_neu(1, 15) = dirpr(5)
     do idire=11, 12
        direc_neu(idire-10, idire)   = dirpr(5)
        direc_neu(idire-10, idire+3) = dirpr(2)
     end do
     direc_neu(3, 10)  = dirpr(5)
     direc_neu(3, 13)  = dirpr(2)

     do idire=16, 18
        direc_neu(idire-15, idire)   = dirpr(1)
        direc_neu(idire-15, idire+3) = dirpr(1)
     end do
     do idire=16, 17
        direc_neu(idire-14, idire)   = dirpr(3)
        direc_neu(idire-14, idire+3) = dirpr(4)
     end do
     direc_neu(1, 18)  = dirpr(3)
     direc_neu(1, 21) = dirpr(4)
     do idire=17, 18
        direc_neu(idire-16, idire)   = dirpr(4)
        direc_neu(idire-16, idire+3) = dirpr(3)
     end do
     direc_neu(3, 16)  = dirpr(4)
     direc_neu(3, 19)  = dirpr(3) 


     do idire =1,3
        weigd_neu(idire) = 0.0801404674_rp
     end do
     do idire =4,6
        weigd_neu(idire) = 0.1357133976_rp
     end do
     do idire =7,9
        weigd_neu(idire) = 0.0239769589_rp
     end do
     do idire =10,15
        weigd_neu(idire) = 0.0974977375_rp
     end do
     do idire =16,21
        weigd_neu(idire) = 0.0443862383_rp
     end do
     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =22, 42
        direc_neu(1,idire)= - direc_neu(1,idire-21)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-21)
        end do
        weigd_neu(idire) = weigd_neu(idire-21)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     ! ORDER CHANGE, USING 7 AND 8 HERE INSTEAD, 
     ! SEE SN4 ABOVE FOR EXPLANATION
     !
     do idire =43, 84 
        direc_neu(2,idire)= - direc_neu(2,idire-42)
        direc_neu(1,idire)= direc_neu(1,idire-42)
      !   direc_neu(3,idire)= direc_neu(3,idire-42)    
        direc_neu(3,idire)= - direc_neu(3,idire-42)    
        weigd_neu(idire) = weigd_neu(idire-42)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire = 85, 168
           direc_neu(2,idire)= direc_neu(2,idire-84)
           direc_neu(1,idire)= direc_neu(1,idire-84)
           direc_neu(3,idire)= -direc_neu(3,idire-84)    
           weigd_neu(idire) = weigd_neu(idire-84)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire)= 2.0_rp * weigd_neu(idire)
        end do
     end if


case ( 16 ) 

     !-------------------------------------------------------------------
     !
     ! SN16: S16-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !

     ! creo la cudratura desde el inicio
      ! N8=num_directions_neu/8
     if( ndime == 3 ) then
        N8=num_directions_neu/8
     else if( ndime == 2 ) then
        N8=num_directions_neu/4
     end if
    
   ! la suma de estos pesos da 1, al aplicar al octante no da pi/2  
   !   weid(1) = 0.28532856
   !   weid(2) = 0.18100812 ! 0.18199812
   !   weid(3) = 0.12803470 ! 0.12893470 
   !   weid(4) =  0.10526928 
   !   weid(5) =  0.0919211
   !   weid(6) =  0.08337745
   !   weid(7) =  0.07954919
   !   weid(8) =  0.0455116

   ! la suma de estos pesos en el octante da pi/2, pero hay que intercambiar el 3 con el 8
   ! es decir, usar el 8 seis veces (en lugar de tres veces y usar el 3 tres veces (en lugar de 6)
   !   weid(1) =0.1149    
   !   weid(2) =0.0729    
   !   weid(3) =0.0516    
   !   weid(4) =0.0424    
   !   weid(5) =0.0370    
   !   weid(6) =0.0336    
   !   weid(7) =0.0320    
   !   weid(8) =0.0183

   ! esto sale de los pesos que sumados dan 1, normalizados por 2.6411109124917767
   ! 2.6411109124917767 es lo que da al sumar los pesos en el primer octante
   ! (4.14864732) y dividirlo por pi/2
     weid(1) = 0.10803354_rp
     weid(2) = 0.06853484_rp
     weid(3) = 0.04847759_rp
     weid(4) = 0.03985796_rp
     weid(5) = 0.03480395_rp
     weid(6) = 0.03156908_rp
     weid(7) = 0.03011959_rp
     weid(8) = 0.01723199_rp

     dirpr(1) = 0.14907120_rp
     dirpr(2) = 0.39440532_rp
     dirpr(3) = 0.53748385_rp
     dirpr(4) = 0.64978629_rp
     dirpr(5) = 0.74535599_rp
     dirpr(6) = 0.82999331_rp
     dirpr(7) = 0.90676470_rp
     dirpr(8) = 0.97752522_rp
     
     !
     ! first octant
     !
     do idire=1, 3
        do idime =1,3
           direc_neu(idime,idire) = dirpr(1)   
        end do
        direc_neu(idire, idire)= dirpr(8)
     end do
     do idire=4, 6
        do idime =1,3
           direc_neu(idime,idire) = dirpr(2)   
        end do
        direc_neu(idire-3, idire)= dirpr(6)
     end do
     
     do idire=7, 9
        do idime =1,3
           direc_neu(idime,idire) = dirpr(4)   
        end do
     end do
   !   direc_neu(2, 7)= dirpr(2)
   !   direc_neu(3, 8)= dirpr(2)
   !   direc_neu(1, 9)= dirpr(2)
     direc_neu(3, 7)= dirpr(2)
     direc_neu(1, 8)= dirpr(2)
     direc_neu(2, 9)= dirpr(2)

     do idire=10, 12
        do idime =1,3
           direc_neu(idime,idire) = dirpr(3)   
        end do
     end do
     direc_neu(1, 10)= dirpr(4)
     direc_neu(2, 11)= dirpr(4)
     direc_neu(3, 12)= dirpr(4)
     
    ! (7, 2, 1) (2, 7, 1)  (1, 7, 2) (1, 2, 7)  (2, 1, 7) (7, 1, 2) 
!   idire=13, 18
    idire=13
    direc_neu(1, idire)   = dirpr(7)
    direc_neu(2, idire)   = dirpr(2)
    direc_neu(3, idire)   = dirpr(1)
    
    idire=14
    direc_neu(1, idire)   = dirpr(2)
    direc_neu(2, idire)   = dirpr(7)
    direc_neu(3, idire)   = dirpr(1)
    
    idire=15
    direc_neu(1, idire)   = dirpr(1)
    direc_neu(2, idire)   = dirpr(7)
    direc_neu(3, idire)   = dirpr(2)

    idire=16
    direc_neu(1, idire)   = dirpr(1)
    direc_neu(2, idire)   = dirpr(2)
    direc_neu(3, idire)   = dirpr(7)

    idire=17
    direc_neu(1, idire)   = dirpr(2)
    direc_neu(2, idire)   = dirpr(1)
    direc_neu(3, idire)   = dirpr(7)

    idire=18
    direc_neu(1, idire)   = dirpr(7)
    direc_neu(2, idire)   = dirpr(1)
    direc_neu(3, idire)   = dirpr(2)


!    idire=19, 24
    ! (6, 3, 1) (3, 6, 1)  (1, 6, 3) (1, 3, 6)  (3, 1, 6) (6, 1, 3) 
    idire=19
    direc_neu(1, idire)   = dirpr(6)
    direc_neu(2, idire)   = dirpr(3)
    direc_neu(3, idire)   = dirpr(1)
    
    idire=20
    direc_neu(1, idire)   = dirpr(3)
    direc_neu(2, idire)   = dirpr(6)
    direc_neu(3, idire)   = dirpr(1)
    
    idire=21
    direc_neu(1, idire)   = dirpr(1)
    direc_neu(2, idire)   = dirpr(6)
    direc_neu(3, idire)   = dirpr(3)

    idire=22
    direc_neu(1, idire)   = dirpr(1)
    direc_neu(2, idire)   = dirpr(3)
    direc_neu(3, idire)   = dirpr(6)

    idire=23
    direc_neu(1, idire)   = dirpr(3)
    direc_neu(2, idire)   = dirpr(1)
    direc_neu(3, idire)   = dirpr(6)

    idire=24
    direc_neu(1, idire)   = dirpr(6)
    direc_neu(2, idire)   = dirpr(1)
    direc_neu(3, idire)   = dirpr(3)

!    idire=25, 30
    ! (5, 4, 1) (4, 5, 1)  (1, 5, 4) (1, 4, 5)  (4, 1, 5) (5, 1, 4) 
    idire=25
    direc_neu(1, idire)   = dirpr(5)
    direc_neu(2, idire)   = dirpr(4)
    direc_neu(3, idire)   = dirpr(1)
    
    idire=26
    direc_neu(1, idire)   = dirpr(4)
    direc_neu(2, idire)   = dirpr(5)
    direc_neu(3, idire)   = dirpr(1)
    
    idire=27
    direc_neu(1, idire)   = dirpr(1)
    direc_neu(2, idire)   = dirpr(5)
    direc_neu(3, idire)   = dirpr(4)

    idire=28
    direc_neu(1, idire)   = dirpr(1)
    direc_neu(2, idire)   = dirpr(4)
    direc_neu(3, idire)   = dirpr(5)

    idire=29
    direc_neu(1, idire)   = dirpr(4)
    direc_neu(2, idire)   = dirpr(1)
    direc_neu(3, idire)   = dirpr(5)

    idire=30
    direc_neu(1, idire)   = dirpr(5)
    direc_neu(2, idire)   = dirpr(1)
    direc_neu(3, idire)   = dirpr(4)
 
!    idire=31, 36
    ! (5, 2, 3) (2, 5, 3)  (3, 5, 2) (3, 2, 5)  (2, 3, 5) (5, 3, 2) 
    idire=31
    direc_neu(1, idire)   = dirpr(5)
    direc_neu(2, idire)   = dirpr(2)
    direc_neu(3, idire)   = dirpr(3)
    
    idire=32
    direc_neu(1, idire)   = dirpr(2)
    direc_neu(2, idire)   = dirpr(5)
    direc_neu(3, idire)   = dirpr(3)
    
    idire=33
    direc_neu(1, idire)   = dirpr(3)
    direc_neu(2, idire)   = dirpr(5)
    direc_neu(3, idire)   = dirpr(2)

    idire=34
    direc_neu(1, idire)   = dirpr(3)
    direc_neu(2, idire)   = dirpr(2)
    direc_neu(3, idire)   = dirpr(5)

    idire=35
    direc_neu(1, idire)   = dirpr(2)
    direc_neu(2, idire)   = dirpr(3)
    direc_neu(3, idire)   = dirpr(5)

    idire=36
    direc_neu(1, idire)   = dirpr(5)
    direc_neu(2, idire)   = dirpr(3)
    direc_neu(3, idire)   = dirpr(2)
 
 

     do idire =1,3
        weigd_neu(idire) = weid(8)
     end do
     do idire =4,6
        weigd_neu(idire) =weid(4)
     end do
     do idire =7,9
        weigd_neu(idire) = weid(2)
     end do
     do idire =10,12
        weigd_neu(idire) = weid(1)
     end do
     do idire =13,18
        weigd_neu(idire) = weid(7)
     end do

     do idire =19,24
        weigd_neu(idire) = weid(6)
     end do
     do idire =25,30
        weigd_neu(idire) = weid(5)
     end do

     do idire =31,36
        weigd_neu(idire) = weid(3)
     end do
     

     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =N8+1, 2*N8
        direc_neu(1,idire)= - direc_neu(1,idire-N8)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-N8)
        end do
        weigd_neu(idire) = weigd_neu(idire-N8)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     ! ORDER CHANGE, USING 7 AND 8 HERE INSTEAD, 
     ! SEE SN4 ABOVE FOR EXPLANATION
     !
     do idire =2*N8+1, 4*N8 
        direc_neu(2,idire)= - direc_neu(2,idire-2*N8)
        direc_neu(1,idire)= direc_neu(1,idire-2*N8)
      !   direc_neu(3,idire)= direc_neu(3,idire-2*N8)    
        direc_neu(3,idire)= - direc_neu(3,idire-2*N8)    
        weigd_neu(idire) = weigd_neu(idire-2*N8)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire = 4*N8+1, num_directions_neu
           direc_neu(2,idire)= direc_neu(2,idire-4*N8)
           direc_neu(1,idire)= direc_neu(1,idire-4*N8)
           direc_neu(3,idire)= -direc_neu(3,idire-4*N8)    
           weigd_neu(idire) = weigd_neu(idire-4*N8)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire)= 2.0_rp * weigd_neu(idire)
        end do
     end if

   case ( 64 ) 

      N8=num_directions_neu/8
      call LDFE(N8)



  end select


end subroutine neu_snordi

subroutine neu_ordico()
  !-----------------------------------------------------------------------
  !constructs ordinates for an icosahedron
  !defines direc_neu and weigd_neu
  !
  !-----------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_parame, only : pi
  use def_domain, only : ndime
  use def_neutro
  use def_master
  implicit none

  integer(ip) :: idire, idime, ipoin, i, kfl_normx_rad!, istat
  real(rp)    :: vertp(3,12), auxve(12), vertr(3, 12)!,  auxdi(3,40)
  real(rp)    :: goldr, normv, rdumm
  external :: vecuni, neu_icodir

  kfl_normx_rad = 0
  !  
  ! construction of icosahedron
  !
  !
  ! ICOSAHEDRON VERTEXS  vertp(idime,ipoin)  ipoin =1, 12
  ! golden ratio
  !
  goldr = 0.5_rp*(1.0_rp+sqrt(5.0_rp))
  goldr = 1.0_rp/goldr
  auxve = 0.0_rp

  do ipoin=1,4
     vertp(1,ipoin)   = 0.0_rp
     vertp(3,ipoin+4) = 0.0_rp
     vertp(2,ipoin+8) = 0.0_rp
  end do
  vertp(2, 1) =  1.0_rp
  vertp(2, 2) =  1.0_rp
  vertp(2, 3) = -1.0_rp
  vertp(2, 4) = -1.0_rp

  vertp(3, 1) =  goldr
  vertp(3, 2) = -goldr
  vertp(3, 3) =  goldr
  vertp(3, 4) = -goldr

  vertp(1, 5) =  1.0_rp
  vertp(1, 6) =  1.0_rp
  vertp(1, 7) = -1.0_rp
  vertp(1, 8) = -1.0_rp

  vertp(2, 5) =  goldr
  vertp(2, 6) = -goldr
  vertp(2, 7) =  goldr
  vertp(2, 8) = -goldr

  vertp(1, 9) =  goldr
  vertp(1,10) =  goldr
  vertp(1,11) = -goldr
  vertp(1,12) = -goldr

  vertp(3, 9) =  1.0_rp
  vertp(3,10) = -1.0_rp
  vertp(3,11) =  1.0_rp
  vertp(3,12) = -1.0_rp  

  do idire =1, 12
     call vecuni(3_ip,vertp(1,idire),normv)  
  end do

  !--------------  
  ! ROTATION MATRIX: Rotation around x vector to have symmetry (casi)respect to xy plane
  !--------------
  !    |1.0           0.0                 0.0            |      
  !    |                                                 |      
  !    |0.0  gr/sqrt(gr*gr+1)    -1.0/sqrt(gr*gr+1)      |
  ! R= |                                                 |
  !    |0.0  1.0/sqrt(gr*gr+1)   gr/sqrt(gr*gr+1)        |
  !    |                                                 |
  !

  rdumm = 1.0_rp/sqrt(goldr*goldr+1.0_rp)

  vertr (1, :) = vertp(1,:)
  vertr (2, :) = rdumm*goldr*vertp(2,:) - rdumm*vertp(3,:)
  vertr (3, :) = rdumm*vertp(2,:) + goldr*rdumm*vertp(3,:)

  if( kfl_normx_rad == 1 ) then
     auxve (:)    = vertr (3,:)
     vertr (3, :) = vertr (2,:) 
     vertr (2, :) = vertr (1,:)
     vertr (1, :) = auxve (:)
  end if
  !
  !
  ! I HAVE 3 VERTEXS PER ICOSAHEDRON FACE
  !
  if( num_directions_neu == 20 ) then

     do idime = 1,3 !direction vectors pointing to the center of icosahedron faces

        direc_neu(idime, 1) = vertp(idime,9)+ vertp(idime,11)+ vertp(idime,1)
        direc_neu(idime, 2) = vertp(idime,9)+ vertp(idime,1)+ vertp(idime,5)
        direc_neu(idime, 3) = vertp(idime,1)+ vertp(idime,5)+ vertp(idime,2)
        direc_neu(idime, 4) = vertp(idime,1)+ vertp(idime,7)+ vertp(idime,2)
        direc_neu(idime, 5) = vertp(idime,1)+ vertp(idime,7)+ vertp(idime,11)
        direc_neu(idime, 6) = vertp(idime,9)+ vertp(idime,11)+ vertp(idime,3)
        direc_neu(idime, 7) = vertp(idime,9)+ vertp(idime,3)+ vertp(idime,6)
        direc_neu(idime, 8) = vertp(idime,3)+ vertp(idime,6)+ vertp(idime,4)
        direc_neu(idime, 9) = vertp(idime,3)+ vertp(idime,4)+ vertp(idime,8)
        direc_neu(idime, 10) = vertp(idime,3)+ vertp(idime,8)+ vertp(idime,11)
        direc_neu(idime, 11) = vertp(idime,10)+ vertp(idime,2)+ vertp(idime,5)
        direc_neu(idime, 12) = vertp(idime,10)+ vertp(idime,6)+ vertp(idime,5)
        direc_neu(idime, 13) = vertp(idime,10)+ vertp(idime,4)+ vertp(idime,6)
        direc_neu(idime, 14) = vertp(idime,10)+ vertp(idime,2)+ vertp(idime,12)
        direc_neu(idime, 15) = vertp(idime,10)+ vertp(idime,4)+ vertp(idime,12)
        direc_neu(idime, 16) = vertp(idime,9)+ vertp(idime,6)+ vertp(idime,5)
        direc_neu(idime, 17) = vertp(idime,11)+ vertp(idime,8)+ vertp(idime,7)
        direc_neu(idime, 18) = vertp(idime,8)+ vertp(idime,4)+ vertp(idime,12)
        direc_neu(idime, 19) = vertp(idime,12)+ vertp(idime,7)+ vertp(idime,2)
        direc_neu(idime, 20) = vertp(idime,12)+ vertp(idime,7)+ vertp(idime,8)     
     end do
     do idire =1, num_directions_neu
        call vecuni(3_ip,direc_neu(:,idire),normv) 
     end do

  else

      i =0
     call neu_icodir(vertp(:,9), vertp(:,11), vertp(:,1), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+ 3*i +5))
     i= i +1
     call neu_icodir(vertp(:,9), vertp(:,1), vertp(:,5), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,1), vertp(:,5), vertp(:,2), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,1), vertp(:,7), vertp(:,2), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,1), vertp(:,7), vertp(:,11), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,9), vertp(:,11), vertp(:,3), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,9), vertp(:,3), vertp(:,6), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,3), vertp(:,6), vertp(:,4), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,3), vertp(:,4), vertp(:,8), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,3), vertp(:,8), vertp(:,11), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,10), vertp(:,2), vertp(:,5), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,10), vertp(:,6), vertp(:,5), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,10), vertp(:,4), vertp(:,6), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,10), vertp(:,2), vertp(:,12), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,10), vertp(:,4), vertp(:,12), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,9), vertp(:,6), vertp(:,5), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,11), vertp(:,8), vertp(:,7), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,8), vertp(:,4), vertp(:,12), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,12), vertp(:,7), vertp(:,2), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
     i= i +1
     call neu_icodir(vertp(:,12), vertp(:,7), vertp(:,8), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
                        direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))


   !    i =0
   !   call neu_icodir(vertp(1,9:11), vertp(1,11:13), vertp(1,1:3), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                       direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+ 3*i +5))
   !   i= i +1
   !   call neu_icodir(vertp(1,9:11), vertp(1,1:4), vertp(1,5:8), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,1:3), vertp(1,5:7), vertp(1,2:4), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,1:3), vertp(1,7:9), vertp(1,2:4), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,1:3), vertp(1,7:9), vertp(1,11:13), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,9:11), vertp(1,11:13), vertp(1,3:5), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,9:11), vertp(1,3:5), vertp(1,6:8), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,3:5), vertp(1,6:8), vertp(1,4:6), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,3:5), vertp(1,4:6), vertp(1,8:10), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,3:5), vertp(1,8:10), vertp(1,11:13), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,10:12), vertp(1,2:4), vertp(1,5:7), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,10:12), vertp(1,6:8), vertp(1,5:7), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,10:12), vertp(1,4:6), vertp(1,6:8), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,10:12), vertp(1,2:4), vertp(1,12:14), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,10:12), vertp(1,4:6), vertp(1,12:14), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,9:11), vertp(1,6:8), vertp(1,5:7), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,11:13), vertp(1,8:10), vertp(1,7:9), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,8:10), vertp(1,4:6), vertp(1,12:14), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,12:14), vertp(1,7:9), vertp(1,2:4), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))
   !   i= i +1
   !   call neu_icodir(vertp(1,12:14), vertp(1,7:9), vertp(1,8:10), direc_neu(1,i+1:i+3), direc_neu(1,20+3*i+1:20+3*i+3), &
   !                   direc_neu(1,20+ 3*i +2:20+ 3*i +4), direc_neu(1,20 + 3*i +3:20+3*i+5))


  end if
  !
  ! STORES ONLY HALF DIRECTIONS FOR 2 DIMENSIONAL CASE (z>0) 
  !
  rdumm = 0.0_rp
  if( ndime == 2 .and. num_directions_neu /= 20 ) then !only directions with z component gt 0.0, save half directions 
     call runend('NEU_DIRECTIONS: NOT CODED')
  end if

  do idire = 1,num_directions_neu
     weigd_neu(idire) = 4-0_rp*pi/real(num_directions_neu,rp)
  end do

end subroutine neu_ordico

subroutine neu_icodir(vert1, vert2, vert3, direm, dire1, dire2, dire3)
  !-----------------------------------------------------------------------
  !subroutine that obtains 4 directions in an icosahedron face
  !vert1, vert2, vert3 are the 3 vertex coordinates
  !direm will be the central direction
  !dire1, dire2, dire3 will be the other 3        /\
  !                                              /--\
  !                                             / \/ \
  !                                            /------\
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  implicit none
  real(rp), intent(in)  :: vert1(3), vert2(3), vert3(3)
  real(rp), intent(out) :: direm(3), dire1(3), dire2(3), dire3(3)
  real(rp)              :: vert12(3), vert13(3), vert23(3)
  real(rp)              :: normv
  integer(ip)           :: idime
  external :: vecuni
  !    do idime = 1, 3
  !     direm(idime) = vert1(idime)+ vert2(idime)+ vert3(idime)
  !     vert12(idime) = vert1(idime) + vert2(idime)
  !     vert13(idime) = vert1(idime) + vert3(idime)
  !     vert23(idime) = vert2(idime) + vert3(idime)
  !    end do
  do idime = 1, 3
     direm(idime)  = vert1(idime)+ vert2(idime)+ vert3(idime)
     vert12(idime) = vert1(idime) + vert2(idime)
     vert13(idime) = vert1(idime) + vert3(idime)
     vert23(idime) = vert2(idime) + vert3(idime)
  end do

  call vecuni(3_ip,direm, normv)  
  call vecuni(3_ip,vert12,normv)  
  call vecuni(3_ip,vert13,normv)  
  call vecuni(3_ip,vert23,normv) 

  do idime = 1,3
     dire1(idime) = vert1(idime) + vert12(idime)+ vert13(idime)
     dire2(idime) = vert2(idime) + vert12(idime)+ vert23(idime)
     dire3(idime) = vert3(idime) + vert13(idime)+ vert23(idime)
  end do

  call vecuni(3_ip,dire1,normv)  
  call vecuni(3_ip,dire2,normv)  
  call vecuni(3_ip,dire3,normv)  

end subroutine neu_icodir



subroutine LDFE(N8)
 use def_kintyp, only : ip,rp
  use def_parame, only : pi
  !use def_domain, only : ndime
  use def_neutro
  use def_master
  implicit none
  integer :: N8

  integer :: idire,idime

direc_neu(1,1)= 4.38752387852336E-02_rp
direc_neu(2,1)= 4.38752387852336E-02_rp
direc_neu(3,1)= 9.98073106963150E-01_rp
weigd_neu(1) = 1.01011818772561E-02_rp
direc_neu(1,2)= 2.23194054503770E-01_rp
direc_neu(2,2)= 4.34603636088626E-02_rp
direc_neu(3,2)= 9.73804708773352E-01_rp
weigd_neu(2) = 1.47114091091454E-02_rp
direc_neu(1,3)= 4.34603636088619E-02_rp
direc_neu(2,3)= 2.23194054503772E-01_rp
direc_neu(3,3)= 9.73804708773351E-01_rp
weigd_neu(3) = 1.47114091091449E-02_rp

direc_neu(1,4)= 9.90147542976669E-02_rp
direc_neu(2,4)= 9.90147542976669E-02_rp
direc_neu(3,4)= 9.90147542976674E-01_rp
weigd_neu(4) = 1.31319081660192E-02_rp

direc_neu(1,5)= 3.91754724400026E-01_rp
direc_neu(2,5)= 5.13254348711689E-02_rp
direc_neu(3,5)= 9.18636998844235E-01_rp
weigd_neu(5) =  2.01635596440037E-02_rp

direc_neu(1,6)= 6.25662723829225E-01_rp
direc_neu(2,6)= 5.09223830015233E-02_rp
direc_neu(3,6)= 7.78429872833796E-01_rp
weigd_neu(6) = 2.42900177326475E-02_rp

direc_neu(1,7)= 4.50060753113050E-01_rp
direc_neu(2,7)= 2.86782825010254E-01_rp
direc_neu(3,7)= 8.45695530191836E-01_rp
weigd_neu(7) = 3.00096725803229E-02_rp

direc_neu(1,8)= 4.92365963917330E-01_rp
direc_neu(2,8)= 1.23091490979332E-01_rp
direc_neu(3,8)= 8.61640436855329E-01_rp
weigd_neu(8) = 2.51087542402306E-02_rp

direc_neu(1,9)= 5.13254348711655E-02_rp
direc_neu(2,9)= 3.91754724400022E-01_rp
direc_neu(3,9)= 9.18636998844237E-01_rp
weigd_neu(9) = 2.01635596440024E-02_rp

direc_neu(1,10)= 2.86782825010250E-01_rp
direc_neu(2,10)= 4.50060753113052E-01_rp
direc_neu(3,10)= 8.45695530191837E-01_rp
weigd_neu(10) = 3.00096725803234E-02_rp

direc_neu(1,11)= 5.09223830015209E-02_rp
direc_neu(2,11)= 6.25662723829230E-01_rp
direc_neu(3,11)= 7.78429872833792E-01_rp
weigd_neu(11) =  2.42900177326465E-02_rp

direc_neu(1,12)= 1.23091490979332E-01_rp
direc_neu(2,12)= 4.92365963917330E-01_rp
direc_neu(3,12)= 8.61640436855329E-01_rp
weigd_neu(12) = 2.51087542402332E-02_rp

direc_neu(1,13)= 3.35052576156500E-01_rp
direc_neu(2,13)= 3.35052576156500E-01_rp
direc_neu(3,13)= 8.80613162757510E-01_rp
weigd_neu(13) = 2.78797764691254E-02_rp

direc_neu(1,14)= 9.47958080258904E-02_rp
direc_neu(2,14)= 2.86222206910509E-01_rp
direc_neu(3,14)= 9.53462428757418E-01_rp
weigd_neu(14) =  1.90475751591700E-02_rp

direc_neu(1,15)= 2.86222206910502E-01_rp
direc_neu(2,15)= 9.47958080259107E-02_rp
direc_neu(3,15)= 9.53462428757418E-01_rp
weigd_neu(15) = 1.90475751591731E-02_rp

direc_neu(1,16)= 2.35702260395515E-01_rp
direc_neu(2,16)= 2.35702260395515E-01_rp
direc_neu(3,16)= 9.42809041582063E-01_rp
weigd_neu(16) =  2.20620660106520E-02_rp

direc_neu(1,17)= 7.78429872833780E-01_rp
direc_neu(2,17)= 5.09223830015111E-02_rp
direc_neu(3,17)= 6.25662723829246E-01_rp
weigd_neu(17) = 2.42900177326465E-02_rp

direc_neu(1,18)= 9.18636998844232E-01_rp
direc_neu(2,18)= 5.13254348711738E-02_rp
direc_neu(3,18)= 3.91754724400033E-01_rp
weigd_neu(18) =  2.01635596440023E-02_rp

direc_neu(1,19)= 8.45695530191835E-01_rp
direc_neu(2,19)= 2.86782825010258E-01_rp
direc_neu(3,19)= 4.50060753113049E-01_rp
weigd_neu(19) =  3.00096725803224E-02_rp

direc_neu(1,20)= 8.61640436855329E-01_rp
direc_neu(2,20)= 1.23091490979332E-01_rp
direc_neu(3,20)= 4.92365963917330E-01_rp
weigd_neu(20) = 2.51087542402348E-02_rp

direc_neu(1,21)= 9.73804708773349E-01_rp
direc_neu(2,21)= 4.34603636088579E-02_rp
direc_neu(3,21)= 2.23194054503781E-01_rp
weigd_neu(21) = 1.47114091091450E-02_rp

direc_neu(1,22)= 9.98073106963150E-01_rp
direc_neu(2,22)= 4.38752387852414E-02_rp
direc_neu(3,22)= 4.38752387852413E-02_rp
weigd_neu(22) =  1.01011818772558E-02_rp

direc_neu(1,23)= 9.73804708773351E-01_rp
direc_neu(2,23)= 2.23194054503774E-01_rp
direc_neu(3,23)= 4.34603636088611E-02_rp
weigd_neu(23) = 1.47114091091441E-02_rp

direc_neu(1,24)= 9.90147542976674E-01_rp
direc_neu(2,24)= 9.90147542976673E-02_rp
direc_neu(3,24)= 9.90147542976674E-02_rp
weigd_neu(24) = 1.31319081660202E-02_rp

direc_neu(1,25)= 8.45695530191832E-01_rp
direc_neu(2,25)= 4.50060753113044E-01_rp
direc_neu(3,25)= 2.86782825010275E-01_rp
weigd_neu(25) = 3.00096725803228E-02_rp

direc_neu(1,26)= 9.18636998844233E-01_rp
direc_neu(2,26)= 3.91754724400032E-01_rp
direc_neu(3,26)= 5.13254348711725E-02_rp
weigd_neu(26) =  2.01635596440023E-02_rp

direc_neu(1,27)= 7.78429872833787E-01_rp
direc_neu(2,27)= 6.25662723829236E-01_rp
direc_neu(3,27)= 5.09223830015169E-02_rp
weigd_neu(27) = 2.42900177326455E-02_rp

direc_neu(1,28)= 8.61640436855329E-01_rp
direc_neu(2,28)= 4.92365963917330E-01_rp
direc_neu(3,28)= 1.23091490979332E-01_rp
weigd_neu(28) = 2.51087542402358E-02_rp

direc_neu(1,29)= 9.53462428757418E-01_rp
direc_neu(2,29)= 2.86222206910519E-01_rp
direc_neu(3,29)= 9.47958080258602E-02_rp
weigd_neu(29) =  1.90475751591677E-02_rp

direc_neu(1,30)= 8.80613162757525E-01_rp
direc_neu(2,30)= 3.35052576156479E-01_rp
direc_neu(3,30)= 3.35052576156479E-01_rp
weigd_neu(30) = 2.78797764691281E-02_rp

direc_neu(1,31)= 9.53462428757418E-01_rp
direc_neu(2,31)= 9.47958080258961E-02_rp
direc_neu(3,31)= 2.86222206910507E-01_rp
weigd_neu(31) = 1.90475751591710E-02_rp

direc_neu(1,32)= 9.42809041582063E-01_rp
direc_neu(2,32)= 2.35702260395515E-01_rp
direc_neu(3,32)= 2.35702260395515E-01_rp
weigd_neu(32) = 2.20620660106529E-02_rp

direc_neu(1,33)= 5.09223830015081E-02_rp
direc_neu(2,33)= 7.78429872833776E-01_rp
direc_neu(3,33)= 6.25662723829251E-01_rp
weigd_neu(33) = 2.42900177326452E-02_rp

direc_neu(1,34)= 2.86782825010257E-01_rp
direc_neu(2,34)= 8.45695530191835E-01_rp
direc_neu(3,34)= 4.50060753113050E-01_rp
weigd_neu(34) =  3.00096725803224E-02_rp

direc_neu(1,35)= 5.13254348711753E-02_rp
direc_neu(2,35)= 9.18636998844231E-01_rp
direc_neu(3,35)= 3.91754724400035E-01_rp
weigd_neu(35) =  2.01635596440032E-02_rp

direc_neu(1,36)= 1.23091490979332E-01_rp
direc_neu(2,36)= 8.61640436855329E-01_rp
direc_neu(3,36)= 4.92365963917330E-01_rp
weigd_neu(36) = 2.51087542402352E-02_rp

direc_neu(1,37)= 4.50060753113045E-01_rp
direc_neu(2,37)= 8.45695530191833E-01_rp
direc_neu(3,37)= 2.86782825010272E-01_rp
weigd_neu(37) =  3.00096725803231E-02_rp

direc_neu(1,38)= 6.25662723829232E-01_rp
direc_neu(2,38)= 7.78429872833790E-01_rp
direc_neu(3,38)= 5.09223830015193E-02_rp
weigd_neu(38) =  2.42900177326471E-02_rp

direc_neu(1,39)= 3.91754724400031E-01_rp
direc_neu(2,39)= 9.18636998844233E-01_rp
direc_neu(3,39)= 5.13254348711724E-02_rp
weigd_neu(39) =  2.01635596440025E-02_rp

direc_neu(1,40)= 4.92365963917331E-01_rp
direc_neu(2,40)= 8.61640436855328E-01_rp
direc_neu(3,40)= 1.23091490979332E-01_rp
weigd_neu(40) =  2.51087542402323E-02_rp

direc_neu(1,41)= 4.34603636088583E-02_rp
direc_neu(2,41)= 9.73804708773350E-01_rp
direc_neu(3,41)= 2.23194054503780E-01_rp
weigd_neu(41) = 1.47114091091460E-02_rp

direc_neu(1,42)= 2.23194054503777E-01_rp
direc_neu(2,42)= 9.73804708773350E-01_rp
direc_neu(3,42)= 4.34603636088597E-02_rp
weigd_neu(42) = 1.47114091091451E-02_rp

direc_neu(1,43)= 4.38752387852359E-02_rp
direc_neu(2,43)= 9.98073106963150E-01_rp
direc_neu(3,43)= 4.38752387852358E-02_rp
weigd_neu(43) =  1.01011818772555E-02_rp

direc_neu(1,44)= 9.90147542976673E-02_rp
direc_neu(2,44)= 9.90147542976674E-01_rp
direc_neu(3,44)= 9.90147542976674E-02_rp
weigd_neu(44) = 1.31319081660197E-02_rp

direc_neu(1,45)= 2.86222206910522E-01_rp
direc_neu(2,45)= 9.53462428757418E-01_rp
direc_neu(3,45)= 9.47958080258519E-02_rp
weigd_neu(45) = 1.90475751591660E-02_rp

direc_neu(1,46)= 9.47958080258752E-02_rp
direc_neu(2,46)= 9.53462428757418E-01_rp
direc_neu(3,46)= 2.86222206910514E-01_rp
weigd_neu(46) = 1.90475751591669E-02_rp

direc_neu(1,47)= 3.35052576156480E-01_rp
direc_neu(2,47)= 8.80613162757525E-01_rp
direc_neu(3,47)= 3.35052576156480E-01_rp
weigd_neu(47) = 2.78797764691266E-02_rp

direc_neu(1,48)= 2.35702260395515E-01_rp
direc_neu(2,48)= 9.42809041582063E-01_rp
direc_neu(3,48)= 2.35702260395515E-01_rp
weigd_neu(48) = 2.20620660106607E-02_rp

direc_neu(1,49)= 7.02505464432133E-01_rp
direc_neu(2,49)= 7.02505464432132E-01_rp
direc_neu(3,49)= 1.13895324249885E-01_rp
weigd_neu(49) = 2.78705916500711E-02_rp

direc_neu(1,50)= 5.29915360873601E-01_rp
direc_neu(2,50)= 7.69766972496521E-01_rp
direc_neu(3,50)= 3.55877111323192E-01_rp
weigd_neu(50) = 3.51928031147286E-02_rp

direc_neu(1,51)= 7.69766972496521E-01_rp
direc_neu(2,51)= 5.29915360873602E-01_rp
direc_neu(3,51)= 3.55877111323192E-01_rp
weigd_neu(51) = 3.51928031147287E-02_rp

direc_neu(1,52)= 6.80413817439771E-01_rp
direc_neu(2,52)= 6.80413817439771E-01_rp
direc_neu(3,52)= 2.72165526975908E-01_rp
weigd_neu(52) = 3.37687299459803E-02_rp

direc_neu(1,53)= 3.55877111323199E-01_rp
direc_neu(2,53)= 7.69766972496527E-01_rp
direc_neu(3,53)= 5.29915360873588E-01_rp
weigd_neu(53) =  3.51928031147297E-02_rp

direc_neu(1,54)= 1.13895324249909E-01_rp
direc_neu(2,54)= 7.02505464432130E-01_rp
direc_neu(3,54)= 7.02505464432130E-01_rp
weigd_neu(54) =  2.78705916500730E-02_rp

direc_neu(1,55)= 3.55877111323192E-01_rp
direc_neu(2,55)= 5.29915360873602E-01_rp
direc_neu(3,55)= 7.69766972496521E-01_rp
weigd_neu(55) =  3.51928031147292E-02_rp

direc_neu(1,56)= 2.72165526975908E-01_rp
direc_neu(2,56)= 6.80413817439771E-01_rp
direc_neu(3,56)= 6.80413817439771E-01_rp
weigd_neu(56) =  3.37687299459794E-02_rp

direc_neu(1,57)= 7.69766972496528E-01_rp
direc_neu(2,57)= 3.55877111323199E-01_rp
direc_neu(3,57)= 5.29915360873587E-01_rp
weigd_neu(57) = 3.51928031147294E-02_rp

direc_neu(1,58)= 5.29915360873602E-01_rp
direc_neu(2,58)= 3.55877111323191E-01_rp
direc_neu(3,58)= 7.69766972496521E-01_rp
weigd_neu(58) = 3.51928031147287E-02_rp

direc_neu(1,59)= 7.02505464432131E-01_rp
direc_neu(2,59)= 1.13895324249906E-01_rp
direc_neu(3,59)= 7.02505464432131E-01_rp
weigd_neu(59) = 2.78705916500716E-02_rp

direc_neu(1,60)= 6.80413817439771E-01_rp
direc_neu(2,60)= 2.72165526975908E-01_rp
direc_neu(3,60)= 6.80413817439771E-01_rp
weigd_neu(60) = 3.37687299459807E-02_rp

direc_neu(1,61)= 4.85081580043393E-01_rp
direc_neu(2,61)= 4.85081580043393E-01_rp
direc_neu(3,61)= 7.27593101537672E-01_rp
weigd_neu(61) = 3.83614334577098E-02_rp

direc_neu(1,62)= 7.27593101537650E-01_rp
direc_neu(2,62)= 4.85081580043409E-01_rp
direc_neu(3,62)= 4.85081580043409E-01_rp
weigd_neu(62) =  3.83614334577117E-02_rp

direc_neu(1,63)= 4.85081580043406E-01_rp
direc_neu(2,63)= 7.27593101537654E-01_rp
direc_neu(3,63)= 4.85081580043406E-01_rp
weigd_neu(63) = 3.83614334577106E-02_rp

direc_neu(1,64)= 5.77350269189625E-01_rp
direc_neu(2,64)= 5.77350269189625E-01_rp
direc_neu(3,64)= 5.77350269189625E-01_rp
weigd_neu(64) = 4.01265145828269E-02_rp


  ! octant II z>0, y>0, x<0 
     !
     do idire =N8+1, 2*N8
        direc_neu(1,idire)= - direc_neu(1,idire-N8)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-N8)
        end do
        weigd_neu(idire) = weigd_neu(idire-N8)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     !
     do idire =2*N8+1, 4*N8 
        direc_neu(2,idire)= - direc_neu(2,idire-2*N8)
        direc_neu(1,idire)= direc_neu(1,idire-2*N8)
        direc_neu(3,idire)= direc_neu(3,idire-2*N8)    
        weigd_neu(idire) = weigd_neu(idire-2*N8)
     end do
     !
     ! octants with z<0
     !
     do idire = 4*N8+1, num_directions_neu
           direc_neu(2,idire)= direc_neu(2,idire-4*N8)
           direc_neu(1,idire)= direc_neu(1,idire-4*N8)
           direc_neu(3,idire)= -direc_neu(3,idire-4*N8)    
           weigd_neu(idire) = weigd_neu(idire-4*N8)
     end do

end subroutine LDFE


