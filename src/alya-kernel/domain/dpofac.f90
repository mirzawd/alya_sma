!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine dpofac(&
     pnodb,xcoor,bocod,fabox,norma,dista,factor,proje)
  !-----------------------------------------------------------------------
  ! NAME
  !    pofadi
  ! DESCRIPTION
  !    Minimun distance between a point and a face
  !    point:  point
  !    bocod:  triangle coordinates
  !    norma:  nornmalized exterior normal to the triangle
  !    dista:  minimum distance between the point and the triangle
  !    factor: inside or outside
  !    proje:  nearest triangle point projection  
  ! USED BY
  !    shdiib
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only           :  ip,rp
  use def_domain, only           :  ndime
  use def_master, only           :  zeror
  
  implicit   none
  integer(ip),   intent(in)      :: pnodb
  real(rp),      intent(in)      :: xcoor(ndime),norma(ndime)
  real(rp),      intent(in)      :: bocod(ndime,*)
  real(rp),      intent(in)      :: fabox(2,ndime)
  real(rp),      intent(out)     :: dista
  real(rp),      intent(out)     :: factor
  real(rp),      intent(out)     :: proje(ndime)
  integer(ip)                    :: idime,ifoun,ntria
  real(rp)                       :: pladi,plapo(ndime)
  real(rp)                       :: dist1,dist2,dist3,dist4
  real(rp)                       :: proj1(3),proj2(3),proj3(3),proj4(3)
  real(rp)                       :: bari1,bari2
  !
  ! Distance between the point and the plane formed by the face
  ! Dot product between the exterior normal and the vector build from the first triangle vertex
  ! to the point  
  !
  pladi = 0.0_rp
  do idime = 1,ndime
     pladi = pladi + norma(idime) * ( xcoor(idime) - bocod(idime,1) )
  end do
  !
  ! Point projection on the plane
  !  
  plapo(1)     = xcoor(1)     - pladi*norma(1)
  plapo(2)     = xcoor(2)     - pladi*norma(2)
  plapo(ndime) = xcoor(ndime) - pladi*norma(ndime)

  if ( pladi < 0.0_rp ) then
     factor = -1.0_rp
  else
     factor =  1.0_rp
  end if
  pladi = pladi * pladi
  ifoun = 0
 
  if( ndime == 3 ) then
     !
     ! Check if we are in bounding box
     !
     if( plapo(1) >= fabox(1,1) .and. plapo(1) <= fabox(2,1) ) then
        if( plapo(2) >= fabox(1,2) .and. plapo(2) <= fabox(2,2) ) then
           if( plapo(3) >= fabox(1,3) .and. plapo(3) <= fabox(2,3) ) then
              !
              ! Determine if the projection point on plane is inside the triangle
              !
              call instr2(pnodb,plapo,bocod,ifoun,bari1,bari2,ntria)
           end if
        end if
     end if
     !
     ! The projection point is inside the triangle
     !
     if( ifoun == 1 ) then                                
        dista    = pladi        
        proje(1) = plapo(1)
        proje(2) = plapo(2)
        proje(3) = plapo(3)
     else
        if( pnodb == 4 ) then
           call dposeg(xcoor,bocod(1,1),bocod(1,2),dist1,proj1,ifoun)
           call dposeg(xcoor,bocod(1,2),bocod(1,4),dist2,proj2,ifoun)
           call dposeg(xcoor,bocod(1,3),bocod(1,4),dist3,proj3,ifoun)         
           call dposeg(xcoor,bocod(1,1),bocod(1,4),dist4,proj4,ifoun)         

           if( dist1 < dist2 .and. dist1 < dist3 .and. dist1 < dist4 ) then              
              dista = dist1 
              proje(1) = proj1(1)
              proje(2) = proj1(2)
              proje(3) = proj1(3)
           else if( dist2 < dist1 .and. dist2 < dist3 .and. dist2 < dist4 ) then        
              dista = dist2
              proje(1) = proj2(1)
              proje(2) = proj2(2)
              proje(3) = proj2(3)
           else if( dist3 < dist1 .and. dist3 < dist2 .and. dist3 < dist4 ) then        
              dista = dist3 
              proje(1) = proj3(1)
              proje(2) = proj3(2)
              proje(3) = proj3(3)              
           else
              dista = dist4
              proje(1) = proj4(1)
              proje(2) = proj4(2)
              proje(3) = proj4(3)
           end if

        else
           call dposeg(xcoor,bocod(1,1),bocod(1,2),dist1,proj1,ifoun)
           call dposeg(xcoor,bocod(1,2),bocod(1,3),dist2,proj2,ifoun)
           call dposeg(xcoor,bocod(1,1),bocod(1,3),dist3,proj3,ifoun)           
           if( dist1 < dist2 .and. dist1 < dist3 ) then              
              dista = dist1
              proje(1) = proj1(1)
              proje(2) = proj1(2)
              proje(3) = proj1(3)
           else if( dist2 < dist1 .and. dist2 < dist3 ) then
              dista = dist2 
              proje(1) = proj2(1)
              proje(2) = proj2(2)
              proje(3) = proj2(3)
           else
              dista = dist3 
              proje(1) = proj3(1)
              proje(2) = proj3(2)
              proje(3) = proj3(3)
           end if
        end if
     end if

  else if( ndime == 2 ) then

     call dposeg(xcoor,bocod(1,1),bocod(1,2),dist1,proje,ifoun)           
     dista = dist1 

  end if

end subroutine dpofac
