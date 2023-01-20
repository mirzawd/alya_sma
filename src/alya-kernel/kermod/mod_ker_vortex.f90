!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_vortex.f90
!> @author  Hadrien Calmet
!> @date    2018-12-29
!> @brief   Vortex
!> @details Vortex
!-----------------------------------------------------------------------

module mod_ker_vortex
  use def_parame
  use def_master
  use def_domain
  use def_elmgeo, only : element_type
  use def_kintyp
  use def_kermod
  use mod_memory  
  implicit none

  private
  public :: ker_vortex

contains

  subroutine ker_vortex(coord_vort,marked_elements)
    implicit none
    real(rp),    intent(inout), pointer           :: coord_vort(:,:)
    logical(lg), intent(inout), pointer, optional :: marked_elements(:)
    real(rp),                   pointer           :: vered(:,:)
    real(rp)                                      :: dummr(2)

    if( INOTMASTER ) then
       !
       ! Allocate if necessary
       !
       if( .not. associated(coord_vort) ) then
          call memory_alloca(mem_modul(1:2,modul),'COORD_VORT','vortex',coord_vort,ndime,npoin)
       end if
       if( present(marked_elements) ) then
          if( .not. associated(marked_elements) ) then
             call memory_alloca(mem_modul(1:2,modul),'MARKED_ELEMENTS','vortex',marked_elements,nelem)
          end if
       end if

       nullify(vered)
       call memory_alloca(mem_modul(1:2,modul),'VERED','vortex',vered,ndime,npoin)
       call ker_redvel(vered)
       if( kfl_vortx == 1 ) then

          if( present(marked_elements) ) then
             call ker_vortex_1(vered,coord_vort,marked_elements)
          else
             call ker_vortex_1(vered,coord_vort)
          end if

       else if( kfl_vortx == 2 ) then
          call ker_vortex_2(vered,coord_vort)
       end if
       call memory_deallo(mem_modul(1:2,modul),'VERED','vortex',vered)
    else
       call ker_redvel(dummr)
    end if

  end subroutine ker_vortex

  subroutine ker_vortex_1(vered,coord_vort,marked_elements)
    !  
    !  cell based
    !
    !  linearly interpolated each components of the reduced velocity
    !  to find the center,we set vered in equationt to 0
    !  a1 + b1.r + c1.s + d1.t = 0  plan1
    !
    !  a2 + b2.r + c2.s + d2.t = 0  plan2
    !
    !  a3 + b3.r + c3.s + d3.t = 0  plan3
    !
    !
    real(rp),    intent(in),    pointer           :: vered(:,:)
    real(rp),    intent(inout), pointer           :: coord_vort(:,:) 
    logical(lg), intent(inout), pointer, optional :: marked_elements(:)
    integer(ip)                                   :: idime,ielem,counti
    integer(ip)                                   :: pnode,pelty,ipoin1,ipoin2,ipoin3,ipoin4
    integer(ip)                                   :: ierr,nsub,jsub
    real(rp)                                      :: a(3),b(3),c(3),d(3),r,s,t
    real(rp)                                      :: vn1(3),vn2(3),vn3(3)
    real(rp)                                      :: ps,mn1,mn2,para,det
    real(rp)                                      :: ma(2,2),mai(2,2),ms(2),slop
    real(rp)                                      :: coord_pt(3,3)
    real(rp)                                      :: veloc_point(3),veloc_mod
    real(rp)                                      :: qvorti_point

    nvort = 0

    if( INOTMASTER ) then   
       !
       !initialise the qvorti array
       !
       call memory_alloca(mem_modul(1:2,modul),'VORTI','ker_vorte1',vorti,ndime+1_ip,npoin)
       call vortic(1_ip)     
       !
       ! cell based 
       !
       elements:do ielem = 1,nelem
          pelty = ltype(ielem) 
          pnode = nnode(pelty)
          if (pnode==4) nsub=1
          if (pnode==6) nsub=3
          if (pnode==8) nsub=6
          if( present(marked_elements) ) marked_elements(ielem) = .false.

          subelements:do jsub=1,nsub
             if (pnode==4) then
                call tetra (ielem,ipoin1,ipoin2,ipoin3,ipoin4,ierr)
                if (ierr == 1) cycle
             else if(pnode==6) then
                call prism (ielem,ipoin1,ipoin2,ipoin3,ipoin4,ierr,jsub)
                if (ierr == 1) cycle   
             else if(pnode==8) then
                call hexa (ielem,ipoin1,ipoin2,ipoin3,ipoin4,ierr,jsub)
                if (ierr == 1) cycle
             else   
                write(*,*)"error in ker_vorte1"
             end if

             a(1) = vered(1,ipoin1)
             a(2) = vered(2,ipoin1)
             a(3) = vered(3,ipoin1)


             b(1) = - vered(1,ipoin1) + vered(1,ipoin2) 
             b(2) = - vered(2,ipoin1) + vered(2,ipoin2)
             b(3) = - vered(3,ipoin1) + vered(3,ipoin2)

             c(1) = - vered(1,ipoin1) + vered(1,ipoin3) 
             c(2) = - vered(2,ipoin1) + vered(2,ipoin3)
             c(3) = - vered(3,ipoin1) + vered(3,ipoin3)

             d(1) = - vered(1,ipoin1) + vered(1,ipoin4) 
             d(2) = - vered(2,ipoin1) + vered(2,ipoin4)
             d(3) = - vered(3,ipoin1) + vered(3,ipoin4)               
             !
             ! the line of intersection produced by the 2 plans  vn1^vn2 = vn3    
             ! 
             vn1(1)=b(1)
             vn1(2)=c(1)
             vn1(3)=d(1)

             vn2(1)=b(2)
             vn2(2)=c(2)
             vn2(3)=d(2)
             !        
             !  condition of non colinear abs(n1.n2)/||n1||*||n2||<1 if parallel =1
             !       
             ps = abs((vn1(1)*vn2(1)+vn1(2)*vn2(2)+vn1(3)*vn2(3)))       
             mn1 =sqrt( vn1(1)*vn1(1)+vn1(2)*vn1(2)+vn1(3)*vn1(3)) 
             mn2 =sqrt( vn2(1)*vn2(1)+vn2(2)*vn2(2)+vn2(3)*vn2(3))

             if (mn1==0.0_rp.or.mn2==0.0_rp)then
                cycle  
             end if

             para = ps / (mn1*mn2)  

             if (para < 1) then  
                call vecpro(vn1,vn2,vn3,ndime)
             else if (para ==1) then 
                !write(*,*)"2 plans are parallel,next element"
                !write(*,*)"element=",ielem
                cycle 
             end if
             !
             !  finding a point is belowed vn3 and the 2 plans
             !  solve a system of 2 eq 2 unknow
             !  ex r=0
             !  a1+c1.s+d1.t=0
             !  a2+c2.s+d2.t=0
             !  
             ! |c1 d1| |s| = |-a1|  ===>  A . S = R ==> S=A-1 R 
             ! |c2 d2| |t| = |-a2|
             !
             !  mai inverse matrix of A
             !
             ma(1,1) = c(1)
             ma(1,2) = d(1)
             ma(2,1) = c(2)
             ma(2,2) = d(2)

             call invmtx(ma,mai,det,ndime-1)

             ms(1)=mai(1,1)*(-a(1))+mai(1,2)*(-a(2))
             ms(2)=mai(2,1)*(-a(1))+mai(2,2)*(-a(2))
             !
             ! point P (0,s,t)
             !
             !r=0
             !s=ms(1)
             !t=ms(2)        
             !
             !solve y=ax+b => y= k vn3 + P
             !
             counti=0
             !
             ! face1 t=0 => k=-p(3)/vn3(3) 
             ! 
             !if (b(1)==c(1).and.c(1)==0) then 
             ! plan parallel to rs so 
             !           
             !check if the line pass through the face 1
             !
             if (vn3(3)==0.0_rp) then
                cycle
             else
                slop=-ms(2)/vn3(3)
                !write(*,*)"slop face 1",slop
                r= 0.0_rp + slop*vn3(1)
                s=  ms(1) + slop*vn3(2)
                t= 0.0_rp
                if (r >= 0.0_rp.and.s >= 0.0_rp.and. &
                     (s+r)<= 1.0_rp ) then
                   !write(*,*)"intersection with face 1"
                   counti=counti+1
                   do idime=1,ndime
                      coord_pt(idime,counti)=(1.0_rp-r-s-t)*coord(idime,ipoin1) &
                           +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                           +t*coord(idime,ipoin4)
                   end do
                end if
             end if
             !
             ! face2 s=0 => k=-p(2)/vn3(2) 
             !        
             if (vn3(2)==0.0_rp) then
                cycle
             else
                slop=-ms(1)/vn3(2)
                !write(*,*)"slop face 2",slop
                r= 0.0_rp + slop*vn3(1)
                s= 0.0_rp
                t= ms(2)  + slop*vn3(3)
                if (r >= 0.0_rp.and.t >= 0.0_rp.and. &
                     (r+t)<= 1.0_rp ) then
                   !write(*,*)"intersection with face 2"
                   counti=counti+1
                   do idime=1,ndime
                      coord_pt(idime,counti)=(1.0_rp-r-s-t)*coord(idime,ipoin1) &
                           +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                           +t*coord(idime,ipoin4)
                   end do
                end if
             end if
             !
             ! face3 r=0 => k=-p(1)/vn3(1) 
             !        
             if (vn3(1)==0.0_rp) then
                cycle
             else
                slop=0.0_rp/vn3(2)
                !write(*,*)"slop face 3",slop
                r= 0.0_rp 
                s= ms(1)  + slop*vn3(2)
                t= ms(2)  + slop*vn3(3)
                if (s >= 0.0_rp.and.t >= 0.0_rp.and. &
                     (s+t)<= 1.0_rp ) then
                   !write(*,*)"intersection with face 3"
                   counti=counti+1
                   do idime=1,ndime
                      coord_pt(idime,counti)=(1.0_rp-r-s-t)*coord(idime,ipoin1) &
                           +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                           +t*coord(idime,ipoin4)
                   end do
                   if (counti==3) then  
                      write(*,*) "the line pass through the 3 faces !"
                   end if
                end if
             end if
             !
             ! face4 r+s+t=1 => k=[1-p(1)-p(2)-p(3)]/[vn3(1)+vn(2)+vn(3)]
             !        
             if ((vn3(1)+vn3(2)+vn3(3))==0.0_rp) then
                cycle
             else
                slop=( 1 -s -t ) / (vn3(1)+vn3(2)+vn3(3))
                !write(*,*)"slop face 4",slop
                r= 0.0_rp + slop*vn3(1)
                s= ms(1)  + slop*vn3(2)
                t= ms(2)  + slop*vn3(3)
                if (r <= 1.0_rp.and.s <= 1.0_rp.and.  &
                     t <= 1.0_rp .and.                &
                     r >= 0.0_rp.and.s >= 0.0_rp.and. &
                     t >= 0.0_rp ) then
                   !write(*,*)"intersection with face 4"
                   counti=counti+1
                   do idime=1,ndime
                      coord_pt(idime,counti)=(1.0_rp-r-s-t)*coord(idime,ipoin1) &
                           +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                           +t*coord(idime,ipoin4)
                   end do
                   if (counti==3) then  
                      write(*,*) "the line pass through the 3 faces !"
                   else if (counti==4) then
                      write(*,*) "la concha ! algo raro pasa !"
                   end if
                end if
             end if
             !
             ! if counti>2 means that the line pass through the element
             !
             if (counti == 2) then
                if (kfl_vortx_thres == 1 ) then
                   !
                   !  first threshold with velocity criterio ( < thr_veloc )
                   !
                   do idime=1,ndime
                      veloc_point(idime)=(1.0_rp-r-s-t)*veloc(idime,ipoin1,1) &
                           +r*veloc(idime,ipoin2,1)+s*veloc(idime,ipoin3,1)&
                           +t*veloc(idime,ipoin4,1)
                   end do
                   veloc_mod=sqrt(veloc_point(1)**2+veloc_point(2)**2+veloc_point(3)**2)
                   !
                   !  second threshold with qvorticity criterio ( > thr_qvort )
                   !
                   qvorti_point=(1.0_rp-r-s-t)*vorti(ndime+1,ipoin1) &
                        +r*vorti(ndime+1,ipoin2)+s*vorti(ndime+1,ipoin3)&
                        +t*vorti(ndime+1,ipoin4)
                   !
                   !
                   !
                   if ((veloc_mod < thr_veloc ).AND.(qvorti_point > thr_qvort)) then
                      nvort=nvort+1
                      do idime=1,ndime
                         coord_vort(idime,nvort)=(coord_pt(idime,1)+coord_pt(idime,2))/2.0_rp
                      end do
                   end if
                else
                   nvort=nvort+1
                   do idime=1,ndime
                      coord_vort(idime,nvort)=(coord_pt(idime,1)+coord_pt(idime,2))/2.0_rp
                   end do
                end if
                if( present(marked_elements) ) marked_elements(ielem) = .true.
             end if
          end do subelements
       end do elements
       !
       ! Hdf5 bug in case of nvort=0, fix it with one point(0,0,0) 
       !
       if (nvort==0) then
          !write(*,*)'Hdf5 bug fixed in Mino ker_verte1'
          nvort=1
          coord_vort(1,1)=0.0_rp
          coord_vort(2,1)=0.0_rp
          coord_vort(3,1)=0.0_rp        
       end if
       !
       !
       !     
       call memory_deallo(mem_modul(1:2,modul),'VORTI','ker_vorte1',vorti)

    end if

  end subroutine ker_vortex_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ker_vortex_2(vered,coord_vort,marked_elements)
    !  
    !  Face based
    !
    !  linearly interpolated each face of the reduced velocity
    !  to find if intersection, we set vered in equationt to 0
    !
    !    a1      + b1.s + c1.t = 0
    !    a2      + b2.s + c2.t = 0
    !    a3      + b3.s + c3.t = 0
    !
    ! (1-s-t).w1 + w2.s + w3.t = 0 
    !  N1            N2     N3
    !
    implicit none
    real(rp),    intent(in),    pointer           :: vered(:,:)
    real(rp),    intent(inout), pointer           :: coord_vort(:,:)
    logical(lg), intent(inout), pointer, optional :: marked_elements(:)
    integer(ip)                                   :: idime,ielem
    integer(ip)                                   :: pnode,pelty,ipoin1,ipoin2,ipoin3
    real(rp)                                      :: a(3),b(3),c(3),s,t
    integer(ip)                                   :: iface,inode1,inode2,inode3

    nvort = 0

    if( INOTMASTER ) then
       !
       ! face based 
       !
       elements: do ielem = 1,nelem 
          pelty = ltype(ielem) 
          pnode = nnode(pelty)
          if( present(marked_elements) ) marked_elements(ielem) = .false.

          if( pnode == 4 ) then
             !
             face:do iface = 1,element_type(pelty) % number_faces 
                !
                !Shape function of triangle(shape2.f90 in kernel/mathru)
                !    
                inode1 = element_type(pelty) % list_faces(1,iface)
                inode2 = element_type(pelty) % list_faces(2,iface)
                inode3 = element_type(pelty) % list_faces(3,iface)

                ipoin1 = lnods(inode1,ielem) 
                if( lpoty(ipoin1) /= 0 ) then
                   cycle
                end if
                ipoin2 = lnods(inode2,ielem) 
                if( lpoty(ipoin2) /= 0 ) then
                   cycle
                end if
                ipoin3 = lnods(inode3,ielem) 
                if( lpoty(ipoin3) /= 0 ) then
                   cycle
                end if

                a(1) = vered(1,ipoin1)
                a(2) = vered(2,ipoin1)
                a(3) = vered(3,ipoin1)          

                b(1) = -vered(1,ipoin1) + vered(1,ipoin2) 
                b(2) = -vered(2,ipoin1) + vered(2,ipoin2)
                b(3) = -vered(3,ipoin1) + vered(3,ipoin2)

                c(1) = -vered(1,ipoin1) + vered(1,ipoin3) 
                c(2) = -vered(2,ipoin1) + vered(2,ipoin3)
                c(3) = -vered(3,ipoin1) + vered(3,ipoin3)
                !
                ! solve Vr(x,y,z) = 0 (Vr(x)=0 and Vr(y)=0) ===> find r,s  
                !
                if( abs(c(1)) > zeror ) then
                   s =   ( a(2) - ( (a(1)*c(2)) / c(1) ) ) / ( -b(2) + ( (b(1)*c(2) ) / c(1) ) )
                   t = - ( a(1) + b(1) * s ) / c(1)

                   if (s >= 0.0_rp .and. s <= 1.0_rp .and. t >= 0.0_rp .and. t <= 1.0_rp) then  
                      if (kfl_vortx_thres == 1 ) then
                         write(*,*)"face threshold not implemented"
                         return
                      end if
                      !
                      !
                      !
                      nvort = nvort + 1
                      do idime = 1,ndime
                         coord_vort(idime,nvort) = &
                              (1.0_rp-s-t)*coord(idime,ipoin1) &
                              +s*coord(idime,ipoin2)+t*coord(idime,ipoin3)                       
                      end do
                      if( present(marked_elements) ) marked_elements(ielem) = .true.

                   end if
                end if

             end do face
             !
             !
          else   
             write(*,*)"error in face_based ker_vorte2"
             stop
          end if

       end do elements
       !
       ! Hdf5 bug in case of nvort=0, fix it with one point(0,0,0) 
       !
       if (nvort==0) then
          !write(*,*)'Hdf5 bug fixed in Mino ker_verte1'
          nvort=1
          coord_vort(1,1)=0.0_rp
          coord_vort(2,1)=0.0_rp
          coord_vort(3,1)=0.0_rp        
       end if
       !
       !
       !     

    end if

  end subroutine ker_vortex_2

  subroutine tetra (ielem,ipoin1,ipoin2,ipoin3,ipoin4,ierr)
    integer(ip),   intent(in)  :: ielem
    integer(ip),   intent(out) :: ipoin1,ipoin2,ipoin3,ipoin4,ierr

    ierr=0

    ipoin1 = lnods(1,ielem)
    if (lpoty(ipoin1)/= 0 ) then
       ierr=1
    end if
    ipoin2 = lnods(2,ielem)
    if (lpoty(ipoin2)/= 0 ) then
       ierr=1
    end if
    ipoin3 = lnods(3,ielem)
    if (lpoty(ipoin3)/= 0 ) then
       ierr=1 
    end if
    ipoin4 = lnods(4,ielem)
    if (lpoty(ipoin4)/= 0 ) then
       ierr=1 
    end if


  end subroutine tetra

  subroutine prism (ielem,ipoin1,ipoin2,ipoin3,ipoin4,ierr,jsub)
    implicit none
    integer(ip),   intent(in)  :: ielem,jsub
    integer(ip),   intent(out) :: ipoin1,ipoin2,ipoin3,ipoin4,ierr

    ierr=0

    if (jsub==1) then ! sub-element 1 -> 1245
       ipoin1 = lnods(1,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(2,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(4,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(5,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==2) then ! sub-element 2 -> 2345 
       ipoin1 = lnods(2,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(3,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(4,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(5,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==3) then ! sub-element 3 -> 3456 
       ipoin1 = lnods(3,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(4,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(5,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(6,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    end if

  end subroutine prism


  subroutine hexa (ielem,ipoin1,ipoin2,ipoin3,ipoin4,ierr,jsub)
    integer(ip),   intent(in)  :: ielem,jsub
    integer(ip),   intent(out) :: ipoin1,ipoin2,ipoin3,ipoin4,ierr

    ierr=0
    if (jsub==1) then ! sub-element 1 -> 1245
       ipoin1 = lnods(1,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(2,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(4,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(5,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==2) then ! sub-element 2 -> 2345
       ipoin1 = lnods(2,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(3,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(4,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(5,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==3) then ! sub-element 3 -> 3458
       ipoin1 = lnods(3,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(4,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(5,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(8,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==4) then ! sub-element 4 -> 2567
       ipoin1 = lnods(2,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(5,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(6,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(7,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==5) then ! sub-element 5 -> 2357
       ipoin1 = lnods(2,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(3,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(5,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(7,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    else if (jsub==6) then ! sub-element 6 -> 3578
       ipoin1 = lnods(3,ielem)
       if (lpoty(ipoin1)/= 0 ) then
          ierr=1
       end if
       ipoin2 = lnods(5,ielem)
       if (lpoty(ipoin2)/= 0 ) then
          ierr=1
       end if
       ipoin3 = lnods(7,ielem)
       if (lpoty(ipoin3)/= 0 ) then
          ierr=1 
       end if
       ipoin4 = lnods(8,ielem)
       if (lpoty(ipoin4)/= 0 ) then
          ierr=1 
       end if
    end if

  end subroutine hexa

end module mod_ker_vortex
!> @}
