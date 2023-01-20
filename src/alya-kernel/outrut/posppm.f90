!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine posppm()
  !-----------------------------------------------------------------------
  !****f* master/Endste
  ! NAME
  !    Endste
  ! DESCRIPTION
  !    This routine closes a time step.
  ! USES
  !    Nastin
  !    Pressr
  !    Codire
  !    Alefor
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_iofile,         only : iofile
  use mod_elsest,         only : elsest_host_element
  implicit none
  integer(ip)                :: ix,iy,i,ipoin,ielem,inode
  integer(ip)                :: imgwidth,imgheight
  integer(ip), save          :: ipass=0
  real(rp)                   :: rmgwidth,rmgheight
  real(rp),    save          :: xplan_min(3)
  real(rp),    save          :: xplan_max(3)
  real(rp)                   :: a,b,c,d,x0,y0,z0
  real(rp)                   :: xx(3)
  real(rp)                   :: xmin,xmax,xvalu
  real(rp)                   :: shapp(mnode),deriv(mnode*3),coloc(3)
  real(rp),    pointer       :: imgpixels(:,:)
  real(rp),    pointer       :: shape_pixels(:,:,:)
  integer(ip), pointer       :: active_pixels(:,:)
  integer(ip), pointer, save :: element_pixels(:,:)
  integer(ip), pointer       :: marked_nodes(:)
  real(rp)                   :: rgb(3),dista
  character(150)             :: fil_pixel
  character(20)              :: wttim

  if( kfl_pixel > 0 ) then 
     if( mod(ittim,kfl_pixel) == 0 ) then

        nullify(imgpixels)
        if( ipass == 0 ) then
           nullify(marked_nodes)
           nullify(element_pixels)
           nullify(active_pixels)
           nullify(shape_pixels)
        end if

        imgwidth  = number_pixel(1)
        imgheight = number_pixel(2)
        call memory_alloca(memor_dom,'IMGPIXELS','posppm',imgpixels,imgwidth,imgheight)

        !----------------------------------------------------------------------
        !
        ! Construct pixelization of the plane
        !
        !----------------------------------------------------------------------

        if( ipass == 0 ) then
           if( INOTMASTER ) then
              if( ndime == 3 ) then
                 a     =  0.0_rp
                 b     =  0.0_rp
                 c     =  0.0_rp
                 d     = -coord_plane_pixel
                 if( plane_pixel == 1 ) then 
                    a = 1.0_rp
                 else if( plane_pixel == 2 ) then 
                    b = 1.0_rp
                 else 
                    c = 1.0_rp
                 end if
                 !
                 ! Mark the nodes of the elements crossed by the plane
                 !
                 call memory_alloca(memor_dom,'MARKED_NODES','posppm',marked_nodes,npoin)
                 call marpla(a,b,c,d,marked_nodes)
                 xplan_min =  huge(1.0_rp)  
                 xplan_max = -huge(1.0_rp)  
                 !
                 ! Bounding box of the plane
                 !
                 do ipoin = 1,npoin
                    if( marked_nodes(ipoin) > 0 ) then
                       xplan_min(1) = min( xplan_min(1) , coord(1,ipoin) )
                       xplan_min(2) = min( xplan_min(2) , coord(2,ipoin) )
                       xplan_min(3) = min( xplan_min(3) , coord(3,ipoin) )
                       xplan_max(1) = max( xplan_max(1) , coord(1,ipoin) )
                       xplan_max(2) = max( xplan_max(2) , coord(2,ipoin) )
                       xplan_max(3) = max( xplan_max(3) , coord(3,ipoin) )
                    end if
                 end do
                 call memory_deallo(memor_dom,'MARKED_NODES','posppm',marked_nodes)
              else
                 do ipoin = 1,npoin
                    xplan_min(1) = min( xplan_min(1) , coord(1,ipoin) )
                    xplan_min(2) = min( xplan_min(2) , coord(2,ipoin) )
                    xplan_max(1) = max( xplan_max(1) , coord(1,ipoin) )
                    xplan_max(2) = max( xplan_max(2) , coord(2,ipoin) )
                 end do
                 plane_pixel = 3
                 coord_plane_pixel = 0.0_rp
              end if
           end if
           call PAR_MIN(ndime,xplan_min,'IN MY CODE')
           call PAR_MAX(ndime,xplan_max,'IN MY CODE')
        end if

        !----------------------------------------------------------------------
        !
        ! Determine who's in charge of each pixel
        !
        !----------------------------------------------------------------------

        if( ipass == 0 ) then
           call memory_alloca(memor_dom,'ACTIVE_PIXELS' ,'posppm',active_pixels, imgwidth,imgheight) 
           if( INOTMASTER ) then
              call memory_alloca(memor_dom,'ELEMENT_PIXELS','posppm',element_pixels,imgwidth,imgheight) 
              call memory_alloca(memor_dom,'ELEMENT_PIXELS','posppm',shape_pixels,imgwidth,imgheight,mnode) 
              rmgwidth  =  1.0_rp / real(imgwidth,rp)
              rmgheight =  1.0_rp / real(imgheight,rp)
              x0        =  xplan_max(1) - xplan_min(1)
              y0        =  xplan_max(2) - xplan_min(2)
              z0        =  xplan_max(3) - xplan_min(3)
              do iy = 1,imgheight
                 do ix = 1,imgwidth
                    if( plane_pixel == 1 ) then
                       xx(1) = coord_plane_pixel
                       xx(2) = xplan_min(2) + real(ix,rp) * rmgwidth  * y0
                       xx(3) = xplan_min(3) + real(iy,rp) * rmgheight * z0
                    else if( plane_pixel == 2 ) then
                       xx(1) = xplan_min(1) + real(ix,rp) * rmgwidth  * x0
                       xx(2) = coord_plane_pixel 
                       xx(3) = xplan_min(3) + real(iy,rp) * rmgheight * z0
                    else 
                       xx(1) = xplan_min(1) + real(ix,rp) * rmgwidth  * x0
                       xx(2) = xplan_min(2) + real(iy,rp) * rmgheight * y0
                       xx(3) = coord_plane_pixel 
                    end if
                    call elsest_host_element(&
                         ielse,relse,1_ip,meshe(ndivi),xx,ielem,&
                         shapp,deriv,coloc,dista)
                    !call elsest(&
                    !     2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
                    !     lnods,ltype,ltopo,coord,xx,relse,&
                    !     ielem,shapp,deriv,coloc,dummi)
                    if( ielem > 0 ) then
                       active_pixels(ix,iy)        = kfl_paral
                       element_pixels(ix,iy)       = ielem
                       shape_pixels(ix,iy,1:mnode) = shapp(1:mnode)
                    end if
                 end do
              end do
           end if
           call PAR_MAX(active_pixels,'IN MY CODE')
           if( INOTMASTER ) then
              do iy = 1,imgheight
                 do ix = 1,imgwidth
                    if( kfl_paral /= active_pixels(ix,iy) ) element_pixels(ix,iy) = 0
                 end do
              end do
              call memory_deallo(memor_dom,'ACTIVE_PIXELS' ,'posppm',active_pixels) 
           end if
        end if

        !----------------------------------------------------------------------
        !
        ! Interpolate xvalu, min and max
        !
        !----------------------------------------------------------------------

        if( INOTMASTER ) then
           xmin      =  huge(1.0_rp)
           xmax      = -huge(1.0_rp)
           do iy = 1,imgheight
              do ix = 1,imgwidth
                 ielem          = element_pixels(ix,iy) 
                 shapp(1:mnode) = shape_pixels(ix,iy,1:mnode) 
                 xvalu          = 0.0_rp
                 if( ielem > 0 ) then

                    if( variable_pixel == ID_PRESS ) then                      ! Pressure
                       do inode = 1,nnode(ltype(ielem))
                          ipoin = lnods(inode,ielem)
                          xvalu = xvalu + shapp(inode) * press(ipoin,1)
                       end do

                    else if( variable_pixel == ID_VELOC ) then                 ! Velocity
                       xx = 0.0_rp
                       do inode = 1,nnode(ltype(ielem))
                          ipoin = lnods(inode,ielem)
                          xx(1:ndime) = xx(1:ndime) + shapp(inode) * veloc(1:ndime,ipoin,1)
                       end do
                       xvalu = sqrt( xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3) )

                    else if( variable_pixel == ID_TEMPER ) then                ! Temper
                       do inode = 1,nnode(ltype(ielem))
                          ipoin = lnods(inode,ielem)
                          xvalu = xvalu + shapp(inode) * tempe(ipoin,1)
                       end do
                    end if

                    xmin = min(xmin,xvalu)
                    xmax = max(xmax,xvalu)
                 end if

                 imgpixels(ix,iy) = xvalu

              end do
           end do
        end if

        !----------------------------------------------------------------------
        !
        ! Get values at all pixels as well as min and max
        !
        !----------------------------------------------------------------------

        call PAR_MIN(xmin,     'IN MY CODE')
        call PAR_MAX(xmax,     'IN MY CODE')
        call PAR_SUM(imgpixels,'IN MY CODE')

        !----------------------------------------------------------------------
        !
        ! Write file
        !
        !----------------------------------------------------------------------

        if( INOTSLAVE ) then 

           if( ittim < 10 ) then
              write(wttim,'(a,i1)') '0000000',ittim
           else if( ittim < 100 ) then
              write(wttim,'(a,i2)') '000000',ittim
           else if( ittim < 1000 ) then
              write(wttim,'(a,i3)') '00000',ittim
           else if( ittim < 10000 ) then
              write(wttim,'(a,i4)') '0000',ittim
           else if( ittim < 100000 ) then
              write(wttim,'(a,i5)') '000',ittim
           else if( ittim < 1000000 ) then
              write(wttim,'(a,i6)') '00',ittim
           else if( ittim < 10000000 ) then
              write(wttim,'(a,i7)') '0',ittim
           else
              write(wttim,'(i8)') ittim         
           end if

           fil_pixel = adjustl(trim(namda))//'-'//trim(wttim)//'.ppm'  
           call iofile(0_ip,lun_tempo,fil_pixel,'AUTOMATIC DETECTION','unknown','formatted')

           open(unit=lun_tempo,file=trim(fil_pixel),form='formatted',status='unknown')

           write(lun_tempo, '(A)') 'P3'
           write(lun_tempo, '(I4)') int(imgwidth,4)
           write(lun_tempo, '(I4)') int(imgheight,4)
           write(lun_tempo, *) 255_4

           do  iy = imgheight,1,-1
              do  ix = 1,imgwidth
                 !
                 ! Assign color
                 !
                 if( imgpixels(ix,iy) > 0.0_rp ) then
                    xvalu = imgpixels(ix,iy)
                    call varrgb(xvalu,xmin,xmax,rgb(1),rgb(2),rgb(3))
                 else
                    rgb = 0.0_rp
                 end if
                 ! 
                 ! Write RGB
                 !
                 write(lun_tempo,100,advance='no') ( int(rgb(i)*255.0_rp,4_4),i=1,3 )
              end do
              write(lun_tempo,110)
           end do

           close(unit=lun_tempo)

        end if

        call memory_deallo(memor_dom,'IMGPIXELS','posppm',imgpixels)

     end if
  end if

100 FORMAT(3(I4))
110 FORMAT((/))

end subroutine posppm

subroutine varrg2(xvalu,xrange,mid,red,green,blue)
  use def_kintyp, only : rp
  implicit none
  real(rp), intent(in)  :: xvalu,xrange,mid
  real(rp), intent(out) :: red,green,blue
  real(rp)              :: s
!  real(rp)              :: h,l,v

  red = (xvalu-mid)/xrange

  if( red > 0.0_rp ) then
     !above mid = red-green
     blue  = 0.0_rp
     green = 1.0_rp-red
  else 
     !lower half green-blue
     blue  = -red
     green =  1.0_rp-blue
     red   =  0.0_rp
  end if
  !
  ! Sepia conversion
  !
  !red   = ( red * 0.393 + green * 0.769 + blue * 0.189 ) / 1.351;
  !green = ( red * 0.349 + green * 0.686 + blue * 0.168 ) / 1.203;
  !blue  = ( red * 0.272 + green * 0.534 + blue * 0.131 ) / 2.140;
  !
  ! Gray scale
  !
  !red   = red * 0.3
  !green = green * 0.59
  !blue  = blue * 0.11
  !red   = red + green + blue
  !green = red
  !blue  = red
  !
  ! Lighten
  !
  !red   = min(1.0_rp,red   + 25.0_rp/255.0_rp)
  !green = min(1.0_rp,green + 25.0_rp/255.0_rp)
  !blue  = min(1.0_rp,blue  + 25.0_rp/255.0_rp)
  !
  ! Grey-violet
  !
  !red = 0.2125 * red + 0.7154 * green + 0.0721 * blue
  !green = 0.2125 * red + 0.7154 * green + 0.0721 * blue
  !blue = 0.2125 * red + 0.7154 * green + 0.0721 * blue
  !
  ! Desaturate
  !
  !call rgbhls(red*255.0_rp,green*255.0_rp,blue*255.0_rp,h,l,s)
  !s=s/10.0_rp
  !l=l*1.005_rp
  !s=max(0.0_rp,s)
  !s=min(255.0_rp,s*2.5_rp)
  !call hlsrgb(h,l,s,red,green,blue)
  !
  ! Saturate
  !
  !!!call rgbhvs(red*255.0_rp,green*255.0_rp,blue*255.0_rp,h,v,s) ! Put back
  s = 0.0_rp !TODO: Needs correciton, s was left uninitialized
  s=min(255.0_rp,s*2.0_rp)
  !h = min(255.0_rp,12.5 * h) 
  !!!call hvsrgb(h,v,s,red,green,blue) ! Put back
  red   = red   / 255.0_rp
  green = green / 255.0_rp
  blue  = blue  / 255.0_rp

  red   = max(0.0_rp,red)
  red   = min(1.0_rp,red)
  green = max(0.0_rp,green)
  green = min(1.0_rp,green)
  blue  = max(0.0_rp,blue)
  blue  = min(1.0_rp,blue)

end subroutine varrg2

subroutine varrgb(xvalu,xmin,xmax,red,green,blue)
  !
  ! blue-green-yellow-red-magenta
  !
  use def_kintyp, only : rp
  implicit none
  real(rp), intent(in)  :: xvalu,xmin,xmax
  real(rp), intent(out) :: red,green,blue
  real(rp)              :: s

  s = xmin + (xvalu-xmin)/(xmax-xmin)
  s = max(0.0_rp,s)
  s = min(1.0_rp,s)

  if( s < 0.2_rp ) then          ! 1. blue1, green+
     blue  = 1.0_rp
     green = s/0.2_rp
     red   = 0.0_rp
  else if( s < 0.4_rp ) then     ! 2. green1, blue-
     blue  = (0.4_rp-s)/0.2_rp
     green = 1.0_rp
     red   = 0.0_rp
  else if( s < 0.6_rp ) then     ! 3. green1, red+
     blue  = 0.0_rp
     green = 1.0_rp
     red   = (s-0.4_rp)/0.2_rp
  else if( s < 0.8_rp ) then     ! 4. green-, red1
     blue  = 0.0_rp
     green = (0.8_rp-s)/0.2_rp
     red   = 1.0_rp
  else                           ! 5. red, blue+
     blue  = (s-0.8_rp)/0.2_rp
     green = 0.0_rp
     red   = 1.0_rp
  end if

end subroutine varrgb

subroutine marpla(a,b,c,d,marked_nodes)
  !-----------------------------------------------------------------------
  !****f* Domain/marpla
  ! NAME
  !    marpla
  ! DESCRIPTION
  !    Inrestect a plane with an element
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memchk
  use def_parame
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  real(rp),    intent(in)  :: a,b,c,d
  integer(ip), intent(out) :: marked_nodes(*)
  integer(ip)              :: ipoin,jpoin,izdom,ifoun
  real(rp)                 :: xp1(ndime),xp2(ndime)
  real(rp)                 :: t

  if( INOTMASTER ) then
     if( ndime == 3 ) then
        !
        ! Equation of the plane: a x + b y + c z + d = 0
        !
        do ipoin = 1,npoin
           izdom = r_dom(ipoin)-1
           do while( izdom < r_dom(ipoin+1)-1 )

              izdom = izdom + 1
              jpoin = c_dom(izdom)
              ifoun = 0

              if( jpoin < ipoin ) then

                 xp1(1) = coord(1,ipoin)
                 xp1(2) = coord(2,ipoin)
                 xp1(3) = coord(3,ipoin)          
                 xp2(1) = coord(1,jpoin)
                 xp2(2) = coord(2,jpoin)
                 xp2(3) = coord(3,jpoin)   
                 !
                 ! Scalar product t = n.(P1,P2)
                 !
                 t = a*(xp2(1)-xp1(1)) + b*(xp2(2)-xp1(2)) + c*(xp2(3)-xp1(3))  

                 if( t /= 0.0_rp ) then
                    !
                    ! Compute parametric coordinate t on P1-P2
                    !
                    t = (- d - a*xp1(1) - b*xp1(2) - c*xp1(3) ) / t
                    if( t >= 0.0_rp .and. t <= 1.0_rp ) ifoun = 1

                 else if( abs( a * xp1(1) + b * xp1(2) + c * xp1(3) + d ) == 0.0_rp ) then
                    !
                    ! (P1,P2) is parallel to plane: check if P1 on plane
                    !
                    ifoun = 1
                 end if

                 if( ifoun == 1 ) then
                    !
                    ! Mark IPOIN and JPOIN
                    !
                    marked_nodes(ipoin) = 1
                    marked_nodes(jpoin) = 1
                    izdom        = r_dom(ipoin+1)-1                    
                 end if

              end if
           end do
        end do

     end if
     !
     ! Exchange in Parallel
     !
     call PAR_INTERFACE_NODE_EXCHANGE(1_ip,marked_nodes,'SUM')
     do ipoin = 1,npoin
        marked_nodes(ipoin) = min(1_ip,marked_nodes(ipoin))
     end do

  end if

end subroutine marpla
