!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    elebin.f90
!> @author  Guillaume Houzeaux
!> @date    04/04/2014
!> @brief   Bin structure for elements neighbors
!> @details Compute a bin structure for elements neighbors
!>          If the mesh changes, this subroutine is called back
!>          by domarr
!> @}
!----------------------------------------------------------------------

subroutine elebin()
  use def_kintyp
  use def_master
  use def_domain
  use mod_memory
  use def_kermod
  use mod_messages, only : livinf
  implicit none
  integer(ip)                  :: ielel,melel,inode,ipoin
  integer(ip)                  :: ielem,jelem,nx,ny,nz,kelem
  integer(ip)                  :: imin(3),imax(3)
  integer(ip)                  :: ii,jj,kk,nn
!  integer(ip)                  :: ll
  real(rp)                     :: dn(3)
  real(rp)                     :: delta(3)
  real(rp)                     :: comin(3)
  real(rp)                     :: comax(3)
  real(rp),    pointer         :: xboun(:,:,:)
!  integer(ip), pointer         :: lnods_tmp(:,:)

  if( kfl_element_bin /= 0 ) then

     call livinf(0_ip,'COMPUTE BIN STRUCTURE FOR ELEMENT NEIGHBORS',0_ip)

     if( INOTMASTER ) then
        !
        ! Size of the bin in each direction
        !
        element_bin_boxes(1:3) = 4
        element_bin_boxes(1:3) = max(1_ip,element_bin_boxes(1:3))
        !
        ! Deallocate if necessary
        !
        if( associated(element_bin) ) then
           do ielem = 1,nelem
              nullify(element_bin(ielem) % bin_size)
              nullify(element_bin(ielem) % list_elements)
              call memory_deallo(memor_dom,'ELEMENT_BIN'  ,'elebin',element_bin(ielem) % bin_size)
              call memory_deallo(memor_dom,'LIST_ELEMENTS','elebin',element_bin(ielem) % list_elements)
           end do
           deallocate(element_bin)
           nullify(element_bin)
        end if
        !
        ! Nullify and allocate
        !
        nullify(xboun)
        allocate(element_bin(nelem))
        do ielem = 1,nelem
           nullify(element_bin(ielem) % bin_size)
           nullify(element_bin(ielem) % list_elements)
           element_bin(ielem) % boxes(1) = element_bin_boxes(1)
           element_bin(ielem) % boxes(2) = element_bin_boxes(2)
           element_bin(ielem) % boxes(3) = element_bin_boxes(3)
           nx                            = element_bin(ielem) % boxes(1)
           ny                            = element_bin(ielem) % boxes(2)
           nz                            = element_bin(ielem) % boxes(3)
           call memory_alloca(memor_dom,'ELEMENT_BIN'  ,'elebin',element_bin(ielem) % bin_size,nx,ny,nz)
           call memory_alloca(memor_dom,'LIST_ELEMENTS','elebin',element_bin(ielem) % list_elements,nx,ny,nz)
           do kk = 1,nz
              do jj = 1,ny
                 do ii = 1,nx
                    nullify(element_bin(ielem) % list_elements(ii,jj,kk) % l)
                 end do
              end do
           end do
        end do
        !
        ! Initialization
        !
        imin(1:3) = 1_ip
        imax(1:3) = 1_ip
        !
        ! Count number of elements per bin
        !
        do ielem = 1,nelem

           comin(1:3) =  huge(1.0_rp)
           comax(1:3) = -huge(1.0_rp)
           !
           ! Allocate bin
           !
           melel = pelel_2(ielem+1)-1
           ielel = pelel_2(ielem)-1
           kelem = 0

           allocate( xboun(2,3,melel-ielel+1) )

           do while( ielel <= melel )
              !
              ! XBOUN = Element bounding box
              ! COMIN = Neighboring minimum coordinates
              ! COMAX = Neighboring maximum coordinates
              !
              !             1/DELTA(1)
              ! 1/DELTA(2) <--------->
              !       /\   +-+-+-+-+-+ COMAX(3)
              !       ||   | | | | | |
              !       ||   +-+-+-+-+-+
              !       ||   | | | | | |
              !       ||   +-+-+-+-+-+
              !       ||   | | | | | |
              !       \/   +-+-+-+-+-+
              !         COMIN(3)
              !
              if( ielel == pelel_2(ielem)-1 ) then
                 jelem = ielem
              else
                 jelem = lelel_2(ielel)
              end if
              ielel              =  ielel + 1
              kelem              =  kelem + 1
              xboun(1,1:3,kelem) =  huge(1.0_rp)
              xboun(2,1:3,kelem) = -huge(1.0_rp)
              do inode = 1,lnnod(jelem)
                 ipoin = lnods(inode,jelem)
                 xboun(1,1:ndime,kelem) = min( coord(1:ndime,ipoin) , xboun(1,1:ndime,kelem) )
                 xboun(2,1:ndime,kelem) = max( coord(1:ndime,ipoin) , xboun(2,1:ndime,kelem) )
              end do
              delta(1:ndime)         = xboun(2,1:ndime,kelem) - xboun(1,1:ndime,kelem)
              xboun(1,1:ndime,kelem) = xboun(1,1:ndime,kelem) - 0.01_rp*delta(1:ndime) ! Rest a 1% tolerance
              xboun(2,1:ndime,kelem) = xboun(2,1:ndime,kelem) + 0.01_rp*delta(1:ndime) ! Add a 1% tolerance
              comin(1:ndime)         = min( comin(1:ndime) , xboun(1,1:ndime,kelem) )
              comax(1:ndime)         = max( comax(1:ndime) , xboun(2,1:ndime,kelem) )
           end do
           !
           ! Bin min and max coordinates
           !
           element_bin(ielem) % comin = comin
           element_bin(ielem) % comax = comax
           !
           ! Compute number of elements per bin
           !
           delta(1:ndime) = comax(1:ndime)-comin(1:ndime)
           delta(1:ndime) = 1.0_rp / delta(1:ndime)
           dn(1:ndime)    = real(element_bin_boxes(1:ndime),rp) * delta(1:ndime)
           melel          = pelel_2(ielem+1)-1
           ielel          = pelel_2(ielem)-1
           kelem          = 0
           do while( ielel <= melel )
              if( ielel == pelel_2(ielem)-1 ) then
                 jelem = ielem
              else
                 jelem = lelel_2(ielel)
              end if
              ielel         = ielel + 1
              kelem         = kelem + 1
              imin(1:ndime) = int( ( xboun(1,1:ndime,kelem) - comin(1:ndime) ) * dn(1:ndime) , ip )  + 1
              imax(1:ndime) = int( ( xboun(2,1:ndime,kelem) - comin(1:ndime) ) * dn(1:ndime) , ip )  + 1
              imin(1:ndime) = max(imin(1:ndime),1_ip)
              imax(1:ndime) = min(imax(1:ndime),element_bin_boxes(1:ndime))
              do kk = imin(3),imax(3)
                 do jj = imin(2),imax(2)
                    do ii = imin(1),imax(1)
                       element_bin(ielem) % bin_size(ii,jj,kk) = element_bin(ielem) % bin_size(ii,jj,kk) + 1
                    end do
                 end do
              end do
           end do
           !
           ! Allocate list of elements
           !
           do kk = 1,element_bin_boxes(3)
              do jj = 1,element_bin_boxes(2)
                 do ii = 1,element_bin_boxes(1)
                    nullify(element_bin(ielem) % list_elements(ii,jj,kk) % l)
                    nn = element_bin(ielem) % bin_size(ii,jj,kk)
                    call memory_alloca(memor_dom,'ELEMENT_BIN','elebin',element_bin(ielem) % list_elements(ii,jj,kk) % l,nn)
                    element_bin(ielem) % bin_size(ii,jj,kk) = 0
                 end do
              end do
           end do
           !
           ! Fill in bin
           !
           melel = pelel_2(ielem+1)-1
           ielel = pelel_2(ielem)-1
           kelem = 0
           nn    = 0
           do while( ielel <= melel )
              if( ielel == pelel_2(ielem)-1 ) then
                 jelem = ielem
              else
                 jelem = lelel_2(ielel)
              end if
              ielel         = ielel + 1
              kelem         = kelem + 1
              imin(1:ndime) = int( ( xboun(1,1:ndime,kelem) - comin(1:ndime) ) * dn(1:ndime) , ip )  + 1
              imax(1:ndime) = int( ( xboun(2,1:ndime,kelem) - comin(1:ndime) ) * dn(1:ndime) , ip )  + 1
              imin(1:ndime) = max(imin(1:ndime),1_ip)
              imax(1:ndime) = min(imax(1:ndime),element_bin_boxes(1:ndime))
              do kk = imin(3),imax(3)
                 do jj = imin(2),imax(2)
                    do ii = imin(1),imax(1)
                       nn = nn + 1
                       element_bin(ielem) % bin_size(ii,jj,kk) = element_bin(ielem) % bin_size(ii,jj,kk) + 1
                       nn = element_bin(ielem) % bin_size(ii,jj,kk)
                       element_bin(ielem) % list_elements(ii,jj,kk) % l(nn) = jelem
                    end do
                 end do
              end do

           end do

           deallocate( xboun )

!!$     ! POSTPROCESS >
!!$     if( leinv_loc(ielem) == 1 ) then
!!$        print*,'a=',comin(1:2),comax(1:2)
!!$        print*,'b=',leinv_loc(ielem)
!!$        write(100,*) 'MESH ELSEST_BIN dimension 2 Elemtype Quadrilateral Nnode 4'
!!$        write(100,*) 'coordinates'
!!$        nn = 0
!!$        do kk = 1,element_bin_boxes(3)
!!$           do jj = 1,element_bin_boxes(2)
!!$              do ii = 1,element_bin_boxes(1)
!!$                 nn = nn + 1 ; write(100,*) nn , comin(1)+real(ii-1,rp)/dn(1) , comin(2)+real(jj-1,rp)/dn(2)
!!$                 nn = nn + 1 ; write(100,*) nn , comin(1)+real(ii,rp)  /dn(1) , comin(2)+real(jj-1,rp)/dn(2)
!!$                 nn = nn + 1 ; write(100,*) nn , comin(1)+real(ii,rp)  /dn(1) , comin(2)+real(jj,rp)  /dn(2)
!!$                 nn = nn + 1 ; write(100,*) nn , comin(1)+real(ii-1,rp)/dn(1) , comin(2)+real(jj,rp)  /dn(2)
!!$              end do
!!$           end do
!!$        end do
!!$        call memgen(1_ip,npoin,0_ip)
!!$        kelem = melel-ielel+1
!!$        allocate( lnods_tmp(mnode,kelem) )
!!$        melel = pelel_2(ielem+1)-1
!!$        ielel = pelel_2(ielem)-1
!!$        kelem = 0
!!$        do while( ielel <= melel )
!!$           if( ielel == pelel_2(ielem)-1 ) then
!!$              jelem = ielem
!!$           else
!!$              jelem = lelel_2(ielel)
!!$           end if
!!$           ielel = ielel + 1
!!$           do inode = 1,lnnod(jelem)
!!$              ipoin = lnods(inode,jelem)
!!$              if( gisca(ipoin) == 0 ) then
!!$                 nn = nn + 1
!!$                 gisca(ipoin) = nn
!!$              end if
!!$           end do
!!$        end do
!!$        do ipoin = 1,npoin
!!$           if( gisca(ipoin) /= 0 ) then
!!$              nn = gisca(ipoin)
!!$              write(100,*) nn,coord(1:ndime,ipoin)
!!$           end if
!!$        end do
!!$        write(100,*) 'end coordinates'
!!$        write(100,*) 'elements'
!!$        kelem = 0
!!$        nn    = 0
!!$        do kk = 1,element_bin_boxes(3)
!!$           do jj = 1,element_bin_boxes(2)
!!$              do ii = 1,element_bin_boxes(1)
!!$                 kelem = kelem + 1
!!$                 write(100,*) kelem,nn+1,nn+2,nn+3,nn+4
!!$                 nn = nn + 4
!!$              end do
!!$           end do
!!$        end do
!!$        write(100,*) 'end elements'
!!$        write(100,*) 'MESH ELSEST_BIN dimension 2 Elemtype Triangle Nnode 3'
!!$        write(100,*) 'elements'
!!$        do kk = 1,element_bin_boxes(3)
!!$           do jj = 1,element_bin_boxes(2)
!!$              do ii = 1,element_bin_boxes(1)
!!$                 do ll = 1,element_bin(ielem) % bin_size(ii,jj,kk)
!!$                    jelem = element_bin(ielem) % list_elements(ii,jj,kk) % l(ll)
!!$                    nn = maths_mapping_3d_to_1d(element_bin_boxes(1),element_bin_boxes(2),element_bin_boxes(3),ii,jj,kk)
!!$                    write(100,'(5(1x,i5))') kelem,gisca(lnods(1:3,jelem)),nn
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$        write(100,*) 'end elements'
!!$        call memgen(3_ip,npoin,0_ip)
!!$        deallocate(lnods_tmp)
!!$     end if
!!$     ! < POSTPROCESS

        end do

     end if

  end if

end subroutine elebin
