!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_temporal_locality.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute temporal locality
!> @details Temporal locality for assembling an element to
!>          a global matrix
!> @} 
!------------------------------------------------------------------------
subroutine par_temporal_locality()

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER
!  use def_master, only : kfl_paral
  use def_domain, only : lnods,lnnod
  use def_domain, only : r_dom,c_dom
  use def_domain, only : nzdom
  use def_domain, only : lezdo
  use def_kermod, only : kfl_element_to_csr
  use mod_parall, only : num_subd_par
  use mod_parall, only : num_pack_par
  use mod_parall, only : list_elements_par

  implicit none

  integer(ip)            :: ipoin,izdom,jpoin,ielem,jnode
  integer(ip)            :: ipass,inode,jcolu,kelem,itime
  integer(ip)            :: isubd,iwhat,ipack
  integer(ip), parameter :: histo_size=2000
  integer(ip)            :: histo(histo_size),ihisto
  real(rp)               :: xmax_nzdom
  real(rp)               :: xmin_nzdom
  real(rp)               :: xave_nzdom
  real(rp)               :: distance
  integer(ip), pointer   :: last_appear(:)
  integer(ip), pointer   :: num_appear(:)
  real(rp),    pointer   :: xmax(:)
  real(rp),    pointer   :: xmin(:)
  real(rp),    pointer   :: xave(:)

  return

  if( INOTMASTER ) then
     !
     ! WHat to do
     ! IWHAT = 1: compute average distance
     ! IWHAT = 2: Compute all statistics and histogram
     !
     iwhat = 1   
     !
     ! Allocate memory
     !
     nullify(last_appear)
     nullify(num_appear)
     nullify(xmax)
     nullify(xmin)
     nullify(xave)

     allocate(last_appear(nzdom))
     allocate(num_appear(nzdom))
     allocate(xave(nzdom))

     if( iwhat == 2 ) then
        allocate(xmax(nzdom))
        allocate(xmin(nzdom))       
        xmin_nzdom    =  huge(1.0_rp)
        xmax_nzdom    = -huge(1.0_rp)
        xmin(1:nzdom) =  huge(1.0_rp)
        xmax(1:nzdom) = -huge(1.0_rp)
        histo         =  0
     end if

     xave(1:nzdom) =  0.0_rp

     do ipass = 1,iwhat
        !
        ! IPASS = 1: Compute distance
        ! IPASS = 2: Fill in histogram
        !
        last_appear(1:nzdom) = 0
        num_appear(1:nzdom)  = 0
        itime                = 0
        do isubd = 1,num_subd_par
           do ipack = 1,num_pack_par(isubd)
              do kelem = 1,size(list_elements_par(isubd) % packs(ipack) % l)
                 ielem = list_elements_par(isubd) % packs(ipack) % l(kelem)
                 itime = itime + 1
                
                 do inode = 1,lnnod(ielem)
                    ipoin = lnods(inode,ielem)
                    do jnode = 1,lnnod(ielem) 
                       jpoin = lnods(jnode,ielem)
                       !
                       ! Find coefficient of edge IPOIN-JPOIN
                       !
                       if( kfl_element_to_csr == 1 ) then
                          izdom = lezdo(inode,jnode,ielem)
                          jcolu = c_dom(izdom)
                       else
                          izdom = r_dom(ipoin)
                          jcolu = c_dom(izdom)
                          do while( jcolu /= jpoin .and. izdom < r_dom(ipoin+1)-1 )
                             izdom = izdom + 1
                             jcolu = c_dom(izdom)
                          end do
                       end if
                       !
                       ! Fill in statistics for edge IZDOM
                       !
                       num_appear(izdom) = num_appear(izdom) + 1
                       if( num_appear(izdom) >= 2 ) then
                          distance = real(abs(itime-last_appear(izdom)),rp)

                          if( ipass == 1 ) then 
                             if( iwhat == 1 ) then
                                xave(izdom)   = xave(izdom) + distance
                             else
                                xmin(izdom)   = min(xmin(izdom),distance)
                                xmax(izdom)   = max(xmax(izdom),distance)
                                xave(izdom)   = xave(izdom) + distance
                                xmin_nzdom    = min(xmin_nzdom,distance)
                                xmax_nzdom    = max(xmax_nzdom,distance)
                             end if
                          else if( ipass == 2 ) then
                             ihisto        = int( ( distance - xmin_nzdom ) / ( xmax_nzdom-xmin_nzdom ) * real(histo_size,rp) , ip )
                             ihisto        = max(ihisto,1_ip)
                             ihisto        = min(ihisto,histo_size)
                             histo(ihisto) = histo(ihisto) + 1
                          end if

                       end if
                       last_appear(izdom) = itime

                    end do
                 end do
              end do
           end do
        end do
     end do

     do izdom = 1,nzdom
        if( num_appear(izdom) > 1 ) then
           xave(izdom) = xave(izdom) / real(num_appear(izdom)-1,rp)
        else
           xave(izdom) = 0.0_rp
        end if
     end do
     xave_nzdom = sum(xave(1:nzdom)) / real(nzdom,rp)

     if( iwhat == 2 ) then
        xmin_nzdom = minval(xmin(1:nzdom))
        xmax_nzdom = maxval(xmax(1:nzdom))        
        !write(kfl_paral+3000,*) '# ',xmin_nzdom,xmax_nzdom,nzdom
        !do ihisto = 1,histo_size
        !   write(kfl_paral+3000,*) xmin_nzdom + real(ihisto-1,rp)/real(histo_size-1,rp) * ( xmax_nzdom-xmin_nzdom ),histo(ihisto)
        !end do
        !flush(kfl_paral+3000)
     end if

     deallocate(xmin)
     deallocate(xmax)
     deallocate(xave)
     deallocate(last_appear)
     deallocate(num_appear)

  end if

end subroutine par_temporal_locality
