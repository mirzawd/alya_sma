!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_funbou.f90
!> @author  Mariano Vazquez
!> @date    October, 2017-Fix bug only in itask=2 (last time step)
!> @date    November, 2017-Adds smooth step
!> @todo    <GGU> Revision itask 3
!> @brief   Calculation of transient boundary conditions
!> @details This subroutine computes transient boundary conditions
!>          for boundary elements coming from transient data files
!> @}
!-----------------------------------------------------------------------

subroutine sld_funbou(itask,ifixi,xretu)

  use def_kintyp,      only : ip, rp
  use def_master,      only : cutim
  use def_domain,      only : ndime
  use def_solidz,      only : kfl_funbo_sld, kfl_funno_sld
  use def_solidz,      only : kfl_fixno_sld, kfl_bodyf_sld
  use def_solidz,      only : kfl_funty_sld
  use def_solidz,      only : mtloa_sld, tload_sld

  implicit none

  integer(ip), intent(in)   :: itask        !> Boundary type
  integer(ip), intent(in)   :: ifixi        !> Counter iboun, ipoin
  real(rp),    intent(out)  :: xretu(ndime) !> Boundary value
  integer(ip)               :: idata,ifunc,idime
  real(rp)                  :: t1,t2,p1,p2,b,m,tlast,plast

  if (itask == 1_ip) then
     !
     ! on boundary elements
     !
     if (ifixi < 0) then
        ifunc = - ifixi
     else
        ifunc = kfl_funbo_sld(ifixi)
     end if

     !EVA
     ifunc = abs(ifunc)

     do idata = 1,mtloa_sld(ifunc)-1
        ! Tabular times at t1 and t2
        t1 = tload_sld(ifunc)%a(ndime+1,idata)
        t2 = tload_sld(ifunc)%a(ndime+1,idata+1)
        ! Find time interval for current time
        if( cutim >= t1 ) then
           if( cutim <= t2 ) then
              ! Tabular values p1 and p2
              p1 = tload_sld(ifunc)%a(1,idata  )
              p2 = tload_sld(ifunc)%a(1,idata+1)
              if( kfl_funty_sld(2,ifunc) == 1_ip ) then
                 ! Smooth step
                 m = (cutim - t1)/(t2 - t1)
                 xretu(1) = p1 + (p2 - p1)*(m**3)*(10.0_rp - 15.0_rp*m + 6.0_rp*(m**2))
              else
                 ! Linear step (Ramp)
                 m = (p2-p1)/(t2-t1)
                 xretu(1) = p1 + m*(cutim - t1)
              end if
              exit
           end if
        end if
     end do

     ! Correction (last time step) when critical time step is used
     tlast = tload_sld(ifunc)%a(ndime+1,mtloa_sld(ifunc))
     if( cutim > tlast ) then
        plast = tload_sld(ifunc)%a(1,mtloa_sld(ifunc))
        xretu(1) = plast*cutim/tlast
     end if

  else if (itask == 2_ip) then
     !
     ! On Boundary Nodes
     !
     if (ifixi < 0) then
        ifunc = - ifixi
     else
        ifunc = kfl_funno_sld(ifixi)
     end if

     !EVA
     ifunc = abs(ifunc)

     ! Read tabular data and return value
     do idata= 1,mtloa_sld(ifunc)-1
        ! Tabular times at t1 and t2
        t1 = tload_sld(ifunc)%a(ndime+1,idata)
        t2 = tload_sld(ifunc)%a(ndime+1,idata+1)
        ! Find time interval for current time
        if (cutim >= t1) then
           if (cutim <= t2) then
              do idime=1, ndime
                 ! Fixed nodes
                 if (kfl_fixno_sld(idime,ifixi) == 1_ip) then
                    ! Tabular values p1 and p2
                    p1 = tload_sld(ifunc)%a(idime,idata  )
                    p2 = tload_sld(ifunc)%a(idime,idata+1)
                    if (kfl_funty_sld(2,ifunc) == 1_ip ) then
                       ! Smooth step
                       m = (cutim - t1)/(t2 - t1)
                       xretu(idime) = p1 + (p2 - p1)*(m**3.0_rp)*(10.0_rp - 15.0_rp*m + 6.0_rp*(m**2.0_rp))
                    else
                       ! Linear step (Ramp)
                       m  = (p2-p1)/(t2-t1)
                       xretu(idime) = p1 + m*(cutim - t1)
                    end if
                 end if
              end do
              ! Exit when current time is between the interval or match with discrete times
              exit
           end if
        end if
     end do

     ! Correction (last time step) when critical time step is used
     tlast = tload_sld(ifunc)%a(ndime+1,mtloa_sld(ifunc))
     if (cutim > tlast) then
        do idime=1,ndime
           if (kfl_fixno_sld(idime,ifixi) == 1_ip) then
              plast = tload_sld(ifunc)%a(idime,mtloa_sld(ifunc))
              xretu(idime) = plast*cutim/tlast
           end if
        end do
     end if

  else if (itask == 3_ip) then
     !
     ! on every node (body forces)
     !
     ifunc = kfl_bodyf_sld
     if (ifunc == 0) then
        xretu = 0.0_rp
        return
     end if
     do idata= 1,mtloa_sld(ifunc)-1
        if (cutim >= tload_sld(ifunc)%a(ndime+1,idata)) then
           if (cutim < tload_sld(ifunc)%a(ndime+1,idata+1)) then
              t1= tload_sld(ifunc)%a(ndime+1,idata)
              t2= tload_sld(ifunc)%a(ndime+1,idata+1)
              do idime=1,ndime
                 ! if (kfl_fixno_sld(idime,ifixi) == 1) then
                 p1= tload_sld(ifunc)%a(idime,idata  )
                 p2= tload_sld(ifunc)%a(idime,idata+1)
                 m=(p2-p1)/(t2-t1)
                 b=p1-m*t1
                 xretu(idime)=m*cutim+b                              ! return value
                 ! end if
              end do
              exit
           end if
        end if
     end do
  end if

end subroutine sld_funbou
