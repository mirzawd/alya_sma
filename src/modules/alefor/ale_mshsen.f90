!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_mshsen.f90
!> @date    16/11/1966
!> @brief   Updates.
!> @details Updates.
!> @} 
!-----------------------------------------------------------------------
subroutine ale_mshsen
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  
  implicit none
  integer(ip)             :: ipoin,idime,itotn
  
  !
  ! For the adjoint case
  !
  if( INOTMASTER ) then

!       do ipoin=1,npoin
!         if (lninv_loc(ipoin) == 6480) print*, "unkno6480", unkno( (ipoin-1)*ndime+1 : (ipoin-1)*ndime+ndime )
!         if (lninv_loc(ipoin) == 6500) print*, "unkno6500", unkno( (ipoin-1)*ndime+1 : (ipoin-1)*ndime+ndime )
!         if (lninv_loc(ipoin) == 7300) print*, "unkno7300", unkno( (ipoin-1)*ndime+1 : (ipoin-1)*ndime+ndime )
!       enddo
        !-------------------------------------------------------------------
        ! 
        ! Update displacement, mesh velocity and new mesh coordinate
        ! FMALE: 1. does not update coordinate
        !        2. Invert mesh velocity
        !
        !-------------------------------------------------------------------

        do ipoin = 1,npoin
           itotn = (ipoin-1_ip) * ndime ! Before in deform_deform

           do idime = 1,ndime

              if( kfl_fixno_ale(idime,ipoin) == -1 .or. kfl_fixno_ale(idime,ipoin) == 3 ) then

                 !
                 ! FMALE type
                 !
                 ! ---------------------------------------------------------------------------
                 ! Change to use the indexes of coord_ale and dispm in the same way as in the 
                 ! other modules
                 !
                 itotn                    = itotn + 1_ip
                 dispm(idime,ipoin,1)     = unkno(itotn)
                 dispm(idime,ipoin,3)     = dispm(idime,ipoin,1)
!                  coord_ale(idime,ipoin,1) = coord_ale(idime,ipoin,1) + dispm(idime,ipoin,1)
!                  velom(idime,ipoin)       =-dtinv * ( coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,3) )
!                  dispm(idime,ipoin,1)     =  coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,2)
!                  velom(idime,ipoin)       = -dtinv * ( coord_ale(idime,ipoin,1) - coord(idime,ipoin) )
!                  bvess_ale(idime,ipoin)   =  0.0_rp                 
              else
                 !
                 ! Normal type
                 !
                 itotn                    = itotn + 1_ip
                 dispm(idime,ipoin,1)     = unkno(itotn)
                 dispm(idime,ipoin,3)     = dispm(idime,ipoin,1)
!                  coord_ale(idime,ipoin,1) = coord_ale(idime,ipoin,1) + dispm(idime,ipoin,1)
!                  velom(idime,ipoin)       = dtinv * ( coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,3) )
!                  coord(idime,ipoin)       = coord_ale(idime,ipoin,1)
!                  dispm(idime,ipoin,1)   =  coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,2)
!                  velom(idime,ipoin)     =  dtinv * ( coord_ale(idime,ipoin,1) - coord(idime,ipoin) )
              end if

           end do

        end do
     end if

end subroutine ale_mshsen

