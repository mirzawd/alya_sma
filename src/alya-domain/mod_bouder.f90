!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_bouder

  use def_kintyp, only     :  ip,rp
  implicit none
  
  public :: bouder

contains

!-----------------------------------------------------------------------
!****f* Domain/bouder
! NAME
!    bouder
! DESCRIPTION
!    This routine calculates the baloc and eucta.
! USES
!    vecnor
!    vecpro
! USED BY
!    nsm_bouope
!    tem_bouope
! SOURCE
!***
!-----------------------------------------------------------------------

  pure subroutine bouder(nnodb,ndime,ndimb,deriv,bocod,baloc,eucta)

    integer(ip), intent(in)  :: nnodb,ndime,ndimb
    real(rp),    intent(in)  :: deriv(max(1_ip,ndimb),nnodb)
    real(rp),    intent(in)  :: bocod(ndime,nnodb)
    real(rp),    intent(out) :: eucta,baloc(ndime,ndime)
    integer(ip)              :: i,j,k
     
    if( ndime == 1 ) then
       !
       ! 1D
       !
       baloc(1,1) =   1.0_rp
       eucta      =   baloc(1,1)*baloc(1,1)
       eucta      =   sqrt(eucta)
       
    else if( ndime == 2 ) then
       !
       ! 2D
       !
       do i=1,ndime
          do j=1,ndimb
             baloc(i,j)=0.0_rp
             do k=1,nnodb
                baloc(i,j)=baloc(i,j)+bocod(i,k)*deriv(j,k)
             end do
          end do
       end do

       baloc(1,2) =   baloc(2,1)
       baloc(2,2) = - baloc(1,1)
       eucta      =   baloc(1,2) * baloc(1,2) &
            &       + baloc(2,2) * baloc(2,2)
       eucta      =   sqrt(eucta)

    else if( ndime == 3 ) then
       !
       ! 3D
       !
       do i=1,3
          do j=1,2
             baloc(i,j)=0.0_rp
             do k=1,nnodb
                baloc(i,j)=baloc(i,j)+bocod(i,k)*deriv(j,k)
             end do
          end do
       end do
       
       baloc(1,3) =   baloc(2,1) * baloc(3,2) - baloc(3,1) * baloc(2,2)
       baloc(2,3) =   baloc(3,1) * baloc(1,2) - baloc(1,1) * baloc(3,2)
       baloc(3,3) =   baloc(1,1) * baloc(2,2) - baloc(2,1) * baloc(1,2)
       eucta      =   baloc(1,3) * baloc(1,3) &
            &       + baloc(2,3) * baloc(2,3) &
            &       + baloc(3,3) * baloc(3,3)
       eucta      =   sqrt(eucta)

       ! recalculate t1 so that it is orthogonal to t2
       baloc(1,1) =   baloc(2,3) * baloc(3,2) - baloc(3,3) * baloc(2,2)
       baloc(2,1) =   baloc(3,3) * baloc(1,2) - baloc(1,3) * baloc(3,2)
       baloc(3,1) =   baloc(1,3) * baloc(2,2) - baloc(2,3) * baloc(1,2)

    end if

  end subroutine bouder

end module mod_bouder
