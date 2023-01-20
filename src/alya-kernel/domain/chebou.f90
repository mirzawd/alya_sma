!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    chebou.f90
!> @author  Guillaume Houzeaux
!> @date    18/10/2012
!> @brief   Check if normal is outwards 
!!          Check if normal is outwards
!!
!!          \verbatim
!!
!!          o----------o
!!          |          |    cge = element center of gravity
!!          |          |    cgb = boundary center of gravity
!!          |    cge   |    v   = vector(cgb,cge)
!!          |    /\    |
!!          |  v ||    |
!!          o----cgb---o <- iboun
!!               ||
!!               \/ n
!!
!!          \endverbatim
!!
!!          The procedure is the following:
!!          - Compute element center of gravity COCOG
!!          - Compute boundary center of gravity COCOB
!!          - v = COCOG - COCOB
!!          - Compute n.v
!!          - Invert BALOC if necessary
!!
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine chebou(pnode,pnodb,baloc,bocod,elcod,ninve)
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb
  real(rp),    intent(in)    :: bocod(ndime,pnodb),elcod(ndime,pnode)
  real(rp),    intent(inout) :: baloc(ndime)
  integer(ip), intent(inout) :: ninve
  integer(ip)                :: inode,inodb
  real(rp)                   :: produ,cocog(3),cocob(3),dummr,dummb
  !
  ! Coordinates center of gravity
  !
  dummr = 1.0_rp / real(pnode,rp)
  dummb = 1.0_rp / real(pnodb,rp)

  if( ndime == 1 ) then

     !-------------------------------------------------------------------
     !
     ! 1D
     !
     !-------------------------------------------------------------------

     cocog(1)=0.0_rp 
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
     end do
     cocog(1)=cocog(1)*dummr

     produ = baloc(1)*baloc(1)
     if( produ /= 0.0_rp ) baloc(1)=baloc(1)/sqrt(produ)

     produ=(cocog(1)-bocod(1,1))*baloc(1)

     if( produ > 0.0_rp ) then
        baloc(1) = -baloc(1)         ! n=-n                      
     end if

  else if( ndime == 2 ) then

     !-------------------------------------------------------------------
     !
     ! 2D
     !
     !-------------------------------------------------------------------

     cocog(1)=0.0_rp
     cocog(2)=0.0_rp
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
        cocog(2)=cocog(2)+elcod(2,inode)
     end do
     cocog(1)=cocog(1)*dummr
     cocog(2)=cocog(2)*dummr

     cocob(1)=0.0_rp
     cocob(2)=0.0_rp
     do inodb=1,pnodb
        cocob(1)=cocob(1)+bocod(1,inodb)
        cocob(2)=cocob(2)+bocod(2,inodb)
     end do
     cocob(1)=cocob(1)*dummb
     cocob(2)=cocob(2)*dummb

     produ= baloc(1)*baloc(1) &
          + baloc(2)*baloc(2)
     if( produ /= 0.0_rp ) then
        produ    = 1.0_rp/sqrt(produ)
        baloc(1) = produ*baloc(1)
        baloc(2) = produ*baloc(2)
     end if

     produ=(cocog(1)-cocob(1))*baloc(1)&
          +(cocog(2)-cocob(2))*baloc(2)

     if(produ>0.0_rp) then
        baloc(1) = -baloc(1)     ! n =-n                    
        baloc(2) = -baloc(2)     ! n =-n
     end if

  else

     !-------------------------------------------------------------------
     !
     ! 3D
     !
     !-------------------------------------------------------------------

     cocog(1) = 0.0_rp
     cocog(2) = 0.0_rp
     cocog(3) = 0.0_rp
     do inode = 1,pnode
        cocog(1) = cocog(1) + elcod(1,inode)
        cocog(2) = cocog(2) + elcod(2,inode)
        cocog(3) = cocog(3) + elcod(3,inode)
     end do
     cocog(1) = cocog(1) * dummr
     cocog(2) = cocog(2) * dummr
     cocog(3) = cocog(3) * dummr

     cocob(1) = 0.0_rp
     cocob(2) = 0.0_rp
     cocob(3) = 0.0_rp
     do inodb = 1,pnodb
        cocob(1) = cocob(1) + bocod(1,inodb)
        cocob(2) = cocob(2) + bocod(2,inodb)
        cocob(3) = cocob(3) + bocod(3,inodb)
     end do
     cocob(1) = cocob(1) * dummb
     cocob(2) = cocob(2) * dummb
     cocob(3) = cocob(3) * dummb

     produ = baloc(1)*baloc(1) &
          +  baloc(2)*baloc(2) &
          +  baloc(3)*baloc(3)

     if( produ /= 0.0_rp ) then
        produ    = 1.0_rp/sqrt(produ)
        baloc(1) = produ*baloc(1)
        baloc(2) = produ*baloc(2)
        baloc(3) = produ*baloc(3)
     end if

     produ =   ( cocog(1) - cocob(1) ) * baloc(1) &
          &  + ( cocog(2) - cocob(2) ) * baloc(2) &
          &  + ( cocog(3) - cocob(3) ) * baloc(3)

     if( produ > 0.0_rp ) then                   
        ninve    =  ninve + 1
        baloc(1) = -baloc(1)     ! n =-n                  
        baloc(2) = -baloc(2)     ! n =-n                  
        baloc(3) = -baloc(3)     ! n =-n
     end if

  end if

end subroutine chebou
