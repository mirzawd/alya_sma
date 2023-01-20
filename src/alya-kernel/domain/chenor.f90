!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    chenor.f90
!> @author  guillaume
!> @date    2021-02-14
!> @brief   Normalize normal
!> @details Normalize exterior normal
!> @} 
!-----------------------------------------------------------------------

subroutine chenor(pnode,baloc,bocod,elcod)
  
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(in)    :: bocod(ndime,*),elcod(ndime,pnode)
  real(rp),    intent(inout) :: baloc(ndime,ndime)
  integer(ip)                :: inode,idime
  real(rp)                   :: produ,cocog(3),dummr
!  integer(ip), intent(in)    :: pnode
!  real(rp),    intent(in)    :: bocod(ndime,*)
!  real(rp),    intent(in)    :: elcod(ndime,pnode)
!  real(rp),    intent(inout) :: baloc(ndime,ndime)
!  integer(ip)                :: idime
!  real(rp)                   :: produ

  do idime = 1,ndime
     dummr                = sqrt(dot_product(baloc(1:ndime,idime),baloc(1:ndime,idime)))
     baloc(1:ndime,idime) = baloc(1:ndime,idime) / dummr
  end do
  return
  
  !
  ! Coordinates center of gravity
  !
  dummr=1.0_rp/real(pnode,rp)

  if(ndime==1) then
     !
     ! 1D
     !
     cocog(1)=0.0_rp
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
     end do
     cocog(1)=cocog(1)*dummr

     produ=baloc(1,1)*baloc(1,1)
     if(produ/=0.0_rp) baloc(1,1)=baloc(1,1)/sqrt(produ)

     produ=(cocog(1)-bocod(1,1))*baloc(1,1)

     if(produ>0.0_rp) then
        baloc(1,1) = -baloc(1,1)         ! n=-n
     end if

  else if(ndime==2) then
     !
     ! 2D
     !
     cocog(1)=0.0_rp
     cocog(2)=0.0_rp
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
        cocog(2)=cocog(2)+elcod(2,inode)
     end do
     cocog(1)=cocog(1)*dummr
     cocog(2)=cocog(2)*dummr

     produ=      baloc(1,1)*baloc(1,1)
     produ=produ+baloc(2,1)*baloc(2,1)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,1)=produ*baloc(1,1)
        baloc(2,1)=produ*baloc(2,1)
     end if

     produ=      baloc(1,2)*baloc(1,2)
     produ=produ+baloc(2,2)*baloc(2,2)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,2)=produ*baloc(1,2)
        baloc(2,2)=produ*baloc(2,2)
     end if

     produ=(cocog(1)-bocod(1,1))*baloc(1,2)&
          +(cocog(2)-bocod(2,1))*baloc(2,2)

     if(produ>0.0_rp) then
        baloc(1,1) = -baloc(1,1)     ! t1=-t1
        baloc(1,2) = -baloc(1,2)     ! n =-n
        baloc(2,1) = -baloc(2,1)     ! t1=-t1
        baloc(2,2) = -baloc(2,2)     ! n =-n
     end if

  else
     !
     ! 3D
     !
     cocog(1)=0.0_rp
     cocog(2)=0.0_rp
     cocog(3)=0.0_rp
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
        cocog(2)=cocog(2)+elcod(2,inode)
        cocog(3)=cocog(3)+elcod(3,inode)
     end do
     cocog(1)=cocog(1)*dummr
     cocog(2)=cocog(2)*dummr
     cocog(3)=cocog(3)*dummr

     produ =       baloc(1,1)*baloc(1,1)
     produ = produ+baloc(2,1)*baloc(2,1)
     produ = produ+baloc(3,1)*baloc(3,1)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,1)=produ*baloc(1,1)
        baloc(2,1)=produ*baloc(2,1)
        baloc(3,1)=produ*baloc(3,1)
     end if
     
     produ =       baloc(1,2)*baloc(1,2)
     produ = produ+baloc(2,2)*baloc(2,2)
     produ = produ+baloc(3,2)*baloc(3,2)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,2)=produ*baloc(1,2)
        baloc(2,2)=produ*baloc(2,2)
        baloc(3,2)=produ*baloc(3,2)
     end if
     produ =       baloc(1,3)*baloc(1,3)
     produ = produ+baloc(2,3)*baloc(2,3)
     produ = produ+baloc(3,3)*baloc(3,3)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,3)=produ*baloc(1,3)
        baloc(2,3)=produ*baloc(2,3)
        baloc(3,3)=produ*baloc(3,3)
     end if

     produ=(cocog(1)-bocod(1,1))*baloc(1,3)&
          +(cocog(2)-bocod(2,1))*baloc(2,3)&
          +(cocog(3)-bocod(3,1))*baloc(3,3)

     if(produ>0.0_rp) then
        baloc(1,1) = -baloc(1,1)     ! t1=-t1
        baloc(1,3) = -baloc(1,3)     ! n =-n
        baloc(2,1) = -baloc(2,1)     ! t1=-t1
        baloc(2,3) = -baloc(2,3)     ! n =-n
        baloc(3,1) = -baloc(3,1)     ! t1=-t1
        baloc(3,3) = -baloc(3,3)     ! n =-n
     end if

  end if

  !! DMM if( produ > 0.0_rp ) print*,'boundary out'
  return
  !do idime = 1,ndime
  !   baloc(:,idime) = baloc(:,idime) / sqrt(dot_product(baloc(:,idime),baloc(:,idime)))
  !end do

end subroutine chenor
