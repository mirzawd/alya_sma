!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine setext(itask,exwor,ierro)
  !-----------------------------------------------------------------------
  !****f* Domain/setext
  ! NAME
  !    setext
  ! DESCRIPTION
  !    This routine computes exterior normals 
  ! OUTPUT
  !    EXNOR(NDIME,NDIME,NBOPO) ... Local basis at boundary nodes
  !    LPOTY(NPOIN) ............... defined by LPOTY(IPOIN)
  !                                 = 0     If IPOIN is interior
  !                                 = IBOPO Number of boundary point
  ! USED BY
  !    extnor
  !***
  !-----------------------------------------------------------------------
  use def_domain
  use def_master,   only    : lninv_loc,zeror
  use def_coupli,   only    : mcoup
  use mod_maths,    only    : maths_local_orthonormal_basis
  use mod_strings,  only    : integer_to_string
  use mod_messages, only    : messages_live 
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(in)    :: exwor(ndime,npoin)
  integer(ip), intent(inout) :: ierro
  integer(ip)                :: ipoin,ibopo
!  integer(ip)                :: iboun,inodb,ielty,iblty,nlboe
  integer(ip)                :: idime,jerro
!  integer(ip)                :: pnode,pnodb,ielem
  real(rp)                   :: xnorm

  if( itask == 1 ) then

     !-------------------------------------------------------------------
     !
     ! EXNOR and LPOTY
     !
     !-------------------------------------------------------------------

     if( ndime == 1 ) then
        !
        ! 1D
        !
        do ibopo = 1,nbopo
           exnor(1,1,ibopo) = 0.0_rp
        end do
        ibopo = 0
        do ipoin = 1,npoin
           if( abs(lpoty(ipoin)) == 1 ) then
              xnorm            = exwor(1,ipoin)*exwor(1,ipoin)
              ibopo            = ibopo+1
              lpoty(ipoin)     = ibopo
              exnor(1,1,ibopo) = exwor(1,ipoin)/sqrt(xnorm) ! Normalize normal
           end if
        end do

     else if( ndime == 2 ) then
        !
        ! 2D
        !
        do ibopo=1,nbopo
           exnor(1,1,ibopo) = 0.0_rp
           exnor(1,2,ibopo) = 0.0_rp
           exnor(2,1,ibopo) = 0.0_rp
           exnor(2,2,ibopo) = 0.0_rp
        end do
        ibopo = 0
        do ipoin = 1,npoin
           if( abs(lpoty(ipoin)) == 1 ) then
              xnorm            =       exwor(1,ipoin)*exwor(1,ipoin)
              xnorm            = xnorm+exwor(2,ipoin)*exwor(2,ipoin)
              ibopo            = ibopo+1
              lpoty(ipoin)     = ibopo
              xnorm            = 1.0_rp/sqrt(xnorm)        ! Normalize normal
              exnor(1,1,ibopo) = exwor(1,ipoin)*xnorm
              exnor(2,1,ibopo) = exwor(2,ipoin)*xnorm
           end if
        end do

     else
        !
        ! 3D
        !
        do ibopo = 1,nbopo
           exnor(1,1,ibopo) = 0.0_rp
           exnor(1,2,ibopo) = 0.0_rp
           exnor(1,3,ibopo) = 0.0_rp
           exnor(2,1,ibopo) = 0.0_rp
           exnor(2,2,ibopo) = 0.0_rp
           exnor(2,3,ibopo) = 0.0_rp
           exnor(3,1,ibopo) = 0.0_rp
           exnor(3,2,ibopo) = 0.0_rp
           exnor(3,3,ibopo) = 0.0_rp
        end do
        ibopo = 0
        do ipoin = 1,npoin
           if( abs(lpoty(ipoin)) == 1 ) then
              xnorm            =       exwor(1,ipoin)*exwor(1,ipoin)
              xnorm            = xnorm+exwor(2,ipoin)*exwor(2,ipoin)
              xnorm            = xnorm+exwor(3,ipoin)*exwor(3,ipoin)
              ibopo            = ibopo+1
              lpoty(ipoin)     = ibopo
              xnorm            = 1.0_rp/sqrt(xnorm)        ! Normalize normal
              exnor(1,1,ibopo) = exwor(1,ipoin)*xnorm
              exnor(2,1,ibopo) = exwor(2,ipoin)*xnorm 
              exnor(3,1,ibopo) = exwor(3,ipoin)*xnorm 
           end if
        end do

     end if
     !
     ! Look for additional boundary nodes not detected by exwor and
     ! explicitely declared by lnodb: OJO check
     !      
     if( ibopo /= nbopo ) then 
        do ipoin = 1,npoin
           if( lpoty(ipoin) == -1 ) then
              ibopo        = ibopo + 1
              lpoty(ipoin) = ibopo
              xnorm = 0.0_rp
              do idime = 1,ndime
                 xnorm = xnorm + exwor(idime,ipoin)*exwor(idime,ipoin)
              end do
              if( xnorm > 20.0_rp*zeror ) then
                 xnorm = 1.0_rp / sqrt(xnorm)
                 do idime = 1,ndime
                    exnor(idime,1,ibopo) = exwor(idime,ipoin)*xnorm
                 end do
              else
                 !
                 ! Enable boudnary nodes to have zero normal:
                 ! One example are extension nodes of Dodeme
                 !
                 exnor(    1,1,ibopo) = 1.0_rp
                 exnor(    2,1,ibopo) = 0.0_rp
                 exnor(ndime,1,ibopo) = 0.0_rp
                 if( mcoup == 0 ) &
                      call messages_live('SETEXT: NORMAL COULD NOT BE COMPUTED FOR NODE '//integer_to_string(lninv_loc(ipoin))//'; ARBITRARILY SET TO 1,0,0 ','WARNING')
                 continue
                 !call runend('SETEXT: ZERO NORMAL FOUND')
              end if
           end if
        end do
     end if
     !if( ibopo /= nbopo ) then 
     !   do iboun = 1,nboun
     !      iblty = ltypb(iboun)
     !      nlboe = nnode(iblty)+1
     !      ielty = ltype(lboel(nlboe,iboun))
     !      do inodb=1,nnode(iblty)
     !         ipoin=lnodb(inodb,iboun)
     !         if(lpoty(ipoin)==-1) then
     !            ibopo        = ibopo+1
     !            lpoty(ipoin) = ibopo
     !            pnodb        = nnode(iblty)
     !            ielem        = lelbo(iboun)
     !            pnode        = nnode(ltype(ielem))
     !            call extcog(iboun,ielem,pnodb,pnode,exnor(1,1,ibopo))
     !         end if
     !      end do
     !   end do
     !end if

  else

     !-------------------------------------------------------------------
     !
     ! EXNOR only
     !
     !-------------------------------------------------------------------

     if( ndime == 1 ) then
        !
        ! 1D
        !
        do ibopo = 1,nbopo
           exnor(1,1,ibopo) = 0.0_rp
        end do
        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 ) then
              ibopo            = lpoty(ipoin)
              xnorm            = exwor(1,ipoin)*exwor(1,ipoin)
              exnor(1,1,ibopo) = exwor(1,ipoin)/sqrt(xnorm) ! Normalize normal
           end if
        end do

     else if( ndime == 2 ) then
        !
        ! 2D
        !
        do ibopo = 1,nbopo
           exnor(1,1,ibopo) = 0.0_rp
           exnor(1,2,ibopo) = 0.0_rp
           exnor(2,1,ibopo) = 0.0_rp
           exnor(2,2,ibopo) = 0.0_rp
        end do

        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 ) then
              ibopo            = lpoty(ipoin)
              xnorm            =       exwor(1,ipoin)*exwor(1,ipoin)
              xnorm            = xnorm+exwor(2,ipoin)*exwor(2,ipoin)
              xnorm            = 1.0_rp/sqrt(xnorm)        ! Normalize normal
              exnor(1,1,ibopo) = exwor(1,ipoin)*xnorm
              exnor(2,1,ibopo) = exwor(2,ipoin)*xnorm
           end if
        end do

     else
        !
        ! 3D
        !
        do ibopo = 1,nbopo
           exnor(1,1,ibopo) = 0.0_rp
           exnor(1,2,ibopo) = 0.0_rp
           exnor(1,3,ibopo) = 0.0_rp
           exnor(2,1,ibopo) = 0.0_rp
           exnor(2,2,ibopo) = 0.0_rp
           exnor(2,3,ibopo) = 0.0_rp
           exnor(3,1,ibopo) = 0.0_rp
           exnor(3,2,ibopo) = 0.0_rp
           exnor(3,3,ibopo) = 0.0_rp
        end do

        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 )  then
              ibopo            = lpoty(ipoin)
              xnorm            =       exwor(1,ipoin)*exwor(1,ipoin)
              xnorm            = xnorm+exwor(2,ipoin)*exwor(2,ipoin)
              xnorm            = xnorm+exwor(3,ipoin)*exwor(3,ipoin)
              xnorm            = 1.0_rp/sqrt(xnorm)        ! Normalize normal
              exnor(1,1,ibopo) = exwor(1,ipoin)*xnorm
              exnor(2,1,ibopo) = exwor(2,ipoin)*xnorm 
              exnor(3,1,ibopo) = exwor(3,ipoin)*xnorm 
           end if
        end do

     end if

  end if
  !
  ! Construct the tangent vector
  ! 
  do ipoin = 1,npoin
     ibopo = lpoty(ipoin)
     if( ibopo /= 0 ) then
        jerro = 0
        call maths_local_orthonormal_basis(ndime,exnor(:,:,ibopo),jerro)
        ierro = ierro + jerro
        if( jerro /= 0 ) then
           call messages_live('SETEXT: NORMAL COULD NOT BE COMPUTED FOR NODE '//integer_to_string(lninv_loc(ipoin)),'WARNING')
        end if
     end if
  end do

end subroutine setext
 
