!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_ibmfor(ielem,pgaus,gpden,gprhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_ibmfor
  ! NAME 
  !    nsi_ibmfor
  ! DESCRIPTION
  !    F = 1/2 rho * A * Ct * u^2
  !    NS = rho * g + rho * a, such that
  !     +-
  !     | rho * a dV = F =>
  !    -+
  !    a = F / ( rho * V )
  !      = 1/2 rho * A * Ct * u^2 / V
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_elmtyp, only     :  ELCUT
  use def_master, only     :  kfl_coupl,ID_NASTIN,ID_IMMBOU,imbou,letib
  use def_domain, only     :  ndime,lelch,cutel
  implicit none
  integer(ip), intent(in)  :: ielem
  integer(ip), intent(in)  :: pgaus
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(out) :: gprhs(ndime,pgaus)
  integer(ip)              :: igaus,idime,iimbo
  real(rp)                 :: Ct,A,u0,F(3),accel(3)
!  real(rp)                 :: xx(3)

  if( kfl_coupl(ID_NASTIN,ID_IMMBOU) /= 0 ) then

     if( lelch(ielem) == ELCUT .or. letib(ielem) /= 0 ) then
        !
        ! IB
        !     
        if( lelch(ielem) == ELCUT ) then
           iimbo = cutel(ielem) % iimbo
        else
           iimbo = letib(ielem)
        end if

        if( imbou(iimbo) % kfl_coupl > 0 ) then
 
           F        =  0.0_rp
           accel    =  0.0_rp
           Ct       =  0.5_rp
           u0       = 15.0_rp
           A        = 30.0_rp
           accel(1) =  0.5_rp * Ct * A * u0**2 / ( imbou(iimbo) % volum )

           if( lelch(ielem) == ELCUT ) then

              !if( ielem == 2706 ) then
              !   do igaus = 1,pgaus
              !      xx = 0.0_rp
              !      do inode = 1,4
              !         do idime = 1,ndime
              !            xx(idime) = xx(idime) + elcod(idime,inode) * gpsha(inode,igaus)
              !         end do
              !      end do
              !      print*,cutel(ielem) % linou(igaus),xx(1:2)
              !   end do
              !   stop
              !end if

              do igaus = 1,pgaus
                 if( cutel(ielem) % linou(igaus) == -1 ) then
                    do idime = 1,ndime
                       gprhs(idime,igaus) = gprhs(idime,igaus) - gpden(igaus) * accel(idime)
                    end do
                 end if
              end do
              
           else if( iimbo /= 0 ) then
              
              do igaus = 1,pgaus
                 do idime = 1,ndime
                    gprhs(idime,igaus) = gprhs(idime,igaus) - gpden(igaus) * accel(idime)
                 end do
              end do
              
           end if

        end if

     end if

  end if

end subroutine nsi_ibmfor
