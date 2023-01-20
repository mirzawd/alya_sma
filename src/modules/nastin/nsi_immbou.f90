!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_immbou()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_immbou
  ! NAME
  !    nsi_immbou
  ! DESCRIPTION
  !    Compute force on particles
  ! USES
  !    bouder
  !    chenor
  ! USED BY
  !    nsi_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_gradie
  use mod_ker_proper
  use mod_messages, only : livinf
  use mod_bouder
  implicit none

  integer(ip)          :: ipoin,idime,iimbo,imeth
  integer(ip)          :: jdime,izdom,jzdom,jpoin
  real(rp)             :: F1v(3),F1p(3),T1p(3),T1v(3),x(3)  
  real(rp),    pointer :: Fp(:),Fv(:),Tp(:),Tv(:)

  if( nimbo == 0 ) return

  call livinf(59_ip,'COMPUTE FORCES AND TORQUES ON IB',0_ip)
  !
  ! Initialization
  !
  imeth = 1
  do iimbo = 1,nimbo
     Fv => imbou(iimbo) % vforce
     Fp => imbou(iimbo) % pforce
     Tv => imbou(iimbo) % vtorqu
     Tp => imbou(iimbo) % ptorqu
     do idime = 1,3
        Fp(idime)  = 0.0_rp
        Fv(idime)  = 0.0_rp
        Tp(idime)  = 0.0_rp
        Tv(idime)  = 0.0_rp

        F1v(idime) = 0.0_rp
        F1p(idime) = 0.0_rp

        T1p(idime) = 0.0_rp
        T1v(idime) = 0.0_rp

        x(idime)   = 0.0_rp
     end do
  end do
  if (ittim <= 1 ) return
  !
  ! Allocate memory: If gradients are smoothed
  !
  if( INOTMASTER ) then
     !
     ! Loop over boundaries
     !
     do iimbo = 1,nimbo

        Fv     => imbou(iimbo) % vforce
        Fp     => imbou(iimbo) % pforce
        Tv     => imbou(iimbo) % vtorqu
        Tp     => imbou(iimbo) % ptorqu

        !-------------------------------------------------------------
        !
        ! Embedded bodies
        !
        !-------------------------------------------------------------

        do ipoin = 1,npoin
           if( lntib(ipoin) == -iimbo) then
              do idime = 1,ndime
                 F1v(idime) = 0.0_rp
                 F1p(idime) = 0.0_rp 
              end do
              do idime = 1,ndime
                 F1v(idime) = F1v(idime) + intfo_nsi(ipoin) % bu(idime)
                 jzdom = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)
                    jzdom = jzdom + 1
                    do jdime = 1,ndime
                       F1v(idime) = F1v(idime) - intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) * veloc(jdime,jpoin,1) 
                    end do
                    F1p(idime) = F1p(idime) - intfo_nsi(ipoin) % Aup(idime,jzdom) * press(jpoin,1)                           
                 end do
              end do

              x(3) = 0.0_rp   ! so that it has the correct value in the 2d case
              do idime = 1,ndime
                 x(idime) = coord(idime,ipoin) - imbou(iimbo) % posil(idime,1) 
              end do
              call vecpro(x,F1p,T1p,3_ip)                 ! T1 = (X-Xg) x F1  (pressure) 
              call vecpro(x,F1v,T1v,3_ip)                 ! T1 = (X-Xg) x F1  (viscous) 
              !
              ! F, T: Actualize force and torque 
              !
              do idime = 1,3
                 Fv(idime) = Fv(idime) + F1v(idime)
                 Tv(idime) = Tv(idime) + T1v(idime)
                 Fp(idime) = Fp(idime) + F1p(idime)
                 Tp(idime) = Tp(idime) + T1p(idime)
              end do
           end if
        end do
     end do
     
  end if
  !
  ! Reduce sum in Parallel
  !
  call ibmdef(6_ip) 

end subroutine nsi_immbou

