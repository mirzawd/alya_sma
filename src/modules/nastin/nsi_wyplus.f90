!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_wyplus()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_wyplus
  ! NAME 
  !    nsi_wyplus
  ! DESCRIPTION
  !    Compute y+ on the wall
  ! USES
  ! USED BY
  !    nsi_
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_gradie
  use mod_memory,                        only : memory_alloca
  use mod_memory,                        only : memory_deallo
  use mod_memory,                        only : memory_size
  use mod_frivel,                        only : frivel
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE

  implicit none
  integer(ip) :: ibopo,idime,ipoin,izdom,jpoin,iwall
  real(rp)    :: rho,mu,nu,ustar,u,y
  real(rp)    :: n1,n2,n3,s11,s12,s22,s13,s23,s33,stnn
  real(rp)    :: sn1,sn2,sn3,stx,sty,stz
  real(rp)    :: delta_aux,dista

  if( delta_nsi == 0.0_rp .and. kfl_delta /= 1 ) then

     !-------------------------------------------------------------------
     !
     ! Up to the wall
     !
     !-------------------------------------------------------------------
     !
     ! Velocity gradients
     ! 
     call memory_alloca(mem_modul(1:2,modul),'GRADV_NSI','nsi_memall',gradv_nsi,ntens,npoin)
     call graten(veloc,gradv_nsi)
     !
     ! y of first node projected on wall
     ! y+ = yU*/nu = y/nu sqrt(Tau_w/rho) = y/sqrt(nu) sqrt(Tau_w') where Tau_w'
     ! is the tangential traction divided by viscosity
     !     
     iwall = 0
     if( associated(walld) ) then
        if( memory_size(walld) > 0 ) then
           iwall = 1
        end if
     end if

     do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if( ibopo > 0 .and. kfl_fixno_nsi(1,ipoin) == 1 ) then
           !
           ! Nearest node JPOIN off the wall at distance y
           !
           y = 0.0_rp
           if( iwall == 1 ) then
              !
              ! Use WALLD
              !
              izdom = r_dom(ipoin)-1
              do while( izdom < r_dom(ipoin+1)-1 ) 
                 izdom = izdom + 1 
                 jpoin = c_dom(izdom)
                 if( walld(jpoin) /= 0.0_rp ) then
                    y     = walld(jpoin)
                    izdom = r_dom(ipoin+1)
                 end if
              end do
           else if( iwall == 0 ) then
              !
              ! Compute shortest distance to interior node
              !
              izdom = r_dom(ipoin)-1
              y     = 1.0e9_rp
              do while( izdom < r_dom(ipoin+1)-1 ) 
                 izdom = izdom + 1 
                 jpoin = c_dom(izdom)
                 if( lpoty(jpoin) == 0 ) then
                    dista = 0.0_rp
                    do idime = 1,ndime
                       dista = dista + (coord(idime,ipoin)-coord(idime,jpoin))**2
                    end do
                    if( dista < y ) y = dista
                 end if
              end do
              y = sqrt(y)
              if( y >= 1.0e9_rp ) y = 0.0_rp
           end if
           !
           ! Tangential traction
           !
           n1    = exnor(1,1,ibopo)
           n2    = exnor(2,1,ibopo)
           s11   = gradv_nsi(1,ipoin)
           s22   = gradv_nsi(2,ipoin)
           s12   = gradv_nsi(3,ipoin)
           sn1   = s11 * n1 + s12 * n2
           sn2   = s12 * n1 + s22 * n2
           rho   = prope_nsi(1,ipoin)
           mu    = prope_nsi(2,ipoin)
           nu    = mu / rho
           if( ndime == 3 ) then
              n3   = exnor(3,1,ibopo) 
              s33  = gradv_nsi(4,ipoin)
              s13  = gradv_nsi(5,ipoin)
              s23  = gradv_nsi(6,ipoin)
              sn1  = sn1 + s13 * n3
              sn2  = sn2 + s23 * n3
              sn3  = s13 * n1 + s23 * n2 + s33 * n3
              stnn = sn1 * n1 - sn2 * n2 - sn3 * n3
              stx  = sn1 - stnn * n1
              sty  = sn2 - stnn * n2
              stz  = sn3 - stnn * n3
              gesca(ipoin) = mu * sqrt( stx*stx + sty*sty + stz*stz )
           else
              stnn = sn1 * n1 + sn2 * n2
              stx  = sn1 - stnn * n1
              sty  = sn2 - stnn * n2
              gesca(ipoin) = mu * sqrt( stx*stx + sty*sty )
           end if
           !
           ! y+ = y/sqrt(nu) sqrt(Tau_w')
           !
           ustar = sqrt( gesca(ipoin) / rho )
           gesca(ipoin) = y * ustar / nu
        else
           gesca(ipoin) = 0.0_rp
        end if
     end do
     !
     ! Take min value 
     !
     if( iwall == 0 .and. ISLAVE ) then
        do ipoin = npoi1+1,npoi2-1
           gesca(ipoin) = 1.0e9_rp
        end do
        do ipoin = npoi3+1,npoin
           gesca(ipoin) = 1.0e9_rp
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MIN')
     end if

     call memory_deallo(mem_modul(1:2,modul),'GRADV_NSI','nsi_wyplus',gradv_nsi)

  else

     !-------------------------------------------------------------------
     !
     ! Wall law
     !
     !-------------------------------------------------------------------
 
     do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if( ibopo /= 0 ) then
           if( kfl_delta==1) then
              delta_aux = ywalp(ibopo)    ! variable wall distance
           else
              delta_aux = delta_nsi       ! fixed wall distance
           end if
           rho   = prope_nsi(1,ipoin)
           mu    = prope_nsi(2,ipoin)
           nu    = mu / rho
           u     = 0.0_rp
           do idime = 1,ndime
              u = u + veloc(idime,ipoin,1) * veloc(idime,ipoin,1)
           end do
           u = sqrt(u)
           if( kfl_rough > 0 ) rough_dom = rough(ipoin) 
           call frivel(kfl_ustar,delta_aux,rough_dom,u,nu,ustar)
           gesca(ipoin) = delta_aux*ustar/nu
        else
           gesca(ipoin) = 0.0_rp
        end if
     end do
  end if

end subroutine nsi_wyplus

