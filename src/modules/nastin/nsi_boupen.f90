!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_boupen(&
     pevat,pnodb,ndofn,lboel,gbsha,bovel,bovfi,gbvis,&
     gbden,hleng,baloc,wmatr)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_boupen
  ! NAME 
  !    nsi_boupen
  ! DESCRIPTION
  !    Penalize LHS by:
  !          +-
  !    1/eps |  (u.n) v.n dS
  !         -+ S
  !
  !    Units of eps:
  !    1/eps * U * L^2 = rho * U^2 / L * L^3
  !    1/eps = rho * U 
  !    tau   = L / ( rho * U ) =>
  !    eps   = tau / L
  ! USES
  !    vecnor
  !    frivel
  ! USED BY
  !    nsi_bouope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
!  use def_master, only     :  dtinv
  implicit none
  integer(ip), intent(in)  :: pnodb,pevat,ndofn
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbvis
  real(rp),    intent(in)  :: gbden
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: bovel(ndime,pnodb)
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(in)  :: bovfi(ndime,pnodb)
  real(rp),    intent(in)  :: baloc(ndime,ndime)
  real(rp),    intent(out) :: wmatr(pevat,pevat)
  integer(ip)              :: idime,inodb,jevab,ievab,jdime,jnodb
  real(rp)                 :: gbvel(3),fact1,fact2,fact3
  real(rp)                 :: oveps,tau,u,h,rho,mu,eps
  !
  ! GBVEL = velocity
  ! 
  do idime = 1,ndime       
     gbvel(idime) = 0.0_rp
  end do
  do inodb = 1,pnodb
     do idime = 1,ndime            
        gbvel(idime) = gbvel(idime) + gbsha(inodb) * bovel(idime,inodb)
     end do
  end do
  u = 0.0_rp
  do idime = 1,ndime
     u = u + gbvel(idime) * gbvel(idime)
  end do
  !
  ! Characteristic values
  !
  u     = sqrt(u)
  h     = hleng(ndime)
  mu    = gbvis
  rho   = gbden  
  !tau   = 1.0_rp / ( rho*dtinv + 2.0_rp*rho*u/h + 4.0_rp*mu/(h*h) )
  tau   = 1.0_rp / ( 2.0_rp*rho*u/h + 4.0_rp*mu/(h*h) )
!  tau   = 1.0_rp / ( 4.0_rp*mu/(h*h) )
  eps   = tau / h
  oveps = 1000.0_rp / eps
  !
  ! WMATR: Compute matrix
  !
  do inodb = 1,pnodb
     ievab = ( lboel(inodb) - 1 ) * ndofn
     fact1 = oveps * gbsha(inodb)
     do idime = 1,ndime
        ievab = ievab + 1
        fact2 = fact1 * baloc(idime,ndime)
        do jnodb = 1,pnodb
           jevab = ( lboel(jnodb) -1 ) * ndofn
           fact3 = fact2 * gbsha(jnodb)
           do jdime = 1,ndime
              jevab = jevab + 1
              wmatr(ievab,jevab) = wmatr(ievab,jevab) + fact3 * baloc(jdime,ndime) 
           end do
        end do
     end do
  end do

end subroutine nsi_boupen
