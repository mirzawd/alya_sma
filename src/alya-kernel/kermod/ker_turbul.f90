!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_turbul(itask,dwall,velmo,const,xvisc,gvelo,hleng,xnutu,gpden,tanmo)
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_turbul
  ! NAME 
  !    ker_turbul
  ! DESCRIPTION
  !    This routine computes the turbulent viscosity from the gradients of the velocity field 
  !    at gauss points
  !    Turbulence models available:
  !       --> SMAGORINSKY                              (itask == 1)
  !       --> SMAGORINSKY with Van Driest damping      (itask == 2)
  !       --> WALE                                     (itask == 3)
  !       --> SIGMA                                    (itask == 4)
  !
  ! USED BY
  !    mod_ker_proper
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame

  implicit none
  integer(ip),   intent(in) :: itask
  real(rp),      intent(in) :: gvelo(ndime,ndime)
  real(rp),      intent(in) :: hleng(ndime)
  real(rp),      intent(in) :: gpden
  real(rp),      intent(in) :: dwall
  real(rp),      intent(in) :: velmo,tanmo
  real(rp),      intent(in) :: const
  real(rp),      intent(in) :: xvisc

  real(rp),      intent(out):: xnutu
  integer(ip)               :: idime,jdime,kdime,idumi
  real(rp)                  :: xmile,seci4,g2_ij(ndime,ndime), &
       g2_kk,seci5,sd_ij,denom,G__ij(3,3), &
       G_val(3),vdumm(3,3),sigma(3),ynorm,fvaDr,Smagc,coeff
  real(rp)                  :: alpha,Bbeta

  ! Compute turbulent viscosity xnutu   

  xnutu = 0.0_rp

  if (ndime == 2) then
     xmile = sqrt ( hleng(1) * hleng(2) )
  else
     xmile = ( hleng(1) * hleng(2) * hleng(3) )** 0.3333333_rp
  end if


  select case ( itask )
     
  case ( 1_ip )
     !
     ! SMAGORINSKY 
     ! 
     seci4 = 0.0_rp
     do idime = 1,ndime                     ! 2 S_ij : S_ij
        do jdime = 1,ndime         
           seci4 = seci4 + gvelo(idime,jdime) * (gvelo(idime,jdime) + gvelo(jdime,idime))
        end do
     end do
     Smagc = const*xmile*xmile

     xnutu = Smagc*sqrt(seci4)

  case ( 2_ip )
     !
     ! SMAGORINSKY using Van Driest dumping function (Piomelli et al 1987) 
     ! 
     seci4 = 0.0_rp
     do idime = 1,ndime                     ! 2 S_ij : S_ij
        do jdime = 1,ndime         
           seci4 = seci4 + gvelo(idime,jdime) * (gvelo(idime,jdime) + gvelo(jdime,idime))
        end do
     end do
     !
     ! Van Driest dumping function
     !
     ynorm = 0.04_rp * dwall * (velmo / real(ndime,rp) ) * gpden  / xvisc  
     fvaDr = 1.0_rp - exp(- (ynorm * ynorm * ynorm))
     Smagc = const*xmile*xmile*fvaDr

     xnutu = Smagc*sqrt(seci4)


  case ( 3_ip ) 
     !
     ! WALE  - Wall-adapting local eddy-viscosity - Nicoud-Ducros 99
     ! 
     !
     ! Computation of the square of the velocity gradient g2_ij

     g2_ij = 0.0_rp
     g2_kk = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              g2_ij(idime,jdime) = g2_ij(idime,jdime) + gvelo(idime,kdime)*gvelo(kdime,jdime)
           end do
        end do
     end do

     g2_kk = g2_ij(1,1) + g2_ij(2,2) + g2_ij(3,3)

     seci4 = 0.0_rp
     seci5 = 0.0_rp
     sd_ij = 0.0_rp
     denom = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           seci4 = seci4 + 0.5_rp*gvelo(idime,jdime) * &          ! S_ij : S_ij
                (gvelo(idime,jdime) + gvelo(jdime,idime)) 
           sd_ij = 0.5_rp * ( g2_ij(idime,jdime) + g2_ij(jdime,idime) ) 
!!$           if (idime == jdime) then
!!$              sd_ij = sd_ij - 0.3333333_rp*g2_kk    ! -1/3 delta_ij  * gkk**2
!!$           endif
           seci5 = seci5 + sd_ij*sd_ij     ! Sd_ij:Sd_ij
        end do

     end do
     seci5 = seci5 - g2_kk*g2_kk/3.0_rp

     denom = max ( zeror , seci4**2.5_rp + seci5**1.25_rp )

     xnutu = const*xmile*xmile*(seci5**1.5_rp)/denom

  case ( 4_ip )     
     !
     ! Sigma model  - Baya Toda, Nicoud 2010
     ! 
     ! h calculated as vol**(1/ndime)  
     !
     G__ij = 0.0_rp                      ! G = g^T*g = g_ki * g_kj 
     G_val = 0.0_rp                      ! Note gvelo(kdime,idime) is g_ik = du_i/dx_k from  Nicoud

     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              G__ij(idime,jdime) = G__ij(idime,jdime) + gvelo(idime,kdime)*gvelo(jdime,kdime)
           end do
        end do
     end do

     call spcdec(G__ij,G_val,vdumm,idumi,0_ip,'KER_TURBUL')   ! eigenvalues sorted in descending order


     sigma(1) = sqrt(max(G_val(1),0.0_rp))
     sigma(2) = sqrt(max(G_val(2),0.0_rp))
     sigma(3) = sqrt(max(G_val(3),0.0_rp))

     xnutu = const*xmile*xmile

     if ( sigma(1)*sigma(1) > 1.0e-10_rp ) then     ! Avoid divide by zero
        xnutu = xnutu * ( sigma(3) * ( sigma(1) - sigma(2) ) * ( sigma(2) - sigma(3) ) ) / ( sigma(1) * sigma(1) )
     else
        xnutu = 0.0_rp
     end if

  case ( 5_ip )     
     !
     ! Vreman model  - PoF Vol 16 No 10 Oct. 2004
     ! 
     ! h calculated as vol**(1/ndime)  
     !
     if( ndime == 2 ) then
        G__ij = 0.0_rp ! G = g^T*g
        do idime = 1_ip,2_ip
           do jdime = 1_ip,2_ip
              do kdime = 1_ip,2_ip
                 G__ij(idime,jdime) = G__ij(idime,jdime) + (gvelo(idime,kdime) &
                      *gvelo(jdime,kdime)*xmile*xmile)
                 ! esto es simplificacion pq se podria poner el delta ij
                 ! que toque
              end do
           end do
        end do

        alpha = (G__ij(1_ip,1_ip) + G__ij(2_ip,2_ip) + G__ij(3_ip,3_ip))/(xmile*xmile)
        Bbeta =  G__ij(1_ip,1_ip)*G__ij(2_ip,2_ip) + G__ij(2_ip,2_ip)*G__ij(3_ip,3_ip) + G__ij(3_ip,3_ip)*G__ij(1_ip,1_ip) &
             - G__ij(1_ip,2_ip)*G__ij(1_ip,2_ip) - G__ij(2_ip,3_ip)*G__ij(2_ip,3_ip) - G__ij(1_ip,3_ip)*G__ij(1_ip,3_ip)

        xnutu = const 
        if ( alpha > 10.0_rp*zeror )then                         ! Avoid divide by zero
           xnutu = xnutu  * sqrt (max(( Bbeta ) / ( alpha ), zeror))
        else
           xnutu = 0.0_rp
        end if
     else
        G__ij = 0.0_rp ! G = g^T*g
        do idime = 1_ip,3_ip
           do jdime = 1_ip,3_ip
              do kdime = 1_ip,3_ip
                 G__ij(idime,jdime) = G__ij(idime,jdime) + (gvelo(idime,kdime) &
                      *gvelo(jdime,kdime)*xmile*xmile)
                 ! esto es simplificacion pq se podria poner el delta ij
                 ! que toque
              end do
           end do
        end do

        alpha = (G__ij(1_ip,1_ip) + G__ij(2_ip,2_ip) + G__ij(3_ip,3_ip))/(xmile*xmile)
        Bbeta =  G__ij(1_ip,1_ip)*G__ij(2_ip,2_ip) + G__ij(2_ip,2_ip)*G__ij(3_ip,3_ip) + G__ij(3_ip,3_ip)*G__ij(1_ip,1_ip) &
             - G__ij(1_ip,2_ip)*G__ij(1_ip,2_ip) - G__ij(2_ip,3_ip)*G__ij(2_ip,3_ip) - G__ij(1_ip,3_ip)*G__ij(1_ip,3_ip)

        xnutu = const 
        if ( alpha > 10.0_rp*zeror )then                         ! Avoid divide by zero
           xnutu = xnutu  * sqrt (max(( Bbeta ) / ( alpha ), zeror))
        else
           xnutu = 0.0_rp
        end if
     end if

  case ( 6_ip )
     !
     ! Algebraic mixing length model using Van Driest dumping function
     ! \nu_{Tml}=(\kappa y^{+})^{2} |S| [ 1 - exp (-y^{+}/A^{+}) ]^{2}
     ! y^{+} = y * u_{*} / \nu
     ! u_{*} = sqrt ( tanmo / rho )
     ! |S| = sqrt (2 S_ij : S_ij)
     ! 
     seci4 = 0.0_rp
     do idime = 1,ndime                     ! 2 S_ij : S_ij
        do jdime = 1,ndime         
           seci4 = seci4 + gvelo(idime,jdime) * (gvelo(idime,jdime) + gvelo(jdime,idime))
        end do
     end do
     !
     ! Van Driest dumping function
     !
     ynorm = dwall * sqrt ( tanmo / gpden) /  xvisc  ! yplus
     fvaDr =  1.0_rp - exp ( -ynorm / 26.0_rp )
     coeff = ( 0.41_rp * ynorm * fvaDr ) ** 2

     xnutu = coeff * sqrt(seci4)

  end select

end subroutine ker_turbul
