!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bousto(&
     pevat,pnodb,wmatr,tract, bvnat, bovel, lnodb, gbsha,&
     baloc,gbden, lboel, ndofn, iboun)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bousto
  ! NAME 
  !    nsi_bousto
  ! DESCRIPTION
  !
  ! Stable outflow condition with resistance
  ! 
  !    Stable outflow condition with resistance (Alfonsito Santiago)
  !    from Bazilevs et al (doi:10.1016/j.cma.2009.04.015)
  !    This routine computes the surface traction for the NS equations at
  !    a given integration point of an open outflow boundary IBOUN received 
  !    by argument. Implements the algorithm of Bazilevs et al
  !                     +-
  ! sig.n = -p0 n - C * | u.n ds  + beta rho*(u.n)_ u 
  !                    -+S
  !           +-
  !           | u.n if u.n < 0 (inflow)
  ! (u.n)_ = -|
  !           | 0   otherwise  (outflow)
  !           +-
  !
  ! p0   = BVNAT(1)
  ! C    = BVNAT(2)
  ! beta = BVNAT(3)
  ! USES
  ! 
  ! USED BY
  !    nsi_bouope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime
!  use def_domain, only      :  kfl_codbo
!  use def_nastin, only      :  kfl_funtb_nsi
  use def_nastin, only      :  kfl_hydro_nsi, bpess_nsi
!  use def_nastin, only      :  outflow_mass
  use mod_windk,  only      :  mod_windk_interact
  implicit none
  integer(ip), intent(in)    :: pevat
  integer(ip), intent(in)    :: pnodb
  real(rp),    intent(inout) :: wmatr(pevat,pevat)
  real(rp),    intent(inout) :: tract(ndime)
  real(rp),    intent(in)    :: bvnat(5)
  real(rp),    intent(in)    :: bovel(ndime,pnodb)
  integer(ip), intent(in)    :: lnodb(pnodb)
  real(rp),    intent(in)    :: gbsha(pnodb)
  real(rp),    intent(in)    :: baloc(ndime,ndime)
  real(rp),    intent(in)    :: gbden
  integer(ip), intent(in)    :: lboel(pnodb)
  integer(ip), intent(in)    :: ndofn
  integer(ip), intent(in)    :: iboun
  
  ! local variables
  integer(ip)              :: iflow, idime,jdime, inodb, ipoin, ievab, jevab, jnodb
  real(rp)                 :: gbvel(ndime), resistive_term(ndime), gbpre, rhoun, udotn, fact1, one, fact2
!  real(rp)                 :: p_in, rdummi
  integer (ip)             :: kfl_linearization ! linearization flag

  !
  ! linearization flag for term : beta* rho*(u.n)_ u, when evaluating u
  !
  kfl_linearization = 2     ! = 0 explicit, =1, implicit fixed point (Picard), =2 full Newton Raphson
  !  
  if (kfl_linearization == 0) then      ! Explicit  beta* rho*(u^{i}.n)_ u^{i}
     one = 1.0_rp 
  else if (kfl_linearization == 1) then ! Picard beta* rho*(u^{i}.n)_ u^{i+1}
     one = 0.0_rp
  else if(kfl_linearization == 2) then  ! Newton Raphson beta* rho*(u^{i}.n)_ u^{i+1} 
     one =-1.0_rp           ! +beta* rho*(u^{i+1}.n)_ u^{i} - beta* rho*(u^{i}.n)_ u^{i}
  end if
  
  iflow = int(bvnat(4),ip)
  gbvel = 0.0_rp
  do inodb = 1,pnodb
     gbvel(1:ndime) = gbvel(1:ndime) + bovel(1:ndime,inodb) * gbsha(inodb)
  end do
  udotn = dot_product(gbvel(1:ndime),baloc(1:ndime,ndime))
  ! - rho*(u.n)_ u term
  rhoun =  -bvnat(3) * gbden * min(udotn,0.0_rp)
 
 
  resistive_term(1:ndime) = bvnat(5) * baloc(1:ndime,ndime)
  
!!$  if( kfl_funtb_nsi(iboun) .eq. 0) then ! Regular Bazilevs condition
!!$    !
!!$    ! This will never work for high values of baloc without
!!$    ! an implicit treatment of the term or, at least, a relaxation
!!$    ! scheme
!!$    !
!!$    resistive_term(1:ndime) = bvnat(2) * outflow_mass(iflow)* baloc(1:ndime,ndime) !(C*\int u.n ds)n
!!$
!!$    if(bvnat(2).gt.1.0_rp)  call runend('NSI_BOUSTO: regular outflow will never work with high (<1.0) values of C. implict treatment required')
!!$
!!$  elseif( kfl_funtb_nsi(iboun) .eq. 3) then ! External windkessel condition
!!$
!!$    p_in=0.0_rp
!!$  
!!$    rdummi=0.0_rp
!!$    call  mod_windk_interact(rdummi, p_in,'NASTIN', kfl_codbo(iboun), OPT_direction='OU', IS_EXPLICIT=NSI_FRACTIONAL_STEP)
!!$    !print*,'a=',p_in,bvnat(5)
!!$    p_in = bvnat(5)
!!$
!!$    resistive_term(1:ndime) =  p_in * baloc(1:ndime,ndime) !(C*\int u.n ds)n
!!$ else
!!$   call runend('nsi_bousto: type of boundary condition not yet programmed')
!!$ endif
 

  
  if(kfl_hydro_nsi/=0)then
     !
     ! Pressure: sig.n = - p n, with p is nodal + bazylev
     !
     gbpre = 0.0_rp
     do inodb = 1,pnodb
        ipoin = lnodb(inodb)
        gbpre = gbpre + bpess_nsi(1,ipoin,1) * gbsha(inodb)
     end do
     tract(1:ndime) = - gbpre * baloc(1:ndime,ndime)  & ! - p n    
          - resistive_term(1:ndime) & !-(C*\int u.n ds)n
          - one*rhoun* gbvel(1:ndime)  ! + rho*(u.n)_ u

  else
     tract(1:ndime) = - bvnat(1)*baloc(1:ndime,ndime) & ! - p0 n
          - resistive_term(1:ndime) & !-(C*\int u.n ds)n
          - one*rhoun* gbvel(1:ndime)  ! + rho*(u.n)_ u
  end if


  !
  ! implicit assembly of rho*(u.n)_ u term
  !
  if (kfl_linearization.ge.1.and.udotn.lt.0.0_rp) then   ! matrix assembly
     do inodb = 1,pnodb
        fact1 = rhoun * gbsha(inodb)
        ievab = (lboel(inodb)-1)*ndofn
        do jnodb = 1,pnodb
           jevab=(lboel(jnodb)-1)*ndofn 
           do idime =1, ndime             
              wmatr(ievab+idime ,jevab+idime) = wmatr(ievab+idime ,jevab +idime) & 
                   +fact1*gbsha(jnodb)
              if (kfl_linearization.eq.2) then
                 do jdime =1, ndime
                    fact2= - bvnat(3)*gbden*gbvel(idime)*baloc(jdime, ndime)*gbsha(inodb)
                    wmatr(ievab+idime ,jevab+jdime) = wmatr(ievab+idime ,jevab +jdime) & 
                         +fact2*gbsha(jnodb)
                 end do
              end if
           end do
        end do
     end do
  end if

end subroutine nsi_bousto
