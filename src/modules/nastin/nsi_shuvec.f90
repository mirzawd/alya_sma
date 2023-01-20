!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_shuvec(itask,Auu,App,Apu,invAuu,Aup,xx,yy)
  !
  ! Muyltiply the approximate pressure Schur complement
  ! by a vector
  !
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  npoin,ndime,nzdom
  use def_master, only     :  NPOIN_TYPE,INOTMASTER,IMASTER
  use mod_nsi_schur_operations
  implicit none
  integer(ip), intent(in)  :: itask
  real(rp),    intent(in)  :: Auu(*)
  real(rp),    intent(in)  :: App(nzdom)
  real(rp),    intent(in)  :: Apu(ndime,nzdom)
  real(rp),    intent(in)  :: invAuu(*)
  real(rp),    intent(in)  :: Aup(ndime,nzdom)
  real(rp),    intent(in)  :: xx(npoin)
  real(rp),    intent(out) :: yy(npoin)
  integer(ip)              :: idofn
  real(rp),    pointer     :: zz(:),ww(:)
  real(rp)                 :: dummr(2)

  if( 1 == 1 ) then

     if( INOTMASTER ) then
        !
        ! Allocate memory
        !
        allocate( zz(ndime*npoin) ) 
        allocate( ww(ndime*npoin) ) 
        !
        ! y  = ( App - Apu.diag(Auu)^-1.Aup ) xx
        ! 
        call nsi_solini(2_ip)    
        call nsi_aupvec(1_ip,Aup,xx,zz)                               ! zz = Aup xx
        call pararr('SLX',NPOIN_TYPE,ndime*npoin,zz)

        do idofn = 1,npoin*ndime
           zz(idofn) = zz(idofn) * invAuu(idofn)                      ! zz =  ( diag(Auu)^-1.Aup ) xx
        end do

        call nsi_appvec(1_ip,App,xx,yy)                               ! yy = App xx
        call nsi_apuvec(2_ip,Apu,zz,yy)                               ! yy =  ( Apu.diag(Auu)^-1.Aup ) xx
        if( itask == 1 ) call pararr('SLX',NPOIN_TYPE,npoin,yy)
        !
        ! Deallocate
        !
        deallocate( zz ) 
        deallocate( ww ) 

     end if

  else

     if( IMASTER ) then
        allocate( zz(1) , ww(1) )
     end if

     if( INOTMASTER ) then
        !
        ! Allocate memory
        !
        allocate( zz(ndime*npoin) ) 
        allocate( ww(ndime*npoin) ) 
        !
        ! y  = ( App - Apu.diag(Auu)^-1.Aup ) xx
        !
        call nsi_aupvec(1_ip,Aup,xx,ww)                               ! zz = Aup xx

        do idofn = 1,npoin*ndime
           zz(idofn) = 0.0_rp
        end do

     end if

     call nsi_solini(1_ip)    
     call solver(ww,zz,Auu,dummr)                                     ! Solve system
     call nsi_solini(2_ip)    

     if( INOTMASTER ) then
        call nsi_appvec(1_ip,App,xx,yy)                               ! yy = App xx
        call nsi_apuvec(2_ip,Apu,zz,yy)                               ! yy =  ( Apu.diag(Auu)^-1.Aup ) xx
        if( itask == 1 ) call pararr('SLX',NPOIN_TYPE,npoin,yy)
     else
        !
        ! Deallocate
        !
        deallocate( zz ) 
        deallocate( ww ) 
     end if

  end if

end subroutine nsi_shuvec
