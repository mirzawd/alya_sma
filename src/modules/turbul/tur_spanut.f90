!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_spanut(itask,nu,nup,nut)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_spanut
  ! NAME 
  !    tur_spanut
  ! DESCRIPTION
  !    This routine transforms for Spalart-Allmaras model:
  !    ITASK = 1 ... nut to nu'
  !          = 2 ... nu' to nut
  ! USED BY
  !    tur_updbcs.f90
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_turbul, only       :  param_tur,ipara_tur
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(in)    :: nu
  real(rp),    intent(inout) :: nup,nut
  integer(ip)                :: kiter,iprod
  real(rp)                   :: parco,visc3,newn3,fdvfr,devfr,diffd
  real(rp)                   :: fv1,Xto3,X,cv1t3

  cv1t3 = param_tur(4)**3
  iprod = ipara_tur(1)

  select case(itask)

  case(1_ip)

     if( iprod == 2 ) then 

        !----------------------------------------------------------------
        !
        !
        ! High-Re SA model: nu'=nut
        !
        !----------------------------------------------------------------
        
        nup=nuT

     else

        !----------------------------------------------------------------
        !
        ! Iterate to obtain nu' knowing nut
        !
        !----------------------------------------------------------------

        kiter=0                                           ! i=0
        nup=nuT                                           ! nu'^0=nu_T
        parco=1.0_rp                                     
        visc3=nu*nu*nu           
        if(nuT/nu>1.0e-2_rp) then
           !                     
           ! Newton-Raphson
           !            
           do while((parco>=1.0e-6_rp).and.kiter<100)        ! Newt.-Raphs.
              kiter=kiter+1                               ! i=i+1 
              newn3=nup*nup*nup
              fdvfr=newn3*nup-nuT*newn3-visc3*cv1t3*nuT   ! F(nu'^i)
              diffd=4.0_rp*newn3-3.0_rp*nuT*nup*nup       ! DF(nu'^i)
              devfr=-fdvfr/diffd                          ! delta(nu'^0)
              parco=abs(devfr/(nu+nuT))                   ! nu'^i=nu'^i-1
              nup=nup+devfr                               ! +delta(nu')
              if(nup/nu<1.0e-10_rp) kiter=100
           end do
        else
           !
           ! Picard
           !                
           do while((parco>=1.0e-6_rp).and.kiter<100)        ! Picard
              kiter=kiter+1                               ! i=i+1
              newn3=nup*nup*nup
              diffd=nup
              nup=(nuT*(newn3+visc3*cv1t3))**0.25_rp
              parco=abs((diffd-nup)/(nu+nuT))                   
           end do
        end if
     end if

  case(2_ip)

     !----------------------------------------------------------------
     !
     ! Obtain nut knowing nu'
     !
     !----------------------------------------------------------------

     if( iprod == 2 ) then  
        !
        ! High-Re SA model
        !    
        fv1 = 1.0_rp

     else    
        !
        ! Standard SA model
        !
        X    = nup/nu
        Xto3 = X*X*X
        fv1  = Xto3/(Xto3+cv1t3)
     end if
     nuT = nup*fv1

  end select

end subroutine tur_spanut
