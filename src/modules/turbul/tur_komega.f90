!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_komega(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,eledd,&
     gpcar,gppro,gppr2,gpgrv,eltur,gprea,gpdif,gprhs,&
     gpgrd)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_komega
  ! NAME
  !   tur_komega
  ! DESCRIPTION
  !    Compute coefficient of the equation of k-w model
  !
  !    Transport equations
  !    -------------------
  !   
  !       dK         dK                        d +-          dK  -+
  !    rho-- + rho*ui--- = P - (b*)*rho*K*w + ---|(mu+mut/sk)---  |
  !       dt         dxi                      dxi+-          dxi -+
  !
  !       dw         dw       w                  d +-          dw  -+
  !    rho-- + rho*ui--- =  a - P - b*rho*w^2 + ---|(mu+mut/sw)---  |
  !       dt         dxi      K                 dxi+-          dxi -+
  !
  !    with  sk=2, sw=2, b=3/40, a=5/9
  !          Rk=6, Rw=27/10, Rb=8, a0*=b/3, a0=1/10
  !
  !               k             dui
  !    and ReT = ----,  P = tij --- 
  !              nu*w           dxj  
  !
  !         a0*+ReT/Rk       5 a0+ReT/Rw                9  5/18+(ReT/Rb)^4
  !    a* = ---------- , a = - ---------(a*)^-1 , b* = --- ---------------
  !           1+ReT/Rk       9  1+ReT/Rw               100    1+(ReT/Rb)^4
  !
  !    Auxiliary relations
  !    -------------------
  !   
  !    L = k^{1/2}/w, eps = b* w k
  !    mu_t=rho*(alpha*)*k/w
  !
  ! OUTPUT 
  !    GPREA .......... r 
  !    GPDIF .......... k 
  !    GPRHS .......... f 
  !    GPGRD(NDIME) ... grad(k) coefficient
  ! USES
  ! USED BY
  !    tur_elmcoe
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_turbul, only       :  nturb_tur,iunkn_tur,kfl_algor_tur,&
       &                        param_tur,ipara_tur
  implicit none
  integer(ip), intent(in)    :: ndime,pnode
  real(rp),    intent(in)    :: eledd(pnode),gptur(nturb_tur)
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gppro,gppr2
  real(rp),    intent(in)    :: gpgrv(ndime,ndime)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(out)   :: gpdif(kfl_algor_tur)
  real(rp),    intent(out)   :: gpgrd(kfl_algor_tur,ndime)
  real(rp),    intent(inout) :: gprea(kfl_algor_tur,kfl_algor_tur)
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(out)   :: gprhs(kfl_algor_tur)
  real(rp)                   :: gpkin,gpome,ReToRb4
  real(rp)                   :: bs,ReT,ReToRk,ReToRw,dummr
  real(rp)                   :: gradk(3),gradw(3),grakw,xk,fbs,fb,b0,bs0,xw
  real(rp)                   :: S(ndime,ndime),W(ndime,ndime)
  integer(ip)                :: idime,jdime,kdime,inode

  gpkin=gptur(1)
  gpome=gptur(2)

  if( ipara_tur(1) == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Modified k-w (Wilcox 1998)
     !
     !-------------------------------------------------------------------

     if(iunkn_tur==1) then

        gradk = 0.0_rp
        gradw = 0.0_rp
        grakw = 0.0_rp
        do inode = 1,pnode
           do idime = 1,ndime
              gradk(idime) = gradk(idime) + eltur(1,inode) * gpcar(idime,inode)
              gradw(idime) = gradw(idime) + eltur(2,inode) * gpcar(idime,inode)
           end do
        end do
        do idime = 1,ndime
           grakw = grakw + gradk(idime) * gradw(idime)
        end do
        xk = 1.0_rp/(gpome*gpome*gpome)*grakw
        if( xk <= 0.0_rp ) then
           fbs = 1.0_rp
        else
           fbs = ( 1.0_rp + 680.0_rp * xk * xk ) / ( 1.0_rp + 400.0_rp * xk * xk )
        end if
        bs0 = 0.09_rp
        bs  = bs0 * fbs
        gprea(1,1) = gpden*bs*gpome                          ! (b*)*rho*w
        gprhs(1)   = gppro + gpgra                           ! P+G

     else

        do idime = 1,ndime
           do jdime = 1,ndime
              S(idime,jdime) = 0.50_rp * ( gpgrv(jdime,idime) + gpgrv(idime,jdime) )
              W(idime,jdime) = 0.50_rp * ( gpgrv(jdime,idime) - gpgrv(idime,jdime) )
           end do
        end do
        xw = 0.0_rp
         do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 xw = xw + W(idime,jdime)*W(jdime,kdime)*S(kdime,idime)
              end do
           end do
        end do
        bs0 = 0.09_rp
        b0  = 9.0_rp/125.0_rp
        xw  = abs(xw)/( (bs0*gpome)**3.0_rp )
        fb  = ( 1.0_rp + 70.0_rp * xw ) / ( 1.0_rp + 80.0_rp * xw )
        gprhs(1)    = 13.0_rp/25.0_rp*gpden*gppr2
        gprea(1,1)  = gpden*b0*fb*gpome                ! b0*fb*w

     end if

  else if (ipara_tur(1)==0) then

     !-------------------------------------------------------------------
     !
     ! Original k-w
     !
     !-------------------------------------------------------------------

     if(kfl_algor_tur==1) then
        !
        ! Decoupled
        !
        if(iunkn_tur==1) then
           !
           ! K-equation
           !
           if(gpome<1.0e-10_rp) then
              bs      = param_tur(6)
           else
              ReT     = gpden*gpkin/(gpvis*gpome)        
              ReToRk  = ReT/param_tur(7)
              ReToRb4 = (ReT/param_tur(9))**4.0_rp
              bs      = param_tur(6)*(5.0_rp/18.0_rp+ReToRb4)&
                   &    /(1.0_rp+ReToRb4)
           end if
!!!bs=0.09_rp ! OJO
           gprea(1,1) = gpden*bs*gpome                          ! (b*)*rho*w
           gprhs(1)   = gppro + gpgra                           ! P+G

        else if(iunkn_tur==2) then
           !
           ! w-equation
           !
           dummr     = gpden*param_tur(3)*gppr2                 ! rho*alpha*P'
           if(gpome>1.0e-10_rp) then
              ReT    = gpden*gpkin/(gpvis*gpome)                ! a*
              ReToRw = ReT/param_tur(8)
              dummr  = dummr*(param_tur(11)+ReToRw)&
                   &         /(1.0_rp      +ReToRw)
           end if
!!!dummr=gpden*5.0_rp/9.0_rp*gppr2 ! OJO
           gprhs(1)    = dummr +  0.0_rp*gpgra                 ! rho*a*P'*(a*) + G
           gprea(1,1)  = gpden*param_tur(4)*gpome              ! b*rho*w

        end if

     else
        !
        ! Coupled: k-equation
        !
        gprhs(1)   = gpgra     
        if(gpome<1.0e-10_rp) then
           bs      = param_tur(6)
        else
           ReT     = gpden*gpkin/(gpvis*gpome)        
           ReToRk  = ReT/param_tur(7)
           ReToRb4 = (ReT/param_tur(9))**4.0_rp
           bs      = param_tur(6)*(5.0_rp/18.0_rp+ReToRb4)/(1.0_rp+ReToRb4)
        end if
        gprea(1,1) = gprea(1,1) + gpden*bs*gpome
        gprea(1,2) = gprea(1,2) + gpden*bs*gpkin
        gprhs(1)   = gprhs(1)   + gppro + gpden*bs*gpome*gpkin
        !
        ! Coupled: w-equation
        !
        gprhs(2) = 0.0_rp*gpgra             ! What should it be?
        dummr    = gpden*param_tur(3)*gppr2 ! alpha*rho*P'
        if(gpome>1.0e-10_rp) then
           ReT    = gpden*gpkin/(gpvis*gpome)
           ReToRw = ReT/param_tur(8)
           dummr  = dummr*(param_tur(11)+ReToRw)&
                &         /(1.0_rp      +ReToRw)
        end if
        gprhs(2)    = gprhs(2)   + dummr
        gprea(2,2)  = gprea(2,2) + gpden*param_tur(4)*gpome

     end if
  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_komega
