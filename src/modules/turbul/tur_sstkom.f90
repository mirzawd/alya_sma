!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_sstkom(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,eledd,&
     eltur,gpcar,gppro,gpwal,gppr2,gprea,gpdif,gprhs,&
     gpgrd,gpgrv,hleng,fddes,gddes,sstf1,sasso)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_sstkom
  ! NAME
  !   tur_sstkom
  ! DESCRIPTION
  !    Compute coefficient of the equation of Mentter SST model - SST-2003
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
  use def_turbul, only       :  nturb_tur,iunkn_tur,kfl_clipp_tur,&
       &                        param_tur,kfl_ddesm_tur, kfl_sasim_tur, kfl_rotat_tur
  implicit none
  integer(ip), intent(in)    :: ndime,pnode
  real(rp),    intent(in)    :: eledd(pnode),gptur(nturb_tur)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(out)   :: fddes, gddes
  real(rp),    intent(out)   :: sstf1, sasso
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gppro
  real(rp),    intent(in)    :: gpwal,gppr2
  real(rp),    intent(out)   :: gpdif
  real(rp),    intent(out)   :: gpgrd(ndime)
  real(rp),    intent(inout) :: gprea
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(out)   :: gprhs
  real(rp),    intent(in)    :: gpgrv(ndime,ndime),hleng(ndime)
  integer(ip)                :: idime,inode,jdime
  real(rp)                   :: gpkin,gpome,F1,CD
  real(rp)                   :: gradk(3),gradw(3),grakw,c1,c2,c3,F11
  real(rp)                   :: sigk1,sigw1,beta1,sigk2,sigw2,beta,gama
  real(rp)                   :: beta2,betas,gama1,gama2,P,arg1,xsmal
  real(rp)                   :: grU_directo, grU_cruzado, Crc, SS, WW, Ri, F4

  gpkin = max(0.0_rp, gptur(1))
  gpome = gptur(2)
  sigk1 = 0.85_rp !+ 0.000001_rp              
  sigw1 = 0.5_rp  !+ 0.000001_rp              
  beta1 = 0.075_rp
  sigk2 = 1.0_rp
  sigw2 = 0.856_rp            
  beta2 = 0.0828_rp
  betas = 0.09_rp
  xsmal = 1.0e-10_rp
  gama1 = 5.0_rp/9.0_rp
  gama2 = 0.44_rp
  Crc = 1.4_rp  ! Magic number for rotating reference frame

  fddes = 0.0_rp
  gddes = 0.0_rp
  sasso = 0.0_rp

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
  if( ( gpome > xsmal ) .or. ( kfl_clipp_tur > 1 ) ) then
     CD = max( 1.0e-10_rp , 2.0_rp*gpden*sigw2*grakw/gpome )
  else
     !     CD = 1.0e-10_rp   ! error
     CD = max( 1.0e-10_rp , 2.0_rp*gpden*sigw2*grakw/xsmal )
  end if
  !
  ! F1 function. k-w model corresponds to F1 = 1
  !
  if( ( gpome > xsmal ) .or. ( kfl_clipp_tur > 1 ) ) then
     c1   = sqrt( gpkin ) / ( betas * gpome * gpwal )
     c2   = 500.0_rp * gpvis / ( gpden * gpwal * gpwal * gpome )
     c3   = 4.0_rp * gpden * sigw2 * gpkin / ( CD * gpwal * gpwal )
     arg1 = min ( max (  c1 , c2 ) , c3 )
  else
     arg1 = 4.0_rp * gpden * sigw2 * gpkin / ( CD * gpwal * gpwal )
  end if
  F1 = tanh( arg1**4.0_rp )
  sstf1 = F1

  F11          = 1.0_rp - F1
  param_tur(1) = F1 * sigk1 + F11 * sigk2 ! sig_k
  param_tur(2) = F1 * sigw1 + F11 * sigw2 ! sig_w
  gama         = F1 * gama1 + F11 * gama2 ! gama
  beta         = F1 * beta1 + F11 * beta2 ! beta

  param_tur(1) = 1.0_rp / param_tur(1)
  param_tur(2) = 1.0_rp / param_tur(2)

  ! Pure k-w model
  !param_tur(1) = 1.0_rp / sigk1
  !param_tur(2) = 1.0_rp / sigw1
  !beta = beta1
  !gama = gama1
  !F1   = 1.0_rp
  !P = gppro

  if ( kfl_ddesm_tur >= 1 .and.  iunkn_tur == 1) then    ! DDES model flag
     call tur_deslsc(ndime,gptur,gpden,gpvis,gpmut,&
                     gpwal,gpgrv,hleng,fddes,gddes, betas)
  else if( kfl_sasim_tur == 1 .and. iunkn_tur == 2) then
      call tur_sasimu(ndime,pnode,gptur,gpden,gradk,gradw,gppr2,sasso)
  end if

  if( iunkn_tur == 1 ) then
     !
     ! K-equation
     !
     P     = min( gppro , 10.0_rp*gpden*betas*gpome*gpkin )
     gprea = gpden * betas * gpome                   ! (b*)*rho*w
     gprhs = P + gpgra                               ! P+G                  
     if ( kfl_ddesm_tur >= 1 ) then                  ! DDES model flag
        gprea = gprea * fddes 
     end if
  else if( iunkn_tur == 2 ) then
     !
     ! w-equation
     !
     P   = min( gppr2 , 10.0_rp*gpden*betas*gpome*gpkin/gpmut )

     if (kfl_rotat_tur == 1) then !Rotating reference frame
        grU_directo=0.0_rp
        grU_cruzado=0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              grU_directo = grU_directo + gpgrv(idime,jdime)*gpgrv(idime,jdime)
              grU_cruzado = grU_cruzado + gpgrv(idime,jdime)*gpgrv(jdime,idime)
           enddo
        enddo
        ! SS = sqrt(2 Sij Sij)    WW = sqrt(2 Wij Wij)
        ! Where Sij = 0.5(grad_j(u_i)+grad_i(u_j)
        !       Wij = 0.5(grad_j(u_i)-grad_i(u_j)
        SS = sqrt( abs(grU_directo + grU_cruzado) ) ! We need to take abs value because small numerical errors can
        WW = sqrt( abs(grU_directo - grU_cruzado) ) ! introduce random NaNs 
        if (SS .gt. 0.0_rp) then
           Ri = (WW/SS) * (WW/SS -1.0_rp)
        else
           Ri = 0.0_rp
        endif
        F4 = (1.0_rp + Crc * Ri )
        if (F4 .gt. 0.0_rp) then
            F4 = 1.0_rp / F4
        else
            F4 = 1.0_rp ! Not really sure what to do with infinite correction...
        endif
        gprea  = F4 * gpden * beta * gpome   ! F4*b*rho*w  see http://turbmodels.larc.nasa.gov/sst.html
     else
        gprea  = gpden * beta * gpome   ! b*rho*w 
     endif
     
     if(kfl_sasim_tur == 1) then
        gprhs  =  gprhs + gama*gpden*P + sasso
     else
        gprhs  =  gprhs + 0.0_rp*gpgra + gama*gpden*P
     end if
     !
     ! Cross-diffusion term - since this term involves dividing by omega the numerical implementation when omega tends to zero should be analysed further
     !
     if( kfl_clipp_tur > 1 ) then     ! if the values of omega have been clipped 'correctly' we can divide by gpomega without problem 
        gprhs  = gprhs + 2.0_rp*F11*gpden*sigw2/gpome*grakw
     else
        gprhs  = gprhs + 2.0_rp*F11*gpden*sigw2*grakw/max(gpome,xsmal)     ! The previous approach Guillaume had chossen was to totaly eliminate this term when gpome<xsmal - Now we prefer to fix it to the 
                                                                           ! value that would correspond gpome=xsmal. In this way we avoid the function to keep growing but we do not introduce a discontinuity.
     end if

  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_sstkom
