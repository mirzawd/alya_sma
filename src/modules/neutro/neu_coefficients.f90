!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_coefficients.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Coefficients
!> @details Here we calculate the source (rhs) term (including the scattered radiation from other directions), &
!>          the current velocity (direction vector) and the total cross-section.
!> @} 
!-----------------------------------------------------------------------
subroutine neu_coefficients(pnode,pgaus,pmate,gpsha,gpabs,gpsca,elrad,elcod,elunk,gpadv,gprea,gprhs,gpunk)
  use def_kintyp, only     :  ip,rp 
  use def_parame, only     :  pi,Stefan_Boltzmann,in4pi
  use def_domain, only     :  ndime
  use def_neutro, only     :  ADR_neu
  use def_neutro, only     :  num_energies_neu,num_directions_neu
  use def_neutro, only     :  current_energy_neu
  use def_neutro, only     :  current_direction_neu
  use def_neutro, only     :  direc_neu
  use def_neutro, only     :  weigd_neu
  use def_neutro, only     :  num_legendre_neu,phi_neu,tita_neu,scatt_neu!,scattering_neu,ener_weigd_neu
  use def_neutro, only     :  funsource!,efectivos_neu_in,efectivos_neu_out
  use def_neutro, only     :  efectivos_neu_sparse, n_efectivos_neu_sparse_acum!, n_efectivos_neu_sparse

  implicit none

  integer(ip), intent(in)  :: pnode 
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: pmate
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpabs(pgaus)
  real(rp),    intent(in)  :: gpsca(pgaus)
  real(rp),    intent(in)  :: elrad(num_energies_neu,num_directions_neu,pnode)
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: elunk(pnode,ADR_neu % ntime)
  real(rp),    intent(out) :: gpadv(ndime,pnode)
  real(rp),    intent(out) :: gprea(pnode) 
  real(rp),    intent(out) :: gprhs(pnode)
  integer(ip)              :: inode,itime!,leff
  integer(ip)              :: igaus,iener,idire,jdire,enerG,en,Lcoef,mcoef!,Kcoef
  real(rp)                 :: gprad(num_energies_neu,num_directions_neu) !> This variable will only be used to calculate the 
                                                                         ! radiation scattered from other directions inthe RHS term
  real(rp)                 :: gpunk(pgaus,ADR_neu % ntime),gpsou(pgaus) !> Values for unkown at Gauss points for each time step
  
  real(rp) :: xmed,ymed,armonic!,Termino_resta
  real(rp), external :: Armo_source_enmater_neu, funharmonic

  gpsou = 0.0_rp !> We keep the source in case we wanted it for the future
  idire = current_direction_neu !> We set the index for the direction into the current direction we consider
  gpunk = 0.0_rp 
  gprea = 0.0_rp
  enerG = current_energy_neu
  
  do igaus = 1,pgaus !>Iterating over the Gauss points 
     !
     ! Source term
     !
     if(funsource(pmate)==-1) then
       
       xmed=0.0_rp
       ymed=0.0_rp
       do inode=1,pnode
           xmed=xmed + elcod(1,inode)/real(pnode, rp)
           ymed=ymed + elcod(2,inode)/real(pnode, rp)
       enddo
       do inode=1,pnode
         !  gpsou(igaus) = gpsou(igaus) + Armo_source_enmater_neu(xmed,ymed,direc_neu(1,idire),gpabs(igaus),gpsca(igaus))*gpsha(inode,igaus)
          gpsou(igaus) = gpsou(igaus) + Armo_source_enmater_neu(xmed,direc_neu(1,idire),gpabs(igaus),gpsca(igaus))*gpsha(inode,igaus)
       enddo

     else
       
       do inode=1,pnode
          gpsou(igaus) =   gpsou(igaus) + funsource(pmate)* gpsha(inode,igaus)
       enddo

     endif

     gprad = 0.0_rp 
     do inode = 1,pnode !> Iterating over the number of nodes of the element considered
        !> We add the contribution of each node of the element 
        ! opt 1: (radiation*shape at the Gauss point)
        gprad(1:num_energies_neu,1:num_directions_neu) = gprad(1:num_energies_neu,1:num_directions_neu) &
             + gpsha(inode,igaus) * (elrad(1:num_energies_neu,1:num_directions_neu,inode)) 

      !   ! opt 2: (radiation/(# of element nodes))
      !   gprad(1:num_energies_neu,1:num_directions_neu) = gprad(1:num_energies_neu,1:num_directions_neu) &
      !        + (elrad(1:num_energies_neu,1:num_directions_neu,inode))/real(pnode)
        do itime = 1,ADR_neu % ntime
           gpunk(igaus,itime) = gpunk(igaus,itime) + gpsha(inode,igaus) * elunk(inode,itime) !> We add the contribution of 
                                                                                                 ! each element to the unkown for each time
        end do
     end do

     gpadv(1:ndime,igaus) = direc_neu(1:ndime,idire) !> The "velocity" is set equal to the direction we consider at the moment
    
     ! Resta unicamente la direccion y energia actuales
   !   Termino_resta=0.0
   !   do lcoef=0,num_legendre_neu
   !       do mcoef=-lcoef,Lcoef
   !          armonic = funharmonic(lcoef, mcoef, phi_neu(idire),tita_neu(idire),phi_neu(idire),tita_neu(idire)) 
   !          Termino_resta=Termino_resta + scatt_neu(pmate,enerG,enerG,lcoef+1) *armonic* weigd_neu(idire)*in4pi !&
                                       ! *scattering_neu(idire,idire) !*ener_weigd_neu(enerG)
   !       enddo
   !   enddo

     gprea(igaus)         = gpabs(igaus) !- (Termino_resta)

     gprhs(igaus)         = gpsou(igaus) !>The right-hand side of the equation is the source term

   !   do en = num_energies_neu,1,-1 
   !   do en = efectivos_neu_in(enerG),efectivos_neu_out(enerG),-1 
     do iener=n_efectivos_neu_sparse_acum(pmate, enerG+1)-1, n_efectivos_neu_sparse_acum(pmate,enerG), -1
         en = efectivos_neu_sparse(pmate, iener)
         do jdire = 1,num_directions_neu
            ! if (en/=enerG .or. idire/=jdire) then
               do lcoef=0,num_legendre_neu
                  do mcoef=-Lcoef,Lcoef
                     armonic = funharmonic(Lcoef,mcoef,phi_neu(idire),tita_neu(idire),phi_neu(jdire),tita_neu(jdire)) 
                     
                     !> We add to RHS the contributions from scattered rad from other directions
                     gprhs(igaus) = gprhs(igaus) + &
                                       scatt_neu(pmate,en,enerG,Lcoef+1) * armonic * &
                                       weigd_neu(jdire) * gprad(en,jdire) * in4pi !*ener_weigd_neu(en) 
                  enddo
               enddo
            ! endif
         enddo
     end do

!     do jdire = 1,num_directions_neu
!        if( jdire /= idire ) then
            !> We add to RHS the contributions from scattered rad from other directions
!           gprhs(igaus) = gprhs(igaus) + gpsca(igaus) * in4pi * weigd_neu(jdire) * scattering_neu(idire,jdire) * gprad(1,jdire) 
!        end if
!     end do



  end do

end subroutine neu_coefficients

! real(8) function Armo_source_enmater_neu(xmed,ymed,ang_neu,st,sc)
real(rp) function Armo_source_enmater_neu(xmed,ang_neu,st,sc)
use def_kintyp_basic, only : rp
implicit none
real(rp) :: xmed,ang_neu,st,sc!,ymed

   Armo_source_enmater_neu = ang_neu*(10.0_rp-2*xmed) - st * (xmed-11.0_rp) * (xmed+1.0_rp) + sc * (xmed-11.0_rp)*(xmed+1.0_rp)
    
end function Armo_source_enmater_neu


real(rp) function funharmonic(LL,m,phi_neu,tita_neu,phi_neu2,tita_neu2)
use def_parame, only : pi
use def_kintyp, only     :  rp,ip
implicit none
integer(ip) :: LL,m
real(rp) :: phi_neu,tita_neu,phi_neu2,tita_neu2
!
integer(ip) :: ml,ml_m     
real(rp) :: funtri,funtri_p,norma,Ylm,Ylm_p, cos_tita, cos_tita2
! real(8) :: PI=3.141529
integer(ip), external :: nfac
real(rp), external :: FUNLEG

      ! ml= m - LL -1 ! ii-L-1
      ml=m

      cos_tita = cos(tita_neu)
      cos_tita2 = cos(tita_neu2)
      ml_m=abs(ml)
      ! if(ml>=0) then
      if(ml>0) then
           funtri   = sqrt(2._rp)*cos(real(ml, rp)*phi_neu)
           funtri_p = sqrt(2._rp)*cos(real(ml, rp)*phi_neu2)
      elseif (ml==0) then
           funtri   = 1.0_rp
           funtri_p = 1.0_rp
      else
           funtri   = sqrt(2._rp)*sin(real(ml_m, rp)*phi_neu)
           funtri_p = sqrt(2._rp)*sin(real(ml_m, rp)*phi_neu2)
      endif

      Norma = sqrt( (real(nfac(LL-ml_m), rp)/real(nfac(LL+ml_m), rp) ) ) ! *sqrt((2*LL+1)/(4*PI)) 

      Ylm   = Norma*FUNLEG(cos_tita ,LL,ml_m)* funtri 
      Ylm_p = Norma*FUNLEG(cos_tita2,LL,ml_m)* funtri_p 
       
       funharmonic = Ylm*Ylm_p

end function funharmonic

real(rp) FUNCTION FUNLEG(X,N,M)
!C
!C       ========================================================
!C       Purpose: Compute associated Legendre functions Pmn(x)
!C                and Pmn'(x) for a given order
!C       Input :  x --- Argument of Pmn(x)
!C                m --- Order of Pmn(x),  m = 0,1,2,...,n
!C                n --- Degree of Pmn(x), n = 0,1,2,...,N
!C       Output:  PM(n) --- Pmn(x)
!C                PD(n) --- Pmn'(x)
!C       ========================================================
!C
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
use def_kintyp_basic, only : rp, ip
implicit none
real(rp) :: X
integer(ip) ::n,m
!local
integer(ip) :: k
real(rp) :: X0,PM0,PMK,PM1,PM2
        
     FUNLEG = 0.0_rp

     IF (DABS(X).EQ.1.0D0) THEN
              IF (M.EQ.0) THEN
                 FUNLEG=1.0D0
                 IF (X.LT.0.0_rp) THEN
                    FUNLEG=(-1.0_rp)**M*FUNLEG
                 ENDIF
              ENDIF
           RETURN
        ENDIF
        X0=DABS(1.0D0-X*X)
        PM0=1.0D0
        PMK=PM0
        DO K=1,M
           PMK=(2.0D0*real(K,rp)-1.0D0)*DSQRT(X0)*PM0
           PM0=PMK
        ENDDO
        PM1=(2.0D0*real(M,rp)+1.0D0)*X*PM0
        
     IF(N.EQ.M) THEN 
       FUNLEG= PMK
        ELSEIF(N.EQ.M+1) THEN
          FUNLEG=PM1
        ELSE

          DO K=M+2,N
             PM2=((2.0D0*real(K,rp)-1.0D0)*X*PM1-(real(K+M,rp)-1.0D0)*PMK)/real(K-M,rp)
             FUNLEG=PM2
             PMK=PM1
             PM1=PM2
          ENDDO

     ENDIF
        
end function funleg

INTEGER(ip) FUNCTION NFAC(N)
use def_kintyp, only     :  ip!,rp
implicit none
integer(ip) :: N
! local
integer(ip) :: kk

NFAC = 1

do KK=1,N
     NFAC=NFAC*KK
enddo

end function nfac
