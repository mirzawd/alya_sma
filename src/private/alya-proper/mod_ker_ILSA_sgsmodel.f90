!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ker_ILSA_sgsmodel.f90
!> @author  Oriol Lehmkuhl
!> @date    04/01/2018
!> @brief   ILSA sgs model
!> @details ILSA sgs model
!> see Dynamic subfilter-scale stress model for large-eddy simulations,
!>     A. Rouhi, U. Piomelli and B.J. Geurts, PHYSICAL REVIEW FLUIDS 1, 044401 (2016) 
!> @} 
!-----------------------------------------------------------------------
module mod_ker_ILSA_sgsmodel
  
  use def_kintyp, only : ip,rp
  use def_master, only : mem_modul
  use def_master, only : modul
  use def_master, only : dtinv
  use def_domain, only : ndime,mnode
  use def_domain, only : nelem,mgaus
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo

  implicit none

  real(rp),   pointer, save :: ck_square(:,:)   ! ck^2
  real(rp),   pointer, save :: ave_kres(:,:)    ! time average resolved kinetic energy 
  real(rp),   pointer, save :: ave_etot(:,:)    ! time average total dissipation 

  real(rp),   pointer, save :: ave_x1(:,:)      ! time average x1 eq(11) 
  real(rp),   pointer, save :: ave_x2(:,:)      ! time average x2 eq(11) 
  real(rp),   pointer, save :: ave_x3(:,:)      ! time average x3 eq(11) 
  real(rp),   pointer, save :: ave_x4(:,:)      ! time average x4 eq(11) 

  real(rp),   pointer, save :: ave_sij(:,:,:,:) ! time average  sij
  real(rp),   pointer, save :: ave_vel(:,:,:)   ! time average  vel

  real(rp),   pointer, save :: gpvsfs(:,:)      ! sfs viscosity
  real(rp),   pointer, save :: gplest(:,:)      ! Lest
  real(rp),   pointer, save :: gpstrain(:,:)    ! |strain|
  real(rp),   pointer, save :: eliti(:)         ! local averaging time

  real(rp),   pointer, save :: delta_f(:,:)     ! Estimated delta of the filter
  real(rp),   pointer, save :: ratio_l_df(:,:)  ! ratio between L/Delta_f
  real(rp),   pointer, save :: ratio_df_h(:,:)   ! ratio between Delta_f/h
  real(rp),   pointer, save :: ave_grad(:,:)    ! average ratio between (Grad·Grad)
  real(rp),   pointer, save :: ave_rij(:,:)     ! average ratio between (Rij·Rij)
  real(rp),   pointer, save :: d_stau(:,:)       ! dynamic_stau

  private

  public :: ker_ILSA_sgs_viscosity
  public :: ker_ILSA_temporal_average
  public :: ker_ILSA_solve_ck
  public :: ker_ILSA_get_ratio_ldf
  public :: ker_ILSA_get_ratio_dfh
  public :: ker_ILSA_get_lest
  public :: ker_ILSA_get_deltaf
  public :: ker_ILSA_get_ck

  real(rp), save :: s_tau
  real(rp), save :: dist
!  real(rp), save :: dfilter
  real(rp), save :: T_ave ! a de venir de fora
  real(rp), save :: gzero

contains 


  subroutine ker_ILSA_allocate()

    integer(ip), save :: ipass = 0
    integer(ip)       :: igaus,idime,jdime,ielem

    !s_tau = 0.022_rp
    gzero = 1.0e-6_rp ! zeror
    !gzero = 1.0e-10_rp ! zeror
  

    if( ipass == 0 ) then

       ipass = 1

       nullify(ck_square)
       nullify(ave_kres)
       nullify(ave_etot)

       nullify(ave_x1)
       nullify(ave_x2)
       nullify(ave_x3)
       nullify(ave_x4)

       nullify(ave_sij)
       nullify(ave_vel)

       nullify(gpvsfs)
       nullify(gplest)
       nullify(gpstrain)

       nullify(delta_f)
       nullify(ratio_l_df)
       nullify(ratio_df_h)
       nullify(ave_grad)
       nullify(ave_rij)
       nullify(d_stau)

       call memory_alloca(mem_modul(1:2,modul),'CKSQUARE','ker_ILSA_allocate',ck_square,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEKRES' ,'ker_ILSA_allocate',ave_kres,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEETOT' ,'ker_ILSA_allocate',ave_etot,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX1'   ,'ker_ILSA_allocate',ave_x1,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX2'   ,'ker_ILSA_allocate',ave_x2,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX3'   ,'ker_ILSA_allocate',ave_x3,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX4'   ,'ker_ILSA_allocate',ave_x4,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVGRAD'  ,'ker_ILSA_allocate',ave_grad,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVRIJ'   ,'ker_ILSA_allocate',ave_rij,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'GPVSFS'  ,'ker_ILSA_allocate',gpvsfs,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'GPLEST'  ,'ker_ILSA_allocate',gplest,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'GPSTRAIN','ker_ILSA_allocate',gpstrain,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEVEL'  ,'ker_ILSA_allocate',ave_vel,nelem,mgaus,ndime)
       call memory_alloca(mem_modul(1:2,modul),'AVESIJ'  ,'ker_ILSA_allocate',ave_sij,nelem,mgaus,ndime,ndime)
       call memory_alloca(mem_modul(1:2,modul),'ELITI'   ,'ker_ILSA_allocate',eliti,nelem)
       call memory_alloca(mem_modul(1:2,modul),'DELTAF'  ,'ker_ILSA_allocate',delta_f,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'RATIOLDF','ker_ILSA_allocate',ratio_l_df,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'RATIODFH','ker_ILSA_allocate',ratio_df_h,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'DSTAU'   ,'ker_ILSA_allocate',d_stau,nelem,mgaus)
       
       do ielem = 1,nelem
          eliti    (ielem) = 0.0_rp
          do igaus = 1,mgaus
             ck_square (ielem,igaus) = 0.0_rp
             ave_kres  (ielem,igaus) = 0.0_rp
             ave_etot  (ielem,igaus) = 0.0_rp  
             ave_x1    (ielem,igaus) = 0.0_rp  
             ave_x2    (ielem,igaus) = 0.0_rp 
             ave_x3    (ielem,igaus) = 0.0_rp   
             ave_x4    (ielem,igaus) = 0.0_rp    
             gpvsfs    (ielem,igaus) = 0.0_rp
             gplest    (ielem,igaus) = 0.0_rp
             gpstrain  (ielem,igaus) = 0.0_rp
             delta_f   (ielem,igaus) = 0.0_rp
             ratio_l_df(ielem,igaus) = 0.0_rp
             ratio_df_h(ielem,igaus) = 0.0_rp
             ave_grad  (ielem,igaus) = 0.0_rp    
             ave_rij   (ielem,igaus) = 0.0_rp    
             d_stau    (ielem,igaus) = 0.0_rp
             do idime = 1,ndime
                ave_vel(ielem,igaus,idime) = 0.0_rp
                do jdime = 1,ndime
                   ave_sij(ielem,igaus,idime,jdime) = 0.0_rp
                end do
             end do
          end do
       end do

    end if

 end subroutine ker_ILSA_allocate

 subroutine ker_ILSA_sgs_viscosity(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem,const1,const2,const3,hleng,gpwal2)

     integer(ip),  intent(in)    :: pgaus,igaui,igauf
     integer(ip),  intent(in)    :: ielem
     real(rp),     intent(in)    :: gpvel(ndime,pgaus)
     real(rp),     intent(in)    :: gpgve(ndime,ndime,pgaus)
     real(rp),     intent(in)    :: gpvis(pgaus)
     real(rp),     intent(in)    :: gpwal2(pgaus)
     real(rp),     intent(inout) :: gpmut(pgaus)
     real(rp),     intent(in)  :: const1,const2,const3
     real(rp),     intent(in)  :: hleng(3)
     integer(ip)   :: igaus

     T_ave = const1 
     s_tau = const2 
     dist  = const3
     !dfilter = const2

     call ker_ILSA_allocate()
     call ker_ILSA_temporal_average(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem,hleng,gpwal2)
     call ker_ILSA_solve_ck(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem)

     ! loop para calcular la mu_sgs

     do igaus = igaui,igauf
        gpmut(igaus)        = ck_square(ielem,igaus)*(gplest(ielem,igaus)**2_ip)*gpstrain(ielem,igaus)
        gpvsfs(ielem,igaus) = gpmut(igaus)
     end do

  end subroutine ker_ILSA_sgs_viscosity   

  subroutine ker_ILSA_temporal_average(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem,hleng,gpwal2)

     integer(ip),  intent(in)    :: pgaus,igaui,igauf
     integer(ip),  intent(in)    :: ielem
     real(rp),     intent(in)    :: gpvel(ndime,pgaus)
     real(rp),     intent(in)    :: gpgve(ndime,ndime,pgaus)
     real(rp),     intent(in)    :: gpvis(pgaus)
     real(rp),     intent(in)    :: gpwal2(pgaus)
     real(rp),     intent(inout) :: gpmut(pgaus)
     real(rp),     intent(in)    :: hleng(3)
     integer(ip)                 :: igaus,idime,jdime,kdime
     real(rp)                    :: gpkres(pgaus),gpvelf(ndime,pgaus)
     real(rp)                    :: gpsij(ndime,ndime,pgaus),gpepst(pgaus),gprij(ndime,ndime,pgaus)
     real(rp)                    :: ave_epsilon, dt, gpsijf(ndime,ndime,pgaus),gpx1(pgaus),gpx2(pgaus),gpx3(pgaus)
!     real(rp)                    :: gpgrad(ndime,ndime,pgaus)
     real(rp)                    :: G__ij(3,3)
     real(rp)                    :: op,alpha,Bbeta

     dt = 1.0_rp/(dtinv+gzero)
     ave_epsilon = dt/(T_ave+dt+gzero)
     if(eliti(ielem) > T_ave) then
        eliti(ielem) = 0.0_rp
     else
        eliti(ielem) = eliti(ielem) + dt
     end if

     eliti(ielem) = 1.0_rp-ave_epsilon

     do igaus = igaui,igauf
     ! get u and sij at gauss points
     ! gpgvel evaluated on turbu no se si duplicar
        do idime = 1,ndime
           do jdime = 1,ndime
              gpsij(idime,jdime,igaus) = 0.5_rp*(gpgve(jdime,idime,igaus)+gpgve(idime,jdime,igaus))
           end do
        end do

        gpstrain(ielem,igaus) = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
            gpstrain(ielem,igaus) = gpstrain(ielem,igaus) + (gpsij(idime,jdime,igaus)*gpsij(idime,jdime,igaus))
           end do
        end do
        gpstrain(ielem,igaus) = sqrt(2.0_rp*gpstrain(ielem,igaus))

     ! u mean and Sij mean
        do idime = 1,ndime
           ave_vel(ielem,igaus,idime) = (ave_epsilon*gpvel(idime,igaus) + eliti(ielem)*ave_vel(ielem,igaus,idime))!/(dt+eliti(ielem))
           do jdime = 1,ndime
              ave_sij(ielem,igaus,idime,jdime) = (ave_epsilon*gpsij(idime,jdime,igaus) + eliti(ielem)*ave_sij(ielem,igaus,idime,jdime))!/(dt+eliti(ielem))
           end do
        end do

     ! k res, epsilon total and R^a_{ij}

        ! fluctuating fields
        gpkres(igaus) = 0.0_rp
        do idime = 1,ndime
           gpvelf(idime,igaus) = gpvel(idime,igaus)-ave_vel(ielem,igaus,idime)
           gpkres(igaus)       = gpkres(igaus) + 0.5_rp*gpvelf(idime,igaus)*gpvelf(idime,igaus)
           do jdime = 1,ndime
               gpsijf(idime,jdime,igaus) = gpsij(idime,jdime,igaus)-ave_sij(ielem,igaus,idime,jdime) 
           end do
        end do

        ! average k and epsilon
        ave_kres(ielem,igaus) = (ave_epsilon*gpkres(igaus) + eliti(ielem)*ave_kres(ielem,igaus))!/(dt+eliti(ielem))

        gpepst(igaus) = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              gpepst(igaus) = gpepst(igaus) + 2.0_rp*(gpvis(igaus)+gpvsfs(ielem,igaus))*gpsijf(idime,jdime,igaus)*gpsijf(idime,jdime,igaus) 
           end do
        end do
        ave_etot(ielem,igaus) = (ave_epsilon*gpepst(igaus) + eliti(ielem)*ave_etot(ielem,igaus))!/(dt+eliti(ielem))

        ! estimated integral length scale
        gplest(ielem,igaus) = sqrt(ave_kres(ielem,igaus)**3_ip)/(gzero+ave_etot(ielem,igaus))

        ! R^a_{ij}
        gprij(1:ndime,1:ndime,igaus) =  0.0_rp

        do idime = 1,ndime
           gprij(idime,idime,igaus) =  - (1.0_rp/3.0_rp)*gpvelf(idime,igaus)*gpvelf(idime,igaus)
           do jdime = 1,ndime
              gprij(idime,jdime,igaus) = gprij(idime,jdime,igaus) + gpvelf(idime,igaus)*gpvelf(jdime,igaus)
           end do
        end do

       G__ij = 0.0_rp ! G = g^T*g
       alpha = 0.0_rp
       do idime = 1_ip,3_ip
          do jdime = 1_ip,3_ip
             do kdime = 1_ip,3_ip
                G__ij(idime,jdime) = G__ij(idime,jdime) + gpgve(kdime,jdime,igaus)*gpgve(kdime,idime,igaus)
             end do
             alpha = alpha + gpgve(jdime,idime,igaus)*gpgve(idime,jdime,igaus)
          end do
       end do

       Bbeta =  G__ij(1_ip,1_ip)*G__ij(2_ip,2_ip) + G__ij(2_ip,2_ip)*G__ij(3_ip,3_ip) + G__ij(3_ip,3_ip)*G__ij(1_ip,1_ip) &
          - G__ij(1_ip,2_ip)*G__ij(1_ip,2_ip) - G__ij(2_ip,3_ip)*G__ij(2_ip,3_ip) - G__ij(1_ip,3_ip)*G__ij(1_ip,3_ip)

       if ( alpha > gzero )then                         ! Avoid divide by zero
          op = sqrt (max(( Bbeta ) / ( alpha ), 0.0_rp))
       else
          op = 0.0_rp
       end if

      !  op = gpstrain(ielem,igaus)

        !gpx1(igaus) = 2.0_rp*(gpstrain(ielem,igaus)**2_ip)
        gpx1(igaus) = 2.0_rp*(gplest(ielem,igaus)**4_ip)*(gpstrain(ielem,igaus)**2_ip)*(op**2_ip)
        gpx2(igaus) = 0.0_rp
        gpx3(igaus) = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              !gpx2(igaus) = gpx2(igaus) + 4.0_rp*(gpsij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
              gpx2(igaus) = gpx2(igaus) + 4.0_rp*(gplest(ielem,igaus)**2_ip)*op*(gpsij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
              gpx3(igaus) = gpx3(igaus) + (gprij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
           end do
        end do

        ave_x1(ielem,igaus) = (ave_epsilon*gpx1(igaus) + eliti(ielem)*ave_x1(ielem,igaus))!/(dt+eliti(ielem))
        ave_x2(ielem,igaus) = (ave_epsilon*gpx2(igaus) + eliti(ielem)*ave_x2(ielem,igaus))!/(dt+eliti(ielem))
        ave_x3(ielem,igaus) = (ave_epsilon*gpx3(igaus) + eliti(ielem)*ave_x3(ielem,igaus))!/(dt+eliti(ielem))

        !ave_x1(ielem,igaus) = gpx1(igaus)
        !ave_x2(ielem,igaus) = gpx2(igaus)
        !ave_x3(ielem,igaus) = gpx3(igaus)

        gpstrain(ielem,igaus) = op

        ! Estimated Delta_f and ratio L/Delta

        !gpgrad(1:ndime,1:ndime,igaus) =  0.0_rp
        !gprij(1:ndime,1:ndime,igaus) =  0.0_rp
        !do idime = 1,ndime
        !   do jdime = 1,ndime
        !         gpgrad(idime,jdime,igaus) = gpgrad(idime,jdime,igaus) + gpgve(jdime,idime,igaus)*gpgve(idime,jdime,igaus)
        !      gprij(idime,jdime,igaus) = gprij(idime,jdime,igaus) + gpvelf(idime,igaus)*gpvelf(jdime,igaus) - 2.0_rp*gpvsfs(ielem,igaus)*gpsij(idime,jdime,igaus)
        !   end do
        !end do

        !gpx2(igaus) = 0.0_rp
        !gpx3(igaus) = 0.0_rp
        !do idime = 1,ndime
        !   do jdime = 1,ndime
        !      gpx2(igaus) = gpx2(igaus) + (gpgrad(idime,jdime,igaus)*gpgrad(idime,jdime,igaus))
        !      gpx3(igaus) = gpx3(igaus) + (gprij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
        !   end do
        !end do
        !ave_grad(ielem,igaus) = (dt*(gpx2(igaus)) + eliti(ielem)*ave_grad(ielem,igaus))/(dt+eliti(ielem))
        !ave_rij(ielem,igaus)  = (dt*(gpx3(igaus)) + eliti(ielem)*ave_rij(ielem,igaus))/(dt+eliti(ielem))

        !if(ave_rij(ielem,igaus)>gzero) then
        !    d_stau(ielem,igaus) = min(2.0_rp*dfilter**2*sqrt(ave_grad(ielem,igaus)/ave_rij(ielem,igaus)), 0.16)
        !else
        !    d_stau(ielem,igaus) = 0.02_rp
        !end if
        if(gpwal2(igaus) <dist) then
           d_stau(ielem,igaus) = 0.99_rp-(0.99_rp-s_tau)*(gpwal2(igaus)/dist)
           !d_stau(ielem,igaus) = 0.99_rp-(0.99_rp-s_tau)*tanh(3.0_rp*gpwal2(igaus)/dist)
           !d_stau(ielem,igaus) = s_tau**(gpwal2(igaus)/dist)
           !d_stau(ielem,igaus) = 0.8_rp
           !d_stau(ielem,igaus) = 1.07_rp-(1.07_rp-s_tau)*cosh(PI*gpwal2(igaus)/(dist))/cos(PI)
           !d_stau(ielem,igaus) = 0.99_rp + ((s_tau-0.99_rp)/dist)*(gpwal2(igaus))
           !d_stau(ielem,igaus) = 0.98_rp 
        else
           d_stau(ielem,igaus) = s_tau
        end if

     end do
  end subroutine ker_ILSA_temporal_average  

  subroutine ker_ILSA_solve_ck(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem)

     integer(ip),  intent(in)  :: pgaus,igaui,igauf
     integer(ip),  intent(in)  :: ielem
     real(rp),     intent(in) :: gpvel(ndime,pgaus)
     real(rp),     intent(in) :: gpvis(pgaus)
     real(rp),     intent(in) :: gpgve(ndime,ndime,pgaus)
     real(rp),     intent(inout) :: gpmut(pgaus)
     real(rp) :: a,b,c,d,x1,x2
     integer(ip)   :: igaus

     a = 0.0_rp
     b = 0.0_rp
     c = 0.0_rp

     ! ax^2 + bx + c = 0

     do igaus = igaui,igauf

        ! d = ave_x2(ielem,igaus)**2_ip + 4.0_rp*((1.0/s_tau)**2_ip - 1.0_rp)*ave_x1(ielem,igaus)*ave_x3(ielem,igaus) 

        ! if(abs(ave_x1(ielem,igaus)) > gzero) then
        !    ck_square(ielem,igaus) = (sqrt(d)-ave_x2(ielem,igaus))/(2.0_rp*ave_x1(ielem,igaus)*((1.0_rp/s_tau)**2_ip - 1.0_rp))
        ! else
        !    if(abs(ave_x2(ielem,igaus))>gzero) then
        !       ck_square(ielem,igaus) = ave_x3(ielem,igaus)/ave_x2(ielem,igaus) 
        !    else
        !       ck_square(ielem,igaus) = 0.0_rp
        !    end if
        ! end if
        ! ck_square(ielem,igaus) = max(ck_square(ielem,igaus), 0.0_rp)

        a = ave_x1(ielem,igaus)*(1.0_rp - (1.0_rp/d_stau(ielem,igaus))**2_ip)
        b = -ave_x2(ielem,igaus)
        c = ave_x3(ielem,igaus)

        !if (abs(a)<gzero) then
        !   if(abs(c)<gzero) then
        !      x1 = 0.0_rp
        !      x2 = 0.0_rp
        !   else
        !      x1 = 2.0_rp*c/(-b + sqrt(b*b - 4.0_rp*a*c))
        !      x2 = 2.0_rp*c/(-b - sqrt(b*b - 4.0_rp*a*c))
        !   end if
        !else
        !   x1 = (-b + sqrt((b*b)-(4.0_rp*a*c))) / ((2.0_rp*a))
        !   x2 = (-b - sqrt((b*b)-(4.0_rp*a*c))) / ((2.0_rp*a))
        !end if
        !d = (b*b)/(b*b - 4.0_rp*a*c)
        !if(abs(d)<10.0_rp) then
        !   if(abs(x1)>abs(x2)) then
        !      x2 = (c/a)/x1
        !   else 
        !      x1 = (c/a)/x2
        !   end if
        !end if

        if(b> 0.0_rp) then
           d = -0.5_rp*(b + sqrt(b*b - 4.0_rp*a*c))
        else 
           d = -0.5_rp*(b - sqrt(b*b - 4.0_rp*a*c))
        end if

        if(abs(d)> gzero) then
           if(abs(a)>gzero) then 
              x1 = d / a
              x2 = c / d   
           else
              x1 = 0.0_rp
              x2 = c / d
           end if
        else
           x1 = 0.0_rp
           x2 = 0.0_rp 
        end if

        ck_square(ielem,igaus) = max(max(x1,x2), 0.0_rp)
        
     end do
     
  end subroutine ker_ILSA_solve_ck 

  pure subroutine ker_ILSA_get_ratio_ldf(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

    do igaus = igaui,igauf
       gpratio(igaus) = ave_etot(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_ratio_ldf 

 pure subroutine ker_ILSA_get_ratio_dfh(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

   do igaus = igaui,igauf
       gpratio(igaus) = ave_etot(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_ratio_dfh 

 pure subroutine ker_ILSA_get_lest(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

    do igaus = igaui,igauf
       gpratio(igaus) = gplest(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_lest 

 pure subroutine ker_ILSA_get_deltaf(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

   do igaus = igaui,igauf
       gpratio(igaus) = ave_kres(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_deltaf 

 pure subroutine ker_ILSA_get_ck(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

   do igaus = igaui,igauf
       gpratio(igaus) = sqrt(ck_square(ielem,igaus))
    end do
 end subroutine ker_ILSA_get_ck 

end module mod_ker_ILSA_sgsmodel
