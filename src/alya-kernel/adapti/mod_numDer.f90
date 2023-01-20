!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!***************************************************************
!*
!*              Module for numerical differentiation
!* 
!***************************************************************
MODULE mod_numDer

  use def_kintyp_basic,       only: ip,rp,lg
  
  IMPLICIT NONE
  SAVE
  
  public :: numericalDerivative_centered,numericalHessian_centered
  public :: numericalDerivative_forward, numericalHessian_forward
  
CONTAINS
  !
  !
  !
  function numericalDerivative_centered(fun,dim,x) result(df)
    !***************************************************
    !*
    !*   Computes the distortion and its derivative on the objective node of the SRF mesh
    !*
    !***************************************************
    implicit none
    !
    ! input-output variables
    interface
      function fun(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function fun
    end interface
    
    integer(ip),intent(in)  :: dim
    real(rp),intent(in)   :: x(dim)
    real(rp) :: df(dim)
    
!    real(rp) :: df_bis(dim)
    
    integer(ip) :: idim
    real(rp)    :: eps3, typ, ei, temp
    real (rp)   :: x_min(dim),x_plus(dim)
    real(rp)    :: f_plus, f_min

    eps3 = 4.64e-6_rp ! 1.0e-4 ! 4.64e-6
    typ  = 0.5_rp

    do idim=1,dim

      ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
      temp = x(idim) +  ei
      ei   = temp -  x(idim)

      x_min=x
      x_min(idim)=x_min(idim)-ei
      f_min= fun(dim,x_min)

      x_plus=x
      x_plus(idim)=x_plus(idim)+ei
      f_plus= fun(dim,x_plus)

      df(idim) = (f_plus-f_min)/(2.0_rp*ei)

    end do
    
!     df_bis = numericalDerivative_harcoded23(fun,dim,x)
!     if( norm2(df-df_bis)>1e-12) then
!       print*,df
!       print*,df_bis
!       call runend("error in df numer")
!     end if

  end function numericalDerivative_centered
  !
  !
  !
  function numericalHessian_centered(fun,dim,x) result(H)
    !
    implicit none
    !
    ! input-output variables
    interface
      function fun(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function fun
    end interface

    integer(ip),intent(in)  :: dim
    real(rp),   intent(in)   :: x(dim)

    real(rp) :: H(dim,dim)
!    real(rp) :: H_bis(dim,dim)

    integer(ip) :: idim, jdim
    real(rp)    :: eps3, typ, ei, ej, temp
    real (rp)   :: x_aux(dim)
    real(rp)    :: fPiPj,fPiMj,fMiPj,fMiMj,f
!    real(rp) :: kkaux

    eps3 = 1.0e-4_rp ! 4.64e-6!
    typ  = 0.5_rp

    f = fun(dim,x)

    do idim=1,int(size(x),ip)

      ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
      temp = x(idim) +  ei
      ei   = temp -  x(idim)
      
      do jdim=idim,int(size(x),ip)
        
        ej   = eps3*max( abs(x(jdim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
        temp = x(jdim) +  ej
        ej   = temp -  x(jdim)
       
        x_aux       = x
        x_aux(idim) = x_aux(idim)+ei
        x_aux(jdim) = x_aux(jdim)+ej
        fPiPj       = fun(dim,x_aux)

        x_aux       = x
        x_aux(idim) = x_aux(idim)-ei
        x_aux(jdim) = x_aux(jdim)-ej
        fMiMj       = fun(dim,x_aux)
        
        if(idim.eq.jdim) then
          fPiMj = f
          fMiPj = f
        else
          x_aux       = x
          x_aux(idim) = x_aux(idim)+ei
          x_aux(jdim) = x_aux(jdim)-ej
          fPiMj       = fun(dim,x_aux)

          x_aux       = x
          x_aux(idim) = x_aux(idim)-ei
          x_aux(jdim) = x_aux(jdim)+ej
          fMiPj       = fun(dim,x_aux)
        end if

        H(idim,jdim) = ((fPiPj+fMiMj)-(fPiMj+fMiPj))/(4_rp*ei*ej)
        H(jdim,idim) = H(idim,jdim)
      end do
    end do
    
!     ! TODO: remove this check
!     H_bis = numericalHessian_harcoded23(fun,dim,x)
!     if( norm2(H-H_bis)>1e-12_rp) then
!       call runend("error in H numer")
!     end if

    return
  end function numericalHessian_centered
  !
  !
  !
  subroutine numericalDerivative_forward(fun,dim,x,f,df)
    !
    use mod_debugTools, only: out_performance, deb_numer
    implicit none
    !
    interface
      function fun(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function fun
    end interface
    
    integer(ip),intent(in)   :: dim
    real(rp),   intent(in)   :: x(dim)
    real(rp),   intent(out)  :: df(dim), f
    
    real(rp) :: df_bis(dim)
    
    integer(ip) :: idim
    real(rp)    :: eps3, typ, ei, temp
    real (rp)   :: x_plus(dim)
    real(rp)    :: f_plus
    real(rp) :: t0,t1
    !
    if(out_performance) call cpu_time(t0)

    eps3 = 4.64e-6_rp ! 1.0e-4 ! 4.64e-6
    typ  = 0.5_rp

    f = fun(dim,x)
    
    do idim=1,dim

      ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
      temp = x(idim) +  ei
      ei   = temp -  x(idim)

      x_plus=x
      x_plus(idim)=x_plus(idim)+ei
      f_plus= fun(dim,x_plus)

      df(idim) = (f_plus-f)/ei
    end do
  
    if(out_performance) then
      call cpu_time(t1)
      deb_numer = deb_numer + (t1-t0)
    end if
    
  end subroutine numericalDerivative_forward
  !
  !
  !
  subroutine numericalHessian_forward(fun,dim,x,f,df,H)
    !
    use mod_debugTools, only: out_performance, deb_numer
    implicit none
    !
    interface
      function fun(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function fun
    end interface

    integer(ip),intent(in)  :: dim
    real(rp),   intent(in)  :: x(dim)
    real(rp),   intent(out) :: H(dim,dim), df(dim), f
    !
    integer(ip) :: idim, jdim
    real(rp)    :: eps3, typ, ei, ej, temp
    real (rp)   :: x_aux(dim),x_plus(dim)
    real(rp)    :: fPiPj
    real(rp) :: kkaux
    real(rp) :: fp(dim)
    real(rp) :: t0,t1
    !
    if(out_performance) call cpu_time(t0)

    eps3 = 1.0e-4_rp ! 4.64e-6!
    typ  = 0.5_rp

    f = fun(dim,x)

!     do idim=1,dim
!       ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
!       temp = x(idim) +  ei
!       ei   = temp -  x(idim)
!
!       x_plus=x
!       x_plus(idim)=x_plus(idim)+ei
!       fp(idim)= fun(dim,x_plus)
!
!       df(idim) = (fp(idim)-f)/ei
!     end do

    do idim=1,int(size(x),ip)

      ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0_rp,x(1))  ! sign is not necessary (centered)
      temp = x(idim) +  ei
      ei   = temp -  x(idim)
      
      do jdim=idim,int(size(x),ip)
        
        ej   = eps3*max( abs(x(jdim)),  typ )!*sign(1.0_rp,x(1))  ! sign is not necessary (centered)
        temp = x(jdim) +  ej
        ej   = temp -  x(jdim)
        
        x_aux       = x
        x_aux(jdim) = x_aux(jdim)+ej
        if(idim==1) then ! this if can be removed if previous loop is uncommented
          fp(jdim)= fun(dim,x_aux)
          df(jdim) = (fp(jdim)-f)/ej
        end if
        x_aux(idim) = x_aux(idim)+ei
        fPiPj       = fun(dim,x_aux)

        H(idim,jdim) = ( fPiPj - fp(idim) - fp(jdim) + f )/( ei*ej )
        H(jdim,idim) = H(idim,jdim)
      end do
    end do
    
  
    if(out_performance) then
      call cpu_time(t1)
      deb_numer = deb_numer + (t1-t0)
    end if
    
    
    return
  end subroutine numericalHessian_forward
  !
  !
  !
  function numericalDerivative_harcoded23(fun,dim,x) result(df)
    !***************************************************
    !*
    !*   Computes the distortion and its derivative on the objective node of the SRF mesh
    !*
    !***************************************************
    implicit none
    !
    ! input-output variables
    integer(ip),intent(in)  :: dim
    real(rp),intent(in)   :: x(dim)

    interface
      function fun(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function fun
    end interface

    real(rp) :: df(dim)
    !
    ! inside variables
    !real(rp)    :: f
    real (rp)  :: typ, ex, ey, ez, temp, eps3
    real (rp)  :: x_min(dim),x_plus(dim), f_plus, f_min,fx,fy,fz

    if(dim>3) then
      call runend('numerical derivative implmeneted only for dim<=3')
    end if

    !f = fun(dim,x)

    eps3 = 4.64e-6_rp ! 1.0e-4 ! 4.64e-6
    typ  = 0.5_rp
    ex   = eps3*max( abs(x(1)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
    temp = x(1) +  ex
    ex   = temp -  x(1)
    ey   = eps3*max( abs(x(2)),  typ)!*sign(1.0,x(2))  ! sign is not necessary (centered)
    temp = x(2) +  ey
    ey   = temp - x(2)
    if(dim>2) then
      ez   = eps3*max(abs(x(3)),   typ )!*sign(1.0,x(3))  ! sign is not necessary (centered)
      temp = x(3) +  ez
      ez   = temp - x(3)
    end if

    ! Derivative with respect x
    x_min   =x
    x_min(1)  =x_min(1)-ex
    f_min   = fun(dim,x_min)

    x_plus  =x
    x_plus(1)  =x_plus(1)+ex
    f_plus   = fun(dim,x_plus)

    fx=(f_plus-f_min)/(2*ex)

!     print*,"ex: ",ex
!     print*,"f_plus: ",f_plus
!     print*,"f_min: ",f_min

    ! Derivative with respect y
    x_min   =x
    x_min(2)  =x_min(2)-ey
    f_min   = fun(dim,x_min)

    x_plus  =x
    x_plus(2)  =x_plus(2)+ey
    f_plus   = fun(dim,x_plus)

     fy=(f_plus-f_min)/(2*ey)

!      print*,"f_plus: ",f_plus
!      print*,"f_min: ",f_min

     if(dim>2) then
       ! Derivative with respect z
       x_min   =x
       x_min(3)  =x_min(3)-ez
       f_min   = fun(dim,x_min)

       x_plus  =x
       x_plus(3)  =x_plus(3)+ez
       f_plus   = fun(dim,x_plus)

        fz=(f_plus-f_min)/(2*ez)

       ! Compose the differential vector
       df(:)=(/ fx, fy, fz /);
     else
       df(:)=(/fx,fy/)
     end if
    return
  end function numericalDerivative_harcoded23
  !
  !
  !
  function numericalHessian_harcoded23(fun,dim,x) result(H)
    !***************************************************
    !*
    !*   Computes the distortion and its derivative on the objective node of the SRF mesh
    !*
    !***************************************************
    implicit none
    !
    ! input-output variables
    integer(ip),intent(in)  :: dim
    real(rp),intent(in)   :: x(dim)

    interface
      function fun(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function fun
    end interface

    real(rp) :: H(dim,dim)
    !
    ! inside variables
    real (rp)  :: f
    real (rp)  :: typ, ex, ey, ez, temp, eps3
    real (rp)  :: x_aux(dim)
    real (rp)  :: fminx,fplusx,fminy,fplusy,fplusz,fminz
    real (rp)  :: fMxPy,fPxMy,fMxPz,fPxMz,fMyPz,fPyMz
    real (rp)  :: fPxPy,fMxMy,fPxPz,fMxMz,fMyMz,fPyPz

    if(dim>3) then
      call runend('numerical derivative implmeneted only for dim<=3')
    end if

    f = fun(dim,x)

    eps3 = 1.0e-4_rp !4.64e-6 !1.0e-4
    typ  = 0.5_rp
    ex   = eps3*max( abs(x(1)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
    temp = x(1) +  ex
    ex   = temp -  x(1)
    ey   = eps3*max( abs(x(2)),  typ)!*sign(1.0,x(2))  ! sign is not necessary (centered)
    temp = x(2) +  ey
    ey   = temp - x(2)
    if(dim>2) then
      ez   = eps3*max( abs(x(3)),   typ )!*sign(1.0,x(3))  ! sign is not necessary (centered)
      temp = x(3) +  ez
      ez   = temp - x(3)
    end if
    
    x_aux   =x
    x_aux(1)  =x(1)-(2_rp*ex)
    fminx   = fun(dim,x_aux)
    x_aux(1)  =x(1)+(2_rp*ex)
    fplusx   = fun(dim,x_aux)

    x_aux   =x
    x_aux(2)  =x(2)-(2_rp*ey)
    fminy   = fun(dim,x_aux)
    x_aux(2)  =x(2)+(2_rp*ey)
    fplusy  = fun(dim,x_aux)

    x_aux   =x
    x_aux(1)  =x(1)+ex
    x_aux(2)  =x(2)-ey
    fPxMy   = fun(dim,x_aux)
    x_aux(1)  =x(1)-ex
    x_aux(2)  =x(2)+ey
    fMxPy   = fun(dim,x_aux)
    x_aux(1)  =x(1)+ex
    x_aux(2)  =x(2)+ey
    fPxPy   = fun(dim,x_aux)
    x_aux(1)  =x(1)-ex
    x_aux(2)  =x(2)-ey
    fMxMy   = fun(dim,x_aux)

    !H(1,1) = (fplusx-2*f+fminx)/(4*ex*ex)
    !H(2,2) = (fplusy-2*f+fminy)/(4*ey*ey)
    !H(1,2) = (fPxPy-fPxMy-fMxPy+fMxMy)/(4*ex*ey)
    !H(2,1) = H(1,2)
    H(1,1) = ((fplusx+fminx)-(2_rp*f))/(4_rp*ex*ex)
    H(2,2) = ((fplusy+fminy)-(2_rp*f))/(4_rp*ey*ey)
    H(1,2) = ((fPxPy+fMxMy)-(fPxMy+fMxPy))/(4_rp*ex*ey)
    H(2,1) = H(1,2)
    
    
!     print*,'in num der'
!     print*,"(fplusy+fminy):         ",(fplusy+fminy)
!     print*,"(fplusy):         ",(fplusy)
!     print*,"(fminy):         ",(fminy)
    !print*,"2*f:                    ",2_rp*f
    !print*,"((fplusy+fminy)-(2*f)): ",((fplusy+fminy)-(2_rp*f))
!     print*,"(4*ey*ey):              ",(4*ey*ey)

!     print*,"fplusx: ",fplusx
!     print*,"f:      ",f
!     print*,"fminx:  ",fminx
!     print*,4*ex*ex
!     print*,"fplusx+fminx: ",fplusx+fminx
!     print*,"-2*f: ",-2*f
!     print*,"2tens: ",fplusx-2*f
!     print*,"3tens: ",fplusx-f-f
!     print*,(fplusx-f)
!     print*,(fplusx-2_rp*f)
!     print*,(fplusx-2_rp*f+fminx)
!     print*,(fplusx+fminx)-(2_rp*f)
!     print*,(4_rp*ex*ex)
!     print*,(fplusx-f-f+fminx)/(4_rp*ex*ex)
!     print*,(fplusx-2_rp*f+fminx)/(4_rp*ex*ex)
!     print*,(fplusx-2*f+fminx)/(4*ex*ex)

    if(dim>2) then
      x_aux   =x
      x_aux(3)  =x(3)-2*ez
      fminz   = fun(dim,x_aux)
      x_aux(3)  =x(3)+2*ez
      fplusz  = fun(dim,x_aux)

      x_aux   =x
      x_aux(1)  =x(1)+ex
      x_aux(3)  =x(3)-ez
      fPxMz  = fun(dim,x_aux)
      x_aux(1)  =x(1)-ex
      x_aux(3)  =x(3)+ez
      fMxPz  = fun(dim,x_aux)
      x_aux(1)  =x(1)+ex
      x_aux(3)  =x(3)+ez
      fPxPz   = fun(dim,x_aux)
      x_aux(1)  =x(1)-ex
      x_aux(3)  =x(3)-ez
      fMxMz   = fun(dim,x_aux)

      x_aux   =x
      x_aux(2)  =x(2)+ey
      x_aux(3)  =x(3)-ez
      fPyMz  = fun(dim,x_aux)
      x_aux(2)  =x(2)-ey
      x_aux(3)  =x(3)+ez
      fMyPz   = fun(dim,x_aux)
      x_aux(2)  =x(2)+ey
      x_aux(3)  =x(3)+ez
      fPyPz   = fun(dim,x_aux)
      x_aux(2)  =x(2)-ey
      x_aux(3)  =x(3)-ez
      fMyMz   = fun(dim,x_aux)

      H(3,3) = ((fplusz+fminz)-(2*f))/(4*ez*ez)
      H(1,3) = ((fPxPz+fMxMz)-(fPxMz+fMxPz))/(4*ex*ez)
      H(2,3) = ((fPyPz+fMyMz)-(fPyMz+fMyPz))/(4*ey*ez)
!       H(1,3) = ((fPxPz+fMxMz)-(fPyMz+fMyPz))/(4*ex*ez)
!       H(2,3) = ((fPyPz+fMyMz)-(fPxMz+fMxPz))/(4*ey*ez)
      H(3,1) = H(1,3)
      H(3,2) = H(2,3)
    end if


    return
  end function numericalHessian_harcoded23
  !
  !
  !
END MODULE mod_numDer
