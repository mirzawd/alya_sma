!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_exacso.f90
!> @author  Guillaume Houzeaux
!> @date    18/04/2013
!> @brief   Manufactured solutions
!> @details Impose the exact RHS using a manufactured solution
!>          \verbatim
!>          - ITASK = 1 ... Compute exact solution and derivatives
!>                  = 2 ... Compute RHS at Gauss point
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------

subroutine sld_exacso(itask,gpcod,u,dudx,d2udx2,Cijkl,gprhs)

  use def_parame 
  use def_master, only       :  cutim
  use def_domain, only       :  ndime
  use def_solidz, only       :  kfl_exacs_sld
  implicit none

  integer(ip), intent(in)    :: itask                                !< Task
  real(rp),    intent(in)    :: gpcod(ndime)                         !< Coordinate
  real(rp),    intent(in)    :: Cijkl(ndime,ndime,ndime,ndime)       !< Cijkl
  real(rp),    intent(out)   :: u(ndime)                             !< Displacement: u(i) = ui
  real(rp),    intent(out)   :: dudx(ndime,ndime)                    !< dudx(i,j) = duj/dxi
  real(rp),    intent(out)   :: d2udx2(ndime,ndime,ndime)            !< d2udx2(i,j,k) = d^2uk/dxidxj
  real(rp),    intent(inout) :: gprhs(ndime)                         !< Force term to be assembled
  integer(ip)                :: ii,jj,kk,ll,mm,nn
  real(rp)                   :: tkron(3,3),x,y,z
  real(rp)                   :: Fnm,Fnl,Fik,term1,term2
  real(rp)                   :: term3,Eml,d2udt2(3),rho,dPijdxj
  real(rp)                   :: t
  !
  ! Initializations
  !
  ! u(i)          = ui
  ! dudx(i,j)     = duj/dxi
  ! d2udx2(i,j,k) = d^2uk/dxidxj
  ! d2udt2(i)     = d^2ui/dt^2
  ! tkron(i,j)    = delta_ij
  !
  x      = gpcod(1)
  y      = gpcod(2)
  z      = gpcod(ndime)
  t      = cutim
  u      = 0.0_rp
  rho    = 1.0_rp
  dudx   = 0.0_rp
  d2udx2 = 0.0_rp
  d2udt2 = 0.0_rp
  tkron  = 0.0_rp
  do ii = 1,ndime
     tkron(ii,ii) = 1.0_rp
  end do
  !
  ! Obtain unknowns and derivatives according to the exact solution
  !
  if( kfl_exacs_sld == 1 ) then
     !
     ! ux = a
     ! uy = b
     ! 
     u(1) = 1.0_rp * 1.0e-4_rp
     u(2) = 2.0_rp * 1.0e-4_rp

  else if( kfl_exacs_sld == 2 ) then
     !
     ! ux =  x+2y
     ! uy = 3x+4y
     ! 
     u(1)      = ( 1.0_rp*x + 2.0_rp*y ) * 1.0e-4_rp
     u(2)      = ( 3.0_rp*x + 4.0_rp*y ) * 1.0e-4_rp
     dudx(1,1) = ( 1.0_rp              ) * 1.0e-4_rp
     dudx(2,1) = ( 2.0_rp              ) * 1.0e-4_rp
     dudx(1,2) = ( 3.0_rp              ) * 1.0e-4_rp
     dudx(2,2) = ( 4.0_rp              ) * 1.0e-4_rp

  else if( kfl_exacs_sld == 3 ) then
     !
     ! ux = x^3 * y^4
     ! uy = x^3 * y^3
     !
     u(1)          =  x**3        * y**4
     dudx(1,1)     =  3.0_rp*x**2 * y**4
     dudx(2,1)     =  4.0_rp*x**3 * y**3
     d2udx2(1,1,1) =  6.0_rp*x    * y**4
     d2udx2(2,2,1) = 12.0_rp*x**3 * y**2
     d2udx2(1,2,1) = 12.0_rp*x**2 * y**3    
     d2udx2(2,1,1) =  d2udx2(1,2,1)             

     u(2)          =  x**3        * y**3
     dudx(1,2)     =  3.0_rp*x**2 * y**3
     dudx(2,2)     =  3.0_rp*x**3 * y**2
     d2udx2(1,1,2) =  6.0_rp*x    * y**3
     d2udx2(2,2,2) =  6.0_rp*x**3 * y
     d2udx2(1,2,2) =  9.0_rp*x**2 * y**2    
     d2udx2(2,1,2) =  d2udx2(1,2,2)             

  else if( kfl_exacs_sld == 4 ) then
     !
     ! ux = ( x+2y  ) * 10-4
     ! uy = ( 3x+4y ) * 10-4
     ! 
     u(1)      = 1.0e-4_rp*x + 2.0e-4_rp*y
     u(2)      = 3.0e-4_rp*x + 4.0e-4_rp*y
     dudx(1,1) = 1.0e-4_rp
     dudx(2,1) = 2.0e-4_rp
     dudx(1,2) = 3.0e-4_rp
     dudx(2,2) = 4.0e-4_rp

  else if( kfl_exacs_sld == 5 ) then
     !
     ! ux = (  x+2y ) * 3t
     ! uy = ( 3x+4y ) * 4t
     ! 
     u(1)      = ( 1.0_rp*x + 2.0_rp*y ) * 3.0_rp * t
     u(2)      = ( 3.0_rp*x + 4.0_rp*y ) * 4.0_rp * t
     dudx(1,1) = 1.0_rp * 3.0_rp * t
     dudx(2,1) = 2.0_rp * 3.0_rp * t
     dudx(1,2) = 3.0_rp * 4.0_rp * t
     dudx(2,2) = 4.0_rp * 4.0_rp * t

     d2udt2(1) = ( 1.0_rp*x + 2.0_rp*y ) * 3.0_rp 
     d2udt2(2) = ( 3.0_rp*x + 4.0_rp*y ) * 4.0_rp

  end if

  !----------------------------------------------------------------------
  !
  ! Compute RHS
  ! ELRHS <=  ELRHS - GPRHS * GPSHA * GPVOL (in sld_elmexa)
  !
  !----------------------------------------------------------------------

  if( itask == 2 ) then
     !
     ! We solve = rho * d^2u/dt^2 - dPij/dXj = rho * d^2u_exa/dt^2 - dPij_exa/dXj
     ! GPRHS = - rho * d^2u_exa/dt^2 + dPij_exa/dXj 
     !
     ! Xi  := reference coordinates
     ! xi  := current coordinates
     ! ui  := xi - Xi = displacement
     !
     ! Pij := Fik Skj (1st Piola-Kirchhoff stress tensor)
     ! Skj := Ckjml Eml
     ! Fik := dui/dxk + tkron_ik := dxi/dXk
     ! Eml := 1/2( Fnm Fnl - tkron_ml )
     !     := 1/2( dum/dXl + dul/dXm + duk/dXl duk/dXm )
     !
     ! Pij  = Fik Ckjml Eml
     !
     ! dPij/dXj = Ckjml( dFik/dXj Eml + Fik dEml/dXj )
     !          = Ckjml( d^2 ui/dXkdXj Eml + Fik/2( dFnm/dXj Fnl + Fnm dFnl/dXj ) ]
     !          = Ckjml( d^2 ui/dXkdXj Eml + Fik/2 ( Fnl d^2un/dXmdXj + Fnm d^2un/dXldXj ) ]
     !          = term2 + term1    
     !
     do ii = 1,ndime

        term1 = 0.0_rp
        term2 = 0.0_rp

        do jj = 1,ndime
           do kk = 1,ndime
              do ll = 1,ndime
                 do mm = 1,ndime
                    Eml   = 0.0_rp
                    term3 = 0.0_rp
                    do nn = 1,ndime
                       Fik   = dudx(kk,ii) + tkron(ii,kk)
                       Fnm   = dudx(mm,nn) + tkron(nn,mm)
                       Fnl   = dudx(ll,nn) + tkron(nn,ll)
                       Eml   = Eml + Fnm*Fnl 
                       term3 = term3 + Fnl*d2udx2(jj,mm,nn) + Fnm*d2udx2(jj,ll,nn)
                    end do
                    Eml   = 0.5_rp * ( Eml - tkron(mm,ll) ) 
                    term1 = term1 + cijkl(kk,jj,mm,ll) * 0.5_rp * Fik * term3
                    term2 = term2 + cijkl(kk,jj,mm,ll) * Eml * d2udx2(jj,kk,ii)
                 end do
              end do
           end do
        end do

        dPijdxj   =   term1 + term2
        gprhs(ii) =  -rho * d2udt2(ii) + dPijdxj

     end do

  end if

end subroutine sld_exacso
