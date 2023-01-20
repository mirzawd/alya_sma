!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_solite
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_solite
  ! NAME 
  !    tur_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the turbulence equations.
  ! USES
  !    tur_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    tur_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_gradie
  use def_kermod, only : kfl_adj_prob
  
  implicit none
  real(rp)          :: cpu_refe1,cpu_refe2
  real(rp), pointer :: untu3(:)

  call cputim(cpu_refe1)
  !
  ! Set up the solver parameters for the turbulence equations
  !
  call tur_inisol(1_ip)
  !
  ! k-eps-phi-f: grad(phi) needed to compute Laplacian
  ! 
  if( kfl_grphi_tur==1 .and. iunkn_tur==4 .and. INOTMASTER ) then 
     untu3 => untur(3,:,1)
     call gradie(untu3,grphi_tur)
  end if
  !
  ! Update inner iteration counter and write headings in solver file
  !
  if(kfl_algor_tur==1) then
     if(iunkn_tur==1.and.itera_tur==1) then
        itinn(modul) = itinn(modul) + 1
        ittot_tur    = ittot_tur + 1
     end if
  else
     itinn(modul) = itinn(modul) + 1
     ittot_tur    = ittot_tur + 1
  end if
  !
  ! Update boundary conditions
  !
  call tur_updbcs(TUR_BEFORE_INNER_ITERATION)
  !
  ! Construct the system matrix and right-hand-side
  !
  call tur_matrix()
  !
  ! Solve the algebraic system
  !
  call tur_updunk(ITASK_INNITE)                       ! Initial condition for UNKNO 
!  if( INOTMASTER )     call tur_relat(amatr, rhsid)
  call solver(rhsid,unkno,amatr,pmatr)          ! Solve system for UNKNO 
  if (kfl_adj_prob == 0) call tur_clippi()

  call cputim(cpu_refe2)
  cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 

end subroutine tur_solite
!
! ****************************************************************************
! ****************************************************************************
!
subroutine tur_relat(Auu, bu)
  use def_kintyp
  use def_domain
  use def_master
  use def_turbul

  implicit none
  real(rp),    intent(inout)   :: Auu(nzdom)
  real(rp),    intent(inout)   :: bu( npoin)
  integer(ip)                :: ipoin,jpoin,izdom
  real(rp)                   :: xdiag, facto
  
  facto = max(1.05_rp**(-ittim), 0.2_rp)
!  facto = 1.0_rp
  print *, facto, 'turbu'
  do ipoin = 1,npoin

     izdom = r_dom(ipoin) - 1
     jpoin = 0
     do while( jpoin /= ipoin )
        izdom = izdom + 1
        jpoin = c_dom(izdom)
     end do 
     
     xdiag = Auu(izdom)*facto
     Auu(izdom) =  Auu(izdom) + xdiag
     bu(ipoin)  =  bu(ipoin)  + xdiag*untur(iunkn_tur, ipoin,3) !*unkno(ipoin)


  end do
end subroutine tur_relat
