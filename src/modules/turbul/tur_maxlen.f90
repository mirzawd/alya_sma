!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_maxlen()

  !*****************************************************
  !tur_maxlen:
  ! Used by tur_endste
  ! This routine calculates de maximum mixing length using Mellor Yamada
  ! Max mixing length = int (sqrt(k)z dz) /int (sqrt(k)dz) 
  !*****************************************************
  ! use of global variables and subroutines
  use def_parame
  use def_master,         only : rhsid, amatr,INOTMASTER, unkno
  use def_master,         only : solve, mem_modul, modul, untur
  use def_domain,         only : meshe, elmar, coord, npoin,ndime, walld
  use mod_ADR,            only : ADR_assemble_convective, ADR_assemble_extension
  use mod_memory,         only : memory_alloca, memory_deallo
  use def_solver,         only : solve_sol
  use mod_solver,         only : solver_solve
  use mod_solver,         only : solver_initialize_matrix_and_rhs
  use def_turbul,         only : kfl_fixno_tur, tur_max_mixlen
  use def_turbul,         only : numer, denom, z_wall, gravi_tur, kfl_cotem_tur
  use def_kermod,         only : ndivi
  use mod_communications, only : PAR_MAX, PAR_MIN
  use mod_ker_regularization, only : regul_k, kfl_regularization
  implicit none  
  ! local variables
  integer(ip)             :: ipoin, idime
  integer(ip), save       :: ipass=0
  real(rp),    save       :: z_top, z_min         ! maximum height
  integer(ip), pointer    :: kfl_fixno_tmp(:,:)
  real(rp),    pointer    :: bvess_tmp(:,:)
  real(rp) ,   pointer    :: conve_tmp(:,:)
  !  real(rp),    pointer    :: numer(:), denom(:), z_wall(:)
  real(rp),    pointer    :: force(:)

  integer(ip)            ::  vcoor , nloop, iloop
  vcoor = 3_ip ! vertical coordinate

  if (kfl_cotem_tur .gt. 0_ip) then !if temper coupling, set vertical
     !coordinate in terms of gravity
     do idime =1, ndime
        if (abs(gravi_tur(idime)).gt.0.5_rp) vcoor = idime    
     end do
  end if


  ! calculate z_top of domain (only once)
  if(ipass==0) then ! to be done only once ! calculates maximum height
     z_top = -1.0e25_rp
     z_min = 1e20_rp
     do ipoin=1, npoin
        z_top = max(coord(vcoor,ipoin),z_top)
        z_min = min(coord(vcoor,ipoin),z_min)
     end do
     call PAR_MAX(z_top) 
     call PAR_MIN(z_min)
     !
     ! Allocate structures 
     !  
     nullify(z_wall)
     nullify(numer)
     nullify(denom)
     call memory_alloca(mem_modul(1:2,modul),'Z_WALL ','tur_endste',z_wall,npoin)
     call memory_alloca(mem_modul(1:2,modul),'NUMER','tur_endste',numer ,npoin)
     call memory_alloca(mem_modul(1:2,modul),'DENOM','tur_endste',denom ,npoin)
  end if ! ipass ==0 ! first step

  !
  ! Allocate structures
  !
  nullify(kfl_fixno_tmp)
  nullify(bvess_tmp)
  nullify(conve_tmp)
  nullify(force)
  call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_TMP','tur_endste',kfl_fixno_tmp,1_ip ,npoin)
  call memory_alloca(mem_modul(1:2,modul),'BVESS_TMP','tur_endste',bvess_tmp    ,1_ip ,npoin)
  call memory_alloca(mem_modul(1:2,modul),'CONVE_TMP','tur_endste',conve_tmp ,ndime,npoin)
  call memory_alloca(mem_modul(1:2,modul),'FORCE','tur_endste',force ,npoin)
  !
  ! load solver  structures
  !
  solve_sol                => solve(3:)
  solve_sol(1) % kfl_fixno => kfl_fixno_tmp 
  solve_sol(1) % bvess     => bvess_tmp
  solve_sol(1) % kfl_iffix =  1_ip  

  nloop= 1 ! use the solver only once
  if(ipass==0) then ! loads vertical wall distance z_wall solving  e_z\nabla z_wall= 1
     nloop = 3
     do ipoin =1,npoin
        conve_tmp(1:3,ipoin) = 0.0_rp ! convection in vertical direction
        conve_tmp(vcoor,ipoin) = 1.0_rp ! from bottom to top
        force(ipoin) = 1.0_rp
        kfl_fixno_tmp(1,ipoin) =  0_ip
        unkno(ipoin) = coord(3, ipoin)
        if(kfl_fixno_tur(1,ipoin,2) == 3) then  ! b.c. at WALL
           bvess_tmp(1,ipoin) = walld(ipoin)  !coord(3, ipoin)
           kfl_fixno_tmp(1,ipoin) =  1_ip
           unkno(ipoin) = bvess_tmp(1,ipoin)
        end if
     end do
     do iloop =1, nloop ! solve twice, and any more,to obtain better convergence
        call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)  
        if (INOTMASTER) &
             call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid)
        call solver_solve(solve_sol,amatr,rhsid,unkno)
     end do

     !  end if
     do ipoin = 1,npoin
        z_wall(ipoin) = unkno(ipoin)
     end do
  end if
  !
  ! integrate from wall to top
  !
  do ipoin = 1,npoin
     conve_tmp(1:3,ipoin) = 0.0_rp ! vertical convection  in vertical direction
     conve_tmp(vcoor,ipoin) = 1.0_rp ! from bottom to top
     kfl_fixno_tmp(1,ipoin) = 0      ! boundary condition free
     unkno(ipoin) = numer(ipoin) ! 0.0_rp
     !  numer(ipoin) = 0.0_rp  ! Unknown sqrt(k) z dz 
     !  denom(ipoin) = 0.0_rp  ! Unknown sqrt(k) dz
     if(kfl_fixno_tur(1,ipoin,2) == 3_ip) then ! if wall law for eps, then fix b.c.
        kfl_fixno_tmp(1,ipoin) =  1_ip
        bvess_tmp(1,ipoin)     =  0.0_rp  ! boundary value
     end if
  end do
  if (kfl_regularization) then        
     do ipoin =1, npoin
        force(ipoin) = sqrt(regul_k(untur(1,ipoin,1)))*z_wall(ipoin)
     end do
  else
     do ipoin = 1,npoin
        force(ipoin) = sqrt(untur(1,ipoin,1))*z_wall(ipoin) !  coord(3,1:npoin)
     end do
  end if

  ! assemble matrix and rhs for numer
  !  do iloop =1, nloop
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)
  if (INOTMASTER) &
       call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid) 
  call solver_solve(solve_sol,amatr,rhsid,unkno)
  !  end do
  do ipoin = 1,npoin
     numer (ipoin) =  unkno(ipoin)
  end do

  ! assemble matrix and rhs  for denom

  if (kfl_regularization) then
     do ipoin=1, npoin
        force(ipoin) = sqrt(regul_k(untur(1,ipoin,1)))
     end do
  else
     do ipoin = 1,npoin
        force(ipoin) = sqrt(untur(1,ipoin,1))
     end do
  end if
  do ipoin = 1,npoin
     unkno(ipoin) = denom(ipoin) !0.0_rp
  end do
  !  do iloop =1, nloop
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)
  if (INOTMASTER) &
       call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid)
  call solver_solve(solve_sol,amatr,rhsid,unkno)
  !  end do
  do ipoin = 1,npoin
     denom (ipoin)= unkno(ipoin)
  end do

  ! topography height ! integral from botton with rhs 0, only done the first time

  if (ipass==0) then
     do ipoin =1, npoin
        if (kfl_fixno_tur(1,ipoin,2) /= 3_ip) then ! not wall, denom .gt.0
           tur_max_mixlen(ipoin) = 0.075_rp*(numer(ipoin)/denom(ipoin) )
        else
           tur_max_mixlen(ipoin) = 0.0_rp
        end if
     end do
  end if

  do ipoin = 1,npoin
     unkno(ipoin) = tur_max_mixlen(ipoin)
  end do
  !     if (kfl_fixno_tur(1,ipoin,2) /= 3) & ! initial condition, denom .gt. 0
  !          tur_max_mixlen(ipoin) = 0.075_rp*(numer(ipoin)/denom(ipoin) ) !-z_wall(ipoin))     



  !
  ! project solution from the top vertically to all domain
  !
  do ipoin = 1,npoin
     conve_tmp(1:3,ipoin) = 0.0_rp
     conve_tmp(vcoor,ipoin) = - 1.0_rp  ! vertical convection from top to bottom
  end do
  do ipoin =1, npoin
     kfl_fixno_tmp(1,ipoin) =  0
     force(ipoin) = 0.0_rp
     if (abs(coord(vcoor, ipoin)-z_top)<0.1_rp ) then ! prescribe top of domain
        kfl_fixno_tmp(1,ipoin) =  1_ip
        unkno(ipoin)           =   0.075_rp*numer(ipoin)/denom(ipoin)
        bvess_tmp(1,ipoin)     =  unkno(ipoin) !tur_max_mixlen(ipoin)
     end if
  end do
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)  
  if (INOTMASTER) &
       call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid)
  call solver_solve(solve_sol,amatr,rhsid,unkno)

  do ipoin = 1,npoin
     tur_max_mixlen(ipoin) =  unkno(ipoin)
  end do

  !  if (INOTMASTER) tur_max_mixlen(1:npoin) =  z_wall(1:npoin)
  !
  ! Deallocate structures 
  !  
  nullify(solve_sol(1) % kfl_fixno) 
  nullify(solve_sol(1) % bvess    )
  if (INOTMASTER) then
     deallocate (kfl_fixno_tmp)
     deallocate (bvess_tmp)
     !  deallocate (z_wall)
     deallocate (conve_tmp)
     !  deallocate (numer)
     !  deallocate (denom)
     deallocate (force)
  end if


  ipass = 1
end subroutine tur_maxlen

