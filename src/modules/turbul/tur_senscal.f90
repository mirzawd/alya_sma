!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> tur_senscal.f90
!> @file tur_senscal.f90 
!> @fn tur_senscal 
!> This subroutine calculates the sensitivities using discrete adjoint
!>

subroutine tur_senscal

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod, only : kfl_ndvars_opt,sens,kfl_dvar_type,sens_mesh,kfl_calc_sens,kfl_nwall
  use def_kermod, only : kfl_fixno_walld_ker
  use def_domain
  use def_turbul
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE,PAR_SUM
  use mod_memory
  
  implicit none
  integer(ip) :: idesvar,ipoin,jdime,jpoin,idime
  real(rp)    :: resdiff_des(npoin)
  real(rp)    :: inner,sens_con
  real(rp),    pointer     :: sens_dim(:)
  
  
  if( kfl_calc_sens == 1 ) then
      
    if (kfl_dvar_type == 4) then
      
      ! initialization
      inner = 0.0_rp
      do ipoin=1,npoin
        resdiff_des(ipoin)=0.0_rp
      end do
      ! start sensitivities calculations
      do idesvar = 1,kfl_ndvars_opt
        do ipoin=1,npoin
          resdiff_des(ipoin)=resdiff_tur(idesvar,ipoin)
        end do
        do ipoin=1,npoin
          sens(idesvar) = sens(idesvar) + resdiff_des(ipoin)*untur(idesvar,ipoin,1)
          ! sens(idesvar) = sens(idesvar) + resdiff_des(ipoin)*resdiff_des(ipoin)
       enddo
       call PAR_SUM(sens(idesvar))
      end do
!       print*,'sens tur', sens
      
    !
    ! Parametric design variables
    !
    elseif (kfl_dvar_type == 6) then
    
      sens_con = 0.0_rp
      if( INOTMASTER) then
        sens_con = sens_mesh(1,1)
      endif
      call PAR_SUM(sens_con, 'IN MY CODE' )
      if( INOTMASTER) then
      else 
        print *, sens_con
      endif
            
    elseif (kfl_dvar_type == 5) then
      !
      ! Geometric design variables: calculate dR/dD and sens_mesh
      !
      nullify(sens_wall)
      if( INOTMASTER) then        
        call memory_alloca(mem_modul(1:2,modul),'SENS_WALL','tur_memall',sens_wall,ndime,kfl_nwall)
      else
        call memory_alloca(mem_modul(1:2,modul),'SENS_WALL','tur_memall',sens_wall,ndime,1_ip)
      end if
      sens_wall = 0.0_rp
      if( INOTMASTER) then
        iunkn_tur = 1
        call tur_elmop2(30_ip)   
      end if
      if( INOTMASTER) then
        call PAR_SUM(sens_wall, 'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
      endif
      if( INOTMASTER) then
        do ipoin = 1,npoi1
          if( kfl_fixno_walld_ker(1,ipoin) == 1 ) then
            do idime = 1,ndime
!               sens_mesh(idime,ipoin) = sens_mesh(idime,ipoin) + sens_wall(idime,walli(ipoin))
              sens_mesh(idime,ipoin) = sens_mesh(idime,ipoin) + sens_wall(idime,wallo(ipoin))
            enddo
          endif
        enddo
        do ipoin = npoi2,npoi3
          if( kfl_fixno_walld_ker(1,ipoin) == 1 ) then
            do idime = 1,ndime
!               sens_mesh(idime,ipoin) = sens_mesh(idime,ipoin) + sens_wall(idime,walli(ipoin))
              sens_mesh(idime,ipoin) = sens_mesh(idime,ipoin) + sens_wall(idime,wallo(ipoin))
            enddo
          endif
        enddo
      endif
            
      nullify(sens_dim)
      allocate( sens_dim(npoin) )
      
      if( INOTMASTER) then
      do jdime = 1,ndime
        do jpoin = 1,npoin  
          sens_dim(jpoin) = 0.0_rp
	      sens_dim(jpoin) = sens_mesh(jdime,jpoin) 
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(sens_dim,'SUM','IN MY CODE')
        sens_mesh(jdime,:) = sens_dim(:)     
      end do !jdime
      endif
      
      
      deallocate( sens_dim )
      call memory_deallo(mem_modul(1:2,modul),'SENS_WALL','tur_memall',sens_wall)

!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 63094 ) print*, "NACA3d63094 TTT", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 165511) print*, "NACA3d165511TTT", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 113354) print*, "NACA3d113354TTT", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 84477 ) print*, "NACA3d84477 TTT", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 104967) print*, "NACA3d104967TTT", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 146834) print*, "NACA3d146834TTT", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 91896 ) print*, "NACA3d91896 TTT", sens_mesh(:,ipoin )
!         enddo
!       endif
      
      
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 178) print*, "TUR-NASA30P178", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 466) print*, "TUR-NASA30P466", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 475) print*, "TUR-NASA30P475", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 486) print*, "TUR-NASA30P486", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 502) print*, "TUR-NASA30P502", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 514) print*, "TUR-NASA30P514", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 782) print*, "TUR-NASA30P782", sens_mesh(:,ipoin )
!         enddo
!       endif
      

!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 316 ) print*, "ONERA316 ", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 1706) print*, "ONERA1706", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 2315) print*, "ONERA2315", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 3148) print*, "ONERA3148", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4541) print*, "ONERA4541", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4769) print*, "ONERA4769", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5037) print*, "ONERA5037", sens_mesh(:,ipoin )
!         enddo
!       endif
      
            
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 3788) print*, "RAMP3788", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4260) print*, "RAMP4260", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4594) print*, "RAMP4594", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5304) print*, "RAMP5304", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 14500) print*, "RAMP14500", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 3643 ) print*, "RAMP3643 ", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 3367 ) print*, "RAMP3367 ", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4098 ) print*, "RAMP4098 ", sens_mesh(:,ipoin )
!         enddo
!       endif
      
      
      
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 4160) print*, "NACA4160ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4875) print*, "NACA4875ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5525) print*, "NACA5525ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5590) print*, "NACA5590ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 9035) print*, "NACA9035ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5915) print*, "NACA5915ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4720) print*, "NACA5720ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5460) print*, "NACA5460ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5330) print*, "NACA5330ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5200) print*, "NACA5200ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 8840) print*, "NACA8840ttt", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 10465) print*, "NACA10465ttt", sens_mesh(:,ipoin )
!         enddo
!       endif
            
    endif !kfl_dvar_type
    
  endif !kfl_calc_sens
  
  
    
  
end subroutine tur_senscal
