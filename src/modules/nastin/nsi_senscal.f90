!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> nsi_senscal.f90
!> @file nsi_senscal.f90 
!> @fn nsi_senscal 
!> This subroutine calculates the sensitivities using discrete adjoint
!>

subroutine nsi_senscal

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod, only : kfl_ndvars_opt,sens,sens_mesh,kfl_dvar_type,kfl_calc_sens
  use def_domain
  use def_nastin
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE,PAR_SUM
!  use mod_parall,         only : PAR_MY_CODE_RANK
  implicit none
  integer(ip) :: idesvar,ipoin,idime,ind,jpoin,jdime
  real(rp),    pointer     :: sens_dim(:)
  
  if( kfl_stead_nsi == 1 ) then
  
    !  
    ! start sensitivities calculations
    !
    if (kfl_dvar_type == 2) then 
    
      sens = 0.0_rp
      
      do idesvar = 1,kfl_ndvars_opt
        do ipoin=1,npoin
          do idime = 1,ndime
            ind = (ipoin-1)*ndime + idime
            sens(idesvar) = sens(idesvar) + resdiff_nsi(idesvar,ind)*veloc(idime,ipoin,1)
          end do
          sens(idesvar) = sens(idesvar) + resdiff_nsi(idesvar,npoin*ndime + ipoin)*press(ipoin,1)
        end do
        sens(idesvar) = sens(idesvar) + dcost_dx_nsi(idesvar)
        call PAR_SUM(sens(idesvar))
      end do
      print*,'sens',sens
   
    elseif (kfl_dvar_type == 5)  then
      !
      ! Geometric design variables: calculate dR/dD and sens_mesh
      !
  
!       sens_mesh = 0.0_rp
      kfl_calc_sens = 1_ip
      !
      ! calculate dR/dD
      !
      if( INOTMASTER) then
        call nsi_elmope_omp(30_ip)   
      end if
            
      nullify(sens_dim)
      allocate( sens_dim(npoin) )
      
      if( INOTMASTER) then
        do jdime = 1,ndime
          do jpoin = 1,npoin  
            idesvar = (jpoin-1)*ndime + jdime
            sens_dim(jpoin) = 0.0_rp
            sens_dim(jpoin) = sens_mesh(jdime,jpoin) + dcost_dx_nsi(idesvar)
          end do
          call PAR_INTERFACE_NODE_EXCHANGE(sens_dim,'SUM','IN MY CODE')
          sens_mesh(jdime,:) = sens_dim(:)
        end do !jdime
      end if

      deallocate( sens_dim )
      
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 63094 ) print*, "NACA3d63094 AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 165511) print*, "NACA3d165511AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 113354) print*, "NACA3d113354AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 84477 ) print*, "NACA3d84477 AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 104967) print*, "NACA3d104967AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 146834) print*, "NACA3d146834AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 91896 ) print*, "NACA3d91896 AAA", sens_mesh(:,ipoin )
!         enddo
!       endif

!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 178) print*, "NSI-NASA30P178", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 466) print*, "NSI-NASA30P466", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 475) print*, "NSI-NASA30P475", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 486) print*, "NSI-NASA30P486", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 502) print*, "NSI-NASA30P502", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 514) print*, "NSI-NASA30P514", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 782) print*, "NSI-NASA30P782", sens_mesh(:,ipoin )
!         enddo
!       endif
      
      
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 316 ) print*, "ONERANNN316 ", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 1706) print*, "ONERANNN1706", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 2315) print*, "ONERANNN2315", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 3148) print*, "ONERANNN3148", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4541) print*, "ONERANNN4541", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4769) print*, "ONERANNN4769", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5037) print*, "ONERANNN5037", sens_mesh(:,ipoin )
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
!           if (lninv_loc(ipoin) == 4160) print*, "NACA4160AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 4875) print*, "NACA4875AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5525) print*, "NACA5525AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5590) print*, "NACA5590AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 9035) print*, "NACA9035AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5915) print*, "NACA5915AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5720) print*, "NACA5720AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5460) print*, "NACA5460AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5330) print*, "NACA5330AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 5200) print*, "NACA5200AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 8840) print*, "NACA8840AAA", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 10465) print*, "NACA10465AAA", sens_mesh(:,ipoin )
!         enddo
!       endif
                  
      
    endif  !kfl_dvar_type
    
  endif !kfl_stead_nsi
    
  
!       do jdime = 1,ndime
!         do jpoin = 1,npoin  
!           idesvar = (jpoin-1)*ndime + jdime
!           sens_dim(jpoin) = 0.0_rp
!           do ipoin=1,npoin
!             do idime = 1,ndime
! 	      ind = (ipoin-1)*ndime + idime
! 	      sens_dim(jpoin) = sens_dim(jpoin) + resdiff_nsi(idesvar,ind)*veloc(idime,ipoin,1)
! 	    enddo
! 	    sens_dim(jpoin) = sens_dim(jpoin) + resdiff_nsi(idesvar,npoin*ndime + ipoin)*press(ipoin,1)
! 	  end do
! 	  sens_dim(jpoin) = sens_dim(jpoin) + dcost_dx_nsi(idesvar)
!         end do
!         call PAR_INTERFACE_NODE_EXCHANGE(sens_dim,'SUM','IN MY CODE')
!         sens_mesh(jdime,:) = sens_dim(:)      
!       end do !jdime
!
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 4379) print*, "NLR2D4379", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 7470) print*, "NLR2D7470", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 7150) print*, "NLR2D7150", sens_mesh(:,ipoin )
!         enddo
!       endif

!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 9670) print*, "cyl3d9670", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 9720) print*, "cyl3d9720", sens_mesh(:,ipoin )
!         enddo
!       endif      

!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 22698) print*, "NLR3d22698", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 23102) print*, "NLR3d23102", sens_mesh(:,ipoin )
!         enddo
!       endif      
! 
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 100203) print*, "NLR3d100203", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 90641) print*, "NLR3d90641", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 177201) print*, "NLR3d177201", sens_mesh(:,ipoin )
!         enddo
!       endif      
      
      
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 79254) print*, "onera1000 79254", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 57058) print*, "onera1000 57058", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 121043) print*, "onera1000 121043", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 122370) print*, "onera1000 122370", sens_mesh(:,ipoin )
!         enddo
!       endif  
!       
!       ! onera
!       if (INOTMASTER) then
!         do ipoin=1,npoin
!           if (lninv_loc(ipoin) == 37488) print*, "onera37488", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 53510) print*, "onera53510", sens_mesh(:,ipoin )
!           if (lninv_loc(ipoin) == 96527) print*, "onera96527", sens_mesh(:,ipoin )
!         enddo
!       endif      


!     ipoin = 1
!     print*, "ssssssssss", sens( 1 : 10 )

!     if (INOTMASTER) then
!       do ipoin=1,npoin
!         if (lninv_loc(ipoin) == 2) print*, "sens 2", sens_mesh(:,ipoin )
!         if (lninv_loc(ipoin) == 3) print*, "sens 3", sens_mesh(:,ipoin )
!         if (lninv_loc(ipoin) == 50) print*, "sens50", sens_mesh(:,ipoin )
!       enddo
!     endif
        
!     if (INOTMASTER) then
!       do ipoin=1,npoin
!         if (lninv_loc(ipoin) == 6480) print*, "senss6480", sens_mesh(:,ipoin )
!         if (lninv_loc(ipoin) == 6500) print*, "senss6500", sens_mesh(:,ipoin )
!         if (lninv_loc(ipoin) == 7300) print*, "senss7300", sens_mesh(:,ipoin )
!       enddo
!     endif


  
  
!     do ipoin=npoi1,npoin
!       numeracion global , PAR_MY_CODE_RANK, 
!     enddo
!     
!     if (INOTMASTER) then
!     do ipoin=npoi1,npoin
!       print*, PAR_MY_CODE_RANK,  lninv_loc(ipoin)
!     enddo
!     endif
  
  
  
  
  
end subroutine nsi_senscal
