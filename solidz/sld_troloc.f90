!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_troloc.f90
!> @author  Mariano Vazquez
!> @brief   Compute local rotation tensor (only for membrane-like aneurysms with prisms)
!> @details Compute local rotation tensor (only for membrane-like aneurysms with prisms)
!> @} 
!-----------------------------------------------------------------------
subroutine sld_troloc(itrol,G_ij)
  use def_domain
  use def_master
  use def_solidz

  implicit none
  integer(ip)  :: ielem, idime, jdime, inode, ipoin,itrol
  real(rp)     :: vecau(ndime), tenau(ndime,ndime), vemod, G_ij(ndime,ndime), vdumm(ndime,ndime)

  ! it is not required to perform an rhsmod afterwards because there is no assembly


  if ( itrol == 0) then
     
     do ielem=1,nelem
        tenau(1,ndime) = coord(1,lnods(2,ielem)) - coord(1,lnods(5,ielem))
        tenau(2,ndime) = coord(2,lnods(2,ielem)) - coord(2,lnods(5,ielem))
        tenau(3,ndime) = coord(3,lnods(2,ielem)) - coord(3,lnods(5,ielem))
        vemod= sqrt(tenau(1,ndime)*tenau(1,ndime) + tenau(2,ndime)*tenau(2,ndime) + tenau(3,ndime)*tenau(3,ndime))
        tenau(1,ndime) = tenau(1,ndime)/vemod 
        tenau(2,ndime) = tenau(2,ndime)/vemod 
        tenau(3,ndime) = tenau(3,ndime)/vemod

        vecau(1) = coord(1,lnods(2,ielem)) - coord(1,lnods(1,ielem))
        vecau(2) = coord(2,lnods(2,ielem)) - coord(2,lnods(1,ielem))
        vecau(3) = coord(3,lnods(2,ielem)) - coord(3,lnods(1,ielem))
        vemod= sqrt(vecau(1)*vecau(1) + vecau(2)*vecau(2) + vecau(3)*vecau(3))
        vecau(1) = vecau(1)/vemod 
        vecau(2) = vecau(2)/vemod 
        vecau(3) = vecau(3)/vemod

        tenau(1,1) = vecau(2) * tenau(3,ndime) - vecau(3) * tenau(2,ndime)
        tenau(2,1) = vecau(3) * tenau(1,ndime) - vecau(1) * tenau(3,ndime)
        tenau(3,1) = vecau(1) * tenau(2,ndime) - vecau(2) * tenau(1,ndime)
        
        tenau(1,2) = tenau(2,1) * tenau(3,ndime) - tenau(3,1) * tenau(2,ndime)
        tenau(2,2) = tenau(3,1) * tenau(1,ndime) - tenau(1,1) * tenau(3,ndime)
        tenau(3,2) = tenau(1,1) * tenau(2,ndime) - tenau(2,1) * tenau(1,ndime)

        do inode= 1,6
           ipoin= lnods(inode,ielem)
           do idime= 1,ndime
              do jdime= 1,ndime
                 roloc_sld(idime,jdime,ipoin)= tenau(idime,jdime)
              end do
           end do
        end do        

     end do
  
  else
     
     ipoin= itrol
     
     !rotate G_ij
     vdumm(1,1) = &
          G_ij(1,1) * roloc_sld(1,1,ipoin) + &
          G_ij(1,2) * roloc_sld(2,1,ipoin) + &
          G_ij(1,3) * roloc_sld(3,1,ipoin) 
     vdumm(1,2) = &
          G_ij(1,1) * roloc_sld(1,2,ipoin) + &
          G_ij(1,2) * roloc_sld(2,2,ipoin) + &
          G_ij(1,3) * roloc_sld(3,2,ipoin) 
     vdumm(1,3) = &
          G_ij(1,1) * roloc_sld(1,3,ipoin) + &
          G_ij(1,2) * roloc_sld(2,3,ipoin) + &
          G_ij(1,3) * roloc_sld(3,3,ipoin) 
     vdumm(2,1) = &
          G_ij(2,1) * roloc_sld(1,1,ipoin) + &
          G_ij(2,2) * roloc_sld(2,1,ipoin) + &
          G_ij(2,3) * roloc_sld(3,1,ipoin) 
     vdumm(2,2) = &
          G_ij(2,1) * roloc_sld(1,2,ipoin) + &
          G_ij(2,2) * roloc_sld(2,2,ipoin) + &
          G_ij(2,3) * roloc_sld(3,2,ipoin) 
     vdumm(2,3) = &
          G_ij(2,1) * roloc_sld(1,3,ipoin) + &
          G_ij(2,2) * roloc_sld(2,3,ipoin) + &
          G_ij(2,3) * roloc_sld(3,3,ipoin) 
     vdumm(3,1) = &
          G_ij(3,1) * roloc_sld(1,1,ipoin) + &
          G_ij(3,2) * roloc_sld(2,1,ipoin) + &
          G_ij(3,3) * roloc_sld(3,1,ipoin) 
     vdumm(3,2) = &
          G_ij(3,1) * roloc_sld(1,2,ipoin) + &
          G_ij(3,2) * roloc_sld(2,2,ipoin) + &
          G_ij(3,3) * roloc_sld(3,2,ipoin) 
     vdumm(3,3) = &
          G_ij(3,1) * roloc_sld(1,3,ipoin) + &
          G_ij(3,2) * roloc_sld(2,3,ipoin) + &
          G_ij(3,3) * roloc_sld(3,3,ipoin) 
     
     G_ij(1,1) = &
          roloc_sld(1,1,ipoin) * vdumm(1,1) + &
          roloc_sld(2,1,ipoin) * vdumm(2,1) + &
          roloc_sld(3,1,ipoin) * vdumm(3,1) 
     G_ij(1,2) = &
          roloc_sld(1,1,ipoin) * vdumm(1,2) + &
          roloc_sld(2,1,ipoin) * vdumm(2,2) + &
          roloc_sld(3,1,ipoin) * vdumm(3,2) 
     G_ij(1,3) = 0.0_rp
     
     G_ij(2,1) = &
          roloc_sld(1,2,ipoin) * vdumm(1,1) + &
          roloc_sld(2,2,ipoin) * vdumm(2,1) + &
          roloc_sld(3,2,ipoin) * vdumm(3,1) 
     G_ij(2,2) = &
          roloc_sld(1,2,ipoin) * vdumm(1,2) + &
          roloc_sld(2,2,ipoin) * vdumm(2,2) + &
          roloc_sld(3,2,ipoin) * vdumm(3,2) 
     G_ij(2,3) = 0.0_rp
     
     G_ij(3,1) = 0.0_rp
     G_ij(3,2) = 0.0_rp
     G_ij(3,3) = 0.0_rp
     
  end if

end subroutine sld_troloc
