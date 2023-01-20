!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_dipeps.f90
!> @author  Mariano Vazquez
!> @date    17/04/2017
!> @brief   Compute displacement perturbation considering BC constraints
!> @details Compute displacement perturbation considering BC constraints
!> @}
!-----------------------------------------------------------------------
subroutine sld_dipeps()
  use def_parame
  use def_master
  use def_domain
  use def_solidz

  implicit none
  integer(ip) :: ipoin,idime,ibopo
  real(rp)    :: vmodu,vdirec(ndime),vepsex(ndime),dummy_matrix(ndime,ndime)

  if( INOTMASTER ) then

     do ipoin = 1,npoin
        vepsex = epsex_sld
        ibopo = lpoty(ipoin)
        if (ibopo > 0) then
           if( kfl_fixno_sld(1,ipoin) == 2 .or. kfl_fixno_sld(1,ipoin) == 3 ) then
              !
              ! Rotate vepsex
              !
              call sld_rotsys(3_ip, &
                   1_ip,1_ip,ndime,ndime,&
                   dummy_matrix,vepsex,&
                   jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))

              !
              ! Correct on local basis
              !
              vepsex(1) = 0.0_rp
              if (kfl_fixno_sld(    2,ipoin) == 1) vepsex(    2) = 0.0_rp
              if (kfl_fixno_sld(ndime,ipoin) == 1) vepsex(ndime) = 0.0_rp

              !
              ! Rotate back the corrected vepsex
              !
              call sld_rotsys(2_ip, &
                   1_ip,1_ip,ndime,ndime,&
                   dummy_matrix,vepsex,&
                   jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))

           else
              !
              ! Check and correct on cartesian basis
              !
              if (kfl_fixno_sld(    1,ipoin) == 1) then
                 vepsex(    1) = 0.0_rp
              end if
              if (kfl_fixno_sld(    2,ipoin) == 1) then
                 vepsex(    2) = 0.0_rp
              end if
              if (kfl_fixno_sld(ndime,ipoin) == 1) then
                 vepsex(    ndime) = 0.0_rp
              end if
           end if
        end if

        !
        ! Normalize last displacements
        !
        vdirec(1:ndime)= displ(1:ndime,ipoin,ITER_K)
        vmodu= vdirec(1)*vdirec(1) + vdirec(2)*vdirec(2)
        if (ndime == 3_ip) vmodu= vmodu + vdirec(ndime)*vdirec(ndime)
        vmodu= sqrt(vmodu)
        if (vmodu > 1.0e-12_rp) then
           !
           ! Compute displacement perturbation for node ipoin in the direction of the
           ! last displacement
           !
           vdirec(1:ndime) = vdirec(1:ndime) / vmodu
           do idime= 1,ndime
              disep_sld(idime,ipoin)= vdirec(idime) * vepsex(idime)
           end do
        else
           !
           ! No last displacement, use directly vepsex
           !
           do idime= 1,ndime
              disep_sld(idime,ipoin)= real(idime,rp)*vepsex(idime)
           end do
        end if
     end do

  end if

end subroutine sld_dipeps
