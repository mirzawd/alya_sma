!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_vect_maths.f90
!> @author  Adria Quintanas-Corominas
!> @date    November 2021
!> @brief   Toolbox with some maths operations for the vectorised assembly
!> @details Toolbox with some maths operations for the vectorised assembly
!------------------------------------------------------------------

module mod_sld_vect_maths

#include "def_solidz_vect_dim.inc" 
   use def_kintyp, only  : ip, rp, lg
   use def_domain, only  : ndime

   implicit none

   contains

      pure function vmath_TRACE(A) result(tr)

         real(rp),    intent(in)  :: A(:,:,:)
         real(rp)                 :: tr(DVS)
         integer(ip)              :: ii

         tr(1:DVS) = 0.0_rp
         do ii = 1, size(A,dim=2)
            tr(1:DVS) = tr(1:DVS) + A(1:DVS,ii,ii)
         enddo

      end function vmath_TRACE


      pure function vmath_MxV(M,V) result(R)

         real(rp),    intent(in)  :: M(:,:,:)
         real(rp),    intent(in)  :: V(:,:)
         real(rp)                 :: R(DVS,size(M,dim=2))
         integer(ip)              :: ii, jj

         R(1:DVS,:) = 0.0_rp
         do jj = 1, size(M,dim=3)
            do ii = 1, size(M,dim=2)
               R(1:DVS,ii) = R(1:DVS,ii) + M(1:DVS,ii,jj) * V(1:DVS,jj)
            enddo
         enddo

      end function vmath_MxV

      pure function MxV_voigt(M,V) result(R)

         real(rp),    intent(in)  :: M(DVS,ZVOIG,ZVOIG)
         real(rp),    intent(in)  :: V(DVS,ZVOIG)
         real(rp)                 :: R(DVS,ZVOIG)
         integer(ip)              :: ii, jj

         R(1:DVS,:) = 0.0_rp
         do jj = 1, ZVOIG
            do ii = 1, ZVOIG
               R(1:DVS,ii) = R(1:DVS,ii) + M(1:DVS,ii,jj) * V(1:DVS,jj)
            enddo
         enddo

      end function MxV_voigt


      pure function vmath_TxV(T,V) result(R)

         real(rp),    intent(in)  :: T(:,:,:)
         real(rp),    intent(in)  :: V(:,:)
         real(rp)                 :: R(DVS,size(T,dim=2))
         integer(ip)              :: ii, jj

         R(1:DVS,:) = 0.0_rp
         do jj = 1, size(T,dim=2)
            do ii = 1, size(T,dim=3)
               R(1:DVS,ii) = R(1:DVS,ii) + T(1:DVS,jj,ii) * V(1:DVS,jj)
            enddo
         enddo

      end function vmath_TxV


      pure function TxV_voigt(T,V) result(R)

         real(rp),    intent(in)  :: T(DVS,ZVOIG,ZVOIG)
         real(rp),    intent(in)  :: V(DVS,ZVOIG)
         real(rp)                 :: R(DVS,ZVOIG)
         integer(ip)              :: ii, jj

         R(1:DVS,:) = 0.0_rp
         do jj = 1, ZVOIG
            do ii = 1, ZVOIG
               R(1:DVS,ii) = R(1:DVS,ii) + T(1:DVS,jj,ii) * V(1:DVS,jj)
            enddo
         enddo

      end function TxV_voigt


      pure function vmath_MxT_v2(A,ai,aj,B,bi,bj) result(C)

         integer(ip), intent(in)  :: ai, aj
         integer(ip), intent(in)  :: bi, bj
         real(rp),    intent(in)  :: A(DVS,ai,aj)
         real(rp),    intent(in)  :: B(DVS,bi,bj)
         real(rp)                 :: C(DVS,ai,bi)
         integer(ip)              :: ii, jj, kk

         do kk = 1, aj
            C(1:DVS,:,:) = 0.0_rp
            do jj = 1, bi
               do ii = 1, ai
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,jj,kk)
               enddo
            enddo
         enddo

      end function vmath_MXT_v2


      pure function vmath_VxMxV(a,M,b) result(r)

         real(rp),    intent(in)  :: a(:,:)
         real(rp),    intent(in)  :: M(:,:,:)
         real(rp),    intent(in)  :: b(:,:)
         real(rp)                 :: r(DVS) 
         integer(ip)              :: ii, jj

         r(1:DVS) = 0.0_rp
         do jj = 1, size(b,dim=2)
            do ii = 1, size(a,dim=2)
               r(1:DVS) = r(1:DVS) + a(1:DVS,ii) * M(1:DVS,ii,jj) * b(1:DVS,jj)
            enddo
         enddo

      end function vmath_VxMxV


      pure function vmath_AxB(A,B) result(C)

         real(rp),    intent(in)  :: A(:,:,:)
         real(rp),    intent(in)  :: B(:,:,:)
         real(rp)                 :: C(DVS,size(A,dim=2),size(B,dim=3))
         integer(ip)              :: ii, jj, kk

         do jj = 1, size(B,dim=3)
            do ii = 1, size(A,dim=2)
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, size(A,dim=3)
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,kk,jj)
               enddo
            enddo
         enddo

      end function vmath_AXB


      pure function AxB_dime(A,B) result(C)

         real(rp),    intent(in)  :: A(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)  :: B(DVS,ZDIME,ZDIME)
         real(rp)                 :: C(DVS,ZDIME,ZDIME)
         integer(ip)              :: ii, jj, kk

         do jj = 1, ZDIME
            do ii = 1, ZDIME
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, ZDIME
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,kk,jj)
               enddo
            enddo
         enddo

      end function AxB_dime

      
      pure function SxAxB_dime(S,A,B) result(C)

         real(rp),    intent(in)  :: S(DVS)
         real(rp),    intent(in)  :: A(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)  :: B(DVS,ZDIME,ZDIME)
         real(rp)                 :: C(DVS,ZDIME,ZDIME)
         integer(ip)              :: ii, jj, kk

         do jj = 1, ZDIME
            do ii = 1, ZDIME
              C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, ZDIME
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,kk,jj)
               enddo
               C(1:DVS,ii,jj) =  S(1:DVS)*C(1:DVS,ii,jj)
            enddo
         enddo

      end function SxAxB_dime

      
      pure function AxB_voigt(A,B) result(C)

         real(rp),    intent(in)  :: A(DVS,ZVOIG,ZVOIG)
         real(rp),    intent(in)  :: B(DVS,ZVOIG,ZVOIG)
         real(rp)                 :: C(DVS,ZVOIG,ZVOIG)
         integer(ip)              :: ii, jj, kk

         do jj = 1, ZVOIG
            do ii = 1, ZVOIG
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, ZVOIG
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,kk,jj)
               enddo
            enddo
         enddo

      end function AxB_voigt


      pure function vmath_MxT(A,B) result(C)

         real(rp),    intent(in)  :: A(:,:,:)
         real(rp),    intent(in)  :: B(:,:,:)
         real(rp)                 :: C(DVS,size(A,dim=2),size(B,dim=2))
         integer(ip)              :: ii, jj, kk

!         do jj = 1, size(B,dim=2)
!            do ii = 1, size(A,dim=2)
!               C(1:DVS,ii,jj) = 0.0_rp
!            enddo
!         enddo

         do jj = 1, size(B,dim=2)
            do ii = 1, size(A,dim=2)
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, size(A,dim=3)
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,jj,kk)
               enddo
            enddo
         enddo

      end function vmath_MXT


      pure function MxT_ndime(A,B) result(C)
         ! Only for NDIME
         real(rp),    intent(in)  :: A(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)  :: B(DVS,ZDIME,ZDIME)
         real(rp)                 :: C(DVS,ZDIME,ZDIME)
         integer(ip)              :: ii, jj, kk

         do jj = 1, ZDIME
            do ii = 1, ZDIME
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, ZDIME
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,ii,kk) * B(1:DVS,jj,kk)
               enddo
            enddo
         enddo

      end function MXT_ndime


      pure function vmath_TxM(A,B) result(C)

         real(rp),    intent(in)  :: A(:,:,:)
         real(rp),    intent(in)  :: B(:,:,:)
         real(rp)                 :: C(DVS,size(A,dim=3),size(B,dim=3))
         integer(ip)              :: ii, jj, kk

         do jj = 1, size(B,dim=3)
            do ii = 1, size(A,dim=3)
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, size(A,dim=2)
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,kk,ii) * B(1:DVS,kk,jj)
               enddo
            enddo
         enddo

      end function vmath_TXM


      pure function TxM_ndime(A,B) result(C)

         real(rp),    intent(in)  :: A(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)  :: B(DVS,ZDIME,ZDIME)
         real(rp)                 :: C(DVS,ZDIME,ZDIME)
         integer(ip)              :: ii, jj, kk

         do jj = 1, ZDIME 
            do ii = 1, ZDIME 
               C(1:DVS,ii,jj) = 0.0_rp
               do kk = 1, ZDIME 
                  C(1:DVS,ii,jj) = C(1:DVS,ii,jj) + A(1:DVS,kk,ii) * B(1:DVS,kk,jj)
               enddo
            enddo
         enddo

      end function TXM_ndime


      pure function vmath_OUT_VxV(u,v) result(R)

         real(rp),    intent(in)  :: u(:,:)
         real(rp),    intent(in)  :: v(:,:)
         real(rp)                 :: R(DVS,size(u,dim=2),size(v,dim=2))
         integer(ip)              :: ii, jj

         do jj = 1, size(v,dim=2)
            do ii = 1, size(u,dim=2)
               R(1:DVS,ii,jj) = u(1:DVS,ii) * v(1:DVS,jj)
            enddo
         enddo
      end function vmath_OUT_VxV


      pure function vmath_OUT_MxM(A,B) result(R)

         real(rp),    intent(in)  :: A(:,:,:), B(:,:,:)
         real(rp)                 :: R(DVS,size(A,dim=2),size(A,dim=3),size(B,dim=2),size(B,dim=3))
         integer(ip)              :: ii, jj, ll, kk

         do ll = 1, size(B,dim=3)
            do kk = 1, size(B,dim=2)
               do jj = 1, size(A,dim=3)
                  do ii = 1, size(A,dim=2)
                     R(1:DVS,ii,jj,kk,ll) = A(1:DVS,ii,jj) * B(1:DVS,kk,ll)
                  enddo
               enddo
            enddo
         enddo

      end function vmath_OUT_MxM


      pure function vmath_DET(ndime,A) result(det)

         integer(ip), intent(in)  :: ndime
         real(rp),    intent(in)  :: A(DVS,ndime,ndime)
         real(rp)                 :: det(DVS)
         real(rp)                 :: t1(DVS), t2(DVS), t3(DVS)

         if( ndime == 2_ip )then

            det(1:DVS) = A(1:DVS,1,1) * A(1:DVS,2,2) - A(1:DVS,2,1) * A(1:DVS,1,2)
            
         else

            t1(1:DVS)  = A(1:DVS,2,2) * A(1:DVS,3,3) - A(1:DVS,3,2) * A(1:DVS,2,3)
            t2(1:DVS)  =-A(1:DVS,2,1) * A(1:DVS,3,3) + A(1:DVS,3,1) * A(1:DVS,2,3)
            t3(1:DVS)  = A(1:DVS,2,1) * A(1:DVS,3,2) - A(1:DVS,3,1) * A(1:DVS,2,2)
            det(1:DVS) = A(1:DVS,1,1) * t1(1:DVS) + A(1:DVS,1,2) * t2(1:DVS) + A(1:DVS,1,3) * t3(1:DVS)

         endif

      end function vmath_DET


      pure function vmath_INV(ndime,A) result(B)

         integer(ip), intent(in)  :: ndime
         real(rp),    intent(in)  :: A(DVS,ndime,ndime)
         real(rp)                 :: B(DVS,ndime,ndime)
         integer(ip)              :: iv
         real(rp)                 :: det(DVS), denom(DVS), t1(DVS), t2(DVS), t3(DVS)

         if( ndime == 2_ip )then

            det(1:DVS) = A(1:DVS,1,1) * A(1:DVS,2,2) - A(1:DVS,2,1) * A(1:DVS,1,2)
            do iv = 1, DVS
               if( abs(det(iv)) > 0.0_rp )then
                  denom(iv) = 1.0_rp / det(iv)
               else
                  denom(iv) = 0.0_rp
               endif
            enddo
            B(1:DVS,1,1) =  A(1:DVS,2,2) * denom(1:DVS)
            B(1:DVS,2,2) =  A(1:DVS,1,1) * denom(1:DVS)
            B(1:DVS,2,1) = -A(1:DVS,2,1) * denom(1:DVS)
            B(1:DVS,1,2) = -A(1:DVS,1,2) * denom(1:DVS)

         else

            t1(1:DVS)  = A(1:DVS,2,2) * A(1:DVS,3,3) - A(1:DVS,3,2) * A(1:DVS,2,3)
            t2(1:DVS)  =-A(1:DVS,2,1) * A(1:DVS,3,3) + A(1:DVS,3,1) * A(1:DVS,2,3)
            t3(1:DVS)  = A(1:DVS,2,1) * A(1:DVS,3,2) - A(1:DVS,3,1) * A(1:DVS,2,2)
            det(1:DVS) = A(1:DVS,1,1) * t1(1:DVS) + A(1:DVS,1,2) * t2(1:DVS) + A(1:DVS,1,3) * t3(1:DVS)
            do iv = 1, DVS
               if( abs(det(iv)) > 0.0_rp )then
                  denom(iv) = 1.0_rp / det(iv)
               else
                  denom(iv) = 0.0_rp
               endif
            enddo
            B(1:DVS,1,1) = t1(1:DVS) * denom(1:DVS)
            B(1:DVS,2,1) = t2(1:DVS) * denom(1:DVS)
            B(1:DVS,3,1) = t3(1:DVS) * denom(1:DVS)
            B(1:DVS,1,2) = (-A(1:DVS,1,2) * A(1:DVS,3,3) + A(1:DVS,3,2) * A(1:DVS,1,3)) * denom(1:DVS)
            B(1:DVS,2,2) = ( A(1:DVS,1,1) * A(1:DVS,3,3) - A(1:DVS,3,1) * A(1:DVS,1,3)) * denom(1:DVS)
            B(1:DVS,3,2) = (-A(1:DVS,1,1) * A(1:DVS,3,2) + A(1:DVS,1,2) * A(1:DVS,3,1)) * denom(1:DVS)
            B(1:DVS,1,3) = ( A(1:DVS,1,2) * A(1:DVS,2,3) - A(1:DVS,2,2) * A(1:DVS,1,3)) * denom(1:DVS)
            B(1:DVS,2,3) = (-A(1:DVS,1,1) * A(1:DVS,2,3) + A(1:DVS,2,1) * A(1:DVS,1,3)) * denom(1:DVS)
            B(1:DVS,3,3) = ( A(1:DVS,1,1) * A(1:DVS,2,2) - A(1:DVS,2,1) * A(1:DVS,1,2)) * denom(1:DVS)

         endif

      end function vmath_INV

end module mod_sld_vect_maths
!> @}
