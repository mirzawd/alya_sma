!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_residual_monolithic.f90
!> @author  Guillaume Houzeaux
!> @brief   Navier-Stokes reisual
!> @details Compute the residual of the momentum and continuity 
!>          equations for a monolithic scheme
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_residual_monolithic(A,b,u,p,resid_mom,resid_con)
  use def_kintyp, only         : ip,rp
  use def_master, only         : INOTMASTER,zeror
  use def_domain, only         : ndime,npoin,r_dom,c_dom
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  real(rp),    intent(in)    :: A(ndime+1,ndime+1,*) 
  real(rp),    intent(in)    :: b(ndime+1,*)
  real(rp),    intent(in)    :: u(ndime,*)
  real(rp),    intent(in)    :: p(*)
  real(rp),    intent(out)   :: resid_mom
  real(rp),    intent(out)   :: resid_con
  integer(ip)                :: ndim1,idime,jdime,jpoin,ipoin,izdom
  real(rp),    allocatable   :: rm(:,:),bm(:,:)
  real(rp),    allocatable   :: rc(:),bc(:)
  real(rp)                   :: denom_mom,denom_con
  integer(ip)                :: npoin1

  npoin1 = max(1_ip,npoin)
  allocate( rm(ndime,npoin1), bm(ndime,npoin1) )
  allocate( rc(npoin1), bc(npoin1) )

  if( INOTMASTER ) then

     ndim1 = ndime + 1
     !$OMP PARALLEL DO SCHEDULE (STATIC)                           &
     !$OMP DEFAULT  ( NONE )                                       &
     !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime, jdime )          &           
     !$OMP SHARED   ( r_dom, c_dom, A, u, p, rm, rc, ndim1,        & 
#ifndef NDIMEPAR
     !$OMP            ndime,                                       &
#endif
     !$OMP            npoin ) 
     do ipoin = 1,npoin
        rm(1:ndime,ipoin) = 0.0_rp
        rc(ipoin)         = 0.0_rp
        do izdom = r_dom(ipoin),r_dom(ipoin+1)-1       
           jpoin = c_dom(izdom)
           do idime = 1,ndime
              do jdime = 1,ndime
                 rm(idime,ipoin) = rm(idime,ipoin) - A(jdime,idime,izdom) * u(jdime,jpoin)      
              end do
              rc(ipoin)       = rc(ipoin)       - A(idime,ndim1,izdom) * u(idime,jpoin)
              rm(idime,ipoin) = rm(idime,ipoin) - A(ndim1,idime,izdom) * p(jpoin)      
           end do
           rc(ipoin) = rc(ipoin) - A(ndim1,ndim1,izdom) * p(jpoin)
        end do
     end do 
     !$OMP END PARALLEL DO 

     call PAR_INTERFACE_NODE_EXCHANGE(ndime,rm,'SUM','IN MY CODE')
     call PAR_INTERFACE_NODE_EXCHANGE(1_ip, rc,'SUM','IN MY CODE')

     do ipoin = 1,npoin
        bm(1:ndime,ipoin) = b(1:ndime,ipoin)
        bc(ipoin)         = b(ndim1,ipoin)
        rm(1:ndime,ipoin) = b(1:ndime,ipoin) + rm(1:ndime,ipoin)
        rc(ipoin)         = b(ndim1,ipoin)   + rc(ipoin)
     end do
     
  end if
  
  call norm2x(ndime,bm,denom_mom)
  call norm2x( 1_ip,bc,denom_con)
  call norm2x(ndime,rm,resid_mom)
  call norm2x( 1_ip,rc,resid_con)

  resid_mom = resid_mom / ( denom_mom + zeror )
  resid_con = resid_con / ( denom_con + zeror )

  deallocate( rm, bm, rc, bc )

end subroutine nsi_residual_monolithic
