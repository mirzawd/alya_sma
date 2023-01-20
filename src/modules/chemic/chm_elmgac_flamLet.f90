!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_elmgac_flamLet(&
     pnode,lnods,elcod,elcon,elvel,elmas)
  !------------------------------------------------------------------------
  !****f* Chemic/chm_elmgac_cfi
  ! NAME
  !    chm_elmgac_cfi
  ! DESCRIPTION
  !    Gather operations for the combustion models
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,coord
  use def_master, only     :  conce,&
                              massk,advec
  use def_chemic, only     :  kfl_advec_chm,&
       &                      nclas_chm,&
       &                      ADR_chm,&
       &                      kfl_lookg_chm
  use mod_ADR,    only     :  BDF

  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  real(rp),    intent(out) :: elcon(pnode,nclas_chm,*)
  real(rp),    intent(out) :: elvel(ndime,pnode)
  real(rp),    intent(out) :: elmas(pnode,nclas_chm)                ! Mass source terms
  integer(ip)              :: inode,ipoin,idime,itime,iclas

  !
  ! Concentration and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     do iclas=1,nclas_chm
        elcon(inode,iclas,1) = conce(ipoin,iclas,1)
     end do
     do idime=1,ndime
        elcod(idime,inode)   = coord(idime,ipoin)
     end do
  end do

  !
  ! Time integration
  !
  if( ADR_chm(1) % kfl_time_integration == 1 ) then
     do iclas = 1,nclas_chm
        do inode = 1,pnode
           ipoin = lnods(inode)
           elcon(inode,iclas,2) = conce(ipoin,iclas,3)
        end do
     end do

     if( ADR_chm(1) % kfl_time_scheme == BDF ) then
        do iclas = 1,nclas_chm
           do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elcon(inode,iclas,itime) = conce(ipoin,iclas,itime+1)
              end do
           end do
        end do
     end if

  end if

  !
  ! Advection
  !
  if( kfl_advec_chm /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elvel(idime,inode) = advec(idime,ipoin,1)
        end do
     end do
  else
     elvel = 0.0_rp
  end if

  !
  ! Mass source terms coefficients
  !
  if (kfl_lookg_chm == 0) then
     do iclas = 1, nclas_chm
        do inode = 1,pnode
           ipoin = lnods(inode)
           elmas(inode,iclas) = massk(ipoin,iclas)
        enddo
     end do
  else
     elmas = 0.0_rp
  endif

end subroutine chm_elmgac_flamLet
