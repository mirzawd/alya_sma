!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_exampl()
  !-----------------------------------------------------------------------
  !****f* Parall/par_exampl
  ! NAME
  !    par_exampl
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_domain
  use def_master
  use def_parall
  implicit none
  integer(ip) :: ix,iy,iz,ipx,ipy,ipz,ielem

  npart_par=nxexp*nyexp*nzexp
  npart=npart_par
  ielem=0
  do iy=1,nyexa
     ipy=int(real(iy-1)/real(nyexa)*real(nyexp))
     do iz=1,nzexa
        ipz=int(real(iz-1)/real(nzexa)*real(nzexp))
        do ix=1,nxexa
           ielem=ielem+1
           ipx=int(real(ix-1)/real(nxexa)*real(nxexp))
           lepar_par(ielem)=ipy*nxexp*nzexp+ipz*nxexp+ipx+1
        end do
     end do
  end do

end subroutine par_exampl
