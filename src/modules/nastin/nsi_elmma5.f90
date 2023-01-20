!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmma5(&
     pnode,pgaus,pevat,gpden,&
     gpvol,gpsha,dtinv_loc,elcmm)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmma5
  ! NAME 
  !    nsi_elmma5
  ! DESCRIPTION
  !    Compute consistent element mass matrix - created from nsi_elmma4
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
!  use def_nastin, only       :  pabdf_nsi
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: dtinv_loc
  real(rp),    intent(out)   :: elcmm(pnode*ndime,pnode*ndime)
  
  integer(ip)                :: inode,jnode,idofn,jdofn
  integer(ip)                :: idof1,idof2,idof3,jdof1,jdof2,jdof3
  integer(ip)                :: igaus
  real(rp)                   :: fact4,fact5,fact8

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  do idofn = 1,pevat
     do jdofn = 1,pevat
        elcmm(jdofn,idofn) = 0.0_rp
     end do
  end do
  !
  ! Mass matrix only Galerkin ( * rho * pabdf(1) / dt)
  !
  do igaus = 1,pgaus

     fact8 = gpden(igaus)*dtinv_loc ! Oriol test consitent matrix + RK 
 !    fact8 = pabdf_nsi(1) * gpden(igaus) * dtinv_loc

     do inode = 1,pnode

        idof1 = ndime*(inode-1)+1
        idof2 = ndime*(inode-1)+2
        if (ndime==3) idof3 = ndime*(inode-1)+3

        fact4 = gpsha(inode,igaus) * gpvol(igaus)

        do jnode = 1,pnode    

           jdof1 = ndime*(jnode-1)+1
           jdof2 = ndime*(jnode-1)+2
           if (ndime==3) jdof3 = ndime*(jnode-1)+3

           fact5 = fact4 * ( fact8 * gpsha(jnode,igaus) )           ! ( rho/dt N_j ) Ni

           elcmm(idof1,jdof1) = elcmm(idof1,jdof1) + fact5 
           elcmm(idof2,jdof2) = elcmm(idof2,jdof2) + fact5
           if (ndime==3) elcmm(idof3,jdof3) = elcmm(idof3,jdof3) + fact5

        end do
     end do
  end do

end subroutine nsi_elmma5

