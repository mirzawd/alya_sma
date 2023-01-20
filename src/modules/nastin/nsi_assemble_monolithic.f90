!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_assemble_monolithic(&
     pnode,pevat,lnods,elauu,elaup,elapp,elapu,&
     elrbu,elrbp,A,b)
  !-----------------------------------------------------------------------
  !****f* mathru/assma3
  ! NAME 
  !    assma3
  ! DESCRIPTION
  !    Assembly an elemental matrix ELMAT in global matrix AMATR
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  nzdom,ndime,r_sol,c_sol
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elauu(pevat,pevat)
  real(rp),    intent(in)    :: elaup(pevat,pnode)
  real(rp),    intent(in)    :: elapp(pnode,pnode)
  real(rp),    intent(in)    :: elapu(pnode,pevat)
  real(rp),    intent(in)    :: elrbu(ndime,pnode)
  real(rp),    intent(in)    :: elrbp(pnode)
  real(rp),    intent(inout) :: A(ndime+1,ndime+1,nzdom)  
  real(rp),    intent(inout) :: b(ndime+1,*)
  integer(ip)                :: ndofn,inode,jnode,iposi,jposi
  integer(ip)                :: idime,jdime,izsol,jpoin,ipoin,jcolu

  ndofn = ndime + 1

  do inode = 1,pnode
     ipoin = lnods(inode)
     !
     ! Matrix
     !
     do jnode = 1,pnode
        jpoin = lnods(jnode)
        izsol = r_sol(ipoin)
        jcolu = c_sol(izsol)
        do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
           izsol = izsol + 1
           jcolu = c_sol(izsol)
        end do
        if( jcolu == jpoin ) then
           do idime = 1,ndime
              iposi = (inode-1) * ndime + idime 
              do jdime = 1,ndime
                 jposi                = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                 !$OMP ATOMIC
#endif
                 A(jdime,idime,izsol) = A(jdime,idime,izsol) + elauu(iposi,jposi)   ! Auu
              end do
#ifdef NO_COLORING
              !$OMP ATOMIC
#endif
              A(ndofn,idime,izsol) = A(ndofn,idime,izsol) + elaup(iposi,jnode)      ! Aup
           end do
           do jdime=1,ndime
              jposi = (jnode-1) * ndime + jdime
#ifdef NO_COLORING
              !$OMP ATOMIC
#endif
              A(jdime,ndofn,izsol) = A(jdime,ndofn,izsol) + elapu(inode,jposi)      ! Apu                 
           end do
#ifdef NO_COLORING
           !$OMP ATOMIC
#endif
           A(ndofn,ndofn,izsol) = A(ndofn,ndofn,izsol) + elapp(inode,jnode)         ! App
        end if
     end do
     !
     ! RHS
     !
     do idime = 1,ndime
#ifdef NO_COLORING
        !$OMP ATOMIC
#endif
        b(idime,ipoin) = b(idime,ipoin) + elrbu(idime,inode)
     end do
#ifdef NO_COLORING
     !$OMP ATOMIC
#endif
     b(ndofn,ipoin)   = b(ndofn,ipoin)   + elrbp(inode)
  end do


end subroutine nsi_assemble_monolithic
