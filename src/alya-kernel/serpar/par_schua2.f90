!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_schua2(&
     itask,nbvar,zi,pp,aib,abi,abb,r_dom_aib,c_dom_aib,&
     r_dom_abb,c_dom_abb,r_dom_abi,c_dom_abi,zz)
  !------------------------------------------------------------------------
  !****f* Parall/par_schuax
  ! NAME 
  !    par_schuax
  ! DESCRIPTION
  !    Compute Sx where S is the Schur complement
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_master, only    :  npoi1,parr1
  use def_master, only    :  INOTMASTER
  use def_domain, only    :  npoin,npoin_bb
  use mod_memchk
  implicit none
  !
  integer(ip), intent(in)          :: itask
  integer(ip), intent(in)          :: nbvar
  real(rp),    intent(inout)       :: zi(*)
  real(rp),    intent(in)          :: pp(*)

  real(rp),    intent(in)          :: aib(nbvar,nbvar,*)
  real(rp),    intent(in)          :: abi(nbvar,nbvar,*)
  real(rp),    intent(in)          :: abb(nbvar,nbvar,*)

  integer(ip), intent(in)          :: c_dom_aib(*),r_dom_aib(*)
  integer(ip), intent(in)          :: c_dom_abi(*),r_dom_abi(*)
  integer(ip), intent(in)          :: c_dom_abb(*),r_dom_abb(*) 
  real(rp),    intent(out), target :: zz(npoin_bb)
  integer(ip)                      :: ii,jj,kk,col

  if( INOTMASTER ) then

     if( itask == 1 ) then
        !
        ! zi = Aib p
        !
        do ii = 1,npoi1
           zi(ii) = 0.0_rp
           do jj = r_dom_aib(ii),r_dom_aib(ii+1)-1
              col    = c_dom_aib(jj) 
              zi(ii) = zi(ii) + aib(1,1,jj) * pp(col)
           end do
        end do

     else if ( itask == 2 ) then
        !
        ! zz = Abb p
        !
        do ii = npoi1+1,npoin
           kk     = ii - npoi1
           zz(kk) = 0.0_rp
           do jj = r_dom_abb(kk),r_dom_abb(kk+1)-1
              col    = c_dom_abb(jj) 
              zz(kk) = zz(kk) + abb(1,1,jj) * pp(col)
           end do
        end do
        !
        ! zz = Abb p - Abi ( Aii^{-1} Aib p )
        !
        do ii = npoi1+1,npoin
           kk = ii - npoi1
           do jj = r_dom_abi(kk),r_dom_abi(kk+1)-1
              col    = c_dom_abi(jj) 
              zz(kk) = zz(kk) - abi(1,1,jj) * zi(col)
           end do
        end do
        
        parr1 => zz
        call par_slesch()         
        !call par_pararr('SCH',0_ip,npoin-npoi1,zz)
     end if

  end if

end subroutine par_schua2
