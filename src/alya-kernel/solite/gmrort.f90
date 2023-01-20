!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine gmrort(nbvar,nbnodes,jj,jj1,kryl,hh)
  !------------------------------------------------------------------------
  !****f* solite/gmrort
  ! NAME 
  !    prodxy
  ! DESCRIPTION
  !    This routine computes the scalar product between the
  !    Krylov vectors before orthogonalization
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr
  implicit none
  integer(ip), intent(in)           :: nbvar,nbnodes,jj,jj1
  real(rp),    intent(in)           :: kryl(nbvar*nbnodes,*)
  real(rp),    target,  intent(out) :: hh(jj)
  integer(ip)                       :: kk,ii,mm1,mm2,mm3

  do ii=1,jj           
     hh(ii) = 0.0_rp
  end do

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do ii=1,jj              
        do kk=1,nbvar*nbnodes
           hh(ii) = hh(ii) + kryl(kk,ii) * kryl(kk,jj1)
        end do
     end do
     
  else if(kfl_paral>=1) then
     !
     ! Parallel: Slaves
     !
     mm1=nbvar*npoi1
     mm2=(npoi2-1)*nbvar+1
     mm3=npoi3*nbvar
     do ii=1,jj              
        do kk=1,mm1
           hh(ii) = hh(ii) + kryl(kk,ii) * kryl(kk,jj1)
        end do
        do kk=mm2,mm3
           hh(ii) = hh(ii) + kryl(kk,ii) * kryl(kk,jj1)
        end do
     end do
  end if

  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     nparr =  jj
     parre => hh
     call par_operat(3_ip)   
  end if

end subroutine gmrort

 
