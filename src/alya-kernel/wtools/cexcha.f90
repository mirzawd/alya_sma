!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine cexcha(nsize,chava)
  !-----------------------------------------------------------------------
  !****f* cexcha
  ! NAME
  !    cexcha
  ! DESCRIPTION
  !    This routine exchange character data individually
  ! USES
  ! USED BY
  !    *_sendat
  !    *_parall
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  implicit none
  integer(ip), intent(in) :: nsize
  character(*)            :: chava
  integer(ip)             :: iposi,parch_size

  iposi = nparc
  nparc = nparc + nsize
  parch_size = len(parch)

  if( parii == 2 ) then

     if( kfl_ptask == 1 ) then

        if( INOTSLAVE ) then
           if( iposi+nsize > parch_size ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
           parch(iposi+1:iposi+nsize) = chava(1:nsize)
        else
           if( iposi+nsize > parch_size ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
           chava(1:nsize) = parch(iposi+1:iposi+nsize)
        end if

     else if( kfl_ptask == 0 ) then

        if( INOTSLAVE ) then
           if( iposi+nsize > parch_size ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
           parch(iposi+1:iposi+nsize) = chava(1:nsize)
        end if

     else if( kfl_ptask == 2 ) then

        if( iposi+nsize > parch_size ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED ')
        chava(1:nsize) = parch(iposi+1:iposi+nsize)

     end if

  end if

end subroutine cexcha
