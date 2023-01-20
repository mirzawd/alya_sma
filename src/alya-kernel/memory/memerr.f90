!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine memerr(itask,vanam,vacal,istat)
  !-----------------------------------------------------------------------
  !****f* memory/memerr
  ! NAME 
  !    memctr
  ! DESCRIPTION
  !    This routine ends Alya when an error has been found
  !    allocating or deallocating memory.
  ! USES
  !    runend
  ! USED BY
  !    Mod_memchk
  !***
  !-----------------------------------------------------------------------
  use def_master
  use mod_memory, only : lbytm
  implicit none
  integer(ip),   intent(in) :: itask
  integer(4),    intent(in) :: istat
  integer(ip)               :: ibyte
  real(rp)                  :: rbyte
  character(6)              :: lbyte
  character*(*), intent(in) :: vanam,vacal
  character(200)            :: wmess
  character(20)             :: wmes2,wmes3

  if(itask==0) then
     !
     ! Allocation
     !
     if(lbytm>=1024*1024*1024) then
        rbyte=1024.0_rp*1024.0_rp*1024.0_rp
        lbyte='Gbytes'
     else if(lbytm>=1024*1024) then 
        rbyte=1024.0_rp*1024.0_rp
        lbyte='Mbytes'     
     else if(lbytm>=1024) then 
        rbyte=1024.0_rp
        lbyte='kbytes'          
     else  
        rbyte=1.0_rp
        lbyte=' bytes'     
     end if 
     ibyte=int(real(lbytm,rp)/rbyte,KIND=ip)
     wmes2=intost(ibyte)
     wmes3=intost(int(istat,ip))
     wmess=trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE ALLOCATED.'&
          //' RUN TIME ERROR: '//trim(wmes3)
     call runend(trim(wmess))

  else if(itask==1) then
     !
     ! Reallocation
     !
     call runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE REALLOCATED')

  else
     !
     ! Deallocation
     !
     call runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE DEALLOCATED')

  end if

end subroutine memerr
