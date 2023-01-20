!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine memgen(itask,ndim1,ndim2)
  !------------------------------------------------------------------------
  !****f* memory/memgen
  ! NAME 
  !    memgen
  ! DESCRIPTION
  !    Allocate and deallocate memory for generic scalar and vector arrays:
  !    ITASK=0 ... Allocate memory
  !    ITASK=2 ... Deallocate memory
  ! USES
  ! USED BY
  !    nsi_outvar
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use mod_memchk
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask,ndim1,ndim2
  integer(4)              :: istat

  select case(itask)

  case(-1_ip)
     !
     ! Deallocate all
     !
     call memory_deallo(memke,'GESCA','memgen',gesca)
     call memory_deallo(memke,'GEVEC','memgen',gevec)
     call memory_deallo(memke,'GETEN','memgen',geten)
     call memory_deallo(memke,'GISCA','memgen',gisca)
     call memory_deallo(memke,'GIVEC','memgen',givec)

  case(0_ip)
     !
     ! Allocate memory for real 
     !
     if(ndim1>0.and.ndim2==0) then
        if ( iasca /= 0_ip ) call runend ('MEMGEN: iasca this would lead to a memory loss')
        iasca = 1
        call memory_alloca(memke,'GESCA','memgen',gesca,ndim1)
     else if(ndim1>0.and.ndim2>0) then
        if ( iavec /= 0_ip ) call runend ('MEMGEN: iavec this would lead to a memory loss')
        iavec = 1
        call memory_alloca(memke,'GEVEC','memgen',gevec,ndim1,ndim2)    
     else if(ndim1<0.and.ndim2>0) then
        if ( iaten /= 0_ip ) call runend ('MEMGEN: iaten this would lead to a memory loss')
        iaten = 1
        call memory_alloca(memke,'GETEN','memgen',geten,-ndim1,-ndim1,ndim2)        
     end if

  case(1_ip)
     !
     ! Allocate memory for integer
     !
     if(ndim1/=0.and.ndim2==0) then
        call memory_alloca(memke,'GISCA','memgen',gisca,ndim1)
     else if(ndim1/=0.and.ndim2/=0) then
        call memory_alloca(memke,'GIVEC','memgen',givec,ndim1,ndim2)     
     end if

  case(2_ip)
     !
     ! Deallocate memory for real
     !
     if(ndim1>0.and.ndim2==0) then
        if ( iasca == 0_ip ) call runend ('MEMGEN: iasca thying to deallocate smothing that is not allocated')
        call memory_deallo(memke,'GESCA','memgen',gesca)
        iasca = 0_ip
     else if(ndim1>0.and.ndim2>0) then 
        if ( iavec == 0_ip ) call runend ('MEMGEN: iavec thying to deallocate smothing that is not allocated')
        call memory_deallo(memke,'GEVEC','memgen',gevec)
        iavec = 0_ip
     else if(ndim1<0.and.ndim2>0) then 
        if ( iaten == 0_ip ) call runend ('MEMGEN: iaten thying to deallocate smothing that is not allocated')
        call memory_deallo(memke,'GETEN','memgen',geten)
        iaten = 0_ip
     end if

  case(3_ip)
     !
     ! Deallocate memory for integer
     !
     if(ndim1/=0.and.ndim2==0) then
        call memory_deallo(memke,'GISCA','memgen',gisca)
     else if(ndim1/=0.and.ndim2/=0) then 
        call memory_deallo(memke,'GIVEC','memgen',givec)        
     end if

  case(4_ip)
     !
     ! Allocate memory for r3p 
     !
     call runend('MEMGEN: NOT PROGRAMMED')

  case(5_ip)
     !
     ! Deallocate memory for r3p 
     !
     call runend('MEMGEN: NOT PROGRAMMED')

  case(6_ip)
     !
     ! Allocate memory for complex 
     !
     if(ndim1>0.and.ndim2==0) then
        if(associated(gescx)) gescx => null()
        allocate(gescx(ndim1),stat=istat)     
        call memchk(zero,istat,memke,'GESCX','memgen',gescx)
        if(istat/=0) call memerr(zero,'GESCX','memgen',0_ip)
     else if(ndim1>0.and.ndim2>0) then      
        if(associated(gevex)) gevex => null()
        allocate(gevex(ndim1,ndim2),stat=istat)     
        call memchk(zero,istat,memke,'GEVEX','memgen',gevex)        
        if(istat/=0) call memerr(zero,'GEVEX','memgen',0_ip)
     end if

  case(7_ip)
     !
     ! Deallocate memory for complex 
     !
     if(ndim1>0.and.ndim2==0) then
        deallocate(gescx,stat=istat)
        call memchk(two,istat,memke,'GESCX','memgen',gescx)
        if(istat/=0) call memerr(two,'GESCX','memgen',0_ip)
     else if(ndim1>0.and.ndim2>0) then 
        deallocate(gevex,stat=istat)
        call memchk(two,istat,memke,'GEVEX','memgen',gevex)
        if(istat/=0) call memerr(two,'GEVEX','memgen',0_ip)         
     end if

  end select

end subroutine memgen
