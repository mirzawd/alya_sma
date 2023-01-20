!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




#ifdef NINJA
module pinnedmemorygpu
  integer(8) :: c_amatr,c_rhsid,c_unkno  
end module pinnedmemorygpu
#endif

subroutine memunk(itask)

  !-----------------------------------------------------------------------
  !****f* master/memunk
  ! NAME
  !    Turnon
  ! DESCRIPTION
  !    This routine allocates memory for all the unknowns of the problem
  !    When using Parall, Master allocates minimum memory
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
#ifdef NINJA
  use pinnedmemorygpu
#endif
  
  use def_parame
  use def_master 
  use def_kermod
  use def_domain
  use def_solver
  use mod_memory
  use mod_domain, only : domain_memory_allocate
#ifdef NINJA
  use iso_c_binding
#endif
  implicit none
  integer(ip), intent(in) :: itask

#ifdef NINJA
  integer(8)  :: d_temp,sz
  type(c_ptr) :: cptr
#endif
  !
  ! Deallocate before allocating
  !
  if( itask == 2 ) then
  end if
  !
  ! Algebraic solver
  !
#ifdef NINJA
  sz = nzmat
  sz = sz * 8

!  call gpumallochost(d_temp,sz)
!  call int2cptr(d_temp,cptr)
!  call c_f_pointer(cptr,amatr,[nzmat])
!  c_amatr = d_temp
  call memory_alloca(memma,'AMATR','memunk',amatr,max(1_ip,nzmat))

  sz = nzrhs
  sz = sz * 8

!  call gpumallochost(d_temp,sz)
!  call int2cptr(d_temp,cptr)
!  call c_f_pointer(cptr,rhsid,[nzrhs])
!  c_rhsid = d_temp
!  call gpumallochost(d_temp,sz)
!  call int2cptr(d_temp,cptr)
!  call c_f_pointer(cptr,unkno,[nzrhs])
!  c_unkno = d_temp
  call memory_alloca(memma,'RHSID','memunk',rhsid,max(1_ip,nzrhs))
  call memory_alloca(memma,'UNKNO','memunk',unkno,max(1_ip,nzrhs))

#else
  
  call memory_alloca(memma,'AMATR','memunk',amatr,max(1_ip,nzmat))
  call memory_alloca(memma,'RHSID','memunk',rhsid,max(1_ip,nzrhs))
  call memory_alloca(memma,'UNKNO','memunk',unkno,max(1_ip,nzrhs))

#endif

  call memory_alloca(memma,'BMATR','memunk',bmatr,max(1_ip,nzmbt))
  call memory_alloca(memma,'PMATR','memunk',pmatr,max(1_ip,nzpre))
  call memory_alloca(memma,'ERRES','memunk',erres,max(1_ip,nzerr))
  !
  ! Eigen solver
  !
  call memory_alloca(memma,'EIGEN','memunk',eigen,max(1_ip,neige))
  call memory_alloca(memma,'EIGVA','memunk',eigva,max(1_ip,neiva))
  !
  ! Algebraic complex solver 
  !
  call memory_alloca(memma,'AMATX','memunk',amatx,max(1_ip,nzmax))
  call memory_alloca(memma,'RHSIX','memunk',rhsix,max(1_ip,nzrhx))
  call memory_alloca(memma,'UNKNX','memunk',unknx,max(1_ip,nzrhx))
  call memory_alloca(memma,'PMATX','memunk',pmatx,max(1_ip,nzprx))
  !
  ! lumpped mass matrix - for use in dual time step preconditioner
  !
   if( INOTMASTER ) call domain_memory_allocate('LUMMA') ! not sure if this is the optimal place
   
end subroutine memunk
