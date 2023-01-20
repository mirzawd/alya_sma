!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine bcsrax_schur(itask,nbnodes,nbvar,ndofn_A3,A1,A2,invA3,A4,invdiag,ja,ia,xx,yy) 
  !----------------------------------------------------------------------
  !****f* mathru/bcsrax
  ! NAME 
  !     bcsrax
  ! DESCRIPTION
  !     Multiply a non symmetric matrix stored in BCSR by a vector
  !     YY = A XX 
  ! INPUT
  !    NBNODES .... Number of equations
  !    NBVAR ...... Number of variables
  !    AN ......... Matrix
  !    JA ......... List of elements
  !    IA ......... Pointer to list of elements
  !    XX ......... Vector
  ! OUTPUT
  !    YY ......... result vector
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only             :  ip,rp
  use def_master, only             :  INOTMASTER,kfl_async,IPARALL
  use def_master, only             :  NPOIN_TYPE
  use def_solver, only             :  solve_sol

#ifdef MYTIMING
  use def_master, only             :  time_bcsrax
#endif
  implicit none
  integer(ip), intent(in)          :: itask,nbnodes,nbvar,ndofn_A3
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(inout)       :: xx(*)
  real(rp),    intent(out), target :: yy(*)
  real(rp),    intent(in)          :: A1(*)
  real(rp),    intent(in)          :: A2(ndofn_A3,*) 
  real(rp),    intent(in)          :: invA3(ndofn_A3,nbnodes)
  real(rp),    intent(in)          :: A4(ndofn_A3,*)
  real(rp),    intent(in)          :: invdiag(*)
  integer(ip)                      :: ii,jj,kk,izdom
  real(rp),    pointer             :: zz(:,:)

#ifdef MYTIMING
  real(rp)    :: time_begin, time_end                   ! Counters for timing

  ! Start timing counter
  call cputim(time_begin)
#endif

  if( IPARALL .and. kfl_async == 1 ) then

     call runend('BCSRAX_SCHUR: NOT CODED')

  else if( INOTMASTER ) then

     if( nbvar == 1 ) then
        !
        ! NBVAR=1
        !
!        do ii = 1,nbnodes
!           yy(ii) = 0.0_rp
!           if( solve_sol(1) % kfl_fixno(1,ii) /= 0 ) then
!              do izdom = ia(ii),ia(ii+1)-1
!                 jj = ja(izdom)
!                 if( ii == jj ) yy(ii) = xx(ii) * A1(izdom)
!              end do
!           else
!              do izdom = ia(ii),ia(ii+1)-1
!                 jj = ja(izdom)
!                 if( solve_sol(1) % kfl_fixno(1,jj) == 0 ) then
!                    yy(ii) = yy(ii) + 2.0_rp * A1(izdom) * xx(jj)
!                 end if
!              end do
!           end if
!        end do
!!        call bcsrax( 1_ip, nbnodes, nbvar, lapla_nsi, ja, ia, xx, yy )  
!        return

        allocate( zz(ndofn_A3,nbnodes) )
        do ii = 1,nbnodes
           do kk = 1,ndofn_A3
              zz(kk,ii) = 0
           end do
           do izdom = ia(ii),ia(ii+1)-1       
              jj = ja(izdom)
              do kk = 1,ndofn_A3
                 !zz(kk,ii) = zz(kk,ii) + A4(kk,izdom)*xx(jj)      
                 zz(kk,jj) = zz(kk,jj) + A2(kk,izdom) * xx(ii) 
              end do
           end do
        end do

        call pararr('SLX',NPOIN_TYPE,ndofn_A3*nbnodes,zz)      
  
        do ii = 1,nbnodes
           do kk = 1,ndofn_A3
              zz(kk,ii) = zz(kk,ii) * invA3(kk,ii)
           end do
        end do

        do ii = 1,nbnodes
           yy(ii) = 0.0_rp
           do izdom = ia(ii),ia(ii+1)-1
              jj = ja(izdom)
              if( solve_sol(1) % kfl_fixno(1,jj) == 0 )  then
                 yy(ii) = yy(ii) + A1(izdom) * xx(jj)
                 do kk = 1,ndofn_A3           
                    yy(ii) = yy(ii) - A2(kk,izdom) * zz(kk,jj)       
                 end do
              end if
           end do
        end do
        
        if( itask == 1 ) call pararr('SLX',NPOIN_TYPE,nbnodes,yy)

        do ii = 1,nbnodes
           if( solve_sol(1) % kfl_fixno(1,ii) /= 0 ) then
              yy(ii) = xx(ii) / invdiag(ii)
           end if
        end do

        deallocate( zz )

     else 

        call runend('BCSRAX_SCHUR: NOT CODED')

     end if

  end if

end subroutine bcsrax_schur
