!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine limite(nbvar,ia,ja,an)  
  !-----------------------------------------------------------------------
  ! 
  ! Algebraic limiter
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER
  use def_solver, only       :  solve_sol
  use def_domain, only       :  npoin
  implicit none
  integer(ip), intent(in)    :: nbvar 
  real(rp),    intent(inout) :: an(nbvar,nbvar,*)
  integer(ip), intent(in)    :: ja(*),ia(*)
  integer(ip)                :: ii,jj,kk,ll,nn

  if( INOTMASTER .and. solve_sol(1)%kfl_limit == 1 ) then

     if( nbvar == 1 ) then
        !
        ! NBVAR = 1
        !
        do ii= 1, npoin 
           jj = ia(ii)
           ll = -1
           do while ( jj < ia(ii+1) .and. ll == -1 )
              if(ja(jj)==ii) then
                 ll = jj
              end if
              jj = jj + 1
           end do
           jj = ia(ii)
           do while ( jj < ia(ii+1) )
              if( jj /= ll .and. an(1,1,jj) > 0.0_rp ) then
                 an(1,1,ll) = an(1,1,ll) + an(1,1,jj)
                 an(1,1,jj) = 0.0_rp
              end if
              jj = jj + 1
           end do
        end do

     else
        !
        ! NBVAR > 1
        !
        do ii= 1, npoin 
           
           jj = ia(ii)
           ll = -1
           do while ( jj < ia(ii+1) .and. ll == -1 )
              if(ja(jj)==ii) then
                 ll = jj
              end if
              jj = jj + 1
           end do

           jj = ia(ii)
           do while ( jj < ia(ii+1) )
              if( jj == ll ) then
                 do kk = 1,nbvar
                    do nn = 1,nbvar
                       if( an(kk,nn,ll) > 0.0_rp .and. kk /= nn ) then
                          an(nn,nn,ll) = an(nn,nn,ll) + an(kk,nn,ll)
                          an(kk,nn,ll) = 0.0_rp
                       end if
                    end do
                 end do
              else
                 do kk = 1,nbvar
                    do nn = 1,nbvar
                       if( an(kk,nn,jj) > 0.0_rp ) then
                          an(nn,nn,ll) = an(nn,nn,ll) + an(kk,nn,jj)
                          an(kk,nn,jj) = 0.0_rp
                       end if
                    end do
                 end do
              end if
              jj = jj + 1
           end do

        end do

     end if
  end if

end subroutine limite
