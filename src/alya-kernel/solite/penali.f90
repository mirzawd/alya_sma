!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine penali(nbvar,npopo,ia,ja,an,bb,xx)
  use def_kintyp, only       :  ip,rp
  use def_solver, only       :  solve_sol
  use def_master, only       :  IMASTER,NPOIN_TYPE
  implicit none
  integer(ip), intent(in)    :: nbvar,npopo
  integer(ip), intent(in)    :: ia(*)
  integer(ip), intent(in)    :: ja(*)
  real(rp),    intent(inout) :: an(nbvar,nbvar,*)
  real(rp),    intent(inout) :: bb(nbvar,*)
  real(rp),    intent(in)    :: xx(nbvar,*)
  real(rp),    pointer       :: diag(:,:)
  integer(ip)                :: ii,jj,kk,ll
  real(rp)                   :: penal

  if( IMASTER ) return

  if( solve_sol(1) % kfl_penal == 1 ) then

     if( solve_sol(1) % kfl_assem == 1 ) then

        solve_sol(1) % kfl_assem = 2
        penal = solve_sol(1) % penal
        allocate(diag(nbvar,npopo))

        if( solve_sol(1) % kfl_symme == 1 ) then
           !
           ! Symmetric graph
           !
           if( nbvar == 1 ) then
              do ii= 1, npopo
                 ll = ia(ii+1)-1
                 diag(1,ii) = an(1,1,ll)
              end do
           else
              do ii= 1, npopo
                 ll = ia(ii+1)-1
                 do kk= 1, nbvar
                    diag(kk,ii) = an(kk,kk,ll)
                 end do
              end do
           end if

        else
           !
           ! Unsymmetric graph
           !
           if( nbvar == 1 ) then
              do ii= 1, npopo 
                 jj = ia(ii)
                 ll = -1
                 do while (jj< ia(ii+1) .and. ll ==-1)
                    if(ja(jj)==ii) then
                       ll = jj
                    end if
                    jj = jj+1
                 end do
                 if(ll/=-1) then
                    diag(1,ii)=an(1,1,ll)
                 end if
              end do
           else
              do ii= 1, npopo 
                 jj = ia(ii)
                 ll = -1
                 do while (jj< ia(ii+1) .and. ll ==-1)
                    if(ja(jj)==ii) then
                       ll = jj
                    end if
                    jj = jj+1
                 end do
                 if(ll/=-1) then
                    do kk= 1, nbvar
                       diag(kk,ii)=an(kk,kk,ll)
                    end do
                 end if
              end do
           end if

        end if

        do ii =  1,npopo
           jj =  ia(ii)
           ll = -1
           do while (jj< ia(ii+1) .and. ll ==-1)
              if(ja(jj)==ii) then
                 ll = jj
              end if
              jj = jj + 1
           end do
           !
           ! A = A + diag(A) * eps
           !
           if( ll /= -1 ) then
              do kk = 1,nbvar
                 an(kk,kk,ll) = an(kk,kk,ll) + diag(kk,ii) * penal
              end do
           end if
        end do
        !
        ! b = b + diag(A) * eps * u^i-1
        !
        call pararr('SLX',NPOIN_TYPE,nbvar*npopo,diag)
        do ii =  1,npopo
           jj =  ia(ii)
           ll = -1
           do while (jj< ia(ii+1) .and. ll ==-1)
              if(ja(jj)==ii) then
                 ll = jj
              end if
              jj = jj + 1
           end do
           if( ll /= -1 ) then
              do kk = 1,nbvar
                 bb(kk,ii) = bb(kk,ii) + diag(kk,ii) * penal * xx(kk,ii)
              end do
           end if
        end do
        !
        ! Deallocate memory
        !
        deallocate(diag)

     end if

  end if

end subroutine penali
