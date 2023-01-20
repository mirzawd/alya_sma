!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine diagon(npopo,nbvar,kfl_symme,kfl_full_rows,ia,ja,an,wa1)
  !-----------------------------------------------------------------------
  !
  ! Compute the diagonal
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER,NPOIN_TYPE
  implicit none
  integer(ip), intent(in)    :: npopo,nbvar,kfl_symme,kfl_full_rows
  integer(ip), intent(in)    :: ia(*),ja(*)
  real(rp),    intent(in)    :: an(nbvar,nbvar,*)
  real(rp),    intent(inout) :: wa1(*)
  integer(ip)                :: ii,jj,kk,ll

  if( INOTMASTER ) then

     if( kfl_symme == 1 ) then
        !
        ! Symmetric graph
        !
        if( nbvar == 1 ) then
           do ii= 1, npopo
              ll = ia(ii+1)-1
              wa1(ii) = an(1,1,ll)
           end do
        else
           do ii= 1, npopo
              ll = ia(ii+1)-1
              jj = (ii-1) * nbvar 
              do kk= 1, nbvar
                 wa1(jj+kk) = an(kk,kk,ll)
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
                 wa1(ii)=an(1,1,ll)
              else
                 wa1(ii)=0.0_rp
                 !print*,'no diagonal=',ii
                 !stop
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
                 jj = (ii-1) * nbvar
                 do kk= 1, nbvar
                    wa1(jj+kk)=an(kk,kk,ll)
                 end do
              else
                 jj = (ii-1) * nbvar
                 do kk= 1, nbvar
                    wa1(jj+kk)=0.0_rp
                 end do
              end if
           end do
        end if

     end if
     !
     ! Periodicity and Parallelization
     !
     if( kfl_full_rows == 0 ) call pararr('SLX',NPOIN_TYPE,nbvar*npopo,wa1)
     !call pararr2('SLX',NPOIN_TYPE,nbvar*npopo,wa1)
     !wa1(1:npopo) = 1.0_rp

  end if

end subroutine diagon
