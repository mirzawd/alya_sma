!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine diagox(npopo,nbvar,kfl_symme,ia,ja,an,wa1)

  !-----------------------------------------------------------------------
  ! This routine computes the diagonal of a matrix
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER,NPOIN_TYPE

  implicit none

  integer(ip), intent(in)    :: npopo,nbvar,kfl_symme
  integer(ip), intent(in)    :: ia(*),ja(*)
  complex(rp), intent(in)    :: an(nbvar,nbvar,*)
  complex(rp), intent(inout) :: wa1(*)

  integer(ip)                :: ii,jj,kk,ll

  if( INOTMASTER ) then

     if( kfl_symme == 1 ) then
 
        !Symmetric graph
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
        
        !Unsymmetric graph
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
              end if
           end do
        end if

     end if
     
     !Periodicity and Parallelization
     call pararx('SLX',NPOIN_TYPE,npopo*nbvar,wa1)
     !call rhsmox(nbvar,wa1) 

  end if

end subroutine diagox
