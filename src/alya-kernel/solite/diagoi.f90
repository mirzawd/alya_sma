!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{                                                                   
!> @file    diagoi.f90
!> @author  Mariano Vazquez
!> @date    12/01/2017
!> @brief   Compute indices for the diagonal of amatr
!> @details Compute indices for the diagonal of amatr
!> @}
!-----------------------------------------------------------------------
subroutine diagoi(npopo,nbvar,kfl_symme,ia,ja,iwa1)
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER,NPOIN_TYPE

  implicit none

  integer(ip), intent(in)    :: npopo,nbvar,kfl_symme
  integer(ip), intent(in)    :: ia(*),ja(*)
  integer(ip), intent(inout) :: iwa1(*)

  integer(ip)                :: ii,jj,kk,ll

  if( INOTMASTER ) then

     if( kfl_symme == 1 ) then
 
        !Symmetric graph
        if( nbvar == 1 ) then
           do ii= 1, npopo
              ll = ia(ii+1)-1
              iwa1(ii) = ll 
           end do
        else
           do ii= 1, npopo
              ll = ia(ii+1)-1
              jj = (ii-1) * nbvar 
              do kk= 1, nbvar
                 iwa1(jj+kk) = ll
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
                 iwa1(ii)= ll
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
                    iwa1(jj+kk)= ll
                 end do
              end if
           end do
        end if

     end if
     
  end if

end subroutine diagoi
