!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine assgr2(ngrou,npopo,nskyl,ndofn,ia,ja,an,askyl)
  !
  ! ASKYL: Factorize group matrix
  !
  use def_kintyp, only               :  ip,rp
  use def_solver, only               :  solve_sol
  use mod_memchk
  implicit none
  integer(ip), intent(in)            :: ngrou,npopo,nskyl,ndofn
  integer(ip), intent(in)            :: ia(*),ja(*)
  real(rp),    intent(in)            :: an(ndofn,ndofn,*)
  real(rp),    intent(inout), target :: askyl(ndofn,ndofn,*)
  integer(ip), pointer               :: lgrou(:)
  integer(ip)                        :: igrou,jgrou,ipoin,izdom,jpoin,izgro
  integer(ip)                        :: ii,jj
  integer(ip), pointer               :: iagro(:),jagro(:)

  !----------------------------------------------------------------
  !
  ! Fill in sparse matrix ASKYL 
  !       
  !----------------------------------------------------------------
  
  lgrou => solve_sol(1) % lgrou
  iagro => solve_sol(1) % iagro
  jagro => solve_sol(1) % jagro
  
  do ipoin = 1,npopo
     if( lgrou(ipoin) > 0 ) then
        igrou = lgrou(ipoin)
        do izdom = ia(ipoin),ia(ipoin+1)-1
           jpoin = ja(izdom)
           if( lgrou(jpoin) > 0 ) then
              jgrou = lgrou(jpoin)
              izgro = iagro(igrou)
              iifzgro1: do while( jagro(izgro) /= jgrou )
                 izgro = izgro + 1
              end do iifzgro1
              do ii = 1,ndofn
                 do jj = 1,ndofn
                    askyl(jj,ii,izgro) = askyl(jj,ii,izgro) + an(jj,ii,izdom)
                 end do
              end do
           end if
        end do
     end if
  end do
end subroutine assgr2
