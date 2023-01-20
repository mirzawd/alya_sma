!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




module mod_mag_wrdata

!##################################################################
! Output parallel_data_structures with matrix_market branch
!
! Open .dat file and add the following lines:
!
! PARALL_SERVICE: On
! OUTPUT
!   NODE_COMMUNICATION_ARRAYS:    On
!   EDGE_COMMUNICATION_ARRAYS:    On
!   FULL_ROWS_COMMUNICATION_ARRAYS:    On
! END_OUTPUT
!
! Adding the following lines in .ker.dat is not needed:
!
! OUTPUT_&_POST_PROCESS
!   ON_LAST_MESH:    On
!   EDGE_GRAPH:    On
!   STEPS = 1
! END_OUTPUT_&_POST_PROCESS
!
!##################################################################

  use def_magnet
  use def_master, only: solve_sol, amatr, rhsid, kfl_paral
  use mod_strings, only: integer_to_string
  use def_domain, only: meshe, ndime
  use def_kermod, only: ndivi
  use mod_mag_nedele
  use mod_mag_lagran
  
  implicit none

contains


  !##############################################################
  subroutine mag_wrsign()
  
    integer(ip) :: ielem, eledof, pelty, idof

    if (kfl_paral /= 0_ip) then

      open(unit=ioun6_mag, file='magnet_signs-'//integer_to_string(kfl_paral)//'.dat')

      do ielem = 1_ip, meshe(ndivi) % nelem

        pelty = meshe(ndivi) % ltype(ielem)
        eledof = mag_edgdof(pelty)

        do idof = 1_ip, eledof - 1_ip
          write(ioun6_mag, '(F8.3)', advance='no') elesig_mag(idof, ielem)
        end do

        if (eledof > 0_ip) write(ioun6_mag, '(F8.3)') elesig_mag(eledof, ielem)

      end do

      close(ioun6_mag)

    end if

  end subroutine mag_wrsign
  !##############################################################


  !##############################################################
  subroutine mag_wreled()

    integer(ip) :: ielem, eledof, pelty, idof

    if (kfl_paral /= 0_ip) then

      open(unit=ioun6_mag, file='magnet_eledg-'//integer_to_string(kfl_paral)//'.dat')

      do ielem = 1_ip, meshe(ndivi) % nelem

        pelty = meshe(ndivi) % ltype(ielem)
        eledof = mag_edgdof(pelty)

        do idof = 1_ip, eledof - 1_ip
          write(ioun6_mag, '(I8)', advance='no') eleedg_mag(idof, ielem)
        end do

        if (eledof > 0_ip) write(ioun6_mag, '(I8)') eleedg_mag(eledof, ielem)

      end do

      close(ioun6_mag)

    end if

  end subroutine mag_wreled
  !##############################################################


  !##############################################################
  subroutine mag_wrmatr()

     integer(ip) :: ii, jj, jj1, jj2

     if (kfl_paral /= 0_ip) then

      open(unit=ioun6_mag, file='magnet_aa-'//integer_to_string(kfl_paral)//'.dat')
      open(unit=ioun7_mag, file='magnet_ja-'//integer_to_string(kfl_paral)//'.dat')

      do ii = 1_ip, solve_sol(1) % nzmat
        write(ioun6_mag, '(e15.8)') amatr(ii)
        write(ioun7_mag, '(I10)') solve_sol(1) % ja(ii)
      end do

      close(ioun6_mag)
      close(ioun7_mag)

      open(unit=ioun6_mag, file='magnet_rhs-'//integer_to_string(kfl_paral)//'.dat')
      open(unit=ioun7_mag, file='magnet_ia-'//integer_to_string(kfl_paral)//'.dat')

      do ii = 1_ip, solve_sol(1) % nzrhs
        write(ioun6_mag, '(e15.8)') rhsid(ii)
        write(ioun7_mag, '(I10)') solve_sol(1) % ia(ii)
      end do
      write(ioun7_mag, '(I10)') solve_sol(1) % ia(solve_sol(1) % nzrhs + 1_ip)        

      close(ioun6_mag)
      close(ioun7_mag)

      open(unit=ioun6_mag, file='magnet_ii-'//integer_to_string(kfl_paral)//'.dat')
      open(unit=ioun7_mag, file='magnet_jj-'//integer_to_string(kfl_paral)//'.dat')               
      
      do ii = 1_ip, solve_sol(1) % nzrhs
        jj1 = solve_sol(1) % ia(ii)
        jj2 = solve_sol(1) % ia(ii + 1_ip) - 1_ip
        do jj = jj1, jj2
          write(ioun6_mag, '(I10)') ii
          write(ioun7_mag, '(I10)') solve_sol(1) % ja(jj)
        end do
      end do

      close(ioun6_mag)
      close(ioun7_mag)

    end if

  end subroutine mag_wrmatr
  !##############################################################


  !##############################################################
  subroutine mag_wrnedg()

!    integer(ip) :: &
!        node_to_nedge(meshe(ndivi) % npoin),    &
!        r_node_to_edge(meshe(ndivi) % npoin + 1_ip),    &
!        c_node_to_edge(meshe(ndivi) % nedge * 2_ip),    &
!        i_node_to_edge(meshe(ndivi) % npoin + 1_ip),    &
!        edge_to_nedge(meshe(ndivi) % nedge),    &
!        r_edge_to_edge(meshe(ndivi) % nedge + 1_ip),    &
!        i_edge_to_edge(meshe(ndivi) % nedge + 1_ip)
!
!    integer(ip), allocatable :: &
!        c_edge_to_edge(:)

    integer(ip) :: iedge,    &
                   jedge,    &
                   ipoin,    &
                   jpoin
!                   ki,    &
!                   kj,    &
!                   ne2e
!    !
!    ! Get number of edges per node:
!    ! node_to_nedge has size npoin
!    !
!    node_to_nedge = 0_ip
!    do iedge = 1, meshe(ndivi) % nedge
!      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
!      jpoin = meshe(ndivi) % edge_to_node(2, iedge)
!
!      node_to_nedge(ipoin) = node_to_nedge(ipoin) + 1_ip
!      node_to_nedge(jpoin) = node_to_nedge(jpoin) + 1_ip
!    end do
!    !
!    !
!    ! Get edges per node:
!    ! r_node_to_node has size npoin+1
!    ! c_node_to_node has size 2*nedge
!    ! i_node_to_node has size npoin+1
!    !
!    r_node_to_edge(1) = 1_ip
!    do ipoin = 1, meshe(ndivi) % npoin
!      r_node_to_edge(ipoin + 1) = r_node_to_edge(ipoin) + node_to_nedge(ipoin)
!    end do
!    i_node_to_edge = r_node_to_edge
!    !
!    do iedge = 1, meshe(ndivi) % nedge
!      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
!      jpoin = meshe(ndivi) % edge_to_node(2, iedge)
!
!      ki = i_node_to_edge(ipoin)
!      c_node_to_edge(ki) = iedge
!      i_node_to_edge(ipoin) = ki + 1_ip
!
!      kj = i_node_to_edge(jpoin)
!      c_node_to_edge(kj) = iedge
!      i_node_to_edge(jpoin) = kj + 1_ip
!    end do
!    !
!    !
!    ! Get number of edges per edge
!    ! edge_to_nedge has size nedge
!    !
!    ne2e = 0
!    do iedge = 1, meshe(ndivi) % nedge
!      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
!      jpoin = meshe(ndivi) % edge_to_node(2, iedge)
!
!      edge_to_nedge(iedge) = r_node_to_edge(ipoin + 1) - r_node_to_edge(ipoin) - 1_ip
!      edge_to_nedge(iedge) = edge_to_nedge(iedge) + r_node_to_edge(jpoin + 1) - r_node_to_edge(jpoin) - 1_ip
!
!      ne2e = ne2e + edge_to_nedge(iedge)
!    end do
!    !
!    allocate(c_edge_to_edge(ne2e))
!    !
!    !
!    ! Get edges per edge:
!    ! r_edge_to_edge has size nedge+1
!    ! c_edge_to_edge has size ne2e
!    ! i_edge_to_edge has size nedge+1
!    !
!    r_edge_to_edge(1) = 1_ip
!    do iedge = 1, meshe(ndivi) % nedge
!      r_edge_to_edge(iedge + 1) = r_edge_to_edge(iedge) + edge_to_nedge(iedge)
!    end do
!    i_edge_to_edge = r_edge_to_edge
!    !
!    do iedge = 1, meshe(ndivi) % nedge
!      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
!      jpoin = meshe(ndivi) % edge_to_node(2, iedge)
!
!      do jedge = r_node_to_edge(ipoin), r_node_to_edge(ipoin+1) - 1_ip
!        if (c_node_to_edge(jedge) /= iedge) then
!          c_edge_to_edge(i_edge_to_edge(iedge)) = c_node_to_edge(jedge)
!          i_edge_to_edge(iedge) = i_edge_to_edge(iedge) + 1_ip 
!        end if
!      end do
!
!      do jedge = r_node_to_edge(jpoin), r_node_to_edge(jpoin+1) - 1_ip
!        if (c_node_to_edge(jedge) /= iedge) then
!          c_edge_to_edge(i_edge_to_edge(iedge)) = c_node_to_edge(jedge)
!          i_edge_to_edge(iedge) = i_edge_to_edge(iedge) + 1_ip
!        end if
!      end do
!    end do
!    !
!    !
    if (kfl_paral /= 0_ip) then

      open(unit=ioun6_mag, file='magnet_iee-'//integer_to_string(kfl_paral)//'.dat')
      open(unit=ioun7_mag, file='magnet_jee-'//integer_to_string(kfl_paral)//'.dat')

      do iedge = 1_ip, meshe(ndivi) % nedge + 1_ip
        write(ioun6_mag, '(I10)') r_edge_to_edge_mag(iedge)
      end do

      do jedge = 1_ip, ne2e_mag
        write(ioun7_mag, '(I10)') c_edge_to_edge_mag(jedge)
      end do

      close(ioun6_mag)
      close(ioun7_mag)


      open(unit=ioun6_mag, file='magnet_ine-'//integer_to_string(kfl_paral)//'.dat')
      open(unit=ioun7_mag, file='magnet_jne-'//integer_to_string(kfl_paral)//'.dat')

      do ipoin = 1_ip, meshe(ndivi) % npoin + 1_ip
        write(ioun6_mag, '(I10)') r_node_to_edge_mag(ipoin)
      end do

      do jpoin = 1_ip, meshe(ndivi) % nedge * 2_ip
        write(ioun7_mag, '(I10)') c_edge_to_edge_mag(jpoin)
      end do

      close(ioun6_mag)
      close(ioun7_mag)

    end if
    !
!    deallocate(c_edge_to_edge_mag)
    !
  end subroutine mag_wrnedg
  !##############################################################


  !##############################################################
  subroutine mag_wredgn()
    !
    implicit none
    !
    integer(ip) :: &
      iedge
    !
    if (kfl_paral /= 0_ip) then
      !
      open(unit=ioun6_mag, file='magnet_ien-'//integer_to_string(kfl_paral)//'.dat')
      !
      do iedge = 1_ip, meshe(ndivi) % nedge
        write(ioun6_mag, '(2I10)') edgnod_mag(1:2, iedge)
      end do
      !
      close(ioun6_mag)
      !  
    end if
    !
  end subroutine mag_wredgn
  !##############################################################


  !##############################################################
  subroutine mag_wrdata()
    !
    call mag_wrsign()
    call mag_wreled()
    call mag_wrmatr()
    call mag_wrnedg()
    call mag_wredgn()
    !
  end subroutine mag_wrdata
  !##############################################################

end module mod_mag_wrdata
