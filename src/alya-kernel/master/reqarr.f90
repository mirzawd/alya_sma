!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine reqarr()
  !-----------------------------------------------------------------------
  !****f* master/reqarr
  ! NAME
  !    reqarr
  ! DESCRIPTION
  !    This routine determines the flags for the list of required arrays.
  !    They are:
  !
  !    KFL_LFACG ... LFACG
  !    KFL_LELBF ... LELBF
  !    KFL_LELP2 ... PELPO_2, LELPO_2
  !    KFL_LELE2 ... PELEL_2, LELEL_2
  !    KFL_SYMGR ... R_SYM, C_SYM
  !    KFL_SCHUR ... NZDOM_*, R_DOM_*, C_DOM_*. * = AII,AIB,ABI,ABB,PREC
  !    KFL_AIIPR ... NZDOM_AII, R_DOM_AII, C_DOM_AII
  !
  !    There can be some dependences. For example, if faces are needed
  !    (KFL_FFACG=1) then extended graph is needed (KFL_LELP2=1)
  !
  ! USES
  ! USED BY
  !    Turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_inpout
  use def_domain
  use def_kermod
  implicit none
  integer(ip) :: ivari,imodu 
  !
  ! List of faces: LFACG
  !
  if( kfl_modul(ID_PARTIS) /= 0 ) then
     kfl_lface = max(kfl_lface,1_ip)
  end if
  !
  ! List of element boundary faces: LELBF
  !
  if( kfl_modul(ID_PARTIS) /= 0 ) then
     kfl_lelbf = max(kfl_lelbf,1_ip)
  end if
  !
  ! Extended node-element: graph: LELPO_2, PELPO_2
  !
  if( kfl_modul(ID_IMMBOU) /= 0 ) then
     kfl_lelp2 = max(kfl_lelp2,1_ip)
  end if
  if( kfl_lface /= 0 ) then
     kfl_lelp2 = max(kfl_lelp2,1_ip)     
  end if
  if( kfl_hangi == 1 ) then
     kfl_lelp2 = max(1_ip,kfl_lelp2)
  end if
  !
  ! Extended element-element: graph: LELEL_2, PELEL_2
  !
  if( kfl_modul(ID_PARTIS) /= 0 ) then
     kfl_lele2 = max(kfl_lele2,1_ip)
  end if
  if( necnt > 0 ) then
     kfl_lele2 = max(kfl_lele2,1_ip)
  end if
  if( kfl_hangi == 1 ) then
     kfl_lele2 = max(1_ip,kfl_lele2)
  end if
  !
  ! Element bin
  !
  !if( kfl_modul(ID_PARTIS) /= 0 ) then
  !   kfl_element_bin = max(kfl_element_bin,1_ip)
  !end if
  !
  ! Symmetric graph: R_SYM, C_SYM
  !
  do imodu = 1,mmodu
     if( kfl_modul(imodu) /= 0 ) then
        if( associated(momod(imodu) % solve) ) then
           do ivari = 1,size( momod(imodu) % solve,KIND=ip)
              if( momod(imodu) % solve(ivari) % kfl_symme == 1 ) then
                 kfl_symgr = max(kfl_symgr,1_ip)
              end if
           end do
        end if
     end if
  end do
  !
  ! NZDOM_*, R_DOM, C_DOM_*. Schur graphs for * = Aii, Aib, Abi, Abb
  !
  do imodu = 1,mmodu
     if( kfl_modul(imodu) /= 0 ) then
        if( associated(momod(imodu) % solve) ) then
           do ivari = 1,size( momod(imodu) % solve,KIND=ip )
              if( momod(imodu) % solve(ivari) % kfl_schur == 1 ) then
                 kfl_schur = max(kfl_schur,1_ip)
              end if
              if(      momod(imodu) % solve(ivari) % kfl_preco == 10 ) then
                 kfl_schur = max(kfl_schur,3_ip)
              else if( momod(imodu) % solve(ivari) % kfl_preco == 11 ) then
                 kfl_schur = max(kfl_schur,2_ip)
              else if( momod(imodu) % solve(ivari) % kfl_preco == 12 ) then
                 kfl_schur = max(kfl_schur,2_ip)
              else if( momod(imodu) % solve(ivari) % kfl_preco == 13 ) then
                 kfl_aiipr = max(kfl_aiipr,1_ip)
              end if
           end do
        end if
     end if
  end do

end subroutine reqarr
