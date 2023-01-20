!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_nulptr()

  use def_magnet

  implicit none

  nullify(edglen_mag)
  nullify(elevol_mag)
  nullify(elesig_mag)
  nullify(edgsig_mag)
  nullify(locboucen_mag)
  nullify(globoucen_mag)
  nullify(locboutan_mag)
  nullify(globoutan_mag)
  
  nullify(eleedg_mag)
  nullify(diredg_mag)
  nullify(edgnod_mag)
  nullify(bouedg_mag)

  nullify(edgdir_mag)
  nullify(edgbou_mag)

  nullify(He_mag)
  nullify(Ha_mag)
  nullify(Hp_mag)
  nullify(dH_mag)
  nullify(bhsid)
  nullify(fhsid)
  nullify(hhsid)

  nullify(Hc_mag)
  nullify(Hn_mag)
  nullify(Jcz_mag)
  nullify(Jnz_mag)
  nullify(Bc_mag)
  nullify(Bn_mag)
  nullify(Jc_mag)
  nullify(Jn_mag)
  nullify(Fn_mag)
  nullify(Fc_mag)
  nullify(Gc_mag)
  nullify(Gn_mag)

  nullify(Hgp_mag)
  nullify(Jgp_mag)
  nullify(Bgp_mag)
  nullify(Fgp_mag)
  nullify(Ggp_mag)

      
  nullify(Hpav_mag)
  nullify(Hxtl_mag)
  nullify(Hsfl_mag)
  nullify(Hsfg_mag)
  nullify(Hbs_mag)
  
  nullify(proc_nbedg)

  nullify(kfl_fixno_mag)
  nullify(kfl_fixbo_mag)
  nullify(bvess_mag)
  nullify(bvnat_mag)
  nullify(edgind_mag)
  nullify(nodind_mag)

  nullify(tecod_mag)
  nullify(tbcod_mag)

  nullify(bdfode_mag % y1)
  nullify(bdfode_mag % y2)
  nullify(bdfode_mag % y3)
  nullify(bdfode_mag % y4)
  nullify(bdfode_mag % y5)
  nullify(bdfode_mag % y6)
  nullify(bdfode_mag % ypi)
  nullify(bdfode_mag % ype)

end subroutine mag_nulptr
