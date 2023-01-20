  !$OMP  PARALLEL DO SCHEDULE (GUIDED) & 
  !$OMP  DEFAULT (NONE)  &
  !$OMP  PRIVATE (       &
  !$OMP           elauu, &
  !$OMP           elaup, &
  !$OMP           elapp, &
  !$OMP           elapu, &
  !$OMP           elrbu, &
  !$OMP           elrbp, &
  !$OMP           elmap, &
  !$OMP           elvel, &
  !$OMP           elpre, &
  !$OMP           elfle, &
  !$OMP           elcod, &
  !$OMP           ielem, &
  !$OMP           igaus, &
  !$OMP           pnode, &
  !$OMP           pgaus, &
  !$OMP           pevat, &
  !$OMP           pelty, &
  !$OMP           plapl, &
  !$OMP           porde, &
  !$OMP           pmate, &
  !$OMP           ptopo, &
  !$OMP           gpcar, &
  !$OMP           gphes, &
  !$OMP           gplap, &
  !$OMP           gpvol, &
  !$OMP           gpvis, &
  !$OMP           gpgvi, &
  !$OMP           gpmut, &
  !$OMP           grvis, &
  !$OMP           gppor, &
  !$OMP           gpden, &
  !$OMP           gpfle, &
  !$OMP           gpst1, &
  !$OMP           gpst2, &
  !$OMP           gpsp1, &
  !$OMP           gpsp2, &
  !$OMP           gptt1, &
  !$OMP           gptt2, &
  !$OMP           gpadv, &
  !$OMP           gprhs, &
  !$OMP           gpsgs, &
  !$OMP           gpsgt, &
  !$OMP           gppre, &
  !$OMP           gpvel, &
  !$OMP           gpgve, &
  !$OMP           gpgpr, &
  !$OMP           gptem, &
  !$OMP           gpsgi, &
  !$OMP           tragl, &
  !$OMP           chave, &
  !$OMP           chale, &
  !$OMP           hleng, &
  !$OMP           dummr, &
  !$OMP           rmomu, &
  !$OMP           rcont, &
  !$OMP           wgrgr, &
  !$OMP           wgrvi, &
  !$OMP           p1vec, &
  !$OMP           p2vec, &
  !$OMP           p2sca) &
  !$OMP   SHARED ( kfl_ellen_nsi, elmar, hnatu, kfl_advec_nsi, nunkn_nsi,  &
  !$OMP            kfl_savco_nsi, npoin, resgs_nsi, itask,                 &
  !$OMP            lorde, solve_nsi, ivari_nsi, llapl, ltype, nnode,       &
  !$OMP            ngaus, kfl_sgste_nsi, tesgs, lnods, kfl_sgsti_nsi,      &
  !$OMP            kfl_sgsco_nsi, ndime, vesgs, rhsid, amatr, lmate,       &
  !$OMP            nmatf, kfl_prtre_nsi, nmate, kfl_colev_nsi, lapla_nsi,  &
  !$OMP            kfl_cotur_nsi, poauu_nsi, poaup_nsi, poapu_nsi,         &
  !$OMP            poapp_nsi, ndbgs_nsi, kfl_ellsh_nsi, kfl_predi_nsi,     &
  !$OMP            nelem, ltopo ) 
