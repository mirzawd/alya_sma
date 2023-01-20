!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_zgrad_lookup()
   !------------------------------------------------------------------------
   ! lookup and integration of heat release rate
   !------------------------------------------------------------------------    
   use def_chemic,             only : kfl_zg_fw_chm, zgrad_gp
   use def_master,             only : INOTEMPTY
   use def_kermod,             only : lookup_fw

   external                  :: chm_post_gp_lookup

   if ( INOTEMPTY ) then
      !
      ! Lookup from HRR table
      !
      call chm_post_gp_lookup(zgrad_gp, lookup_fw( kfl_zg_fw_chm ))

   endif

end subroutine chm_zgrad_lookup
