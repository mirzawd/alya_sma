!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_iniunk()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_iniunk
  ! NAME 
  !    lev_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the wave amplitude
  ! USED BY
  !    lev_begste
  !***
  !-----------------------------------------------------------------------
  use      def_parame, only       :  ip,rp
  use      def_master, only       :  fleve,kfl_rstar,IPARALL,NPARR,PARRE,NPARI,PARIN,INOTSLAVE,INOTMASTER
  use      def_domain
  use      def_levels
  use      mod_memchk
  use def_kermod,     only : kfl_cutel
  use mod_ker_elmcut, only : ker_elmcut_free_surface
  implicit none
  integer(ip)              :: ipoin,i
  real(rp)                 :: t
  !
  ! Load initial conditions for the level set function
  !
  if ( kfl_rstar == 0 ) then
     ! the part that was here is now in lev_inibcs so that ker_updpro(ITASK_INIUNK) has fleve
  else     
     !
     ! Read restart file
     !
     call lev_restar(1_ip)    
  end if

  if( INOTMASTER ) then
     if(kfl_conbc_lev==0)then
        do ipoin = 1,npoin
           if(kfl_funno_lev(ipoin)==0)fleve(ipoin,1) = fleve(ipoin,3)
        end do
     else
        do ipoin = 1,npoin
           fleve(ipoin,1) = fleve(ipoin,3)
        end do
     end if
  end if
  !
  ! Optional initial redistanciation
  !
  if( inred_lev == 1 ) then
     if(       tyred_lev == 1 ) then
        call lev_redist()
     else if( tyred_lev == -1 ) then
        call lev_redist_generalized_distance()
     else if( tyred_lev == -2 ) then
        call lev_redist_geometrical_distance()
     else if( tyred_lev >   1 ) then
        call lev_redieq()
     end if
  end if
  !
  ! Update unknowns
  !
  call lev_updunk(11_ip)
  !
  ! Postprocess Gauges
  !
  if(npp_gauge_lev==1) then

     do i = 1, npp_nbgau_lev
        valga_lev(i) = 0.0_rp
        findg_lev(i) = 0
     end do

     call lev_outgau() 

     if( IPARALL ) then

        nparr =  ngaug_lev
        parre => valga_lev
        call par_operat(3_ip)

        npari =  ngaug_lev
        parin => findg_lev
        call par_operat(3_ip)

     end if

     do i = 1, npp_nbgau_lev
        if(findg_lev(i)/=0) valga_lev(i) = valga_lev(i) / real(findg_lev(i),rp)
     enddo

     t=0.0_rp

     if( INOTSLAVE ) then
        write(lun_gauge_lev,1) t, (valga_lev(i),i=1,npp_nbgau_lev)
        flush(lun_gauge_lev)
     endif

  endif
  !
  ! Cut elements
  !
  if( kfl_cutel == 1 ) call ker_elmcut_free_surface()

  ! 
  ! Output initial liquid phase volume
  !

  if(kfl_corvo_lev>0_ip) then
     call lev_calvol()
     volrf_lev=volit_lev
     if( INOTSLAVE ) then   
        write(lun_volum_lev, *) 0_rp,' ',volit_lev  
        !print*,' vol ref ',volrf_lev
     endif
  endif

1 format(e13.6,10(1x,e12.6))

end subroutine lev_iniunk
