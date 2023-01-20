!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_endste
  !-----------------------------------------------------------------------
  !****f* Levels/lev_endste
  ! NAME 
  !    lev_endste
  ! DESCRIPTION
  !    This routine ends a time step of the level set convection.
  ! USES
  !    lev_cvgunk
  !    lev_updunk
  !    lev_output
  ! USED BY
  !    Levels
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain

  implicit none
  integer(ip)  :: i
  real(rp)     :: xintf,phi1,phi2,x1,x2,a,g,xadim,tadim
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if( kfl_stead_lev == 0 ) then
     call lev_cvgunk(three)
     call lev_updunk( five)
  end if
  !
  ! Postprocess Gauges
  !
  if(npp_gauge_lev==1) then

     do i = 1, npp_nbgau_lev
        valga_lev(i) = 0.0
        findg_lev(i) = 0
     enddo

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
        if(findg_lev(i)/=0) valga_lev(i) = valga_lev(i) /findg_lev(i)
     enddo

     if( INOTSLAVE ) then
        write(lun_gauge_lev,1) cutim, (valga_lev(i),i=1,npp_nbgau_lev)
        flush(lun_gauge_lev)
     endif

  endif

  !
  ! Write restart file
  !
  call lev_restar(two)

  if(kfl_paral==-1_ip) then
     if(kfl_inlev_lev==8) then
        a=0.25_rp
        g=9.81_rp
        do i=1,ncapt_lev-1
           phi1=fleve(capin_lev(i),1)
           phi2=fleve(capin_lev(i+1),1)
           if(phi1*phi2<0.0_rp) then
              x1=coord(1,capin_lev(i))
              x2=coord(1,capin_lev(i+1))
              xintf=x1+(x2-x1)*phi1/(phi1-phi2)
              xadim=xintf/a
              tadim=cutim*sqrt(2*g/a)
              write(lun_capte_lev, *) tadim,' ',xadim
           endif
        enddo
     endif
  endif
  !
  ! If not steady, go on
  !
  if(kfl_stead_lev==0) kfl_gotim = 1

  call lev_calvol


1 format(e12.6,10(1x,e12.6))

end subroutine lev_endste
