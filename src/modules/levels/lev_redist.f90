!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_redist
  !-----------------------------------------------------------------------
  !****f* Levels/lev_redist
  ! NAME 
  !    lev_redist
  ! DESCRIPTION
  !    Compute the level set function redistanciation geometrically
  ! USES
  !    
  ! USED BY
  !    lev_iniunk
  !    lev_updunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use def_solver
  use mod_ker_updpro,     only : ker_updpro
  use mod_messages,       only : livinf
  use mod_communications, only : PAR_BARRIER
  use mod_communications, only : PAR_MIN
  use mod_run_config,     only : run_config

  implicit none
  integer(ip)             :: idime,icomp
  integer(ip)             :: ipoin
  integer(ip)             :: compp(ndime)
  integer(ip), save       :: ipass=0
  integer(ip)             :: prfin

  real(rp)                :: dista, time1, time
  real(rp)                :: d,dx,dy,ratio,px,py,pz
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,xa,ya
  real(rp)                :: timr0,timr1

  ! Initial verification of the redistanciation necessity

  call livinf(59_ip,'GEOMETRICAL REDISTANTIATION',0_ip)
  if( INOTSLAVE ) call cputim(time)

  !
  ! Initializes the distance
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        dista_lev(ipoin) = 1.0e+06
     end do
     prfin = 0_ip
  end if
  !
  ! COORD_LEV + LNODB_LEV: First phase = discrete interface construction
  !
  if(run_config%timing) call PAR_BARRIER()
  call cputim(timr0)
  call lev_conint()
  if(run_config%timing) call PAR_BARRIER()
  call cputim(timr1)
  if( INOTSLAVE .and. run_config%timing) then
     print*,'cpu_time: lev_conint',timr1-timr0
  end if
  
  if( kfl_geodi_lev == 1) then
  !if( 1 == 1 ) then
     !
     ! New geometrical distance
     !
     call lev_cristo()

  else
     !
     ! Last step : distance is computed from discrete interface data
     !
     if( INOTMASTER ) then

        if(ndime==2) then
           !
           ! Compute the distance 2d
           !
           do ipoin=1,npoin

              px=coord(1,ipoin)
              py=coord(2,ipoin)

              do icomp=1,nboun_lev

                 do idime=1,ndime
                    compp(idime)=lnodb_lev(idime,icomp)
                 end do

                 x1=coord_lev(1,compp(1))
                 y1=coord_lev(2,compp(1))
                 x2=coord_lev(1,compp(2))
                 y2=coord_lev(2,compp(2))
                 dx=x2-x1
                 dy=y2-y1

                 ratio = ((px-x1)*dx+(py-y1)*dy)/(dx*dx+dy*dy)

                 if( ratio < 0.0_rp ) then 
                    d = sqrt((px-x1)*(px-x1)+(py-y1)*(py-y1))
                 else if(ratio>1) then 
                    d = sqrt((px-x2)*(px-x2)+(py-y2)*(py-y2))
                 else
                    xa = (1.0_rp-ratio)*x1+ratio*x2
                    ya = (1.0_rp-ratio)*y1+ratio*y2
                    d  = sqrt((px-xa)*(px-xa)+(py-ya)*(py-ya))
                 endif

                 dista            = dista_lev(ipoin)
                 dista            = min(dista,d)
                 dista_lev(ipoin) = dista

              end do

           end do

        else if(ndime==3) then

           do ipoin=1,npoin

              px=coord(1,ipoin)
              py=coord(2,ipoin)
              pz=coord(3,ipoin)


              do icomp=1,nboun_lev

                 do idime=1,ndime
                    compp(idime)=lnodb_lev(idime,icomp)
                 end do

                 x1=coord_lev(1,compp(1))
                 y1=coord_lev(2,compp(1))
                 z1=coord_lev(3,compp(1))

                 x2=coord_lev(1,compp(2))
                 y2=coord_lev(2,compp(2))
                 z2=coord_lev(3,compp(2))

                 x3=coord_lev(1,compp(3))
                 y3=coord_lev(2,compp(3))
                 z3=coord_lev(3,compp(3))

                 call dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

                 dista=dista_lev(ipoin)
                 dista=min(dista,d)
                 dista_lev(ipoin)=dista


              end do


           end do
        endif
     end if
  end if
  !
  ! Deallocate memory
  !
  call lev_memall(8_ip)
  !
  !    Updating of redistanced level set function
  !
  if( INOTMASTER ) then
     if (kfl_zonal_lev==0)then
        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)<1) then
              if(fleve(ipoin,1)>=0.0_rp) then
                 fleve(ipoin,1)= dista_lev(ipoin)
              else
                 fleve(ipoin,1)=-dista_lev(ipoin)
              endif
           end if
        end do
     else
        call lev_zonalr()
     end if
  end if

  if( INOTMASTER ) prfin = 1_ip

  call PAR_MIN(prfin)

  if( INOTSLAVE ) then
     if(prfin==1_ip) then
        if(ipass==0) then
           if(kfl_rstar/=2) then
              ipass=1
              write(lun_cored_lev,100)
           end if
        end if
        call cputim(time1)
        write(lun_cored_lev,101) ittim,time1-time
        flush(lun_cored_lev)
     endif
  endif

100 format('# --|  Redistanciation duration ' ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time step         2. Redistanciation duration  ')
101 format(4x,i9,2x,e12.6)

end subroutine lev_redist

subroutine lev_cristo()
  use def_elmtyp
  use def_master
  use def_domain
  use mod_kdtree
  use def_levels
  use mod_memory
  use mod_communications, only : PAR_BARRIER
  use mod_run_config,     only : run_config
  implicit none
  real(rp)             :: xcoor(3)
!  real(rp)             :: bobox(3,2)
!  real(rp),    pointer :: fabox(:,:,:) => null()
!  real(rp),    pointer :: sabox(:,:,:) => null()
!  integer(ip), pointer :: blink(:)     => null()
!  integer(ip), pointer :: stru2(:)     => null()
!  type(netyp), pointer :: lnele(:)     => null()
!  real(rp),    pointer :: ldist(:)     => null()
  integer(ip), pointer :: ltypb_lev(:) => null()
  type(typ_kdtree)     :: kdtree_lev
  integer(ip)          :: mnodb_lev,iboun
!  integer(ip)          :: mnodb_lev,iboun,kboun
  real(rp)             :: chkdi
!  real(rp)             :: norma(ndime)
  real(rp)             :: proje(ndime)
  real(rp)             :: time0,time1,time2
  integer(ip)          :: ipoin
  !
  ! Initialize variables
  !
  if( INOTMASTER ) then
     chkdi     = 1.0e10_rp
     mnodb_lev = ndime
     call memory_alloca(memor_dom,'LTYPB_LEV','lev_cristo',ltypb_lev,nboun_lev)
     if( ndime == 2 ) then
        do iboun = 1,nboun_lev
           ltypb_lev(iboun) = BAR02
        end do
     else
        do iboun = 1,nboun_lev
           ltypb_lev(iboun) = TRI03
        end do
     end if
  end if
  ! 
  ! Compute kd-tree structure
  !
  if(run_config%timing) call PAR_BARRIER()
  call cputim(time0)
  if( INOTMASTER ) then
     call kdtree_construct(&
          nboun_lev,npoin_lev,lnodb_lev,ltypb_lev,coord_lev,kdtree_lev)
     !call kdtree(& 
     !1_ip,mnodb_lev,npoin_lev,nboun_lev,coord_lev,&
     !lnodb_lev,ltypb_lev,fabox,bobox,sabox,blink,&
     !stru2,ldist,lnele)
  end if
  if(run_config%timing) call PAR_BARRIER()
  call cputim(time1)
  !
  ! Compute shortest distance
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        xcoor(    1) = coord(    1,ipoin)
        xcoor(    2) = coord(    2,ipoin)
        xcoor(ndime) = coord(ndime,ipoin)
        call kdtree_nearest_boundary(xcoor,kdtree_lev,iboun,dista_lev(ipoin),proje)
        !call dpopar(&
        !     ndime,xcoor,npoin_lev,mnodb_lev,nboun_lev,chkdi,ltypb_lev,&
        !     lnodb_lev,coord_lev,dista_lev(ipoin),norma,proje,kboun,&
        !     fabox,sabox,blink,stru2,ldist,lnele)

        dista_lev(ipoin) = abs(dista_lev(ipoin))
     end do
  end if
  if(run_config%timing) call PAR_BARRIER()
  call cputim(time2)
  !
  ! Deallocate memory
  !
  if( INOTMASTER ) then
     call kdtree_deallocate(kdtree_lev)
     !call kdtree(&
     !     2_ip,mnodb_lev,npoin_lev,nboun_lev,coord_lev,&
     !     lnodb_lev,ltypb_lev,fabox,bobox,sabox,blink,&
     !     stru2,ldist,lnele)
     call memory_deallo(memor_dom,'LTYPB_LEV','lev_cristo',ltypb_lev)
  end if

  if( INOTSLAVE .and.run_config%timing) then
     print*,'KD-TREE CARAS =',nboun_lev
     print*,'KD-TREE TIME 1=',time1-time0
     print*,'KD-TREE TIME 2=',time2-time1
  end if

end subroutine lev_cristo
