!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_redisc()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_redisc
  ! NAME 
  !    lev_redisc
  ! DESCRIPTION
  !    Compute the level set function redistanciation geometrically only on cut nodes
  ! USES
  !    
  ! USED BY
  !    lev_redieq
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_levels
  use def_domain
  use def_solver
  use mod_memchk
  use mod_communications
  use mod_messages,       only : livinf
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  integer(ip)             :: idime,icomp,ii,jj
  integer(ip)             :: ipoin
  integer(ip)             :: nlevp,nlevn,ielem,inode,pnode,pelty
  integer(ip)             :: compp(ndime)
  integer(ip), save       :: ipass=0
  integer(ip)             :: prfin
  integer(4)              :: istat

  real(rp)                :: dista, time1, time
  real(rp)                :: d,dx,dy,ratio,px,py,pz
  real(rp)                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,xa,ya
  real(rp)                :: ellev(mnode)

  ! Initial verification of the redistanciation necessity

  call livinf(59_ip,'GEOMETRICAL REDIST. IN CUT NODES',0_ip)
  if( INOTSLAVE ) call cputim(time)
  !
  ! Identification of cut nodes
  !
  if( INOTMASTER ) then
     call memgen(1_ip,npoin,0_ip)  ! allocate gisca
     do ipoin = 1,npoin
        gisca(ipoin) = 0
     end do
     do ielem=1,nelem
        nlevp=0
        nlevn=0
        pelty=ltype(ielem)
        if( lelch(ielem) == ELFEM ) then
           pnode=nnode(pelty)
           ! Gather operations
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode)=fleve(ipoin,1)
              if(ellev(inode)>0.0_rp) nlevp=nlevp+1
              if(ellev(inode)<0.0_rp) nlevn=nlevn+1
           end do
           if (min(nlevp,nlevn)>0)then
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 gisca(ipoin) = gisca(ipoin) + 1
              end do
           end if
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
  end if
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
  call lev_conint()
  !
  ! Last step : distance is computed from discrete interface data
  !
  if( INOTMASTER ) then

     if(ndime==2) then
        !
        ! Compute the distance 2d
        !
        do ipoin=1,npoin
           if(gisca(ipoin)>0) then

              px=coord(1,ipoin)
              py=coord(2,ipoin)

              do ii = pefpo(ipoin),pefpo(ipoin+1)-1
                 ielem = lefpo(ii)                             ! elements related to point ipoin
                 do jj = psbel_lev(ielem),psbel_lev(ielem+1)-1
                    icomp = lsbel_lev (jj)                     ! surfaces belonging to element

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
              !              write(kfl_paral+700,*)ipoin,dista   !hhh
           end if
        end do

     else if(ndime==3) then

        do ipoin=1,npoin
           if(gisca(ipoin)>0) then

              px=coord(1,ipoin)
              py=coord(2,ipoin)
              pz=coord(3,ipoin)


              do ii = pefpo(ipoin),pefpo(ipoin+1)-1
                 ielem = lefpo(ii)                             ! elements related to point ipoin
                 do jj = psbel_lev(ielem),psbel_lev(ielem+1)-1
                    icomp = lsbel_lev (jj)                     ! surfaces belonging to element

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

           end if
        end do
     endif
     !
     ! Deallocate memory
     !
     call lev_memall(8_ip)
     !
     !    Updating of redistanced level set function
     !
     if (kfl_zonal_lev==0)then
        do ipoin=1,npoin
           if(gisca(ipoin)>0)then
              if(kfl_fixno_lev(1,ipoin)<1) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              end if
           end if
        end do
     elseif ( (kfl_zonal_lev>=3).or.(kfl_zonal_lev<=6) )then
        call lev_zonalr
     else
        call runend('ZONAL 1 and 2 not ready in lev_redisc')
     end if

  endif

  if( INOTMASTER ) prfin = 1_ip

  if( IPARALL ) then
     call PAR_MIN(prfin)
  end if

  if(kfl_paral<=0) then
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

  if( INOTMASTER ) then
     do ipoin=1,npoin
        icupt_lev(ipoin)=gisca(ipoin)
     end do
     call memgen(3,npoin,0)  ! deallocate gisca
  end if

  if( INOTMASTER ) then

     !
     ! Deallocate LSBEL_LEV that has been allocated in lev_conint
     !
     call memchk(two,istat,mem_modul(1:2,modul),'LSBEL_LEV','lev_redisc',lsbel_lev)
     deallocate(lsbel_lev, stat=istat)
     if(istat/=0) call memerr(two,'LSBEL_LEV','lev_redisc',0_ip)

  end if


100 format('# --|  Redistanciation duration ' ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time step         2. Redistanciation duration  ')
101 format(4x,i9,2x,e12.6)

end subroutine lev_redisc

