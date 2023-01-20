!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_wallav(itask)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_wallav
  ! NAME 
  !    nsi_wallav
  ! DESCRIPTION
  !    Calculate average velocity at the boundaries:
  !    U(n+1) = e*u(n+1) + (1-e)*U(n)
  !    where:
  !          e = dt/T is the weighting function
  !          T is the averaging period
  !          U(n) is the averaged velocity at time-step n
  !          u(n) is the instantaneous velocity at time-step n
  !
  !    ITASK = 1 ... Calculates time-averaged velocity at the boundary gauss points  (velav_ker)
  !            2 ... Calculates time-averaged velocity at the nodes (for postprocessing) (avupo_nsi)
  !            3 ... Calculates initial value of time-averaged velocity (for gauss points)
  !
  !            Beware I belive the nodal value should be also initialized
  !
  ! USES
  ! USED BY
  !    nsi_averag
  !    nsi_outtan
  ! OUTPUTS
  !    velav_ker(ndime, ngaub, nboun)
  !    avupo_nsi(ndime, npoin)
  !***
  !-----------------------------------------------------------------------
  use def_master,             only : ip,rp,dtime,solve,veloc,inotmaster,kfl_rstar
!  use def_master,             only : ittim
  use def_domain,             only : ltypb,lelbo,ltype,ngaus,lmate,nelem,nboun,nnode,lnods,coord,lnodb,kfl_codbo
  use def_domain,             only : ndimb,elmar,ndime,mnode,mnodb,npoin,lpoty
  use def_nastin,             only : kfl_fixbo_nsi,ivari_nsi,kfl_wlare_nsi
  use def_kermod,             only : kfl_waexl_ker, velel_ker, lexlo_ker, velav_ker, tpeav_ker
  use def_kermod,             only : kfl_noslw_ker,kfl_nswpo_ker,avupo_ker,lnsw_exch,kfl_fixbo_nsw_ker,kfl_boexc_ker
  use mod_interpolation,      only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_calc_avta1_nsw_ker, only : calc_avta1_nsw_ker
  use mod_nsi_nsw_averages,   only : nsi_nsw_averages
  use mod_wall_exchange,      only : ker_waexlo_getval
  use def_kintyp,             only : lg
  use mod_bouder

  implicit none

  integer(ip),intent(in)    :: itask

  real(rp)                  :: eucta
  real(rp)                  :: baloc(ndime,ndime)
  real(rp)                  :: elvel(ndime,mnode)
  real(rp)                  :: elcod(ndime,mnode)
  real(rp)                  :: bovel(ndime,mnodb)
  real(rp)                  :: bocod(ndime,mnodb)
  real(rp)                  :: vewal(ndime), dt
  real(rp)                  :: tvelo(ndime)
  real(rp)                  :: avwei   ! e=dt/T, weighting function for averaging
!  real(rp)                  :: ux,uy

  integer(ip)               :: ipoin,igaub,iboun,ielem
  integer(ip)               :: idime,ibopo,inode,inodb
  integer(ip)               :: kdime,pblty,pnodb,pelty
  integer(ip)               :: pnode,pgaub,pevat,pmate
!  integer(ip)               :: nav   ! nav=N first time steps where e=1/ittim
  logical(lg)               :: lauxi


  !  if ( kfl_timco /= 0 ) then
  !     call runend('TIME AVERAGING FOR THE WALL LAW ONLY READY FOR PRESCRIBED dt')
  !  else
  ! nav = int(tpeav_ker/dtime) ! nav=N first steps where I perform standard averaging  ! this may not be teh best option fon non-constant time step - REVISE
  !  end if

  ! if ( itask < 3_ip ) then
  !    if ( ittim .le. nav ) then
  !       if (ittim == 0_ip)  then       ! to avoid divide by 0  -- check if it is the best option
  !          avwei = 1.0_rp              ! For the first N time-steps e=1/ittim
  !       else
  !          avwei = 1.0_rp/real(ittim,rp) ! For the first N time-steps e=1/ittim
  !       end if
  !    else
  !       avwei = dtime/tpeav_ker    ! e=dt/T  where T is the time-period for averaging
  !    end if
  !end if

  avwei = 1.0_rp-(dtime/(dtime+tpeav_ker))
  dt = dtime/(dtime+tpeav_ker)

  if ( itask == 1_ip ) then

     !-----------------------------------------------------------------------------
     !
     ! Calculate average velocity at the boundary gauss points for nsi_bouwal
     ! in the case of no slip wall law this is used in ker_nsw_visc (comes in avelavv)
     !
     !-----------------------------------------------------------------------------
     
     if (kfl_noslw_ker == 1_ip ) then
        elements0: do ielem = 1,nelem
           if ( lnsw_exch(ielem)%nbogp /= 0_ip ) then
              lnsw_exch(ielem)%vel_aux(1:ndime) = 0.0_rp
           end if
        end do elements0
     end if

     boundaries: do iboun = 1,nboun
        lauxi = .false.        ! to allocate kfl_fixbo_nsw_ker only if kfl_noslw_ker /= 0
        if( kfl_noslw_ker /= 0 ) then
           if( kfl_fixbo_nsw_ker(iboun) == 1_ip ) lauxi = .true. ! no slip wall law 
        end if

        if(  kfl_fixbo_nsi(iboun) ==  3 .or. &          ! Wall law
             & kfl_fixbo_nsi(iboun) == 13 .or. &        ! Wall law + open pressure
             & kfl_fixbo_nsi(iboun) == 18 .or. &        ! u.n in weak form
             & lauxi ) then                             ! no slip wall law

           !
           ! Element properties and dimensions
           !
           pblty = ltypb(iboun) 
           pnodb = nnode(pblty)
           ielem = lelbo(iboun)
           pelty = ltype(ielem)

           if( pelty > 0 ) then

              pnode = nnode(pelty)
              pgaub = ngaus(pblty) 
              pevat = solve(ivari_nsi) % ndofn * pnode
              pmate = lmate(ielem)

              if( pmate /= -1 ) then
                 !
                 ! Gather operations: ELVEL, ELCOD, BOVEL
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
                    elcod(1:ndime,inode) = coord(1:ndime,ipoin)    
                 end do

                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    bocod(1:ndime,inodb) = coord(1:ndime,ipoin)    
                    bovel(1:ndime,inodb) = veloc(1:ndime,ipoin,1)
                 end do

                 gauss_points: do igaub = 1,pgaub

                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                         bocod,baloc,eucta)                                   ! and Jacobian
                    call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                    if ( kfl_waexl_ker == 0_ip .and. kfl_noslw_ker == 0_ip ) then  ! normal behaviour

                       do idime = 1,ndime                                     ! Velocity
                          vewal(idime) = 0.0_rp
                       end do
                       do inodb = 1,pnodb
                          do idime = 1,ndime
                             vewal(idime) = vewal(idime) &
                                  + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                          end do
                       end do

                    else  ! using velocity from exchange location -- also used for no slip wall law

                       if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
                          do idime = 1,ndime
                             vewal(idime) = velel_ker(idime,lexlo_ker(igaub,iboun))   ! Velocity - this comes as an input to this subroutine
                          end do
                       end if

                    end if

                    do idime = 1,ndime
                       tvelo(idime) = vewal(idime)
                       do kdime = 1,ndime
                          tvelo(idime) = tvelo(idime)   &
                               - baloc(idime,ndime) &
                               * baloc(kdime,ndime) * vewal(kdime)
                       end do
                    end do
                    ! averaged tangent velocity
                    ! first option for no slip wall was obtain the average (over all associated boundary gauss points) and then the average in time   
                    !                    if ( kfl_fixbo_nsw_ker(iboun) == 1_ip) then ! no slip wall law - an elemental value is obtained averaging over all boundary gauss points relaed to the element
                    !                       lnsw_exch(ielem)%vel_aux(1:ndime) = lnsw_exch(ielem)%vel_aux(1:ndime) + tvelo(1:ndime)                       
                    !                    else    ! rest of cases normal behaviour
                    !                       velav_ker(1:ndime,igaub,iboun) = avwei*tvelo(1:ndime)+(1.0_rp-avwei)*velav_ker(1:ndime,igaub,iboun)
                    !                    end if
                    !Second option first average in time and then for no slip wall law average in over all associated boundary gauss points
                    ! they should be the same but the second option uses  velav_ker as the exchange location case thus it si preferable
                    ! This second choice also alters the call to call nsi_nsw_averaging   before I was using nsi_nsw_averaging(0_ip   since teh time averaging had not been performed but now I change to
                    ! nsi_nsw_averaging(1_ip  - Actaully this is better because it eliminates the need for this option and I can therefore eliminate it
                    if (kfl_noslw_ker == 0_ip) then
                       velav_ker(1:ndime,igaub,iboun) = dt*tvelo(1:ndime)+avwei*velav_ker(1:ndime,igaub,iboun)
                       !velav_ker(1:ndime,igaub,iboun) = avwei*tvelo(1:ndime)+(1.0_rp-avwei)*velav_ker(1:ndime,igaub,iboun)
                    else 
                       velav_ker(1:ndime,igaub,iboun) = tvelo(1:ndime)
                    endif

                 end do gauss_points
              end if
           end if
        end if
     end do boundaries
     call nsi_nsw_averaging(0_ip,avwei,dt)     ! obtains lnsw_exch(ielem)%velav

     call calc_avta1_nsw_ker(avwei)
     if( kfl_noslw_ker /= 0 ) call nsi_nsw_averages(avwei)

  else if ( itask == 2_ip ) then

     !-----------------------------------------------------------------------------
     !
     ! Calculate average velocity at the boundary nodes for post-processing
     !
     !-----------------------------------------------------------------------------

     if (kfl_noslw_ker == 0_ip) then
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if (ibopo > 0_ip) then
              avupo_ker(1:ndime,ipoin) = dt*veloc(1:ndime,ipoin,1)+avwei*avupo_ker(1:ndime,ipoin)
              !avupo_ker(1:ndime,ipoin) = avwei*veloc(1:ndime,ipoin,1)+(1.0_rp-avwei)*avupo_ker(1:ndime,ipoin)
           end if
        end do
     else
        do ipoin=1,npoin
           if (kfl_nswpo_ker(ipoin) > 0_ip) then
              avupo_ker(1:ndime,ipoin) = dt*veloc(1:ndime,ipoin,1)+avwei*avupo_ker(1:ndime,ipoin)
              !avupo_ker(1:ndime,ipoin) = avwei*veloc(1:ndime,ipoin,1)+(1.0_rp-avwei)*avupo_ker(1:ndime,ipoin)
           end if
        end do
     end if

  else if ( itask == 3_ip ) then

     !-----------------------------------------------------------------------------
     !
     ! Calculate initial value of average velocity
     !
     !-----------------------------------------------------------------------------
     if ( kfl_waexl_ker == 1_ip ) then
        call ker_waexlo_getval(veloc,velel_ker) 
     end if

     if ( INOTMASTER )then

        if ( kfl_wlare_nsi == 0_ip )then ! No restart file found

           if (kfl_noslw_ker == 1_ip ) then
              elements03: do ielem = 1,nelem
                 if ( lnsw_exch(ielem)%nbogp /= 0_ip ) then
                    lnsw_exch(ielem)%vel_aux(1:ndime) = 0.0_rp
                 end if
              end do elements03
           end if

           boundaries3: do iboun = 1,nboun
              lauxi = .false.        ! to allocate kfl_fixbo_nsw_ker only if kfl_noslw_ker /= 0
              if( kfl_noslw_ker /= 0 ) then
                 if( kfl_fixbo_nsw_ker(iboun) == 1_ip ) lauxi = .true. ! no slip wall law 
              end if
              if(    kfl_fixbo_nsi(iboun) ==  3 .or. &        ! Wall law
                   & kfl_fixbo_nsi(iboun) == 13 .or. &        ! Wall law + open pressure
                   & kfl_fixbo_nsi(iboun) == 18 .or. &        ! u.n in weak form
                   & lauxi) then   ! no slip wall law
                 !
                 ! Element properties and dimensions
                 !
                 pblty = ltypb(iboun) 
                 pnodb = nnode(pblty)
                 ielem = lelbo(iboun)
                 pelty = ltype(ielem)

                 if( pelty > 0 ) then

                    pnode = nnode(pelty)
                    pgaub = ngaus(pblty) 
                    pevat = solve(ivari_nsi) % ndofn * pnode
                    pmate = lmate(ielem)

                    if( pmate /= -1 ) then
                       !
                       ! Gather operations: ELVEL, ELCOD, BOVEL
                       !
                       do inode = 1,pnode
                          ipoin = lnods(inode,ielem)
                          elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
                          elcod(1:ndime,inode) = coord(1:ndime,ipoin)    
                       end do

                       do inodb = 1,pnodb
                          ipoin = lnodb(inodb,iboun)
                          bocod(1:ndime,inodb) = coord(1:ndime,ipoin)    
                          bovel(1:ndime,inodb) = veloc(1:ndime,ipoin,1)
                       end do

                       do igaub = 1,pgaub

                          call bouder(&
                               pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                               bocod,baloc,eucta)                                   ! and Jacobian
                          call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                          if ( kfl_waexl_ker == 0_ip ) then
                             !
                             ! normal behaviour
                             !
                             do idime = 1,ndime                                     ! Velocity
                                vewal(idime) = 0.0_rp
                             end do
                             do inodb = 1,pnodb
                                do idime = 1,ndime
                                   vewal(idime) = vewal(idime) &
                                        + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                                end do
                             end do

                          else
                             !
                             ! using velocity from exchange location
                             !
                             if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
                                do idime = 1,ndime
                                   vewal(idime) = velel_ker(idime,lexlo_ker(igaub,iboun))   ! Velocity - this comes as an input to this subroutine
                                end do
                             end if

                          end if

                          do idime = 1,ndime
                             tvelo(idime) = vewal(idime)
                             do kdime = 1,ndime
                                tvelo(idime) = tvelo(idime)   &
                                     - baloc(idime,ndime)     &
                                     * baloc(kdime,ndime) * vewal(kdime)
                             end do
                          end do

                          velav_ker(1:ndime,igaub,iboun) = tvelo(1:ndime)  

                          !  no longer used -- ! related to L184 above
                          !                       if ( kfl_fixbo_nsw_ker(iboun) == 1_ip) then ! no slip wall law - an elemental value is obtained averaging over all boundary gauss points relaed to the element
                          !                          lnsw_exch(ielem)%vel_aux(1:ndime) = lnsw_exch(ielem)%vel_aux(1:ndime) + tvelo(1:ndime)                       
                          !                       else    ! rest of cases normal behaviour
                          !                          velav_ker(1:ndime,igaub,iboun) = tvelo(1:ndime)
                          !                       end if

                       end do

                    end if

                 end if

              end if

           end do boundaries3
        end if
        call nsi_nsw_averaging(1_ip,avwei,dt)  ! Obtains lnsw_exch(ielem)%velav(1:ndime) -- here avwei is not used due to 1_ip - this must be done with or without restart

     end if

  else if ( itask == 4_ip ) then

     !-----------------------------------------------------------------------------
     !
     ! Calculate initial value of average velocity at the boundary nodes for post-processing
     !
     !-----------------------------------------------------------------------------
     if (kfl_rstar == 0) then
        if (kfl_noslw_ker == 0_ip) then
           do ipoin=1,npoin
              ibopo=lpoty(ipoin)
              if (ibopo > 0_ip) then
                 avupo_ker(1:ndime,ipoin) =  veloc(1:ndime,ipoin,1)
              end if
           end do
        else
           do ipoin=1,npoin
              if (kfl_nswpo_ker(ipoin) > 0_ip) then
                 avupo_ker(1:ndime,ipoin) =  veloc(1:ndime,ipoin,1)
              end if
           end do
        end if
     end if

  end if
  
end subroutine nsi_wallav


subroutine nsi_nsw_averaging(kfl_initilization,avwei,dt)
  use def_master,             only : ip,rp
  use def_kermod,             only : kfl_noslw_ker,lnsw_exch,kfl_fixbo_nsw_ker,velav_ker
  use def_domain,             only : ltypb,lelbo,ltype,ngaus,lmate,nelem,nboun,ndime

  implicit none


  integer(ip),intent(in)    :: kfl_initilization
  real(rp),intent(in)       :: avwei,dt   ! e=dt/T, weighting function for averaging
  integer(ip)               :: iboun,ielem,pblty,pelty,pmate,pgaub,igaub
  !
  ! Averaging for no slip wall law
  ! first add all velav_ker in lnsw_exch(ielem)%vel_aux(1:ndime)
  ! then divide by number of boundary gauss points associated to element
  !
  if (kfl_noslw_ker == 1_ip ) then
     boundaries3b: do iboun = 1,nboun
        if(  kfl_fixbo_nsw_ker(iboun) == 1_ip) then   ! no slip wall law
           !
           ! Element properties and dimensions
           !
           pblty = ltypb(iboun) 
           ielem = lelbo(iboun)
           pelty = ltype(ielem)
           if( pelty > 0 ) then
              pgaub = ngaus(pblty) 
              pmate = lmate(ielem)
              lnsw_exch(ielem)%vel_aux(1:ndime) = 0.0_rp
              if( pmate /= -1 ) then
                 do igaub = 1,pgaub
                    lnsw_exch(ielem)%vel_aux(1:ndime) = lnsw_exch(ielem)%vel_aux(1:ndime) + velav_ker(1:ndime,igaub,iboun)
                 end do
              end if
           end if
        end if
     end do boundaries3b



     if (kfl_initilization==1_ip) then
        elements3: do ielem = 1,nelem
           if ( lnsw_exch(ielem)%nbogp /= 0_ip ) then
              lnsw_exch(ielem)%velav(1:ndime) = lnsw_exch(ielem)%fact * lnsw_exch(ielem)%vel_aux(1:ndime)
           end if
        end do elements3
     else    ! This option is no longer used -- eliminate it
        elements: do ielem = 1,nelem
           if ( lnsw_exch(ielem)%nbogp /= 0_ip ) then
              lnsw_exch(ielem)%velav(1:ndime) = dt * lnsw_exch(ielem)%fact * lnsw_exch(ielem)%vel_aux(1:ndime) +avwei*lnsw_exch(ielem)%velav(1:ndime)
              lnsw_exch(ielem)%vel_aux(1:ndime) = lnsw_exch(ielem)%fact * lnsw_exch(ielem)%vel_aux(1:ndime)
           end if
        end do elements
     end if

  end if

end subroutine nsi_nsw_averaging
