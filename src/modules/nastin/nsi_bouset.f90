!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouset(kbsec,kbset)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouset
  ! NAME 
  !    nsi_bouset
  ! DESCRIPTION
  !    This routine computes variables on a boundary set W.
  !    The variable are: 
  !     0       SETSU: surface            =  int_W  = meas(W)      
  !     2       SETMP: mean pressure      =  int_W  P/meas(W)             
  !     3       SETMA: mass               =  int_W  rho*u.n                 
  !     3 -> 5  SETFV: viscous force      =  int_W  [2 mu E(u)].n          
  !     6 -> 8  SETFP: pressure force     =  int_W  [-pI].n                 
  !     9 -> 11 SETTV: viscous torque     =  int_W  (r-rc) x ( [2 mu E(u)].n ) 
  !    12 -> 14 SETTP: pressure torque    =  int_W  (r-rc) x ( [-pI].n ) 
  !    15       SETYP: mean y+            =  int_W  y+/meas(W)                
  !    16       SETMV: mean velocity      =  int_W  u.n/meas(W)               
  !    17 -> 19 SETWV: viscous wet force  =  int_Ww [2 mu E(u)].n  
  !    20 -> 22 SETWP: pressure wet force =  int_Ww [-pI].n     
  !    23       SETWS: wet surface        =  int_Ww 
  !    24 -> 26 SETIN: Internal force     =  Sum_{i in W} ( Auu.u + Aup.p - bu ) |_i
  !    27 -> 30 SETRE: Reattachment       =  
  !    31       SEDOT: u.n                =  int_W  u.n
  !    32 -> 34 SEPOR: Porous force     =  Sum_{i in W} ( bupor_nsi ) |_i
  !
  !    Force and torque are exerted by the solid on the fluid.
  !    Values of SETMP and SETYP are averaged further on in nsi_outset
  !
  ! USES 
  !    bouder
  !    chenor
  ! USED BY
  !    nsi_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_nsi_frixgb,     only : nsi_ml_ustar_all,ustars     ! machine learning
  use mod_bouder
  implicit none

  integer(ip), intent(in)  :: kbsec,kbset
  real(rp),    pointer     :: setsu(:),setmp(:)
  real(rp),    pointer     :: setfv(:),setfp(:)
  real(rp),    pointer     :: settv(:),settp(:)
  real(rp),    pointer     :: setyp(:),setmv(:)
  real(rp),    pointer     :: setwv(:),setwp(:)
  real(rp),    pointer     :: setws(:),setin(:)
  real(rp),    pointer     :: setre(:),setma(:)
  real(rp),    pointer     :: sedot(:),sepor(:)
  real(rp),    pointer     :: selin(:)
  real(rp),    pointer     :: tau_wall(:,:)
  real(rp),    pointer     :: setac(:)
  integer(ip)              :: ielem,inode,ipoin,igaus,idime,nn,kfl_force
  integer(ip)              :: pnode,pgaus,iboun,igaub,inodb,dummi
  integer(ip)              :: pelty,pblty,pnodb,pgaub,pmate,jdime,imate
  integer(ip)              :: jpoin,izdom,jzdom,ipoin1,ipoin2
  integer(ip)              :: rank_max_owner
  real(rp)                 :: baloc(ndime,ndime),bopre(mnodb)
  real(rp)                 :: bovel(ndime,mnodb,2),elvel(ndime,mnode)
  real(rp)                 :: bovfi(ndime,mnodb)     
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                 :: cartb(ndime,mnode),gpcar(ndime,mnode,mgaus)
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime) 
  real(rp)                 :: shapp(mnode),gbcoo(3),grave(ndime,ndime)
  real(rp)                 :: gblev(mgaus)
  real(rp)                 :: gbpre(max(mgaus,mgaub)),hleng(3),gbgve(ndime,ndime,mgaus)
  real(rp)                 :: gbsur,eucta,gbden(max(mgaus,mgaub)),gbvis(max(mgaus,mgaub))
  real(rp)                 :: gbvel(ndime,max(mgaus,mgaub),2),tauwa,velno,xfact
  real(rp)                 :: gbvdt(ndime,max(mgaus,mgaub))            ! tangent component of velocity - prescribed velocity.
  real(rp)                 :: tmpfv(3),tmpfp(3)
  real(rp)                 :: tmptv(3),tmptp(3),setdv(3),setdp(3)
  real(rp)                 :: dummr(3),nu,ustar,tragl(9),veaux(3)
  real(rp)                 :: fauxi(3),fauxn,taudot,taudotmax
  real(rp)                 :: tanfo(3)
  real(rp)                 :: velex(3),veave(3)
  real(rp)                 :: norfo,rho_air,gbvis_aux, roughness,udotn
  logical(lg)              :: lauxi

  if( INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! Initialization
     !
     !----------------------------------------------------------------------

     nullify(tau_wall)

     nn    =  postp(1) % nvabs + 1
     setsu => vbset( nn:nn , kbset ) ! Surface    
     setmp => vbset(  1: 1 , kbset ) ! Mean pressure     
     setma => vbset(  2: 2 , kbset ) ! Mass              
     setfv => vbset(  3: 5 , kbset ) ! Viscous force     
     setfp => vbset(  6: 8 , kbset ) ! Pressure force    
     settv => vbset(  9:11 , kbset ) ! Viscous torque    
     settp => vbset( 12:14 , kbset ) ! Pressure torque   
     setyp => vbset( 15:15 , kbset ) ! Mean y+           
     setmv => vbset( 16:16 , kbset ) ! Mean velocity     
     setwv => vbset( 17:19 , kbset ) ! Viscous wet force 
     setwp => vbset( 20:22 , kbset ) ! Pressure wet force
     setws => vbset( 23:23 , kbset ) ! Wet surface       
     setin => vbset( 24:26 , kbset ) ! Internal force
     setre => vbset( 27:30 , kbset ) ! Reattachment
     sedot => vbset( 31:31 , kbset ) ! udotn
     sepor => vbset( 32:34 , kbset ) ! Porous force
     setac => vbset( 35:35 , kbset ) ! Acceleration
     selin => vbset( 36:38 , kbset ) ! Linear momentum
     setsu =  0.0_rp
     setmp =  0.0_rp
     setma =  0.0_rp
     setfv =  0.0_rp
     setfp =  0.0_rp
     settv =  0.0_rp
     settp =  0.0_rp
     setyp =  0.0_rp
     setmv =  0.0_rp
     setwv =  0.0_rp
     setwp =  0.0_rp
     setws =  0.0_rp
     setin =  0.0_rp
     setdv =  0.0_rp
     setdp =  0.0_rp
     setre =  0.0_rp
     sedot =  0.0_rp
     sepor =  0.0_rp
     setac =  0.0_rp
     selin = 0.0_rp
     !
     ! Others
     !
     gbgve = 0.0_rp   ! Needed for Smagorinsky
     if(  postp(1) % npp_setsb( 3) /= 0 .or. &
          postp(1) % npp_setsb( 9) /= 0 .or. &
          postp(1) % npp_setsb(15) /= 0 ) then
        kfl_force = 1 ! Force is needed to compute force!, torque and wet force
     else
        kfl_force = 0
     end if
     !
     ! Reattachment
     !
     if( postp(1) % npp_setsb(27) /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'TAU_WALL','nsi_bouset',tau_wall,ndime+1,npoin)
     end if 

     boundaries: do iboun = 1,nboun

        if( lbset(iboun) == kbsec ) then

           !----------------------------------------------------------------
           !
           ! Element properties and dimensions and gather
           !
           !----------------------------------------------------------------

           pblty = ltypb(iboun) 
           pnodb = nnode(pblty)
           pgaub = ngaus(pblty)
           pmate = 1

           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              if (kfl_adj_prob == 0) then
                bopre(inodb) = press(ipoin,1)
              else
                bopre(inodb) = press_forw(ipoin,1)
              endif
              do idime = 1,ndime
                 bocod(idime,inodb) = coord(idime,ipoin)
                 if (kfl_adj_prob == 0) then
                   bovel(idime,inodb,1) = veloc(idime,ipoin,1)
                 else
                   bovel(idime,inodb,1) = veloc_forw(idime,ipoin,1)
                 endif
              end do              
              if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then 
                 do idime = 1,ndime
                    bovfi(idime,inodb) = velom(idime,ipoin)    ! see comments in nsi_bouope
                 end do
              else
                 do idime = 1,ndime
                    bovfi(idime,inodb) = 0.0_rp
                 end do
              end if
           end do
           if( postp(1) % npp_setsb(35) /= 0 ) then
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 bovel(1:ndime,inodb,2) = veloc(1:ndime,ipoin,3)
              end do
           end if
           
           ielem = lelbo(iboun)
           pelty = ltype(ielem)
           if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    if (kfl_adj_prob == 0) then
                      elvel(idime,inode) = veloc(idime,ipoin,1)
                    else
                      elvel(idime,inode) = veloc_forw(idime,ipoin,1)
                    endif
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              !
              ! Element length HLENG
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                   hnatu(pelty),hleng)

              !----------------------------------------------------------------
              !
              ! Values at Gauss points: GBPRE, GBVEL, GBLEV, GBTEM
              !
              !----------------------------------------------------------------
!               do inodb = 1,pnodb
!                  ipoin = lnodb(inodb,iboun)
!                    if (kfl_adj_prob == 0 .and. ipoin == 600 ) bopre(inodb) = bopre(inodb) + 0.00000001_rp
!               end do
              do igaub = 1,pgaub
                 gbpre(igaub) = 0.0_rp
                 do idime = 1,ndime
                    gbvel(idime,igaub,1) = 0.0_rp
                    gbvdt(idime,igaub)   = 0.0_rp
                 end do
                 !              do inodb = 1,pnodb
                 !                 ipoin = lnodb(inodb,iboun)
                 !                 gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * press(ipoin,1)
                 !                 do idime = 1,ndime
                 !                    gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * veloc(idime,ipoin,1)
                 !                 end do
                 !              end do
                 ! OJO no veo porqeu guillaume lo ha hecho con press y veloc cuando ya ha calculado bopre y bovel
                 ! yo prefiero hacerlo con bopre bovel  y para la differencia de veloc con bovfi(porque este último este
                 ! en eje global, lo he rotado más arriba
                 do inodb = 1,pnodb
                    gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * bopre(inodb)
                    do idime = 1,ndime
                       gbvel(idime,igaub,1) = gbvel(idime,igaub,1) + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb,1)
                       gbvdt(idime,igaub)   = gbvdt(idime,igaub)   + elmar(pblty)%shape(inodb,igaub) * ( bovel(idime,inodb,1)  &
                            &    - bovfi(idime,inodb) )  ! for the moment it includes normal component (substracted later)
                    end do
                 end do                 
              end do
              if( postp(1) % npp_setsb(35) /= 0 ) then
                 do igaub = 1,pgaub
                    gbvel(:,igaub,2) = 0.0_rp
                    do inodb = 1,pnodb
                       gbvel(1:ndime,igaub,2) = gbvel(1:ndime,igaub,2) + elmar(pblty)%shape(inodb,igaub) * bovel(1:ndime,inodb,2)
                    end do
                 end do
              end if
              if( kfl_hydro_nsi /= 0 ) then
                 do igaub = 1,pgaub
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       gbpre(igaub) = gbpre(igaub) &
                            + elmar(pblty)%shape(inodb,igaub) * bpess_nsi(1,ipoin,1)
                    end do
                 end do
              end if
              if( kfl_colev_nsi /= 0 ) then               ! Level set coupling: needs level set: GBLEV
                 do igaub = 1,pgaub
                    gblev(igaub) = 0.0_rp
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       gblev(igaub) = gblev(igaub) + elmar(pblty)%shape(inodb,igaub) * fleve(ipoin,kfl_colev_nsi)
                    end do
                 end do
              end if
              !
              ! Properties: mu and rho 
              !
              call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
              call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis)

              if( kfl_force == 1 ) then
                 do igaus = 1,pgaus
                    call elmder(&      
                         pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&            ! Cartesian derivative
                         elcod,gpcar(1,1,igaus),dummr,xjacm,xjaci)              ! and Jacobian
                 end do
              end if

              !----------------------------------------------------------------
              !
              ! Loop over Gauss points
              !
              !----------------------------------------------------------------

              gauss_points: do igaub = 1,pgaub

                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                      bocod,baloc,eucta)  
                 gbsur = elmar(pblty)%weigp(igaub)*eucta
                 setsu = setsu + gbsur
                 call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                 !if(kbsec==1) then
                 !   print*,'a=',iboun,igaub,gbsur,pnodb,pblty,setsu
                 !   if(igaub==9)stop
                 !end if
                 !
                 ! Velocity gradients grad(u): GRAVE(i,j) = duj/dxi
                 !
                 if( kfl_force == 1 ) then
                    call cartbo(&
                         1_ip,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                         elmar(pelty)%shaga,gpcar,elmar(pelty)%shape,&
                         shapp,cartb,pnodb,pnode,pgaus)
                    do idime = 1,ndime
                       do jdime = 1,ndime
                          grave(idime,jdime) = 0.0_rp
                          do inode = 1,pnode
                             grave(idime,jdime) = grave(idime,jdime) &
                                  + cartb(idime,inode) * elvel(jdime,inode)
                          end do
                       end do
                    end do
                 end if
                 
                 !-------------------------------------------------------------
                 !
                 ! rho * du/dt
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(35) /= 0 ) then
                    dummr(1:ndime) = gbvel(1:ndime,igaub,2)-gbvel(1:ndime,igaub,1)
                    setac = setac + gbden(igaub) * sqrt(dot_product(dummr(1:ndime),dummr(1:ndime))) * dtinv_nsi * gbsur
                 end if
                 
                 !-------------------------------------------------------------
                 !
                 ! Mean pressure
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(1) /= 0 ) then
                    setmp = setmp + gbpre(igaub) * gbsur
                 end if

                 !-------------------------------------------------------------
                 !
                 ! Mass
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(2) /= 0 ) then
                    xfact = gbden(igaub) * gbsur
                    do idime = 1,ndime
                       setma = setma + xfact * gbvel(idime,igaub,1) * baloc(idime,ndime)
                    end do
                 end if

                 !-------------------------------------------------------------
                 !
                 ! u.n
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(31) /= 0 ) then
                    sedot = sedot + gbsur * dot_product(gbvel(1:ndime,igaub,1),baloc(1:ndime,ndime))
                 end if

                 !-------------------------------------------------------------
                 !
                 ! rho * u * u.n
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(36)/=0 ) then
                    udotn        = gbden(igaub) * dot_product(gbvel(1:ndime,igaub,1),baloc(1:ndime,ndime))
                    selin(    1) = selin(    1) + gbsur * udotn * gbvel(    1,igaub,1)
                    selin(    2) = selin(    2) + gbsur * udotn * gbvel(    2,igaub,1)
                    if( ndime == 3 ) selin(ndime) = selin(ndime) + gbsur * udotn * gbvel(ndime,igaub,1)
                 end if

                 if( kfl_force == 1 ) then

                    tmpfv(1) = 0.0_rp              
                    tmpfp(1) = 0.0_rp
                    tmpfv(2) = 0.0_rp              
                    tmpfp(2) = 0.0_rp
                    tmpfv(3) = 0.0_rp              
                    tmpfp(3) = 0.0_rp

                    !----------------------------------------------------------
                    !
                    ! Pressure force: Fp = -p.n
                    !
                    !----------------------------------------------------------

                    do idime = 1,ndime  
                       tmpfp(idime) = -gbpre(igaub)*baloc(idime,ndime)
                    end do

                    !----------------------------------------------------------
                    !
                    ! Viscous force: Fv
                    !
                    !----------------------------------------------------------

                    if( kfl_fixbo_nsi(iboun) == 3 ) then
                       !
                       ! Law of the wall: F = - rho * (u*)^2 * (u_tan-u_fix_tan)/|u_tan-u_fix_tan|
                       !
                       if( kfl_rough > 0 ) then ! variable
                          roughness = 0.0_rp
                          do inodb = 1,pnodb
                             ipoin = lnodb(inodb,iboun)
                             roughness = roughness + rough(ipoin) * elmar(pblty)%shape(inodb,igaub)
                          end do
                       else    ! constant value or not value
                          roughness = rough_dom
                       end if

                       if ( kfl_waexl_ker == 1_ip ) then !if exchange location for wall law
                          velex(1:ndime) = velel_ker(1:ndime,lexlo_ker(igaub,iboun))
                          gbvdt(1:ndime,igaub) = velex(1:ndime)
                       else
                          velex(1:ndime) = 0.0_rp
                       end if

                       if ( kfl_wlaav_ker == 1_ip ) then !if exchange location for wall law
                          veave(1:ndime) = velav_ker(1:ndime,igaub,iboun)
                       else
                          veave(1:ndime) = 0.0_rp
                       end if

		               if (kfl_mlwm_ker == 1_ip) then
                         if(.not.allocated(ustars)) then 
                           call nsi_bouope_all(1_ip)
                         end if
                         ustar = ustars(igaub,iboun)
                       end if
                       call nsi_bouwal(&                        
                            2_ip,1_ip,pnodb,dummi,iboun,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                            bovel,bovfi,dummr,gbvis(igaub),gbden(igaub),baloc,ustar,dummr,roughness, dummr(1), &
                            velex,veave, igaub,lelbo(iboun))

                       do idime = 1,ndime              ! Substract normal component from gbvdt
                          veaux(idime) = gbvdt(idime,igaub)
                       end do
                       do idime = 1,ndime     
                          do jdime = 1,ndime
                             gbvdt(idime,igaub) = gbvdt(idime,igaub)   &
                                  - baloc(idime,ndime) &
                                  * baloc(jdime,ndime) * veaux(jdime)
                          end do
                       end do

                       call vecnor(gbvdt(1,igaub),ndime,velno,2_ip)
                       if( velno == 0.0_rp ) velno = 1.0_rp
                       tauwa = gbden(igaub) * ustar * ustar
                       do idime = 1,ndime
                          tmpfv(idime) = - tauwa * gbvdt(idime,igaub) / velno
                       end do
                       !
                       ! Add the normal component of the traction vector  F
                       ! (F.n)n = (n.sig.n)n = {n.(mu * [ grad(u) + grad(u)^t ] . n )}n
                       ! because the wall law only deals with tangential component
                       !
                       do idime = 1,ndime   ! obtain F
                          fauxi(idime) = 0.0_rp
                          do jdime = 1,ndime
                             fauxi(idime) = fauxi(idime)   &
                                  + gbvis(igaub) * ( grave(jdime,idime) + grave(idime,jdime) )  &
                                  * baloc(jdime,ndime)
                          end do
                       end do
                       fauxn = 0.0_rp ! F.n
                       do idime = 1,ndime
                          fauxn = fauxi(idime)* baloc(idime,ndime)
                       end do
                       do idime = 1,ndime    ! add (F.n)n
                          tmpfv(idime) = tmpfv(idime)    &
                               + fauxn * baloc(idime,ndime)
                       end do

                    else
                       lauxi = .false.
                       if( kfl_noslw_ker /= 0 ) then
                          if( kfl_fixbo_nsw_ker(iboun) == 1_ip ) lauxi = .true. ! no slip wall law uses the variational force not the one from gradients
                       end if
                       
                       if(.not.lauxi) then   ! gradient based traction   
                          !
                          ! No-slip: F = sig.n = mu * [ grad(u) + grad(u)^t ] . n  ! actually all non-wall law cases (Slip etc)
                          !
                          tauwa = 0.0_rp
                          gbvis_aux = gbvis(igaub)
                          if( kfl_noslw_ker /= 0 ) gbvis_aux = gbvis_aux + el_nsw_visc(ielem)  !this will be no longer needed because I will obtain it variationally
                          do idime = 1,ndime
                             do jdime = 1,ndime
                                tmpfv(idime) = tmpfv(idime)&
                                     + gbvis_aux * ( grave(jdime,idime) + grave(idime,jdime) )&
                                     * baloc(jdime,ndime)
                             end do
                             tauwa = tauwa + tmpfv(idime) * tmpfv(idime)
                          end do
                          tauwa = sqrt(tauwa)
                          ustar = sqrt(tauwa/gbden(igaub))
                       end if

                    end if

                    setfv(1) = setfv(1) + gbsur * tmpfv(1)
                    setfp(1) = setfp(1) + gbsur * tmpfp(1)
                    setfv(2) = setfv(2) + gbsur * tmpfv(2)
                    setfp(2) = setfp(2) + gbsur * tmpfp(2)
                    setfv(3) = setfv(3) + gbsur * tmpfv(3)
                    setfp(3) = setfp(3) + gbsur * tmpfp(3)
                    
                    !----------------------------------------------------------
                    !
                    ! Torque: Tv and Tp
                    !
                    !----------------------------------------------------------

                    if( postp(1) % npp_setsb(9) /= 0 ) then
                       do idime = 1,3
                          tmptv(idime) = 0.0_rp
                          tmptp(idime) = 0.0_rp
                          gbcoo(idime) = 0.0_rp
                       end do
                       do inodb = 1,pnodb
                          ipoin = lnodb(inodb,iboun)
                          do idime = 1,ndime
                             gbcoo(idime) = gbcoo(idime) &
                                  + elmar(pblty)%shape(inodb,igaub) &
                                  * coord(idime,ipoin)
                          end do
                       end do
                       do idime = 1,ndime
                          gbcoo(idime) = gbcoo(idime) - postp(1) % pabse(idime,9)
                          tmptv(idime) = gbcoo(idime)
                          tmptp(idime) = gbcoo(idime)
                       end do
                       call vecpro(tmptv,tmpfv,tmptv,3)       ! Viscous torque:  (r-rc) x Fv
                       call vecpro(tmptp,tmpfp,tmptp,3)       ! Pressure torque: (r-rc) x Fp
                       settv(1) = settv(1) + gbsur * tmptv(1)
                       settp(1) = settp(1) + gbsur * tmptp(1)
                       settv(2) = settv(2) + gbsur * tmptv(2)
                       settp(2) = settp(2) + gbsur * tmptp(2)
                       settv(3) = settv(3) + gbsur * tmptv(3)
                       settp(3) = settp(3) + gbsur * tmptp(3)
                    end if

                 end if

                 !-------------------------------------------------------------
                 !
                 ! Wet force and/or wet surface
                 !
                 !-------------------------------------------------------------

                 if( kfl_colev_nsi /= 0 .and. ( postp(1) % npp_setsb(17) /= 0 .or. postp(1) % npp_setsb(23) /= 0 ) ) then
                    imate = 1
                    if( nmate > 1 ) then
                       call runend('NSI_ANALYTICAL_PRESSURE: DO SOMETHING WHEN HAVING MORE THAN 1 MATERIAL')
                    end if
                    rho_air = densi_ker % rlaws(2,imate)
                    if( gblev(igaub) >= 0.0_rp ) then ! We are wet (in water)
                       !                    setwv(1) = setwv(1) + gbsur * tmpfv(1)
                       !                    setwp(1) = setwp(1) + gbsur * tmpfp(1)
                       !                    setwv(2) = setwv(2) + gbsur * tmpfv(2)
                       !                    setwp(2) = setwp(2) + gbsur * tmpfp(2)
                       !                    setwv(3) = setwv(3) + gbsur * tmpfv(3)
                       !                    setwp(3) = setwp(3) + gbsur * tmpfp(3)
                       setws    = setws    + gbsur
                    else ! We are dry (in air)                           ! approximation for the forces in the 'air'
                       setdv(1) = setdv(1) + gbsur * tmpfv(1) * rho_air / gbden(igaub)
                       setdp(1) = setdp(1) + gbsur * tmpfp(1) * rho_air / gbden(igaub)
                       setdv(2) = setdv(2) + gbsur * tmpfv(2) * rho_air / gbden(igaub)
                       setdp(2) = setdp(2) + gbsur * tmpfp(2) * rho_air / gbden(igaub)
                       setdv(3) = setdv(3) + gbsur * tmpfv(3) * rho_air / gbden(igaub)
                       setdp(3) = setdp(3) + gbsur * tmpfp(3) * rho_air / gbden(igaub)
                    end if
                 end if

                 !-------------------------------------------------------------
                 !
                 ! Mean y+ (ustar was calculated previously)
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(15) /= 0 ) then
                    nu    = gbvis(igaub) / gbden(igaub)
                    if (kfl_delta == 1) then
                       setyp = setyp + ywalb(iboun) * ustar / nu * gbsur
                    else
                       setyp = setyp + delta_nsi * ustar / nu * gbsur
                    end if
                 end if

                 !-------------------------------------------------------------
                 !
                 ! Mean normal velocity
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(16) /= 0 ) then
                    do idime = 1,ndime
                       setmv = setmv + gbsur * gbvel(idime,igaub,1) * baloc(idime,ndime)
                    end do
                 end if

                 !-------------------------------------------------------------
                 !
                 ! Reattachment
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(27) /= 0 ) then

                    ipoin = lnodb(1,iboun)
                    if( sum(kfl_fixno_nsi(1:ndime,ipoin)) /= ndime .and. kfl_fixbo_nsi(iboun) /= 3 ) then
                       !
                       ! I guess we have slip condition here
                       ! We then look for a change of sign in the tangential velocity
                       !
                       norfo = dot_product(gbvel(1:ndime,igaub,1),baloc(1:ndime,ndime))
                       tanfo(1:ndime) = gbvel(1:ndime,igaub,1) - baloc(1:ndime,ndime) * norfo                       
                    else
                       !
                       ! I guess we have a wall law here or a no-slip condition
                       ! We then look for the change of sign of the wall shear stress
                       ! Wall shear stress = sig.n - (n.sig.n) n
                       !
                       norfo = dot_product(tmpfv(1:ndime),baloc(1:ndime,ndime))
                       tanfo(1:ndime) = tmpfv(1:ndime) - baloc(1:ndime,ndime) * norfo
                    end if
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       tau_wall(1:ndime,ipoin) = tau_wall(1:ndime,ipoin) &
                            + tanfo(1:ndime) * elmar(pblty) % shape(inodb,igaub) * gbsur
                       tau_wall(ndime+1,ipoin) = tau_wall(ndime+1,ipoin) &
                            + elmar(pblty) % shape(inodb,igaub) * gbsur
                    end do
                 end if

              end do gauss_points

           end if

        end if

     end do boundaries
          
  end if
  
  
  
  !
  ! Reattachment
  !
  if( postp(1) % npp_setsb(27) /= 0 ) then
     if( INOTMASTER ) then
        !
        ! Project wall shear stress
        !
        call PAR_INTERFACE_NODE_EXCHANGE(tau_wall,'SUM','IN MY CODE') 
        do ipoin = 1,npoin
           if( tau_wall(ndime+1,ipoin) > zeror ) & 
                tau_wall(1:ndime,ipoin) = tau_wall(1:ndime,ipoin) / tau_wall(ndime+1,ipoin)
        end do
        taudotmax = -1.0_rp
        ipoin1    = 0
        ipoin2    = 0
        do ipoin = 1,npoin
           if( tau_wall(ndime+1,ipoin) > zeror ) then
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin  = c_dom(izdom)
                 if( jpoin /= ipoin .and. tau_wall(ndime+1,jpoin) > zeror ) then
                    taudot = dot_product(tau_wall(1:ndime,ipoin),tau_wall(1:ndime,jpoin))
                    if( taudot < -zeror ) then
                       !
                       ! tau_wall changes direction
                       ! Compute maximum gradient
                       !
                       taudot = 0.0_rp
                       do idime = 1,ndime
                          taudot = taudot + (tau_wall(idime,ipoin)-tau_wall(idime,jpoin))**2
                       end do
                       taudot = sqrt(taudot) / (0.5_rp*(tau_wall(ndime+1,ipoin)+tau_wall(ndime+1,jpoin)))
                       !
                       ! Keep the value with highest gradient
                       !
                       if( taudot > taudotmax ) then
                          ipoin1    = ipoin
                          ipoin2    = jpoin
                          taudotmax = taudot
                       end if
                    end if
                 end if
              end do
           end if
        end do
     end if
     !
     ! Check which worker has the maximum gradient value
     !
     call PAR_MAX(taudotmax,'IN MY CODE',rank_max_owner)
     
     if( kfl_paral == rank_max_owner .and. INOTMASTER ) then
        if( ipoin1 /= 0 .and. ipoin2 /= 0 ) then
           setre(1:ndime) =  0.5_rp*(coord(1:ndime,ipoin1)+coord(1:ndime,ipoin2))
           setre(ndime+1) =  taudotmax
        else
           setre(1:ndime) =  0.0_rp
           setre(ndime+1) = -1.0_rp
        end if
     end if
     call memory_deallo(mem_modul(1:2,modul),'TAU_WALL','nsi_bouset',tau_wall)
  end if

  if( INOTMASTER ) then
     !
     ! Level set
     !
     if( kfl_colev_nsi /= 0 .and. ( postp(1) % npp_setsb(17) /= 0 .or. postp(1) % npp_setsb(23) /= 0 ) ) then
        setwv(1) = setfv(1) - setdv(1)
        setwv(2) = setfv(2) - setdv(2)
        setwv(3) = setfv(3) - setdv(3)
        setwp(1) = setfp(1) - setdp(1)
        setwp(2) = setfp(2) - setdp(2)
        setwp(3) = setfp(3) - setdp(3)
     end if
     !
     ! Internal force
     !
     if(  postp(1) % npp_setsb(24) /= 0 .or. &
          postp(1) % npp_setsb(25) /= 0 .or. &
          postp(1) % npp_setsb(26) /= 0 ) then

        call memgen(1_ip,npoin,0_ip)
        gisca = 0_ip
        do iboun = 1,nboun
           if( lbset(iboun) == kbsec ) then
              lauxi = .false.        ! to allocate kfl_fixbo_nsw_ker only if kfl_noslw_ker /= 0
              if( kfl_noslw_ker /= 0 ) then
                 if( kfl_fixbo_nsw_ker(iboun) == 1_ip ) lauxi = .true. ! no slip wall law uses the variational force not the one from gradients
              end if
              do inodb = 1,lnnob(iboun)
                 ipoin = lnodb(inodb,iboun)            
                 if( lauxi ) then
                    gisca(ipoin) = 2_ip
                 else
                    gisca(ipoin) = 1_ip
                 end if
              end do
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
        do ipoin = 1,npoin
           if( gisca(ipoin) == 2 ) then ! related to kfl_fixbo_nsw_ker(iboun) == 1_ip - no slip wall law
              if (ipoin <=  npoin_own)  setin(1:ndime) = setin(1:ndime) + vafor_nsi(1:ndime,ipoin) !only own nodes to avoid repeating at interfase
           else if( gisca(ipoin) /= 0) then ! original behaviour
              if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then
                 do idime = 1,ndime
                    setin(idime) = setin(idime) - intfo_nsi(ipoin) % bu(idime)
                    jzdom = 0
                    do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                       jpoin = c_dom(izdom)
                       jzdom = jzdom + 1
                       do jdime = 1,ndime
                          setin(idime) = setin(idime) + intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) * veloc(jdime,jpoin,1) 
                       end do
                       setin(idime) = setin(idime) + intfo_nsi(ipoin) % Aup(idime,jzdom) * press(jpoin,1)                           
                    end do
                 end do
              end if
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
  end if

  if( INOTMASTER ) then
     !
     ! Porous force
     !
     if(  postp(1) % npp_setsb(32) /= 0 .or. &
          postp(1) % npp_setsb(33) /= 0 .or. &
          postp(1) % npp_setsb(34) /= 0 ) then
        sepor(1:ndime) = porfo_nsi(1:ndime)
     end if
  end if


  !write(6,*) '|------DEBUGDEBUGDEBUGDEBUG. kbsec: ', kbsec, 'kbset: ',kbset, 'setmv: ',setmv

end subroutine nsi_bouset
