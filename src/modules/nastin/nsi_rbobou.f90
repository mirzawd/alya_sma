!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_rbobou()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_rbobou
  ! NAME
  !    nsi_rbobou
  ! DESCRIPTION
  !    Compute force on rigid bodies
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
  use mod_gradie
  use mod_ker_proper
  use mod_messages, only : livinf
  use mod_bouder

  implicit none
  real(rp)             :: elvel(ndime,mnode)
  real(rp)             :: elcod(ndime,mnode)

  real(rp)             :: bocod(ndime,mnodb)
  real(rp)             :: bovel(ndime,mnodb)
  real(rp)             :: bovfi(ndime,mnodb)     
  real(rp)             :: bopre(mnodb)
  real(rp)             :: baloc(ndime,ndime)

  real(rp)             :: gbden(mgaus)
  real(rp)             :: gbvis(mgaus)
  real(rp)             :: gbpre(mgaus)
  real(rp)             :: gbvel(ndime,mgaus)
  real(rp)             :: gbvdt(ndime,mgaus)   ! tangent component of velocity - prescribed velocity.
  real(rp)             :: gbcoo(3)
  real(rp)             :: grave(3,3)
  real(rp)             :: cartb(ndime,mnode)
  real(rp)             :: gpcar(ndime,mnode,mgaus)

  integer(ip)          :: ielem,inode,ipoin,idime,iimbo,imeth,kgaus
  integer(ip)          :: pnode,iboun,igaub,inodb,dummi
  integer(ip)          :: pelty,pblty,pnodb,pgaub,pmate,igaus,pgaus
  integer(ip)          :: jdime,kboun,izdom,jzdom,jpoin,kauxi
  real(rp)             :: xjaci(9),xjacm(9),hleng(3),velno,tauwa,ustar
  real(rp)             :: gbsur,eucta,dummr
  real(rp)             :: tragl(9),F1v(3),F1p(3),T1p(3),T1v(3) 
  real(rp)             :: veaux(3),fauxi(3),fauxn, roughness
  real(rp),    pointer :: vpfor(:,:),vptor(:,:)
  real(rp),    pointer :: Fp(:),Fv(:),Tp(:),Tv(:)
  real(rp),    pointer :: pp_pf(:,:),pp_vf(:,:),pp_pt(:,:),pp_vt(:,:)

  integer(ip)          :: linau(99)   ! I understand the max num of sets is 99
  integer(ip)          :: ii,iiaux

  if( nrbod == 0 ) return  ! I understand that if there kfl_rigid_ale=0 =>nrbod=0  so it aonly enters here for kfl_rigid_ale/=0

!  COMO Ahoar estoy usando nrbod  en lugar de nimbo estas lineas deberian ser alpedo
!  if (.not.associated(rbbou)) stop 'return al comienzo de nsi_rbobou'
!  if (.not.associated(rbbou)) return ! needed to overcome the case when nrbod /= 0 due to ibm
!                                     ! the alternative would be to use a new variable ej nrbbo instead of nrbod  ! esto lo correg'i!!!!!!!

  call livinf(59_ip,'COMPUTE FORCES AND TORQUES ON IB',0_ip)

  !
  ! Initialization
  !
  imeth = 1
  do iimbo = 1,nrbod

     !
     ! auxiliar to postprocess forces separatelly by sets
     !
     do ii = 1,99
        linau(ii) = 0
     end do
     do ii = 1,rbbou(iimbo) % nrbse
        linau( rbbou(iimbo) % lrbse(ii) ) = ii
     end do

     Fv => rbbou(iimbo) % vforce
     Fp => rbbou(iimbo) % pforce
     Tv => rbbou(iimbo) % vtorqu
     Tp => rbbou(iimbo) % ptorqu

     pp_pf => rbbou(iimbo) % pp_pf   ! to postprocess forces separatelly by sets
     pp_vf => rbbou(iimbo) % pp_vf
     pp_pt => rbbou(iimbo) % pp_pt
     pp_vt => rbbou(iimbo) % pp_vt

     do idime = 1,3
        Fp(idime)  = 0.0_rp
        Fv(idime)  = 0.0_rp
        Tp(idime)  = 0.0_rp
        Tv(idime)  = 0.0_rp
        do ii =1,10
           pp_pf(idime,ii)   = 0.0_rp
           pp_vf(idime,ii)   = 0.0_rp
           pp_pt(idime,ii)   = 0.0_rp
           pp_vt(idime,ii)   = 0.0_rp
        end do
     end do
  end do

  if (ittim <= 1 ) return
  !
  ! Allocate memory: If gradients are smoothed
  !
  if( INOTMASTER ) then
     if( imeth == 1 .and. kfl_intfo_nsi == 0 ) then
        call memgen(0_ip,ntens,npoin)
        call graten(veloc,gevec)
     else if( kfl_intfo_nsi >= 1 ) then
        call memgen(1_ip,npoin,0_ip)
     end if
  end if
  !
  ! Loop over boundaries
  !
  if( INOTMASTER ) then

     kgaus = 0
     do iimbo = 1,nrbod

        Fv     => rbbou(iimbo) % vforce
        Fp     => rbbou(iimbo) % pforce
        Tv     => rbbou(iimbo) % vtorqu
        Tp     => rbbou(iimbo) % ptorqu

        pp_pf => rbbou(iimbo) % pp_pf   ! to postprocess forces separatelly by sets
        pp_vf => rbbou(iimbo) % pp_vf
        pp_pt => rbbou(iimbo) % pp_pt
        pp_vt => rbbou(iimbo) % pp_vt

        !-------------------------------------------------------------
        !
        ! Body fitted bodies
        !
        !-------------------------------------------------------------

        if( kfl_intfo_nsi >= 1 ) then
           !
           ! Force is computed using nodal force Fi: F = Sum_i Fi 
           ! Think how to implement forces by sets - what to do in nodes at the intersection
           !
           do kboun = 1,rbbou(iimbo) % nboib
              iboun = rbbou(iimbo) % lbinv(kboun)
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 gisca(ipoin) = 1
              end do
           end do
           do ipoin = 1,npoin
              if( gisca(ipoin) == 1 ) then
                 do idime = 1,ndime
                    F1v(idime) = intfo_nsi(ipoin) % bu(idime)
                    F1p(idime) = 0.0_rp
                    T1v(idime) = 0.0_rp
                    T1p(idime) = 0.0_rp
                    jzdom = 0
                    do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                       jpoin = c_dom(izdom)
                       jzdom = jzdom + 1
                       do jdime = 1,ndime
                          F1v(idime) = F1v(idime) - intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) * veloc(jdime,jpoin,1) 
                       end do
                       F1p(idime) = F1p(idime) - intfo_nsi(ipoin) % Aup(idime,jzdom) * press(jpoin,1)                           
                    end do
                 end do
                 if( kfl_hydro_gravity_nsi /= 0 ) then  
                    call runend('OUPS')
                    do idime = 1,ndime
                       jzdom = 0
                       do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                          jpoin = c_dom(izdom)
                          jzdom = jzdom + 1
                          F1p(idime) = F1p(idime) - intfo_nsi(ipoin) % Aup(idime,jzdom) * bpess_nsi(1,jpoin,1)         
                       end do
                    end do
                 end if
                 gisca(ipoin) = 0
                 !
                 ! obtain Torque
                 !
                 do idime = 1,ndime
                    veaux(idime) = coord(idime,ipoin) - rbbou(iimbo) % posil(idime,1)
                 end do
                 call vecpro(veaux,F1v,T1v,3_ip)       ! Viscous torque:  (r-rc) x Fv
                 call vecpro(veaux,F1p,T1p,3_ip)       ! Pressure torque: (r-rc) x Fp

                 iiaux = linau (lbpse(ipoin))          ! to postprocess forces separatelly by sets

                 do idime = 1,3

                    Fv(idime) =  Fv(idime) + F1v(idime)
                    Fp(idime) =  Fp(idime) + F1p(idime)
                    Tv(idime) =  Tv(idime) + T1v(idime)
                    Tp(idime) =  Tp(idime) + T1p(idime)

                    pp_vf(idime,iiaux) =  pp_vf(idime,iiaux) + F1v(idime)
                    pp_pf(idime,iiaux) =  pp_pf(idime,iiaux) + F1p(idime)
                    pp_vt(idime,iiaux) =  pp_vt(idime,iiaux) + T1v(idime)
                    pp_pt(idime,iiaux) =  pp_pt(idime,iiaux) + T1p(idime)

                 end do

              end if
           end do
           !
           ! For wall law we need to add the force due to the wall law as is done in the case (INTEG) because as the  
           ! tangential direction is free the residual formulation calculates no force in this direction. 
           ! We reapeat what is done in the INTEG case but only for wall law
           !
           kauxi = 0_ip
           boundaries3: do kboun = 1,rbbou(iimbo) % nboib
              iboun = rbbou(iimbo) % lbinv(kboun)
              if( kfl_fixbo_nsi(iboun) == 3 ) kauxi = kauxi + 1_ip
           end do boundaries3
           if (kauxi>0) then
              !
              ! Only wall law force is computed using traction: F = int_S sig.n ds
              !             
              boundaries4: do kboun = 1,rbbou(iimbo) % nboib

                 iboun = rbbou(iimbo) % lbinv(kboun)
                 pblty = ltypb(iboun) 
                 pnodb = nnode(pblty)
                 pgaub = ngaus(pblty)
                 ielem = lelbo(iboun)
                 pelty = ltype(ielem)
                 pmate = 1

                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    bopre(inodb) = press(ipoin,1)
                    do idime = 1,ndime
                       bocod(idime,inodb) = coord(idime,ipoin)
                       bovel(idime,inodb) = veloc(idime,ipoin,1)
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

                 if( pelty > 0 ) then

                    pnode = nnode(pelty)
                    pgaus = ngaus(pelty)
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       do idime = 1,ndime
                          elvel(idime,inode) = veloc(idime,ipoin,1)
                          elcod(idime,inode) = coord(idime,ipoin)
                       end do
                    end do
                    !
                    ! Element length HLENG
                    !
                    call elmlen(&
                         ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                         hnatu(pelty),hleng)
                    !
                    ! Values at Gauss points: GBPRE, GBVEL, GBLEV, GBTEM
                    !
                    do igaub = 1,pgaub
                       gbpre(igaub) = 0.0_rp
                       do idime = 1,ndime
                          gbvel(idime,igaub) = 0.0_rp
                          gbvdt(idime,igaub) = 0.0_rp
                       end do
                       !                          do inodb = 1,pnodb
                       !                             ipoin = lnodb(inodb,iboun)
                       !                             gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * press(ipoin,1)
                       !                             do idime = 1,ndime
                       !                                gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * veloc(idime,ipoin,1)
                       !                             end do
                       !                          end do
                       ! OJO no veo porqeu guillaume lo ha hecho con press y veloc cuando ya ha calculado bopre y bovel
                       ! yo prefiero hacerlo con bopre bovel  y para la differencia de veloc con bovfi(porque este último est en eje global, lo he rotado
                       ! más arriba
                       do inodb = 1,pnodb
                          gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * bopre(inodb)
                          do idime = 1,ndime
                             gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                             gbvdt(idime,igaub) = gbvdt(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * ( bovel(idime,inodb)  &
                                  &    - bovfi(idime,inodb) )  ! for the moment it includes normal component (substracted later)
                          end do
                       end do
                    end do
                    if( kfl_hydro_gravity_nsi /= 0 ) then
                       call runend('UPS 2')
                       do igaub = 1,pgaub
                          do inodb = 1,pnodb
                             ipoin = lnodb(inodb,iboun)
                             gbpre(igaub) = gbpre(igaub) &
                                  + elmar(pblty)%shape(inodb,igaub) * bpess_nsi(1,ipoin,1)
                          end do
                       end do
                    end if
                    !
                    ! Properties: mu and rho 
                    !
                    call ker_proper('DENSI','PGAUB',dummi,iboun,gbden,pnode,pgaub,elmar(pblty)%shape)
                    call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis,pnode,pgaub,elmar(pblty)%shape)       
                    !
                    ! Cartesian derivatives
                    !
                    do igaus = 1,pgaus
                       call elmder(&
                            pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&          ! Cartesian derivative
                            elcod,gpcar(1,1,igaus),dummr,xjacm,xjaci)            ! and Jacobian
                    end do
                    !
                    ! Loop over Gauss points
                    !
                    gauss_points: do igaub = 1,pgaub

                       call bouder(&
                            pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                            bocod,baloc,eucta)                                   ! and Jacobian
                       gbsur = elmar(pblty)%weigp(igaub)*eucta 
                       call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                       do idime = 1,3  
                          F1v(idime) = 0.0_rp              
                          T1v(idime) = 0.0_rp
                       end do
                       !
                       ! Viscous force: Fv
                       !
                       if( kfl_fixbo_nsi(iboun) == 3 ) then
                          !
                          ! Law of the wall: F = - rho * (u*)^2 * (u_tan-u_fix_tan)/|u_tan-u_fix_tan|
                          !
                          if( kfl_rough > 0 ) then
                             roughness = 0.0_rp
                             do inodb = 1,pnodb
                                ipoin = lnodb(inodb,iboun)
                                roughness = roughness + rough(ipoin) * elmar(pblty)%shape(inodb,igaub)
                             end do
                          else
                             roughness = rough_dom
                          end if
                          call nsi_bouwal(&                        
                               2_ip,1_ip,pnodb,dummi,iboun,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                               bovel,bovfi,dummr,gbvis(igaub),gbden(igaub),baloc,ustar,dummr,roughness, dummr,dummr,dummr, igaub,lelbo(iboun))

                          do idime = 1,ndime              ! Substract normal componenet from gbvdt
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
                          F1v(3) = 0.0_rp   ! so that in the 2d case it has the correct value
                          do idime = 1,ndime
                             F1v(idime) = - tauwa * gbvdt(idime,igaub) / velno
                          end do

                       end if
                       !
                       ! Torque: Tv
                       !
                       do idime = 1,3
                          gbcoo(idime) = 0.0_rp
                       end do
                       do inodb = 1,pnodb
                          ipoin = lnodb(inodb,iboun)
                          do idime = 1,ndime
                             gbcoo(idime) = gbcoo(idime) &
                                  + elmar(pblty) % shape(inodb,igaub) &
                                  * coord(idime,ipoin)
                          end do
                       end do
                       do idime = 1,ndime
                          gbcoo(idime) = gbcoo(idime) - rbbou(iimbo) % posil(idime,1)
                       end do
                       call vecpro(gbcoo,F1v,T1v,3_ip)       ! Viscous torque:  (r-rc) x Fv

                       do idime = 1,3
                          Fv(idime) =  Fv(idime) - gbsur * F1v(idime)
                          Tv(idime) =  Tv(idime) - gbsur * T1v(idime)
                       end do

                    end do gauss_points

                 end if

              end do boundaries4
           end if
           ! if( ittim <= 1 ) then
           !    Fv = 0.0_rp
           !    Fp = 0.0_rp
           ! end if
           !print*,'force a=',Fv(2),Fp(2)
           !   Fv = 0.0_rp
           !   Fp = 0.0_rp

        else
           !
           ! Forces is computed using traction: F = int_S sig.n ds
           !            
 !          print *,"DEBUG: EN EL ELSE DEL NSI_RIGID: "

           boundaries2: do kboun = 1,rbbou(iimbo) % nboib

              iboun = rbbou(iimbo) % lbinv(kboun)
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              pgaub = ngaus(pblty)
              ielem = lelbo(iboun)
              pelty = ltype(ielem)
              pmate = 1
              iiaux = linau (lbset(iboun))  ! to postprocess forces separatelly by sets

              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 bopre(inodb) = press(ipoin,1)
                 do idime = 1,ndime
                    bocod(idime,inodb) = coord(idime,ipoin)
                    bovel(idime,inodb) = veloc(idime,ipoin,1)
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

              if( pelty > 0 ) then
                 pnode = nnode(pelty)
                 pgaus = ngaus(pelty)
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1,ndime
                       elvel(idime,inode) = veloc(idime,ipoin,1)
                       elcod(idime,inode) = coord(idime,ipoin)
                    end do
                 end do
                 !
                 ! Element length HLENG
                 !
                 call elmlen(&
                      ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                      hnatu(pelty),hleng)
                 !
                 ! Values at Gauss points: GBPRE, GBVEL, GBLEV, GBTEM
                 !
                 do igaub = 1,pgaub
                    gbpre(igaub) = 0.0_rp
                    do idime = 1,ndime
                       gbvel(idime,igaub) = 0.0_rp
                       gbvdt(idime,igaub) = 0.0_rp
                    end do
                    !                       do inodb = 1,pnodb
                    !                          ipoin = lnodb(inodb,iboun)
                    !                          gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * press(ipoin,1)
                    !                          do idime = 1,ndime
                    !                             gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * veloc(idime,ipoin,1)
                    !                          end do
                    !                       end do
                    ! OJO no veo porqeu guillaume lo ha hecho con press y veloc cuando ya ha calculado bopre y bovel
                    ! yo prefiero hacerlo con bopre, bovel  y para la differencia de veloc con bovfi(porque este último esta 
                    ! en eje global, lo he rotado más arriba)
                    ! podría haber puesto al diff en gbvel directamente pero no porque se usa en nsi_elmpro
                    do inodb = 1,pnodb
                       gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * bopre(inodb)
                       do idime = 1,ndime
                          gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                          gbvdt(idime,igaub) = gbvdt(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * ( bovel(idime,inodb)  &
                               &  - bovfi(idime,inodb) )  ! for the moment it includes normal component (substracted later)
                       end do
                    end do
                 end do
                 if( kfl_hydro_gravity_nsi /= 0 ) then
                    call runend('OUPS 3')
                    do igaub = 1,pgaub
                       do inodb = 1,pnodb
                          ipoin = lnodb(inodb,iboun)
                          gbpre(igaub) = gbpre(igaub) &
                               + elmar(pblty)%shape(inodb,igaub) * bpess_nsi(1,ipoin,1)
                       end do
                    end do
                 end if
                 !
                 ! Properties: mu and rho 
                 !
                 call ker_proper('DENSI','PGAUB',dummi,iboun,gbden,pnode,pgaub,elmar(pblty)%shape)
                 call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis,pnode,pgaub,elmar(pblty)%shape)  
                 !
                 ! Cartesian derivatives
                 !
                 do igaus = 1,pgaus
                    call elmder(&
                         pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&          ! Cartesian derivative
                         elcod,gpcar(1,1,igaus),dummr,xjacm,xjaci)            ! and Jacobian
                 end do
                 !
                 ! Loop over Gauss points
                 !
                 gauss_points1: do igaub = 1,pgaub

                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                         bocod,baloc,eucta)                                   ! and Jacobian
                    gbsur = elmar(pblty)%weigp(igaub)*eucta 
                    call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
                    !
                    ! Velocity gradients grad(u): GRAVE
                    !
                    call cartbo(&
                         1_ip,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                         elmar(pelty)%shaga,gpcar,elmar(pelty)%shape,&
                         dummr,cartb,pnodb,pnode,pgaus)
                    do idime = 1,ndime
                       do jdime = 1,ndime
                          grave(idime,jdime) = 0.0_rp
                          do inode = 1,pnode
                             grave(idime,jdime) = grave(idime,jdime) &
                                  + cartb(idime,inode) * elvel(jdime,inode)
                          end do
                       end do
                    end do

                    do idime = 1,3  
                       F1v(idime) = 0.0_rp              
                       F1p(idime) = 0.0_rp
                       T1v(idime) = 0.0_rp
                       T1p(idime) = 0.0_rp                       
                    end do
                    !
                    ! Pressure force: Fp = -p.n
                    !
                    do idime = 1,ndime  
                       F1p(idime) = -gbpre(igaub)*baloc(idime,ndime)
                    end do
                    !
                    ! Viscous force: Fv
                    !
                    if( kfl_fixbo_nsi(iboun) == 3 ) then
                       !
                       ! Law of the wall: F = - rho * (u*)^2 * (u_tan-u_fix_tan)/|u_tan-u_fix_tan|
                       !
                       if( kfl_rough > 0 ) then
                          roughness = 0.0_rp
                          do inodb = 1,pnodb
                             ipoin = lnodb(inodb,iboun)
                             roughness = roughness + rough(ipoin) * elmar(pblty)%shape(inodb,igaub)
                          end do
                       else
                          roughness = rough_dom
                       end if
                       call nsi_bouwal(&                        
                            2_ip,1_ip,pnodb,dummi,iboun,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                            bovel,bovfi,dummr,gbvis(igaub),gbden(igaub),baloc,ustar,dummr,roughness, dummr,dummr,dummr, igaub,lelbo(iboun))

                       do idime = 1,ndime              ! Substract normal componenet from gbvdt
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
                          F1v(idime) = - tauwa * gbvdt(idime,igaub) / velno
                       end do
                       !
                       ! Add the normal component of the traction vector  F
                       ! (F.n)n = (n.sig.n)n = {n.(mu * [ grad(u) + grad(u)^t ] . n )}n
                       ! because the wall law only deals with tangential component
                       !
                       do idime = 1,ndime   ! obtain F
                          fauxi(idime) = 0.0_rp
                          do jdime = 1,ndime
                             fauxi(idime) = fauxi(idime)&
                                  + gbvis(igaub) * ( grave(jdime,idime) + grave(idime,jdime) )&
                                  * baloc(jdime,ndime)
                          end do
                       end do
                       fauxn = 0.0_rp ! F.n
                       do idime = 1,ndime
                          fauxn = fauxi(idime)* baloc(idime,ndime)
                       end do
                       do idime = 1,ndime    ! add (F.n)n
                          F1v(idime) = F1v(idime)&
                               + fauxn * baloc(idime,ndime)
                       end do

                    else
                       !
                       ! No-slip: F = sig.n = mu * [ grad(u) + grad(u)^t ] . n ! actually all non-wall law cases (Slip etc)
                       !
                       do idime = 1,ndime
                          do jdime = 1,ndime
                             F1v(idime) = F1v(idime)&
                                  + gbvis(igaub) * ( grave(jdime,idime) + grave(idime,jdime) )&
                                  * baloc(jdime,ndime)
                          end do
                       end do

                    end if
                    !
                    ! Torque: Tv and Tp
                    !
                    do idime = 1,3
                       gbcoo(idime) = 0.0_rp
                    end do
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       do idime = 1,ndime
                          gbcoo(idime) = gbcoo(idime) &
                               + elmar(pblty) % shape(inodb,igaub) &
                               * coord(idime,ipoin)
                       end do
                    end do
                    do idime = 1,ndime
                       gbcoo(idime) = gbcoo(idime) - rbbou(iimbo) % posil(idime,1)
                    end do
                    call vecpro(gbcoo,F1v,T1v,3_ip)       ! Viscous torque:  (r-rc) x Fv
                    call vecpro(gbcoo,F1p,T1p,3_ip)       ! Pressure torque: (r-rc) x Fp

                    do idime = 1,3

                       Fv(idime) =  Fv(idime) - gbsur * F1v(idime)
                       Fp(idime) =  Fp(idime) - gbsur * F1p(idime)
                       Tv(idime) =  Tv(idime) - gbsur * T1v(idime)
                       Tp(idime) =  Tp(idime) - gbsur * T1p(idime)

                       pp_vf(idime,iiaux) =  pp_vf(idime,iiaux) - gbsur * F1v(idime)
                       pp_pf(idime,iiaux) =  pp_pf(idime,iiaux) - gbsur * F1p(idime)
                       pp_vt(idime,iiaux) =  pp_vt(idime,iiaux) - gbsur * T1v(idime)
                       pp_pt(idime,iiaux) =  pp_pt(idime,iiaux) - gbsur * T1p(idime)

                    end do
                 end do gauss_points1
              end if
           end do boundaries2
        end if
     end do
  end if


  !write(99,*) Fv(2),Fp(2)
  !do iimbo = 1,nrbod
  !   Fv => rbbou(iimbo) % vforce
  !   Fp => rbbou(iimbo) % pforce
  !   Tv => rbbou(iimbo) % vtorqu
  !   Tp => rbbou(iimbo) % ptorqu
  !   do idime = 1,3
  !!      Fp(idime)  = 0.0_rp
  !!      Fv(idime)  = 0.0_rp
  !!      Tp(idime)  = 0.0_rp
  !!      Tv(idime)  = 0.0_rp
  !   end do
  !end do
  !
  ! Add buoyancy: OJO, poner en el bucle elemental mas arriba
  !
  !do iimbo = 1,nrbod
  !   Fp  => rbbou(iimbo) % pforce
  !   do idime = 1,ndime
  !      Fp(idime) = Fp(idime) - grnor * gravi(idime) * rbbou(iimbo) % volum * densi_nsi(1,1)
  !   end do
  !end do
  !
  ! Deallocate
  !
  if( INOTMASTER ) then
     if( imeth == 1 .and. kfl_intfo_nsi == 0 ) then
        call memgen(2_ip,ntens,npoin)
     else if( kfl_intfo_nsi >= 1 ) then
        call memgen(3_ip,npoin,0_ip)
     end if
  end if
  !
  ! Reduce sum in Parallel
  !
  call ibmdef(7_ip) 

  do iimbo = 1,nrbod

     Fv => rbbou(iimbo) % vforce
     Fp => rbbou(iimbo) % pforce
     Tv => rbbou(iimbo) % vtorqu
     Tp => rbbou(iimbo) % ptorqu
     vpfor => rbbou(iimbo) % vpfor
     vptor => rbbou(iimbo) % vptor

     do idime = 1,3
        vpfor(idime,1)  = Fv(idime) + Fp(idime)
        vptor(idime,1)  = Tv(idime) + Tp(idime)
        !        print *,"DEBUG: RB FORCES IN NASTIN VPFOR: ", vpfor(idime,1)

     end do

     ! print*,"DEBUG: fuerza en nastin ", vpfor

  end do




end subroutine nsi_rbobou
