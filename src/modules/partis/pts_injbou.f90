!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_injbou()

  use def_kintyp,        only : ip,rp
  use def_master
  use def_domain
  use def_partis
  use mod_memory
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE 
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_pts_injection,  only : PTS_INJ_GEO_CIRCLE    
  use mod_bouder

  implicit none
  integer(ip)          :: iinj,ipoin,ielem,inode
  integer(ip)          :: iboun,inodb,pnodb,igaub,pelty,pgaub,pgaus,pblty,pnode
  integer(ip)          :: nnodes_corona
  integer(ip), pointer :: codno_pts(:),aux(:),lnodes_corona(:)
  real(rp)             :: center_gravity(ndime) , normal_surface(ndime), surface, radius, dista(ndime)
  real(rp)             :: gbsur,gbcoo(ndime),eucta
  real(rp)             :: bocod(ndime,mnodb),elcod(ndime,mnode),baloc(ndime,ndime)

  nullify(codno_pts)
  nullify(aux)
  nullify(lnodes_corona)

  do iinj = 1,pts_minj
           

     if( codbo_pts(iinj) /= 0 ) then
        nnodes_corona  = 0_ip
        center_gravity = 0.0_rp
        normal_surface = 0.0_rp
        surface        = 0.0_rp
        radius         = 0.0_rp

        !
        ! Definir la corona
        !
        if( INOTEMPTY ) then
           call memory_alloca(mem_modul(1:2,modul),'CODNO_PTS','pts_injbou',codno_pts,npoin)               
           call memory_alloca(mem_modul(1:2,modul),'AUX',      'pts_injbou',aux,npoin)               
           codno_pts(1:npoin) = 0_ip
           aux(1:npoin)       = 0_ip
           !
           ! Mark boundary nodes which do not belong to our target
           !
           do iboun = 1,nboun
              if( kfl_codbo(iboun) /= codbo_pts(iinj) ) then
                 pnodb = nnode(ltypb(iboun))
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    codno_pts(ipoin) = kfl_codbo(iboun)
                 end do
              end if
           end do

           call PAR_INTERFACE_NODE_EXCHANGE(codno_pts,'MAX')

           do iboun = 1,nboun
              if( kfl_codbo(iboun) == codbo_pts(iinj) ) then
                 pnodb = nnode(ltypb(iboun))
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    if(codno_pts(ipoin)== 0) then
                       codno_pts(ipoin) = kfl_codbo(iboun)
                    else if (codno_pts(ipoin) == kfl_codbo(iboun)) then
                       codno_pts(ipoin) = kfl_codbo(iboun)
                    else 
                       !nodo de corona
                       nnodes_corona = nnodes_corona + 1_ip
                       aux(ipoin)    = ipoin
                    end if
                 end do
              end if
           end do

           call memory_alloca(mem_modul(1:2,modul),'LNODES_CORONA','pts_injbou',lnodes_corona,nnodes_corona)
           
           nnodes_corona = 0_ip
           do ipoin = 1,npoin
              if( aux(ipoin) /= 0 ) then
                 nnodes_corona = nnodes_corona + 1_ip
                 lnodes_corona(nnodes_corona) = aux(ipoin)
              end if
           end do
           
           call memory_deallo(mem_modul(1:2,modul),'CODNO_PTS','pts_injbou',codno_pts)
           call memory_deallo(mem_modul(1:2,modul),'AUX',      'pts_injbou',aux)
           !
           ! Definir el centro de gravedad de la superficie total
           !
           do iboun = 1,nboun
              pblty = ltypb(iboun) 
              if( pblty > 0 .and. kfl_codbo(iboun) == codbo_pts(iinj) ) then
                 pnodb = lnnob(iboun)
                 pgaub = ngaus(pblty)                 
                 ielem = lelbo(iboun)
                 pelty = ltype(ielem)
                 if( pelty > 0 ) then
                    pnode = nnode(pelty)
                    pgaus = ngaus(pelty)
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                    end do
                 end if
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
                 end do

                 gauss_points: do igaub = 1,pgaub
                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                         bocod,baloc,eucta)                                   ! and Jacobian
                    gbsur = elmar(pblty)%weigp(igaub)*eucta
                    gbcoo = 0.0_rp
                    do inodb = 1,pnodb
                       gbcoo(1:ndime) = gbcoo(1:ndime) &
                            + elmar(pblty)%shape(inodb,igaub) &
                            * bocod(1:ndime,inodb)
                    end do
                    call chenor(pnode,baloc,bocod,elcod)
                    
                    center_gravity(1:ndime) = center_gravity(1:ndime) + gbcoo(1:ndime)       * gbsur
                    normal_surface(1:ndime) = normal_surface(1:ndime) + baloc(1:ndime,ndime) * gbsur
                    surface                 = surface                 + gbsur
                    
                 end do gauss_points

              end if
           end do

        end if
        !
        ! Parallelization
        !
        call PAR_SUM(ndime,center_gravity)
        call PAR_SUM(ndime,normal_surface)
        call PAR_SUM(surface)
        !
        ! Radius
        !        
        do ipoin = 1,nnodes_corona
           dista(1:ndime) = center_gravity(1:ndime) - coord(1:ndime,lnodes_corona(nnodes_corona))
           radius = max(0.0_rp,sqrt(dot_product(dista,dista)))
        end do
        !
        ! Parallelization
        !
        call PAR_MAX(radius)

        if( surface /= 0.0_rp ) then
           center_gravity = center_gravity / surface
           normal_surface = normal_surface / surface
        else
           call runend('PTS_INJBOU: NULL SURFACE')
        end if

        injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_CIRCLE
        injection_pts(iinj) % geo_coord    = center_gravity
        injection_pts(iinj) % geo_rad      = radius
        injection_pts(iinj) % geo_normal   = normal_surface
     end if
  end do
  
end subroutine pts_injbou
