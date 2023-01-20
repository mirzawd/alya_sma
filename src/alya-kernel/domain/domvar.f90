!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    domvar.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Define some domain variables
!> @details Define some domain variables
!>          \verbatim
!>          LPOIZ(IZONE) % L(KPOIN) ....... =  0 if IPOIN does not belong to zone IZONE
!>                                  ....... /= 0 if IPOIN does
!>          LNOCH(NPOIN) .................. Node characterstics, compute from the element
!>                                          characteristics.
!>                                          = NOFEM ... Node is a normal node
!>                                          = NOHOL ... Node is a hole node
!>          LMAST(NPOIN) .................. List of masters
!>                                          LMAST(JPOIN) = IPOIN if JPOIN is a slave and IPOIN master
!>                                          = 0 if JPOIN is not a slave
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

subroutine domvar(itask)

  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_kermod,              only : kfl_posdi
  use def_kermod,              only : kfl_lapla
  use def_kermod,              only : kfl_rotation_axe,ndivi
  use def_kermod,              only : rotation_angle
  use def_kermod,              only : rotation_axis
  use def_kermod,              only : kfl_difun
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_memory,              only : memory_size
  use mod_parall,              only : mcode
  use mod_parall,              only : msubd
  use mod_parall,              only : mzone
  use mod_parall,              only : mcolo
  use mod_communications,      only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,      only : PAR_GHOST_BOUNDARY_EXCHANGE
  use mod_communications,      only : PAR_SUM
  use mod_communications,      only : PAR_MAX
  use mod_communications,      only : PAR_BROADCAST
  use mod_elmgeo,              only : elmgeo_element_type_initialization
  use mod_elmgeo,              only : elmgeo_number_nodes
  use mod_reaset,              only : reaset_set_connectivity
  use mod_maths_basic,         only : maths_rotation_matrix_2D
  use mod_maths_basic,         only : maths_rotation_matrix_3D
  use def_domain,              only : htable_lninv_loc
  use mod_materials,           only : materials_from_boundaries
  use mod_htable,              only : HtableMaxPrimeNumber
  use mod_htable,              only : hash_t
  use mod_htable,              only : htaini
  use mod_htable,              only : htaadd
  use mod_htable,              only : htades
  use mod_messages,            only : livinf
  use mod_domain,              only : domain_memory_allocate
  use mod_messages,            only : messages_live
  use mod_strings,             only : integer_to_string
  use def_elmgeo,              only : element_type
  use mod_boundary_conditions, only : boundary_conditions_number_boundary_codes
  use mod_mpio_config,         only : mpio_config
  use mod_std

  implicit none

  integer(ip),            intent(in) :: itask
  integer(ip)                        :: ielem,pelty,iboun,pblty,dummi
  integer(ip)                        :: inode,ipoin,jdime,idime,isubd
  integer(ip)                        :: ihole_aux,inodb,iperi,pnode
  integer(ip)                        :: ipoin_master,ipoin_slave,ifiel
  integer(ip)                        :: i1,i2,i3,n1,n2,n3
  integer(ip),            pointer    :: list_boundary_nodes(:)
  character(20)                      :: chnod
  real(rp)                           :: rotation_matrix(3,3),theta,xx(3)

  nullify(list_boundary_nodes)

  select case ( itask )

  case ( 0_ip )

     !-------------------------------------------------------------------
     !
     ! Materials, subomains, cohesive elements, number of groups
     ! NGROU_DOM should be broadcast if the master has found a number
     ! of groups different from that prescribed
     !
     !-------------------------------------------------------------------

     nmate     = 0
     nsubd     = 0
     nperi     = 0
     kfl_elcoh = 0
     call materials_from_boundaries()
     if( INOTMASTER .and. nelem > 0 ) then
        nmate     = maxval(lmate)
        nsubd     = maxval(lesub)
        nperi     = maxval(lmast)
        kfl_elcoh = min(1_ip,count(lelch == ELCOH,KIND=ip))
        kfl_elint = min(1_ip,count(lelch == ELINT,KIND=ip))
     end if
     call PAR_MAX      (nmate)
     call PAR_MAX      (nsubd)
     call PAR_MAX      (nperi)
     call PAR_MAX      (kfl_elcoh)
     call PAR_MAX      (kfl_elint)
     call PAR_BROADCAST(ngrou_dom)
     if( nelem == 0 ) then
        IEMPTY    = .true.
        INOTEMPTY = .false.
     else
        IEMPTY    = .false.
        INOTEMPTY = .true.
     end if
     mzone     = nzone
     msubd     = nsubd
     call PAR_MAX(mzone,'IN THE WORLD')
     call PAR_MAX(msubd,'IN THE WORLD')
     mcolo = (mcode+1)*(mzone+1)*(msubd+1)-1
     ncodb = boundary_conditions_number_boundary_codes()
     
     !-------------------------------------------------------------------
     !
     ! If some meshes are moving, activate kfl_domar
     !
     !-------------------------------------------------------------------

     if( kfl_modul(ID_ALEFOR) /= 0 .or. kfl_difun /= 0 ) then
        kfl_domar       = 1
        kfl_domar_world = 1
     end if
     call PAR_MAX(kfl_domar)
     call PAR_MAX(kfl_domar_world,'IN THE WORLD')     

  case ( 1_ip )
     !
     ! If we export the mesh, desactivate mesh multiplication
     !
     if( mpio_config%output%post_process%export_only ) then
        ndivi            = 0
        kfl_posdi        = 0
        trans(1)         = 0.0_rp
        trans(2)         = 0.0_rp
        trans(3)         = 0.0_rp
        xscal(1)         = 1.0_rp
        xscal(2)         = 1.0_rp
        xscal(3)         = 1.0_rp
        kfl_rotation_axe = 0
     end if

     !-------------------------------------------------------------------
     !
     ! Mesh dependent arrays
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        !
        ! LNLEV, LELEV, LBLEV: level of mesh multiplication
        !
        call memory_alloca(memor_dom,'LNLEV','domvar',lnlev,npoin,'IDENTITY')
        call memory_alloca(memor_dom,'LELEV','domvar',lelev,nelem,'IDENTITY')
        call memory_alloca(memor_dom,'LBLEV','domvar',lblev,nboun,'IDENTITY')
        !
        ! COORD: translation
        !
        if( trans(1) /= 0.0_rp .or. trans(2) /= 0.0_rp .or. trans(3) /= 0.0_rp ) then
           do ipoin = 1,npoin
              coord(1:ndime,ipoin) = trans(1:ndime) + coord(1:ndime,ipoin)
           end do
        end if
        !
        ! COORD: scale factor
        !
        if( xscal(1) /= 1.0_rp .or. xscal(2) /= 1.0_rp .or. xscal(3) /= 1.0_rp ) then
           do ipoin = 1,npoin
              coord(1:ndime,ipoin) = xscal(1:ndime) * coord(1:ndime,ipoin)
           end do
        end if
        !
        ! COORD: rotation
        !
        if( kfl_rotation_axe /= 0 ) then
           rotation_matrix = 0.0_rp
           theta           = rotation_angle/180.0_rp*pi
           if( ndime == 2 ) then
              rotation_matrix(1:2,1:2) = maths_rotation_matrix_2D(theta)
           else
              rotation_matrix = maths_rotation_matrix_3D(theta,v=rotation_axis)
           end if
           do ipoin = 1,npoin
              xx(1:ndime)          = coord(1:ndime,ipoin)
              coord(1:ndime,ipoin) = 0.0_rp
              do idime = 1,ndime
                 do jdime = 1,ndime
                    coord(idime,ipoin) = coord(idime,ipoin) + rotation_matrix(idime,jdime) * xx(jdime)
                 end do
              end do
           end do
        end if
        !
        ! Fields scaling
        !
        if( kfl_xscal_fields /= 0 ) then
           do ifiel = 1,nfiel
              n1 = memory_size(xfiel(ifiel) % a,1_ip)
              n2 = memory_size(xfiel(ifiel) % a,2_ip)
              n3 = memory_size(xfiel(ifiel) % a,3_ip)
              do i3 = 1,n3
                 do i2 = 1,n2
                    do i1 = 1,n1
                       xfiel(ifiel) % a(i1,i2,i3) = xscal_fields(ifiel) * xfiel(ifiel) % a(i1,i2,i3)  
                    end do
                 end do
              end do
           end do
        end if
     end if

  case ( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Mesh dependent variables: this subroutine must be called whenever
     ! new elements are created (for example after mesh multiplication)
     !
     !-------------------------------------------------------------------

     !
     ! Number of elements per type
     !
     lnuty = 0
     do ielem = 1,nelem
        pelty = abs(ltype(ielem))
        lnuty(pelty) = lnuty(pelty) + 1
     end do
     !
     ! Number of boundaries per type
     !
     do iboun = 1,nboun
        pblty = abs(ltypb(iboun))
        lnuty(pblty) = lnuty(pblty) + 1
     end do
     !
     ! Check if high order element exist
     !
     do ielem = 1,nelem
        pelty = abs(ltype(ielem))
        if( pelty <= element_end ) then
           if( lorde(pelty) > 1 ) kfl_horde = 1
        end if
     end do
     !
     ! Contact/hole/cohesive elements/interface elements
     !
     do ielem = 1,nelem
        if( lelch(ielem) == ELCNT ) then
           ltype(ielem) = -abs(ltype(ielem))
        else if( lelch(ielem) == ELHOL ) then
           ltype(ielem) = -abs(ltype(ielem))
        end if
     end do
     call PAR_MAX(kfl_horde,'IN MY CODE')
     !
     ! Copy of original mesh size for geometrical halo dimensions
     !
     if( ISEQUEN ) then
        npoi1        = npoin
        npoi2        = npoin+1
        npoi3        = npoin
        npoin_par(1) = npoin
        nelem_par(1) = nelem
        nboun_par(1) = nboun
        npoin_total  = npoin
        nelem_total  = nelem
        nboun_total  = nboun
     end if
     nelem_2    = nelem
     nboun_2    = nboun
     npoin_2    = npoin
     npoin_own  = npoi3
     npoin_halo = npoin_own
     !
     ! Fields dimensions
     !
     do ifiel = 1,nfiel
        if( kfl_field(1,ifiel) /= 0 ) then
           if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
              kfl_field(5,ifiel) = nelem
           else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
              kfl_field(5,ifiel) = npoin
           else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
              kfl_field(5,ifiel) = nboun
           else
              call runend('OUTDOM: UNDEFINED FIELD TYPE')
           end if
        end if
     end do
     !
     ! Cancel Laplacian if required
     !
     if( kfl_lapla == 0 ) llapl = 0

  case ( 3_ip )

     !-------------------------------------------------------------------
     !
     ! Must be called after renelm() and mesh multiplication
     !
     !-------------------------------------------------------------------
     !
     ! Htable for LNINV_LOC
     !
     if( INOTEMPTY ) then
        call htades( htable_lninv_loc, memor_opt=memor_dom  )
        if( nperi == -2 ) then
           call htaini( htable_lninv_loc, npoin, lidson=.true., AUTOMATIC_SIZE=.true.,memor_opt=memor_dom,REPEATED_ELEMENTS=.true.)
        else
           call htaini( htable_lninv_loc, npoin, lidson=.true., AUTOMATIC_SIZE=.true.,memor_opt=memor_dom)
        end if
        call htaadd( htable_lninv_loc, lninv_loc, memor_opt=memor_dom )
        if(htable_lninv_loc % nelem /= npoin ) then
           call runend('DOMVAR: TROUBLES WITH LNINV_LOC')
        end if
     end if
     !
     ! Set connectivity
     !
     call reaset_set_connectivity()
     !
     ! Define LNOCH according to LELCH and detect holes automatically
     !
     if( INOTMASTER .and. npoin > 0 ) then
        ihole_aux= 0
        call memgen(1_ip,npoin,0_ip)
        do ielem = 1,nelem
           pnode = elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
           if( lelch(ielem) /= ELHOL ) then
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 gisca(ipoin) = 1
              end do
           else
              ihole_aux= ihole_aux+1
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
        do ipoin = 1,npoin
           if( gisca(ipoin) == 0 ) then
              lnoch(ipoin) = NOHOL
              if (ihole_aux == 0) then
                 chnod = intost(ipoin)
                 call livinf(10000_ip, &
                      'NODE NUMBER ( '//adjustl(trim(chnod)) &
                      // ' HAS NO CONNECTIVITY DEFINED)',0_ip)
              end if
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
     !
     ! List of boundary nodes
     ! NBONO .......... Number of boundary nodes
     ! LBONO(NBONO) ... List of boudanry nodes
     !
     if( INOTMASTER .and. npoin > 0 ) then
        call memory_alloca(memor_dom,'LIST_BOUNDARY_NODES','memgeo',list_boundary_nodes,npoin)
        nbono = 0
        do ipoin = 1,npoin
           list_boundary_nodes(ipoin) = 0
        end do
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              if ( ipoin>npoin ) then
                 call runend( "Domvar.f90: Array index of bounds when trying to access point ID " &
                      //trim(intost(ipoin))//" for boundary element "//trim(intost(iboun))// &
                      ". Max point ID "//trim(intost(npoin))//&
                      ". LNODB likely was not read correctly, see BOUNDARIES..END_BOUNDARIES section in dom.dat for errors." )
              end if
              if( list_boundary_nodes(ipoin) == 0 ) then
                 nbono = nbono + 1
                 list_boundary_nodes(ipoin) = 1
              end if
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(list_boundary_nodes,'MAX','IN MY CODE')
        nbono = 0
        do ipoin = 1,npoin
           nbono = nbono + list_boundary_nodes(ipoin)
        end do
        call domain_memory_allocate('LBONO')
        nbono = 0
        do ipoin = 1,npoin
           if( list_boundary_nodes(ipoin) /= 0 ) then
              nbono = nbono + 1
              lbono(nbono) = ipoin
           end if
        end do
        call memory_deallo(memor_dom,'LIST_BOUNDARY_NODES','memgeo',list_boundary_nodes)
     end if
     !
     ! List of periodicity couples using LMAST
     ! LPERI(1,:) = master
     ! LPERI(2,:) = slave
     !
     ! LMAST(JPOIN) = IPOIN, IPOIN is the master of slave JPOIN
     !
     if( INOTEMPTY .and. new_periodicity == 0 ) then
        nperi = 0
        nperi = count( lmast(1:npoin) /= 0 ,KIND=ip)
        if( nperi > 0 ) then
           call domain_memory_allocate('LPERI')
           iperi = 0
           do ipoin_slave = 1,npoin
              ipoin_master = lmast(ipoin_slave)
              if( ipoin_master > 0 ) then
                 iperi = iperi + 1
                 lperi(1,iperi) = ipoin_master
                 lperi(2,iperi) = ipoin_slave
              end if
           end do
        end if
     end if
     !
     ! Check LESUB
     !
     call memgen(1_ip,npoin,0_ip)
     dummi = 0
     do ielem = 1,nelem
        isubd = lesub(ielem)
        do inode = 1,element_type(abs(ltype(ielem))) % number_nodes
           ipoin = lnods(inode,ielem)
           if( gisca(ipoin) == 0 .or. gisca(ipoin) == isubd ) then
              gisca(ipoin) = lesub(ielem)
           else
              dummi = lninv_loc(ipoin)
           end if
        end do
     end do
     call memgen(3_ip,npoin,0_ip)
     call PAR_MAX(dummi)
     if( dummi > 0 ) call runend('DOMVAR: NODE '//integer_to_string(dummi)//' BELONGS TO TWO SUBDOMAINS')

  case ( 4_ip )
     !
     ! LNNOD: number of nodes per element. Treat special case of virtual elements
     !
     call domain_memory_allocate('LNNOD')
     do ielem = 1,nelem
        lnnod(ielem) = elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
     end do
     meshe(ndivi) % lnnod => lnnod
     !
     ! LNNOB: Number of nodes per boundary
     !
     call domain_memory_allocate('LNNOB')
     do iboun = 1,nboun
        lnnob(iboun) = nnode(abs(ltypb(iboun)))
     end do
     meshe(ndivi) % lnnob => lnnob
     !
     ! LGAUS: Number of Gauss points per element
     !
     call domain_memory_allocate('LGAUS')
     do ielem = 1,nelem
        pelty = abs(ltype(ielem))
        if( pelty <= element_end ) lgaus(ielem) = ngaus(pelty)
     end do
     meshe(ndivi) % lgaus => lgaus
     !
     ! Empty subdomain
     !
     if( nelem == 0 ) then
        IEMPTY    = .true.
        INOTEMPTY = .false.
     else
        IEMPTY    = .false.
        INOTEMPTY = .true.
     end if

  end select

end subroutine domvar

