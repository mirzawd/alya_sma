!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_outvar(ivari,imesh)
  !------------------------------------------------------------------------
  !****f* Master/ker_output
  ! NAME
  !    ker_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    ker_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_elmtyp
  use def_kermod
  use def_domain
  use def_coupli
  use mod_ker_proper
  use mod_ker_vortex
  use mod_memory,         only : memory_size
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_parall,         only : par_omp_num_colors
  use mod_parall,         only : par_omp_ia_colors
  use mod_parall,         only : par_omp_ja_colors
  use mod_parall,         only : par_omp_num_threads
  use mod_parall,         only : par_omp_nelem_chunk
  use mod_parall,         only : commd
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_SUM
  use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
  use mod_renumbering,    only : renumbering_elements
  use mod_graphs,         only : graphs_eleele
  use mod_graphs,         only : graphs_dealep
  use mod_gradie,         only : gradie
  use mod_projec,         only : projec_elements_to_nodes
  use mod_projec,         only : projec_boundaries_to_nodes
  use mod_ker_nsw_visc2,  only : ker_nod_nsw_visc_0
  use mod_outvar,         only : outvar
  use mod_par_tools,      only : par_tools_comm_to_array
  use mod_arrays,         only : arrays
  use mod_arrays,         only : arrays_name
  use mod_arrays,         only : arrays_number
  use mod_postpr

  use mod_matrix,         only : matrix_CSR_parallel_SpMV
  use mod_matrices,       only : matrices_laplacian
  use mod_biofibers,      only : biofib_point_nodal_fibers

  !!$ use mod_AMR
  !!$ use def_AMR

  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ibopo,ipoin,idime,ielpo,isubd,nenti,kdime
  integer(ip)             :: ielem,icont,iline,jpoin,iboun,inode,ithre
  integer(ip)             :: dummi,inodb,imate,incnt,kpoin,istack,kmesh
  integer(ip)             :: nstack,izdom,icoup,icolo,kelem,ii,ifiel,ienti
  integer(ip)             :: igaub,pblty,pgaub,igaus,jdime,icoun
  real(rp)                :: rutim,dummr,rmate
  character(5)            :: wopos(3)
  integer(ip), pointer    :: lstack(:)
  integer(ip), pointer    :: list_colors(:)
  real(rp),    pointer    :: velbo_ker(:,:,:)
  character(20)           :: wfiel

  real(rp), pointer :: lapla(:)
  real(rp), pointer :: aux(:)
  !!$ real(rp)          :: auxr(ndime)
  !!$ integer(ip)       :: pnode

!  real(rp), pointer :: hh_opt_elem(:), hh_opt_node(:)

  nullify(lapla)
  nullify(aux)

  !!$ nullify(hh_opt_elem)
  !!$ nullify(hh_opt_node)

  if( ivari == 0 ) return
  !
  ! Define postprocess variable
  !
  rutim = cutim
  
  select case ( arrays_name(ivari) )  

  case ( 'EXNOR' )
     !
     ! EXNOR
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = exnor(idime,1,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 'ERRNO' )
     !
     ! ERRNO
     !
     if( INOTEMPTY ) then
        allocate(lapla(nzdom))
        call matrices_laplacian(lapla)
        call memgen(0_ip,npoin,0_ip)
        allocate(aux(npoin))
        do ipoin = 1,npoin
           aux(ipoin) = sqrt(dot_product(veloc(1:ndime,ipoin,1),veloc(1:ndime,ipoin,1)))
        end do
        call matrix_CSR_parallel_SpMV(npoin,1_ip,1_ip,r_dom,c_dom,lapla,aux,gesca)
        deallocate(aux,lapla)
     end if

     !!$if (INOTMASTER) then

     !!$    call memgen(0_ip,npoin,0_ip)

     !!$    allocate(hh_opt_elem(nelem),hh_opt_node(npoin))

     !!$    do ielem = 1,nelem
     !!$       auxr = 0.0_rp
     !!$       pnode     = lnnod(ielem)
     !!$       do inode = 1,pnode
     !!$          ipoin                = lnods(inode,ielem)
     !!$          auxr(1:ndime)        = auxr(1:ndime) + ((0.1_rp/real(pnode)) * coord(1:ndime,ipoin))
     !!$       end do
     !!$       auxr(1:ndime) = amrp0(1:ndime)-auxr(1:ndime)
     !!$       hh_opt_elem(ielem) = sqrt(dot_product(auxr,auxr))-0.45_rp
     !!$    end do
     !!$    do ielem = 1,nelem
     !!$       pnode     = lnnod(ielem)
     !!$       do inode = 1,pnode
     !!$          ipoin                = lnods(inode,ielem)
     !!$          hh_opt_node(ipoin)  = hh_opt_elem(ielem)
     !!$       end do
     !!$    end do
     !!$    gesca = hh_opt_node

     !!$    deallocate(hh_opt_elem,hh_opt_node)
     !!$end if

     !!$if (associated(hh_opt))then
     !!$   if ( INOTMASTER ) then
     !!$       call memgen(0_ip,npoin,0_ip)
     !!$       allocate(hh_opt_node(npoin))
     !!$       !print*, 'hh_opt_elem := ',hh_opt
     !!$       do ielem = 1,nelem
     !!$          pnode     = lnnod(ielem)
     !!$          do inode = 1,pnode
     !!$             ipoin                = lnods(inode,ielem)
     !!$             hh_opt_node(ipoin)  = hh_opt(ielem)
     !!$          end do
     !!$       end do
     !!$       gesca = hh_opt_node

     !!$       deallocate(hh_opt_node)
     !!$   end if
     !!$end if

  case ( 'LPOIN' )
     !
     ! Geometrical local basis and type of point: LPOIN, SKCOS
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,1,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 'SKCO1' )
     !
     ! Geometrical local basis SKCOS(:,1,:)
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,1,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 'SKCO2' )
     !
     ! Geometrical local basis SKCOS(:,2,:)
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,2,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 'SKCO3' )
     !
     ! Geometrical local basis SKCOS(:,NDIME,:)
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              do idime = 1,ndime
                 gevec(idime,ipoin) = skcos(idime,ndime,ibopo)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case ( 'DISPM' )
     !
     ! DISPM: Mesh displacement
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        if( associated(dispm) ) then
           do ipoin = 1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin) = dispm(idime,ipoin,1)
              end do
           end do
        else
           do ipoin = 1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end do
        end if
        !if( nimbo > 0 ) then
        !   kpoin = npoin
        !   do iimbo = 1,nimbo
        !      do ipoib = 1,imbou(iimbo)%npoib
        !         kpoin = kpoin + 1
        !         do idime = 1,ndime
        !            gevec(idime,kpoin) = imbou(iimbo)%cooi2(idime,ipoib) - imbou(iimbo)%cooin(idime,ipoib)
        !         end do
        !      end do
        !   end do
        !end if
     end if
     !if( kfl_outib == 3 ) kfl_outib = 4

  case ( 'LPOTY' )
     !
     ! Boundary points: LPOTY
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           if(lpoty(ipoin)/=0) gesca(ipoin)=1.0_rp
        end do
     end if

  case ( 'PONUM' )
     !
     ! PONUM: Point Numbering
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lninv_loc(ipoin),rp)
        end do
     end if

  case ( 'CODNO' )
     !
     ! Boundary codes: KFL_CODNO
     !
     if( kfl_icodn > 0 ) then
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime = 1,ndime
                 gevec(idime,ipoin)=real(kfl_codno(idime,ipoin),rp)
              end do
           end do
        end if
     else
        return
     end if

  case ( 'YWALP' )
     !
     ! YWALP: Wall distance
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo /= 0 ) gesca(ipoin)=ywalp(ibopo)
        end do
     end if

  case ( 'LNTIB' )
     !
     ! LNTIB: node type
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin)=real(lntib(ipoin),rp)
        end do
     end if

  case ( 'HOLES' )
     !
     ! HOLES: with node type
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           if( lntib(ipoin) > 0 ) then
              gesca(ipoin) = 1.0_rp
           end if
        end do
     end if

  case ( 'DENSI' )
     !
     ! DENSI: Density
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        densi_ker % kfl_nedsm = 1
        call ker_proper('DENSI','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'VISCO' )
     !
     ! VISCO: viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('VISCO','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'POROS' )
     !
     ! POROS : porosity
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        poros_ker % kfl_nedsm = 1
        call ker_proper('POROS','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'CONDU' )
     !
     ! CONDU: Conductivity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('CONDU','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'SPECI' )
     !
     ! SPECI: Specific Heat
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('SPHEA','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'GROUP' )
     !
     ! LGROU_DOM
     !
     if( ngrou_dom <= 0 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lgrou_dom(ipoin),rp)
        end do
     end if

  case ( 'MASSM' )
     !
     ! MASSM
     !
     gesca => vmass

  case ( 'MASSC' )
     !
     ! MASSC
     !
     gesca => vmasc

  case ( 'GEONO' )
     !
     ! KFL_GEONO
     !
     if( kfl_geome == 0 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo /= 0 ) gesca(ipoin)=real(kfl_geono(ibopo),rp)
        end do
     end if

  case ( 'SUBDO' )
     !
     ! Subdomains
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoi1
           gesca(ipoin) = real(kfl_paral,rp)
        end do
        do ipoin = npoi1+1,npoin
           gesca(ipoin) = -1.0_rp
        end do
        ! tuve que descomentariar estas lineas sino daba mal
        ! porque estan comentadas   - habra que agregar if iparall????
        do ipoin = npoi2,npoi3
           gesca(ipoin) = real(kfl_paral,rp)
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 'WALLD' )
     !
     ! WALLD
     !
     if( kfl_walld == 0 ) return
     gesca => walld

  case ( 'ROUGH' )
     !
     ! ROUGH
     !
     if( kfl_rough < 1 ) return
     gesca => rough

  case( 'KEKET' )
     !
     ! LINEL: Linelets of preconditioner CG
     !
     if( INOTMASTER ) then
        icont = 0
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        do iline = 1,solve(2) % nline
           icont = icont+1
           do ipoin = solve(2) % lline(iline),solve(2) % lline(iline+1)-1
              jpoin = solve(2) % lrenup(ipoin)
              rhsid(jpoin) = real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case ( 'CODBO' )
     !
     ! Boundary codes: CODBO
     !
     if( INOTMASTER ) then
        call memgen(zero,nboun,zero)
        do iboun = 1,nboun
           gesca(iboun)=real(kfl_codbo(iboun),rp)
        end do
     end if

  case ( 'MATER' )
     !
     ! MATER: Material Numbering (Elemental)
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(lmate(ielem),rp)
        end do
     end if

  case ( 'LMATN' )
     !
     ! Nodal material (maximum)
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do imate = 1,nmate
           rmate = real(imate,rp)
           do ii = 1,memory_size(lmatn(imate)%l)
              ipoin = lmatn(imate)%l(ii)
              gesca(ipoin) = max(rmate,gesca(ipoin))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')
     end if

  case ( 'COMMU' )
     !
     ! COMMU: communication arrays
     !
     if( ISEQUEN ) return
     call par_tools_comm_to_array(PAR_COMM_MY_CODE_ARRAY(1),npoin,givec,'GIVEC')
     ii = memory_size(givec,1_ip)
     call memgen(zero,ii,npoin)

     do ipoin = 1,npoin
        gevec(:,ipoin) = real(givec(:,ipoin),rp)
     end do
     call memgen(3_ip,ii,npoin)
     
  case ( 'LELEV' )
     !
     ! LELEV
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(lelev(ielem),rp)
        end do
     end if

  case ( 'ELNUM' )
     !
     ! ELNUM: Global Element Number
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(leinv_loc(ielem),rp)
        end do
     end if

  case ( 'CODBB' )
     !
     ! Boundary codes: CODBB
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do iboun = 1,nboun
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              gesca(ipoin)=real(kfl_codbo(iboun),rp)
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')
     end if

  case ( 'LETIB' )
     !
     ! LETIB
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = max(real(letib(ielem),rp),gesca(ipoin))
           end do
        end do
     end if

  case ( 'LTYP2' )
     !
     ! LTYPE
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = real(ltype(ielem),rp)
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MIN','IN MY CODE')
     end if

  case ( 'LELC2' )
     !
     ! LELCH
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              dummr = real(lelch(ielem),rp)
              gesca(ipoin) = max(gesca(ipoin),dummr)
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 'CONTA' )
     !
     ! Contact
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do incnt = 1,nncnt
           ipoin = lncnt(1,incnt)
           jpoin = lncnt(2,incnt)
           gesca(ipoin) = 1.0_rp
           gesca(jpoin) = 2.0_rp
        end do
     end if

  case ( 'RDOM ' )
     !
     ! R_DOM
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(r_dom(ipoin+1)-r_dom(ipoin),rp)
        end do
     end if

  case ( 'VORTX' )
     !
     ! VORTX: extraction of vortex core
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
     else
        call memgen(zero,1_ip,1_ip)
     end if
     call ker_vortex(gevec)
     if( IMASTER ) call memgen(two,1_ip,1_ip)

  case ( 'LESET' )
     !
     ! Element sets
     !
     if( neset < 1 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = max(gesca(ipoin),real(leset(ielem),rp))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')
     end if

  case ( 'DISMM' )
     !
     ! DISPL_KER
     !
     if( kfl_suppo == 0 ) return
     gevec => displ_ker

  case ( 'VELOC' )
     !
     ! VELOC
     !
     if( INOTMASTER ) then
        if( associated(veloc) ) then
           gevec => veloc(:,:,1)
        else if( associated(advec) ) then
           gevec => advec(:,:,1)
        end if
     end if

  case ( 'LBSET' )
     !
     ! LBSET: boundary sets
     !
     if( nbset < 1 ) return
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        gesca = 0.0_rp
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              gesca(ipoin) = max(gesca(ipoin),real(lbset(iboun),rp))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 'LMESH' )
     !
     ! Connected meshes
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        allocate(lstack(npoin))
        do ipoin = 1,npoin
           lstack(ipoin) = 0
        end do

        kmesh        = 0
        kpoin        = 0

        do while( kpoin /= npoin )

           kmesh = kmesh + 1
           ipoin = 1
           do while( gesca(ipoin) /= 0.0_rp )
              ipoin = ipoin + 1
           end do
           nstack       = 1
           lstack(1)    = ipoin
           gesca(ipoin) = real(kmesh,rp)
           istack       = 0
           kpoin        = kpoin + 1

           do
              if( istack == nstack ) exit
              istack = istack + 1
              ipoin  = lstack(istack)
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( gesca(jpoin) == 0.0_rp ) then
                    gesca(jpoin)   = real(kmesh,rp)
                    nstack         = nstack + 1
                    lstack(nstack) = jpoin
                    kpoin          = kpoin + 1
                 end if
              end do
           end do
        end do
        deallocate(lstack)

     end if

  case ( 'TURBU' )
     !
     ! TURBU: Turbulent viscosity
     !
     if( INOTMASTER ) then
        turmu_ker % kfl_nedsm = 1
        call memgen(zero,npoin,zero)
        call ker_proper('TURBU','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'LNSUB' )
     !
     ! LNSUB: Element subdomain
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ielem = 1,nelem
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              gesca(ipoin) = max(gesca(ipoin),real(lesub(ielem),rp))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')
     end if

  case ( 'WETNO' )
     !
     ! Wet nodes
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do icoup = 1,mcoup
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              if( ipoin /= 0 ) gesca(ipoin) = real(icoup,rp)
           end do
        end do
     end if

  case ( 'LESUB' )
     !
     ! LESUB: Element Subdomains
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = max(gesca(ielem),real(lesub(ielem),rp))
        end do
     end if

  case ( 'CANOP' )
     !
     ! CANOPY HEIGHT
     !
     if( kfl_canhe < 1 ) return
     gesca => canhe

  case ( 'HEIGH' )
     !
     ! HEIGHT OVER TERRAIN
     !
     if( kfl_heiov < 1 ) return
     gesca => heiov

  case ( 'BEATR' )
     !
     ! BEATR: Element subdomain + characteristic + material
     !
     if( INOTMASTER ) then
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           icont = 100 * lmate(ielem) + 10 * lesub(ielem) + lelch(ielem)
           gesca(ielem) = max(gesca(ielem),real(icont,rp))
        end do
     end if

  case ( 'COLOR' )
     !
     ! COLORING of the openmp stratgy
     !
     if( INOTMASTER ) then
        allocate(list_colors(nelem))
        call memgen(zero,nelem,zero)
        do icolo = 1,par_omp_num_colors
           do kelem = par_omp_ia_colors(icolo),par_omp_ia_colors(icolo+1)-1
              ielem = par_omp_ja_colors(kelem)
              list_colors(ielem) = icolo
           end do
        end do
        do ielem = 1,nelem
           gesca(ielem) = real(list_colors(ielem),rp)
        end do
        !do ipoin = 1,npoin
        !   gesca(ipoin) = huge(1.0_rp)
        !end do
        !do ielem = 1,nelem
        !   do inode = 1,lnnod(ielem)
        !      ipoin = lnods(inode,ielem)
        !      gesca(ipoin) = min(gesca(ipoin),real(list_colors(ielem),rp))
        !   end do
        !end do
        !call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MIN','IN MY CODE')
        deallocate(list_colors)
     end if

  case ( 'WALLN' )
     !
     ! WALLN
     !
     gevec => walln

  case ( 'LOCAL' )
     !
     ! Local numbering
     ! what was here originally did not work well
     ! now I leave 2 option that both work well
     !
     if(1==1) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) = real(ipoin,rp)
           end do
        end if
     else if(1==2) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              ! esto da mal
              !gesca(ipoin) = real(ipoin,rp)
           end do
           do ipoin = 1,meshe(ndivi)%npoi1
              gesca(ipoin) = real(ipoin,rp)
           end do
           do ipoin = meshe(ndivi)%npoi1+1 , meshe(ndivi)%npoin
              gesca(ipoin) = real(PAR_COMM_MY_CODE_ARRAY(1) % node_number_in_owner(ipoin-meshe(ndivi)%npoi1),rp)
           end do
        end if
     else
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoi1
              gesca(ipoin) = real(ipoin,rp)
           end do
           do ipoin = npoi1+1,npoin
              gesca(ipoin) = -1.0_rp
           end do
           do ipoin = npoi2,npoi3
              gesca(ipoin) = real(ipoin,rp)
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
        end if
     end if

  case ( 'ADVEC' )
     !
     ! Advection
     !
     gevec => advec(:,:,1)

  case ( 'RENEL' )
     !
     ! RENEL: local element renumbering
     !
     call memgen(0_ip,nelem,zero)
     do ielem = 1,nelem
        gesca(ielem) = real(ielem,rp)
     end do

  case ( 'RENPO' )
     !
     ! RENPO: local node renumbering
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(ipoin,rp)
        end do
     end if

  case ( 'OMPSS' )
     !
     ! OMPSS element subdomains
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do isubd = 1, size(ompss_domains)
           do kelem = 1,size(ompss_domains(isubd)%elements)
              ielem = ompss_domains(isubd)%elements(kelem)
              gesca(ielem) = real(isubd,rp)
           end do
        end do
        !
        ! Mark separator
        !
        ielpo = 0
!!$           if( .not. associated(lelel) ) then
!!$              call graphs_eleele(&
!!$                   nelem,npoin,mnode,mepoi,lnods,lnnod,&
!!$                   pelpo,lelpo,nedg1,medg1,pelel,lelel)
!!$           end if
!!$           do ielem = 1,nelem
!!$              do ielpo = pelel(ielem),pelel(ielem+1)-1
!!$                 kelem = lelel(ielpo)
!!$                 if( gesca(ielem) /= gesca(kelem) .and. gesca(ielem) > 0.0_rp .and. gesca(kelem) > 0.0_rp ) then
!!$                    if( gesca(ielem) < gesca(kelem) ) then
!!$                       gesca(ielem) = -abs(gesca(ielem))
!!$                    else
!!$                       gesca(kelem) = -abs(gesca(kelem))
!!$                    end if
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$           ipoin = 1
!!$           dummi = 0
!!$           do while( ipoin > 0 )
!!$              dummi = dummi + 1
!!$              print*,'passs=',dummi
!!$              ipoin = 0
!!$              do ielem = 1,nelem
!!$                 if( gesca(ielem) < 0.0_rp ) then
!!$                    inode = 0
!!$                    do ielpo = pelel(ielem),pelel(ielem+1)-1
!!$                       kelem = lelel(ielpo)
!!$                       if( abs(gesca(ielem)) /= gesca(kelem) .or. gesca(kelem) > 0.0_rp ) then
!!$                          inode = 1
!!$                       end if
!!$                    end do
!!$                    if( inode == 0 ) then
!!$                       ipoin = ipoin + 1
!!$                       gesca(ielem) = abs(gesca(ielem))
!!$                    end if
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$           gesca = max(gesca,0.0_rp)

        if( ielpo /= 0 ) then
           call graphs_dealep(pelel,lelel,memor=memor_dom)
        end if

     end if

  case ( 'OPENM' )
     !
     ! OPENMP element chunks
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        ielem = 0
        element_loop: do
           do ithre = 1,par_omp_num_threads
              do kelem = 1,par_omp_nelem_chunk
                 ielem = ielem + 1
                 if( ielem > nelem ) then
                    exit element_loop
                 else
                    gesca(ielem) = real(ithre,rp)
                 end if
              end do
           end do
        end do element_loop
     end if

  case ( 'GAUSS' )
     !
     ! LGAUS: number of Gauss points
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,zero)
        do ielem = 1,nelem
           gesca(ielem) = real(lgaus(ielem),rp)
        end do
     end if

  case ( 'WALLO' )
     !
     ! WALLD
     !
     if( kfl_walld == 2 .or. kfl_walld == 3) then
        do ipoin = 1,npoin
           gesca(ipoin) = real(wallo(ipoin),rp)
        end do
     endif

  case ( 'MIXIN' )
     !
     ! VISCA: air viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('MIXIN','NPOIN',dummi,dummi,gesca)
     end if

  case ( 'FIELD' )
     !
     ! FIELD
     !
     ! only if no prealoading is enabled
     do ifiel = 1,nfiel

        if (kfl_field(6,ifiel) /= 1 ) then !not on demand


            wfiel =  intost(ifiel)
            if(ifiel<10) then
               wopos(1)=postp(1)%wopos(1,70)(1:3)//'0'//trim(wfiel(1:2))
            else
               wopos(1)=postp(1)%wopos(1,70)(1:3)//trim(wfiel(1:2))
            end if
            if( kfl_field(1,ifiel) > 0 ) then
               if(      kfl_field(2,ifiel) == NPOIN_TYPE ) then
                  nenti    = npoin
                  wopos(3) = 'NPOIN'
               else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                  nenti    = nboun
                  wopos(3) = 'NBOUN'
               else if( kfl_field(2,ifiel) == NELEM_TYPE ) then
                  nenti    = nelem
                  wopos(3) = 'NELEM'
               end if
               kdime = kfl_field(1,ifiel)
               if(      kdime == 1     ) then
                  wopos(2) = 'SCALA'
               else if( kdime == ndime ) then
                  wopos(2) = 'VECTO'
               else
                  kdime = 0
               end if
               
               if( kdime == 1 ) then
                  if( INOTMASTER ) then
                     call memgen(zero,nenti,0_ip)
                     do ienti = 1,nenti
                        gesca(ienti) = xfiel(ifiel) % a(1,ienti,1)
                     end do
                  end if
                  call postpr(gesca,wopos,ittim,rutim)
                  if( INOTMASTER ) call memgen(two,nenti,0_ip)
                  
               else if( kdime == ndime ) then
                  
                  if( INOTMASTER ) then
                     call memgen(zero,ndime,nenti)
                     do ienti = 1,nenti
                        gevec(1:ndime,ienti) = xfiel(ifiel) % a(1:ndime,ienti,1)
                     end do
                  end if
                  call postpr(gevec,wopos,ittim,rutim)
                  if( INOTMASTER ) call memgen(two,ndime,nenti)
                  
               end if
               
            end if
        end if !not on demand
     end do
     
     return

  case ( 'VELOM' )
     !
     ! VELOM
     !
     if( INOTMASTER ) then 
        call memgen(zero,ndime,npoin)
        if( associated(velom) ) then
           gevec(1:ndime,1:npoin) = velom(1:ndime,1:npoin)
        else
           gevec = 0.0_rp
        end if
     end if
     
  case ( 'OMPSB' )
     !
     ! OMPSS boundary subdomains
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nboun,zero)
        do isubd = 1, size(ompss_boundaries)
           do kelem = 1,size(ompss_boundaries(isubd)%elements)
              ielem = ompss_boundaries(isubd)%elements(kelem)
              gesca(ielem) = real(isubd,rp)
           end do
        end do
     end if
     
  case ( 'NSWVI' )
     !
     ! NSWVI: nodal projection of no slip wall  viscosity
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        if( kfl_noslw_ker /= 0 ) then
           call ker_nod_nsw_visc_0
           gesca = 0.0_rp
           call projec_elements_to_nodes(el_nsw_visc,gesca)
        end if
     end if

  case ( 'NUMBE' )
     !
     ! NUMBERING
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(ipoin,rp)
        end do
     end if
     
  case ( 'VELAV' )
     !
     ! VELAV
     !
     call arrays(arrays_number('VELAV'),'POSTPROCESS',velav_ker,MESH_ID=imesh)
     return
     
  case ( 'MPIRA' )
     !
     ! MPI RANK
     !
    if( INOTMASTER ) then
       call memgen(0_ip,nelem,zero)
       do ielem = 1,nelem
          gesca(ielem) = real(kfl_paral,rp)
       end do
    end if

  case ( 'VELAE' )
    !
    ! VELAE: VELAV on elements
    !
    call memgen(0_ip,ndime,nelem)
    do iboun = 1,nboun
       ielem = lelbo(iboun)
       gevec(1:ndime,ielem) = velav_ker(1:ndime,1,iboun)
    end do

  case ( 'LAD  ' )
     !
     ! CANOPY LAD
     !
     if( kfl_canla < 1 ) return
     gesca => canla
      
  case ( 'TEMPE' )
    !
    ! TEMPER
    !
    gesca => tempe(:,1)

  case ( 'BFIBL' )
      !
      ! BFIBL: Bio-fibers longitudinal-fiber direction
      !
     call biofib_point_nodal_fibers( gevec, 'LONGITUDINAL', 'CURRENT' )

  case ( 'BFIBS' )
     !
     ! BFIBS: Bio-fibers transversal-Sheet direction
     !
     call biofib_point_nodal_fibers( gevec, 'SHEET', 'CURRENT' )

  case ( 'BFIBN' )
     !
     ! BFIBN: Bio-fibers transversal-Nheet direction
     !
     call biofib_point_nodal_fibers( gevec, 'NORMAL', 'CURRENT' )
    
  case ( 'INTNO' )
     !
     ! Boundary
     !
     call memgen(0_ip,npoin,0_ip)
     do ii = 1,commd % bound_dim
        ipoin = commd % bound_perm(ii)
        if( ipoin <= npoin_own ) & 
             gesca(ipoin) = 1.0_rp
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM')
    
  case ( 'VELNO' )
     !
     ! VELNO_KER
     !
     ii = 0
     if( memory_size(lexlo_ker) > 0 ) ii = 1
     call PAR_SUM(ii)
     if( ii == 0 ) return
     
     nullify(velbo_ker)
     call memory_alloca(mem_modul(1:2,modul),'VELBO_KER','ker_memall',velbo_ker,ndime,mgaub,nboun)
     call memgen(0_ip,ndime,npoin)
     do iboun = 1,nboun
        pblty = abs(ltypb(iboun))
        pgaub = ngaus(pblty)                                                     ! Number of Gauss points
        do igaub = 1,pgaub
           if( lexlo_ker(igaub,iboun) /= 0 ) then
              velbo_ker(1:ndime,igaub,iboun) = velel_ker(1:ndime,lexlo_ker(igaub,iboun))
           end if
        end do
     end do
     call projec_boundaries_to_nodes(velbo_ker,meshe(ndivi),gevec)
     call memory_deallo(mem_modul(1:2,modul),'VELBO_KER','ker_memall',velbo_ker)
    
  case ( 'FRING' )
     !
     ! FRING: Fringe nodes in wetnode structure
     !
     if( INOTMASTER ) then
        do icoup = 1,mcoup
           if ( coupling_type(icoup) % code_target == current_code ) then
              call memgen(zero,npoin,zero)
              if ( coupling_type(icoup) % wet % kfl_get_fringe ) then
                 do ipoin = 1,npoin
                    gesca(ipoin) = real(coupling_type(icoup) % wet % kfl_fringe_wetnodes(ipoin),rp)
                 end do
              else
                 do ipoin = 1,npoin
                    gesca(ipoin) = 0_ip
                 end do
              end if
           end if
        end do
     end if
    
  case ( 'DAVID' )
     !
     ! DAVID: Characteristic nodes ### debugging ###
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lnoch(ipoin),rp)
        end do
     end if

  case ( 'GRWAL' )
     !
     ! Walld gradient
     !
     call memgen(0_ip,ndime,npoin)
     call gradie(walld,gevec)

  case ( 'ANIPO' )
     !
     ! ANIPO
     ! 
     call memgen(0_ip,ndime*ndime,nelem)
     do ielem = 1,nelem
        igaus = 1
        icoun = 0
        do jdime = 1,ndime
           do idime = 1,ndime
              icoun = icoun + 1
              ii    = (idime-1)*ndime+jdime
              gevec(icoun,ielem) = anipo_ker % value_ielem(ielem) % a(ii)
           end do
        end do
     end do
     
  case ( 'COMMI' )
     !
     ! COMMI: inverse communication arrays
     !
     if( ISEQUEN ) return
     call par_tools_comm_to_array(PAR_COMM_MY_CODE_ARRAY(1),npoin,givec,'GIVEC',GATHER=.false.)
     ii = memory_size(givec,1_ip)
     call memgen(zero,ii,npoin)

     do ipoin = 1,npoin
        gevec(:,ipoin) = real(givec(:,ipoin),rp)
     end do
     call memgen(3_ip,ii,npoin)
     

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine ker_outvar
