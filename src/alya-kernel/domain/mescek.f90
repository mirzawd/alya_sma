!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mescek(itask)
  !-----------------------------------------------------------------------
  !****f* domain/mescek
  ! NAME
  !    mescek
  ! DESCRIPTION
  !    This routine checks if the mesh is correct
  ! OUTPUT
  !    VODOM ... Total domain volume
  !    VOAVE ... Averaged domain volume
  !    VOMIN ... Minimum element volume
  !    VOMAX ... Maximum element volume
  !    ELMIN ... Element with minimum volume
  !    ELMAX ... Element with maximum volume
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_domain
  use def_master
  use def_kermod
  use mod_memchk
  use mod_communications_global, only : PAR_SUM
  use mod_communications_global, only : PAR_MIN
  use mod_communications_global, only : PAR_MAX
  use mod_communications_global, only : PAR_BROADCAST
  use mod_outfor,                only : outfor
  use mod_messages,              only : messages_live
  use mod_iofile,                only : iofile_flush_unit
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_deallo
  use mod_meshes,                only : meshes_check_mesh
  use mod_mesh_type,             only : mesh_type_update_last_mesh
  use mod_bouder,                only : bouder
  use def_elmgeo,                only : element_type
  use mod_strings,               only : integer_to_string
  use mod_strings,               only : string_continue
  use mod_communications,        only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_live_info_config,      only : live_info_config
  use mod_memory,                only : memory_size
  use, intrinsic :: IEEE_ARITHMETIC, only : ieee_is_finite, ieee_is_nan

  implicit none

  integer(ip), intent(in)       :: itask
  integer(ip)                   :: ielem,inode,pgaus,pnode,pelty,igaus
  integer(ip)                   :: iboun,inodb,pgaub,pnodb,pblty,igaub
  integer(ip)                   :: i,j,k, ifiel
  integer(ip)                   :: kstar,kzero,keror,ieror,ipoin,npoit,idime
  integer(ip)                   :: neror(7),jpoin,iblty,nelet,pmate
  integer(ip)                   :: kpoin
  integer(ip)                   :: kfl_paral_min,kfl_paral_max,ipoin1,ipoin2
  real(rp)                      :: detjm,volum,surfa
  real(rp)                      :: cartd(ndime,mnode)
  real(rp)                      :: elcod(ndime,mnode),baloc(ndime,ndime)
  real(rp)                      :: xjaci(ndime,ndime),xjacm(ndime,ndime)
  character(10)                 :: mess1,mess2
  logical(lg),      pointer     :: touch(:)
  character(len=:), allocatable :: message
  integer(ip)                   :: ierro_element_wrong_node
  integer(ip)                   :: ierro_floating_node
  integer(ip)                   :: ierro_boundary_wrong_node
  integer(ip)                   :: ierro_local_node
  integer(ip)                   :: ierro_lelbo
  integer(ip)                   :: isize
  integer(ip),      pointer     :: touch_nodes(:)
  integer(ip)                   :: ierro_all(5)
  logical(lg)                   :: ifoun
#if defined(PNODE_VALUE) || defined(PGAUS_VALUE)
  integer(ip)                   :: ierro,ieaux
#endif

  neror = 0
  nullify(touch)

  select case ( itask )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Element Jacobian
     !
     !-------------------------------------------------------------------

     call messages_live('CHECK ELEMENT ORDERING')
     if( associated(vomat) ) then
        call memory_deallo(memor_dom,'VOMAT','mescek' ,vomat)
     end if
     call memory_alloca(memor_dom,'VOMAT','mescek' ,vomat,nmate)

     if( INOTMASTER ) then
        !
        ! Check element ordering by computing Jacobian sign
        !
        vomin    = huge(1.0_rp)
        vomax    = 0.0_rp
        vodom    = 0.0_rp
        elmax    = 1
        do ielem = 1,nelem
           pelty = ltype(ielem)
           if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              pmate = lmate(ielem)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              volum = 0.0_rp

              if( pelty == SHELL ) then
                 !
                 ! SHELL
                 !
                 call bouder(&
                      pnode,ndime,ndimb,elmar(pelty) % dercg,&
                      elcod,baloc,detjm)
                 volum = volum + elmar(pelty) % weicg * detjm
              else if( pelty == BAR3D ) then
                 !
                 ! BAR3D
                 !
                 detjm  = 1.0_rp
                 ipoin1 = lnods(1,ielem)
                 ipoin2 = lnods(2,ielem)
                 volum  = volum + sqrt(dot_product(coord(:,ipoin1)-coord(:,ipoin2),coord(:,ipoin1)-coord(:,ipoin2)))
              else
                 !
                 ! Other elements
                 !
                 do igaus = 1,pgaus
                    call jacobi(&
                         ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                         xjacm,xjaci,cartd,detjm)
                    volum = volum + elmar(pelty)%weigp(igaus)*detjm
                 end do
                 !call jacobi(&
                 !     ndime,pnode,elcod,elmar(pelty) % dercg,&
                 !     xjacm,xjaci,cartd,detjm)
              end if

              if( detjm <= 0.0_rp .and. lelch(ielem) /= ELINT ) then
                 neror(7) = neror(7)+1
                 mess1    = intost(1_ip)
                 mess2    = intost(ielem)
                 !block
                 !  use mod_strings
                 !  call meshe(ndivi) % output (filename='mesh-jacobian-'//integer_to_string(kfl_paral))
                 !end block
                 !print*,'subdo=',kfl_paral
                 !print*,'lnods=',lnods(:,ielem)
                 !do inode = 1,lnnod(ielem)
                 !   ipoin = lnods(inode,ielem)
                 !   print*,'coord=',ipoin,coord(:,ipoin)
                 !end do
                 call outfor(-2_ip,live_info_config%lun_livei,&
                      'JACOBIAN AT GAUSS POINT '//trim(mess1)&
                      //' OF ELEMENT '//trim(mess2)//' IS ZERO OR NEGATIVE')
                 call runend('ERROR IN MESH')
              end if

              vomat(pmate) = vomat(pmate) + volum   ! Volume per material        
              vodom        = vodom + volum          ! Total volume
              if( volum >= vomax ) then
                 vomax = volum                      ! Maximum volume
                 elmax = ielem
              end if
              if( volum <= vomin) then
                 vomin = volum                      ! Minimum volume
                 elmin = ielem
              end if
           end if
        end do
     end if

     if( INOTMASTER ) nelet = nelem
     call PAR_SUM(neror(7))
     call PAR_SUM(nelet)
     call PAR_SUM(vodom)
     call PAR_SUM(vomat)
     call PAR_MIN(vomin,'IN MY CODE',kfl_paral_min)
     call PAR_MAX(vomax,'IN MY CODE',kfl_paral_max)
     voave = vodom/real(nelet,rp)               ! Averaged volume

     if( INOTMASTER ) then
        if( kfl_paral_min /= kfl_paral ) then
           elmin = 0
        else
           elmin = leinv_loc(elmin)
        end if
        if( kfl_paral_max /= kfl_paral ) then
           elmax = 0
        else
           elmax = leinv_loc(elmax)
        end if
     end if
     call PAR_MAX(elmin)
     call PAR_MAX(elmax)

     if( neror(7) /= 0 ) call runend('MESCEK: ELEMENT(S) WITH NEGATIVE JACOBIAN')

     call messages_live('CHECK BOUNDARY ORDERING')

     !-------------------------------------------------------------------
     !
     ! Boundary Jacobian
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        !
        ! Check boundary ordering by computing Jacobian sign
        !
        surfa = 0.0_rp
        do iboun = 1,nboun
           pblty = ltypb(iboun)
           ielem = lelbo(iboun)
           if( pblty > 0 .and. pblty /= POINT ) then
              pnodb = nnode(pblty)
              pgaub = ngaus(pblty)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 do idime = 1,ndime
                    elcod(idime,inodb) = coord(idime,ipoin)
                 end do
              end do
              do igaub = 1,pgaub
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                      elcod,baloc,detjm)                                   ! and Jacobian
                 !call bouder(&
                 !     pnodb,ndime,ndimb,elmar(pblty)%dercg,&               ! Cartesian derivative
                 !     elcod,baloc,detjm)                                   ! and Jacobian
                 surfa = surfa + elmar(pblty)%weigp(igaub) * detjm

                 if( detjm <= 0.0_rp .and. lelch(ielem) /= ELINT ) then
                    neror(7) = neror(7) + 1
                    mess1    = intost(1_ip)
                    mess2    = intost(iboun)
                    call outfor(-2_ip,live_info_config%lun_livei,&
                         'JACOBIAN AT GAUSS POINT '//trim(mess1)&
                         //' OF BOUNDARY '//trim(mess2)//' IS ZERO OR NEGATIVE')
                 end if
              end do
           end if
        end do

     end if

     call PAR_SUM(neror(7))
     if( neror(7) /= 0 ) call runend('MESCEK: BOUNDARY(IES) WITH NEGATIVE JACOBIAN')

     if( INOTMASTER ) then

        if( kfl_chege == 1 ) then
           !
           ! Checks against two identical nonzero nodal coordinates
           !
           npoit=npoin-1
           npoit=0   ! >>> Not activated. Too much work
           do ipoin=1,npoit
              do jpoin=ipoin+1,npoin
                 keror=0
                 do idime=1,ndime
                    if(coord(idime,ipoin)==coord(idime,jpoin)) keror=keror+1
                 end do
                 if (keror==ndime) then
                    neror(1)=neror(1)+1
                    mess1=intost(ipoin)
                    mess2=intost(jpoin)
                    call outfor(1_ip,lun_outpu,&
                         ' IDENTICAL COORDINATES HAVE BEEN FOUND FOR POINTS NUMBER '&
                         //trim(mess1)//','//trim(mess2))
                 end if
              end do
           end do
           !
           ! Check if all the nodes belong to an element
           !
           call messages_live('CHECK IF ALL NODES BELONG TO AN ELEMENT')
           call memory_alloca(memor_dom,'TOUCH','mescek' ,touch,npoin)
           touch=.false.
           do ielem=1,nelem
              do inode = 1,nnode(ltype(ielem))
                 ipoin=lnods(inode,ielem)
                 touch(ipoin)=.true.
              end do
           end do
           do ipoin=1,npoin
              if(.not.touch(ipoin)) then
                 mess1=intost(ipoin)
                 call outfor(1_ip,lun_outpu,'NODE '//trim(mess1)&
                      //' DOES NOT BELONG TO ANY ELEMENT')
                 neror(6)=neror(6)+1
              end if
           end do
           call memory_deallo(memor_dom,'TOUCH','mescek' ,touch)
           !
           ! Checks for any repetition of a node number within an element
           !
           npoit=npoin
           npoit=0   ! >>> Not activated. Too much work
           do ipoin=1,npoit
              !
              ! Seek first,last and intermediate appearances of node ipoin
              ! & calculate increase or decrease in frontwidth at each element stage
              !
              kstar=0
              do ielem=1,nelem
                 kzero=0
                 do inode=1,lnnod(ielem)
                    if(lnods(inode,ielem)==ipoin) then
                       kzero=kzero+1
                       if(kzero>1) then
                          neror(3)=neror(3)+1
                          mess1=intost(ipoin)
                          mess2=intost(ielem)
                          call outfor(1_ip,lun_outpu,&
                               'NODE '//trim(mess1)&
                               //' APPEARS MORE THAN ONCE IN THE LIST OF'&
                               // ' NODAL CONNECTIONS OF ELEMENT NUMBER '//trim(mess2))
                       end if
                       if(kstar==0) kstar=ielem
                    end if
                 end do
              end do
              if(kstar==0) then
                 !
                 ! Checks if this is an unused node & if it has non-zero coordinates
                 !
                 mess1=intost(ipoin)
                 call outfor(1_ip,lun_outpu,&
                      'NODE '//trim(mess1)//' NEVER APPEARS IN THE ELEMENT CONNECTIVITY')
                 neror(4)=neror(4)+1
              end if
           end do
        end if
     end if

  case ( 2 )
     !
     ! Let us check the mesh just in case!
     !
     if( kfl_chege /= 0 ) then
        call messages_live('CHECKING THE MESH CONNECTIVITY JUST IN CASE...')

        nullify(touch_nodes)

        call memory_alloca(memor_dom,'TOUCH_NODES','meshes_check_element_connectivity',touch_nodes,npoin)

        ierro_element_wrong_node  = 0
        ierro_floating_node       = 0
        ierro_boundary_wrong_node = 0
        ierro_local_node          = 0
        ierro_lelbo               = 0
        !
        ! Types
        !
        !
        ! Element connectivity 
        !
        do ielem = 1,nelem
           do inode = 1,element_type(abs(ltype(ielem))) % number_nodes 
              ipoin = lnods(inode,ielem)
              if( ipoin < 1 .or. ipoin > npoin ) then
                 ierro_element_wrong_node = leinv_loc(ielem)
              else
                 touch_nodes(ipoin) = 1 
              end if
           end do
        end do
        !
        ! Boundary connectivity
        !
        loop_iboun: do iboun = 1,nboun
           ielem = lelbo(iboun)
           if( ielem < 1 .or. ielem > nelem ) then
              ierro_lelbo = lbinv_loc(iboun)
              exit loop_iboun
           end if
           pnode = element_type(abs(ltype(ielem))) % number_nodes
           do inodb = 1,element_type(abs(ltypb(iboun))) % number_nodes
              ipoin = lnodb(inodb,iboun)
              if( ipoin == 0 ) then
                 ierro_boundary_wrong_node = 1
              else
                 ifoun = .false.
                 loop_inode: do inode = 1,pnode
                    if( lnods(inode,ielem) == ipoin ) then
                       ifoun = .true.
                       exit loop_inode
                    end if
                 end do loop_inode
                 if( .not. ifoun ) then
                    ierro_local_node = lbinv_loc(iboun)
                 end if
              end if
           end do
        end do loop_iboun

        call PAR_INTERFACE_NODE_EXCHANGE(touch_nodes,'SUM')

        do ipoin = 1,npoin
           if( touch_nodes(ipoin) == 0 ) then
              ierro_floating_node = lninv_loc(ipoin)
           end if
        end do

        call memory_deallo(memor_dom,'TOUCH_NODES','meshes_check_element_connectivity',touch_nodes)

        ierro_all(1) = ierro_element_wrong_node
        ierro_all(2) = ierro_floating_node
        ierro_all(3) = ierro_boundary_wrong_node 
        ierro_all(4) = ierro_local_node   
        ierro_all(5) = ierro_lelbo
        isize        = size(ierro_all,KIND=ip)

        call PAR_MAX(isize,ierro_all)

        ierro_element_wrong_node  = ierro_all(1)
        ierro_floating_node       = ierro_all(2)
        ierro_boundary_wrong_node = ierro_all(3) 
        ierro_local_node          = ierro_all(4) 
        ierro_lelbo               = ierro_all(5) 

        if( ierro_floating_node       /= 0 ) call string_continue(message,'; ','THERE ARE NODES WITHOUT ELEMENTS, CHECK NODE '//integer_to_string(ierro_floating_node))
        if( ierro_local_node          /= 0 ) call string_continue(message,'; ','WRONG BOUNDARY CONNECTIVITY')
        if( ierro_element_wrong_node  /= 0 ) call string_continue(message,'; ','SOME ELEMENTS HAVE NODES OUT OF RANGE, CHECK ELEMENT '//integer_to_string(ierro_element_wrong_node))
        if( ierro_boundary_wrong_node /= 0 ) call string_continue(message,'; ','SOME BOUNDARY ELEMENTS HAVE NULL NODES')
        if( ierro_lelbo               /= 0 ) call string_continue(message,'; ','WRONG BOUNDARY/ELEMENT CONNECTIVITY FOR BOUNDARY '//integer_to_string(ierro_lelbo))


        do ifiel=1,nfiel
          do k=1,memory_size( xfiel(ifiel) % a, 3_ip )
            do j=1,memory_size( xfiel(ifiel) % a, 2_ip )
               do i=1,memory_size( xfiel(ifiel) % a, 1_ip )
                  if (.not. ieee_is_finite(xfiel(ifiel) % a(i,j,k) ) ) &
                     call string_continue(message,'; ',"FIELD "//trim(intost(ifiel))//"("//&
                        trim(intost(i))//","//trim(intost(j))//","//trim(intost(k))//&
                        ") IS INFINITE" )
                  if ( ieee_is_nan(xfiel(ifiel) % a(i,j,k) ) ) &
                     call string_continue(message,'; ',"FIELD "//trim(intost(ifiel))//"("//&
                        trim(intost(i))//","//trim(intost(j))//","//trim(intost(k))//&
                        ") IS NAN" )
               end do
            end do
          end do
        end do


        if( allocated(message) ) call runend('CHECKING MESH: '//message)        
     end if

  case ( 8 )
     !
     ! Some checks according to optimization options
     !
#ifdef PNODE_VALUE
     ierro = 0
     if( IMASTER ) then
        do pelty = element_num_ini(ndime),min(ubound(nnode,1_ip),element_num_end(ndime))
           pnode = nnode(pelty)
           if( lexis(pelty) /= 0 .and. pnode /= int(PNODE_VALUE,ip) ) then
              ierro = 1
              ieaux = pnode
           end if
        end do
     end if
     call PAR_BROADCAST(ierro)
     if( ierro /= 0 ) call runend('MESCEK: -DPNODE_VALUE OPTIONS NOT COMPATIBLE WITH YOUR MESH, FOUND PNODE= '//trim(intost(ieaux)))
#endif
#ifdef PGAUS_VALUE
     ierro = 0
     if( IMASTER ) then
        do pelty = element_num_ini(ndime),min(ubound(ngaus,1_ip),element_num_end(ndime))
           pgaus = ngaus(pelty)
           if( lexis(pelty) /= 0 .and. pgaus /= int(PGAUS_VALUE,ip) ) then
              ierro = 1
              ieaux = pgaus
           end if
        end do
     end if
     call PAR_BROADCAST(ierro)
     if( ierro /= 0 ) call runend('MESCEK: -DPGAUS_VALUE OPTIONS NOT COMPATIBLE WITH YOUR MESH, FOUND PGAUS= '//trim(intost(ieaux)))
#endif
     return

  end select
  !
  ! Verifies if any mess errors have been detected
  !
  keror=0
  do ieror=1,7
     keror=keror+neror(ieror)
  end do
  if(keror/=0) then
     mess1=intost(keror)
     if(keror==1) then
        call runend('1 ERROR IN MESH DEFINITION HAS BEEN DETECTED')
     else
        call runend(trim(mess1)//' ERRORS IN MESH DEFINITION HAVE BEEN DETECTED')
     end if
  end if

end subroutine mescek

