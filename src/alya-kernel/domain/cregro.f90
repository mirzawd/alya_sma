!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine cregro()
  !-----------------------------------------------------------------------
  !****f* domain/cregro
  ! NAME
  !    cregro
  ! DESCRIPTION
  !    This routine computes the groups for the deflated CG.
  !    Starting nodes are those node ipoin imposed, that is when 
  !    LIMPO(ipoin)>0
  !
  !    TO DO 1: RENUMBER GROUPS TO MINIMIZE SLYKINE
  !    TO DO 2: PERFORM A GROUP SMOOTHING USING A CG ONTO GROUPS TO SOLVE
  !             -Lapl(g) + r g = g' r >> 1 in the iterior
  !                                 r -> 0 near the boundary
  !
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use def_elmtyp
  use mod_postpr
  use mod_couplings,      only : COU_PUT_VALUE_ON_TARGET
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_SUM
  use mod_par_tools,      only : PAR_GLOBAL_TO_LOCAL_NODE
  use mod_messages,       only : messages_live
  implicit none
  integer(ip)          :: ipoin,ibopo,nqueu,npogr,igrou,jgrou,iqueu
  integer(ip)          :: kgrou,kpoin,izdom,jpoin,npomi,ngrou
  integer(ip)          :: nfron,nfold,jqueu,ifront,nfnew,icheck
  integer(ip)          :: ndofn,ierro,nmarkt,nmark,inode
  integer(ip)          :: npoin_all,ngrou_par,ii,nimpo
  real(rp)             :: f,g,h
  integer(ip), pointer :: lngro(:)
  integer(ip), pointer :: lqueu(:)
  integer(ip), pointer :: lfron(:)
  integer(ip), pointer :: lgrou(:)
  integer(ip), pointer :: limpo(:)
  integer(ip), pointer :: gigro(:)
  logical(lg), pointer :: lmark(:)
  logical(lg)          :: lggro

  nullify(lngro)
  nullify(lqueu)
  nullify(lfron)
  nullify(lgrou)
  nullify(limpo)
  nullify(gigro)
  nullify(lmark)

  kpoin = 0
  ierro = 0
  lggro = .false.

  if( solve_sol(1) % ngrou /= 0 .or. ngrou_dom /= 0 ) lggro = .true.

  if( lggro ) then

     !----------------------------------------------------------------------
     !
     ! Allocate memory and impose boundary conditions
     !
     !----------------------------------------------------------------------

     call messages_live(namod(modul)//': COMPUTE GROUPS OF DEFLATION FOR '//trim(solve_sol(1) % wprob))
     
     if( ngrou_dom == 0 ) then
        npoin_all = 0
        if( ISEQUEN ) then
           npoin_all = npoin
        else if( INOTMASTER ) then
           npoin_all = npoi1 + (npoi3-npoi2+1)
        end if
        call PAR_SUM(npoin_all)
        !
        ! Automatic number of groups
        !
        if( solve_sol(1) % ngrou < 0 ) then
           f = real(npoin_all,rp) / 2000.0_rp + 10.0_rp
           g = 50.0_rp * log(real(npoin_all,rp))
           call mixing(1_ip,h,1.0e-6_rp,1.0e6_rp,real(npoin_all,rp))
           solve_sol(1) % ngrou = int( g*h+(1.0_rp-h)*f ,ip)
        end if
        if( ISEQUEN ) then
           ngrou_par = solve_sol(1) % ngrou
        else
           ngrou_par = solve_sol(1) % ngrou / npart
           ngrou_par = max(1_ip,ngrou_par)
        end if
     else
        solve_sol(1) % ngrou = ngrou_dom
     end if

     ierro = 0
     if( INOTMASTER ) then
        !
        ! LGROU: allocate memory
        !
        limpo => solve_sol(1) % limpo
        call memory_alloca(memit,'SOLVE % LGROU','cregro',solve_sol(1) % lgrou,npoin)
        lgrou => solve_sol(1) % lgrou
        !
        ! Assign group from LGROU_DOM
        !
        if( ngrou_dom /= 0 ) then
           do ipoin = 1,npoin
              solve_sol(1) % lgrou(ipoin) = lgrou_dom(ipoin)
           end do
        end if
     end if                                  

     if( INOTMASTER ) then

        if( IEMPTY ) then                           ! All nodes

           continue
           
        else if(  solve_sol(1) % ifbop == 1 ) then  ! Boundary nodes

           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo /= 0 ) then
                 if( limpo(ibopo) > 0 ) then
                    lgrou(ipoin) = -1
                 end if
              end if
           end do

        else if( .not. associated(limpo) ) then    ! Choose one node
          
           kpoin = huge(1_ip)
           lpoty_loop1: do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo /= 0 ) then
                 kpoin = ipoin
                 exit lpoty_loop1
              end if
           end do lpoty_loop1
           if( kpoin >=1 .and. kpoin <= npoin ) kpoin = lninv_loc(kpoin)
           call PAR_MIN(kpoin,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
           if( kpoin > npoin ) call runend('CREGRO: COULD NOT FIND NODE')
           ipoin = PAR_GLOBAL_TO_LOCAL_NODE(kpoin)
           if( ipoin > 0 ) lgrou(ipoin) = -1
           
        else if( memory_size(limpo) == npoin ) then       ! All nodes

           do ipoin = 1,npoin
              if( limpo(ipoin) > 0 ) then
                 lgrou(ipoin) = -1
              end if
           end do

        else if( memory_size(limpo) == 1 ) then          ! Only one node

           if( limpo(1) /= 0 ) then
              write(*,*) 'CREGRO: LIMPO NOT CODED'
           end if

           nodes: do ipoin = 1,npoin
              if( limpo(1) == ipoin ) then
                 lgrou(ipoin) = -1
                 exit nodes
              end if
           end do nodes

        end if
        !
        ! If there is a coupling, impose boundary condition
        !
        call COU_PUT_VALUE_ON_TARGET(-1_ip,lgrou)        
        !
        ! Impose periodic slaves
        !
        do ipoin = 1,npoin
           if(lmast(ipoin)/=0) lgrou(ipoin) = -2
        end do
        !
        ! Impose holes
        !
        do ipoin = 1,npoin
           if( lnoch(ipoin) == NOHOL ) lgrou(ipoin) = -2
        end do
     end if
     !
     ! Count imposed nodes
     !
     if( INOTMASTER .and. ngrou_dom == 0 ) then
        nmarkt = npoin
        do ipoin = 1,npoin
           if( lgrou(ipoin) == -1 ) nmarkt = nmarkt - 1
        end do       
        !
        ! Fill in LQUEU, LMARK, LFRON
        !
        call memory_alloca(memit,'LQUEU','cregro',lqueu,npoin)
        call memory_alloca(memit,'LMARK','cregro',lmark,npoin)
        call memory_alloca(memit,'LFRON','cregro',lfron,npoin)
        do ipoin = 1,npoin
           lmark(ipoin) = .false.
        end do
        !
        ! Compute points per group and minimal threshold
        !
        !npogr  = max(npoin/solve_sol(1) % ngrou,1_ip)
        npogr  = max(npoin/ngrou_par,1_ip)
        npomi  = npogr/3
        !
        ! Find first non marked point
        !
        if( solve_sol(1) % ifbop == 1 ) then

           loop1: do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo /= 0 ) then
                 if( limpo(ibopo) > 0 ) then
                    do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                       jpoin = c_dom(izdom)
                       if( lgrou(jpoin) == 0 ) then                   
                          kpoin = jpoin
                          exit loop1
                       endif
                    enddo
                 end if
              end if
           end do loop1

        else if( memory_size(limpo) == npoin ) then

           loop2: do ipoin = 1,npoin
              if( limpo(ipoin) > 0 ) then
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if( lgrou(jpoin) == 0 ) then                   
                       kpoin = jpoin
                       exit loop2
                    endif
                 enddo
              end if
           end do loop2

        else if( memory_size(limpo) == 1 ) then

           loop3: do ipoin = 1,npoin
              if( limpo(1) == ipoin ) then
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if( lgrou(jpoin) == 0 ) then                   
                       kpoin = jpoin
                       exit loop3
                    endif
                 enddo
              end if
           end do loop3

        end if
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Check errors
  !
  !----------------------------------------------------------------------

  !if( ISEQUEN .and. solve_sol(1) % ngrou /= 0 ) then    
  !   if( ngrou_par /= 0 ) then    
  !      if( kpoin == 0 ) ierro = 1
  !      call errors(2_ip,ierro,iwarn,'CREGRO: COULD NOT FIND IMPOSED NODE')
  !   end if
  !end if

  !----------------------------------------------------------------------
  !
  ! Decide who starts
  !
  !----------------------------------------------------------------------

  !npart = 2
  !if( IPARALL ) then
  !   allocate( lkpoi(npart+1) )
  !   lkpoi(kfl_paral+1) = kpoin
  !   call parari('SUM',1_ip,npart+1,lkpoi)
  !   ii = 2
  !   do while( ii <= npart )
  !      if( lkpoi(ii) > 0 ) then
  !         jj = ii -1
  !         ii = npart
  !      end if
  !      ii = ii + 1
  !   end do
  !   if( jj /= kfl_paral ) kpoin = 0
  !end if

  !if( ISLAVE ) then
  !   pard1 = kpoin
  !   call Parall(800_ip)
  !   kpoin = pard1
  !end if

  if( ngrou_dom == 0 ) then

     !----------------------------------------------------------------------
     !
     ! If not already found, look for a starting node
     !
     !----------------------------------------------------------------------

     if( ISLAVE .and. kpoin == 0 ) then
        ipoin = 0
        do while( ipoin < npoin )
           ipoin = ipoin + 1
           if( lpoty(ipoin) /= 0 .and. ipoin > npoi1 ) then
              kpoin = ipoin
              ipoin = npoin
           end if
        end do
        if( kpoin == 0 ) then
           ipoin = 0
           do while( ipoin < npoin )
              ipoin = ipoin + 1
              if( ipoin > npoi1 ) then
                 kpoin = ipoin
                 ipoin = npoin
              end if
           end do
        end if
     end if

     !----------------------------------------------------------------------
     !
     ! Construct groups
     !
     !----------------------------------------------------------------------

     igrou = 0

     if( lggro .and. INOTMASTER ) then
        !
        !  Initialize main loop
        !
        nmark        = 0
        igrou        = 1
        lqueu(1)     = kpoin
        lgrou(kpoin) = igrou
        nmark        = nmark + 1
        nqueu        = 1
        iqueu        = 1
        nfron        = 0

        open_loop: do

           !
           ! Loop on queue
           ! 

           neigh:do 
              ipoin = lqueu(iqueu)
              !
              ! Loop on neighbors
              ! 
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 !
                 ! Has the group of this point been assigned before ?
                 !
                 if( lgrou(jpoin) == 0 ) then
                    !
                    ! Did we reach the maximum number of points per group
                    !
                    if( nqueu == npogr ) exit neigh
                    lgrou(jpoin) = igrou
                    nmark        = nmark + 1    
                    nqueu        = nqueu + 1
                    lqueu(nqueu) = jpoin
                    !
                    ! Does this point belong to the current front
                    !
                    if( lmark(jpoin) .eqv. .false. ) lmark(jpoin) = .true.
                 endif
              end do
              !
              ! Did we exhaust the queue?
              !        
              if( iqueu == nqueu )exit neigh
              iqueu = iqueu+1

           end do neigh
           !
           ! Add nodes to the front
           !
           nfold = nfron+1

           do jqueu = iqueu,nqueu
              ipoin = lqueu(jqueu)
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( lgrou(jpoin) == 0 ) then
                    if( lmark(jpoin) .eqv. .false. ) then
                       nfron = nfron+1
                       lfron(nfron) = jpoin
                       lmark(jpoin) = .true.
                    end if
                 end if
              end do
           end do
           !
           ! Clean up the last points
           !
           do ifront = nfold,nfron
              ipoin = lfron(ifront)
              lmark(ipoin) = .false.
           enddo
           !
           ! Compress the front
           !
           nfnew = 0
           do ifront = 1,nfron
              ipoin = lfron(ifront)
              if( lmark(ipoin) .eqv. .false. ) then
                 nfnew = nfnew + 1
                 lfron(nfnew) = lfron(ifront)
              endif
           enddo
           nfron = nfnew
           !
           ! Special case: Do we have enough nodes   
           !
           if( nqueu < npomi ) then

              !
              !   Find a neighboring group
              !

              glue: do iqueu = 1,nqueu
                 ipoin = lqueu(iqueu)
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    jgrou = lgrou(jpoin)
                    if( jgrou /= igrou .and. jgrou /= -1 ) then
                       exit glue
                    endif
                 end do
              end do glue

              do ipoin = 1,npoin
                 kgrou = lgrou(ipoin)
                 if( kgrou == igrou ) then
                    lgrou(ipoin) = jgrou
                    !nmark        = nmark + 1    
                 endif
              enddo

              nqueu = 0
              igrou = igrou-1

           endif
           !
           ! Is the front empty --> go home
           !
           !
           !     Is the front empty --> go home
           !
           !if(nfron==0) exit open_loop
           if( nfron == 0 ) then
              !
              !     Did we mark all the points
              !
              if( nmark == nmarkt ) then
                 !
                 !     Is the front empty --> go home
                 !
                 exit open_loop

              else
                 !
                 !     Find a non marked point
                 !
                 icheck=0_ip
                 do ipoin=1,npoin
                    if(lgrou(ipoin)==0)then
                       lfron(1)=ipoin
                       icheck=1_ip
                       exit
                    endif
                 enddo
                 if(icheck==0)then
                    call runend('CREGRO: nmark/=nmarkt and point not found')
                 endif
              endif

           endif
           !
           ! Find new seed
           !
           jpoin=lfron(1)
           !
           ! Initialize new round
           !
           igrou        = igrou+1
           lqueu(1)     = jpoin
           lgrou(jpoin) = igrou
           nmark        = nmark + 1    
           nqueu        = 1
           iqueu        = 1 
           lmark(jpoin) = .true.

        end do open_loop

     end if
     !
     ! Redefine ngrou
     !
     ngrou_par = igrou
     ngrou     = ngrou_par
     !
     ! Parall
     !
     if( IPARALL ) then

        call memory_alloca(memit,'LNGRO','cregro',lngro,npart+1_ip)        
        lngro(kfl_paral+1) = ngrou
        call PAR_SUM(npart+1_ip,lngro)
        lngro(1) = 0
        ngrou    = 0
        do ii = 2,npart+1
           ngrou     = ngrou     + lngro(ii)
           lngro(ii) = lngro(ii) + lngro(ii-1)
        end do

        if( INOTMASTER ) then
           do ipoin = 1,npoin
              if( solve_sol(1) % lgrou(ipoin) > 0 ) then
                 solve_sol(1) % lgrou(ipoin) = solve_sol(1) % lgrou(ipoin) + lngro(kfl_paral)
              end if
           end do
        end if
        call memory_deallo(memit,'LNGRO','cregro',lngro)        

     end if

     solve_sol(1) % ngrou = ngrou
     ndofn = solve_sol(1) % ndofn

     if( ISLAVE ) then
        pari1 => solve_sol(1) % lgrou
        call par_cregrp(1_ip) ! Uniquely define goups on interface nodes
     end if

  else if( ngrou_dom /= 0 ) then
     !
     ! Groups are defined by LGROU_DOM
     !
     ngrou = solve_sol(1) % ngrou 
     ndofn = solve_sol(1) % ndofn

  end if

  if( INOTMASTER .and. ngrou_dom == 0 ) then
     !
     ! LQUEU, LMARK, LFRON: Deallocate memory
     !
     call memory_deallo(memit,'LNGRO','cregro',lqueu)    
     call memory_deallo(memit,'LMARK','cregro',lmark)    
     call memory_deallo(memit,'LFRON','cregro',lfron)    

  end if

  !----------------------------------------------------------------------
  !
  ! Check if no node has been imposed
  !
  !----------------------------------------------------------------------
  !
  ! NIMPO= # nodes imposed in my zone (apart from holes)
  !
  nimpo = 0
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        if( lgrou(ipoin) == -1 ) nimpo = nimpo + 1
     end do
  end if
  call PAR_SUM(nimpo)
  !
  ! Detect node to impose
  !
  if( nimpo == 0 ) then

     if( INOTMASTER ) then 
        inode = huge(ip)
        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 .and. lgrou(ipoin) > 0 ) then
              if( lninv_loc(ipoin) < inode ) inode = lninv_loc(ipoin)
           end if
        end do
     end if
     call PAR_MIN(inode)
     if( INOTMASTER ) then 
        jpoin = 0 
        ipoin = 0
        do while( ipoin < npoin )
           ipoin = ipoin + 1
           if( lninv_loc(ipoin) == inode ) then
              jpoin = ipoin
              ipoin = npoin
           end if
        end do
        if( jpoin > 0 ) then
           lgrou(jpoin) = -1 
        end if 
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Give '0' flag to prescribed (0) / hole (-2) / periodic (-2) nodes
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     do ipoin = 1,npoin
        igrou = lgrou(ipoin)  
        if( igrou <= -1 ) lgrou(ipoin) = 0
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Fix nodes, take off the empty groups, and renumber them 
  !
  !----------------------------------------------------------------------

  call memory_alloca(memit,'GIGRO','cregro',gigro,ngrou)

  if( INOTMASTER ) then
     do ipoin = 1,npoin
        igrou = lgrou(ipoin)  
        if( igrou > ngrou ) then
           print*,'PROBLEM IN CREGRO: ',kfl_paral,igrou,lninv_loc(ipoin)
           call runend('CREGRO: TOO FEW GROUPS!')
        end if
        if( igrou > 0 ) gigro(igrou) = gigro(igrou) + 1
     end do
  end if

  call PAR_SUM(ngrou,gigro)

  jgrou = 0
  do igrou = 1,ngrou
     if( gigro(igrou) > 0 ) then
        jgrou = jgrou + 1
        gigro(igrou) = jgrou
     end if
  end do

  if( INOTMASTER ) then
     do ipoin = 1,npoin
        igrou = lgrou(ipoin)
        if( igrou > 0 ) lgrou(ipoin) = gigro(igrou)
     end do
  end if

  ngrou = jgrou
  solve_sol(1) % ngrou = ngrou
  call memory_deallo(memit,'GIGRO','cregro',gigro)

  !----------------------------------------------------------------------
  !
  ! Compute matrix assembly format sparse/skyline
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_defas == 0 ) then
     call skygro()
  else if( solve_sol(1) % kfl_defas == 1 ) then
     call csrgro()
  else if( solve_sol(1) % kfl_defas == 2 ) then 
     call runend('CREGRO: ASSEMBLY NOT CODED')
  end if

  !----------------------------------------------------------------------
  !
  ! Parallelization
  !
  !----------------------------------------------------------------------

  if( IPARALL ) then
     !
     ! Compute parallelization strategy (only if allgatherv is used)
     !
     if( solve_sol(1) % kfl_gathe == 1 ) call par_gatgro()


  end if

end subroutine cregro

subroutine possla(itask,dummi,dummr,dummv)
  use def_domain
  use def_kintyp
  use mod_postpr
  use def_master
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), intent(in) :: dummi(*)
  real(rp),    intent(in) :: dummr(*)
  real(rp),    intent(in) :: dummv(ndime,*)
  integer(ip)             :: ipoin,ipois,ieles,kfl_markm
  character(5)            :: wopos(2)

  if( IMASTER ) return

  lun_outpu_dom = 90
  lun_postp     = 91
  kfl_markm     =  0
  if( itask < 0 ) kfl_markm = -1
  open( unit=lun_outpu_dom, file='subd'//trim(intost(kfl_paral))//'.post.msh',status='unknown')

  ipois = 0
  ieles = 0
  if(kfl_paral==2) then
     !ipois=9833
     !ieles=19200
  end if
!  call gidele(&
!       iesta_dom,iesto_dom,mnode,npoin,nelem,lun_outpu_dom,lexis,&
!       ltype,lnods,coord,ltype,ieles,ipois,'VOL',0_ip)

  !call geogid()
 
  if( itask > 0 ) then
     open( unit=lun_postp,     file='subd'//trim(intost(kfl_paral))//'.post.res',status='unknown')
     write(lun_postp,'(a)') 'GiD Post Results File 1.0'
     write(lun_postp,'(a)') ' '
     wopos(1) = 'PRUEB'
     if( itask == 3 ) then
        write(lun_postp,2) wopos(1),'PRUEBA',0.0_rp,'Vector'
        write(lun_postp,3) wopos(1)//'_X',wopos(1)//'_Y',wopos(1)//'_Z'
     else
        write(lun_postp,2) wopos(1),'PRUEBA',0.0_rp,'Scalar'
        write(lun_postp,3) wopos(1)
     end if
     write(lun_postp,1) 'Values'
     if( itask == 1 ) then
        do ipoin = 1,npoin
           write(lun_postp,4) ipoin+ipois,dummi(ipoin)
        end do
     else if( itask == 2 ) then
        do ipoin = 1,npoin
           write(lun_postp,5) ipoin+ipois,dummr(ipoin)
        end do
     else if( itask == 3 ) then
        do ipoin = 1,npoin
           write(lun_postp,5) ipoin+ipois,dummv(1:ndime,ipoin)
        end do
     end if
     write(lun_postp,1) 'End values'
     close(91)
  end if

  close(90)
  !
  ! GiD formats
  !
1 format(a)
2 format('Result ',a,' ',a,' ',e14.8,' ',a,' OnNodes')
3 format('ComponentNames ',a)
4 format(i9, 3(1x,i7))
5 format(i9, 3(1x,e16.8E3))

end subroutine possla
