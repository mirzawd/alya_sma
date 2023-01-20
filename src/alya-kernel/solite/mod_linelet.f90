!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Preconditioner
!> @{
!> @name    Linelet preconditioner
!> @file    mod_linelet.f90
!> @author  Guillaume Houzeaux
!> @date    23/06/2017
!> @brief   ToolBox for linelet
!> @details ToolBox for linelet
!>
!------------------------------------------------------------------------

module mod_linelet

  use def_parame
  use def_master
  use def_domain
  use def_solver
  use mod_graphs,         only : graphs_number_to_linked_list
  use mod_elmgeo,         only : element_type
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_communications, only : PAR_SUM
  implicit none

  private

  public :: linelet_initialization
  public :: linelet_factorization
  public :: linelet_solution

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-26
  !> @brief   This routine computes the linelets preconditioner
  !> @details Variables computed
  !>          NLINE .......... Number of linelets
  !>          LLINE(NLINE) ... Starting point of line NLINE
  !>          LMARK(IPOIN) ... Line of IPOIN
  !>          LRENU(IPOIN)
  !>          LRENUP(IPOIN)
  !> 
  !-----------------------------------------------------------------------

  subroutine linelet_initialization()

    integer(ip)          :: ipoin,icont,izdom,jpoin,nline,iline,nlin1,ipoi1
    integer(ip)          :: mpop2,ielem,inode,lsize,jelem,ilisn,iedgg,nlelp
    integer(ip)          :: lnods_loc(mnode),jnode,pelty
    integer(ip)          :: imeth
    real(rp)             :: tole2
    integer(ip), pointer :: lmark(:)
    integer(ip), pointer :: c_new(:)
    integer(ip), pointer :: r_new(:)
    integer(ip), pointer :: lline(:)
    integer(ip), pointer :: ledgl(:)
    real(rp),    pointer :: redge(:)
    integer(ip), pointer :: nepoi_new(:)
    integer(ip), pointer :: pelpo_new(:)
    integer(ip), pointer :: lelpo_new(:)
    logical(lg), pointer :: mask_node(:)

    nline = 0 ! Number of linelets
    nlin1 = 0 ! Number of linelets crossed by boundaries
    ipoi1 = 0 ! Percentage of linelet nodes

    nullify(lmark)
    nullify(c_new)
    nullify(r_new)
    nullify(lline)
    nullify(ledgl)
    nullify(redge)
    nullify(nepoi_new)
    nullify(pelpo_new)
    nullify(lelpo_new)
    nullify(mask_node)

    if( INOTMASTER ) then

       if( solve_sol(1) % kfl_linty < 0 ) then

          !----------------------------------------------------------------
          !
          ! Prescribed linelets through a field
          !
          !----------------------------------------------------------------

          call prelin()

       else

          !----------------------------------------------------------------
          !
          ! Find linelets
          !
          !----------------------------------------------------------------
          !
          ! Tolerance (squared) on the anisotropy
          ! Ratio = edge_max^2/edge^2 > 100^2
          !
          tole2 = solve_sol(1)%toler**2 
          !
          ! Graph method
          ! = 1 ... based on geometrical edges
          ! = 2 ... based on matrix graph
          !
          imeth = 1

          if( imeth == 1 ) then
             !
             ! Create a mask
             !
             if( solve_sol(1) % kfl_linty == 3 ) then
                call memory_alloca(memit,'MASK_NODE','crelin',mask_node,npoin)
                do ielem = 1,nelem
                   pelty = ltype(ielem)
                   if( element_type(pelty) % topology == 0 .or. element_type(pelty) % topology == 2 ) then
                      do inode = 1,element_type(pelty) % number_nodes
                         ipoin = lnods(inode,ielem)
                         mask_node(ipoin) = .true.
                      end do
                   end if
                end do
             end if
             !
             ! Allocate memory for NEPOI_NEW and compute it
             !
             call memory_alloca(memit,'NEPOI_NEW','crelin',nepoi_new,npoin)
             mpop2 = 0
             do ielem = 1,nelem
                mpop2 = mpop2 + lnnod(ielem)*lnnod(ielem)
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   nepoi_new(ipoin) = nepoi_new(ipoin) + 1
                end do
             end do
             !
             ! Allocate memory for PELPO_NEW and compute it
             !
             call memory_alloca(memit,'PELPO_NEW','crelin',pelpo_new,npoin+1)
             pelpo_new(1) = 1
             do ipoin = 1,npoin
                pelpo_new(ipoin+1) = pelpo_new(ipoin) + nepoi_new(ipoin)
             end do
             !
             ! Allocate memory for LELPO_NEW and construct the list
             !
             nlelp = pelpo_new(npoin+1)
             call memory_alloca(memit,'LELPO_NEW','crelin',lelpo_new,nlelp)
             do ielem = 1,nelem
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   lelpo_new(pelpo_new(ipoin)) = ielem
                   pelpo_new(ipoin) = pelpo_new(ipoin) + 1
                end do
             end do
             !
             ! Recompute PELPO_NEW and maximum number of element neighbors MEPOI
             !
             pelpo_new(1) =  1
             do ipoin = 1,npoin
                pelpo_new(ipoin+1) = pelpo_new(ipoin) + nepoi_new(ipoin)
             end do
             !
             ! Allocate memory
             !
             call memory_alloca(memit,'LEDGL','crelin',ledgl,mpop2)
             call memory_alloca(memit,'R_NEW','crelin',r_new,npoin+1)    
             !
             ! Construct the array of indexes
             !
             r_new(1) = 1
             if( solve_sol(1) % kfl_linty == 3 ) then
                do ipoin = 1,npoin                
                   lsize = 0
                   if( mask_node(ipoin) ) then
                      do ielem = pelpo_new(ipoin),pelpo_new(ipoin+1)-1
                         jelem = lelpo_new(ielem)
                         inode = 0
                         do jnode = 1,lnnod(jelem)
                            jpoin = lnods(jnode,jelem)
                            if( mask_node(jpoin) ) then
                               inode = inode + 1
                               lnods_loc(inode) = jpoin
                            end if
                         end do
                         if( inode > 0 ) then
                            call mergl5( ipoin, ledgl(r_new(ipoin)), lsize, lnods_loc,inode )
                         end if
                      end do
                   end if
                   r_new(ipoin+1) = r_new(ipoin) + lsize
                end do
             else
                do ipoin = 1,npoin
                   lsize = 0
                   do ielem = pelpo_new(ipoin),pelpo_new(ipoin+1)-1
                      jelem = lelpo_new(ielem)
                      call mergl5( ipoin, ledgl(r_new(ipoin)), lsize, lnods(:,jelem), lnnod(jelem) )
                   end do
                   r_new(ipoin+1) = r_new(ipoin) + lsize
                end do
             end if

             nedge = r_new(npoin+1) - 1           
             !
             ! Allocate new graph
             !
             call memory_alloca(memit,'C_NEW','crelin',c_new,nedge)
             !
             ! Fill in edge graph
             !    
             iedgg = 0
             do ipoin = 1,npoin
                do ilisn = 1,r_new(ipoin+1)-r_new(ipoin)
                   iedgg          = iedgg + 1
                   jpoin          = ledgl(iedgg)
                   c_new(iedgg)   = jpoin
                end do
             end do
             !
             ! Deallocate memory
             !
             call memory_deallo(memit,'LEDGL    ','crelin',ledgl)    
             call memory_deallo(memit,'LELPO_new','crelin',lelpo_new)    
             call memory_deallo(memit,'PELPO_new','crelin',pelpo_new)    
             call memory_deallo(memit,'NEPOI_new','crelin',nepoi_new)   

          else
             !
             ! Count the number of new internal edges without the diagonal
             !
             icont = 0_ip
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jpoin = c_dom(izdom)
                   if( jpoin /= ipoin ) icont = icont+1
                end do
             end do
             nedge = icont
             !
             ! Allocate new graph
             !
             call memory_alloca(memit,'R_NEW','crelin',r_new,npoin+1)    
             call memory_alloca(memit,'C_NEW','crelin',c_new,nedge)
             !
             ! Copy the local graph and take out the diagonal
             !
             r_new(1) = 1
             icont    = 0
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jpoin = c_dom(izdom)
                   if( jpoin /= ipoin )then
                      icont = icont+1
                      c_new(icont) = jpoin
                   end if
                end do
                r_new(ipoin+1) = icont+1
             end do
          end if
          !
          ! Allocate the array of the length of the new edges
          !
          call memory_alloca(memit,'REDGE','crelin',redge,nedge)
          !
          ! Create source points
          ! 
          call memory_alloca(memit,'LMARK','crelin',lmark,npoin)

          call inilin(nedge,tole2,redge,r_new,c_new,lmark)
          !
          ! Cancel
          !
          if( associated(mask_node) ) then
             do ipoin = 1,npoin
                if( .not. mask_node(ipoin) ) lmark(ipoin) = 0
             end do
             call memory_deallo(memit,'MASK_NODE','crelin',mask_node)
          end if
          !
          ! NLINE: number of linelets
          !
          nline = 0_ip
          do ipoin = 1,npoin
             if( lmark(ipoin) == 1 ) nline = nline + 1
          end do

          solve_sol(1) % nline = nline
          iline = nline
          !
          ! Allocate lline with nline+1 for after
          !
          call memory_alloca(memit,'LLINE','crelin',solve_sol(1)%lline,nline+1_ip)
          lline => solve_sol(1) % lline
          !
          ! Mark lline and lmark
          !
          nline = 0
          do ipoin = 1,npoin
             if( lmark(ipoin) == 1 ) then
                nline        = nline + 1
                lline(nline) = ipoin      ! Starting point IPOIN of line NLINE
                lmark(ipoin) = nline      ! Line NLINE to which IPOIN is associated
             end if
          end do
          !
          ! Grow linelets from source points
          !
          call grolin(nedge,tole2,redge,r_new,c_new,lmark)
          !
          ! Renumber linelets
          !
          call renlin(r_new,c_new,lmark)          
          !
          ! Deallocate memory
          !
          call memory_deallo(memit,'LMARK','crelin',lmark)   
          call memory_deallo(memit,'REDGE','crelin',redge)   
          call memory_deallo(memit,'C_NEW','crelin',c_new)   
          call memory_deallo(memit,'R_NEW','crelin',r_new)   

       end if
       !
       ! Allocate lpntr
       !
       call memory_alloca(memit,'SOLVE % LPNTR','crelin',solve_sol(1)%lpntr,solve_sol(1)%npntr)
       call memory_alloca(memit,'SOLVE % TRIMA','crelin',solve_sol(1)%trima,solve_sol(1)%npntr*solve_sol(1)%ndofn)

    end if
    !
    ! Output information
    !
    if( INOTMASTER ) then
       nline = solve_sol(1) % nline
       nlin1 = solve_sol(1) % nlin1
       ipoi1 = solve_sol(1) % npoin
    end if
    if( IPARALL ) then
       call PAR_SUM(nline)
       call PAR_SUM(nlin1)
       call PAR_SUM(ipoi1)
       if( IMASTER ) then
          solve_sol(1) % nline = nline
          solve_sol(1) % nlin1 = nlin1        
          solve_sol(1) % npoin = ipoi1
       end if
    end if
    if( INOTSLAVE ) then
       solve_sol(1) % npoin = int( real(ipoi1,rp)/real(npoin_total,rp) * 100.0_rp , ip ) 
    end if

  end subroutine linelet_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-09-24
  !> @brief   Initial linelet points
  !> @details This sub creates the initial points for the linelets
  !> 
  !-----------------------------------------------------------------------

  subroutine inilin(nedge,tole2,redge,r_new,c_new,lmark)

    integer(ip),          intent(in)    :: nedge
    real(rp),             intent(in)    :: tole2
    real(rp),    pointer, intent(inout) :: redge(:)
    integer(ip), pointer, intent(in)    :: r_new(:)
    integer(ip), pointer, intent(inout) :: c_new(:)
    integer(ip), pointer, intent(inout) :: lmark(:)
    integer(ip)                         :: ipoin,izdom,iloc,jloc,jpoin,kpoin
    integer(ip)                         :: idime,neigh
    real(rp)                            :: rloc,rdife,rilen,rjlen,ratio
    !
    ! Loop on internal points
    !
    !if( .not. associated(lmark) ) print*,'a=',nedge
    do ipoin = 1,npoin 
       do izdom = r_new(ipoin),r_new(ipoin+1)-1
          jpoin = c_new(izdom)
          !
          ! Compute the edge length 
          !
          rloc = 0.0_rp
          do idime = 1,ndime 
             rdife = coord(idime,ipoin)-coord(idime,jpoin) 
             rloc = rloc+rdife*rdife
          end do
          redge(izdom) = rloc
       end do
       !
       ! Order the edges (brute force)
       !
       do iloc = r_new(ipoin),r_new(ipoin+1)-2
          rilen = redge(iloc) 
          do jloc = iloc+1,r_new(ipoin+1)-1
             rjlen = redge(jloc)
             if( rilen > rjlen ) then
                redge(iloc) = rjlen   
                redge(jloc) = rilen
                rilen       = rjlen   
                neigh       = c_new(iloc)
                c_new(iloc) = c_new(jloc)
                c_new(jloc) = neigh
             end if
          end do
       end do
    end do
    !
    !   Create the origin points 
    !
    if( solve_sol(1) % kfl_linty == -999 ) then

       do ipoin = 1,npoin
          if( solve_sol(1) % limli(ipoin) == 1 ) then
             lmark(ipoin) = 1
          end if
       end do

    else if( solve_sol(1) % kfl_linty == 0 ) then
       !
       ! Mark any interior or boundary node
       !
       do ipoin = 1,npoin
          if( r_new(ipoin+1)-1 >= r_new(ipoin) ) then
             ratio = redge(r_new(ipoin+1)-1)/redge(r_new(ipoin))
             if( ratio > tole2 ) lmark(ipoin) = 1
          end if
       end do

    else if( solve_sol(1) % kfl_linty == 1 .or. solve_sol(1) % kfl_linty == 3 ) then
       !
       ! Mark only from boundary 
       !
       do ipoin = 1,npoin
          if( lpoty(ipoin) > 0 .and. r_new(ipoin+1)-1 >= r_new(ipoin) ) then
             ratio = redge(r_new(ipoin+1)-1)/redge(r_new(ipoin)) 
             if( ratio > tole2 ) lmark(ipoin) = 1                    
          end if
       end do

    else if( solve_sol(1) % kfl_linty == 2 ) then
       !
       ! Where value is prescribed
       !
       do ipoin = 1,npoin
          if( lpoty(ipoin) > 0 .and. solve_sol(1) % limli(ipoin) /= 0 .and. r_new(ipoin+1)-1 >= r_new(ipoin) ) then
             ratio = redge(r_new(ipoin+1)-1)/redge(r_new(ipoin)) 
             if( ratio > tole2 ) lmark(ipoin) = 1                    
          end if
       end do

    end if
    !
    ! Check for repeated sources on the same linelet
    !
    do ipoin = 1,npoin
       !
       ! Has the point been marked for potential source
       ! 
       if( lmark(ipoin) == 1 ) then
          jpoin = c_new(r_new(ipoin))
          kpoin = c_new(r_new(jpoin))
          if( ipoin /= kpoin ) then
             !
             ! kpoin should be the source
             !
             lmark(ipoin) = 0_ip 
          else
             !
             ! Must check one order further
             ! 
             rilen = redge(r_new(ipoin)+1)
             rjlen = redge(r_new(jpoin)+1)
             if( rilen < rjlen ) then
                lmark(ipoin) = 0_ip
             else
                lmark(jpoin) = 0_ip
             end if
          end if
       end if
    end do
    !
    ! Take off interior nodes
    !
    if( solve_sol(1) % kfl_linty /= 0 ) then
       do ipoin = 1,npoin
          if( lpoty(ipoin) == 0 ) lmark(ipoin) = 0_ip
       end do
    end if

  end subroutine inilin

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-09-24
  !> @brief   Grows linelets
  !> @details This sub grows the linelets from the seed points 
  !> 
  !-----------------------------------------------------------------------

  subroutine grolin(nedge,tole2,redge,r_new,c_new,lmark)

    integer(ip),          intent(in)    :: nedge
    real(rp),             intent(in)    :: tole2
    real(rp),    pointer, intent(in)    :: redge(:)
    integer(ip), pointer, intent(in)    :: r_new(:)
    integer(ip), pointer, intent(in)    :: c_new(:)
    integer(ip), pointer, intent(inout) :: lmark(:)
    integer(ip)                         :: iline,izdom,ipoin,jpoin,iprog,nline,kpoin
    real(rp)                            :: ratio,rlmax
    integer(ip), pointer                :: lline(:)
    integer(ip), pointer                :: lline_sav(:)
    integer(ip), pointer                :: nlin1(:)

    nullify(lline)
    nullify(lline_sav)
    nullify(nlin1)
    !
    !   Main loop
    !
    lline => solve_sol(1) % lline
    nline =  solve_sol(1) % nline
    iprog =  1

    if( nline > 0 ) then
       allocate(lline_sav(nline))
       do iline = 1,nline
          lline_sav(iline) = lline(iline)
       end do
    end if
    
    open_loop: do

       if( iprog /= 1 ) exit open_loop
       iprog = 0

       do iline = 1,nline
          kpoin = lline(iline)

          if( kpoin > 0 ) then
             ipoin        =  kpoin
             lline(iline) = -kpoin
             rlmax        =  redge(r_new(ipoin+1)-1)
             !
             !   Loop on neighbors
             !
             neigh:do izdom = r_new(ipoin),r_new(ipoin+1)-1
                jpoin = c_new(izdom) 
                !
                ! Already used?
                !
                if( lmark(jpoin) == 0 ) then
                   !
                   !  Good aspect ratio?
                   !
                   ratio = rlmax/redge(izdom) 
                   if( ratio > tole2 ) then      
                      !
                      ! New point found
                      !
                      lmark(jpoin) = iline
                      lline(iline) = jpoin
                      iprog        = 1
                      exit neigh
                   end if
                end if
             end do neigh
          end if
       end do

    end do open_loop
    !
    ! Start linelets from original point
    !
    do iline = 1,nline
       lline(iline) = lline_sav(iline) ! -lline(iline)
    end do
    if( nline > 0 ) deallocate(lline_sav)
    !
    !   Wake not yet ready .......
    !

    !
    ! Count number of linelets with at least one boundary node
    !
    call parlin(lmark)

  end subroutine grolin

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-09-24
  !> @brief   Linelets with at least one boundary node
  !> @details Count number of linelets with at least one boundary node
  !> 
  !-----------------------------------------------------------------------

  subroutine parlin(lmark)

    integer(ip), pointer, intent(inout) :: lmark(:)
    integer(ip)                         :: iline,nline,ipoin
    integer(ip), pointer                :: nlin1(:)
    
    nline = solve_sol(1) % nline
    solve_sol(1) % nlin1 = 0
    solve_sol(1) % npoin = 0
    if( nline > 0 ) then
       allocate( nlin1(nline) )
       do iline = 1,nline
          nlin1(iline) = 0
       end do
       do ipoin = npoi1+1,npoin
          iline = lmark(ipoin)
          if( iline > 0 ) nlin1(iline) = 1
       end do
       do iline = 1,nline
          solve_sol(1) % nlin1 = solve_sol(1) % nlin1 + nlin1(iline)
       end do

       do ipoin = 1,npoin
          if( lmark(ipoin) > 0 ) solve_sol(1) % npoin = solve_sol(1) % npoin + 1
       end do

       deallocate( nlin1 )

    end if

  end subroutine parlin

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-09-24
  !> @brief   Linelet renumbering
  !> @details This subroutine renumbers the points respect to the linelets 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine renlin(r_new,c_new,lmark)
    
    integer(ip), pointer, intent(in) :: r_new(:)
    integer(ip), pointer, intent(in) :: c_new(:)
    integer(ip), pointer, intent(in) :: lmark(:)
    integer(ip)                      :: nrenu,ipoin,iline
    integer(ip)                      :: onlin,izdom,jpoin,nline  
    integer(ip), pointer             :: lrenu(:)
    integer(ip), pointer             :: lrenup(:)
    integer(ip), pointer             :: lline(:)
    !
    ! Allocate
    !
    call memory_alloca(memit,'SOLVE % LRENU' ,'renlin',solve_sol(1) % lrenu ,npoin)
    call memory_alloca(memit,'SOLVE % LRENUP','renlin',solve_sol(1) % lrenup,npoin)
    ! 
    ! Initialize renumbering counter
    !
    nrenu  =  0_ip
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    lline  => solve_sol(1) % lline
    nline  =  solve_sol(1) % nline 
    solve_sol(1) % npntr = 0
    !
    ! Loop on lines
    !
    ipoin    = lline(1)
    lline(1) = 1

    do iline = 1,nline
       
       nrenu                = nrenu + 1
       lrenu(ipoin)         = nrenu
       lrenup(nrenu)        = ipoin
       onlin                = 1_ip
       solve_sol(1) % npntr = solve_sol(1) % npntr + 1  !count for trima

       open_loop: do
          !
          ! Did we exhaust the linelet
          !          
          if( onlin /= 1 ) exit open_loop
          onlin = 0_ip
          ! 
          ! Loop on neighbors
          !      
          neigh: do izdom = r_new(ipoin),r_new(ipoin+1)-1
             jpoin = c_new(izdom)
             !
             ! Has the point been renumbered already?  
             ! 
             if( lrenu(jpoin) == 0 ) then
                !
                ! Does this point belong to the same linelet
                !                      
                if( lmark(jpoin) == iline ) then
                   !
                   ! Renumber the point
                   !
                   nrenu                = nrenu + 1
                   lrenu(jpoin)         = nrenu
                   lrenup(nrenu)        = jpoin
                   ipoin                = jpoin
                   onlin                = 1_ip 
                   solve_sol(1) % npntr = solve_sol(1) % npntr + 2 !count for trima  
                   exit neigh 
                end if
             end if
          end do neigh
       end do open_loop
       !
       ! Next first point and remember the pointer in LLINE
       !
       ipoin          = lline(iline+1)
       lline(iline+1) = nrenu+1
    end do
    !
    ! Remember the points in linelet
    !
    solve_sol(1) % nlpntr = nrenu
    !
    ! Renumber the points not belonging to any linelet
    !
    do ipoin = 1,npoin
       if( lrenu(ipoin) == 0 ) then
          nrenu                = nrenu+1
          lrenu(ipoin)         = nrenu
          lrenup(nrenu)        = ipoin
          solve_sol(1) % npntr = solve_sol(1) % npntr + 1  !count for trima 
       end if
    end do

  end subroutine renlin

  subroutine matlin()
    !
    !     This subroutine computes the pointer to the graph for fast factorization
    !     if the matrix changes in time
    !
    use def_kintyp
    use def_domain, only    :  r_sym,c_sym,npoin
    use def_solver, only    :  solve_sol
    implicit none
    integer(ip)             :: icont,iline,isto0,isto1,ipoin
    integer(ip)             :: ia,ib,jp,izdom,jpoin,kpoin,nline
    integer(ip), pointer    :: lrenu(:),lrenup(:),lpntr(:),lline(:)

    icont  =  0_ip
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    lpntr  => solve_sol(1) % lpntr
    lline  => solve_sol(1) % lline
    nline  =  solve_sol(1) % nline 
    !
    ! Loop on linelets
    !
    do iline=1,nline
       !
       !     Loop on points of iline
       ! 
       isto0 = lline(iline)
       isto1 = lline(iline+1)-1
       !
       !     First point of the linelet
       !
       ipoin = lrenup(isto0)
       isto0 = isto0+1
       !
       !     Pointer to the diagonal
       !
       icont        = icont+1
       lpntr(icont) = r_sym(ipoin+1)-1 
       !      
       do jp = isto0,isto1
          jpoin = lrenup(jp)
          !
          !     Find the place in the matrix for edge IPOIN-JPOIN
          !
          !
          !     Pointer to the previous in the linelet
          !
          if( ipoin > jpoin )then
             ia = ipoin
             ib = jpoin
          else
             ia = jpoin
             ib = ipoin   
          end if

          neigh: do izdom = r_sym(ia),r_sym(ia+1)-1
             kpoin = c_sym(izdom)
             if( kpoin == ib ) exit neigh
          end do neigh
          icont        = icont + 1
          lpntr(icont) = izdom
          !
          !     Pointer to the diagonal
          !
          icont        = icont + 1
          lpntr(icont) = r_sym(jpoin+1)-1 
          !
          !     Goto next edge
          !
          ipoin        = jpoin
       end do
    end do
    !
    !     Points not in linelets
    !
    do ipoin = solve_sol(1) % nlpntr+1,npoin
       icont = icont + 1
       lpntr(icont) = r_sym(ipoin+1)-1
    end do

  end subroutine matlin

  subroutine linelet_factorization(ndofn,an,invdiag)
    !
    ! This subroutine performs a cholesky factorization of the tridiagonal linelet matrix  
    !
    implicit none 
    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(in)    :: an(*),invdiag(*)

    if( IMASTER ) return

    if( ndofn == 1 ) then
       call facli1(an,invdiag)
    else
       call facli2(ndofn,an,invdiag,solve_sol(1)%trima)
    end if

  end subroutine linelet_factorization

  subroutine facli1(an,invdiag)
    !
    ! This subroutine performs a cholesky factorization of the tridiagonal linelet matrix  
    !

    real(rp),    intent(in)    :: an(*),invdiag(*)
    integer(ip)                :: ipntr,icont,ipoin,ibopo,ipnt,nline
    integer(ip)                :: iline,isto0,isto,isto1,iplace,npntr
    integer(ip)                :: iperi,jpoin
    real(rp)                   :: pivot,extra
    integer(ip), pointer       :: limli(:),lpntr(:),lrenu(:),lrenup(:)
    integer(ip), pointer       :: lmark(:),lline(:)
    real(rp),    pointer       :: trima(:)
    logical(lg), pointer       :: gmark(:)

    nullify(lmark)
    nullify(gmark)
    !
    ! Only the slaves
    !
    if( IMASTER ) return
    !
    ! Factorize matrix
    !
    if( solve_sol(1)%kfl_factl == 0 ) then
       solve_sol(1)%kfl_factl = 1
       if( solve_sol(1)%kfl_symme == 1 ) then
          call matlin()
       else
          call matliu()
       end if
    end if

    nline  =  solve_sol(1) % nline
    npntr  =  solve_sol(1) % npntr
    trima  => solve_sol(1) % trima
    lpntr  => solve_sol(1) % lpntr
    limli  => solve_sol(1) % limli
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    lline  => solve_sol(1) % lline
    !
    ! Fill in matrix
    !
    do ipntr = 1,npntr
       trima(ipntr) = an(lpntr(ipntr))
    end do
    !
    ! Correspondance in trimat, mark < 0 the first point of each linelet  
    !
    call memory_alloca(memit,'LMARK','factlin',lmark,npoin)
    call memory_alloca(memit,'GMARK','factlin',gmark,npoin)
    do ipoin=1,npoin
       gmark(ipoin) = .false.
    end do
    icont = 0
    do iline = 1,nline
       isto0 = lline(iline)
       isto1 = lline(iline+1)-1
       !
       ! First point  
       ! 
       icont        =  icont + 1
       ipoin        =  lrenup(isto0)
       lmark(ipoin) = -icont
       do isto = isto0+1,isto1
          icont        = icont+2
          ipoin        = lrenup(isto)
          lmark(ipoin) = icont
       end do
       gmark(ipoin) = .true.
    end do
    !
    ! Imposed points: IPLACE < 0 is the first node of the linelet
    !  
    if( memory_size(limli) == nbopo .and. nbopo /= 0 ) then
       call runend('LINELET MAY NOT WORK')
       do ipoin = 1,npoin
          ibopo = lpoty(ipoin)
          if( ibopo /= 0 ) then
             if( limli(ibopo) > 0 )then
                iplace = lmark(ipoin)
                if( iplace < 0 ) then
                   trima(-iplace)   = 1.0_rp/invdiag(ipoin)          ! diagonal
                   trima(-iplace+1) = 0.0_rp                         ! extra-diagonal
                else if( iplace > 0 ) then
                   trima(iplace)   = 1.0_rp/invdiag(ipoin)           ! diagonal
                   trima(iplace-1) = 0.0_rp                          ! extra-diagonal
                   if( .not. gmark(ipoin) ) trima(iplace+1) = 0.0_rp ! extra-diagonal
                end if
             end if
          end if
       end do
    else
       do ipoin = 1,npoin
          if( limli(ipoin) > 0 ) then
             iplace = lmark(ipoin)
             if( iplace < 0 ) then
                trima(-iplace)   = 1.0_rp/invdiag(ipoin)             ! diagonal
                trima(-iplace+1) = 0.0_rp                            ! extra-diagonal
             else if( iplace > 0 ) then
                trima(iplace)   = 1.0_rp/invdiag(ipoin)              ! diagonal
                trima(iplace-1) = 0.0_rp                             ! extra-diagonal
                if( .not. gmark(ipoin) ) trima(iplace+1) = 0.0_rp    ! extra-diagonal
             end if
          end if
       end do
    end if
    !
    ! Slave neighbor nodes: Impose diagonal
    ! 
    do ipoin = npoi1+1,npoin
       iplace = lmark(ipoin)
       if( iplace < 0 ) then
          trima(-iplace)   = 1.0_rp/invdiag(ipoin)                   ! diagonal
          trima(-iplace+1) = 0.0_rp                                  ! extra-diagonal
       else if( iplace > 0 ) then
          trima(iplace)   = 1.0_rp/invdiag(ipoin)                    ! diagonal
          trima(iplace-1) = 0.0_rp                                   ! extra-diagonal
          if( .not. gmark(ipoin) ) trima(iplace+1) = 0.0_rp          ! extra-diagonal
       end if
    end do
    !
    ! Periodicity
    !
    do iperi = 1,nperi
       ipoin = lperi(1,iperi)
       jpoin = lperi(2,iperi)
       if( ipoin > 0 ) then
          iplace = lmark(ipoin)
          if( iplace < 0 ) then
             trima(-iplace)   = 1.0_rp/invdiag(ipoin)                   ! diagonal
             trima(-iplace+1) = 0.0_rp                                  ! extra-diagonal
          else if( iplace > 0 ) then
             trima(iplace)   = 1.0_rp/invdiag(ipoin)                    ! diagonal
             trima(iplace-1) = 0.0_rp                                   ! extra-diagonal
             if( .not. gmark(ipoin) ) trima(iplace+1) = 0.0_rp          ! extra-diagonal
          end if
       end if
       if( jpoin > 0 ) then
          ipoin = jpoin
          iplace = lmark(ipoin)
          if( iplace < 0 ) then
             trima(-iplace)   = 1.0_rp                                  ! diagonal
             trima(-iplace+1) = 0.0_rp                                  ! extra-diagonal
          else if( iplace > 0 ) then
             trima(iplace)   = 1.0_rp                                   ! diagonal
             trima(iplace-1) = 0.0_rp                                   ! extra-diagonal
             if( .not. gmark(ipoin) ) trima(iplace+1) = 0.0_rp          ! extra-diagonal
          end if
       end if
    end do
    !
    ! Factorization of the matrix
    ! 
    icont = 0_ip

    do iline = 1,nline

       icont        = icont+1 
       trima(icont) = sqrt(abs(trima(icont)))
       pivot        = trima(icont)

       do ipnt = lline(iline)+1,lline(iline+1)-1
          !
          ! Extra diagonal
          !
          icont        = icont+1
          trima(icont) = trima(icont)/pivot
          extra        = trima(icont)
          !
          ! Diagonal
          !         
          icont        = icont+1
          trima(icont) = sqrt(abs(trima(icont)-extra*extra))
          pivot        = trima(icont)

       end do
    end do
    !
    ! Invert the diagonal
    !
    do ipoin = solve_sol(1)%nlpntr+1,npoin
       icont        = icont + 1
       trima(icont) = invdiag(ipoin)
    end do
    !
    ! Deallocate LMARK and GMARK
    !
    call memory_deallo(memit,'LMARK','factlin',lmark)
    call memory_deallo(memit,'GMARK','factlin',gmark)

  end subroutine facli1

  subroutine linelet_solution(ndofn,rhsid,wa2,invdiag)
    !-----------------------------------------------------------------------
    !****f* domain/sollin
    ! NAME
    !    sollin
    ! DESCRIPTION
    !    This sub solves the forward/backward substitution of the 
    !    tridiagonal matrix
    ! INPUT
    !    RHSID: Right-hand side to be modified
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(in)    :: invdiag(*)
    real(rp),    intent(inout) :: rhsid(*)
    real(rp),    intent(out)   :: wa2(*)

    if( IMASTER ) return

    if( ndofn == 1 ) then
       call solli1(rhsid,wa2,invdiag)
    else
       call solli2(ndofn,rhsid,wa2,invdiag,solve_sol(1)%trima)
    end if

  end subroutine linelet_solution

  subroutine solli1(rhsid,wa2,invdiag)
    !-----------------------------------------------------------------------
    !****f* domain/sollin
    ! NAME
    !    sollin
    ! DESCRIPTION
    !    This sub solves the forward/backward substitution of the 
    !    tridiagonal matrix
    ! INPUT
    !    RHSID: Right-hand side to be modified
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------

    real(rp),    intent(in)    :: invdiag(*)
    real(rp),    intent(inout) :: rhsid(*)
    real(rp),    intent(out)   :: wa2(*)
    integer(ip)                :: ipoin,ncont,nlpntr,nline,npntr,iline,iptdi
    real(rp)                   :: rval,extra
    integer(ip), pointer       :: lrenu(:),lrenup(:),lline(:)
    real(rp),    pointer       :: trima(:)
    !
    ! Only the slaves
    !
    if( IMASTER ) return
    nline  =  solve_sol(1) % nline
    npntr  =  solve_sol(1) % npntr
    trima  => solve_sol(1) % trima
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    nlpntr =  solve_sol(1) % nlpntr
    lline  => solve_sol(1) % lline
    !
    ! Start solving
    !
    do ipoin = 1,npoin
       wa2(ipoin) = rhsid(lrenup(ipoin))
    end do
    !
    ! Beginning of the diagonal part in trima
    !
    iptdi = npntr-(npoin-nlpntr)+1
    ! 
    ! Forward substitution
    !       
    ncont = 0

    do iline = 1,nline
       ipoin      = lline(iline)
       ncont      = ncont+1
       wa2(ipoin) = wa2(ipoin)/trima(ncont)
       rval       = wa2(ipoin)
       do ipoin = lline(iline)+1,lline(iline+1)-1
          ncont      = ncont+1
          extra      = trima(ncont)
          ncont      = ncont+1
          wa2(ipoin) = (wa2(ipoin)-extra*rval)/trima(ncont)
          rval       = wa2(ipoin)
       end do
    end do
    ! 
    ! Backward substitution
    !
    ncont = iptdi

    do iline = nline,1,-1

       ipoin      = lline(iline+1) - 1 
       ncont      = ncont - 1 
       wa2(ipoin) = wa2(ipoin)/trima(ncont)
       rval       = wa2(ipoin)

       do ipoin = lline(iline+1)-2,lline(iline),-1
          ncont      = ncont - 1
          extra      = trima(ncont)
          ncont      = ncont - 1
          wa2(ipoin) = (wa2(ipoin)-extra*rval)/trima(ncont)
          rval       = wa2(ipoin)
       end do
    end do
    !
    ! Points not in linelets
    !
    ncont = iptdi-1

    do ipoin = nlpntr+1,npoin
       wa2(ipoin) = wa2(ipoin) * invdiag(lrenup(ipoin))
    end do
    !
    ! Copy in input 
    !
    do ipoin = 1,npoin
       rhsid(ipoin) = wa2(lrenu(ipoin))
    end do

  end subroutine solli1

  subroutine matliu()
    !
    !     This subroutine computes the pointer to the graph for fast factorization
    !     if the matrix changes in time
    !
    use def_kintyp
    use def_domain, only    :  r_dom,c_dom,npoin
    use def_solver, only    :  solve_sol
    implicit none
    integer(ip)             :: icont,iline,isto0,isto1,ipoin,jzdom
    integer(ip)             :: ia,ib,jp,izdom,jpoin,kpoin,nline
    integer(ip), pointer    :: lrenu(:)
    integer(ip), pointer    :: lrenup(:)
    integer(ip), pointer    :: lpntr(:)
    integer(ip), pointer    :: lline(:)

    icont  =  0_ip
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    lpntr  => solve_sol(1) % lpntr
    lline  => solve_sol(1) % lline
    nline  =  solve_sol(1) % nline 
    !
    !     Loop on linelets
    !
    do iline = 1,nline
       !
       !     Loop on points of iline
       ! 
       isto0 = lline(iline)
       isto1 = lline(iline+1)-1
       !
       !     First point of the linelet
       !
       ipoin = lrenup(isto0)
       isto0 = isto0+1
       !
       !     Pointer to the diagonal
       !
       icont        = icont+1
       do_izdom: do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
          if( c_dom(izdom) == ipoin ) then
             lpntr(icont) = izdom
             exit do_izdom
          end if
       end do do_izdom
       !      
       do jp=isto0,isto1
          jpoin=lrenup(jp)
          !
          !     Find the place in the matrix 
          !
          !
          !     Pointer to the previous in the linelet
          !
          if( ipoin > jpoin )then
             ia = ipoin
             ib = jpoin
          else
             ia = jpoin
             ib = ipoin   
          end if

          neigh: do izdom = r_dom(ia),r_dom(ia+1)-1
             kpoin = c_dom(izdom)
             if( kpoin == ib ) exit neigh
          end do neigh
          icont        = icont+1
          lpntr(icont) = izdom
          !
          !     Pointer to the diagonal
          !
          icont        = icont+1

          do_jzdom: do jzdom = r_dom(jpoin),r_dom(jpoin+1)-1
             if( c_dom(jzdom) == jpoin ) then
                lpntr(icont) = jzdom
                exit do_jzdom
             end if
          end do do_jzdom
          ipoin = jpoin
       end do
    end do
    !
    !     Points not in linelets
    !
    do ipoin = solve_sol(1)%nlpntr+1,npoin
       icont = icont+1
       do_kzdom: do jzdom = r_dom(ipoin),r_dom(ipoin+1)-1
          if( c_dom(jzdom) == ipoin ) then
             lpntr(icont) = jzdom
             exit do_kzdom
          end if
       end do do_kzdom
    end do

  end subroutine matliu


  subroutine facli2(ndofn,an,invdiag,trima)
    !
    ! This subroutine performs a cholesky factorization of the tridiagonal linelet matrix  
    !
    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(in)    :: an(ndofn,ndofn,*),invdiag(ndofn,*)
    real(rp),    intent(out)   :: trima(ndofn,*)
    integer(ip)                :: ipntr,icont,ipoin,ibopo,ipnt,nline,idofn
    integer(ip)                :: iline,isto0,isto,isto1,iplace,npntr
    real(rp)                   :: pivot(ndofn),extra(ndofn)
    integer(ip), pointer       :: limli(:),lpntr(:),lrenu(:),lrenup(:)
    integer(ip), pointer       :: lmark(:),lline(:)
    logical(lg), pointer       :: gmark(:)

    nullify(limli)
    nullify(lpntr)
    nullify(lrenu)
    nullify(lrenup)
    nullify(lmark)
    nullify(lline)
    nullify(gmark)
    !
    ! Only the slaves
    !
    if( IMASTER ) return
    !
    ! Factorize matrix
    !
    if( solve_sol(1)%kfl_factl == 0 ) then
       solve_sol(1)%kfl_factl = 1
       if( solve_sol(1)%kfl_symme == 1 ) then
          call matlin()
       else
          call matliu()
       end if
    end if

    nline  =  solve_sol(1) % nline
    npntr  =  solve_sol(1) % npntr
    lpntr  => solve_sol(1) % lpntr
    limli  => solve_sol(1) % limli
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    lline  => solve_sol(1) % lline
    !
    ! Fill in matrix
    !
    do ipntr = 1,npntr
       do idofn = 1,ndofn
          trima(idofn,ipntr) = an(idofn,idofn,lpntr(ipntr))
       end do
    end do
    !
    ! Correspondance in trimat, mark < 0 the first point of each linelet  
    !
    call memory_alloca(memit,'LMARK','factlin',lmark,npoin)
    call memory_alloca(memit,'GMARK','factlin',gmark,npoin)

    do ipoin=1,npoin
       gmark(ipoin) = .false.
    end do

    icont = 0
    do iline = 1,nline
       isto0 = lline(iline)
       isto1 = lline(iline+1)-1
       !
       ! First point  
       ! 
       icont        = icont+1
       ipoin        = lrenup(isto0)
       lmark(ipoin) = -icont
       do isto = isto0+1,isto1
          icont        = icont+2
          ipoin        = lrenup(isto)
          lmark(ipoin) = icont
       end do
       gmark(ipoin) = .true.
    end do
    !
    ! Imposed points
    !  
    if( memory_size(limli) == nbopo ) then
       call runend('LINELET MAY NOT WORK')
       do ipoin = 1,npoin
          ibopo = lpoty(ipoin)
          if( ibopo /= 0 ) then
             if( limli(ibopo) > 0 )then
                iplace = lmark(ipoin)
                if( iplace < 0 ) then
                   do idofn = 1,ndofn
                      trima(idofn,-iplace)   = 1.0_rp/invdiag(idofn,ipoin)       ! diagonal
                      trima(idofn,-iplace+1) = 0.0_rp   
                   end do
                else if( iplace > 0) then
                   do idofn = 1,ndofn
                      trima(idofn,iplace)   = 1.0_rp/invdiag(idofn,ipoin)        ! diagonal
                      trima(idofn,iplace-1) = 0.0_rp                             ! extra-diagonal
                      if( .not. gmark(ipoin) )  trima(idofn,iplace+1) = 0.0_rp   ! extra-diagonal
                   end do
                end if
             end if
          end if
       end do
    else
       do ipoin = 1,npoin
          if( limli(ipoin) > 0 ) then
             iplace = lmark(ipoin)
             if( iplace < 0 ) then
                do idofn = 1,ndofn
                   trima(idofn,-iplace)   = 1.0_rp / invdiag(idofn,ipoin)    ! diagonal
                   trima(idofn,-iplace+1) = 0.0_rp   
                end do
             else if( iplace > 0 ) then
                do idofn = 1,ndofn
                   trima(idofn, iplace)   = 1.0_rp / invdiag(idofn,ipoin)    ! diagonal
                   trima(idofn, iplace-1) = 0.0_rp                           ! extra-diagonal
                   if( .not. gmark(ipoin) ) trima(idofn, iplace+1) = 0.0_rp  ! extra-diagonal
                end do
             end if
          end if
       end do
    end if
    !
    ! Slave neighbor nodes: Impose diagonal
    ! 
    do ipoin = npoi1+1,npoin
       iplace = lmark(ipoin)
       if( iplace < 0 ) then
          do idofn = 1,ndofn
             trima(idofn,-iplace)   = 1.0_rp/invdiag(idofn,ipoin)          ! diagonal
             trima(idofn,-iplace+1) = 0.0_rp   
          end do
       else if( iplace > 0 ) then
          do idofn = 1,ndofn
             trima(idofn,iplace)   = 1.0_rp/invdiag(idofn,ipoin)           ! diagonal
             trima(idofn,iplace-1) = 0.0_rp                                ! extra-diagonal
             if( .not. gmark(ipoin) ) trima(idofn,iplace+1) = 0.0_rp       ! extra-diagonal
          end do
       end if
    end do
    !
    ! Factorization of the matrix
    ! 
    icont = 0_ip

    do iline = 1,nline

       icont        = icont+1 

       do idofn = 1,ndofn
          trima(idofn,icont) = sqrt(abs(trima(idofn,icont)))
          pivot(idofn)       = trima(idofn,icont)
       end do

       do ipnt = lline(iline)+1,lline(iline+1)-1
          !
          ! Extra diagonal
          !
          icont = icont+1         
          do idofn = 1,ndofn
             trima(idofn,icont) = trima(idofn,icont)/pivot(idofn)
             extra(idofn)       = trima(idofn,icont)
          end do
          !
          ! Diagonal
          !         
          icont = icont+1
          do idofn = 1,ndofn
             trima(idofn,icont) = sqrt(abs(trima(idofn,icont)-extra(idofn)*extra(idofn)))
             pivot(idofn)       = trima(idofn,icont)
          end do

       end do
    end do
    !
    ! Invert the diagonal
    !
    do ipoin = solve_sol(1)%nlpntr+1,npoin
       icont = icont + 1
       do idofn = 1,ndofn
          trima(idofn,icont) = invdiag(idofn,ipoin)
       end do
    end do
    !
    ! Deallocate LMARK and GMARK
    !
    call memory_deallo(memit,'LMARK','factlin',lmark)
    call memory_deallo(memit,'GMARK','factlin',gmark)

  end subroutine facli2

  subroutine solli2(ndofn,rhsid,wa2,invdiag,trima)
    !-----------------------------------------------------------------------
    !****f* domain/sollin
    ! NAME
    !    sollin
    ! DESCRIPTION
    !    This sub solves the forward/backward substitution of the 
    !    tridiagonal matrix
    ! INPUT
    !    RHSID: Right-hand side to be modified
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(in)    :: invdiag(ndofn,*)
    real(rp),    intent(inout) :: rhsid(ndofn,*)
    real(rp),    intent(out)   :: wa2(ndofn,*)
    real(rp),    intent(in)    :: trima(ndofn,*)
    integer(ip)                :: ipoin,ncont,nlpntr,nline,npntr
    integer(ip)                :: iline,iptdi,jpoin,idofn
    real(rp)                   :: rval(ndofn),extra(ndofn)
    integer(ip), pointer       :: lrenu(:),lrenup(:),lline(:)
    !
    ! Only the slaves
    !
    if( IMASTER ) return
    nline  =  solve_sol(1) % nline
    npntr  =  solve_sol(1) % npntr
    lrenu  => solve_sol(1) % lrenu
    lrenup => solve_sol(1) % lrenup
    nlpntr =  solve_sol(1) % nlpntr
    lline  => solve_sol(1) % lline
    !
    ! Start solving
    !
    do ipoin = 1,npoin
       jpoin = lrenup(ipoin)
       do idofn = 1,ndofn
          wa2(idofn,ipoin) = rhsid(idofn,jpoin)
       end do
    end do
    !
    ! Beginning of the diagonal part in trima
    !
    iptdi = npntr-(npoin-nlpntr)+1
    ! 
    ! Forward substitution
    !       
    ncont = 0

    do iline = 1,nline
       ipoin = lline(iline)
       ncont = ncont + 1
       do idofn = 1,ndofn
          wa2(idofn,ipoin) = wa2(idofn,ipoin)/trima(idofn,ncont)
          rval(idofn)      = wa2(idofn,ipoin)
       end do
       do ipoin = lline(iline)+1,lline(iline+1)-1
          ncont = ncont + 1
          do idofn = 1,ndofn
             extra(idofn) = trima(idofn,ncont)
          end do
          ncont = ncont + 1
          do idofn = 1,ndofn
             wa2(idofn,ipoin) = (wa2(idofn,ipoin)-extra(idofn)*rval(idofn))/trima(idofn,ncont)
             rval(idofn)      = wa2(idofn,ipoin)
          end do
       end do
    end do
    ! 
    ! Backward substitution
    !
    ncont = iptdi

    do iline = nline,1,-1

       ipoin      = lline(iline+1) - 1 
       ncont      = ncont - 1 
       do idofn = 1,ndofn
          wa2(idofn,ipoin) = wa2(idofn,ipoin)/trima(idofn,ncont)
          rval(idofn)      = wa2(idofn,ipoin)
       end do

       do ipoin = lline(iline+1)-2,lline(iline),-1
          ncont = ncont - 1
          do idofn = 1,ndofn
             extra(idofn) = trima(idofn,ncont)
          end do
          ncont = ncont - 1
          do idofn = 1,ndofn
             wa2(idofn,ipoin) = (wa2(idofn,ipoin)-extra(idofn)*rval(idofn))/trima(idofn,ncont)
             rval(idofn)      = wa2(idofn,ipoin)
          end do
       end do
    end do
    !
    ! Points not in linelets
    !
    ncont = iptdi-1

    do ipoin = nlpntr+1,npoin
       jpoin = lrenup(ipoin)
       do idofn = 1,ndofn
          wa2(idofn,ipoin) = wa2(idofn,ipoin) * invdiag(idofn,jpoin)
       end do
    end do
    !
    ! Copy in input 
    !
    do ipoin = 1,npoin
       jpoin = lrenu(ipoin)
       do idofn = 1,ndofn
          rhsid(idofn,ipoin) = wa2(idofn,jpoin)
       end do
    end do

  end subroutine solli2

  subroutine linsgs(istar,nbnodes,nbvar,lrenup,an,dd,xx,bb) 
    !----------------------------------------------------------------------
    !****f* mathru/bcsrax
    ! NAME 
    !     bcsrax
    ! DESCRIPTION
    !     Multiply a non symmetric matrix stored in BCSR by a vector
    !     XX = A BB 
    ! INPUT
    !    NBNODES .... Number of equations
    !    NBVAR ...... Number of variables
    !    AN ......... Matrix
    !    JA ......... List of elements
    !    IA ......... Pointer to list of elements
    !    BB ......... Vector
    ! OUTPUT
    !    XX ......... result vector
    ! USES
    ! USED BY
    !***
    !----------------------------------------------------------------------

    integer(ip), intent(in)          :: istar,nbnodes,nbvar
    integer(ip), intent(in)          :: lrenup(*)
    real(rp),    intent(in)          :: an(nbvar,nbvar,*)
    real(rp),    intent(in)          :: dd(nbvar,*)
    real(rp),    intent(in)          :: bb(nbvar,*)
    real(rp),    intent(inout)       :: xx(nbvar,*)
    integer(ip)                      :: ii,jj,kk,ll,col,mm
    real(rp),    pointer             :: yy(:,:)

    if( INOTMASTER ) then

       allocate(yy(nbvar,nbnodes))
       !
       ! NBVAR = whatever
       !
       !*OMP   PARALLEL DO SCHEDULE (GUIDED)               & 
       !*OMP   DEFAULT (NONE)                              &
       !*OMP   PRIVATE ( ii, jj, kk, ll, col, raux)        &
       !*OMP   SHARED ( nbnodes, nbvar, bb, xx, r_dom, c_dom, an)

       !
       ! (L+D) D^-1 (U+D) x = b
       ! 1. (L+D) y = b
       ! 2. z = D y
       ! 3. (U+D) x = z
       !
       do ii = npoi1+1,nbnodes
          do kk = 1,nbvar
             xx(kk,ii) = bb(kk,ii) * dd(kk,ii)
             yy(kk,ii) = xx(kk,ii)
          end do
       end do

       if( nbvar == 1 ) then
          !
          ! 1. (L+D) y = b
          !           
          do kk = istar,npoi1
             !
             ! x_i^{k+1} = 1/a_ii * b_i
             !
             ii = lrenup(kk)

             yy(1,ii) = bb(1,ii) * dd(1,ii)

             do jj  = r_dom(ii),r_dom(ii+1)-1
                col = c_dom(jj) 
                if( col < ii ) then
                   !
                   ! j<i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^{k+1}
                   !
                   yy(1,ii) = yy(1,ii) - dd(1,ii) * an(1,1,jj) * yy(1,col) 
                end if
             end do
          end do
          !
          ! 2. z = D y
          !              
          do kk = 1,npoi1
             ii = lrenup(kk)
             yy(1,ii) = yy(1,ii) / dd(1,ii)
          end do
          !
          ! 3. (U+D) x = z
          !
          do kk = npoi1,istar,-1
             !
             ! x_i^{k+1} = 1/a_ii * b_i
             !
             ii = lrenup(kk)

             xx(1,ii) = yy(1,ii) * dd(1,ii)

             do jj  = r_dom(ii),r_dom(ii+1)-1
                col = c_dom(jj) 
                if( col > ii ) then
                   !
                   ! j>i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^k
                   !
                   xx(1,ii) = xx(1,ii) - dd(1,ii) * an(1,1,jj) * xx(1,col) 
                end if
             end do
          end do

       else
          !
          ! 1. (L+D) y = b
          !           
          do mm = istar,npoi1
             !
             ! x_i^{k+1} = 1/a_ii * b_i
             !
             ii = lrenup(mm)
             do kk = 1,nbvar
                yy(kk,ii) = bb(kk,ii) * dd(kk,ii)
             end do

             do jj  = r_dom(ii),r_dom(ii+1)-1
                col = c_dom(jj) 
                if( col < ii ) then
                   !
                   ! j<i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^{k+1}
                   !
                   do kk = 1,nbvar
                      do ll = 1,nbvar
                         yy(kk,ii) = yy(kk,ii) - dd(kk,ii) * an(ll,kk,jj) * yy(ll,col) 
                      end do
                   end do
                end if
             end do
          end do
          !
          ! 2. z = D y
          !              
          do mm = istar,npoi1
             ii = lrenup(mm)
             do kk = 1,nbvar
                yy(kk,ii) = yy(kk,ii) / dd(kk,ii)
             end do
          end do
          !
          ! 3. (U+D) x = z
          !
          do mm = npoi1,istar,-1
             !
             ! x_i^{k+1} = 1/a_ii * b_i
             !           
             ii = lrenup(mm)
             do kk = 1,nbvar
                xx(kk,ii) = yy(kk,ii) * dd(kk,ii)
             end do

             do jj  = r_dom(ii),r_dom(ii+1)-1
                col = c_dom(jj) 
                if( col > ii ) then
                   !
                   ! j>i: x_i^{k+1} = x_i^{k+1} - 1/a_ii * a_ij * x_j^k
                   !
                   do kk = 1,nbvar
                      do ll = 1,nbvar
                         xx(kk,ii) = xx(kk,ii) - dd(kk,ii) * an(ll,kk,jj) * xx(ll,col) 
                      end do
                   end do
                end if
             end do
          end do

       end if

       deallocate(yy)

    end if

  end subroutine linsgs

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-26
  !> @brief   Linelet
  !> @details Linelet is given by a field
  !> 
  !-----------------------------------------------------------------------

  subroutine prelin()
    integer(ip)          :: iline,ifiel,izdom,nz
    integer(ip)          :: ipoin,jpoin,nline
    integer(ip), pointer :: lmark(:)
    integer(ip), pointer :: lposi(:)
    logical(lg)          :: negative = .false.
    ! 
    ! Initialize renumbering counter
    !
    nullify(lmark,lposi)

    if( INOTMASTER ) then

       !-------------------------------------------------------------------
       !
       ! LMARK(IPOIN) = Linelet node
       !
       !-------------------------------------------------------------------

       ifiel = -solve_sol(1) % kfl_linty
       call memory_alloca(memit,'LPOSI','prelin',lposi,npoin)
       call memory_alloca(memit,'LMARK','prelin',lmark,npoin)
       do ipoin = 1,npoin
          iline = int(xfiel(ifiel) % a(1,ipoin,1),ip)
          if( iline > 0 ) then
             lmark(ipoin) = abs(iline)
             lposi(iline) = 1
          end if
       end do

       !-------------------------------------------------------------------
       !
       ! Count and renumber linelets (global to local)
       !
       !-------------------------------------------------------------------

       nline = 0
       do iline = 1,npoin
          if( lposi(iline) > 0 ) then
             nline = nline + 1
             lposi(iline) = nline
          end if
       end do
       solve_sol(1) % nline = nline
       do ipoin = 1,npoin
          iline = lmark(ipoin)
          if( iline > 0 ) lmark(ipoin) = lposi(iline)
       end do

       solve_sol(1) % nline = nline
       call memory_alloca(memit,'SOLVE % LLINE','prelin',solve_sol(1) % lline,nline+1_ip)

       !-------------------------------------------------------------------
       !
       ! Loop for the seed
       !
       !-------------------------------------------------------------------

       negative = .false.
       do ipoin = 1,npoin
          iline = lmark(ipoin)
          if( iline < 0 ) then
             solve_sol(1) % lline(-iline) = ipoin
             negative = .true.
          end if
       end do

       if( .not. negative ) then  
          do ipoin = 1,npoin
             iline = lmark(ipoin)
             nz = 0
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                if( lmark(jpoin) == iline ) nz = nz + 1
             end do
             if( nz == 2 ) solve_sol(1) % lline(iline) = ipoin
          end do
       end if

       call parlin(lmark)
       call renlin(r_dom,c_dom,lmark) 
       call memory_deallo(memit,'LPOSI','prelin',lposi)   
       call memory_deallo(memit,'LMARK','prelin',lmark)   

    end if

  end subroutine prelin

end module mod_linelet
!> @}
