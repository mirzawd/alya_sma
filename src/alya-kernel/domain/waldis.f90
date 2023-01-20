!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine waldis(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/waldis
  ! NAME
  !    waldis
  ! DESCRIPTION
  !    This routine calculates the wall distance when it is variable depending on the element size (kfl_delta==1).
  !    Beware: for the moment it calculates the wall distance on all boundaries not just boundaries with wall_law.
  ! OUTPUT
  !    ywalb ... Boundary wall distance
  !    ywalp ... Nodal wall distance, projection of ywalb. It is used mainly for postprocesing but also in tur_walgen
  !    to extend the wall distance to the whole domain this extension of the wall distance is added to the distance from the boundary to each 
  !    node to obtain an approximation to the distance from the actual object boundary to the point. 
  !    It would be better to project using only wall_law boundaries not all boundaries but as we don not have the classification of the boundaries here 
  !    (kerner/domain) we have decided to leave it like this for the moment. The aproximation would be poor where a wall law boundary meets an oulet for example.
  !    Think also about wall boundaries.
  !    In the parallel case there might be subdomains without boundaries. In that case one option might be to return. But despite there are no boundaries there might be a
  !    bopo therefore I leave it anyway so that the bopo receives the correct value thanks to rhsmod.
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_memchk
  use def_domain
  use mod_messages,       only : livinf
  use mod_domain,         only : domain_memory_allocate
  use mod_bouder
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iboun,iblty,pnodb,ielem,pnode,inode,inodb,ipoin,idime,pblty,ibopo,ichkm
  real(rp)                :: bonor(ndime),elcog(ndime),bocog(ndime)
  real(rp)                :: bocod(ndime,mnodb),bomas(mnodb),borhs(mnodb)
  real(rp)                :: eucta,baloc(ndime,ndime)
  real(rp)                :: aux
  character(15)           :: chnod

  if( kfl_delta /= 1 ) return

  call livinf(0_ip,'COMPUTE WALL DISTANCE',0_ip)

  if( INOTMASTER ) then

     if( itask == 1 ) call domain_memory_allocate('VARIABLE WALL DISTANCE') ! YWALB(NBOUN) and YWALP(NBOPO)

     if ( kfl_noslw_ker /= 0_ip)  ywale = huge(1.0_rp)  ! Init
     
     do iboun = 1,nboun
        iblty = ltypb(iboun)
        pnodb = nnode(iblty)
        ielem = lelbo(iboun)
        pnode = nnode(ltype(ielem))
        call extcog(iboun,ielem,pnodb,pnode,bonor)   ! obtain bonor: the normalized normal to the boundary 
        !
        !coordinates of element cog & boundary cog
        !
        do idime = 1,ndime
           elcog(idime)   = 0.0_rp
           bocog(idime)   = 0.0_rp
        end do
        do inode = 1,pnode  ! gather elcod
           ipoin = lnods(inode,ielem)
           do idime = 1,ndime
              elcog(idime)   = elcog(idime) + coord(idime,ipoin)
           end do
        end do
        do inodb = 1,pnodb  ! gather elcod
           ipoin = lnodb(inodb,iboun)
           do idime = 1,ndime
              bocog(idime)   = bocog(idime) + coord(idime,ipoin)
           end do
        end do
        do idime = 1,ndime
           elcog(idime)   = elcog(idime)/real(pnode,rp)
           bocog(idime)   = bocog(idime)/real(pnodb,rp)
        end do
        ! Obtain ywalb = abs ( (cogbo-cogel).bonor )
        aux = 0.0_rp      
        do idime = 1,ndime
           aux   = aux + bonor(idime) * (elcog(idime)-bocog(idime)) 
        end do
        ywalb(iboun) = delmu_dom * abs(aux)

        if ( kfl_noslw_ker /= 0_ip)  ywale(ielem) = min(ywale(ielem),ywalb(iboun))  ! another option might have been average - perhaps take into account only boundaries type 23
        
     end do

     !
     ! Proyect ywalb to obtain ywalp 
     !

     call memgen(zero,2_ip,max(1_ip,npoin))  
     !! allocates gevec  of size(2,nbopo) , the first one I will use it for the boundary mass matrix and the second one for the rhs
     ! I need gevec to be of size 2,nbopo and not nbopo,2 so that I can use rhsmod
     ! I do not see any need to use pararr('SLX'  If I am using rhsmod instead , discuss with guillaume 
     !
     ! Loop over boundaries to contruct lumped boundary mass matrix and RHS
     !
     if( ( kfl_naxis == 1 ).or.( kfl_spher == 1 )) call runend('waldis/massmb: axis and spher not ready')
     boundaries: do iboun=1,nboun

        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)

        if( ndime == 1 ) then
           do inodb=1,pnodb
              ipoin          = lnodb(inodb,iboun)
              bocod(1,inodb) = coord(1,ipoin) 
           end do
        else if( ndime == 2 ) then
           do inodb=1,pnodb
              ipoin          = lnodb(inodb,iboun)
              bocod(1,inodb) = coord(1,ipoin)
              bocod(2,inodb) = coord(2,ipoin)
           end do
        else
           do inodb = 1,pnodb
              ipoin          = lnodb(inodb,iboun)
              bocod(1,inodb) = coord(1,ipoin)
              bocod(2,inodb) = coord(2,ipoin)
              bocod(3,inodb) = coord(3,ipoin)
           end do
        end if
        !
        ! Loop over Gauss points (which are nodes)
        !
        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deric(:,:,inodb),&
                bocod,baloc,eucta)    ! we don't need baloc here, just eucta 
           aux   = elmar(pblty)%weigc(inodb)*eucta
           bomas(inodb) = aux
           borhs(inodb) = aux * ywalb(iboun)
        end do
        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           ibopo = lpoty(ipoin)
           gevec(1,ipoin) = gevec(1,ipoin) + bomas(inodb)  ! boundary mass matrix
           gevec(2,ipoin) = gevec(2,ipoin) + borhs(inodb)  ! boundary rhs
        end do

     end do boundaries
     !
     ! Modify RHS due to periodicity and Parall service
     !
     call rhsmod(2_ip,gevec)
     !
     ! Hanging nodes
     ! 
     if( nhang > 0 ) call runend('MASSMB: NOT PROGRAMMED')
     !
     ! Loop over boundary nodes to control zero-volume points
     !
     ichkm = 0
     aux   = epsilon(1.0_rp)
     boundary_nodes: do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if( ibopo > 0 ) then
           if( gevec(1,ipoin) < aux ) then
              ichkm = 1
              chnod = intost(ibopo)
              call livinf(10000_ip, &
                   'BOUNDARY NODE NUMBER ( '//adjustl(trim(chnod)) &
                   // ' HAS ZERO MASS)',0_ip) 
              if(1==1) then
                 chnod = intost(ipoin)
                 print*,'ipoin,kfl_paral,ipoin_glob',ipoin,kfl_paral,lninv_loc(ipoin)
                 call livinf(10000_ip, &
                      'CORRESPONDS TO POINT ( '//adjustl(trim(chnod)),0_ip)
              end if
           end if
        end if
     end do boundary_nodes

     if( ichkm == 1 ) call runend('WALDIS: ZERO MASS BOUNDARY NODES FOUND')
     
     do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if( ibopo > 0 ) ywalp(ibopo) = gevec(2,ipoin)/gevec(1,ipoin)
     end do

     call memgen(two,2_ip,max(1_ip,npoin))  ! deallocates gevec

     !
     ! Check that no boundary has ywalp == 0, that would give divide by 0 in tur_updbcs w=U*/(sqrt(beta*)*kap*y)
     !
     do ibopo = 1,nbopo
        if (ywalp(ibopo)<1.0e-12_rp) then
           do ipoin=1,npoin
              if (lpoty(ipoin)==ibopo) print*,'ibopo,ipoin,ywalp(ibopo)',ibopo,ipoin,ywalp(ibopo)
           end do
        end if
     end do

  end if   ! not master  ! perhaps better use a return if master

end subroutine waldis



!  real(rp)    :: auxve(ndime),elno2,elno1
!  real(rp)    :: tragl(ndime,ndime),hleng(ndime)  ! Element characteristics
!        guillaume says that this way of calculatin element length in a certain direction does not work fine. 
!        Instead I will calculate distance in the normal dierction as ywalb = abs ( (cogbo-cogel).bonor)
!        I will conserve it in case I want to go back to it
!        !
!        ! HLENG and TRAGL at center of gravity
!        !
!        call elmlen(&
!             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
!             hnatu(pelty),hleng)
!        !
!        ! Compute the length in the direction normal to the boundary
!        !
!        call mbvab1(auxve,tragl,bonor,ndime,ndime,elno2,elno1)   ! borrowed from elmchl in velocity direction
!        if(elno2>1.0e-16.and.elno1>1.0e-16) then
!           ywalb(iboun) = hnatu*elno1/elno2
!        else
!           ywalb(iboun) = hleng(ndime)  ! minimum elememt length , find better option   
!        end if
!        ywalb(iboun) = ywalb(iboun)/2.0_rp
