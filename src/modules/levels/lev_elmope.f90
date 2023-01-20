!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmope(order)

  !------------------------------------------------------------------------
  !****f* Levels/lev_elmope
  ! NAME 
  !    lev_elmope
  ! DESCRIPTION
  !    ORDER=1:
  ! Level Set convection equation
  !      1. Compute elemental matrix 
  !      2. Assemble them
  !    ORDER=2:
  !      Update the subgrid scale
  !      1. Compute rhs
  ! 
  ! USES
  !    lev_elmdir
  !    lev_elmmat 
  ! USED BY
  !    lev_matrix
  !------------------------------------------------------------------------

  use      def_parame
  use      def_elmtyp
  use      def_master
  use      def_kermod
  use      def_domain
  use      def_levels
  use      def_solver
  implicit none

  integer(ip), intent(in) :: order                          ! 
  integer(ip)             :: ielem,igaus,idime,itime        ! Indices and dimensions
  integer(ip)             :: inode,jnode,ipoin
  integer(ip)             :: pelty,pnode,pgaus,plapl,porde
  real(rp)                :: elmat(mnode,mnode)             ! Element matrices
  real(rp)                :: elrhs(mnode)
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode,ncomp_lev),elle0(mnode)
  real(rp)                :: elvel(ndime,mnode)

  real(rp)                :: xjaci(ndime,ndime) 
  real(rp)                :: xjacm(ndime,ndime) 


  real(rp)                :: tragl(ndime,ndime),chave(ndime,2)     ! Stabilization
  real(rp)                :: chale(2)
  real(rp)                :: gpnve,hleng(ndime),taust  

  real(rp)                :: gpcar(ndime,mnode,mgaus)                ! dNk/dxj
  real(rp)                :: gpvel(3),gpgrl(ndime)                   ! u,grad(phi)
  real(rp)                :: gpvol,gpdet                             ! |J|*w,|J|
  real(rp)                :: gplev(ncomp_lev),gple0                  ! phi 
  real(rp)                :: fact1,fact2,fact3
  real(rp)                :: tleve,dtn,gpdiv,dummr

  real(rp)                :: thire,almti,almtj
  integer(ip)             :: islim


  select case(order)
     !
     ! Explicit treatment case 
     !
  case(1_ip)
     !
     ! Left-hand side
     !
     dtn   = 1.0_rp/dtinv_lev

     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty=ltype(ielem)
        if( pelty > 0 ) then
           pnode=nnode(pelty)
           pgaus=ngaus(pelty)
           plapl=0
           porde=lorde(pelty)
           do inode=1,pnode
              do jnode=1,pnode
                 elmat(inode,jnode)=0.0_rp
              end do
           end do
           !
           ! Gather operations
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode,1)=fleve(ipoin,1)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do

           call lev_velfun(kfl_advec_lev,ndime,pnode,lnods(1,ielem),elcod,elvel)
           !
           ! Mesh velocity
           !     
           if( kfl_coupl(ID_LEVELS,ID_ALEFOR) >= 1 ) then
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elvel(idime,inode) = elvel(idime,inode) - velom(idime,ipoin)
                 end do
              end do
           end if          
           !
           ! hleng and tragl at center of gravity
           !
           call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                hnatu(pelty),hleng)
           !
           ! Compute the characteristic length CHALE
           !
           call elmchl(&
                tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                porde,hnatu(pelty),kfl_advec_lev,kfl_ellen_lev)

           do igaus=1,pgaus
              !
              ! Cartesian derivatives, and volume: GPCAR, PGVOL
              !
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
              gpvol=elmar(pelty)%weigp(igaus)*gpdet

              gpnve=0.0_rp
              !
              ! Gauss point values
              !
              do idime=1,ndime
                 gpvel(idime)=0.0_rp
                 gpgrl(idime)=0.0_rp
              end do

              gplev(1) = 0.0_rp
              gplev(2) = 0.0_rp

              do inode=1,pnode
                 ipoin=lnods(inode,ielem)              
                 do idime=1,ndime
                    gpvel(idime)=gpvel(idime)+elmar(pelty)%shape(inode,igaus)*elvel(idime,inode)
                    gpgrl(idime)=gpgrl(idime)+gpcar(idime,inode,igaus)*ellev(inode,1)
                 end do
                 if(kfl_tiacc_lev==1) then
                    if(itinn(modul)<=1) then   
                       gplev(1)=gplev(1)+elmar(pelty)%shape(inode,igaus)*fleve(ipoin,3)
                       gplev(2)=gplev(2)+elmar(pelty)%shape(inode,igaus)*fleve(ipoin,4)
                    else
                       gplev(1)=gplev(1)+elmar(pelty)%shape(inode,igaus)*fleve(ipoin,1)
                       gplev(2)=gplev(2)+elmar(pelty)%shape(inode,igaus)*fleve(ipoin,3)
                    endif
                 else if(kfl_tiacc_lev==2.and.supgp_lev==1_rp) then
                    gplev(1)=gplev(1)+elmar(pelty)%shape(inode,igaus)*fleve(ipoin,4)
                    gplev(2)=gplev(2)+elmar(pelty)%shape(inode,igaus)*fleve(ipoin,3)
                 endif
              end do
              !
              ! Compute SUPG stabilization term
              !
              call vecnor(gpvel,ndime,gpnve,2_ip) 
              if(gpnve==0.0) then
                 taust=0.0
              else 
                 taust=chale(1)/(2.0_rp*gpnve)
              endif

              !
              ! Extension elements
              !
              if( lelch(ielem) == ELEXT ) then
                 call runend('NOT CODED')
                 call elmext(&
                      4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
                      dummr,elrhs,ellev)
              end if


              do inode=1,pnode
                 ipoin=lnods(inode,ielem)              
                 tleve=0.0_rp
                 do idime=1,ndime
                    tleve=tleve+gpcar(idime,inode,igaus)*gpvel(idime)
                 end do
                 tleve=taust*tleve

                 rhsid(ipoin)=rhsid(ipoin)-(gplev(1)-gplev(2))*tleve*gpvol*supgp_lev

                 tleve=tleve+elmar(pelty)%shape(inode,igaus)
                 do idime=1,ndime
                    rhsid(ipoin)=rhsid(ipoin)-gpgrl(idime)*gpvel(idime)*dtn*tleve*gpvol
                 end do
              end do

           end do

        end if

     end do


  case(2_ip)
     !
     ! Left-hand side Implicit Treatment
     !
     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)
        pelty = abs(pelty)

        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = 0
           porde = lorde(pelty)
           !
           ! Initialization
           !
           call lev_elmini(pnode,elmat,elrhs)
           !
           ! Gather operations
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode,1)=fleve(ipoin,1)
              ellev(inode,2)=fleve(ipoin,3)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do

           if(kfl_tisch_lev==2) then
              do itime = 2,kfl_tiacc_lev
                 do inode=1,pnode
                    ipoin=lnods(inode,ielem)
                    ellev(inode,itime+1) = fleve(ipoin,itime+2)
                 enddo
              enddo
           endif

           call lev_velfun(kfl_advec_lev,ndime,pnode,lnods(1,ielem),elcod,elvel)
           !
           ! Mesh velocity
           !     
           if( kfl_coupl(ID_LEVELS,ID_ALEFOR) >= 1 ) then
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elvel(idime,inode) = elvel(idime,inode) - velom(idime,ipoin)
                 end do
              end do
           end if
           !
           ! hleng and tragl at center of gravity
           !
           call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                hnatu(pelty),hleng)
           !
           ! Compute the characteristic length CHALE
           !
           call elmchl(tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                porde,hnatu(pelty),kfl_advec_lev,kfl_ellen_lev)

           do igaus=1,pgaus
              !
              ! Cartesian derivatives, and volume: GPCAR, PGVOL
              !
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
              gpvol = elmar(pelty)%weigp(igaus)*gpdet

              fact1 = 0.0_rp
              fact2 = 0.0_rp
              fact3 = 0.0_rp
              gpnve = 0.0_rp
              gpvel = 0.0_rp
              gplev = 0.0_rp
              !
              ! Gauss point values
              !
              do inode = 1,pnode
                 gplev(2)       = gplev(2)       + elmar(pelty) % shape(inode,igaus) * ellev(inode,2)   
                 gpvel(1:ndime) = gpvel(1:ndime) + elmar(pelty) % shape(inode,igaus) * elvel(1:ndime,inode)
              end do
              if( kfl_tisch_lev == 2 ) then
                 do itime = 2,kfl_tiacc_lev
                    do inode = 1,pnode
                       gplev(itime+1) = gplev(itime+1) + elmar(pelty) % shape(inode,igaus) * ellev(inode,itime+1)   
                    end do
                 end do
              endif
              !
              ! Compute SUPG stabilization term
              !
              call vecnor(gpvel,ndime,gpnve,2_ip) 
              if( gpnve < zeror ) then
                 taust = 0.0_rp
              else 
                 taust = chale(1)/(2.0_rp*gpnve)
              endif

              if( kfl_timco /= 2 ) then
                 if( kfl_tiacc_lev == 1 ) then
                    fact1 = gpvol * dtinv_lev
                 else if( kfl_tiacc_lev == 2 ) then
                    !
                    ! Crank-Nicolson method 
                    !        
                    if( kfl_tisch_lev == 1 ) then 
                       fact1 = gpvol * dtinv_lev
                       !
                       ! BDF scheme
                       !
                    else if( kfl_tisch_lev == 2 )then
                       fact1 = gpvol * dtinv_lev
                    end if
                 end if
              else
                 if( taust /= 0.0_rp ) then
                    fact1 = gpvol/(taust*safet_lev)  ! check if it is ok for CN & bdf2
                 else
                    fact1 = gpvol*dtinv_lev
                 end if
              end if
              gpdiv = 0.0_rp
              !
              ! Assembling in element matrix
              !
              call lev_elmmat(&
                   pnode,elmar(pelty)%shape(1,igaus),gpcar(1,1,igaus),&
                   gpvel,gplev,gpvol,taust,fact1,elmat,elrhs)

           end do
           !
           ! Dirichlet condition
           !
           call lev_elmdir(pnode,lnods(1,ielem),elmat,elrhs)
           !
           ! Extension elements
           !
           if( lelch(ielem) == ELEXT ) then
              call elmext(&
                   4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
                   dummr,elrhs,ellev)
           end if
           !
           ! Assembling in domain matrices
           !
           call assrhs(&
                solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
           call assmat(&
                solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
                ielem,lnods(1,ielem),elmat,amatr)
        end if
     end do

  case(3_ip)
     !
     ! Implicit Treatment for redistanciation with Sussman equation
     !
     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = 0
           porde = lorde(pelty)
           !
           ! Initialization
           !
           call lev_elmini(pnode,elmat,elrhs)
           !
           ! Gather operations
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode,1)=fleve(ipoin,1)
              elle0(inode)=flev0_lev(ipoin)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
                 elvel(idime,inode)=norml_lev(idime,ipoin)
              end do
           end do
           !
           ! Mesh velocity
           !     
           if( kfl_coupl(ID_LEVELS,ID_ALEFOR) >= 1 ) then
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elvel(idime,inode) = elvel(idime,inode) - velom(idime,ipoin)
                 end do
              end do
           end if          
           !
           ! hleng and tragl at center of gravity
           !
           call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                hnatu(pelty),hleng)
           !
           ! Compute the characteristic length CHALE
           !
           call elmchl(tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                porde,hnatu(pelty),kfl_advec_lev,kfl_ellen_lev)

           do igaus=1,pgaus
              !
              ! Cartesian derivatives, and volume: GPCAR, PGVOL
              !
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
              gpvol=elmar(pelty)%weigp(igaus)*gpdet

              fact1    = 0.0_rp
              fact2    = 0.0_rp
              fact3    = 0.0_rp
              gpnve    = 0.0_rp
              gpvel(1) = 0.0_rp
              gpvel(2) = 0.0_rp
              gplev(2) = 0.0_rp
              gple0    = 0.0_rp
              !
              ! Gauss point values
              !
              if(ndime==2) then
                 do inode=1,pnode
                    gplev(2)=gplev(2)&
                         +elmar(pelty)%shape(inode,igaus)*ellev(inode,1)   
                    gple0 = gple0&
                         +elmar(pelty)%shape(inode,igaus)*elle0(inode)   
                    gpvel(1)=gpvel(1)&
                         +elmar(pelty)%shape(inode,igaus)*elvel(1,inode)
                    gpvel(2)=gpvel(2)&
                         +elmar(pelty)%shape(inode,igaus)*elvel(2,inode)
                 end do
              else
                 gpvel(3)=0.0_rp
                 do inode=1,pnode
                    gplev(2)=gplev(2)&
                         +elmar(pelty)%shape(inode,igaus)*ellev(inode,1)   
                    gple0 = gple0&
                         +elmar(pelty)%shape(inode,igaus)*elle0(inode)   
                    gpvel(1)=gpvel(1)&
                         +elmar(pelty)%shape(inode,igaus)*elvel(1,inode)
                    gpvel(2)=gpvel(2)&
                         +elmar(pelty)%shape(inode,igaus)*elvel(2,inode)
                    gpvel(3)=gpvel(3)&
                         +elmar(pelty)%shape(inode,igaus)*elvel(3,inode)
                 end do
              end if
              !
              ! Compute SUPG stabilization term
              !
              call vecnor(gpvel,ndime,gpnve,2_ip) 
              if( gpnve == 0.0_rp ) then
                 taust = 0.0_rp
              else 
                 taust = chale(1)/(2.0_rp*gpnve)
              endif

              fact1 = gpvol*dtred_lev
              gpdiv = 0.0_rp

              ! Assembling in element matrice
              call lev_elmma2(&
                   pnode,elmar(pelty)%shape(1,igaus),gpcar(1,1,igaus),&
                   gpvel,gplev,gple0,gpvol,taust,fact1,elmat,elrhs)

           end do
           !
           ! Decide elements not to assemble when kfl_locre_lev>1 , also modify icupt_lev
           !
           if(kfl_locre_lev>0) then
              islim=0    ! the element belong to the limit that devides nodes that are redistanced from those that are not 
              thire = 10.0_rp*thicl  ! I should add thire in def_levels
              do inode = 1,pnode-1
                 almti = abs(elle0(inode))-thire   ! Abs Level Minus Thire in node Inode
                 if (almti==0.0_rp) almti = -0.0001_rp*thire  ! solution to avoid nodes exactly at thire
                 do jnode = inode,pnode
                    almtj = abs(elle0(inode))-thire   ! Abs Level Minus Thire in node Inode
                    if (almtj==0.0_rp) almtj = -0.0001_rp*thire  ! solution to avoid nodes exactly at thire
                    if( almti * almtj < 0.0_rp) islim=1
                 end do
              end do
              if(islim==1)  elmat=0.0_rp  ! setting elmat to 0 I eliminate the contribution from elements in the boundary
              !
              ! Modify icupt_lev 
              !
              if(kfl_locre_lev>0) then
                 do inode = 1,pnode
                    if ( abs(elle0(inode)) >= thire ) icupt_lev(lnods(inode,ielem))=1
                 end do
              end if
           end if
           !
           ! Extension elements
           !
           if( lelch(ielem) == ELEXT ) then
              call elmext(&
                   4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
                   dummr,elrhs,ellev)
           end if
           !
           ! Dirichlet condition
           !
           call lev_elmdi2(pnode,lnods(1,ielem),elmat,elrhs)
           !
           ! Assembling in domain matrices
           !
           if( solve_sol(1)%kfl_algso == 9 ) then
              call assric(&
                   solve_sol(1)%ndofn,pnode,lnods(1,ielem),elrhs,elmat,ellev,rhsid)    
           else
              call assrhs(&
                   solve_sol(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
              call assmat(&
                   solve_sol(1)%ndofn,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
                   ielem,lnods(1,ielem),elmat,amatr)
           end if

        end if

     end do

  end select

end subroutine lev_elmope

!------------------------------------------------------------------------
! NOTES
! 
! Governing equation:
!
! L(phi) = d phi / dt + u.grad(phi) =0
!
! The equation is stabilized using the SUPG model:
!
! (D_k \phi^{n+1}_h /dt+u^{n+1}_h.grad(phi^{n+1}_h),v_h)
!  +(tau^{n+1} u^{n+1}_h .grad(v_h) , D_k \phi^{n+1}_h + u^{n+1}_h . \grad(phi^{n+1})_h )=0
!
! with
!                     h_e                          
! tau   =  ----------------------  
!               2|u_e|
!where h_e is the element length in the direction of the flow and
! |u_e| the velocity norm of element e.
!
!***
!------------------------------------------------------------------------
