!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_distuc()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_distuc
  ! NAME 
  !    lev_distuc
  ! DESCRIPTION
  !    Compute the generalized distance to the cut nodes using the same strategy useg by tur_walgen perhaps we should unify
  !    Compared to tur_walgen I only leave the option corresponding to imeth=0 
  !    Poisson equation:
  !    1. Solve Lapl(f)=-1, with f=0 on wall
  !    2. d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !    See the following references:
  !    P.G. Tucker, Differential equation-based wall distance computation for
  !         DES and RANS, J. Comp. Phys. 190 (2003) 229-248.
  !    P.G. Tucker, Int. J. Numer. Fluids 33 (2000) 869.
  !    P.G. Tucker, Appl. Math. Model. 22 (1998) 293.
  !
  !    It outputs the results directly into fleve(ipoin,1) 
  !
  ! USES
  ! USED BY
  !    lev_redieq
  !*** 
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  use def_solver
  use mod_gradie
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: ielem,igaus,idime,inode,ipoin,pnode,pgaus,pelty
  real(rp)    :: elmal(mnode,mnode),elrhs(mnode),elcod(ndime,mnode)
  real(rp)    :: gpcar(ndime,mnode,mgaus),gpvol(mgaus),gpdet,fact1
  real(rp)    :: xjaci(9),xjacm(9),fact2
  
  !
  ! First Extend distance on the wall to the rest of the domian smoothly to be used instead of delta_tur that is not constant
  ! Solve system: Lapl(g)=0, with g = ywalp on nodes with icupt_lev == 1
  ! For the moment I use the same solver used to calculate wall distance (next step) and I only leave the option corresponding to imeth=0 
  !
  !
  ! Initialize solver
  !
  call livinf(59_ip,'EXTEND DISTANCE FROM CUT NODES',modul)
  call lev_inisol(3_ip)

  if( INOTMASTER ) then
     call inisol()
     do ipoin = 1,npoin
        unkno(ipoin) = 0.0_rp
     end do
     do ielem = 1,nelem
        !
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        !
        ! Gather operations: ELCOD
        !
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        !
        ! 1st order Cartesian derivatives GPCAR and GPVOL=dV=|J|*wgx
        !
        do igaus = 1,pgaus     
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
           gpvol(igaus) = elmar(pelty)%weigp(igaus)*gpdet  
        end do
        !
        ! Compute element matrix ELMAL and assemble LAPLA_LEV
        !
        call lev_elmlap(&
             two,pnode,pgaus,lnods(1,ielem),gpcar,elmar(pelty)%shape,&
             gpvol,elmal,elrhs)
        call assmat(&
             1_ip,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
             ielem,lnods(1,ielem),elmal,amatr)
        call assrhs(&
             1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
     end do
  end if
  !
  ! Solve system: Lapl(g)=0, with g = flev0_lev on nodes with icupt_lev > 0
  !
  call solver(rhsid,unkno,amatr,pmatr)
 
  if( INOTMASTER ) then
     call memgen(zero,npoin,zero)         ! allocates gesca
     do ipoin = 1,npoin
        gesca(ipoin) = unkno(ipoin)       ! saves result is gesca for later use
     end do
  end if

  ! Second part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !
  ! Initialize solver
  !
  call livinf(59_ip,'COMPUTE DISTANCE TO CUT NODES',modul)
  call lev_inisol(3_ip)

  if( INOTMASTER ) then
     call inisol()
     do ipoin = 1,npoin
        unkno(ipoin) = 0.0_rp
     end do
     do ielem = 1,nelem
        !
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        !
        ! Gather operations: ELCOD
        !
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        !
        ! 1st order Cartesian derivatives GPCAR and GPVOL=dV=|J|*wgx
        !
        do igaus = 1,pgaus     
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
           gpvol(igaus) = elmar(pelty)%weigp(igaus)*gpdet  
        end do
        !
        ! Compute element matrix ELMAL and assemble LAPLA_TUR
        !
        call lev_elmlap(&
             one,pnode,pgaus,lnods(1,ielem),gpcar,elmar(pelty)%shape,&
             gpvol,elmal,elrhs)
        call assmat(&
             1_ip,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
             ielem,lnods(1,ielem),elmal,amatr)
        call assrhs(&
             1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
     end do
  end if
  !
  ! Solve system: Lapl(f)=-1, with f=0 on wall
  !
  call solver(rhsid,unkno,amatr,pmatr) 
  !
  ! Compute wall distance: d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !
  if( INOTMASTER ) then
     call memgen(zero,ndime,npoin)
     call gradie(unkno,gevec)
     do ipoin = 1,npoin
        fact1 = 0.0_rp 
        do idime = 1,ndime
           fact1 = fact1 + gevec(idime,ipoin) * gevec(idime,ipoin)
        end do
        fact2 = fact1 + 2.0_rp*max(unkno(ipoin),0.0_rp)
        if( fact2 < 0.0_rp ) then
           call runend('WRONG DISTANCE CALCULATION IN LEV_DISTUC')
        else
           if(flev0_lev(ipoin)>0.0_rp) then
              fleve(ipoin,1) = sqrt(fact2) - sqrt(fact1) + gesca(ipoin)
           else
              fleve(ipoin,1) = -( sqrt(fact2) - sqrt(fact1) + gesca(ipoin) )
           end if
        end if
     end do
     call memgen(two,ndime,npoin)
  end if

  if ( INOTMASTER ) call memgen(two,npoin,zero)         ! deallocates gesca

end subroutine lev_distuc
