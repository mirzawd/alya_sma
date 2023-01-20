!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_solite
  !-----------------------------------------------------------------------
  !****f* Levels/lev_solite
  ! NAME 
  !    lev_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the temperature equations.
  ! USES
  !    lev_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    lev_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_levels
  use mod_postpr
  use def_solver
  implicit none
  integer(ip) :: ipoin,ielem
  !
  ! Update inner iteration counter
  !
  itinn(modul) = itinn(modul) + 1

  do ielem = 1,nelem
     lelch(ielem) = ELFEM
     ltype(ielem) = abs(ltype(ielem))
  end do
  do ipoin = 1,npoin
     lnoch(ipoin) = NOFEM
  end do
  !
  ! Construct the system matrix and right-hand-side
  !
  if( INOTMASTER ) call lev_matrix()
  !
  ! Solve the algebraic system.
  !
  if(kfl_timet_lev==1) then
     if( INOTMASTER ) call lev_explic()
  else
     call solver(rhsid,unkno,amatr,pmatr)
  end if
  !
  ! Cut elements
  !
  !call lev_cut_elements()
  
end subroutine lev_solite


subroutine lev_cut_elements()
  use def_kintyp,         only : ip,rp
  use def_master
  use def_domain
  use def_elmtyp
  use mod_cutele
  use def_kermod,         only : kfl_cutel 
  use mod_ker_updpro,     only : ker_updpro
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  integer(ip) :: ielem,inode,ipoin
  integer(ip) :: knode,pnode,jpoin
  integer(ip) :: jnode,pelty
  real(rp)    :: rsign,norma(3)
  real(rp)    :: fi,fj,t,xx(ndime,mnode)
  real(rp)    :: pcoor(3),dummr

  if( INOTMASTER .and. kfl_cutel /= 0 ) then ! TESTEO
     !
     ! Initialize hole and cut elements
     !
     do ielem = 1,nelem
        if( lelch(ielem) == ELHOL ) then
           lelch(ielem) = ELFEM
           ltype(ielem) = abs(ltype(ielem))
        else if( lelch(ielem) == ELCUT ) then
           lelch(ielem) = ELFEM       
        end if
     end do
     !
     ! Detect elements to be cut
     !
     do ielem = 1,nelem
        ipoin = lnods(1,ielem)
        rsign = fleve(ipoin,1)
        loop_inode: do inode = 2,lnnod(ielem)
           ipoin = lnods(inode,ielem)
           if( rsign*fleve(ipoin,1) < 0.0_rp ) then
              lelch(ielem) = ELCUT
              exit loop_inode
           end if
        end do loop_inode
     end do
     ! 
     ! Cut elements
     !
     do ielem = 1,nelem
        if( lelch(ielem) == ELCUT ) then       
           knode = 0
           pnode = lnnod(ielem)
           pelty = ltype(ielem)
           do inode = 1,pnode
              jnode = lenex(inode,pelty)
              ipoin = lnods(inode,ielem)
              jpoin = lnods(jnode,ielem)
              fi    = fleve(ipoin,1)
              fj    = fleve(jpoin,1)
              t     = -fi / ( fj - fi + zeror)
              if( t >= 0.0_rp .and. t <= 1.0_rp ) then
                 knode = knode + 1
                 xx(1:ndime,knode) = coord(1:ndime,ipoin) + t * ( coord(1:ndime,jpoin) - coord(1:ndime,ipoin) )
              end if
           end do
           if( knode /= 2 ) then
              call runend('PROBLEM CUTTING ELEMENT')
           else
              norma(1) = - ( xx(2,2) - xx(2,1) ) 
              norma(2) =   ( xx(1,2) - xx(1,1) ) 
              call vecuni(ndime,norma,dummr)
              pcoor = 0.0_rp
              do inode = 1,knode
                 pcoor(1:ndime) = pcoor(1:ndime) + xx(1:ndime,knode)
              end do
              pcoor(1:ndime) = pcoor(1:ndime) / real(knode,rp)
              call onecut(-1_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
           end if
        else
           knode = 0
           pnode = lnnod(ielem)
           pelty = ltype(ielem)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              if( fleve(ipoin,1) < 0.0_rp ) knode = knode + 1
           end do
           if( knode == pnode ) then
              lelch(ielem) = ELHOL
              ltype(ielem) = -ltype(ielem)
           end if
        end if
     end do
     !
     ! Update properties, mainly for cut elements
     !
     call ker_updpro()
     !
     ! Detect hole nodes
     !
     !
     ! Define LNOCH according to LELCH and detect holes automatically
     !
     call memgen(1_ip,npoin_2,0_ip)
     do ipoin = 1,npoin
        lnoch(ipoin) = NOFEM
     end do
     do ielem = 1,nelem
        if( lelch(ielem) /= ELHOL ) then
           do inode = 1,nnode(abs(ltype(ielem)))
              ipoin = lnods(inode,ielem)
              gisca(ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
     do ipoin = 1,npoin
        if( gisca(ipoin) == 0 ) then
           lnoch(ipoin) = NOHOL
        end if
     end do
     call memgen(3_ip,npoin_2,0_ip)

  end if

end subroutine lev_cut_elements
