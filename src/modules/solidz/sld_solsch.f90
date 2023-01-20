!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_solshe()
  use def_kintyp,    only : ip,rp
  use def_elmtyp,    only : SHELL,BAR3D  
  use def_domain,    only : lnods,coord,nelem,lexis
  use def_domain,    only : ltype,npoin,ndime,nnode
  use def_master,    only : solve,rhsid,amatr,unkno
  use def_master,    only : solve_sol,displ,gevec
  use mod_solver,    only : solver_solve
  use mod_sld_shell, only : STIFCSTDKT 
  use mod_sld_shell, only : STIFPOR3D
  use mod_postpr,    only : postpr_right_now
  implicit none
  integer(ip) :: ielem,inode,ipoin,pelty
  integer(ip) :: pevat,pnode,idime,idofn
  real(rp)    :: elmat(18,18)
  !real(rp)    :: elmat(12,12)
  real(rp)    :: elrhs(6)
  real(rp)    :: elcod(3,3)
  integer(ip) :: EST
  real(rp)    :: E,POI,T
  real(rp)    :: JT,IY,IZ,AREA,L

  if( lexis(SHELL) == 1 .or. lexis(BAR3D) == 1 ) then

     solve_sol => solve(2:)
     call inisol()

     amatr = 0.0_rp
     EST   = 1
     E     = 3.0e+7_rp
     POI   = 0.3_rp
     T     = 0.5_rp
     elrhs = 0.0_rp

     E     = 2.8e10_rp  
     POI   = 0.2_rp
     JT    = 2e-3_rp
     IY    = 2.0833e-3_rp
     IZ    = 1.0_rp
     AREA  = 0.1_rp

     do ielem = 1,nelem

        pelty = ltype(ielem)

        if( pelty == SHELL .or. pelty == BAR3D ) then

           pnode = nnode(pelty)
           pevat = 6 * pnode
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(inode,1:3) = coord(1:3,ipoin)
           end do
           elmat = 0.0_rp

           if( ltype(ielem) == SHELL ) then
              call STIFCSTDKT(elcod,E,POI,T,EST,elmat)
           else if( ltype(ielem) == BAR3D ) then
              L = 0.0_rp
              do idime = 1,2
                 L = L + (elcod(1,idime)-elcod(2,idime))**2
              end do
              L = sqrt(L)
              call STIFPOR3D(L,E,POI,AREA,JT,IY,IZ,elmat)
           end if

           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              if(  ltype(ielem) == SHELL .and. &
                   ( ipoin == 1 .or. &
                   ipoin == 3 .or. &
                   ipoin == 4 .or. &
                   ipoin == 11 .or. &
                   ipoin == 12 .or. &
                   ipoin == 16 .or. &
                   ipoin == 28 .or. &
                   ipoin == 31 .or. &
                   ipoin == 35 .or. &
                   ipoin == 46 ) ) then
                 do idime = 1,6
                    idofn = (inode-1)*6 + idime
                    elmat(idofn,:)     = 0.0_rp
                    elmat(idofn,idofn) = 1.0_rp
                 end do
              else if(  ltype(ielem) == BAR3D .and. ipoin == 1 ) then
                 do idime = 1,3
                    idofn = (inode-1)*6 + idime
                    elmat(idofn,:)     = 0.0_rp
                    elmat(idofn,idofn) = 1.0_rp
                 end do
              else if(  ltype(ielem) == BAR3D .and. ipoin == 101 ) then
                 do idime = 3,3
                    idofn = (inode-1)*6 + idime
                    elmat(idofn,:)     = 0.0_rp
                    elmat(idofn,idofn) = 1.0_rp
                 end do
              end if
           end do

           call assmat(&
                solve(2) % ndofn,pnode,pevat,solve(2) % nunkn,&
                solve(2) % kfl_algso,ielem,lnods(1,ielem),elmat,amatr)

        end if
     end do

     do ipoin = 1,6*npoin
        rhsid(ipoin) = 0.0_rp
        unkno(ipoin) = 0.0_rp
     end do
!print*,'a'
     ipoin = 154; rhsid( (ipoin-1)*6 + 2 ) = -5000.0_rp
     ipoin = 160; rhsid( (ipoin-1)*6 + 2 ) = -5000.0_rp

     !ipoin = 50; rhsid( (ipoin-1)*6 + 3 ) = -6000.0_rp

!print*,'b'
     call solver_solve(solve_sol,amatr,rhsid,unkno)  
!print*,'c'
    
     do ipoin = 1,npoin
        do idime = 1,ndime
           displ(idime,ipoin,1) = unkno( (ipoin-1)*6+idime )
        end do
     end do 
     gevec => displ(:,:,1) 
     call postpr_right_now('XXXXX','VECTO','NPOIN',gevec)
stop
     
  end if

end subroutine sld_solshe
 
