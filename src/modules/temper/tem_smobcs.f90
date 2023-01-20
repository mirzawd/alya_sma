!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_smobcs()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_smobcs
  ! NAME 
  !    tem_smobcs
  ! DESCRIPTION
  !    This routine smoothes the boundary conditions
  ! USED BY
  !    tem_iniunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
!  integer(ip) :: ipoin

!  if( INOTMASTER ) then
!     !
!     ! Initialize solver
!     !
!     call tem_inisol()
!     if( INOTMASTER ) then
!        call inisol()
!        do ipoin = 1,npoin
!           unkno(ipoin) = 0.0_rp
!        end do        
!        do ielem = 1,nelem
!           pelty = ltype(ielem)
!           pnode = nnode(pelty)
!           pgaus = ngaus(pelty)
!           do inode = 1,pnode
!              ipoin = lnods(inode,ielem)
!              do idime = 1,ndime
!                 elcod(idime,inode) = coord(idime,ipoin)
!              end do
!           end do
!           do igaus = 1,pgaus     
!              call elmder(&
!                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
!                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
!              gpvol(igaus)=elmar(pelty)%weigp(igaus)*gpdet  
!           end do
!           !
!           ! Compute element matrix ELMAL and assemble LAPLA_TUR
!           !
!           call tur_elmlap(&
!                pnode,pgaus,lnods(1,ielem),gpcar,elmar(pelty)%shape,&
!                gpvol,elmal,elrhs)
!           call assmat(&
!                1_ip,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
!                lnods(1,ielem),elmal,amatr)
!           call assrhs(1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
!        end do
!     end if
!     !
!     ! Solve system: Lapl(f)=-1, with f=0 on wall
!     !
!     call solver(rhsid,unkno,amatr,pmatr) 
!     !
!     ! Compute wall distance: d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
!     !
!     if( INOTMASTER ) then
!        call memgen(zero,ndime,npoin)
!        call gradie(unkno,gevec)
!        do ipoin=1,npoin
!           fact1=0.0_rp 
!           do idime=1,ndime
!              fact1=fact1+gevec(idime,ipoin)*gevec(idime,ipoin)
!           end do       
!           fact2=fact1+2.0_rp*max(unkno(ipoin),0.0_rp)
!           if(fact2<0.0_rp) then
!              call runend('WRONG DISTANCE TO THE WALL')
!           else
!              walld_tur(ipoin)=sqrt(fact2)-sqrt(fact1)+delta_tur
!           end if
!        end do
!        call memgen(two,ndime,npoin)
!     end if
!  end if

end subroutine tem_smobcs
