!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_assemble_schur.f90
!> @author  Guillaume Houzeaux
!> @brief   Matrix and RHS assembly
!> @details Assembly of matrix and RHS:
!>          1. Auu, Aup, Apu, App, bu, bp
!>          2. Only App 
!>          3. CMM consitent mass matrix 
!> @} 
!------------------------------------------------------------------------

subroutine nsi_assemble_schur(&
     itask,pnode,pevat,ielem,lnods,elauu,elaup,elapp,elapu,&
     elrbu,elrbp,Auu,Aup,App,Apu,bu,bp)
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  nzdom,ndime,r_sol
  use def_domain, only       :  r_sym,c_sym,c_sol,lezdo
  use def_master, only       :  solve
  use def_kermod, only       :  kfl_element_to_csr
  implicit none
  integer(ip), intent(in)    :: itask,pnode,pevat,ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elauu(pevat,pevat)
  real(rp),    intent(in)    :: elaup(pevat,pnode)
  real(rp),    intent(in)    :: elapp(pnode,pnode)
  real(rp),    intent(in)    :: elapu(pnode,pevat)
  real(rp),    intent(in)    :: elrbu(ndime,pnode)
  real(rp),    intent(in)    :: elrbp(pnode)
  real(rp),    intent(inout) :: Auu(ndime,ndime,nzdom) 
  real(rp),    intent(inout) :: Aup(ndime,nzdom)
  real(rp),    intent(inout) :: Apu(ndime,nzdom)
  real(rp),    intent(inout) :: App(solve(2)%nzmat)
  real(rp),    intent(inout) :: bu(ndime,*),bp(*) 
  integer(ip)                :: ndofn,inode,jnode,iposi,jposi,izsym 
  integer(ip)                :: idime,jdime,izsol,jpoin,ipoin,jcolu

  ndofn = ndime + 1

  select case(itask)

  case (1_ip)

     !----------------------------------------------------------------
     !
     ! Assemble the following matrices:
     !
     ! Auu <= Auu^(e)
     ! Aup <= Aup^(e)
     ! Apu <= Apu^(e)
     ! App <= App^(e)
     ! bu  <= bu^(e)
     ! bp  <= bp^(e)
     !
     !----------------------------------------------------------------

     if( solve(2) % kfl_symme == 1 ) then
        !
        ! App requires symmetric assembly
        !
        do inode = 1,pnode
           ipoin = lnods(inode)
           do jnode = 1,pnode
              jpoin = lnods(jnode)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jcolu == jpoin ) then
                 do idime=1,ndime
                    iposi=(inode-1)*ndime+idime
                    do jdime=1,ndime
                       jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                            &                   + elauu(iposi,jposi)          ! Auu
                    end do
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode)  ! Aup
                 end do
                 iposi = iposi+1
                 do jdime=1,ndime
                    jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu                 
                 end do
              end if
              if( jpoin <= ipoin ) then 
                 izsym = r_sym(ipoin)
                 jcolu = c_sym(izsym)
                 do while( jcolu /= jpoin .and. izsym < r_sym(ipoin+1)-1 )
                    izsym = izsym + 1
                    jcolu = c_sym(izsym)
                 end do
                 if( jcolu == jpoin ) then
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    App(izsym) = App(izsym) + elapp(inode,jnode)              ! App
                 end if
              end if
           end do
        end do
     else
        !
        ! App requires assymmetric assembly
        !
        if( kfl_element_to_csr == 1 ) then
           do inode = 1,pnode
              do jnode = 1,pnode
                 izsol = lezdo(inode,jnode,ielem)
                 do idime = 1,ndime
                    iposi = (inode-1)*ndime+idime
                    do jdime = 1,ndime
                       jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol)&
                            +                   elauu(iposi,jposi)           ! Auu
                    end do
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode) ! Aup
                 end do
                 do jdime=1,ndime
                    jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu                 
                 end do
#ifdef NO_COLORING
                 !$OMP ATOMIC
#endif
                 App(izsol) = App(izsol) + elapp(inode,jnode)                 ! App
              end do
           end do
        else
           do inode = 1,pnode
              ipoin = lnods(inode)
              do jnode = 1,pnode
                 jpoin = lnods(jnode)
                 izsol = r_sol(ipoin)
                 jcolu = c_sol(izsol)
                 do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                    izsol = izsol + 1
                    jcolu = c_sol(izsol)
                 end do
                 if( jcolu == jpoin ) then
                    do idime = 1,ndime
                       iposi = (inode-1)*ndime+idime
                       do jdime = 1,ndime
                          jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                          !$OMP ATOMIC
#endif
                          Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol)&
                               +                   elauu(iposi,jposi)           ! Auu
                       end do
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       Aup(idime,izsol) = Aup(idime,izsol) + elaup(iposi,jnode) ! Aup
                    end do

                    do jdime=1,ndime
                       jposi = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       Apu(jdime,izsol) = Apu(jdime,izsol) + elapu(inode,jposi)  ! Apu                 
                    end do
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    App(izsol) = App(izsol) + elapp(inode,jnode)                 ! App
                 end if
              end do
           end do
        end if
     end if
     !
     ! RHS
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
#ifdef NO_COLORING
           !$OMP ATOMIC
#endif
           bu(idime,ipoin) = bu(idime,ipoin) + elrbu(idime,inode)
        end do
#ifdef NO_COLORING
        !$OMP ATOMIC
#endif
        bp(ipoin) = bp(ipoin) + elrbp(inode)
     end do

  case (2_ip)

     !----------------------------------------------------------------
     !
     ! Assemble only (App is pressure Schur preconditioner):
     ! App <= App^(e)
     !
     !----------------------------------------------------------------

     if( solve(2) % kfl_symme == 1 ) then
        !
        ! App requires symmetric assembly
        !
        do inode = 1,pnode
           ipoin = lnods(inode)
           do jnode = 1,pnode
              jpoin = lnods(jnode)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jpoin <= ipoin ) then 
                 izsym = r_sym(ipoin)
                 jcolu = c_sym(izsym)
                 do while( jcolu /= jpoin .and. izsym < r_sym(ipoin+1)-1 )
                    izsym = izsym + 1
                    jcolu = c_sym(izsym)
                 end do
                 if( jcolu == jpoin ) then
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    App(izsym) = App(izsym) + elapp(inode,jnode) ! App
                 end if
              end if
           end do
        end do
     else
        !
        ! App requires asymmetric assembly
        !
        if( kfl_element_to_csr == 1 ) then
           do inode = 1,pnode
              do jnode = 1,pnode
                 izsol = lezdo(inode,jnode,ielem)
#ifdef NO_COLORING
                 !$OMP ATOMIC
#endif
                 App(izsol) = App(izsol) + elapp(inode,jnode) ! App
              end do
           end do
        else
           do inode = 1,pnode
              ipoin = lnods(inode)
              do jnode = 1,pnode
                 jpoin = lnods(jnode)
                 izsol = r_sol(ipoin)
                 jcolu = c_sol(izsol)
                 do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                    izsol = izsol + 1
                    jcolu = c_sol(izsol)
                 end do
                 if( jcolu == jpoin ) then
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    App(izsol) = App(izsol) + elapp(inode,jnode) ! App
                 end if
              end do
           end do
        end if

     end if

  case (3_ip)

     !----------------------------------------------------------------
     !
     ! Assemble only Auu <= Auu^(e)  - Actually I will use it to assemble the consistent mass matrix that has the sam shape
     !
     !----------------------------------------------------------------

     do inode = 1,pnode
        ipoin = lnods(inode)
        do jnode = 1,pnode
           jpoin = lnods(jnode)
           izsol = r_sol(ipoin)
           jcolu = c_sol(izsol)
           do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
              izsol = izsol + 1
              jcolu = c_sol(izsol)
           end do
           if( jcolu == jpoin ) then
              do idime=1,ndime
                 iposi=(inode-1)*ndime+idime
                 do jdime=1,ndime
                    jposi                  = (jnode-1)*ndime+jdime
#ifdef NO_COLORING
                    !$OMP ATOMIC
#endif
                    Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) &
                         &                   + elauu(iposi,jposi)          ! Auu
                 end do
              end do

           end if
        end do
     end do

  end select

end subroutine nsi_assemble_schur
