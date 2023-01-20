!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_assmat(&
     itask,pnode,pevat,lnods,elmat,Auu,Aup,Apu,App)
  !-----------------------------------------------------------------------
  !****f* mathru/assmat
  ! NAME 
  !    assmat
  ! DESCRIPTION
  !    Assembly an elemental matrix ELMAT in global matrix AMATR
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  nzdom,ndime,r_sol,c_sol
  use def_domain, only       :  r_sym,c_sym,nzsym
  use def_solver, only       :  nzdom_aii,nzdom_aib,nzdom_abi
  use def_master, only       :  solve
  implicit none
  integer(ip), intent(in)    :: itask,pnode,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elmat(pevat,pevat)
  real(rp),    intent(inout) :: Auu(ndime,ndime,nzdom) 
  real(rp),    intent(inout) :: Aup(ndime,nzdom)
  real(rp),    intent(inout) :: Apu(ndime,nzdom)
  real(rp),    intent(inout) :: App(solve(2)%nzmat)
  integer(ip)                :: ndofn,inode,jnode,iposi,jposi,izsym 
  integer(ip)                :: idime,jdime,izsol,jpoin,ipoin,jcolu
  integer(ip)                :: poaii,poaib,poabi,poabb

  ndofn = ndime+1

  select case ( itask )

  case ( -1_ip )
     !
     ! Assemble Schur complement preconditioner (renamed here App)
     !
     if( solve(2) % kfl_schur == 1 ) then
        poaii = 1
        poaib = poaii + nzdom_aii
        poabi = poaib + nzdom_aib 
        poabb = poabi + nzdom_abi
       
        call csrshu(elmat,App(poaii),App(poaib),App(poabi),App(poabb),&
             1_ip,pnode,pnode,lnods,2_ip)
     else
        if( solve(2) % nzmat == nzsym ) then
           call csrase(elmat,App,1_ip,pnode,pnode,lnods,3_ip)
        else
           call csrase(elmat,App,1_ip,pnode,pnode,lnods,2_ip)
        end if
     end if

  case ( 1_ip )

     !solve(1) % kfl_assem = 1 ! Momentum   equation is assembled
     !solve(2) % kfl_assem = 1 ! Continuity equation is assembled

     if( solve(2) % nzmat == nzsym ) then
        !
        ! App requires symmetric assembly
        !
        do inode = 1,pnode
           ipoin = lnods(inode)
           do jnode = 1,pnode
              jpoin = lnods(jnode)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while(jcolu/=jpoin)
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              do idime=1,ndime
                 iposi=(inode-1)*ndofn+idime
                 do jdime=1,ndime
                    jposi                  = (jnode-1)*ndofn+jdime
                    !$OMP ATOMIC
                    Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) + elmat(iposi,jposi)  ! Auu
                 end do
                 jposi            = jposi+1
                 !$OMP ATOMIC
                 Aup(idime,izsol) = Aup(idime,izsol) + elmat(iposi,jposi)                 ! Aup
              end do
              iposi = iposi+1
              do jdime=1,ndime
                 jposi = (jnode-1)*ndofn+jdime
                 !$OMP ATOMIC
                 Apu(jdime,izsol) = Apu(jdime,izsol) + elmat(iposi,jposi)                 ! Apu                 
              end do
              if(jpoin<=ipoin) then 
                 izsym = r_sym(ipoin)
                 jcolu = c_sym(izsym)
                 do while(jcolu/=jpoin)
                    izsym = izsym + 1
                    jcolu = c_sym(izsym)
                 end do
                 jposi      = jnode*ndofn
                 iposi      = inode*ndofn
                 !$OMP ATOMIC
                 App(izsym) = App(izsym) + elmat(iposi,jposi)                             ! App
              end if
           end do
        end do
     else
        !
        ! App requires asymmetric assembly
        !
        do inode = 1,pnode
           ipoin = lnods(inode)
           do jnode = 1,pnode
              jpoin = lnods(jnode)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while(jcolu/=jpoin .and. izsol < r_sol(ipoin+1)-1 )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jcolu == jpoin ) then
                 iposi=(inode-1)*ndofn
                 do idime=1,ndime
                    iposi=iposi+1
                    do jdime=1,ndime
                       jposi                  = (jnode-1)*ndofn+jdime
                       !$OMP ATOMIC
                       Auu(jdime,idime,izsol) = Auu(jdime,idime,izsol) + elmat(iposi,jposi) ! Auu
                    end do
                    jposi            = jposi+1
                    !$OMP ATOMIC
                    Aup(idime,izsol) = Aup(idime,izsol) + elmat(iposi,jposi) ! Aup
                 end do
                 iposi = iposi+1
                 jposi = (jnode-1)*ndofn
                 do jdime=1,ndime
                    jposi = jposi + 1
                    !$OMP ATOMIC
                    Apu(jdime,izsol) = Apu(jdime,izsol) + elmat(iposi,jposi)  ! Apu                 
                 end do
                 iposi      = ndofn*inode
                 jposi      = ndofn*jnode
                 !$OMP ATOMIC                 
                 App(izsol) = App(izsol) + elmat(iposi,jposi)                 ! App
              end if
           end do
        end do
     end if

  end select

end subroutine nsi_assmat
