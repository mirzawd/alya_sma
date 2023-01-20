!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_redsus
  !-----------------------------------------------------------------------
  !****f* Levels/lev_redsus
  ! NAME
  !    lev_redieq
  ! DESCRIPTION
  !    Compute the level set function redistanciation with Sussman equation the non cut nodes
  ! USES
  !
  ! USED BY
  !    lev_redieq
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_elmtyp
  use def_kermod
  use def_levels
  use def_domain
  use def_solver
  use mod_gradie
  use mod_communications
  implicit none
  integer(ip)             :: ielem,idime
  integer(ip)             :: inode,ipoin,it
  integer(ip)             :: pelty,pnode
  integer(ip), save       :: ipass=0
  real(rp)                :: elcod(ndime,mnode),elvel(ndime,mnode)
  real(rp)                :: dtcri,rdinv,gpdet,h,c,v,time,time1
  real(rp)                :: gpcar(ndime,mnode),xjacm(9),xjaci(9)
  real(rp)                :: norm, signl,dtmin
  real(rp),    target     :: resid

  !----------------------------------------------------------------------
  !
  ! Solve equation
  !
  !----------------------------------------------------------------------
  call cputim(time)

  do it = 1,nbitr_lev
     !
     ! grad(phi)
     !
     call grasca(fleve,norml_lev)
     !
     ! sign(phi) * grad(phi) / |grad(phi)|
     !
     if(INOTMASTER) then

        do ipoin = 1,npoin

           norm = 0.0_rp
           do idime = 1,ndime
              norm = norm + norml_lev(idime,ipoin)*norml_lev(idime,ipoin)
           enddo
           norm = sqrt(norm)

           if( flev0_lev(ipoin) == 0.0_rp ) then
              signl = 1.0_rp
           else
              signl = flev0_lev(ipoin)/abs(flev0_lev(ipoin))
           endif

           if( norm /= 0.0_rp ) then
              norm = signl / norm
              do idime = 1,ndime
                 norml_lev(idime,ipoin) = norml_lev(idime,ipoin)*norm
              end do
           end if
        end do
        !
        ! Get critical time step
        !
        c        = 0.0_rp
        v        = 0.0_rp
        rdinv    = 1.0_rp/real(ndime)
        dtmin    = 1e+06
        dtcri    = 1e+06

        do ielem = 1,nelem
           pelty = ltype(ielem)

           if( lelch(ielem) /= ELHOL ) then

              pnode = nnode(pelty)

              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                    elvel(idime,inode) = norml_lev(idime,ipoin)
                 end do
              end do
              do inode = 1,pnode
                 v = 0.0_rp
                 do idime = 1,ndime
                    v = v + elvel(idime,inode) * elvel(idime,inode)
                 end do
                 c = max( sqrt(v) , c )
              end do

              call elmder(pnode,ndime,elmar(pelty)%dercg,elcod,gpcar,gpdet,xjacm,xjaci)
              h = (elmar(pelty)%weicg*gpdet)**rdinv    ! h
              if ( c > 1.e-10_rp ) then
                 dtcri = h/c                           ! dt=h/c
              end if
              dtmin = min( dtmin , dtcri )

           end if

        end do
     end if

     if( IPARALL ) then
        call PAR_MIN(dtmin)
     endif
     dtred_lev = 1.0_rp/(dtmin*1.0_rp)

     !
     ! Assemble and solve system
     !
     call lev_inisol(2_ip)
     call inisol()

     if( INOTMASTER ) call lev_elmope(3_ip)

     call solver(rhsid,unkno,amatr,pmatr)
     !
     ! Compute and output residual
     !
     call residu(2_ip,one,one,unkno,fleve,one,one,one,1.0_rp,resid)

     if( INOTSLAVE ) then
        if( ipass == 0 .and. kfl_rstar /= 2 ) then
           ipass = 1
           write(lun_cored_lev,100)
        end if
        call cputim(time1)
        write(lun_cored_lev,101) ittim,it,resid,time1-time,dtmin
        flush(lun_cored_lev)
     end if
     !
     ! Update solution: FLEVE <= UNKNO
     !
     call lev_updunk(3_ip)

  end do

  if( INOTMASTER ) then
     if(kfl_locre_lev==1) then
        do ipoin=1,npoin
           if(abs(flev0_lev(ipoin))>10.0_rp*thicl) then
              fleve(ipoin,1) = flev0_lev(ipoin)
           end if
        end do
     end if
  end if


100 format('# --| Convergence of redistanciation '       ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time step  2. Inner Iteration   '  ,/,&
       & '# --| 3. Norm       4. Redistanciation duration  5. Time step')
101 format(4x,i9,2x,i9,6(2x,e12.6))

end subroutine lev_redsus
