!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_walgen(itask)
  !-----------------------------------------------------------------------
  !****f* domain/ker_walgen
  ! NAME 
  !    ker_walgen
  ! DESCRIPTION
  !    Compute the generalized distance to the wall via a 
  !     IF (kfl_walld == 1)
  !       Poisson equation:
  !       1. Solve Lapl(f)=-1, with f=0 on wall
  !       2. d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !       See the following references:
  !       P.G. Tucker, Differential equation-based wall distance computation for
  !         DES and RANS, J. Comp. Phys. 190 (2003) 229-248.
  !       P.G. Tucker, Int. J. Numer. Fluids 33 (2000) 869.
  !       P.G. Tucker, Appl. Math. Model. 22 (1998) 293.
  !
  !     IF (kfl_walld == 2)
  !        Search the mimimum distance 
  ! USES
  ! USED BY
  !    Domain
  !*** 
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_gradie
  use mod_memory
  use mod_solver
  use mod_messages,      only : messages_live
  use mod_moduls_conf,   only : moduls_set_current_solver
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem,igaus,idime,inode,ipoin,pnode,pgaus,pelty,imeth
  integer(ip)             :: iiter,niter,izmat,jzmat,jpoin
  real(rp)                :: elmal(mnode,mnode)
  real(rp)                :: elrhs(mnode)
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: gpcar(ndime,mnode,mgaus)
  real(rp)                :: gpvol(mgaus),gpdet,fact1
  real(rp)                :: xjaci(9)
  real(rp)                :: xjacm(9)
  real(rp)                :: fact2,eps,xnorm,time1,time2
  real(rp),    pointer    :: auxwd(:)
  integer(ip), pointer    :: auxwi(:)
  integer(ip)             :: npoin_walld,i
  real(rp)                :: distance_min,distance2
  real(rp)                :: pcoor(ndime),wcoor(ndime)

  nullify(auxwd)

  if( kfl_walld == 1 ) then

     call moduls_set_current_solver(ID_KERMOD,'WALL_DISTANCE')
     
     if( solve_sol(1) % kfl_algso /= -999 ) then

        !-------------------------------------------------------------------
        !
        ! Extend distance on the wall to the rest of the domian smoothly to be used instead of delta_tur that is not constant
        ! Solve system: Lapl(g)=0, with g = ywalp on nodes with kfl_fixno_tur == 3 or 4
        ! For the moment I use the same solver used to calculate wall distance (next step) and I only leave the option corresponding to imeth=0 
        !     
        !-------------------------------------------------------------------

        if( kfl_delta == 1 ) then
           !
           ! Initialize solver
           !
           call messages_live('KERMOD: EXTEND DISTANCE FROM THE WALL')
           call cputim(time1)

           if( INOTMASTER ) then

              call inisol()
              if( itask == 1 ) then
                 do ipoin = 1,npoin
                    unkno(ipoin) = 0.0_rp
                 end do
              else
                 do ipoin = 1,npoin
                    unkno(ipoin) = uwal2_ker(ipoin)
                 end do
              end if

              do ielem = 1,nelem

                 if( lelch(ielem) /= ELHOL ) then
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
                          !                            if (lninv_loc(ipoin) == 4160 .and. idime == 2) elcod(idime,inode) = elcod(idime,inode) - 0.000001_rp
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
                    ! Compute element matrix ELMAL and assemble LAPLA
                    !
                    call elmlap(&
                         two,pnode,pgaus,lnods(1,ielem),lelch(ielem),gpcar,&
                         elmar(pelty)%shape,gpvol,elmal,elrhs)
                    call assmat(&
                         1_ip,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
                         ielem,lnods(1,ielem),elmal,amatr)
                    call assrhs(&
                         1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
                 end if
              end do
           end if
           call cputim(time2) 
           cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time2 - time1
           !
           ! Solve system: Lapl(g)=0, with g = ywalp on nodes with gisca = 1
           !           
           call solver(rhsid,unkno,amatr,pmatr) 

           if( INOTMASTER ) then 

              call memory_alloca(mem_modul(1:2,modul),'AUXWD','ker_walgen',auxwd,npoin)

              do ipoin = 1,npoin
                 auxwd(ipoin) = unkno(ipoin)
                 uwal2_ker(ipoin) = unkno(ipoin)     ! Actually now that I store uwal2_ker I could eliminate auxwd 
              end do
           end if

        end if !kfl_delta

        !-------------------------------------------------------------------
        !
        ! Initialize solver
        !
        !-------------------------------------------------------------------

        call messages_live('KERMOD: COMPUTE DISTANCE TO THE WALL')
        call cputim(time1)

        if( INOTMASTER ) then

           call inisol()
           if( itask == 1 ) then
              do ipoin = 1,npoin
                 unkno(ipoin) = 0.0_rp
              end do
           else
              do ipoin = 1,npoin
                 unkno(ipoin) = uwall_ker(ipoin)
              end do
           end if
           
           do ielem = 1,nelem

              if( lelch(ielem) /= ELHOL ) then
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
!                         if (lninv_loc(ipoin) == 4160 .and. idime == 2) elcod(idime,inode) = elcod(idime,inode) - 0.000001_rp
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
                 ! Compute element matrix ELMAL and assemble LAPLA
                 ! 
                 call elmlap(&
                      one,pnode,pgaus,lnods(1,ielem),lelch(ielem),gpcar,&
                      elmar(pelty)%shape,gpvol,elmal,elrhs)                 
                 call assmat(&
                      1_ip,pnode,pnode,npoin,solve_sol(1)%kfl_algso,&
                      ielem,lnods(1,ielem),elmal,amatr)
                 call assrhs(&
                      1_ip,pnode,lnods(1,ielem),elrhs,rhsid)
              end if
           end do
        end if

        call cputim(time2) 
        cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time2 - time1

        imeth = 0
        eps   = 1.0_rp

        if( imeth == 0 ) then
           niter = 1
        else
           niter = 100
           jpoin = npoin
           jzmat = solve_sol(1)%nzmat
           do ipoin = 1,npoin
              jpoin = jpoin + 1
              rhsid(jpoin) = rhsid(ipoin)
           end do
           do izmat = 1,solve_sol(1)%nzmat
              jzmat = jzmat + 1
              amatr(jzmat) = amatr(izmat) 
           end do
           do ipoin = 1,npoin
              if( kfl_fixno_walld_ker(1,ipoin) == 1 ) then
                 continue
              else
                 call csrdia(ipoin,solve_sol(1)%kfl_symme,izmat)
                 amatr(izmat) = amatr(izmat) + eps * vmasc(ipoin)
              end if
           end do
        end if

        iiter = 0
        do while( iiter < niter )
           iiter = iiter + 1
           !
           ! Penalize equation
           !
           if( imeth == 1 ) then
              do ipoin = 1,npoin
                 if( kfl_fixno_walld_ker(1,ipoin) == 1 ) then
                    continue
                 else
                    jpoin = npoin + ipoin
                    rhsid(ipoin) = rhsid(jpoin) + eps * vmasc(ipoin) * unkno(ipoin) 
                 end if
              end do
           end if
           !
           ! Solve system: Lapl(f)=-1, with f=0 on wall
           !
           call solver(rhsid,unkno,amatr,pmatr)           
           !
           ! Check convergence
           !
           if( imeth == 1 ) then
              call resnor(1_ip,solve_sol(1)%kfl_symme,solve_sol(1)%kfl_full_rows,rhsid(npoin+1),&
                   unkno,amatr(solve_sol(1)%nzmat+1),xnorm)
              !write(88,*)
           end if

        end do

        !-------------------------------------------------------------------
        !
        ! Compute wall distance: d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
        !
        !-------------------------------------------------------------------

        if( INOTEMPTY ) then
           call memgen(zero,ndime,npoin)
           call gradie(unkno,gevec)
           do ipoin = 1,npoin
              fact1 = 0.0_rp 
              do idime = 1,ndime
                 fact1 = fact1 + gevec(idime,ipoin) * gevec(idime,ipoin)
              end do
              uwall_ker(ipoin) = unkno(ipoin)
              fact2 = fact1 + 2.0_rp*max(unkno(ipoin),0.0_rp)
              if( fact2 < 0.0_rp ) then
                 call runend('WRONG DISTANCE TO THE WALL')
              else
                 if (kfl_delta==1) then
                    walld(ipoin) = sqrt(fact2) - sqrt(fact1) + auxwd(ipoin)
                    !                 walld(ipoin) = auxwd(ipoin)   ! temp test  to see only auxwd  in walld !!!!!!!!!!!!!!!!!
                 else
                    walld(ipoin) = sqrt(fact2) - sqrt(fact1) + delta_dom
                 end if
                 !walld(ipoin) = unkno(ipoin)
              end if
           end do
           !
           ! Deallocate memory
           !
           call memgen(two,ndime,npoin)

        end if

        if( ( kfl_delta == 1 ) .and.  INOTMASTER ) then         ! deallocates auxwd
           call memory_deallo(mem_modul(1:2,modul),'AUXWD','ker_walgen',auxwd)
        end if

     end if

  elseif (kfl_walld == 2) then
     !-------------------------------------------------------------------
     !
     ! Compute wall distance based on the exhaustive search
     !
     !-------------------------------------------------------------------
     nullify(auxwi)
     call messages_live('KERMOD: COMPUTE DISTANCE TO THE WALL EXHAUSTIVE SEARCH')
     open(10,file='wallo.dat')
     open(11,file='wallcoor.dat')

     npoin_walld = 0
     do ipoin = 1,npoin
        if( kfl_fixno_walld_ker(1,ipoin) == 1 ) then
           npoin_walld = npoin_walld + 1
        endif
     enddo
     call memory_alloca(mem_modul(1:2,modul),'AUXWI','ker_walgen',auxwi,npoin_walld)

     i = 0
     do ipoin = 1,npoin
        if( kfl_fixno_walld_ker(1,ipoin) == 1 ) then
           i = i + 1
           auxwi(i) = ipoin
        endif
     enddo

     do ipoin = 1,npoin

        distance_min = 1.0e10_rp
        pcoor(1:ndime) = coord(1:ndime,ipoin)
        do jpoin = 1,npoin_walld

           wcoor(1:ndime) = coord(1:ndime,auxwi(jpoin))
           distance2 = 0.0_rp
           do idime = 1,ndime
              distance2 = distance2 + (pcoor(idime)-wcoor(idime))**2
           enddo
           if ( sqrt(distance2) < distance_min) then
              !           walli(ipoin) = auxwi(jpoin)
              wallo(ipoin) = jpoin
              !           walli(ipoin) = jpoin
              wallcoor(1:ndime,ipoin) = wcoor(1:ndime)
              distance_min = sqrt(distance2)
           endif

        enddo ! jpoin

        !       write (10,*) ipoin,walli(ipoin)
        write (10,*) ipoin,wallo(ipoin)
        write (11,'(i8,3(1x,e22.14))')ipoin,wallcoor(:,ipoin) 


     enddo !ipoin

     ! deallocates auxwi
     call memory_deallo(mem_modul(1:2,modul),'AUXWI','ker_walgen',auxwi)
     close(10)
     close(11)

  elseif (kfl_walld == 3) then

     call messages_live('KERMOD: COMPUTE DISTANCE TO THE WALL READ FILES')
    !
    ! Read from the field file
    !
    if( INOTMASTER ) then
      do ipoin=1,npoin
!         walli(ipoin) = int(xfiel(-kfl_walld_field(1))%a(1,ipoin))
        wallo(ipoin) = int(xfiel(-kfl_walld_field(1))%a(1,ipoin,1))
        wallcoor(1:ndime,ipoin) = xfiel(-kfl_walld_field(2))%a(1:ndime,ipoin,1)
        pcoor(1:ndime) = coord(1:ndime,ipoin)
        wcoor(1:ndime) = wallcoor(1:ndime,ipoin)
        distance2 = 0.0_rp
        do idime = 1,ndime
          distance2 = distance2 + (pcoor(idime)-wcoor(idime))**2
        enddo
        walld(ipoin) = sqrt(distance2)
      end do
    endif    
    
    
  end if !kfl_walld


end subroutine ker_walgen
