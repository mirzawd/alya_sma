!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_smooth()
  !-----------------------------------------------------------------------
  !****f* domain/ale_smooth
  ! NAME
  !    domain
  ! DESCRIPTION
  !    Laplacian ale_smoothing using Gauss-Seidel
  ! USED BY
  !    Turnon 
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_elmtyp
  use def_parame
  use def_elmtyp
  use def_alefor
  use def_domain
  use def_kermod, only : kfl_element_to_csr
  use mod_solver, only : solver_solve
  use mod_solver, only : solver_initialize_matrix_and_rhs
  use mod_matrix, only : matrix_assemble_element_RHS
  use mod_matrix, only : matrix_assemble_element_matrix_to_CSR
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  implicit none
  integer(ip)           :: iters,ipoin,inode,idime
  integer(ip)           :: ielem,pelty,pnode,idofn
  integer(ip)           :: kfl_recov
!  integer(ip)           :: izdom,jpoin
  real(rp)              :: xjacm(9),w1,w
!  real(rp)              :: ri(3),ro(3),gpdet
  real(rp)              :: elcod(ndime,mnode)
!  real(rp)              :: gpcar(ndime,mnode,mgaus)
  real(rp)              :: hleng(3),asele,resid,denom
  integer(ip), pointer  :: kfl_fixno_sav(:,:)
  integer(ip), pointer  :: kfl_fixrs_sav(:)
  real(rp),    pointer  :: bvess_sav(:,:)
  real(rp),    pointer  :: unkno_sav(:)
  real(rp),    pointer  :: invdiag(:)

  integer(ip)           :: kfl_elfix(ndime,mnode)
  real(rp)              :: elmat(ndime*mnode,ndime*mnode)
  real(rp)              :: elrhs(ndime,mnode)
  real(rp)              :: eless(ndime,mnode)

  w  = resmo_ale
  w1 = (1.0_rp-w)

  nullify(kfl_fixno_sav)
  nullify(kfl_fixrs_sav)
  nullify(bvess_sav)
  nullify(unkno_sav)
  nullify(invdiag)

  call memory_alloca(mem_modul(1:2,modul),'UNKNO_SAV'    ,'ale_smooth',unkno_sav,max(1_ip,ndime*npoin))
  call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_SAV','ale_smooth',kfl_fixno_sav,ndime,max(1_ip,npoin))
  call memory_alloca(mem_modul(1:2,modul),'KFL_FIXRS_SAV','ale_smooth',kfl_fixrs_sav,max(1_ip,npoin))
  call memory_alloca(mem_modul(1:2,modul),'BVESS_SAV'    ,'ale_smooth',bvess_sav,ndime,max(1_ip,npoin))

  kfl_recov = solve_sol(1) % kfl_recov
  solve_sol(1) % kfl_recov = 1

  if( INOTEMPTY ) then

     !-------------------------------------------------------------------
     !
     ! Assemble matrix
     !
     !-------------------------------------------------------------------

     kfl_fixno_sav    = kfl_fixno_ale
     kfl_fixrs_sav    = kfl_fixrs_ale
     bvess_sav        = bvess_ale(:,:,1)
     bvess_ale(:,:,1) = coord_ale(:,:,1)
     
     do ipoin = 1,npoin
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime+idime
           unkno(idofn)     = coord_ale(idime,ipoin,1) 
           unkno_sav(idofn) = unkno(idofn)
        end do
     end do

     do ipoin = 1,npoin
        kfl_fixrs_ale(ipoin) = 0
        if( lpoty(ipoin) /= 0 ) then
           kfl_fixno_ale(1:ndime,ipoin) = 1
        else
           kfl_fixno_ale(1:ndime,ipoin) = 0
        end if
     end do
     !
     ! Prescribe hole nodes
     !
     do ipoin = 1,npoin
        if( lnoch(ipoin) == NOHOL ) then
           kfl_fixno_ale(1:ndime,ipoin) = 1
        end if
     end do
  end if
  !
  ! Assemble matrix
  !
  if( kfl_smoot_ale == 1 ) then

     do iters = 1,nsmoo_ale

        call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)

        do ielem = 1,nelem

           pelty = ltype(ielem)

           if( pelty > 0 ) then

              pnode = nnode(pelty)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = unkno((ipoin-1)*ndime+idime)
                 end do
                 kfl_elfix(1:ndime,inode) = kfl_fixno_ale(1:ndime,ipoin)
                 eless(1:ndime,inode)     = bvess_ale(1:ndime,ipoin,1)
              end do
              call elmlen(&
                   ndime,pnode,elmar(pelty) % dercg,xjacm,elcod,&
                   hnatu(pelty),hleng)
              asele = hleng(1) / hleng(ndime)
              call ale_elmpno(&
                   pnode,pelty,lelch(ielem),kfl_elfix,solve_sol(1) % kfl_iffix,&
                   elcod,eless,asele,elmat,elrhs)
              call matrix_assemble_element_RHS(&
                   solve_sol(1) % ndofn,solve_sol(1) % ndofn,pnode,&
                   lnods(:,ielem),elrhs,rhsid)
              call matrix_assemble_element_matrix_to_CSR(&
                   kfl_element_to_csr,solve_sol(1) % ndofn,pnode,pnode*ndime,&
                   ielem,lnods(:,ielem),elmat,solve_sol(1) % ia,&
                   solve_sol(1) % ja,amatr,lezdo)  
           end if
        end do

        !-------------------------------------------------------------------
        !
        ! Implicit: Solve non-linear system
        !
        !-------------------------------------------------------------------
        
        call solver_solve(solve_sol,amatr,rhsid,unkno)
        

        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn            = (ipoin-1)*ndime+idime
              unkno(idofn)     = w1 * unkno_sav(idofn) + w * unkno(idofn)
              unkno_sav(idofn) = unkno_sav(idofn)-unkno(idofn)
           end do
        end do

        call norm2x(ndime,unkno_sav,resid)
        call norm2x(ndime,unkno,denom)
        !if( INOTSLAVE ) write(90,*) iters,resid/max(denom,zeror)

        do ipoin = 1,npoin*ndime
           unkno_sav(ipoin) = unkno(ipoin)
        end do

     end do

     solve_sol(1) % kfl_recov = kfl_recov 
     !
     ! Unknow should be the displacement to be loaded in DISPM
     !
     do ipoin = 1,npoin
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime+idime
           unkno(idofn) = unkno(idofn)-coord_ale(idime,ipoin,1)
        end do
     end do

  end if

  if( 2 == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Gauss-Seidel
     !
     !-------------------------------------------------------------------

     call runend('ALE_SMOOTH: RECODE')
!!$     if( INOTMASTER ) then
!!$
!!$        call memory_alloca(mem_modul(1:2,modul),'UNKNO_TMP','ale_smooth',unkno_tmp,ndime,npoin)
!!$        call memory_alloca(mem_modul(1:2,modul),'INVDIAG'  ,'ale_smooth',invdiag,npoin)
!!$
!!$        do ipoin = 1,npoin
!!$           izdom = r_dom(ipoin)
!!$           jpoin = c_dom(izdom)
!!$           do while( jpoin /= ipoin )
!!$              izdom = izdom + 1
!!$              jpoin = c_dom(izdom)
!!$           end do
!!$           invdiag(ipoin) = amatr_tmp(izdom)
!!$        end do
!!$        call pararr('SLX',NPOIN_TYPE,npoin,invdiag)
!!$        call pararr('SLX',NPOIN_TYPE,ndime*npoin,rhsid_tmp)
!!$        do ipoin = 1,npoin
!!$           invdiag(ipoin) = 1.0_rp / invdiag(ipoin)
!!$        end do
!!$        !
!!$        ! Scale matrix: Put zero on invdiagnal
!!$        !
!!$        do ipoin = 1,npoin
!!$           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$              jpoin = c_dom(izdom)
!!$              if( ipoin == jpoin ) then
!!$                 amatr_tmp(izdom) = invdiag(ipoin) * amatr_tmp(izdom) 
!!$              else
!!$                 amatr_tmp(izdom) = invdiag(ipoin) * amatr_tmp(izdom) 
!!$              end if
!!$           end do
!!$        end do
!!$        !
!!$        ! Loop over iterations
!!$        !
!!$        do iters = 1,nsmoo_ale
!!$           !
!!$           ! Update interior nodes
!!$           !
!!$           do ipoin = 1,npoi1
!!$              if( kfl_fixno_ale(1,ipoin) == 0 ) then
!!$
!!$                 w = resmo_ale
!!$                 do idime = 1,ndime
!!$                    ri(idime) = coord_ale(idime,ipoin,1) + invdiag(ipoin) * rhsid_tmp(idime,ipoin)
!!$                    ro(idime) = coord_ale(idime,ipoin,1)
!!$                 end do
!!$
!!$                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$                    jpoin = c_dom(izdom)
!!$
!!$                    do idime = 1,ndime
!!$                       ri(idime) = ri(idime) - amatr_tmp(izdom) * coord_ale(idime,jpoin,1) 
!!$                    end do
!!$                 end do
!!$
!!$                 do idime = 1,ndime
!!$                    coord_ale(idime,ipoin,1) = (1.0_rp-w) * ro(idime) + w * ri(idime)
!!$                 end do
!!$
!!$                 if( kfl_smoot_ale == 3 ) then
!!$                    !
!!$                    ! Check Jacobian
!!$                    !
!!$                    do idime = 1,ndime
!!$                       unkno_tmp(idime,ipoin) = coord_ale(idime,ipoin,1)
!!$                    end do
!!$                    izdom = pelpo(ipoin)-1
!!$                    do while( izdom < pelpo(ipoin+1)-1 )
!!$                       izdom = izdom + 1
!!$                       ielem = lelpo(izdom)
!!$                       pelty = ltype(ielem)
!!$                       pnode = nnode(pelty)
!!$                       do inode = 1,pnode
!!$                          jpoin = lnods(inode,ielem)
!!$                          do idime = 1,ndime
!!$                             elcod(idime,inode) = coord_ale(idime,jpoin,1)
!!$                          end do
!!$                       end do
!!$                       call determ(&
!!$                            ndime,pnode,elcod,elmar(pelty) % dercg,xjacm,gpcar,gpdet)
!!$                       if( gpdet <= 0.0_rp ) then
!!$                          w = w / 2.0_rp
!!$                          do idime = 1,ndime
!!$                             coord_ale(idime,ipoin,1) = (1.0_rp-w) * ro(idime) + w * ri(idime)
!!$                          end do
!!$                          izdom = pelpo(ipoin) - 1
!!$                       end if
!!$                    end do
!!$                 end if
!!$              end if
!!$           end do
!!$           !
!!$           ! Update boundary nodes
!!$           !
!!$           do ipoin = npoi1+1,npoin
!!$              if( kfl_fixno_ale(1,ipoin) == 0 ) then
!!$                 do idime = 1,ndime
!!$                    unkno_tmp(idime,ipoin) = 0.0_rp
!!$                 end do
!!$                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$                    jpoin = c_dom(izdom)
!!$                    do idime = 1,ndime
!!$                       unkno_tmp(idime,ipoin) = unkno_tmp(idime,ipoin) - amatr_tmp(izdom) * coord_ale(idime,jpoin,1) 
!!$                    end do
!!$                 end do
!!$              end if
!!$           end do
!!$
!!$           call pararr('SLX',NPOIN_TYPE,ndime*npoin,unkno_tmp)
!!$
!!$           do ipoin = npoi1+1,npoin
!!$              if( kfl_fixno_ale(1,ipoin) == 0 ) then
!!$                 do idime = 1,ndime
!!$                    unkno_tmp(idime,ipoin) = unkno_tmp(idime,ipoin) &
!!$                         + coord_ale(idime,ipoin,1) + invdiag(ipoin) * rhsid_tmp(idime,ipoin)
!!$                    coord_ale(idime,ipoin,1) = (1.0_rp-w) * coord_ale(idime,ipoin,1) + w * unkno_tmp(idime,ipoin)
!!$                 end do
!!$              end if
!!$           end do
!!$
!!$        end do
!!$
!!$        call memory_deallo(mem_modul(1:2,modul),'UNKNO_TMP','ale_smooth',unkno_tmp)
!!$        call memory_deallo(mem_modul(1:2,modul),'INVDIAG'  ,'ale_smooth',invdiag)
!!$  end if

  end if
  !
  ! Deallocate memory
  !
  if( INOTEMPTY ) then
     kfl_fixno_ale = kfl_fixno_sav
     kfl_fixrs_ale = kfl_fixrs_sav
     bvess_ale(:,:,1) = bvess_sav
  end if
  call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_SAV','ale_smooth',kfl_fixno_sav)
  call memory_deallo(mem_modul(1:2,modul),'BVESS_SAV'    ,'ale_smooth',bvess_sav)
  call memory_deallo(mem_modul(1:2,modul),'UNKNO_SAV'    ,'ale_smooth',unkno_sav)

end subroutine ale_smooth

subroutine ale_smooth_fixity(ndofn,iffix,xv,aa,bb,xx)

  use def_domain
  
  implicit none
  integer(ip), intent(in)    :: ndofn
  integer(ip), intent(in)    :: iffix(ndofn,npoin)
  real(rp),    intent(in)    :: xv(ndofn,npoin)
  real(rp),    intent(inout) :: aa(ndofn,ndofn,nzdom)
  real(rp),    intent(out)   :: bb(ndofn,npoin)
  real(rp),    intent(out)   :: xx(ndofn,npoin)
  integer(ip)                :: ii,izdom,jj,idofn
  real(rp)                   :: adiag(ndofn)

  do ii = 1,npoin 
     if( iffix(1,ii) > 0 ) then
        do izdom = r_dom(ii),r_dom(ii+1)-1
           jj = c_dom(izdom)
           if( ii == jj ) then              
              do idofn = 1,ndofn
                 if( aa(idofn,idofn,izdom) /= 0.0_rp ) then
                    adiag(idofn) = aa(idofn,idofn,izdom)
                 else
                    adiag(idofn) = 1.0_rp
                 end if
              end do
              aa(:,:,izdom) = 0.0_rp
              do idofn = 1,ndofn
                 aa(idofn,idofn,izdom) = adiag(idofn)
              end do
           else
              aa(:,:,izdom) = 0.0_rp
           end if
        end do
        bb(:,ii) = adiag(:) * xv(:,ii)
        xx(:,ii) = xv(:,ii)
     end if
  end do

end subroutine ale_smooth_fixity
