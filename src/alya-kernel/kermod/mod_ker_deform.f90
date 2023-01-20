!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ker_deform

  use def_kintyp,         only : ip,rp
  use def_kintyp_solvers, only : soltyp
  use def_domain,         only : elmar,ndime,mnode,mgaus,nelem
  use def_domain,         only : lnods,hnatu,ltype,lnnod,lelch,ngaus,npoin,lnoch
  use def_domain,         only : lezdo,walld
  use def_master,         only : INOTMASTER,solve_sol,lzone,NPOIN_TYPE,kfl_paral
  use def_master,         only : lninv_loc,ittim,momod,modul,IMASTER,cpu_modul,npoi3
  use def_master,         only : mem_modul,CPU_ASSEMBLY
  use def_elmtyp,         only : ELEXT,ELHOL,QUA04,HEX08,PEN06,NODE_CONTACT_FLUID
  use mod_solver,         only : solver_solve
  use def_solver,         only : SOL_MATRIX_HAS_CHANGED
  use mod_matrix,         only : matrix_assemble_element_RHS
  use mod_matrix,         only : matrix_assemble_element_matrix_to_CSR  
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_MIN
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use def_kermod,         only : kfl_adj_prob,sens_mesh
  use def_kermod,         only : kfl_element_to_csr
  use mod_ker_timeline,   only : ker_timeline

  private

  public :: deform_deform

contains

  subroutine deform_deform(&
       ndefo,kfl_dmeth,defor_param,imodu,kfl_fixno,bvess,&
       coord,amatr,unkno,rhsid,solve_in)
    !-----------------------------------------------------------------------
    !****f* domain/ale_deform
    ! NAME
    !    domain
    ! DESCRIPTION
    !    This routines ale_deformes the mesh. Idea: Put a node at the middle
    !    of its neighbors. A Gauss Seidel in 1D would give:
    !
    !    X^{k+1}(i) = 1/2 [  X^k(i-1) + X^k(i+1) ]
    !    X^{k+1}(i) = X^k(i) + 1/2 [ X^k(i-1) -2 X^k(i) + X^k(i+1) ]
    !    X^{k+1}(i) = X^k(i) + 1/2 h^2 [ X^k(i-1) -2 X^k(i) + X^k(i+1) ]/h^2
    !    X^{k+1}(i) = X^k(i) + h^2 [ 0 - Lapl(X) ]
    !
    !    Let X = x0+d. x0 is initial solution. d is displacement.
    !    Solve div[ a * grad(X) ] = 0 with X=x0 on boundary
    ! 
    !    In 1D:
    !
    !    i-1   i   i+1
    !     o----o----o
    !        h    h
    !    +-
    !    |  a* grad(X).grad(v) dx =   a*h* [X(i)-X(i-1)]/h (1/h) 
    !   -+                          + a*h* [X(i+1)-X(i)]/h (-1/h) 
    !                             = -alpha*h* Lapla(X)
    !   => a = h
    !
    ! USED BY
    !    Turnon 
    !***
    !-----------------------------------------------------------------------
    implicit none
    integer(ip),                       intent(in)    :: ndefo            !< # deformaiton steps
    integer(ip),                       intent(in)    :: kfl_dmeth        !< Deformation method
    real(rp),                          intent(in)    :: defor_param      !< Deformation parameter
    integer(ip),                       intent(in)    :: imodu            !< Current module
    integer(ip),  pointer,             intent(in)    :: kfl_fixno(:,:)   !< Fixity
    real(rp),     pointer,             intent(inout) :: bvess(:,:)       !< Boundary value
    real(rp),     pointer,             intent(inout) :: coord(:,:)       !< Coordinates
    real(rp),     pointer, contiguous, intent(inout) :: amatr(:)         !< Matrix
    real(rp),     pointer,             intent(inout) :: unkno(:)         !< Displacement
    real(rp),     pointer,             intent(inout) :: rhsid(:)         !< RHS
    type(soltyp), pointer,             intent(inout) :: solve_in(:)      !< Solver pointer
    integer(ip)                                      :: ielem,igaus,idime,inode
    integer(ip)                                      :: ipoin,pnode,pgaus
    integer(ip)                                      :: itotn,pevat,pelty,iters
    integer(ip)                                      :: kfl_matdi
    integer(ip)                                      :: kfl_elfix(ndime,mnode)
    real(rp)                                         :: elmat(mnode*mnode*ndime*ndime)
    real(rp)                                         :: elrhs(mnode*ndime)
    real(rp)                                         :: elcod(ndime,mnode)
    real(rp)                                         :: elbve(ndime,mnode)
    real(rp)                                         :: gpwal(mgaus)
    real(rp)                                         :: gpcar(ndime,mnode,mgaus)
    real(rp)                                         :: gpvol(mgaus),hleng(3)
    real(rp)                                         :: xjaci(9),xjacm(9),gpdet
    real(rp)                                         :: asmin,asmax,asele
    real(rp)                                         :: vomin,vomax,xfact,voele
    real(rp)                                         :: time1,time2
    real(rp),     pointer                            :: dispm_inc(:,:)

    nullify(dispm_inc)
    solve_sol => solve_in
    kfl_matdi =  solve_sol(1) % kfl_iffix 
    xfact     =  1.0_rp / real(ndefo,rp)
    !
    ! Incremental deformation
    !
    if( ndefo > 1 .and. INOTMASTER ) then 
       call memory_alloca(mem_modul(1:2,modul),'DISPM_INC','deform_deform',dispm_inc,ndime,npoin)
    end if
    
    do iters = 1,ndefo

       !-------------------------------------------------------------------
       !
       ! VOMIN, VOMAX
       !
       !-------------------------------------------------------------------

       call ker_timeline('INI_ASSEMBLY')

       if( INOTMASTER ) then

          vomin = 1e09_rp
          vomax =  0.0_rp
          asmin = 1e09_rp
          asmax =  0.0_rp

          do ielem = 1,nelem
             pelty = ltype(ielem)

             if( pelty > 0 ) then
                pnode = lnnod(ielem)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do
                if( ndefo > 1 ) then
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1:ndime,inode) = elcod(1:ndime,inode) + dispm_inc(1:ndime,ipoin)
                   end do
                end if
                call jacobi(&
                     ndime,pnode,elcod,elmar(pelty) % dercg,&
                     xjacm,xjaci,gpcar,gpdet)
                call elmlen(&
                     ndime,pnode,elmar(pelty) % dercg,xjacm,elcod,&
                     hnatu(pelty),hleng)
                voele = elmar(pelty) % weicg * gpdet
                asele = hleng(1)/hleng(ndime)
                vomin = min(vomin,voele)
                vomax = max(vomax,voele)
                asmin = min(asmin,asele)
                asmax = max(asmax,asele)
             end if
          end do

       end if

       call PAR_MIN(vomin,'IN MY CODE')
       call PAR_MAX(vomax,'IN MY CODE')
       call PAR_MIN(asmin,'IN MY CODE')
       call PAR_MAX(asmax,'IN MY CODE')

       if( vomin < 0.0_rp ) call runend('DEFORM_DEFORM: NEGATIVE JACOBIAN FOUND')

       !-------------------------------------------------------------------
       !
       ! Assemble system
       !
       !-------------------------------------------------------------------

       call cputim(time1)
 
       if( INOTMASTER ) then

          call inisol()

          do ielem = 1,nelem
             pelty = ltype(ielem)

             if( pelty > 0 ) then
                !
                ! Element properties and dimensions
                !
                pnode = lnnod(ielem)
                pgaus = ngaus(pelty)
                pevat = pnode * ndime
                !
                ! Gather operations: ELCOD
                !
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elcod(1:ndime,inode)     = coord(1:ndime,ipoin)
                   elbve(1:ndime,inode)     = bvess(1:ndime,ipoin) * xfact
                   if (kfl_adj_prob == 1) then
                     elbve(1:ndime,inode) = sens_mesh(1:ndime,ipoin)
                   endif
                   kfl_elfix(1:ndime,inode) = kfl_fixno(1:ndime,ipoin)
                end do
                if( ndefo > 1 ) then
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1:ndime,inode) = elcod(1:ndime,inode) + dispm_inc(1:ndime,ipoin)
                   end do
                end if
                if( kfl_dmeth == 8 ) then
                   gpwal = 0.0_rp
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do igaus = 1,pgaus                      
                         gpwal(igaus) = gpwal(igaus) + walld(ipoin) * elmar(pelty) % shape(inode,igaus)
                      end do
                   end do
                end if
                !
                ! 1st order Cartesian derivatives GPCAR and GPVOL=dV=|J|*wgx
                ! VOELE: Element volume
                !
                voele = 0.0_rp
                do igaus = 1,pgaus     
                   call elmder(&
                        pnode,ndime,elmar(pelty) % deriv(1,1,igaus),& 
                        elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
                   gpvol(igaus) = elmar(pelty) % weigp(igaus) * gpdet  
                   voele        = voele + gpvol(igaus)
                end do
                !
                ! HLENG and element aspect ratio ASELE
                !
                call elmlen(&
                     ndime,pnode,elmar(pelty) % dercg,xjacm,elcod,&
                     hnatu(pelty),hleng)
                asele = hleng(1) / hleng(ndime)
                !
                ! Compute element matrix ELMAT 
                !
                call elmsmo(                                               &
                     pnode,pgaus,pelty,pevat,lelch(ielem),kfl_dmeth,       & 
                     kfl_elfix,gpcar,gpvol,elcod,elbve,hleng,voele,asmin,  &
                     asmax,asele,vomin,vomax,kfl_matdi,elmat,elrhs,gpwal,  &
                     defor_param)
                !
                ! Assemble
                !
                call matrix_assemble_element_RHS(&
                     solve_sol(1) % ndofn,solve_sol(1) % ndofn,pnode,&
                     lnods(:,ielem),elrhs,rhsid)
                call matrix_assemble_element_matrix_to_CSR(&
                     kfl_element_to_csr,solve_sol(1) % ndofn,pnode,pnode*solve_sol(1) % ndofn,&
                     ielem,lnods(:,ielem),elmat,solve_sol(1) % ia,&
                     solve_sol(1) % ja,amatr,lezdo)       
             end if
          end do

       end if

       call cputim(time2)
       cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time2 - time1
       call ker_timeline('END_ASSEMBLY')

       !-------------------------------------------------------------------
       ! 
       ! Impose boundary condition and solve system
       !
       !-------------------------------------------------------------------

       do ipoin = 1,npoin
          itotn = (ipoin-1) * ndime

          do idime = 1,ndime
             itotn = itotn + 1
             unkno(itotn) = xfact * bvess(idime,ipoin)
          end do

       end do
       

       solve_sol(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED


       call solver_solve(solve_sol,amatr,rhsid,unkno)
 
       !-------------------------------------------------------------------
       ! 
       ! When using increments, accumulate displacement
       !
       !-------------------------------------------------------------------

       if( ndefo > 1 ) then
          do ipoin = 1,npoin
             do idime = 1,ndime
                dispm_inc(idime,ipoin) = dispm_inc(idime,ipoin) + unkno( (ipoin-1)*ndime+idime )
             end do
          end do
       end if

    end do

    !-------------------------------------------------------------------
    ! 
    ! Put total deformation in UNKNO
    !
    !-------------------------------------------------------------------

    if( ndefo > 1 .and. INOTMASTER ) then
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotn = (ipoin-1) * ndime + idime
             unkno(itotn) = dispm_inc(idime,ipoin) 
          end do
       end do
       call memory_deallo(mem_modul(1:2,modul),'DISPM_INC','deform_deform',dispm_inc)
    end if
    
  end subroutine deform_deform

  subroutine elmsmo(                                  &
       pnode,pgaus,pelty,pevat,lelch,kfl_dmeth,       &
       kfl_elfix,gpcar,gpvol,elcod,elbve,hleng,       &
       voele,asmin,asmax,asele,vomin,vomax,kfl_matdi, &
       elmat,elrhs,gpwal,defor_param)
    !----------------------------------------------------------------------
    !****f* domain/ale_elmsmo
    ! NAME 
    !    ale_elmsmo
    ! DESCRIPTION
    !    Compute the elemental weighted Laplacian matrix
    !    ( alpha * grad X , grad v )
    ! USES
    ! USED BY
    !    smooth
    !***
    !----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: pelty
    integer(ip), intent(in)    :: lelch
    integer(ip), intent(in)    :: kfl_dmeth
    integer(ip), intent(in)    :: kfl_elfix(ndime,*)  
    real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpvol(pgaus)
    real(rp),    intent(in)    :: elcod(ndime,pnode)
    real(rp),    intent(in)    :: elbve(ndime,pnode)
    real(rp),    intent(in)    :: gpwal(pgaus)
    real(rp),    intent(in)    :: defor_param
    real(rp),    intent(in)    :: hleng(ndime)
    real(rp),    intent(in)    :: voele
    real(rp),    intent(in)    :: asmin
    real(rp),    intent(in)    :: asmax
    real(rp),    intent(inout) :: asele 
    real(rp),    intent(in)    :: vomin
    real(rp),    intent(in)    :: vomax
    integer(ip), intent(in)    :: kfl_matdi
    real(rp),    intent(out)   :: elmat(ndime*pnode,ndime*pnode)
    real(rp),    intent(out)   :: elrhs(ndime*pnode)
    integer(ip)                :: inode,jnode,kdime,igaus,idime
    integer(ip)                :: jdime,pevat,idofn,jdofn
    real(rp)                   :: fact1,gpder(ndime,ndime)
    real(rp)                   :: f1,f2,gpdif(3)

    !----------------------------------------------------------------------
    !
    ! Determine stiffness
    !
    !----------------------------------------------------------------------

    if( kfl_dmeth == 1 ) then
       !
       ! Only smoothing
       !
       do idime = 1,ndime
          gpdif(idime) = hleng(idime)
       end do

    else if( kfl_dmeth == 2 ) then
       !
       ! Conserve small elements only
       !
       do idime = 1,ndime
          gpdif(idime) = (1.0_rp+(1.0_rp-vomin/vomax)/(voele/vomax))
       end do

    else if( kfl_dmeth == 4 ) then
       !
       ! Sort hleng: hleng(1)=max; hleng(ndime)=min
       ! HLENG(1)     = Max length
       ! HLENG(NDIME) = Min length
       !     
       do idime = 1,ndime
          gpdif(idime) = asele
       end do

    else if( kfl_dmeth == 5 ) then
       !
       ! Isotropic and uniform Laplacian
       !
       do idime = 1,ndime
          gpdif(idime) = 1.0_rp
       end do

    else if( kfl_dmeth == 6 ) then
       !
       ! Aspect ratio and element volume for HEX and PEN
       !     
       if( pelty /= QUA04 .and. pelty /= HEX08 .and. pelty /= PEN06 ) then
          asele = 1.0_rp
       end if

       do idime = 1,ndime
          !f1           = (1.0_rp+(1.0_rp-asmin/asmax)/(asele/asmax))**1.0_rp
          !f2           = (1.0_rp+(1.0_rp-vomin/vomax)/(voele/vomax))**1.0_rp
          f1           = asele
          f2           = vomax/voele
          gpdif(idime) = f1 * f2
       end do

    else if( kfl_dmeth == 7 ) then
       !
       ! Aspect ratio and element volume for all elements
       !     
       do idime = 1,ndime
          !f1           = (1.0_rp+(1.0_rp-asmin/asmax)/(asele/asmax))**1.0_rp
          !f2           = (1.0_rp+(1.0_rp-vomin/vomax)/(voele/vomax))**1.0_rp
          f1           = asele
          f2           = vomax/voele
          gpdif(idime) = f1 * f2
       end do
       
    else if( kfl_dmeth == 8 ) then
       !
       ! Distance to the wall (averaged over Gauss points)
       !
       gpdif = 0.0_rp
       do igaus = 1,pgaus
          gpdif(1) = gpdif(1) + gpwal(igaus)
       end do
       gpdif(1)       = real(pgaus,rp) / (gpdif(1)**defor_param)
       gpdif(2:ndime) = gpdif(1)

    end if

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    pevat = pnode * ndime
    do jdofn = 1,pevat
       do idofn = 1,pevat
          elmat(idofn,jdofn) = 0.0_rp
       end do
       elrhs(jdofn) = 0.0_rp
    end do

    !----------------------------------------------------------------------
    !
    ! Laplacian matrix: ( alpha * grad X , grad v )
    !
    !----------------------------------------------------------------------

    do igaus = 1,pgaus

       gpder = 0.0_rp
       do inode = 1,pnode
          do jdime = 1,ndime
             do idime = 1,ndime
                gpder(idime,jdime) = gpder(idime,jdime) &
                     + elcod(jdime,inode) * gpcar(idime,inode,igaus)
             end do
          end do
       end do

       do inode = 1,pnode
          do jnode = inode+1,pnode
             fact1 = 0.0_rp
             do kdime = 1,ndime
                fact1 = fact1 + gpcar(kdime,inode,igaus) &
                     &        * gpcar(kdime,jnode,igaus)
             end do
             fact1 = fact1 * gpvol(igaus)
             idofn = ( inode - 1 ) * ndime 
             jdofn = ( jnode - 1 ) * ndime 
             do idime = 1,ndime
                idofn = idofn + 1
                jdofn = jdofn + 1
                elmat(idofn,jdofn) = elmat(idofn,jdofn) + fact1 * gpdif(idime) 
                elmat(jdofn,idofn) = elmat(jdofn,idofn) + fact1 * gpdif(idime) 
             end do
          end do
          fact1 = 0.0_rp
          do kdime = 1,ndime
             fact1 = fact1 + gpcar(kdime,inode,igaus) &
                  &        * gpcar(kdime,inode,igaus)
          end do
          fact1 = fact1 * gpvol(igaus)
          idofn = ( inode - 1 ) * ndime 
          do idime = 1,ndime
             idofn = idofn + 1
             elmat(idofn,idofn) = elmat(idofn,idofn) + fact1 * gpdif(idime) 
          end do

       end do

    end do


    !----------------------------------------------------------------------
    !
    ! Extension elements
    !
    !----------------------------------------------------------------------

    if( lelch == ELEXT ) then
       do idofn = ndime+1,pevat
          do jdofn = 1,pevat
             elmat(idofn,jdofn) = 0.0_rp
          end do
          elrhs(idofn) = 0.0_rp
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Prescribe Laplacian on walls
    !
    !----------------------------------------------------------------------

    if( kfl_matdi == 0 .and. kfl_adj_prob == 0 ) then

       do inode = 1,pnode
          do idime = 1,ndime
             if( kfl_elfix(idime,inode) > 0 ) then
                idofn = ( inode - 1 ) * ndime + idime
                fact1 = elmat(idofn,idofn)
                do jdofn = 1,pevat
                   elrhs(jdofn)       = elrhs(jdofn) - elmat(jdofn,idofn) * elbve(idime,inode)
                   elmat(jdofn,idofn) = 0.0_rp
                   elmat(idofn,jdofn) = 0.0_rp
                end do
                elmat(idofn,idofn) = fact1
                elrhs(idofn)       = fact1 * elbve(idime,inode)
             end if
          end do 
       end do
    
    elseif( kfl_matdi == 0 .and. kfl_adj_prob == 1 ) then
    
       do inode = 1,pnode
          do idime = 1,ndime
             idofn = ( inode - 1 ) * ndime + idime
             elrhs(idofn)       = elbve(idime,inode)
             if( kfl_elfix(idime,inode) > 0 ) then
                fact1 = elmat(idofn,idofn)
                do jdofn = 1,pevat
                   elmat(idofn,jdofn) = 0.0_rp
                end do
                elmat(idofn,idofn) = fact1
                elrhs(idofn)       = fact1 * elbve(idime,inode)
             end if
          end do
       end do
    

    end if

    !----------------------------------------------------------------------
    !
    ! Transpose operation for adjoint
    !
    !----------------------------------------------------------------------
    
    if(kfl_adj_prob == 1 ) then
      elmat = transpose(elmat)
    endif
    
    
    
  end subroutine elmsmo

end module mod_ker_deform

