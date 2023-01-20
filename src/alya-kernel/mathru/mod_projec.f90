!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Projection_Toolbox
!> @{
!> @name    ToolBox for L2 projections
!> @file    mod_projec.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for matrix operations
!> @details ToolBox for matrix operations: fill in, etc.
!------------------------------------------------------------------------

module mod_projec

  use def_kintyp,         only : ip,rp,r3p,lg,elm
  use def_elmtyp,         only : ELEXT
  use def_elmtyp,         only : PYR05
  use def_master,         only : party,pardi,parki,parr2,IMASTER
  use def_master,         only : NPOIN_REAL_2DIM,INOTMASTER,ISLAVE
  use def_master,         only : INOTSLAVE,zeror,gevec,kfl_paral
  use def_kermod,         only : kfl_grpro
  use def_domain,         only : ndime,npoin,nelem,nnode,mnode,ntens
  use def_domain,         only : lnods,ltype,coord,vmass,elmar,vmasc
  use def_domain,         only : lexis,ngaus,kfl_naxis,lnodb,mnodb
  use def_domain,         only : lelbo,mgaus,lelch,nboun,ltypb,ndimb
  use def_domain,         only : lnnob,mgaub,memor_dom,mnode
  use def_domain,         only : npoin_own,lpoty
  use def_domain,         only : mesh_type
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_alloca_min
  use mod_memory,         only : memory_size
  use mod_ker_proper,     only : ker_proper
  use mod_bouder

  implicit none

  private

  integer(ip) :: knode
  integer(ip) :: ipoin,idime,inode,ielem,igaus,jnode
  integer(ip) :: pnode,pelty,pgaus,jdime
  real(rp)    :: gpdet,gpvol
  real(rp)    :: xjaci(9),xjacm(9),xfact
  integer(8)  :: memor(2)

  interface projec_elements_to_nodes
     module procedure projec_elements_to_nodes_r3p_1,&
          &           projec_elements_to_nodes_r3p_2,&
          &           projec_elements_to_nodes_rp1  ,&
          &           projec_elements_to_nodes_rp2
  end interface projec_elements_to_nodes

  interface projec_boundaries_to_nodes
     module procedure projec_boundaries_to_nodes_rp_1,&
          &           projec_boundaries_to_nodes_rp_2,&
          &           projec_boundaries_to_nodes_rp_3,&
          &           projec_boundaries_to_nodes_rp_21
  end interface projec_boundaries_to_nodes

  public :: projec_mass_conservation
  public :: projec_norm_symmetric_gradient_vector
  public :: projec_elements_to_nodes
  public :: projec_boundaries_to_nodes
  public :: calc_massb
  public :: projec_elements_to_nodes_unity_test
  
contains

  subroutine projec_norm_symmetric_gradient_vector(unkno,grunk)

    !----------------------------------------------------------------------
    !
    !> @author  Guillaume Houzeaux
    !> @date    12/06/2013
    !> @brief   Project the norm of a symmetric gradient of a Tensor
    !> @details Given ui, project P=SijSij, with Sij=1/2(dui/dxj+duj/dxi)
    !
    !----------------------------------------------------------------------

    real(rp),   intent(in)             :: unkno(ndime,npoin)
    real(rp),   intent(inout), pointer :: grunk(:)
    real(rp)                           :: elgra(ndime,ndime)
    real(rp)                           :: elunk(ndime,mnode)
    real(rp)                           :: gpcar(ndime,mnode)
    real(rp)                           :: elcod(ndime,mnode)
    real(rp)                           :: S11,S12,S22,S13,S23,S33
    if( INOTMASTER ) then
       !
       ! Allocate memory?
       !
       if( .not. associated(grunk) ) then
          call memory_alloca(memor,'GRUNK','mod_projec',grunk,npoin)
       end if
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          grunk(ipoin) = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       if( kfl_grpro == 0 ) then

          if( ndime == 3 ) then
             !
             ! 3D open rule
             !
             do ielem = 1,nelem
                pelty = ltype(ielem)
                if( pelty > 0 ) then
                   pnode = nnode(pelty)
                   pgaus = ngaus(pelty)
                   !
                   ! Treat extension elements
                   !
                   if( lelch(ielem) == ELEXT ) then
                      knode = 1
                   else
                      knode = pnode
                   end if
                   !
                   ! Gather vectors
                   !
                   do inode = 1,pnode
                      ipoin          = lnods(inode,ielem)
                      elcod(1,inode) = coord(1,ipoin)
                      elcod(2,inode) = coord(2,ipoin)
                      elcod(3,inode) = coord(3,ipoin)
                      elunk(1,inode) = unkno(1,ipoin)
                      elunk(2,inode) = unkno(2,ipoin)
                      elunk(3,inode) = unkno(3,ipoin)
                   end do
                   !
                   ! Loop over Gauss points
                   !
                   do igaus = 1,pgaus
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty) % weigp(igaus) * gpdet
                      !
                      ! Velocity strain rates
                      !
                      elgra(1,1) = 0.0_rp
                      elgra(1,2) = 0.0_rp
                      elgra(1,3) = 0.0_rp
                      elgra(2,1) = 0.0_rp
                      elgra(2,2) = 0.0_rp
                      elgra(2,3) = 0.0_rp
                      elgra(3,1) = 0.0_rp
                      elgra(3,2) = 0.0_rp
                      elgra(3,3) = 0.0_rp
                      do jnode = 1,pnode
                         elgra(1,1) = elgra(1,1) + gpcar(1,jnode) * elunk(1,jnode)
                         elgra(1,2) = elgra(1,2) + gpcar(1,jnode) * elunk(2,jnode)
                         elgra(1,3) = elgra(1,3) + gpcar(1,jnode) * elunk(3,jnode)
                         elgra(2,1) = elgra(2,1) + gpcar(2,jnode) * elunk(1,jnode)
                         elgra(2,2) = elgra(2,2) + gpcar(2,jnode) * elunk(2,jnode)
                         elgra(2,3) = elgra(2,3) + gpcar(2,jnode) * elunk(3,jnode)
                         elgra(3,1) = elgra(3,1) + gpcar(3,jnode) * elunk(1,jnode)
                         elgra(3,2) = elgra(3,2) + gpcar(3,jnode) * elunk(2,jnode)
                         elgra(3,3) = elgra(3,3) + gpcar(3,jnode) * elunk(3,jnode)
                      end do
                      !
                      ! Assembly
                      !
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                         S11   = 0.5_rp * ( elgra(1,1) + elgra(1,1) )
                         S22   = 0.5_rp * ( elgra(2,2) + elgra(2,2) )
                         S12   = 0.5_rp * ( elgra(1,2) + elgra(2,1) )
                         S33   = 0.5_rp * ( elgra(3,3) + elgra(3,3) )
                         S13   = 0.5_rp * ( elgra(1,3) + elgra(3,1) )
                         S23   = 0.5_rp * ( elgra(2,3) + elgra(3,2) )
                         grunk(ipoin) = grunk(ipoin) + xfact * ( &
                              &   S11*S11 + S12*S12 + S13*S13    &
                              & + S12*S12 + S22*S22 + S23*S23    &
                              & + S13*S13 + S23*S23 + S33*S33  )
                      end do
                   end do
                end if
             end do

          else
             !
             ! 2D open rule
             !
             do ielem = 1,nelem
                pelty = ltype(ielem)
                if( pelty > 0 ) then
                   pnode = nnode(pelty)
                   pgaus = ngaus(pelty)
                   !
                   ! Treat extension elements
                   !
                   if( lelch(ielem) == ELEXT ) then
                      knode = 1
                   else
                      knode = pnode
                   end if
                   !
                   ! Gather vectors
                   !
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                      elunk(1:ndime,inode) = unkno(1:ndime,ipoin)
                   end do
                   !
                   ! Loop over Gauss points (which are nodes)
                   !
                   do igaus = 1,pgaus
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty)%weigp(igaus)*gpdet
                      !
                      ! Velocity strain rates
                      !
                      elgra(1,1) = 0.0_rp
                      elgra(1,2) = 0.0_rp
                      elgra(2,1) = 0.0_rp
                      elgra(2,2) = 0.0_rp
                      do jnode = 1,pnode
                         elgra(1,1) = elgra(1,1) + gpcar(1,jnode) * elunk(1,jnode)
                         elgra(1,2) = elgra(1,2) + gpcar(1,jnode) * elunk(2,jnode)
                         elgra(2,1) = elgra(2,1) + gpcar(2,jnode) * elunk(1,jnode)
                         elgra(2,2) = elgra(2,2) + gpcar(2,jnode) * elunk(2,jnode)
                      end do
                      !
                      ! Assembly
                      !
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                         S11   = 0.5_rp * ( elgra(1,1) + elgra(1,1) )
                         S22   = 0.5_rp * ( elgra(2,2) + elgra(2,2) )
                         S12   = 0.5_rp * ( elgra(1,2) + elgra(2,1) )
                         grunk(ipoin) = grunk(ipoin) + xfact * ( &
                              &   S11*S11 + S12*S12              &
                              & + S12*S12 + S22*S22            )
                      end do
                   end do
                end if
             end do
          end if

       else
          !
          ! 3D close rule
          !
          if( ndime == 3 ) then

             do ielem = 1,nelem
                pelty = ltype(ielem)

                if( pelty > 0 ) then
                   pgaus = ngaus(pelty)
                   pnode = nnode(pelty)
                   !
                   ! Treat extension elements
                   !
                   if( lelch(ielem) == ELEXT ) then
                      knode = 1
                   else
                      knode = pnode
                   end if
                   !
                   ! Gather vectors
                   !
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1,inode) = coord(1,ipoin)
                      elcod(2,inode) = coord(2,ipoin)
                      elcod(3,inode) = coord(3,ipoin)
                      elunk(1,inode) = unkno(1,ipoin)
                      elunk(2,inode) = unkno(2,ipoin)
                      elunk(3,inode) = unkno(3,ipoin)
                   end do
                   !
                   ! Loop over Gauss points (which are nodes)
                   !
                   if( pelty /= PYR05 ) then

                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         call elmder(&
                              pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                              elcod,gpcar,gpdet,xjacm,xjaci)
                         gpvol = elmar(pelty)%weigc(inode)*gpdet
                         !
                         ! Velocity strain rates
                         !
                         elgra(1,1) = 0.0_rp
                         elgra(1,2) = 0.0_rp
                         elgra(1,3) = 0.0_rp
                         elgra(2,1) = 0.0_rp
                         elgra(2,2) = 0.0_rp
                         elgra(2,3) = 0.0_rp
                         elgra(3,1) = 0.0_rp
                         elgra(3,2) = 0.0_rp
                         elgra(3,3) = 0.0_rp
                         do jnode = 1,pnode
                            elgra(1,1) = elgra(1,1) + gpcar(1,jnode) * elunk(1,jnode)
                            elgra(1,2) = elgra(1,2) + gpcar(1,jnode) * elunk(2,jnode)
                            elgra(1,3) = elgra(1,3) + gpcar(1,jnode) * elunk(3,jnode)
                            elgra(2,1) = elgra(2,1) + gpcar(2,jnode) * elunk(1,jnode)
                            elgra(2,2) = elgra(2,2) + gpcar(2,jnode) * elunk(2,jnode)
                            elgra(2,3) = elgra(2,3) + gpcar(2,jnode) * elunk(3,jnode)
                            elgra(3,1) = elgra(3,1) + gpcar(3,jnode) * elunk(1,jnode)
                            elgra(3,2) = elgra(3,2) + gpcar(3,jnode) * elunk(2,jnode)
                            elgra(3,3) = elgra(3,3) + gpcar(3,jnode) * elunk(3,jnode)
                         end do
                         !
                         ! Assembly
                         !
                         xfact = gpvol
                         S11   = 0.5_rp * ( elgra(1,1) + elgra(1,1) )
                         S22   = 0.5_rp * ( elgra(2,2) + elgra(2,2) )
                         S12   = 0.5_rp * ( elgra(1,2) + elgra(2,1) )
                         S33   = 0.5_rp * ( elgra(3,3) + elgra(3,3) )
                         S13   = 0.5_rp * ( elgra(1,3) + elgra(3,1) )
                         S23   = 0.5_rp * ( elgra(2,3) + elgra(3,2) )
                         grunk(ipoin) = grunk(ipoin) + xfact * ( &
                              & + S11*S11 + S12*S12 + S13*S13    &
                              & + S12*S12 + S22*S22 + S23*S23    &
                              & + S13*S13 + S23*S23 + S33*S33  )

                      end do

                   else

                      do igaus = 1,pgaus
                         call elmder(&
                              pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                              elcod,gpcar,gpdet,xjacm,xjaci)
                         gpvol = elmar(pelty)%weigp(igaus)*gpdet
                         !
                         ! Velocity strain rates
                         !
                         elgra(1,1) = 0.0_rp
                         elgra(1,2) = 0.0_rp
                         elgra(1,3) = 0.0_rp
                         elgra(2,1) = 0.0_rp
                         elgra(2,2) = 0.0_rp
                         elgra(2,3) = 0.0_rp
                         elgra(3,1) = 0.0_rp
                         elgra(3,2) = 0.0_rp
                         elgra(3,3) = 0.0_rp
                         do jnode = 1,pnode
                            elgra(1,1) = elgra(1,1) + gpcar(1,jnode) * elunk(1,jnode)
                            elgra(1,2) = elgra(1,2) + gpcar(1,jnode) * elunk(2,jnode)
                            elgra(1,3) = elgra(1,3) + gpcar(1,jnode) * elunk(3,jnode)
                            elgra(2,1) = elgra(2,1) + gpcar(2,jnode) * elunk(1,jnode)
                            elgra(2,2) = elgra(2,2) + gpcar(2,jnode) * elunk(2,jnode)
                            elgra(2,3) = elgra(2,3) + gpcar(2,jnode) * elunk(3,jnode)
                            elgra(3,1) = elgra(3,1) + gpcar(3,jnode) * elunk(1,jnode)
                            elgra(3,2) = elgra(3,2) + gpcar(3,jnode) * elunk(2,jnode)
                            elgra(3,3) = elgra(3,3) + gpcar(3,jnode) * elunk(3,jnode)
                         end do
                         !
                         ! Assembly
                         !
                         do inode = 1,knode
                            ipoin = lnods(inode,ielem)
                            xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                            S11   = 0.5_rp * ( elgra(1,1) + elgra(1,1) )
                            S22   = 0.5_rp * ( elgra(2,2) + elgra(2,2) )
                            S12   = 0.5_rp * ( elgra(1,2) + elgra(2,1) )
                            S33   = 0.5_rp * ( elgra(3,3) + elgra(3,3) )
                            S13   = 0.5_rp * ( elgra(1,3) + elgra(3,1) )
                            S23   = 0.5_rp * ( elgra(2,3) + elgra(3,2) )
                            grunk(ipoin) = grunk(ipoin) + xfact * ( &
                                 & + S11*S11 + S12*S12 + S13*S13    &
                                 & + S12*S12 + S22*S22 + S23*S23    &
                                 & + S13*S13 + S23*S23 + S33*S33  )
                         end do
                      end do

                   end if

                end if
             end do

          else
             !
             ! 2D close rule
             !
             do ielem = 1,nelem
                pelty = ltype(ielem)
                if( pelty > 0 ) then
                   pnode = nnode(pelty)
                   pgaus = ngaus(pelty)
                   !
                   ! Treat extension elements
                   !
                   if( lelch(ielem) == ELEXT ) then
                      knode = 1
                   else
                      knode = pnode
                   end if
                   !
                   ! Gather vectors
                   !
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elcod(1:ndime,inode)=coord(1:ndime,ipoin)
                      elunk(1:ndime,inode)=unkno(1:ndime,ipoin)
                   end do
                   !
                   ! Loop over Gauss points (which are nodes)
                   !
                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty)%weigc(inode)*gpdet
                      !
                      ! Velocity strain rates
                      !
                      elgra(1,1) = 0.0_rp
                      elgra(1,2) = 0.0_rp
                      elgra(2,1) = 0.0_rp
                      elgra(2,2) = 0.0_rp
                      do jnode = 1,pnode
                         elgra(1,1) = elgra(1,1) + gpcar(1,jnode) * elunk(1,jnode)
                         elgra(1,2) = elgra(1,2) + gpcar(1,jnode) * elunk(2,jnode)
                         elgra(2,1) = elgra(2,1) + gpcar(2,jnode) * elunk(1,jnode)
                         elgra(2,2) = elgra(2,2) + gpcar(2,jnode) * elunk(2,jnode)
                      end do
                      !
                      ! Assembly
                      !
                      xfact = gpvol
                      S11   = 0.5_rp * ( elgra(1,1) + elgra(1,1) )
                      S22   = 0.5_rp * ( elgra(2,2) + elgra(2,2) )
                      S12   = 0.5_rp * ( elgra(1,2) + elgra(2,1) )
                      grunk(ipoin) = grunk(ipoin) + xfact * ( &
                           & + S11*S11 + S12*S12              &
                           & + S12*S12 + S22*S22            )
                   end do
                end if
             end do

          end if

       end if
       !
       ! Parall and periodicity
       !
       Call rhsmod(1_ip,grunk)
       !
       ! Solve diagonal system
       !
       if( kfl_grpro == 0 ) then
          do ipoin = 1,npoin
             grunk(ipoin) = grunk(ipoin) / vmass(ipoin)
          end do
       else
          do ipoin = 1,npoin
             grunk(ipoin) = grunk(ipoin) / vmasc(ipoin)
          end do
       end if
    end if

  end subroutine projec_norm_symmetric_gradient_vector

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    12/09/2014
  !> @brief   Conserve mass
  !> @details Given a velocity vector u*, find u so that mass is conserved.
  !>          Solve the following constrained system:
  !>          \verbatim
  !>           +-
  !>           | Minimize int_S | u - u* |^2 dS
  !>           | Under    f(u) = r
  !>           +-
  !>                        ------
  !>           \partial W = S U S'
  !>
  !>                  S'
  !>              +--------+
  !>              |        |
  !>           S' |   W    | S
  !>              |        |
  !>              +--------+
  !>                  S'
  !>
  !>           where f(u) = int_S(u.n)ds and
  !>           Total mass: r =  x - int_S'(u.n)ds
  !>           Local mass: r =  x
  !>
  !>           In algebraic sense:
  !>           +-
  !>           | Minimize int_S | Nu - Nu* |^2 dS
  !>           | Under    R^t.u = r
  !>           +-
  !>
  !>           +-         -+ +-      -+   +-    -+
  !>           |  M    -R  | |    u   |   | M u* |  (1)
  !>           |           | |        | = |      |
  !>           | R^t    0  | | lambda |   |  r   |  (2)
  !>           +-         -+ +-      -+   +-    -+
  !>
  !>           which solution is obtained by computing lambda by doing
  !>           R^t * (1) - (2) and substitute the result in (1):
  !>           lambda = (R^t.M^-1.R)^-1 (r - R^t u*)
  !>           u      = u* + (M^-1.R) * lambda
  !>
  !>           To impose Dirichlet, eliminate b.c. in second equation:
  !>           R^t u = r - R^t u* and then put R|_Sd = 0 so that
  !>           first equation is u = u* on Sd
  !>
  !>           Note that:
  !>           - (r-R^tu*) is the constraint such that u=u* if it is fulfilled
  !>           - lambda = 0 if u=u*
  !>           - (R^t.M^-1.R)^-1 is a scalar
  !>          \endverbatim
  !>
  !----------------------------------------------------------------------

  subroutine projec_mass_conservation(unkno,lboun,lpoin,what,xmass_opt,where_opt,given_normal_opt,ERROR)

    real(rp),     intent(inout), pointer           :: unkno(:,:)
    logical(lg),  intent(in),    pointer           :: lboun(:)
    real(rp),     intent(in),            optional  :: given_normal_opt(ndime)
    logical(lg),  intent(in),   pointer, optional  :: lpoin(:)
    character(*), intent(in),            optional  :: what
    real(rp),     intent(in),            optional  :: xmass_opt
    character(*), intent(in),            optional  :: where_opt    
    integer(ip),  intent(out),           optional  :: ERROR
    integer(ip)                                    :: iboun,pblty,pnodb,pgaub,igaub
    integer(ip)                                    :: inodb,ipoin,idime,ndim1,dummi
    real(rp)                                       :: baloc(ndime,ndime),gbsur,gbdet
    real(rp)                                       :: bocod(ndime,mnodb)
    real(rp)                                       :: elcod(ndime,mnode)
    real(rp)                                       :: rr,ss,xmass,lambda,Rtu,rrr
    real(rp)                                       :: xfact,gbden(mgaub)
    real(rp)                                       :: uloc(ndime), ulocversor(ndime), ugivector(ndime), tanvector(ndime)
    real(rp)                                       :: ulocnorm, cosang, ugivenorm,unkno_tmp(1)
    real(rp),                   pointer            :: restr(:,:)
    real(rp),                   pointer            :: massb(:)
    real(rp),                   pointer            :: MinvR(:,:)
    logical(lg),                pointer            :: kpoin(:)
    logical(lg)                                    :: global_mass_normal
    logical(lg)                                    :: mass_local_global
    character(50)                                  :: where_message
    !
    ! By default, restriction is int_S u.n dS = 0
    !
    global_mass_normal = .false.
    mass_local_global  = .false.
    !
    ! Optional arguments
    !
    if( present(what) ) then
       if(      trim(what) == 'LOCAL MASS' ) then
          global_mass_normal = .false.
          mass_local_global  = .true.
       else if( trim(what) == 'TOTAL MASS' .or. trim(what) == 'GLOBAL MASS' ) then
          global_mass_normal = .true.
          mass_local_global  = .true.
       else if( trim(what) == 'LOCAL NORMAL' ) then
          global_mass_normal = .false.
          mass_local_global  = .false.
       else if( trim(what) == 'TOTAL NORMAL' .or. trim(what) == 'GLOBAL NORMAL' ) then
          global_mass_normal = .true.
          mass_local_global  = .false.
       end if
    end if
    !
    ! Prescribed flow rate
    !
    if( present(xmass_opt) ) then
       xmass = xmass_opt
    else
       xmass = 0.0_rp
    end if
    !
    ! Error
    !
    if( present(ERROR) ) ERROR = 0    
    !
    ! Communicator to use
    !
    if( present(where_opt) ) then
       where_message = trim(where_opt)
    else
       where_message = 'IN MY CODE'
    end if

    nullify( restr      )
    nullify( massb      ) 
    nullify( MinvR      )
    nullify( kpoin      )
    gbden = 1.0_rp
    rr    = xmass
    rrr   = 0.0_rp
    !
    ! BOundary condition on IPOIN:
    ! LPOIN(IPOIN) = .FALSE. : Dirichlet is prescribed using input value
    !              = .TRUE.  : Value should be update
    !
    if( IMASTER ) then

       if( .not. present(lpoin) ) then
          call memory_alloca_min(memor_dom,'KPOIN','projec_mass_conservation',kpoin)
       end if
       call memory_alloca_min(memor_dom,'RESTR','projec_mass_conservation',restr)
       call memory_alloca_min(memor_dom,'MASSB','projec_mass_conservation',massb)
       call memory_alloca_min(memor_dom,'MINVR','projec_mass_conservation',MinvR)


    else if( INOTMASTER ) then

       if( present(lpoin) ) then
          kpoin => lpoin
       else if( associated(lboun) ) then
          allocate( kpoin(npoin) )
          kpoin = .false.
          do iboun = 1,nboun
             if( lboun(iboun) ) then
                do inodb = 1,lnnob(iboun)
                   ipoin = lnodb(inodb,iboun)
                   kpoin(ipoin) = .true.
                end do
             end if
          end do
       end if

       ndim1 = ndime + 1
       call memory_alloca(memor_dom,'RESTR','projec_mass_conservation',restr,ndime,npoin)  ! R = restriction operator
       call memory_alloca(memor_dom,'MASSB','projec_mass_conservation',massb,npoin)        ! M = boundary mass
       call memory_alloca(memor_dom,'MINVR','projec_mass_conservation',MinvR,ndime,npoin)  ! M^-1.R
       !
       ! Loop over marked boundaries LBOUN(:) = .TRUE.
       !
       do iboun = 1,nboun

          if( lboun(iboun) .or. global_mass_normal ) then

             pblty = abs(ltypb(iboun))
             pnodb = nnode(pblty)
             pgaub = ngaus(pblty)
             ielem = lelbo(iboun)
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
             end do
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do

             if( mass_local_global ) then
                call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
             end if

             do igaub = 1,pgaub
                call bouder(&
                     pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),&
                     bocod,baloc,gbdet)
                call chenor(pnode,baloc,bocod,elcod)                                                         ! Check normal
                gbsur = elmar(pblty) % weigp(igaub) * gbdet                                                  ! |J| ds
                do inodb = 1,pnodb
                   ipoin = lnodb(inodb,iboun)
                   xfact = gbsur * elmar(pblty) % shape(inodb,igaub)                                         ! Ni |J| ds
                   restr(1:ndime,ipoin) = restr(1:ndime,ipoin) + xfact * baloc(1:ndime,ndime) * gbden(igaub) ! Restriction R
                   massb(ipoin)         = massb(ipoin)         + xfact                                       ! Mass matrix M
                end do
                if( global_mass_normal ) then                                                                ! Remaing mass
                   do inodb = 1,pnodb
                      ipoin = lnodb(inodb,iboun)
                      do idime = 1,ndime
                         rrr = rrr + gbsur * elmar(pblty) % shape(inodb,igaub) * baloc(idime,ndime) &
                              * unkno(idime,ipoin) * gbden(igaub)
                      end do
                   end do
                end if
             end do

          end if

       end do
       !
       ! Parallel exchange and periodicity
       !
       call PAR_INTERFACE_NODE_EXCHANGE(restr,'SUM','IN MY CODE')
       call PAR_INTERFACE_NODE_EXCHANGE(massb,'SUM','IN MY CODE')
    end if
    !
    ! Compute M^-1.R for nodes located on S. Result is zero on prescribed nodes
    !
    do ipoin = 1,npoin
       if( massb(ipoin) > 0.0_rp .and. kpoin(ipoin) ) then
          MinvR(1:ndime,ipoin) = restr(1:ndime,ipoin) / massb(ipoin)
       end if
    end do
    !
    ! rr = r - R^t.u
    ! ss = R^t.M^-1.R
    !
    if( associated(unkno) ) then
       call prodxy(ndime,npoin,restr,unkno,Rtu)
    else
       call prodxy(ndime,npoin,restr,unkno_tmp,Rtu)       
    end if
    call prodxy(ndime,npoin,restr,MinvR,ss)

    if( INOTMASTER ) then
       !
       ! Add remaining mass
       !
       if( global_mass_normal ) then
          call PAR_SUM(rrr,where_message)
          rr = rr + rrr
       end if

       if( ss /= 0.0_rp ) then
          !
          ! Lagrange multiplier LAMBDA = ( r - R^t.u ) / ( R^t.M^-1.R )^-1
          !
          lambda = ( rr - Rtu ) / ss
          !
          ! Actualize solution u = u* + M^-1.R * lambda
          !
          do ipoin = 1,npoin
             if( kpoin(ipoin) ) then
                unkno(1:ndime,ipoin) = unkno(1:ndime,ipoin) + MinvR(1:ndime,ipoin) * lambda

                !
                ! Correct if there is a given_normal:
                !
                ! uloc_i : unkno
                ! ulocnorm= |uloc_i|
                ! ulocversor_i = uloc_i / ulocnorm  : normalized vector unkno 
                ! n_i : (normalized) given_normal_opt
                ! 
                ! cos = ulocversor_i n_i
                ! ugivenorm = |u| / abs(cos)
                ! ugivector_i = ugivenorm * n_i : velocity along the given normal
                !
                ! tanvector_i = ugivector_i - uloc_i : what i have to add to u_i to obtain q_i
                ! 
                ! uloc_i = uloc_i + tanvector_i : final corrected unkno
                !
                if( present(given_normal_opt) ) then

                   uloc(1:ndime)= unkno(1:ndime,ipoin)
                   ulocnorm= 0.0_rp
                   do idime=1,ndime
                      ulocnorm= ulocnorm + uloc(idime)*uloc(idime)
                   end do
                   ulocnorm= sqrt(ulocnorm)
                   ulocversor= 0.0_rp
                   if (ulocnorm > 0.0_rp) then
                      ulocversor= uloc / ulocnorm
                   end if
                   cosang= 0.0_rp
                   do idime=1,ndime
                      cosang= cosang + ulocversor(idime)*given_normal_opt(idime)
                   end do
                   ugivenorm= 0.0_rp
                   if (abs(cosang) > 0.0_rp) ugivenorm= ulocnorm / abs(cosang)
                   ugivector= ugivenorm * given_normal_opt
                   tanvector= ugivector - uloc
                   uloc = uloc + tanvector

                   unkno(1:ndime,ipoin)= uloc(1:ndime)

                end if


             end if
          end do

       else

          if( present(ERROR) ) then
             ERROR = 1
          else
             call runend('PROJEC_MASS_CONSERVATION: WRONG LAGRANGE MULTIPLIER')
          end if

       end if

       if( .not. present(lpoin) ) then
          call memory_deallo(memor_dom,'KPOIN','projec_mass_conservation',kpoin)
       end if
       call memory_deallo(memor_dom,'RESTR','projec_mass_conservation',restr)
       call memory_deallo(memor_dom,'MASSB','projec_mass_conservation',massb)
       call memory_deallo(memor_dom,'MinvR','projec_mass_conservation',MinvR)

    end if

    if( present(ERROR) ) call PAR_MAX(ERROR)


  end subroutine projec_mass_conservation

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    09/12/2015
  !> @brief   Boundary projection
  !> @details Project an element type array of constant values
  !>          on the nodes
  !>
  !----------------------------------------------------------------------

  subroutine projec_elements_to_nodes_rp1(unkno,xvalu)
    implicit none
    real(rp),     pointer, intent(in)    :: unkno(:)
    real(rp),     pointer, intent(inout) :: xvalu(:)
    call projec_elements_to_nodes_rp(1_ip,unkno,xvalu)
  end subroutine projec_elements_to_nodes_rp1

  subroutine projec_elements_to_nodes_rp2(unkno,xvalu)
    implicit none
    real(rp),     pointer, intent(in)    :: unkno(:,:)
    real(rp),     pointer, intent(inout) :: xvalu(:,:)
    integer(ip)                          :: ndofn
    ndofn = memory_size(unkno,1_ip)
    if( memory_size(unkno,1_ip) /= memory_size(xvalu,1_ip) ) call runend('PROJEC_ELEMENTS_TO_NODES: WRONG SIZE')
    call projec_elements_to_nodes_rp(ndofn,unkno,xvalu)
  end subroutine projec_elements_to_nodes_rp2

  subroutine projec_elements_to_nodes_rp(ndofn,unkno,xvalu)

    implicit none
    integer(ip),  intent(in)  :: ndofn
    real(rp),     intent(in)  :: unkno(ndofn,*)
    real(rp),     intent(out) :: xvalu(ndofn,*)
    integer(ip)               :: ipoin,inode,ielem,igaus
    integer(ip)               :: pnode,pelty,knode
    real(rp)                  :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                  :: xjaci(9),xjacm(9)
    real(rp)                  :: elcod(ndime,mnode)
    real(rp)                  :: elunk(ndofn,mnode)

    if( INOTMASTER ) then
       !
       ! Errors
       !
       if( kfl_naxis == 1 ) then
          call runend('MOD_GRADIE: NOT CODED')
       end if
       !
       ! Initialization
       !
       xvalu(1:ndofn,1:npoin) = 0.0_rp
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 ) then
             !
             ! Gather vectors
             !
             pnode                  = nnode(pelty)
             elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
             elunk(1:ndofn,1:pnode) = 0.0_rp
             !
             ! Extension
             !
             if( lelch(ielem) == ELEXT ) then
                knode = 1
             else
                knode = pnode
             end if
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,ngaus(pelty)
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm
                do inode = 1,knode
                   elunk(1:ndofn,inode) = elunk(1:ndofn,inode) + gpvol * elmar(pelty) % shape(inode,igaus)
                end do

             end do gauss_points

             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                xvalu(1:ndofn,ipoin) = xvalu(1:ndofn,ipoin) + elunk(1:ndofn,inode) * unkno(1:ndofn,ielem)
             end do

          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(ndofn,xvalu)
       !
       ! Solve system
       !
       do ipoin = 1,npoin
          xvalu(1:ndofn,ipoin) = xvalu(1:ndofn,ipoin) / vmass(ipoin)
       end do

    end if

  end subroutine projec_elements_to_nodes_rp

  subroutine projec_elements_to_nodes_r3p_1(unkno,xvalu)

    type(r3p),    pointer, intent(in)    :: unkno(:)
    real(rp),     pointer, intent(inout) :: xvalu(:)
    integer(ip)                          :: ipoin,inode,ielem,igaus
    integer(ip)                          :: pnode,pelty,knode
    real(rp)                             :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                             :: xjaci(9),xjacm(9)
    real(rp)                             :: elcod(ndime,mnode)
    real(rp)                             :: elval(mnode)
    real(rp)                             :: gpunk(mgaus)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin = 1,npoin          
          xvalu(ipoin) = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 ) then
             !
             ! Gather 
             !
             pnode                  = nnode(pelty)
             pgaus                  = ngaus(pelty)
             elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
             gpunk(1:pgaus)         = unkno(ielem) % a(1,1:pgaus,1)
             elval                  = 0.0_rp
             !
             ! Extension
             !
             if( lelch(ielem) == ELEXT ) then
                knode = 1
             else
                knode = pnode
             end if
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm
                do inode = 1,pnode
                   elval(inode) = elval(inode) + &
                        gpunk(igaus) * gpvol * elmar(pelty) % shape(inode,igaus)
                end do
            end do gauss_points
            !
            ! Scatter
            !
            do inode = 1,pnode
               ipoin = lnods(inode,ielem)
               xvalu(ipoin) = xvalu(ipoin) +  elval(inode)
            end do
            
         end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(1_ip,xvalu)
       !
       ! Solve system
       !
       do ipoin = 1,npoin
          xvalu(ipoin) = xvalu(ipoin) / vmass(ipoin)
       end do

    end if

  end subroutine projec_elements_to_nodes_r3p_1

  subroutine projec_elements_to_nodes_r3p_2(ndofn,unkno,xvalu)

    integer(ip),           intent(in)    :: ndofn
    type(r3p),    pointer, intent(in)    :: unkno(:)
    real(rp),     pointer, intent(inout) :: xvalu(:,:)
    integer(ip)                          :: ipoin,inode,ielem,igaus
    integer(ip)                          :: pnode,pelty,knode
    real(rp)                             :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                             :: xjaci(9),xjacm(9)
    real(rp)                             :: elcod(ndime,mnode)
    real(rp)                             :: elval(ndofn,mnode)
    real(rp)                             :: gpunk(ndofn,mgaus)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin = 1,npoin          
          xvalu(1:ndofn,ipoin) = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 ) then
             !
             ! Gather 
             !
             pnode                  = nnode(pelty)
             pgaus                  = ngaus(pelty)
             elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
             gpunk(1:ndofn,1:pgaus) = unkno(ielem) % a(1:ndofn,1:pgaus,1)
             elval                  = 0.0_rp
             !
             ! Extension
             !
             if( lelch(ielem) == ELEXT ) then
                knode = 1
             else
                knode = pnode
             end if
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm
                do inode = 1,pnode
                   elval(1:ndofn,inode) = elval(1:ndofn,inode) + &
                        gpunk(1:ndofn,igaus) * gpvol * elmar(pelty) % shape(inode,igaus)
                end do
            end do gauss_points
            !
            ! Scatter
            !
            do inode = 1,pnode
               ipoin = lnods(inode,ielem)
               xvalu(1:ndofn,ipoin) = xvalu(1:ndofn,ipoin) + elval(1:ndofn,inode)
            end do
            
         end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(ndofn,xvalu)
       !
       ! Solve system
       !
       do ipoin = 1,npoin
          xvalu(1:ndofn,ipoin) = xvalu(1:ndofn,ipoin) / vmass(ipoin)
       end do

    end if

  end subroutine projec_elements_to_nodes_r3p_2

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    09/12/2015
  !> @brief   Boundary projection
  !> @details Project an boundary array on the nodes
  !>
  !----------------------------------------------------------------------

  subroutine projec_boundaries_to_nodes_rp_1(xarra_in,meshe_in,elmar_in,xarra_out,boundary_mask)
    real(rp),        intent(in),   pointer            :: xarra_in(:)       !< Value on boundaries
    type(mesh_type), intent(in)                       :: meshe_in          !< Mesh structure
    type(elm),       intent(in),   pointer            :: elmar_in(:)       !< Integration structure
    real(rp),        intent(inout),pointer            :: xarra_out(:)      !< Projected value
    logical(lg),     intent(in),   pointer, optional  :: boundary_mask(:)
    real(rp)                                          :: xarra_tmp(1,1)

    if( INOTMASTER ) then

       if( associated(xarra_in) ) then
          call projec_boundaries_to_nodes_rp(1_ip,xarra_in,meshe_in,elmar_in,xarra_out,boundary_mask)
       else
          call projec_boundaries_to_nodes_rp(1_ip,xarra_tmp,meshe_in,elmar_in,xarra_out,boundary_mask)          
       end if

    end if

  end subroutine projec_boundaries_to_nodes_rp_1

  subroutine projec_boundaries_to_nodes_rp_2(xarra_in,meshe_in,elmar_in,xarra_out,boundary_mask)
    real(rp),        intent(in),   pointer            :: xarra_in(:,:)     !< Value on boundaries
    type(mesh_type), intent(in)                       :: meshe_in          !< Mesh structure
    type(elm),       intent(in),   pointer            :: elmar_in(:)       !< Integration structure
    real(rp),        intent(inout),pointer            :: xarra_out(:,:)    !< Projected value
    logical(lg),     intent(in),   pointer, optional  :: boundary_mask(:)
    integer(ip)                                       :: ndofn
    real(rp)                                          :: xarra_tmp(1,1)

    if( INOTMASTER ) then

       ndofn = memory_size(xarra_in,1_ip)
       if( associated(xarra_in) ) then
          call projec_boundaries_to_nodes_rp(ndofn,xarra_in ,meshe_in,elmar_in,xarra_out,boundary_mask)
       else
          call projec_boundaries_to_nodes_rp(ndofn,xarra_tmp,meshe_in,elmar_in,xarra_out,boundary_mask)
       end if

    end if

  end subroutine projec_boundaries_to_nodes_rp_2

  subroutine projec_boundaries_to_nodes_rp_21(xarra_in,meshe_in,elmar_in,xarra_out,boundary_mask)
    real(rp),        intent(in),   pointer            :: xarra_in(:,:)     !< Value on boundaries
    type(mesh_type), intent(in)                       :: meshe_in          !< Mesh structure
    type(elm),       intent(in),   pointer            :: elmar_in(:)       !< Integration structure
    real(rp),        intent(inout),pointer            :: xarra_out(:)      !< Projected value
    logical(lg),     intent(in)  , pointer, optional  :: boundary_mask(:)
    integer(ip)                                       :: ndofn
    real(rp)                                          :: xarra_tmp(1,1)

    if( INOTMASTER ) then

       ndofn = memory_size(xarra_in,1_ip)
       !if( ndofn /= 1 ) call runend('projec_boundaries_to_nodes_rp_21: DONT KNOW WHAT TO DO')

       if( associated(xarra_in) ) then
          call projec_boundaries_to_nodes_rp(1_ip,xarra_in,meshe_in,elmar_in,xarra_out,boundary_mask)
       else          
          call projec_boundaries_to_nodes_rp(1_ip,xarra_tmp,meshe_in,elmar_in,xarra_out,boundary_mask)
       end if

    end if

  end subroutine projec_boundaries_to_nodes_rp_21

  subroutine projec_boundaries_to_nodes_rp(ndofn,xarra_in,meshe_in,elmar_in,xarra_out,boundary_mask)
    integer(ip),     intent(in)            :: ndofn
    real(rp),        intent(in)            :: xarra_in(ndofn,*)     !< Value on boundaries
    type(mesh_type), intent(in)            :: meshe_in              !< Mesh structure
    type(elm),       intent(in)            :: elmar_in(*)           !< Integration structure
    real(rp),        intent(out)           :: xarra_out(ndofn,*)    !< Projected value
    logical(lg),     intent(in), pointer, optional  :: boundary_mask(:)
    integer(ip)                            :: iboun,pblty,pnodb
    integer(ip)                            :: pgaub,igaub
    integer(ip)                            :: inodb,ipoin
    real(rp)                               :: baloc(ndime,ndime)
    real(rp)                               :: gbsur,gbdet,xvalu(ndofn),epsil
    real(rp)                               :: bocod(ndime,meshe_in % mnode)
    real(rp)                               :: elmas(mnode)
    real(rp)                               :: elxar(ndofn,mnode)
    logical(lg)                            :: compute_iboun
    real(rp),                   pointer    :: massb(:)

    if( INOTMASTER ) then

       nullify( massb )
       allocate( massb(meshe_in % npoin) )
       epsil = epsilon(1.0_rp)

       do ipoin = 1,meshe_in % npoin
          xarra_out(1:ndofn,ipoin) = 0.0_rp
          massb(ipoin)             = 0.0_rp
       end do

       do iboun = 1,meshe_in % nboun

          if( present(boundary_mask) ) then
             compute_iboun = boundary_mask(iboun)
          else
             compute_iboun = .true.
          end if

          if( compute_iboun ) then

             pblty = abs(meshe_in % ltypb(iboun))
             pgaub = elmar_in(pblty) % pgaus
             pnodb = meshe_in % lnnob(iboun)
             xvalu(1:ndofn) = xarra_in(1:ndofn,iboun)

             do inodb = 1,pnodb
                ipoin = meshe_in % lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = meshe_in % coord(1:ndime,ipoin)
             end do
             elmas(1:pnodb)         = 0.0_rp
             elxar(1:ndofn,1:pnodb) = 0.0_rp

             do igaub = 1,pgaub
                call bouder(&
                     pnodb,ndime,ndimb,elmar_in(pblty) % deriv(:,:,igaub),&
                     bocod,baloc,gbdet)
                gbsur = elmar_in(pblty) % weigp(igaub) * gbdet
                do inodb = 1,pnodb
                   elxar(1:ndofn,inodb) = elxar(1:ndofn,inodb) + gbsur * elmar_in(pblty) % shape(inodb,igaub) * xvalu(1:ndofn)
                   elmas(inodb) = elmas(inodb) + gbsur * elmar_in(pblty) % shape(inodb,igaub)
                end do
             end do

             do inodb = 1,pnodb
                ipoin            = meshe_in % lnodb(inodb,iboun)
                xarra_out(1:ndofn,ipoin) = xarra_out(1:ndofn,ipoin) + elxar(1:ndofn,inodb)
                massb(ipoin)     = massb(ipoin)     + elmas(inodb)
             end do
          end if
       end do
       !
       ! Parallel exchange and periodicity
       !
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,xarra_out,'SUM','IN MY CODE')
       call PAR_INTERFACE_NODE_EXCHANGE(1_ip, massb    ,'SUM','IN MY CODE')
       !
       ! Solve diagonal system
       !
       do ipoin = 1,meshe_in % npoin
          if( massb(ipoin) > 0.0_rp) then
             xarra_out(1:ndofn,ipoin) = xarra_out(1:ndofn,ipoin) / massb(ipoin)
          end if
       end do
       deallocate(massb)

    end if

  end subroutine projec_boundaries_to_nodes_rp

  subroutine projec_boundaries_to_nodes_rp_3(xarra_in,meshe_in,xarra_out,boundary_mask)

    real(rp),        intent(in),    pointer           :: xarra_in(:,:,:) !< Value on boundaries
    type(mesh_type), intent(in)                       :: meshe_in                !< Mesh structure
    real(rp),        intent(inout), pointer           :: xarra_out(:,:)       !< Projected value
    logical(lg),     intent(in),    pointer, optional :: boundary_mask(:)
    integer(ip)                                       :: iboun,pblty,pnodb
    integer(ip)                                       :: pgaub,igaub
    integer(ip)                                       :: inodb,ipoin,ndofn
    real(rp)                                          :: baloc(ndime,ndime)
    real(rp)                                          :: gbsur,gbdet,epsil
    real(rp)                                          :: bocod(ndime,meshe_in % mnode)
    real(rp)                                          :: elmas(mnode)
    real(rp),                   pointer               :: elxar(:,:)
    logical(lg)                                       :: compute_iboun
    real(rp),                   pointer               :: massb(:)
        
    ndofn = memory_size(xarra_in,1_ip)
    call PAR_MAX(ndofn)
    
    if( INOTMASTER ) then

       nullify ( massb )
       nullify ( elxar )
       allocate( massb(meshe_in % npoin) )
       allocate( elxar(ndofn,mnode) )
       epsil = epsilon(1.0_rp)

       do ipoin = 1,meshe_in % npoin
          xarra_out(1:ndofn,ipoin) = 0.0_rp
          massb(ipoin)             = 0.0_rp
       end do

       do iboun = 1,meshe_in % nboun

          if( present(boundary_mask) ) then
             compute_iboun = boundary_mask(iboun)
          else
             compute_iboun = .true.
          end if

          if( compute_iboun ) then

             pblty = abs(meshe_in % ltypb(iboun))
             pgaub = elmar(pblty) % pgaus
             pnodb = meshe_in % lnnob(iboun)

             do inodb = 1,pnodb
                ipoin = meshe_in % lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = meshe_in % coord(1:ndime,ipoin)
             end do
             elmas(1:pnodb)         = 0.0_rp
             elxar(1:ndofn,1:pnodb) = 0.0_rp

             do igaub = 1,pgaub
                call bouder(&
                     pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),&
                     bocod,baloc,gbdet)
                gbsur = elmar(pblty) % weigp(igaub) * gbdet
                do inodb = 1,pnodb
                   elxar(1:ndofn,inodb) = elxar(1:ndofn,inodb) + gbsur * elmar(pblty) % shape(inodb,igaub) * xarra_in(1:ndofn,igaub,iboun)
                   elmas(inodb) = elmas(inodb) + gbsur * elmar(pblty) % shape(inodb,igaub)
                end do
             end do

             do inodb = 1,pnodb
                ipoin            = meshe_in % lnodb(inodb,iboun)
                xarra_out(1:ndofn,ipoin) = xarra_out(1:ndofn,ipoin) + elxar(1:ndofn,inodb)
                massb(ipoin)     = massb(ipoin)     + elmas(inodb)
             end do
          end if
       end do
       !
       ! Parallel exchange and periodicity
       !
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,xarra_out,'SUM','IN MY CODE')
       call PAR_INTERFACE_NODE_EXCHANGE(1_ip, massb    ,'SUM','IN MY CODE')
       !
       ! Solve diagonal system
       !
       do ipoin = 1,meshe_in % npoin
          if( massb(ipoin) > 0.0_rp) then
             xarra_out(1:ndofn,ipoin) = xarra_out(1:ndofn,ipoin) / massb(ipoin)
          end if
       end do
       deallocate(massb)
       deallocate(elxar)
       
    end if

  end subroutine projec_boundaries_to_nodes_rp_3

  subroutine calc_massb(meshe_in,elmar_in,massb,boundary_mask)
    type(mesh_type), intent(in)            :: meshe_in                             !< Mesh structure
    type(elm),       intent(in)            :: elmar_in(*)                          !< Integration structure
    real(rp),        intent(out)           :: massb(meshe_in % npoin)              !< boundary mass matrix
    logical(lg),     intent(in), pointer, optional  :: boundary_mask(:)
    integer(ip)                            :: iboun,pblty,pnodb
    integer(ip)                            :: pgaub,igaub
    integer(ip)                            :: inodb,ipoin
    real(rp)                               :: baloc(ndime,ndime)
    real(rp)                               :: gbsur,gbdet,epsil
    real(rp)                               :: bocod(ndime,meshe_in % mnode)
    real(rp)                               :: elmas(mnode)
    logical(lg)                            :: compute_iboun

    if( INOTMASTER ) then

       epsil = epsilon(1.0_rp)

       do ipoin = 1,meshe_in % npoin
          massb(ipoin)             = 0.0_rp
       end do

       do iboun = 1,meshe_in % nboun

          if( present(boundary_mask) ) then
             compute_iboun = boundary_mask(iboun)
          else
             compute_iboun = .true.
          end if

          if( compute_iboun ) then

             pblty = abs(meshe_in % ltypb(iboun))
             pgaub = elmar_in(pblty) % pgaus
             pnodb = meshe_in % lnnob(iboun)

             do inodb = 1,pnodb
                ipoin = meshe_in % lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = meshe_in % coord(1:ndime,ipoin)
             end do
             elmas(1:pnodb)         = 0.0_rp

             do igaub = 1,pgaub
                call bouder(&
                     pnodb,ndime,ndimb,elmar_in(pblty) % deriv(:,:,igaub),&
                     bocod,baloc,gbdet)
                gbsur = elmar_in(pblty) % weigp(igaub) * gbdet
                do inodb = 1,pnodb
                   elmas(inodb) = elmas(inodb) + gbsur * elmar_in(pblty) % shape(inodb,igaub)
                end do
             end do

             do inodb = 1,pnodb
                ipoin            = meshe_in % lnodb(inodb,iboun)
                massb(ipoin)     = massb(ipoin)     + elmas(inodb)
             end do
          end if
       end do
       !
       ! Parallel exchange and periodicity
       !
       ! call PAR_INTERFACE_NODE_EXCHANGE(1_ip, massb    ,'SUM','IN MY CODE')
       ! guillaume me dijo usar rhsmod en lugar de PAR_INTERFACE_NODE_EXCHANGE - me lo arreglo !!!!
       ! dice que en projec_boundaries_to_nodes_rp tambien serái más correcto con rhsmod - pero ahi ambos dan bien aca no
       call rhsmod( 1_ip,massb)         

    end if

  end subroutine calc_massb

  subroutine projec_elements_to_nodes_unity_test()

    integer(ip)            :: ndofn,pgaus,ielem,igaus
    integer(ip)            :: inode,pnode,ipoin
    real(rp)               :: gpcoo(ndime),numer,denom,xx
    type(r3p),    pointer  :: unkno(:)
    type(r3p),    pointer  :: unkn1(:)
    real(rp),     pointer  :: xvalu(:,:)
    real(rp),     pointer  :: unkn2(:)
    real(rp),     pointer  :: xval1(:)
    !
    ! Testing projec_elements_to_nodes_r3p_2
    ! Solution= ( x+2y+3z, 2(x+2y+3z), 3(x+2y+3z) )
    !
    ndofn = ndime
    if( INOTMASTER ) then
       allocate(unkno(nelem))
       allocate(xvalu(ndofn,npoin))
       do ielem = 1,nelem
          pelty = ltype(ielem)
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)         
          allocate(unkno(ielem) % a(ndofn,pgaus,1))
          do igaus = 1,pgaus
             gpcoo = 0.0_rp
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                gpcoo(1:ndime) = gpcoo(1:ndime) + elmar(pelty) % shape(inode,igaus) * coord(1:ndime,ipoin)
             end do
             xx = 0.0_rp
             do idime = 1,ndime
                xx = xx + real(idime,rp) * gpcoo(idime)
             end do
             do idime = 1,ndime
                unkno(ielem) % a(idime,igaus,1) = real(idime,rp) * xx
             end do
          end do
       end do
       call projec_elements_to_nodes(ndofn,unkno,xvalu)
    end if
    numer = 0.0_rp
    denom = 0.0_rp
    do ipoin = 1,npoin_own
       if(lpoty(ipoin)==0) then
          do idime = 1,ndime
             xx = 0.0_rp
             do jdime = 1,ndime
                xx = xx + real(jdime,rp) * coord(jdime,ipoin)
             end do
             xx = real(idime,rp) * xx
             numer = numer + ( xvalu(idime,ipoin) - xx )**2
             denom = denom + ( xx )**2
          end do
       end if
    end do
    call PAR_SUM(numer)
    call PAR_SUM(denom)
    numer = sqrt(numer/denom)
    if( IMASTER ) print*,'Error 1=',numer
    !
    ! Testing projec_elements_to_nodes_r3p_1
    ! Solution = x+2y+3z
    !
    if( INOTMASTER ) then
       allocate(unkn1(nelem))
       allocate(xval1(npoin))
       do ielem = 1,nelem
          pelty = ltype(ielem)
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)         
          allocate(unkn1(ielem) % a(1,pgaus,1))
          do igaus = 1,pgaus
             gpcoo = 0.0_rp
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                gpcoo(1:ndime) = gpcoo(1:ndime) + elmar(pelty) % shape(inode,igaus) * coord(1:ndime,ipoin)
             end do
             xx = 0.0_rp
             do idime = 1,ndime
                xx = xx + real(idime,rp) * gpcoo(idime)
             end do
             unkn1(ielem) % a(1,igaus,1) = xx
          end do
       end do
       call projec_elements_to_nodes(unkn1,xval1)
    end if
    numer = 0.0_rp
    denom = 0.0_rp
    do ipoin = 1,npoin_own
       if(lpoty(ipoin)==0) then
          xx = 0.0_rp
          do idime = 1,ndime
             xx = xx + real(idime,rp) * coord(idime,ipoin) 
          end do
          numer = numer + ( xval1(ipoin) - xx )**2
          denom = denom + ( xx )**2
       end if
    end do
    call PAR_SUM(numer)
    call PAR_SUM(denom)
    numer = sqrt(numer/denom)
    if( IMASTER ) print*,'Error 2=',numer
    !
    ! Testing projec_elements_to_nodes_1
    ! Solution = x+2y+3z
    !
    if( INOTMASTER ) then
       allocate(unkn2(nelem))
       do ielem = 1,nelem
          pelty = ltype(ielem)
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)
          xx    = 0.0_rp
          do inode = 1,pnode
             xx = xx + coord(1,lnods(inode,ielem))
          end do
          unkn2(ielem) = xx / real(pnode,rp)
       end do
       call projec_elements_to_nodes(unkn2,xval1)
    end if
    numer = 0.0_rp
    denom = 0.0_rp
    do ipoin = 1,npoin_own
       if(lpoty(ipoin)==0) then
          xx = coord(1,ipoin)
          numer = numer + ( xval1(ipoin) - xx )**2
          denom = denom + ( xx )**2
       end if
    end do
    call PAR_SUM(numer)
    call PAR_SUM(denom)
    numer = sqrt(numer/denom)
    if( IMASTER ) print*,'Error 3=',numer

    if( INOTMASTER ) then
       deallocate(unkno,xvalu)
       deallocate(unkn1,xval1)
       deallocate(unkn2)
    end if
    
    !call postpr_right_now('LPART','VECTO','NPOIN',xvalu)
    call runend('O.K.!')
    
  end subroutine projec_elements_to_nodes_unity_test
  
end module mod_projec
!> @}
