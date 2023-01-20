!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_gradie
  !------------------------------------------------------------------------
  !****f* Mathru/grasca
  ! NAME 
  !    grasca
  ! DESCRIPTION
  !    This routine computes the gradients of a scalar
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master, only             :  party,pardi,parki,parr2,&
       &                              NPOIN_REAL_2DIM,INOTMASTER,ISLAVE
  use def_kermod, only             :  kfl_grpro
  use def_domain, only             :  ndime,npoin,nelem,nnode,mnode,ntens,&
       &                              lnods,ltype,coord,vmass,elmar,vmasc,&
       &                              lexis,ngaus,kfl_naxis,&
       &                              mgaus,lelch,memor_dom
  use mod_solver, only             :  solver_lumped_mass_system
  use mod_memory, only             :  memory_alloca
  use mod_memory, only             :  memory_deallo
  
  use def_elmtyp
  implicit none

  private
  integer(ip), private             :: knode
  integer(ip), private             :: ipoin,idime,inode,ielem,igaus,jnode
  integer(ip), private             :: pnode,pelty,pgaus,jdime,itens
  real(rp),    private             :: gpdet,gpvol
  real(rp),    private             :: xjaci(9),xjacm(9),xfact

  interface gradie
     module procedure grasca,grasc2,graten,gravec,gravec_RP_3
  end interface

  public :: projec_hessian
  public :: grasca
  public :: gradie
  public :: gravec
  public :: graten
  
contains

  subroutine grasca(unkno,grunk)

    !------------------------------------------------------------------------
    !
    ! Case of a scalar
    !
    !------------------------------------------------------------------------
    implicit none
    real(rp),    intent(in)          :: unkno(npoin)
    real(rp),    intent(out), target :: grunk(ndime,npoin)
    real(rp)                         :: elunk(mnode)
    real(rp)                         :: gpcar(ndime,mnode)
    real(rp)                         :: gpgru(ndime)
    real(rp)                         :: elcod(ndime,mnode)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          do idime = 1,ndime
             grunk(idime,ipoin) = 0.0_rp
          end do
       end do
       !
       ! Loop over elements
       !
       if( kfl_grpro == 0 ) then

          elements: do ielem = 1,nelem
             pelty = ltype(ielem)
             if( pelty > 0 ) then
                pnode = nnode(pelty)
                pgaus = ngaus(pelty)
                !
                ! Gather vectors
                !
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      elcod(idime,inode)=coord(idime,ipoin)
                   end do
                   elunk(inode)=unkno(ipoin)
                end do
                !
                ! Treat extension elements
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if
                !
                ! Loop over Gauss points (which are nodes)
                !
                gauss_points: do igaus = 1,pgaus
                   call elmder(&
                        pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                        elcod,gpcar,gpdet,xjacm,xjaci)
                   gpvol = elmar(pelty)%weigp(igaus)*gpdet
                   !
                   ! Compute gradients
                   !
                   do idime = 1,ndime
                      gpgru(idime) = 0.0_rp
                      do inode = 1,pnode
                         gpgru(idime) = gpgru(idime) + gpcar(idime,inode) * elunk(inode)
                      end do
                   end do

                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      if( kfl_naxis == 1 ) then
                         call runend('MOD_GRADIE: NOT CODED')
                         !if(elcod(1,inode)==0.0_rp) then
                         !   gpvol=gpvol*twopi*1.0e-12
                         !else
                         !   gpvol=gpvol*twopi*elcod(1,inode)
                         !end if
                      end if
                      !
                      ! Assemble gradients
                      !
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                      do idime = 1,ndime               
                         grunk(idime,ipoin) = grunk(idime,ipoin) + xfact * gpgru(idime)
                      end do
                   end do

                end do gauss_points

             end if

          end do elements

       else

          elements2: do ielem = 1,nelem
             pelty = ltype(ielem) 
             if( pelty > 0 ) then
                pnode = nnode(pelty)
                pgaus = ngaus(pelty)
                !
                ! Gather vectors
                !
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      elcod(idime,inode) = coord(idime,ipoin)
                   end do
                   elunk(inode) = unkno(ipoin)
                end do
                !
                ! Treat extension elements
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if
                
                if( pelty /= PYR05 ) then
                   !
                   ! Loop over Gauss points (which are nodes)
                   !
                   gauss_points2: do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty) % weigc(inode) * gpdet
                      if( kfl_naxis == 1 ) then
                         if( elcod(1,inode) == 0.0_rp ) then
                            gpvol = gpvol*twopi*1.0e-12_rp
                         else
                            gpvol = gpvol*twopi*elcod(1,inode)
                         end if
                      end if
                      !
                      ! Assemble Gradients
                      !
                      do idime = 1,ndime
                         gpgru(idime) = 0.0_rp
                         do jnode = 1,pnode
                            gpgru(idime) = gpgru(idime) + gpcar(idime,jnode) * elunk(jnode)
                         end do
                      end do
                      do idime=1,ndime               
                         grunk(idime,ipoin) = grunk(idime,ipoin) + gpvol * gpgru(idime)
                      end do
                      
                   end do gauss_points2

                else
                   !
                   ! Treat PYR05 apart
                   !
                   gauss_points3: do igaus = 1,pgaus
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty)%weigp(igaus)*gpdet
                      !
                      ! Compute gradients
                      !
                      do idime = 1,ndime
                         gpgru(idime) = 0.0_rp
                         do inode = 1,pnode
                            gpgru(idime) = gpgru(idime) + gpcar(idime,inode) * elunk(inode)
                         end do
                      end do
                      !
                      ! Assemble gradients
                      !
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         if( kfl_naxis == 1 ) then
                            call runend('MOD_GRADIE: NOT CODED')
                         end if
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                         do idime = 1,ndime               
                            grunk(idime,ipoin) = grunk(idime,ipoin) + xfact * gpgru(idime)
                         end do
                      end do

                   end do gauss_points3
                end if
             end if
          end do elements2
       end if
       !
       ! Parallelization
       !
       !
       ! Solve diagonal system
       !
       if( kfl_grpro == 0 ) then
          !call solver_lumped_mass_system(ndime,grunk)
          call rhsmod(ndime,grunk)
          do ipoin = 1,npoin
             xfact = 1.0_rp / vmass(ipoin)
             do idime = 1,ndime
                grunk(idime,ipoin) = xfact * grunk(idime,ipoin)
             end do
          end do
       else
          call rhsmod(ndime,grunk)
          do ipoin = 1,npoin
             xfact = 1.0_rp / vmasc(ipoin)
             do idime = 1,ndime
                grunk(idime,ipoin) = xfact * grunk(idime,ipoin)
             end do
          end do
       end if
    end if

  end subroutine grasca

  subroutine grasc2(unkno,grunk,diffu)

    !------------------------------------------------------------------------
    !
    ! Case of a scalar
    !
    !------------------------------------------------------------------------

    implicit none
    real(rp),   intent(in)           :: unkno(npoin),diffu(npoin)
    real(rp),   intent(out), target  :: grunk(ndime,npoin)
    real(rp)                         :: elunk(mnode)
    real(rp)                         :: gpgru(ndime)
    real(rp)                         :: gpcar(ndime,mnode)
    real(rp)                         :: elcod(ndime,mnode)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          do idime = 1,ndime
             grunk(idime,ipoin) = 0.0_rp
          end do
       end do
       !
       ! Loop over elements
       !
       if( kfl_grpro == 0 ) then

          elements: do ielem=1,nelem
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
                   elunk(inode) = unkno(ipoin)
                end do
                !
                ! Loop over Gauss points 
                !
                gauss_points: do igaus = 1,pgaus
                   call elmder(&
                        pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                        elcod,gpcar,gpdet,xjacm,xjaci)
                   gpvol = elmar(pelty)%weigp(igaus)*gpdet
                   !
                   ! Compute gradients
                   !
                   do idime = 1,ndime
                      gpgru(idime) = 0.0_rp
                      do inode = 1,pnode
                         gpgru(idime) = gpgru(idime) + gpcar(idime,inode) * elunk(inode)
                      end do
                   end do
                   !
                   ! Assemble gradients
                   !
                   do inode = 1,knode
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus) 
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime               
                         grunk(idime,ipoin) = grunk(idime,ipoin) + xfact * diffu(ipoin) * gpgru(idime) 
                      end do
                   end do

                end do gauss_points
             end if
          end do elements

       else

          elements2: do ielem = 1,nelem
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
                do inode=1,pnode
                   ipoin=lnods(inode,ielem)
                   elcod(1:ndime,inode)=coord(1:ndime,ipoin)
                   elunk(inode)=unkno(ipoin)
                end do
                !
                ! Loop over Gauss points (which are nodes)
                !
                if( pelty /= PYR05 ) then

                   gauss_points3: do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty)%weigc(inode)*gpdet
                      !
                      ! Gradients
                      !
                      do idime = 1,ndime
                         gpgru(idime) = 0.0_rp
                         do jnode = 1,pnode
                            gpgru(idime) = gpgru(idime) + gpcar(idime,jnode) * elunk(jnode)
                         end do
                      end do
                      xfact = gpvol * diffu(ipoin)
                      do idime = 1,ndime               
                         grunk(idime,ipoin)=grunk(idime,ipoin) + xfact * gpgru(idime)
                      end do
                   end do gauss_points3

                else

                   gauss_points2: do igaus = 1,pgaus
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty)%weigp(igaus)*gpdet
                      !
                      ! Compute gradients
                      !
                      do idime = 1,ndime
                         gpgru(idime) = 0.0_rp
                         do inode = 1,pnode
                            gpgru(idime) = gpgru(idime) + gpcar(idime,inode) * elunk(inode)
                         end do
                      end do
                      !
                      ! Assemble gradients
                      !
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                         do idime = 1,ndime               
                            grunk(idime,ipoin) = grunk(idime,ipoin) + xfact * diffu(ipoin) * gpgru(idime) 
                         end do
                      end do
                   end do gauss_points2

                end if
             end if
          end do elements2

       end if
       !
       ! Periodicity
       !
       call rhsmod(ndime,grunk) 
       !
       ! Solve diagonal system
       !
       if( kfl_grpro == 0 ) then
          do ipoin = 1,npoin
             xfact = 1.0_rp / vmass(ipoin)
             do idime = 1,ndime
                grunk(idime,ipoin) = xfact * grunk(idime,ipoin)
             end do
          end do
       else
          do ipoin = 1,npoin
             xfact = 1.0_rp / vmasc(ipoin)
             do idime = 1,ndime
                grunk(idime,ipoin) = xfact * grunk(idime,ipoin)
             end do
          end do
       end if
    end if

  end subroutine grasc2

  subroutine graten(unkno,grunk)

    !------------------------------------------------------------------------
    !
    ! Case of a tensor
    !
    !------------------------------------------------------------------------
    implicit none
    real(rp),   intent(in)          :: unkno(ndime,npoin)
    real(rp),   intent(out), target :: grunk(ntens,npoin)
    real(rp)                        :: elgra(ndime,ndime)
    real(rp)                        :: stres(ndime,ndime)
    real(rp)                        :: elunk(ndime,mnode)
    real(rp)                        :: gpcar(ndime,mnode)
    real(rp)                        :: elcod(ndime,mnode)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          do itens = 1,ntens
             grunk(itens,ipoin) = 0.0_rp
          end do
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
                      ipoin = lnods(inode,ielem)
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
                      ! OJO: OPTIMZE LOOP. PUT RESULT IN ELEMENTAL ARRAY ELGRU
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus) 
                         grunk(1,ipoin) = grunk(1,ipoin) + xfact * ( elgra(1,1) + elgra(1,1) )
                         grunk(2,ipoin) = grunk(2,ipoin) + xfact * ( elgra(2,2) + elgra(2,2) ) 
                         grunk(3,ipoin) = grunk(3,ipoin) + xfact * ( elgra(1,2) + elgra(2,1) ) 
                         grunk(4,ipoin) = grunk(4,ipoin) + xfact * ( elgra(3,3) + elgra(3,3) ) 
                         grunk(5,ipoin) = grunk(5,ipoin) + xfact * ( elgra(1,3) + elgra(3,1) ) 
                         grunk(6,ipoin) = grunk(6,ipoin) + xfact * ( elgra(2,3) + elgra(3,2) ) 
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
                      do idime=1,ndime
                         do jdime=1,ndime                 
                            elgra(idime,jdime) = &
                                 dot_product(gpcar(idime,1:pnode),elunk(jdime,1:pnode))
                         end do
                      end do
                      do idime=1,ndime
                         do jdime=idime,ndime 
                            stres(idime,jdime) = (elgra(idime,jdime)+elgra(jdime,idime))
                         end do
                      end do
                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                         grunk(1,ipoin) = grunk(1,ipoin) + xfact * stres(1,1)
                         grunk(2,ipoin) = grunk(2,ipoin) + xfact * stres(2,2)
                         grunk(3,ipoin) = grunk(3,ipoin) + xfact * stres(1,2)
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

                         grunk(1,ipoin) = grunk(1,ipoin) + gpvol * ( elgra(1,1) + elgra(1,1) )
                         grunk(2,ipoin) = grunk(2,ipoin) + gpvol * ( elgra(2,2) + elgra(2,2) ) 
                         grunk(3,ipoin) = grunk(3,ipoin) + gpvol * ( elgra(1,2) + elgra(2,1) ) 
                         grunk(4,ipoin) = grunk(4,ipoin) + gpvol * ( elgra(3,3) + elgra(3,3) ) 
                         grunk(5,ipoin) = grunk(5,ipoin) + gpvol * ( elgra(1,3) + elgra(3,1) ) 
                         grunk(6,ipoin) = grunk(6,ipoin) + gpvol * ( elgra(2,3) + elgra(3,2) ) 

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
                         ! OJO: OPTIMZE LOOP. PUT RESULT IN ELEMENTAL ARRAY ELGRU
                         do inode = 1,knode
                            ipoin = lnods(inode,ielem)
                            xfact = gpvol * elmar(pelty) % shape(inode,igaus) 
                            grunk(1,ipoin) = grunk(1,ipoin) + xfact * ( elgra(1,1) + elgra(1,1) )
                            grunk(2,ipoin) = grunk(2,ipoin) + xfact * ( elgra(2,2) + elgra(2,2) ) 
                            grunk(3,ipoin) = grunk(3,ipoin) + xfact * ( elgra(1,2) + elgra(2,1) ) 
                            grunk(4,ipoin) = grunk(4,ipoin) + xfact * ( elgra(3,3) + elgra(3,3) ) 
                            grunk(5,ipoin) = grunk(5,ipoin) + xfact * ( elgra(1,3) + elgra(3,1) ) 
                            grunk(6,ipoin) = grunk(6,ipoin) + xfact * ( elgra(2,3) + elgra(3,2) ) 
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
                      do idime = 1,ndime
                         do jdime = 1,ndime                 
                            elgra(idime,jdime) = &
                                 dot_product(gpcar(idime,1:pnode),elunk(jdime,1:pnode))
                         end do
                      end do
                      do idime = 1,ndime
                         do jdime = idime,ndime 
                            stres(idime,jdime) = (elgra(idime,jdime)+elgra(jdime,idime))
                         end do
                      end do
                      grunk(1,ipoin) = grunk(1,ipoin) + gpvol * stres(1,1)
                      grunk(2,ipoin) = grunk(2,ipoin) + gpvol * stres(2,2)
                      grunk(3,ipoin) = grunk(3,ipoin) + gpvol * stres(1,2)
                   end do

                end if
             end do

          end if

       end if
       !
       ! Periodicity
       !
       call rhsmod(ntens,grunk) 
       !
       ! Solve diagonal system
       !
       if( kfl_grpro == 0 ) then
          if( ndime==1 ) then
             do ipoin=1,npoin
                grunk(1,ipoin) = grunk(1,ipoin) / vmass(ipoin)
             end do
          else if( ndime==2 ) then
             do ipoin=1,npoin
                xfact          = 1.0_rp / vmass(ipoin)
                grunk(1,ipoin) = xfact * grunk(1,ipoin)
                grunk(2,ipoin) = xfact * grunk(2,ipoin)
                grunk(3,ipoin) = xfact * grunk(3,ipoin)
             end do
          else
             do ipoin=1,npoin
                xfact          = 1.0_rp / vmass(ipoin)
                grunk(1,ipoin) = xfact * grunk(1,ipoin)
                grunk(2,ipoin) = xfact * grunk(2,ipoin)
                grunk(3,ipoin) = xfact * grunk(3,ipoin)
                grunk(4,ipoin) = xfact * grunk(4,ipoin)
                grunk(5,ipoin) = xfact * grunk(5,ipoin)
                grunk(6,ipoin) = xfact * grunk(6,ipoin)
             end do
          end if
       else
          if( ndime==1 ) then
             do ipoin=1,npoin
                grunk(1,ipoin) = grunk(1,ipoin) / vmasc(ipoin)
             end do
          else if( ndime==2 ) then
             do ipoin=1,npoin
                xfact          = 1.0_rp / vmasc(ipoin)
                grunk(1,ipoin) = xfact * grunk(1,ipoin)
                grunk(2,ipoin) = xfact * grunk(2,ipoin)
                grunk(3,ipoin) = xfact * grunk(3,ipoin)
             end do
          else
             do ipoin=1,npoin
                xfact          = 1.0_rp / vmasc(ipoin)
                grunk(1,ipoin) = xfact * grunk(1,ipoin)
                grunk(2,ipoin) = xfact * grunk(2,ipoin)
                grunk(3,ipoin) = xfact * grunk(3,ipoin)
                grunk(4,ipoin) = xfact * grunk(4,ipoin)
                grunk(5,ipoin) = xfact * grunk(5,ipoin)
                grunk(6,ipoin) = xfact * grunk(6,ipoin)
             end do
          end if
       end if
    end if

  end subroutine graten

  subroutine gravec_RP_3(unkno,grunk)

    implicit none
    real(rp),   intent(in),  pointer :: unkno(:,:,:)
    real(rp),   intent(inout), pointer :: grunk(:,:)
    
    if( associated(unkno) .and. associated(grunk) ) then
       call gravec(unkno,grunk)
    end if
    
  end subroutine gravec_RP_3
  
  subroutine gravec(unkno,grunk)

    !------------------------------------------------------------------------
    !
    ! Case of a vector
    !
    !------------------------------------------------------------------------

    implicit none
    real(rp),   intent(in)          :: unkno(ndime,npoin)
    real(rp),   intent(out), target :: grunk(ndime,ndime,npoin)
    real(rp)                        :: elgra(ndime,ndime)
    real(rp)                        :: elunk(ndime,mnode)
    real(rp)                        :: gpcar(ndime,mnode)
    real(rp)                        :: elcod(ndime,mnode)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin=1,npoin
          do idime=1,ndime
             do jdime=1,ndime
                grunk(jdime,idime,ipoin)=0.0_rp
             end do
          end do
       end do
       !
       ! Loop over elements
       !
       if( kfl_grpro == 0 ) then
          !
          ! Open rule
          !
          elements: do ielem = 1,nelem
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
                gauss_points: do igaus = 1,pgaus
                   call elmder(&
                        pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                        elcod,gpcar,gpdet,xjacm,xjaci)
                   gpvol = elmar(pelty)%weigp(igaus)*gpdet
                   !
                   ! Vector derivatives
                   !
                   do idime=1,ndime
                      do jdime=1,ndime                 
                         elgra(idime,jdime)=&
                              dot_product(gpcar(idime,1:pnode),elunk(jdime,1:pnode))
                      end do
                   end do

                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            grunk(jdime,idime,ipoin) = grunk(jdime,idime,ipoin) &
                                 + xfact * elgra(jdime,idime)
                         end do
                      end do
                   end do

                end do gauss_points
             end if
          end do elements

       else
          !
          ! Close rule
          !
          elements2: do ielem=1,nelem
             pelty = ltype(ielem) 
             if( pelty > 0 ) then
                pnode=nnode(pelty)
                pgaus=ngaus(pelty)
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
                do inode=1,pnode
                   ipoin=lnods(inode,ielem)
                   elcod(1:ndime,inode)=coord(1:ndime,ipoin)
                   elunk(1:ndime,inode)=unkno(1:ndime,ipoin)
                end do
                !
                ! Loop over Gauss points (which are nodes)
                !
                if( pelty /= PYR05 ) then

                   gauss_points3: do inode = 1,knode
                      ipoin=lnods(inode,ielem)
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol=elmar(pelty)%weigc(inode)*gpdet
                      !
                      ! Vector derivatives
                      !
                      do idime=1,ndime
                         do jdime=1,ndime                 
                            elgra(idime,jdime)=&
                                 dot_product(gpcar(idime,1:pnode),elunk(jdime,1:pnode))
                         end do
                      end do
                      do idime=1,ndime
                         do jdime=1,ndime
                            grunk(jdime,idime,ipoin)=grunk(jdime,idime,ipoin)&
                                 +gpvol*elgra(jdime,idime)
                         end do
                      end do
                   end do gauss_points3

                else

                   gauss_points2: do igaus = 1,pgaus
                      call elmder(&
                           pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                           elcod,gpcar,gpdet,xjacm,xjaci)
                      gpvol = elmar(pelty)%weigp(igaus)*gpdet
                      !
                      ! Vector derivatives
                      !
                      do idime=1,ndime
                         do jdime=1,ndime                 
                            elgra(idime,jdime)=&
                                 dot_product(gpcar(idime,1:pnode),elunk(jdime,1:pnode))
                         end do
                      end do

                      do inode = 1,knode
                         ipoin = lnods(inode,ielem)
                         xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                         do idime = 1,ndime
                            do jdime = 1,ndime
                               grunk(jdime,idime,ipoin) = grunk(jdime,idime,ipoin) &
                                    + xfact * elgra(jdime,idime)
                            end do
                         end do
                      end do
                      
                   end do gauss_points2

                end if
             end if
          end do elements2
       end if
       !
       ! Parall service
       !
       call rhsmod(ndime*ndime,grunk) 
       !
       ! Solve diagonal system
       !
       if( kfl_grpro == 0 ) then
          do ipoin = 1,npoin
             xfact = 1.0_rp / vmass(ipoin)
             do idime = 1,ndime
                do jdime = 1,ndime
                   grunk(jdime,idime,ipoin) = grunk(jdime,idime,ipoin) * xfact
                end do
             end do
          end do
       else
          do ipoin = 1,npoin
             xfact = 1.0_rp / vmasc(ipoin)
             do idime = 1,ndime
                do jdime = 1,ndime
                   grunk(jdime,idime,ipoin) = grunk(jdime,idime,ipoin) * xfact
                end do
             end do
          end do
       end if

    end if

  end subroutine gravec

  subroutine smosca(unkno,smunk)

    !------------------------------------------------------------------------
    !
    ! Case of a scalar
    !
    !------------------------------------------------------------------------
    implicit none
    type(r1p),    intent(in)          :: unkno(nelem)
    real(rp),     intent(out), target :: smunk(npoin)
    real(rp)                          :: gpcar(ndime,mnode)
    real(rp)                          :: elcod(ndime,mnode)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          smunk(ipoin) = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 .and. associated(unkno(ielem)%a) ) then
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
                do idime = 1,ndime
                   elcod(idime,inode) = coord(idime,ipoin)
                end do
             end do
             !
             ! Loop over Gauss points 
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,gpdet,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * gpdet
                !
                ! Assemble
                !
                if( kfl_naxis == 1 ) then
                   call runend('MOD_GRADIE: NOT CODED')
                end if
                do inode = 1,knode
                   ipoin        = lnods(inode,ielem)
                   xfact        = gpvol * elmar(pelty) % shape(inode,igaus)
                   smunk(ipoin) = smunk(ipoin) + xfact * unkno(ielem) % a(igaus)
                end do

             end do gauss_points
          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(1_ip,smunk)
       !
       ! Solve diagonal system
       !
       do ipoin = 1,npoin
          smunk(ipoin) = smunk(ipoin) / vmass(ipoin)
       end do
    end if

  end subroutine smosca

    !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-07-05
  !> @brief   Hessian
  !> @details Compute the hessian of a scalar using a double projection
  !> 
  !-----------------------------------------------------------------------

  subroutine projec_hessian(unkno,hessi)
    
    real(rp),    intent(in),    pointer :: unkno(:)
    real(rp),    intent(inout), pointer :: hessi(:,:,:)
    real(rp),                   pointer :: gradi(:,:)
    integer(ip)                         :: idime,jdime,ipoin
    
    nullify(gradi)
    call memory_alloca(memor_dom,'GRADI','mod_projec',gradi,ndime,npoin)

    if( .not. associated(hessi) ) then
       call memory_alloca(memor_dom,'HESSI','mod_projec',hessi,ndime,ndime,npoin)
    end if

    call grasca(unkno,gradi)
    call gravec(gradi,hessi)
    
    call memory_deallo(memor_dom,'GRADI','mod_projec',gradi)
    !
    ! Symmetrize Hessian
    !
    do ipoin = 1,npoin
       do idime = 1,ndime
          do jdime = idime+1,ndime
             hessi(idime,jdime,ipoin) = 0.5_rp*(hessi(idime,jdime,ipoin)+hessi(jdime,idime,ipoin))
             hessi(jdime,idime,ipoin) = hessi(idime,jdime,ipoin)
          end do
       end do
    end do
    
  end subroutine projec_hessian
  
end module mod_gradie
  
