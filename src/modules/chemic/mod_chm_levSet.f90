!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_levSet

   use def_master,              only : INOTMASTER, ITASK_ENDINN, ITASK_INNITE, solve, mem_modul, modul, dtinv, advec, unkno
   use def_domain,              only : ndime, npoin, npoin
   use mod_gradie,              only : gradie
   use mod_memory,              only : memory_alloca
   use def_kintyp,              only : ip, rp
   use mod_element_integration, only : element_shape_function_derivatives_jacobian

   implicit none

   real(rp),   pointer, save :: rhs_ls(:,:)     ! RHS
   real(rp),   pointer, save :: phi(:)   ! fraction of volume
   real(rp),   pointer, save :: phig(:)   ! fraction of volume
   real(rp),   pointer, save :: phi0(:)  ! fraction of volume 0 in pseudo time
   real(rp),   pointer, save :: r(:,:)   ! residual
   real(rp),   pointer, save :: rg(:,:)   ! residual
   real(rp),   pointer, save :: Mass(:)  ! lumped mass
   real(rp),   pointer, save :: gradPhi(:,:)
   real(rp),   pointer, save :: mass_inv(:)
   real(rp),   pointer, save :: epsil(:)   ! fraction of volume
   real(rp),   pointer, save :: tau(:)

   real (rp), save :: a(5)
   real (rp), save :: b(5,5)
   real (rp), save :: sigma(5)
   real (rp), save :: pseudo_dtinv

! Uncomment if needed of delete
!   type(elm)   :: elmtyp

   private

   public :: chm_element_operations_levSet
   public :: chm_gather_levSet
   public :: chm_pseudo_time_reinit_levSet

contains

   subroutine chm_multi_step_fs_eval_levSet(istep)

      integer(ip) , intent(in) :: istep
      integer(ip) :: ipoin

      external    :: divvec
      external    :: rhsmod

      call chm_element_operations_levSet()

      if( INOTMASTER ) then
         do ipoin = 1,npoin
            r(ipoin,istep) = rhs_ls(ipoin,1)
         end do
      end if

      call rhsmod(1_ip,r(:,istep))

      if(istep == 4_ip) then
         call chm_multi_step_fs_solution_sij_levSet(4_ip,sigma,1.0_rp)
      else
         call chm_multi_step_fs_solution_sij_levSet(istep,b(istep+1,:),a(istep+1))
      endif

   end subroutine chm_multi_step_fs_eval_levSet

   subroutine chm_multi_step_fs_solution_sij_levSet(istep,weight1,weight2)
      use def_chemic,      only : nclas_chm, kfl_fixno_chm

      integer(ip) , intent(in) :: istep
      real(rp) ,    intent(in) :: weight1(5)
      real(rp) ,    intent(in) :: weight2
      integer(ip)       :: ipoin,jstep,kpoin
      real(rp) :: aux2

      if( INOTMASTER ) then

         do ipoin = 1,npoin
            aux2 = 0.0_rp
            do jstep=1,istep
               aux2   = aux2 + r(ipoin,jstep)*weight1(jstep)
            end do
            !rhs_ls(ipoin,1) = rhs_ls(ipoin,2)  + (weight2*aux2)
            phi(ipoin) = phi0(ipoin)  + (weight2*aux2)/(dtinv*Mass(ipoin))
         end do

      end if

      !call matrix_jacobi_CSR(1_ip,npoin,solve(1)%ia,solve(1)%ja,mass_inv,rhs_ls(:,1),phi)


      kpoin = 0
      do ipoin = 1,npoin
         if( kfl_fixno_chm(3,ipoin) > 0 ) then
            kpoin = (ipoin - 1) * nclas_chm + 3
            phi(ipoin) = unkno(kpoin)
         end if
      end do

   end subroutine chm_multi_step_fs_solution_sij_levSet

! Uncomment if needed of delete
!   subroutine chm_multi_step_fs_eval_levSet_galerkin(istep)
!
!      integer(ip) , intent(in) :: istep
!      integer(ip) :: ipoin
!
!      call chm_element_operations_levSet_galerkin()
!
!      if( INOTMASTER ) then
!         do ipoin = 1,npoin
!            rg(ipoin,istep) = rhs_ls(ipoin,1)
!         end do
!      end if
!
!      call rhsmod(1_ip,rg(:,istep))
!
!      if(istep == 4_ip) then
!         call chm_multi_step_fs_solution_sij_levSet_galerkin(4_ip,sigma,1.0_rp)
!      else
!         call chm_multi_step_fs_solution_sij_levSet_galerkin(istep,b(istep+1,:),a(istep+1))
!      endif
!
!   end subroutine chm_multi_step_fs_eval_levSet_galerkin
!
!   subroutine chm_multi_step_fs_solution_sij_levSet_galerkin(istep,weight1,weight2)
!      use def_chemic,      only : nclas_chm, kfl_fixno_chm
!      integer(ip) , intent(in) :: istep
!      real(rp) ,    intent(in) :: weight1(5)
!      real(rp) ,    intent(in) :: weight2
!      integer(ip)       :: ipoin,jstep,kpoin
!      real(rp) :: aux2
!
!      if( INOTMASTER ) then
!
!         do ipoin = 1,npoin
!            aux2 = 0.0_rp
!            do jstep=1,istep
!               aux2   = aux2 + rg(ipoin,jstep)*weight1(jstep)
!            end do
!            phig(ipoin) = phi0(ipoin)  + (weight2*aux2)/(dtinv*Mass(ipoin))
!         end do
!
!      end if
!
!      kpoin = 0
!      do ipoin = 1,npoin
!         if( kfl_fixno_chm(3,ipoin) > 0 ) then
!            kpoin = (ipoin - 1) * nclas_chm + 3
!            phig(ipoin) = unkno(kpoin)
!         end if
!      end do
!
!   end subroutine chm_multi_step_fs_solution_sij_levSet_galerkin

   subroutine chm_levSet_galerkin()

      use def_chemic,  only : nclas_chm, kfl_fixno_chm

      integer(ip) :: ipoin, kpoin

      external    :: rhsmod

      call chm_element_operations_levSet_galerkin()

      do ipoin = 1,npoin
         rg(ipoin,1) = rhs_ls(ipoin,1)
      end do

      call rhsmod(1_ip,rg(:,1))

      do ipoin = 1,npoin
         !phig(ipoin) = phi0(ipoin) + (1.5_rp*rg(ipoin,1) - 0.5_rp*rg(ipoin,2))/ (Mass(ipoin)*dtinv)
         phig(ipoin) = phi0(ipoin) + (rg(ipoin,1))/ (Mass(ipoin)*dtinv)
      end do

      kpoin = 0
      do ipoin = 1,npoin
         if( kfl_fixno_chm(3,ipoin) > 0 ) then
            kpoin = (ipoin - 1) * nclas_chm + 3
            phig(ipoin) = unkno(kpoin)
         end if
         rg(ipoin,2) = rg(ipoin,1)
      end do

   end subroutine chm_levSet_galerkin

   subroutine chm_pseudo_time_reinit_levSet()

      use def_chemic,    only: nclas_chm, kfl_fixno_chm, comax_chm, comin_chm, rtpts_chm

      integer(ip) :: iiter
      integer(ip) :: ipoin, kpoin, iclas_phi
      real(rp) :: mod_grad
      real(rp) :: gradPhi00(ndime,npoin),phi00(npoin),grad(ndime,npoin),gradc(ndime,npoin)
      real(rp) :: compress(npoin)
      integer(ip) :: nn
      integer(ip) :: do_reinit

      external    :: chm_endite
      external    :: divvec
      external    :: rhsmod

      do_reinit = 0_ip
      iclas_phi = 3
      rtpts_chm =  0.0_rp
      comin_chm =  1.0e9_rp
      comax_chm = -1.0e9_rp

      nn = solve(1)%ndofn
      solve(1) % ndofn = 1

      call chm_reinit_allocate_levSet()

      kpoin = 0
      do ipoin = 1,npoin
         phi0(ipoin)   = phi(ipoin)
         kpoin = (ipoin - 1) * nclas_chm + iclas_phi
         phi(ipoin) = unkno(kpoin)
      end do

      call gradie(phi,gradPhi)
      call chm_levSet_galerkin()
      !call chm_levSet_tau()

      !call chm_multi_step_fs_eval_levSet_galerkin(2_ip)      call chm_multi_step_fs_eval_levSet(2_ip)
      call chm_multi_step_fs_eval_levSet(2_ip)
      call gradie(phi,gradPhi)
      call chm_endite(ITASK_INNITE)

      !call chm_multi_step_fs_eval_levSet_galerkin(3_ip)
      call chm_multi_step_fs_eval_levSet(3_ip)
      call gradie(phi,gradPhi)
      call chm_endite(ITASK_INNITE)

      !call chm_multi_step_fs_eval_levSet_galerkin(4_ip)
      call chm_multi_step_fs_eval_levSet(4_ip)
      call chm_endite(ITASK_INNITE)

      ! re init

      if (do_reinit == 1_ip) then
         do ipoin = 1,npoin
            phi(ipoin) = min(phi(ipoin), 1.0_rp)
            phi(ipoin) = max(phi(ipoin), 0.0_rp)
         end do
         r(1:npoin,1:4) = 0.0_rp
         do iiter = 1, 3
            do ipoin = 1,npoin
               phi00(ipoin) = phi(ipoin)
            end do

            call gradie(phi,gradPhi00)

            do ipoin = 1, npoin
               grad(1:ndime,ipoin) = gradPhi00(1:ndime,ipoin)
               mod_grad  = sqrt(dot_product(gradPhi00(1:ndime,ipoin),gradPhi00(1:ndime,ipoin)))
               if(mod_grad > 0.0_rp) then
                  gradPhi00(1:ndime,ipoin) = gradPhi00(1:ndime,ipoin) / mod_grad
               else
                  gradPhi00(1:ndime,ipoin) = 0.0_rp
               end if
            end do
            gradc(1:ndime,1:npoin) = 0.0_rp
            do ipoin=1,npoin
               gradc(1:ndime,ipoin) =grad(1:ndime,ipoin)*epsil(ipoin) - phi(ipoin)*(1.0_rp - phi(ipoin))*gradPhi00(1:ndime,ipoin)
            end do

            call divvec(gradc,compress)
            do ipoin = 1,npoin
               if ( (phi(ipoin)>0.10_rp) .and. (phi(ipoin)<0.90_rp) ) then
                  r(ipoin,1) =  compress(ipoin)
               endif
            end do

            call rhsmod(1_ip,r(:,1))

            do ipoin = 1,npoin
               !phi(ipoin) = phi00(ipoin) + r(ipoin,1)/pseudo_dtinv
               phi(ipoin) = phi00(ipoin) + ((3.0_rp/2.0_rp)*r(ipoin,1) - 0.5_rp*r(ipoin,2) )/pseudo_dtinv
               phi(ipoin) = min(phi(ipoin), 1.0_rp)
               phi(ipoin) = max(phi(ipoin), 0.0_rp)
               r(ipoin,2) = r(ipoin,1)
            end do

            kpoin = 0
            do ipoin = 1,npoin
               if( kfl_fixno_chm(3,ipoin) > 0 ) then
                  kpoin = (ipoin - 1) * nclas_chm + 3
                  phi(ipoin) = unkno(kpoin)
               end if
            end do
            call chm_endite(ITASK_INNITE)
         end do
      end if

      kpoin = 0
      do ipoin = 1,npoin
         kpoin = (ipoin - 1) * nclas_chm + iclas_phi
         unkno(kpoin) = phi(ipoin)
      end do

      call chm_endite(ITASK_ENDINN)

      solve(1)%ndofn = nn

   end subroutine chm_pseudo_time_reinit_levSet

   subroutine chm_reinit_allocate_levSet()

      use def_chemic,         only : nclas_chm
      use mod_communications, only : PAR_MIN
      use def_domain,         only : ntens, mgaus, elmar, nelem, ltype, nnode, ngaus, llapl, lorde, ltopo, lnods, coord, mnode,&
                                     vmass, hnatu

      integer(ip), save :: ipass = 0
      integer(ip)       :: idime
      integer(ip)       :: ielem, inode, ipoin
      integer(ip)       :: pelty,pnode,pgaus
      integer(ip)       :: plapl,porde,ptopo
      integer(ip)       :: kpoin
      real(rp)          :: eps,c_eps,n_eps
      real(rp)          :: delta
      real(rp)          :: hleng(3),tragl(9)
      real(rp)          :: gphes(ntens,mnode,mgaus)
      real(rp)          :: gpvol(mgaus)
      real(rp)          :: gpcar(ndime,mnode,mgaus)
      real(rp)          :: elcod(ndime,mnode)
      real(rp)          :: pseudo_time_step

      external          :: elmlen
      external          :: elmcar

      if( ipass == 0 ) then

         ipass = 1
         nullify(r)
         nullify(rhs_ls)
         nullify(phi)
         nullify(phi0)
         nullify(rg)
         nullify(phig)
         nullify(Mass)

         nullify(epsil)


         if( INOTMASTER ) then
            call memory_alloca(mem_modul(1:2,modul),'R'      ,'chm_reinit_allocate_levSet',r, npoin,5_ip)
            call memory_alloca(mem_modul(1:2,modul),'RG'     ,'chm_reinit_allocate_levSet',rg, npoin,5_ip)
            call memory_alloca(mem_modul(1:2,modul),'RHS'    ,'chm_reinit_allocate_levSet',rhs_ls, npoin,2_ip)
            call memory_alloca(mem_modul(1:2,modul),'PHI'    ,'chm_reinit_allocate_levSet',phi, npoin)
            call memory_alloca(mem_modul(1:2,modul),'PHIG'   ,'chm_reinit_allocate_levSet',phig, npoin)
            call memory_alloca(mem_modul(1:2,modul),'PHI0'   ,'chm_reinit_allocate_levSet',phi0, npoin)
            call memory_alloca(mem_modul(1:2,modul),'MASS'   ,'chm_reinit_allocate_levSet',Mass, npoin)
            call memory_alloca(mem_modul(1:2,modul),'GRADPHI'  ,'chm_reinit_allocate_levSet',gradPhi, ndime, npoin)

            call memory_alloca(mem_modul(1:2,modul),'EPSILON'   ,'chm_reinit_allocate_levSet',epsil, npoin)
            call memory_alloca(mem_modul(1:2,modul),'TAU'   ,'chm_reinit_allocate_levSet',tau, npoin)
            call memory_alloca(mem_modul(1:2,modul),'MASSINV'   ,'chm_reinit_allocate_levSet',mass_inv, solve(1)%nzmat)

         else
            call memory_alloca(mem_modul(1:2,modul),'R'      ,'chm_reinit_allocate_levSet',r, 1_ip,4_ip)
            call memory_alloca(mem_modul(1:2,modul),'RG'     ,'chm_reinit_allocate_levSet',rg, 1_ip,4_ip)
            call memory_alloca(mem_modul(1:2,modul),'RHS'    ,'chm_reinit_allocate_levSet',rhs_ls, 1_ip,2_ip)
            call memory_alloca(mem_modul(1:2,modul),'PHI'    ,'chm_reinit_allocate_levSet',phi, 1_ip)
            call memory_alloca(mem_modul(1:2,modul),'PHIG'   ,'chm_reinit_allocate_levSet',phig, 1_ip)
            call memory_alloca(mem_modul(1:2,modul),'PHI0'   ,'chm_reinit_allocate_levSet',phi0, 1_ip)
            call memory_alloca(mem_modul(1:2,modul),'MASS'   ,'chm_reinit_allocate_levSet',Mass, 1_ip)
            call memory_alloca(mem_modul(1:2,modul),'GRADPHI'  ,'chm_reinit_allocate_levSet',gradPhi, 1_ip, 1_ip)

            call memory_alloca(mem_modul(1:2,modul),'MASSINV'   ,'chm_reinit_allocate_levSet',mass_inv, 1_ip)
            call memory_alloca(mem_modul(1:2,modul),'EPSILON'   ,'chm_reinit_allocate_levSet',epsil, 1_ip)
            call memory_alloca(mem_modul(1:2,modul),'TAU'   ,'chm_reinit_allocate_levSet',tau, 1_ip)

         end if

         r(1:npoin,1:4) = 0.0_rp
         kpoin = 0
         do ipoin = 1,npoin
            kpoin = (ipoin - 1) * nclas_chm + 3
            phi(ipoin) = unkno(kpoin)
            phi0(ipoin)   = phi(ipoin)
         end do

         a = 0.0_rp
         b = 0.0_rp
         sigma = 0.0_rp

         a(3) = 0.5_rp
         a(4) = 1.0_rp

         b(3,2) = 1_rp
         b(4,2) = 1.0_rp/4.0_rp
         b(4,3) = 1.0_rp/4.0_rp

         sigma(2) = 1.0_rp/6.0_rp
         sigma(3) = 1.0_rp/6.0_rp
         sigma(4) = 2.0_rp/3.0_rp

         !mysolve = solve(1)

      ! if(INOTMASTER) then
      !   call element_integration_initialization(elmtyp)
      !   call element_integration(ltype(1),elmtyp,'OPEN GAUSS LEGENDRE')
      !   !call element_integration(ltype(1),elmtyp,'CLOSE GAUSS LEGENDRE')

      !   Mass = 0.0_rp
      !   do ielem = 1, nelem
      !      pelty = ltype(ielem)
      !      pnode = nnode(pelty)
      !      pgaus = ngaus(pelty) ! is closed
      !      plapl = llapl(pelty)
      !      porde = lorde(pelty)
      !      ptopo = ltopo(pelty)
      !      gpvol = 0.0_rp
      !      gpcar = 0.0_rp
      !      gphes = 0.0_rp

      !      call element_shape_function_derivatives_jacobian(&
      !         pnode,pgaus,plapl,elmtyp % weigp,elmtyp %shape,&
      !         elmtyp % deriv,elmtyp % heslo,&
      !         elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)

      !      do inode=1,pnode
      !         ipoin=lnods(inode,ielem)
      !         do idime=1,ndime
      !            elcod(idime,inode)  = coord(idime,ipoin)
      !         end do
      !      end do

      !      call elmcar(&
      !         pnode,pgaus,plapl,elmtyp%weigp,elmtyp%shape,elmtyp%deriv, &
      !         elmtyp%heslo,elcod,gpvol,gpcar,gphes,ielem)

      !      elmat = 0.0_rp
      !      do inode = 1,pnode
      !         ipoin = lnods(inode,ielem)
      !         do jnode=1,pnode
      !            do igaus=1,pgaus
      !               elmat(inode,jnode) = elmat(inode,jnode) + elmtyp%shape(inode,igaus)*elmtyp%shape(jnode,igaus)*gpvol(igaus)
      !            end do
      !            Mass(ipoin) = Mass(ipoin) + elmat(inode,jnode)
      !         end do
      !      end do
      !      call matrix_assemble_element_matrix_to_CSR(1_ip,1_ip,pnode,pnode,&
      !         ielem,lnods(1:pnode,ielem),elmat,solve(1)%ia,solve(1)%ja,mass_inv)
      !   end do
      ! end if

       ! call rhsmod(1_ip,Mass)

       if(INOTMASTER) Mass(1:npoin) = vmass(1:npoin)

       pseudo_time_step = 1.0e10_rp
        do ielem = 1,nelem

           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = llapl(pelty)
           porde = lorde(pelty)
           ptopo = ltopo(pelty)

           gpvol = 0.0_rp
           gpcar = 0.0_rp
           gphes = 0.0_rp


           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)  = coord(idime,ipoin)
              end do
           end do

           call elmlen(&
              ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
              hleng)

           if (ndime == 2 ) then
              delta   = ( hleng(1) * hleng(2) )** 0.5_rp
           else
              delta = ( hleng(1) * hleng(2) * hleng(3) )** 0.3333333_rp
           endif

           c_eps = 0.5_rp
           n_eps = 0.9_rp
           eps   = c_eps * delta **(n_eps)

           pseudo_time_step = min(pseudo_time_step, 0.01_rp*(delta**2)/eps )


           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              epsil(ipoin) = max(epsil(ipoin), eps)
           end do

        end do

        call PAR_MIN(pseudo_time_step, 'IN MY CODE')

        pseudo_dtinv = 1.0_rp/pseudo_time_step

      end if

   end subroutine chm_reinit_allocate_levSet

!  subroutine matrix_jacobi_CSR(n1,n2,ia,ja,an,b,x)
!
!     use mod_memory,  only : memory_alloca, memory_deallo
!     use mod_matrix,  only : matrix_diagonal_CSR
!
!     integer(ip),  intent(in)                    :: n1                !< Starting node
!     integer(ip),  intent(in)                    :: n2                !< Final node
!     integer(ip),  intent(in), pointer           :: ia(:)             !< Matrix graph ia
!     integer(ip),  intent(in), pointer           :: ja(:)             !< Matrix graph ja
!     real(rp),     intent(inout)                 :: an(*)             !< Matrix
!     real(rp),     intent(in)                    :: b(*)              !< rhs
!     real(rp),     intent(inout)                 :: x(*)              !< x
!     integer(ip)                                 :: ii,jj,iz,pp,iter
!     real(rp),    pointer                        :: invdiag(:),r(:)
!
!     nullify(r,invdiag)
!     call memory_alloca(mem_modul(1:2,modul),'INVDIAG','matrix_jacobi_CSR',invdiag, npoin)
!     call memory_alloca(mem_modul(1:2,modul),'R','matrix_jacobi_CSR',r, npoin)
!
!     call matrix_diagonal_CSR(npoin,1_ip,0_ip,ia,ja,an,invdiag)
!     call rhsmod(1_ip,invdiag)
!     invdiag(1:npoin) = 1.0_rp/invdiag(1:npoin)
!
!     do ii = n1,n2
!        do iz = ia(ii),ia(ii+1)-1
!           jj = ja(iz)
!           if(ii == jj) then
!              an(iz) = 0.0_rp
!           end if
!        end do
!     end do
!
!     do iter=1,5
!        call matrix_CSR_SpMV(1_ip,npoin,1_ip,1_ip,ia,ja,an,x,r)
!        x(n1:n2) = (b(n1:n2)-r(n1:n2))*invdiag(n1:n2)
!        call rhsmod(1_ip,x)
!     end do
!
!     call memory_deallo(mem_modul(1:2,modul),'INVDIAG', 'matrix_jacobi_CSR',invdiag)
!     call memory_deallo(mem_modul(1:2,modul),'R', 'matrix_jacobi_CSR',r)
!
!  end subroutine matrix_jacobi_CSR

   subroutine chm_gather_levSet(pnode,lnods,elcon,elcod,elvel,elgra)
     use def_chemic, only     :  nclas_chm
     use def_domain, only     :  ndime,coord

     implicit none

     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: lnods(pnode)
     real(rp),    intent(out) :: elcon(pnode,nclas_chm,*)
     real(rp),    intent(out) :: elcod(ndime,pnode)
     real(rp),    intent(out) :: elvel(ndime,pnode)
     real(rp),    intent(out) :: elgra(ndime,pnode)
     integer(ip)              :: inode,ipoin,iclas,idime

     !
     ! Concentration and coordinates
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        do iclas=1,nclas_chm
           elcon(inode,iclas,1) = phi(ipoin)
           elcon(inode,iclas,2) = phig(ipoin)
           elcon(inode,iclas,3) = phi0(ipoin)
        end do
        do idime=1,ndime
           elcod(idime,inode)  = coord(idime,ipoin)
           elvel(idime,inode)  = advec(idime,ipoin,1)
           elgra(idime,inode)  = gradPhi(idime,ipoin)
        end do
     end do


   end subroutine chm_gather_levSet

   subroutine chm_element_operations_levSet_galerkin()

     use def_domain, only      : ndime,mnode,mgaus,nelem,ltype,nnode,     &
                                 ngaus,llapl,lorde,ltopo,lnods,elmar,     &
                                 ntens
     use def_chemic, only      : nclas_chm,                               &
                                 ncomp_chm

     implicit none
     integer(ip)               :: ielem, igaus, inode, ipoin,jnode
     integer(ip)               :: pelty,pnode,pgaus
     integer(ip)               :: plapl,porde,ptopo

     real(rp)                  :: elcon(mnode,nclas_chm,3)
     real(rp)                  :: elcod(ndime,mnode)
     real(rp)                  :: elrb(mnode)
     real(rp)                  :: elvel(ndime,mnode)
     real(rp)                  :: elgra(ndime,mnode)
     real(rp)                  :: gpvol(mgaus)
     real(rp)                  :: gpcon(mgaus,nclas_chm,ncomp_chm)
     real(rp)                  :: gpcar(ndime,mnode,mgaus)
     real(rp)                  :: gphes(ntens,mnode,mgaus)
     real(rp)                  :: gpder(ndime,mnode,mgaus)
     real(rp)                  :: gpvel(ndime,mgaus)
     real(rp)                  :: gpphi(mgaus)
     real(rp)                  :: gprhs(mgaus)
     real(rp)                  :: gpsha(mnode,mgaus)
     real(rp)                  :: gppr(mgaus)
     real(rp)                  :: gpadv(mnode,mgaus)
     real(rp)                  :: elmat(mnode,mnode)

     external                  :: elmcar

     !
     ! Loop over elements
     !
     rhs_ls = 0.0_rp

     elements: do ielem = 1,nelem

        !
        ! Element dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        !
        ! Initialization variables
        !
        gpcon = 0.0_rp
        gpvol = 0.0_rp
        gpcar = 0.0_rp
        gphes = 0.0_rp
        gprhs = 0.0_rp

        !
        ! Gather all
        !
        call chm_gather_levSet(&
           pnode,lnods(1:pnode,ielem),elcon(1:pnode,:,:),elcod,elvel,elgra)

        !
        ! CHALE, HLENG and TRAGL
        !

        call element_shape_function_derivatives_jacobian(&
           pnode,pgaus,plapl,elmar(pelty)% weigp,elmar(pelty)%shape,&
           elmar(pelty)% deriv,elmar(pelty)% heslo,&
           elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
        call elmcar(&
           pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
           elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)

        gpvel = 0.0_rp
        gpphi = 0.0_rp
        gppr  = 0.0_rp
        do inode=1,pnode
           do igaus=1, pgaus
              gpvel(1:ndime,igaus) = gpvel(1:ndime,igaus) + elvel(1:ndime,inode)*gpsha(inode,igaus)
           end do
        end do

        gpadv = 0.0_rp
        do inode=1,pnode
           do igaus=1, pgaus
              gpadv(inode,igaus) = gpadv(inode,igaus) + dot_product(gpvel(1:ndime,igaus),gpcar(1:ndime,inode,igaus))
           end do
        end do

        elrb(1:pnode) = 0.0_rp
        elmat = 0.0_rp
        do igaus=1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode
                 elmat(inode,jnode) =  elmat(inode,jnode) + gpadv(jnode,igaus)*gpsha(inode,igaus)*gpvol(igaus)
              end do
           end do
        end do

        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do jnode = 1,pnode
              rhs_ls(ipoin,1) = rhs_ls(ipoin,1) - elmat(inode,jnode)*elcon(jnode,1,1)
           end do
        end do

     end do elements

  end subroutine chm_element_operations_levSet_galerkin

   subroutine chm_element_operations_levSet()

     use def_domain,         only : ndime,mnode,mgaus,nelem,ltype,nnode,     &
                                    ngaus,llapl,lorde,ltopo,lnods,elmar,     &
                                    hnatu,ntens
     use def_chemic,         only : ncomp_chm,nclas_chm,                      &
                                    diffu_chm
     use mod_ker_proper,     only : ker_proper
     use def_kermod,         only : turmu_ker
     use mod_communications, only : PAR_MIN, PAR_MAX

     implicit none
     integer(ip)               :: ielem, igaus, inode, ipoin,jnode
     integer(ip)               :: pelty,pnode,pgaus
     integer(ip)               :: plapl,porde,ptopo

     real(rp)                  :: elcon(mnode,nclas_chm,3)
     real(rp)                  :: elcod(ndime,mnode)
     real(rp)                  :: elrb(mnode)
     real(rp)                  :: elvel(ndime,mnode)
     real(rp)                  :: elgra(ndime,mnode)
     real(rp)                  :: gpvol(mgaus)
     real(rp)                  :: gpcon(mgaus,nclas_chm,ncomp_chm)
     real(rp)                  :: gpcar(ndime,mnode,mgaus)
     real(rp)                  :: gphes(ntens,mnode,mgaus)
     real(rp)                  :: gpder(ndime,mnode,mgaus)
     real(rp)                  :: gpvel(ndime,mgaus)
     real(rp)                  :: gpphi(mgaus)
     real(rp)                  :: gprhs(mgaus)
     real(rp)                  :: hleng(3),tragl(9)
     real(rp)                  :: gpsha(mnode,mgaus)
     real(rp)                  :: aux,aux_min,aux_max
     real(rp)                  :: Emax,Emin,elEres(mnode)
     real(rp)                  :: fact(mnode),fact2,gppr(mgaus),fact3
     real(rp)                  :: gpadv(mnode,mgaus)
     real(rp)          :: elmat(mnode,mnode),ellap(mnode,mnode), stab, ve, betae
     real(rp)           :: gpmuturb(mgaus)
     integer(ip)        :: dummi

     external                  :: elmcar
     external                  :: elmlen

     !
     ! Loop over elements
     !
     rhs_ls = 0.0_rp

     elements: do ielem = 1,nelem

        !
        ! Element dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        !
        ! Initialization variables
        !
        gpcon = 0.0_rp
        gpvol = 0.0_rp
        gpcar = 0.0_rp
        gphes = 0.0_rp
        gprhs = 0.0_rp

        !
        ! Gather all
        !
        call chm_gather_levSet(&
           pnode,lnods(1:pnode,ielem),elcon(1:pnode,:,:),elcod,elvel,elgra)

        !
        ! CHALE, HLENG and TRAGL
        !

        call element_shape_function_derivatives_jacobian(&
           pnode,pgaus,plapl,elmar(pelty)% weigp,elmar(pelty)%shape,&
           elmar(pelty)% deriv,elmar(pelty)% heslo,&
           elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
        call elmcar(&
           pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
           elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)
        call elmlen(&
           ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
           hleng)

        gpvel = 0.0_rp
        gpphi = 0.0_rp
        gppr  = 0.0_rp

        ! Galerkin
        do inode=1,pnode
           do igaus=1, pgaus
              gpvel(1:ndime,igaus) = gpvel(1:ndime,igaus) + elvel(1:ndime,inode)*gpsha(inode,igaus)
           end do
        end do

        gpadv = 0.0_rp
        do inode=1,pnode
           do igaus=1, pgaus
              gpadv(inode,igaus) = gpadv(inode,igaus) + dot_product(gpvel(1:ndime,igaus),gpcar(1:ndime,inode,igaus))
           end do
        end do

        elrb(1:pnode) = 0.0_rp
        elmat = 0.0_rp
        do igaus=1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode
                 elmat(inode,jnode) =  elmat(inode,jnode) + gpadv(jnode,igaus)*gpsha(inode,igaus)*gpvol(igaus)
              end do
           end do
        end do

        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do jnode = 1,pnode
              rhs_ls(ipoin,1) = rhs_ls(ipoin,1) - elmat(inode,jnode)*elcon(jnode,1,1)
           end do
        end do

        ! stab
        fact = 0
        elEres = 0.0_rp
        Emin = 1.0e10_rp
        Emax = 1.0e-10_rp
        do inode=1,pnode
           Emax = max(Emax, -log(abs(elcon(inode,1,1)*(1.0_rp-elcon(inode,1,1))) + 1.0e-14_rp))
           Emin = min(Emin, -log(abs(elcon(inode,1,1)*(1.0_rp-elcon(inode,1,1))) + 1.0e-14_rp))
        end do

        do igaus=1,pgaus
           do inode=1,pnode
              fact2 = elcon(inode,1,1)
              fact3  = (-fact2*(1.0_rp-fact2)*(1.0_rp-2.0_rp*fact2))&
                  /(1.0e-14_rp+abs(fact2*(1.0_rp-fact2))*(abs(fact2*(1.0_rp-fact2)+1.0e-14_rp)))
              fact3 = gpvol(igaus)*gpsha(inode,igaus)*fact3/ (Emax-Emin + 1.0e-14_rp)
              elEres(inode) = elEres(inode) + fact3*((elcon(inode,1,2)-elcon(inode,1,1))*dtinv&
                  + dot_product(elvel(1:ndime,inode),elgra(1:ndime,inode)))
           end do
        end do
        ve = 0.0_rp
        betae = 0.0_rp
        stab = 1.0e10_rp
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           ve = max(abs(elEres(inode)/Mass(ipoin)),ve)
          betae = max(dot_product(elvel(1:ndime,inode),elvel(1:ndime,inode)), betae)
        end do

        aux = 1e15_rp
        aux = min(hleng(1),aux)
        aux = min(hleng(2),aux)
        if (ndime == 3) then
           aux = min(hleng(3),aux)
        end if

        ellap = 0.0_rp
        stab = min(0.5_rp*betae*aux,3.0_rp*ve*aux**2)
        !stab = min(0.5_rp*betae*aux,1.5_rp*ve*aux**2)


        ! anti difusion for the interfase
        aux = 0.0_rp
        aux_min = 1.0e10_rp
        aux_max = 1.0e-10_rp
        do inode=1,pnode
            aux = aux + elcon(inode,1,1)
            aux_max = max(aux_max,elcon(inode,1,1))
            aux_min = min(aux_min,elcon(inode,1,1))
        end do
        aux = aux/(1.0_rp*real(pnode,rp))

        stab = stab*max((1.0_rp - max(aux*(1.0_rp-aux),0.0_rp)/(abs(aux_max-aux_min) + 1.0e-10_rp)), 0.0_rp)

        do igaus=1,pgaus
           do jnode=1,pnode
              do inode = 1,pnode
                 ellap(inode,jnode) = ellap(inode,jnode) + stab*dot_product(gpcar(1:ndime,inode,igaus), gpcar(1:ndime,jnode,igaus))&
                     *gpvol(igaus)
              end do
           end do
        end do

        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do jnode = 1,pnode
              rhs_ls(ipoin,1) = rhs_ls(ipoin,1) - ellap(inode,jnode)*elcon(jnode,1,1)
           end do
        end do

        ! SGS clousue

        if(turmu_ker % kfl_exist /= 0_ip) then
           ellap = 0.0_rp

           call ker_proper('TURBU','PGAUS',dummi,ielem,gpmuturb,pnode,pgaus,elmar(pelty)%shape,gpcar)

           do igaus=1,pgaus
              do jnode=1,pnode
                 do inode = 1,pnode
                    ellap(inode,jnode) = ellap(inode,jnode) + (gpmuturb(igaus)/diffu_chm(1,1))&
                        * dot_product(gpcar(1:ndime,inode,igaus), gpcar(1:ndime,jnode,igaus))*gpvol(igaus)
                 end do
              end do
           end do

           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do jnode = 1,pnode
                 rhs_ls(ipoin,1) = rhs_ls(ipoin,1) - ellap(inode,jnode)*elcon(jnode,1,1)
              end do
           end do

        end if

     end do elements

  end subroutine chm_element_operations_levSet

! Uncomment if needed of delete
!  subroutine chm_levSet_tau()
!
!     use def_kintyp, only      : ip,rp
!     use def_domain, only      : ndime,mnode,mgaus,nelem,ltype,nnode,     &
!                                 ngaus,llapl,lorde,ltopo,lnods,elmar,     &
!                                 hnatu,ntens
!     use def_chemic, only      : ncomp_chm,nclas_chm
!
!     use mod_communications, only : PAR_MIN, PAR_MAX
!
!     implicit none
!     integer(ip)               :: ielem, inode, ipoin
!     integer(ip)               :: pelty,pnode,pgaus
!     integer(ip)               :: plapl,porde,ptopo
!
!     real(rp)                  :: elcon(mnode,nclas_chm,3)
!     real(rp)                  :: elcod(ndime,mnode)
!     real(rp)                  :: elvel(ndime,mnode)
!     real(rp)                  :: elgra(ndime,mnode)
!     real(rp)                  :: gpvol(mgaus)
!     real(rp)                  :: gpcon(mgaus,nclas_chm,ncomp_chm)
!     real(rp)                  :: gpcar(ndime,mnode,mgaus)
!     real(rp)                  :: gphes(ntens,mnode,mgaus)
!     real(rp)                  :: gpder(ndime,mnode,mgaus)
!     real(rp)                  :: gprhs(mgaus)
!     real(rp)                  :: hleng(3),tragl(9)
!     real(rp)                  :: gpsha(mnode,mgaus)
!     real(rp)                  :: aux,aux2
!
!     !
!     ! Loop over elements
!     !
!     rhs_ls = 0.0_rp
!
!     elements: do ielem = 1,nelem
!
!        !
!        ! Element dimensions
!        !
!        pelty = ltype(ielem)
!        pnode = nnode(pelty)
!        pgaus = ngaus(pelty)
!        plapl = llapl(pelty)
!        porde = lorde(pelty)
!        ptopo = ltopo(pelty)
!
!        !
!        ! Initialization variables
!        !
!        gpcon = 0.0_rp
!        gpvol = 0.0_rp
!        gpcar = 0.0_rp
!        gphes = 0.0_rp
!        gprhs = 0.0_rp
!
!        !
!        ! Gather all
!        !
!        call chm_gather_levSet(&
!           pnode,lnods(1:pnode,ielem),elcon(1:pnode,:,:),elcod,elvel,elgra)
!
!        !
!        ! CHALE, HLENG and TRAGL
!        !
!
!        call element_shape_function_derivatives_jacobian(&
!           pnode,pgaus,plapl,elmar(pelty)% weigp,elmar(pelty)%shape,&
!           elmar(pelty)% deriv,elmtyp % heslo,&
!           elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
!        call elmcar(&
!           pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
!           elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)
!        call elmlen(&
!           ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
!           hleng)
!
!
!        aux = 1e15_rp
!        aux = min(hleng(1),aux)
!        aux = min(hleng(2),aux)
!        if (ndime == 3) then
!           aux = min(hleng(3),aux)
!        end if
!
!        do inode=1,pnode
!           aux2 = 1e-15_rp
!           aux2 = max(abs(elvel(1,inode)),aux2)
!           aux2 = max(abs(elvel(2,inode)),aux2)
!           if (ndime == 3) then
!              aux2 = max(abs(elvel(3,inode)),aux2)
!           end if
!        end do
!
!        ipoin = lnods(inode,ielem)
!        tau(ipoin) = 1.0_rp/(2.0_rp*aux2/aux + dtinv)
!
!     end do elements
!
!  end subroutine chm_levSet_tau

end module mod_chm_levSet
