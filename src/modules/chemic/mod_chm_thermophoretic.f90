!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_thermophoretic
   use def_kintyp,     only : ip,rp
#ifdef CANTERA
   use cantera,        only : setMassFractions, setState_HP, temperature
#endif

   implicit none


   !
   ! Thermophoretic diffusion
   !
   real(rp),  pointer                    :: &
        Vthermoph_nodes(:,:),               & ! Thermophoretic velocity at nodes
        TermVT(:,:)

   private

   !
   ! Public variables
   !
   public :: Vthermoph_nodes

   !
   ! Public subroutines
   !
   public :: chm_thermophoretic_memory
   public :: chm_thermophoretic_gather
   public :: chm_thermophoretic_elmpre
   public :: chm_thermophoretic_calc_nodal_terms

contains


   subroutine chm_thermophoretic_gather(pnode,neq,iequa_first,lnods,elVtherm)
      implicit none
      integer(ip),           intent(in)  :: pnode
      integer(ip),           intent(in)  :: neq
      integer(ip),           intent(in)  :: iequa_first
      integer(ip),           intent(in)  :: lnods(pnode)
      real(rp),              intent(out) :: elVtherm(pnode,neq)

      integer(ip)                        :: inode,ipoin,iequa_last

      !
      ! Initialization
      !
      iequa_last = iequa_first + neq - 1
      elVtherm(1:pnode,1:neq) = 0.0_rp
      do inode=1,pnode
         ipoin=lnods(inode)
         elVtherm(inode,1:neq) = TermVT(ipoin,iequa_first:iequa_last)
      enddo

   end subroutine chm_thermophoretic_gather

   subroutine chm_thermophoretic_elmpre(pnode,pgaus,neq,gpsha,elVtherm,gpdivVt)
      implicit none
      integer(ip),          intent(in)    :: pnode
      integer(ip),          intent(in)    :: pgaus
      integer(ip),          intent(in)    :: neq
      real(rp),             intent(in)    :: gpsha(pnode,pgaus)
      real(rp),             intent(in)    :: elVtherm(pnode,neq)
      real(rp),             intent(inout) :: gpdivVt(pgaus,neq)

      integer(ip)                         :: inode,igaus

      !
      ! Soot mass frctions
      !
      gpdivVt(1:pgaus,1:neq) = 0.0_rp
      do igaus = 1,pgaus
         do inode = 1,pnode
            gpdivVt(igaus,1:neq) = gpdivVt(igaus,1:neq) &
                                 + gpsha(inode,igaus) * elVtherm(inode,1:neq)
         end do
      end do
   end subroutine chm_thermophoretic_elmpre



   subroutine chm_thermophoretic_calc_nodal_terms(neq,iequa_first)
      use def_domain,     only : ndime
      use def_domain,     only : npoin
      use def_domain,     only : nelem,ltype,ngaus
      use def_chemic,     only : kfl_model_chm
      use def_master,     only : conce,INOTEMPTY
      use def_master,     only : tempe
      use def_master,     only : tempe_gp
      use mod_memory,     only : memory_alloca
      use mod_memory,     only : memory_deallo
      use def_master,     only : mem_modul,modul
#ifdef CANTERA
      use def_chemic,     only : gas_chm,nspec_chm
      use def_master,     only : therm,prthe
#endif
      use mod_gradie,     only : gradie
      use mod_ker_proper, only : ker_proper
      use def_kintyp,     only : r1p

      implicit none
      integer(ip),           intent(in)  :: neq
      integer(ip),           intent(in)  :: iequa_first


      integer(ip)                        :: ipoin,dummi,iclas,iequa_last
      integer(ip)                        :: ielem,pelty,pgaus
      real(rp)                           :: coeff
      real(rp), pointer                  :: auxvec(:,:),   &
                                            prope_mu(:),   &
                                            prope_tem(:)
      type(r1p),pointer                  :: aux_r1p(:)

      character(35), parameter :: vacal = 'chm_thermophoretic_calc_nodal_terms'

      external                           :: smooth
      external                           :: divvec

      iequa_last = iequa_first + neq - 1

      if (INOTEMPTY) then

         nullify ( auxvec )
         nullify ( prope_mu )
         nullify ( prope_tem )

         call memory_alloca(mem_modul(1:2,modul),'AUXVEC',  vacal,auxvec,  ndime,npoin)
         call memory_alloca(mem_modul(1:2,modul),'PROPE_MU',vacal,prope_mu,      npoin)
         if ( kfl_model_chm == 3 .or. (.not. associated(tempe))) then
            call memory_alloca(mem_modul(1:2,modul),'PROPE_TEM',vacal,prope_tem,npoin)
         endif

         call ker_proper('VISCO','NPOIN',dummi,dummi,prope_mu)

         if ( kfl_model_chm == 3) then
#ifdef CANTERA
            do ipoin=1,npoin
               !
               ! Set Yk,H,P at node
               !
               call setMassFractions(gas_chm,conce(ipoin,1:nspec_chm,1))
               call setState_HP(gas_chm,therm(ipoin,1),prthe(1))
               !
               ! Get temperature (might be slightly different from tempe(ipoin,1)
               !
               prope_tem(ipoin) = temperature(gas_chm)
            end do
#endif
         elseif(associated(tempe)) then
            prope_tem => tempe(1:npoin,1)
         elseif(associated(tempe_gp)) then
            nullify(aux_r1p)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P',vacal,aux_r1p,nelem)
            do ielem = 1,nelem
               pelty = ltype(ielem)
               pgaus = ngaus(pelty)
               call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A',vacal,aux_r1p(ielem)%a,pgaus)
               aux_r1p(ielem) % a = tempe_gp(ielem) % a(:,1,1)
            end do
            call smooth (aux_r1p, prope_tem)
            call memory_deallo(mem_modul(1:2,modul),'AUX_R1P',vacal,aux_r1p)
         endif

         !
         ! Temperature gradient
         !
         call gradie(prope_tem,auxvec)

         !
         ! Thermophoretic velocity
         !
         do ipoin=1,npoin
            coeff = (-0.75_rp / (1.0_rp + 3.14156_rp*0.9_rp*0.125_rp)) * prope_mu(ipoin) / prope_tem(ipoin)
            Vthermoph_nodes(1:ndime,ipoin) = coeff*auxvec(1:ndime,ipoin)
         end do

         !
         ! Divergence of mass fraction flux transported with the thermophoretic veocity
         !
         do iclas = iequa_first,iequa_last
            do ipoin=1,npoin
               TermVT(ipoin,iclas)   = 0.0_rp
               auxvec(1:ndime,ipoin) = Vthermoph_nodes(1:ndime,ipoin) * conce(ipoin,iclas,1)
            enddo
            call divvec(auxvec,TermVT(1:npoin,iclas))
         end do

         call memory_deallo(mem_modul(1:2,modul),'AUXVEC',  vacal,auxvec)
         call memory_deallo(mem_modul(1:2,modul),'PROPE_MU',vacal,prope_mu)
         if ( kfl_model_chm == 3 .or. (.not. associated(tempe))) then
            call memory_deallo(mem_modul(1:2,modul),'PROPE_TEM',vacal,prope_tem)
         endif

      end if

   end subroutine chm_thermophoretic_calc_nodal_terms



   subroutine chm_thermophoretic_memory(nvari)
      use mod_memory, only : memory_alloca
      use def_master, only : mem_modul,modul
      use def_domain, only : ndime, npoin
      implicit none
      integer(ip), intent(in) :: nvari

      !
      ! Thermophoretic velocity
      !
      call memory_alloca(mem_modul(1:2,modul),'VTHERMOPH_NODES','chm_thermophoretic_memory',Vthermoph_nodes,ndime,npoin)
      call memory_alloca(mem_modul(1:2,modul),'TERMVT',         'chm_thermophoretic_memory',TermVT,         npoin,nvari)
   end subroutine chm_thermophoretic_memory

end module mod_chm_thermophoretic

