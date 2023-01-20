!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    mod_tem_fast.f90
!> @author  goyarzun
!> @date    2020-08-01
!> @brief   Fast version of temper
!> @details Vectorization reordering applied to temper

module mod_tem_fast

#include "def_vector_size.inc"

#define DEF_VECT 1:VECTOR_SIZE

    implicit none

    private

    public :: tem_chemic_fast
    public :: tem_elmchl_fast
    public :: tem_elmdir_fast
    public :: tem_elmgat_fast
    public :: tem_rhocpdt_fast
    public :: tem_turbul_fast
    public :: tem_velfun_fast
    public :: tem_elmpre_fast
    public :: tem_elmlen_fast
    public :: tem_element_assembly
    public :: mod_tem_element_operations_fast
contains
    subroutine tem_chemic_fast(&
         VECTOR_DIM,pgaus,gprhs,list_elements_p)
      !-----------------------------------------------------------------------
      !****f* Temper/tem_chemic
      ! NAME
      !   tem_radiat
      ! DESCRIPTION
      !    Couple to the heat source from chemical processes
      ! USES
      ! USED BY
      !    tem_elmop2 
      !***
      !-----------------------------------------------------------------------
      use def_kintyp, only       :  ip,rp
      use def_master, only       :  div_enthalpy_transport,kfl_coupl, &
                                    ID_TEMPER,ID_CHEMIC
      implicit none 
      integer(ip), intent(in)    :: VECTOR_DIM,pgaus
      integer(ip), intent(in)    :: list_elements_p(VECTOR_DIM)
      real(rp),    intent(inout) :: gprhs(VECTOR_DIM,pgaus)
      integer(ip)                :: ivect,ielem
    
      !
      ! Coupling with CHEMIC
      !
      if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 .and. associated(div_enthalpy_transport)) then
         !
         ! RHS terms: Divergence of enthalpy transport by diffusion
         !
         do ivect=1,VECTOR_DIM
            ielem=list_elements_p(ivect)
             gprhs(ivect,:) = gprhs(ivect,:) + div_enthalpy_transport(ielem)%a(:,1,1)
         end do
    
      endif
    
    end subroutine tem_chemic_fast

    subroutine tem_elmchl_fast(VECTOR_DIM, &
     tragl,hleng,elcod,elvel,chave,chale,pnode,porde,&
     hnatu,kfl_advec,kfl_ellen)
      !-----------------------------------------------------------------------
      !****f* Domain/elmchl
      ! NAME
      !   elmchl
      ! DESCRIPTION
      !   This routine computes the characteristic element lengths CHALE 
      !   according to a given strategy. CHALE is divided by two for
      !   quadratic elements:
      !   KFL_ELLEN = 0 ... CHALE(1) = Minimum element length
      !                 ... CHALE(2) = Minimum element length
      !   KFL_ELLEN = 1 ... CHALE(1) = Maximum element length
      !                 ... CHALE(2) = Maximum element length
      !   KFL_ELLEN = 2 ... CHALE(1) = Average element length
      !                 ... CHALE(2) = Average element length
      !   KFL_ELLEN = 3 ... IF KFL_ADVEC = 1:
      !                     CHALE(1) = Flow direction
      !                     CHALE(2) = Flow direction
      !                     ELSE IF KFL_ADVEC =0:
      !                     CHALE(1) = Minimum element length
      !                     CHALE(2) = Minimum element length
      !   KFL_ELLEN = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
      !                 ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
      !   KFL_ELLEN = 5 ... CHALE(1) = Length in flow direction
      !                 ... CHALE(2) = Minimum element kength
      ! OUTPUT
      !   CHALE
      ! USES
      ! USED BY
      !***
      !-----------------------------------------------------------------------
      use def_kintyp, only     :  ip,rp
      use def_domain, only     :  ndime
      implicit none
      integer(ip), intent(in)  :: VECTOR_DIM,pnode,porde,kfl_advec,kfl_ellen
      real(rp),    intent(in)  :: hnatu
      real(rp),    intent(out) :: chale(VECTOR_DIM,2)
      real(rp),    intent(in)  :: tragl(VECTOR_DIM,ndime,ndime),hleng(VECTOR_DIM,ndime)
      real(rp),    intent(in)  :: elcod(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(in)  :: elvel(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(out) :: chave(VECTOR_DIM,ndime,2)
      integer(ip)              :: idime,inode,ivect
      real(rp)                 :: elno1,elno2
  
  
    if(kfl_ellen==0) then 
       !
       ! Minimum element length
       !
       chale(DEF_VECT,1)=hleng(DEF_VECT,ndime) 
       chale(DEF_VECT,2)=chale(DEF_VECT,1)
  
    else if(kfl_ellen==1) then   
       !
       ! Maximum element length
       !     
       chale(DEF_VECT,1)=hleng(DEF_VECT,1) 
       chale(DEF_VECT,2)=chale(DEF_VECT,1)
  
    else if(kfl_ellen==2) then 
       !
       ! Average length
       !
       chale(DEF_VECT,1)=0.0_rp
       do idime=1,ndime
          chale(DEF_VECT,1)=chale(DEF_VECT,1)+hleng(DEF_VECT,idime)
       end do
       chale(DEF_VECT,1)=chale(DEF_VECT,1)/real(ndime,rp) 
       chale(DEF_VECT,2)=chale(DEF_VECT,1)
  
    else if(kfl_ellen==3) then 
       !
       ! Length in flow direction
       !
       if(kfl_advec/=0) then 
          !
          ! Characteristic element velocity (average)
          !
          chave=0.0_rp
          do idime=1,ndime
             do inode=1,pnode
                chave(DEF_VECT,idime,1)=chave(DEF_VECT,idime,1)+elvel(DEF_VECT,idime,inode)
             end do
             chave(DEF_VECT,idime,1)=chave(DEF_VECT,idime,1)/real(pnode,rp)
          end do
          !
          ! Characteristic element length u^l = J^(-t) u^g
          !
          do ivect=1,VECTOR_DIM
            call mbvab1(chave(ivect,1,2),tragl(ivect,:,:),chave(ivect,1,1),ndime,ndime,elno2,elno1)
            if(elno2>1.0e-16_rp.and.elno1>1.0e-16_rp) then
               chale(ivect,1)=hnatu*elno1/elno2
            else
               chale(ivect,1)=hleng(ivect,ndime)
            end if
          end do
  
          chale(DEF_VECT,2)=chale(DEF_VECT,1)
          chale(DEF_VECT,2)=hleng(DEF_VECT,ndime)
  
          if (ndime ==3 ) then
             chale(DEF_VECT,2)=(hleng(DEF_VECT,ndime)*hleng(DEF_VECT,2)*hleng(DEF_VECT,1))**(1.0_rp/3.0_rp)
          else if (ndime==2) then
             chale(DEF_VECT,2)=sqrt(hleng(DEF_VECT,2)*hleng(DEF_VECT,1))
          end if
       else
          chale(DEF_VECT,1)=hleng(DEF_VECT,ndime)       
          chale(DEF_VECT,2)=chale(DEF_VECT,1)
       end if
  
    else if(kfl_ellen==4) then 
       !
       ! sqrt(hmin*hmax)
       !
       chale(DEF_VECT,1)=sqrt(hleng(DEF_VECT,1)*hleng(DEF_VECT,ndime))     
       chale(DEF_VECT,2)=chale(DEF_VECT,1)
  
    else if(kfl_ellen==5) then 
       !
       ! Along velocity direction
       !
       do ivect=1,VECTOR_DIM
         call velchl(pnode,elcod(ivect,:,:),elvel(ivect,:,:),chale(ivect,:),hleng(ivect,:))
       end do
  
    else if(kfl_ellen==6) then 
       !
       ! Mixed element length - hmin for tau1, hmax for tau2 - here we only obtain the values for tau1 - tau2 directly in nsi_elmsgs
       !
       chale(DEF_VECT,1)=hleng(DEF_VECT,ndime) 
       chale(DEF_VECT,2)=chale(DEF_VECT,1)
  
    end if
    !
    ! Divide h by 2 for quadratic elements and 3 for cubic elements
    !
    chale(DEF_VECT,1) = chale(DEF_VECT,1)/real(porde,rp)
    chale(DEF_VECT,2) = chale(DEF_VECT,2)/real(porde,rp)
  
  end subroutine tem_elmchl_fast

    subroutine tem_elmdir_fast(VECTOR_DIM,&
         pnode,lnods_loc,elmat,elrhs,list_elements_p) !ielem)
      !------------------------------------------------------------------------
      !****f* Temper/tem_elmdir
      ! NAME 
      !    tem_elmdir
      ! DESCRIPTION
      !    This routine prescribes the boundary conditions for the 
      !    temperature equations. 
      ! USES
      ! USED BY
      !    tem_elmope
      !    tem_bouope
      !------------------------------------------------------------------------
      use def_kintyp, only       :  ip,rp
      use def_temper, only       :  bvess_tem,kfl_fixno_tem
      use def_domain, only       :  nhang,lhang,lelch
      use def_elmtyp, only       :  ELHAN
      implicit none
      integer(ip), intent(in)    :: VECTOR_DIM,pnode
      integer(ip), intent(in)    :: lnods_loc(VECTOR_DIM,pnode)
      real(rp),    intent(inout) :: elmat(VECTOR_DIM,pnode,pnode),elrhs(VECTOR_DIM,pnode)
      integer(ip), intent(in)    :: list_elements_p(VECTOR_DIM)
      integer(ip)                :: inode,ipoin,jnode,jpoin,ihang,ivect,ielem
!      integer(ip)                :: knode
      real(rp)                   :: adiag,xvalu
      
      do ivect = 1, VECTOR_DIM
         ielem=list_elements_p(ivect)
         if( lelch(ielem) == ELHAN ) then
            do inode = 1,pnode
               ipoin = lnods_loc(ivect,inode)
               do ihang = 1,nhang
                  if( lhang(ihang) % l(1) == ipoin ) then
                     do jnode = 1,pnode
                        jpoin = lnods_loc(ivect,jnode)
                        if( lhang(ihang) % l(2) == jpoin .or. lhang(ihang) % l(3) == jpoin ) then
                           elrhs(ivect,jnode) = elrhs(ivect,jnode) + elrhs(ivect,inode)
                           elrhs(ivect,inode) = 0.0_rp
    !                       do knode = 1,pnode
                              elmat(ivect,jnode,:) = elmat(ivect,jnode,:) + elmat(ivect,inode,:)
                              elmat(ivect,inode,:) = 0.0_rp
    !                       end do
       
    !                       do knode = 1,pnode
                              elmat(ivect,:,jnode) = elmat(ivect,:,jnode) + elmat(ivect,:,inode)
                              elmat(ivect,:,inode) = 0.0_rp
    !                       end do
                           elmat(ivect,inode,jnode) = -1.0_rp
                           elmat(ivect,inode,inode) =  1.0_rp
                        end if
                     end do
                  end if
               end do
            end do
         end if
      end do
      !if(ielem==2) then
      !   inode = 1
      !   ipoin = lnods(inode)
      !   do while( ipoin /= 6 )
      !      inode = inode + 1
      !      ipoin = lnods(inode)
      !   end do
      !   jnode = 1
      !   jpoin = lnods(jnode)
      !   do while( jpoin /= 3 )
      !      jnode = jnode + 1
      !      jpoin = lnods(jnode)
      !   end do
      !end if
      !if(ielem==3) then
      !   inode = 1
      !   ipoin = lnods(inode)
      !   do while( ipoin /= 6 )
      !      inode = inode + 1
      !      ipoin = lnods(inode)
      !   end do
      !   jnode = 1
      !   jpoin = lnods(jnode)
      !  do while( jpoin /= 9 )
      !      jnode = jnode + 1
      !      jpoin = lnods(jnode)
      !   end do
      !end if
      !if( ielem == 2 .or. ielem == 3 ) then
      !   elrhs(jnode) = elrhs(jnode) + elrhs(inode)
      !   elrhs(inode) = 0.0_rp
      !   do knode = 1,pnode
      !      elmat(jnode,knode) = elmat(jnode,knode) + elmat(inode,knode)
      !      elmat(inode,knode) = 0.0_rp
      !   end do
      !   do knode = 1,pnode
      !      elmat(knode,jnode) = elmat(knode,jnode) + elmat(knode,inode)
      !      elmat(knode,inode) = 0.0_rp
      !   end do
      !   elmat(inode,jnode) = -1.0_rp
      !   elmat(inode,inode) =  1.0_rp
      !end if
      do ivect = 1, VECTOR_DIM
         do inode = 1,pnode
            ipoin = lnods_loc(ivect,inode)
            if(  kfl_fixno_tem(1,ipoin) == 1 .or.&
                 kfl_fixno_tem(1,ipoin) == 4 .or.&
                 kfl_fixno_tem(1,ipoin) == 5 .or.&
                 kfl_fixno_tem(1,ipoin) == 7 ) then
               adiag = elmat(ivect,inode,inode)
               xvalu = bvess_tem(1,ipoin,1)
    !           do jnode = 1,pnode
                  elmat(ivect,inode,:) = 0.0_rp
                  elrhs(ivect,:)       = elrhs(ivect,:) - elmat(ivect,:,inode) * xvalu
                  elmat(ivect,:,inode) = 0.0_rp
     !          end do
               if( abs(adiag) > 0.0_rp ) then
                  elmat(ivect,inode,inode) = adiag
                  elrhs(ivect,inode)       = adiag * xvalu
               else
                  elmat(ivect,inode,inode) = 1.0_rp
                  elrhs(ivect,inode)       = xvalu           
               end if
            end if
         end do
      end do
    end subroutine tem_elmdir_fast

    subroutine tem_elmgat_fast(&
         VECTOR_DIM,pnode,lnods,eltem,elvel,elcod,eledd,elmsh)
      !------------------------------------------------------------------------
      !****f* Temper/tem_elmgat
      ! NAME 
      !    tem_elmgat
      ! DESCRIPTION
      !    This routine performs the gather operations
      ! USES
      ! USED BY
      !    tem_elmope
      !    tem_bouope
      !***
      !------------------------------------------------------------------------ 
      use def_kintyp, only     :  ip,rp 
      use def_domain, only     :  ndime,coord
      use def_master, only     :  therm,advec,veloc_forw,velom
!      use def_master, only     :  turmu
      use def_temper, only     :  ADR_tem,kfl_advec_tem
      use mod_ADR,    only     :  BDF
      use def_kermod, only     :  kfl_adj_prob
      use mod_ker_proper
      implicit none
      integer(ip), intent(in)  :: pnode 
      integer(ip), intent(in)  :: VECTOR_DIM
      integer(ip), intent(in)  :: lnods(VECTOR_DIM,pnode)
      real(rp),    intent(out) :: eltem(VECTOR_DIM,pnode,*)
      real(rp),    intent(out) :: elcod(VECTOR_DIM,ndime,pnode),elvel(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(out) :: eledd(VECTOR_DIM,pnode)
      real(rp),    intent(out) :: elmsh(VECTOR_DIM,ndime,pnode)
      integer(ip)              :: inode,ipoin,itime,ivect
      !
      ! Current temperature and coordinates
      !
      do ivect=1,VECTOR_DIM
          do inode = 1,pnode
             ipoin = lnods(ivect,inode)
             eltem(ivect,inode,1) = therm(ipoin,1)
             elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
          end do
      end do
      !
      ! Turbulent viscosity mut from RANS   !Podria comentarse
      !
    !  if( kfl_cotur_tem > 0 ) then
    !    do ivect =1,VECTOR_DIM
    !       do inode = 1,pnode
    !          ipoin = lnods(ivect,inode)
    !          eledd(ivect,inode) = turmu(ipoin)
    !       end do
    !    end do
    !  end if
      !
      ! Time integration
      !
      if( ADR_tem % kfl_time_integration /= 0 ) then
         do ivect=1,VECTOR_DIM
            do inode = 1,pnode
               ipoin = lnods(ivect,inode)
               eltem(ivect,inode,2) = therm(ipoin,3)
            end do
         end do
    
         if( ADR_tem % kfl_time_scheme == BDF ) then
            do ivect = 1,VECTOR_DIM
              do itime = 3,ADR_tem % kfl_time_order + 1
                 do inode = 1,pnode 
                    ipoin = lnods(ivect,inode)
                    eltem(ivect,inode,itime) = therm(ipoin,itime+1)
                 end do
              end do
            end do
         end if
      end if
      !
      ! ALE terms
      !
      if( associated(velom) ) then
         do ivect = 1,VECTOR_DIM
           do inode = 1,pnode
              ipoin = lnods(ivect,inode)
              elmsh(ivect,1:ndime,inode)=velom(1:ndime,ipoin)
           end do
         end do
      else
         elmsh = 0.0_rp 
      end if
      !
      ! Coupling with flow equations
      !
      if( kfl_advec_tem == 1 ) then
         do ivect = 1,VECTOR_DIM
           do inode = 1,pnode
              ipoin = lnods(ivect,inode)
              if (kfl_adj_prob == 0_ip) then
                 elvel(ivect,1:ndime,inode)=advec(1:ndime,ipoin,1)
              else
                 elvel(ivect,1:ndime,inode)=veloc_forw(1:ndime,ipoin,1)
              endif
           end do
         end do
    
      else if( kfl_advec_tem >= 2 ) then
         call tem_velfun_fast(VECTOR_DIM,pnode,elcod,elvel)
    
      end if
    
    end subroutine tem_elmgat_fast

    subroutine tem_rhocpdt_fast(VECTOR_DIM,pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden,list_elements_p,dt_rho_cp_tem_loc)
      use def_parame
      use def_elmtyp
      use def_master
      use def_domain
      use def_kermod
      use mod_ker_proper 
      use def_temper
      implicit none
    
      integer(ip), intent(in)  :: VECTOR_DIM
      integer(ip), intent(in)  :: pnode
      integer(ip), intent(in)  :: pgaus
      integer(ip), intent(in)  :: porde
      integer(ip), intent(in)  :: lnods_loc(VECTOR_DIM,pnode)
      real(rp),    intent(in)  :: gpsha(VECTOR_DIM,pnode,pgaus)
      real(rp),    intent(in)  :: gpvol(VECTOR_DIM,pgaus)
      real(rp),    intent(in)  :: gpden(VECTOR_DIM,pgaus)
      integer(ip), intent(in)  :: list_elements_p(VECTOR_DIM)
      real(rp),    intent(out) :: dt_rho_cp_tem_loc(*)
      integer(ip) :: inode,jnode,ipoin,igaus,ielemone,ivect,ielem
      real(rp) :: eldtcprho(VECTOR_DIM,pnode), fact(VECTOR_DIM)
!      real(rp) :: d,T
      real(rp)                :: elmat(VECTOR_DIM,pnode,pnode)
      real(rp)                :: trace(VECTOR_DIM), elmass(VECTOR_DIM)       
      ielemone= list_elements_p(1)
    
      if (kfl_rhs_scal_tem == 0 ) then
    
         if( porde == 1 ) then
            !
            ! Element assembly
            !
            eldtcprho = 0.0_rp
    
            do igaus = 1,pgaus
               fact(DEF_VECT) = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  * dtinv )
               do inode = 1,pnode
                  eldtcprho(DEF_VECT,inode) = eldtcprho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * fact(DEF_VECT)
               end do
            end do
    
            !
            ! Nodal projection
            !
            do ivect = 1, VECTOR_DIM
               ielem = list_elements_p(ivect)
               if(ivect == 1 .or. (ivect > 1 .and. ielem /= ielemone)) then
                 do inode = 1,pnode
                    ipoin = lnods_loc(ivect,inode)
                    dt_rho_cp_tem_loc(ipoin) = dt_rho_cp_tem_loc(ipoin) + eldtcprho(ivect,inode) 
                 end do
               end if
            end do
    
         else
            !
            ! Element assembly
            !
    !@        eldtcprho = 0.0_rp
    
            !  T = 0.0_rp
            !  d = 0.0_rp
            !  do igaus = 1,pgaus
            !     fact = gpvol(igaus) / ( gpden(igaus) * dtinv )
            !     T = T + fact
            !     do inode = 1,pnode
            !        eldtcprho(inode) = eldtcprho(inode) + gpsha(inode,igaus)**2 * fact
            !     end do
            !  end do
            !  do inode = 1,pnode
            !     d = d + eldtcprho(inode)
            !  end do
            !  do inode = 1,pnode
            !     eldtcprho(inode) = eldtcprho(inode) * T / d
            !  end do
    
            elmat(DEF_VECT,:,:) = 0.0_rp
    !@        do inode=1,pnode
    !@           do jnode=1,pnode
    !@              elmat(inode,jnode)=0.0_rp
    !@           end do
    !@        end do
    
            do igaus=1,pgaus
               do inode=1,pnode
                  fact(DEF_VECT)=gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,inode,igaus)/(gpden(DEF_VECT,igaus)*dtinv)
                  do jnode=1,pnode
                     elmat(DEF_VECT,inode,jnode)=elmat(DEF_VECT,inode,jnode) +fact*gpsha(DEF_VECT,jnode,igaus)
                  end do     
               end do
            end do
    
            trace(DEF_VECT)  = 0.0_rp
            elmass(DEF_VECT) = 0.0_rp
            do inode = 1,pnode                       
               trace(DEF_VECT) = trace(DEF_VECT) + elmat(DEF_VECT,inode,inode)
               do jnode = 1,pnode                       
                  elmass(DEF_VECT) = elmass(DEF_VECT) + elmat(DEF_VECT,inode,jnode)
               end do
            end do
    
            !
            ! Nodal projection
            !
            do ivect = 1, VECTOR_DIM
              ielem = list_elements_p(ivect)
              if(ivect == 1 .or. (ivect > 1 .and. ielem /= ielemone)) then
                 do inode = 1,pnode
                    ipoin = lnods_loc(ivect,inode)
                    dt_rho_cp_tem_loc(ipoin) = dt_rho_cp_tem_loc(ipoin) + elmat(ivect,inode,inode)*(elmass(ivect)/trace(ivect)) 
                 end do
               end if
            end do
         end if
      endif
    end subroutine tem_rhocpdt_fast

    subroutine tem_turbul_fast(&
         VECTOR_DIM,pnode,pgaus,gpsha,gpcon,gpsph, &
         gpdif,gpgrd,gpden,gptur)   
      !-----------------------------------------------------------------------
      !****f* Temper/tem_turbul
      ! NAME 
      !    tem_turbul
      ! DESCRIPTION
      !  Coupling with TURBUL
      !    Compute effective diffusion coefficient rho*(D+D_t)
      ! USES
      ! USED BY
      !    tem_elmope_new
      !***
      !-----------------------------------------------------------------------
      use def_kintyp, only      : ip,rp
      use def_domain, only      : ndime
      use def_kermod, only      : turmu_ker
      use def_temper, only      : prtur_tem,kfl_grdif_tem, &
                                  kfl_regim_tem, kfl_rhs_scal_tem
      use def_master, only      : div_enthalpy_transport,kfl_coupl,ID_TEMPER,ID_CHEMIC, &
                                  kfl_htran
    
      implicit none
      integer(ip), intent(in)   :: VECTOR_DIM,pnode,pgaus!,igaui,igauf
      real(rp),    intent(in)   :: gpsha(VECTOR_DIM,pnode,pgaus)
      real(rp),    intent(in)   :: gpcon(VECTOR_DIM,pgaus)
      real(rp),    intent(in)   :: gpsph(VECTOR_DIM,pgaus)
      real(rp),    intent(in)   :: gpden(VECTOR_DIM,pgaus)
      real(rp),    intent(inout):: gptur(VECTOR_DIM,pgaus)
      real(rp),    intent(inout):: gpgrd(VECTOR_DIM,ndime,pgaus)
      real(rp),    intent(out)  :: gpdif(VECTOR_DIM,pgaus)
      !
      ! ENTHALPY EQUATION
      ! 
      if (kfl_regim_tem == 4) then
        
        !
        ! Laminar contribution (rho*D)
        !
        if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 .and. associated(div_enthalpy_transport) .and. kfl_htran == 1) then
           gpdif(DEF_VECT,:) = 0.0_rp   ! In this case div_enthalpy_transport is computed in Chemic and no additional diffusion term is required for the energy eq.
        else
           gpdif(DEF_VECT,:) = gpcon(DEF_VECT,:) / gpsph(DEF_VECT,:)
        end if
        ! 
        ! Turbulent contribution (rho*D_t) for RANS & LES
        !
        if( turmu_ker % kfl_exist /= 0_ip ) then
          !
          ! Compute mu_t for LES
          !
          gptur(DEF_VECT,:) = gptur(DEF_VECT,:) * gpden(DEF_VECT,:)
          !
          ! Effective diffusion coefficient: rho*(D+D_t)
          !
          gpdif(DEF_VECT,:) = gpdif(DEF_VECT,:) + gptur(DEF_VECT,:) / prtur_tem
     
        end if
    
      !
      ! TEMPERATURE EQUATION 
      !
      else
        !
        ! Laminar contribution (k)
        !
        gpdif(DEF_VECT,:) = gpcon(DEF_VECT,:)
        ! 
        ! Turbulent contribution (k_t) for RANS & LES
        !
        if( turmu_ker % kfl_exist /= 0_ip ) then
          !
          ! Compute mu_t for LES
          !
          gptur(DEF_VECT,:) = gptur(DEF_VECT,:) * gpden(DEF_VECT,:)
          !
          ! Effective diffusion coefficient: (k + k_t)
          !
          gpdif(DEF_VECT,:) = gpdif(DEF_VECT,:) + gpsph(DEF_VECT,:)*gptur(DEF_VECT,:)/prtur_tem
    
    
        end if
    
    
      end if
    
      if ( kfl_rhs_scal_tem > 0 ) gpdif(DEF_VECT,:) = gpdif(DEF_VECT,:) / gpden(DEF_VECT,:)  
    
    
      if(kfl_grdif_tem/=0) then
         !
         ! GPGRD=grad(k+kt)
         !
         if( turmu_ker % kfl_exist /= 0_ip ) call runend('TEM_TURBUL: GRADIENT OF TURBULENT VISCOSITY NOT CODED')
    
      endif
    
    end subroutine tem_turbul_fast

    subroutine tem_velfun_fast(VECTOR_DIM,pnode,coord,vefun)
      !-----------------------------------------------------------------------
      !****f* Temper/tem_velfun fast
      ! NAME
      !   tem_velfun
      ! DESCRIPTION
      !   Compute velocity according to the function number 
      ! INPUT
      !   KFL_ADVEC_TEM ... Function
      !   COORD ........... Coordinates
      !   NPOIN ........... Number of nodes or Gauss points
      ! OUTPUT 
      !   VEFUN ........... Velocity
      ! USES
      ! USED BY
      !    tem_elmope
      !***
      !-----------------------------------------------------------------------
      use def_kintyp, only       :  ip,rp 
      use def_temper, only       :  kfl_advec_tem
      use def_master, only       :  veloc
      use def_domain, only       :  ndime
      implicit none
      integer(ip), intent(in)    :: VECTOR_DIM
      integer(ip), intent(in)    :: pnode
      real(rp),    intent(in)    :: coord(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(out)   :: vefun(VECTOR_DIM,ndime,pnode)
      integer(ip)                :: idime,ipoin,ivect
      real(rp)                   :: x,y,theta
    
    
      if(kfl_advec_tem==0) then
        vefun(DEF_VECT,:,:) = 0.0_rp
         
      else if(kfl_advec_tem==1) then
         do ivect = 1,VECTOR_DIM
           do ipoin = 1,pnode
              do idime = 1,ndime
                 vefun(ivect,idime,ipoin)= veloc(idime,ipoin,1)
              end do
           end do
         end do
         
      else if(kfl_advec_tem==2) then
        do ivect = 1,VECTOR_DIM
           do ipoin=1,pnode
              x=coord(ivect,1,ipoin)
              y=coord(ivect,2,ipoin)
              vefun(ivect,1,ipoin)= 0.5_rp*(1.0_rp-x*x)*(1.0_rp+y)
              vefun(ivect,2,ipoin)=-0.5_rp*x*(4.0_rp-(1.0_rp+y)**2)
           end do
        end do
    
      else if(kfl_advec_tem==3) then
          vefun(DEF_VECT,1,:) =  1.0_rp
          vefun(DEF_VECT,2,:) = -1.0_rp
          vefun(DEF_VECT,ndime,:) = -1.0_rp
    
      else if(kfl_advec_tem==4) then
           vefun(DEF_VECT,1,:) = 1.0_rp
           if(ndime>1) vefun(DEF_VECT,2:ndime,:) = 0.0_rp
    
      else if(kfl_advec_tem==5) then
         theta = 1.1071487_rp  ! tan-1(2)
         vefun(DEF_VECT,1,:)=  cos(theta)
         vefun(DEF_VECT,2,:)= -sin(theta)
            
    
      else if(kfl_advec_tem==6) then
          vefun(DEF_VECT,1:2,:) = 1.0_rp/sqrt(2.0_rp)
     
      else if(kfl_advec_tem==7) then
         vefun(DEF_VECT,1,:)       = 1.0_rp
         vefun(DEF_VECT,2:ndime,:) = 0.0_rp
    
      else if(kfl_advec_tem==8) then
         if(ndime == 2 ) then
            vefun(DEF_VECT,:,:) = 1.0_rp
         else
            vefun(DEF_VECT,1,:) = 1.0_rp
            vefun(DEF_VECT,2,:) = 0.0_rp
            vefun(DEF_VECT,3,:) = 1.0_rp
         end if
    
      else if(kfl_advec_tem==9) then
         vefun(DEF_VECT,1,:) = 0.0_rp
         vefun(DEF_VECT,2,:) = 0.0_rp
         vefun(DEF_VECT,ndime,:) = 1.0_rp
       
    
      else if(kfl_advec_tem==10) then
         do ivect = 1, VECTOR_DIM
           do ipoin = 1,pnode
              vefun(ivect,1,ipoin) = -(coord(ivect,2,ipoin)-0.5_rp)
              vefun(ivect,2,ipoin) =  (coord(ivect,1,ipoin)-0.5_rp)
           end do
         end do
    
      else if(kfl_advec_tem==11) then
            
         vefun(DEF_VECT,1:ndime-1,:) = 0.0_rp
         vefun(DEF_VECT,ndime,:)     = 1.0_rp
    
    
      end if
    
    end subroutine tem_velfun_fast

    !-----------------------------------------------------------------------
    !> @addtogroup Tempre Fast
    !> @{
    !> @file    tem_elmpre_fast.f90
    !> @author  houzeaux
    !> @date    2020-04-02
    !> @brief   Compute some Gauss values
    !> @details The variables are:
    !>          GPGRD(NDIME) ... grad(k) coefficient
    !>          GPTEM .......... Temperature of previous iterations and time steps
    !>          GPVEL .......... Advection a
    !>          GPADV .......... Advection term a.grad(Ni)
    !>          GPRGD .......... Thermal conductivity gradient grad(k+kt)
    !> @} 
    !-----------------------------------------------------------------------
    
    subroutine tem_elmpre_fast(&
         VECTOR_DIM,pnode,pgaus,pmate,gpden,gpsph,gpsgv,&
         gpsha,gpcar,gphes,elvel,eltem,elcod,elmsh,&
         gpvel,gptem,gprhs,gpcod,gpgrt,lnods,list_elements)
     
      use def_kintyp,                  only : ip,rp
      use def_master,                  only : dpthe,cutim,velom
      use def_domain,                  only : mnode,ndime,ntens
      use def_domain,                  only : xfiel
      use def_temper,                  only : kfl_advec_tem,ADR_tem
      use def_temper,                  only : kfl_sourc_tem
      use def_temper,                  only : kfl_regim_tem
      use def_temper,                  only : sourc_tem,kfl_rhs_scal_tem
      use def_temper,                  only : kfl_sonum_tem
      use def_temper,                  only : lsour_material_tem
      use def_temper,                  only : xsour_material_tem
      use def_temper,                  only : SOURCE_TERM_SPACE_TIME
      use def_temper,                  only : SOURCE_TERM_FIELD
      use def_temper,                  only : SOURCE_TERM_MATERIAL
      use mod_ADR,                     only : BDF
      use def_kermod,                  only : gasco 
      use mod_ker_space_time_function, only : ker_space_time_function
      
      implicit none
      
      integer(ip), intent(in)    :: VECTOR_DIM
      integer(ip), intent(in)    :: pnode
      integer(ip), intent(in)    :: pgaus
      integer(ip), intent(in)    :: pmate
      integer(ip), intent(in)    :: lnods(VECTOR_DIM,pnode)
      real(rp),    intent(in)    :: gpsph(VECTOR_DIM,pgaus)
      real(rp),    intent(in)    :: gpsgv(VECTOR_DIM,ndime,pgaus)
      real(rp),    intent(in)    :: gpsha(VECTOR_DIM,pnode,pgaus)
      real(rp),    intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
      real(rp),    intent(in)    :: gphes(VECTOR_DIM,ntens,mnode,pgaus)
      real(rp),    intent(in)    :: elvel(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(in)    :: eltem(VECTOR_DIM,pnode,*)
      real(rp),    intent(in)    :: elcod(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(in)    :: elmsh(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(inout) :: gpden(VECTOR_DIM,pgaus)
      real(rp),    intent(out)   :: gpvel(VECTOR_DIM,ndime,pgaus)
      real(rp),    intent(out)   :: gpcod(VECTOR_DIM,ndime,pgaus)
      real(rp),    intent(out)   :: gprhs(VECTOR_DIM,pgaus)
      real(rp),    intent(out)   :: gptem(VECTOR_DIM,pgaus,*)
      real(rp),    intent(out)   :: gpgrt(VECTOR_DIM,ndime,pgaus)
      integer(ip), intent(in)    :: list_elements(VECTOR_DIM)                        !< List of elements
     
      integer(ip)                :: idime,inode,igaus,mxdim,ivect,ielem
      real(rp)                   :: gpmsh(VECTOR_DIM,ndime,pgaus)

      gpmsh = 0.0_rp
      !
      ! Temperature: GPTEM
      !
      if( ADR_tem % kfl_time_integration /= 0 ) then
         do igaus = 1,pgaus
            gptem(DEF_VECT,igaus,2) = 0.0_rp
            do inode = 1,pnode
               gptem(DEF_VECT,igaus,2) = gptem(DEF_VECT,igaus,2) + gpsha(DEF_VECT,inode,igaus) * eltem(DEF_VECT,inode,2)
            end do
         end do
      end if
      !
      ! Density: GPDEN=rho*Cp
      !
      if( kfl_regim_tem == 1 ) then !Este se puede comentar se usan 3y4
         do igaus = 1,pgaus
            gpden(DEF_VECT,igaus) = gpden(DEF_VECT,igaus) * (gpsph(DEF_VECT,igaus)-gasco)     ! rho*Cv=rho*(Cp-R)
         end do
      else if (kfl_regim_tem /= 4) then
         do igaus = 1,pgaus
            gpden(DEF_VECT,igaus) = gpden(DEF_VECT,igaus) * gpsph(DEF_VECT,igaus)             ! rho*Cp
         end do
      end if
      
     !@ if ( kfl_rhs_scal_tem > 0 )  gpden(1:pgaus) = 1.0_rp
      if ( kfl_rhs_scal_tem > 0 )  gpden(DEF_VECT,:) = 1.0_rp
      !
      ! Coordinates
      !
      do igaus = 1,pgaus
         do idime = 1,ndime
            gpcod(DEF_VECT,idime,igaus) = 0.0_rp
         end do
         do inode = 1,pnode
            do idime = 1,ndime
               gpcod(DEF_VECT,idime,igaus) = gpcod(DEF_VECT,idime,igaus)&
                    + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
            end do
         end do
      end do
      !
      ! Velocity GPVEL=a and advection term GPADV=a.grad(Ni)-(k/r,0).grad(Ni)
      !
    
      gpvel(DEF_VECT,:,:) = 0.0_rp
     
    
      if( kfl_advec_tem /= 0 ) then
    
         if( kfl_advec_tem == 1 ) then
            do igaus = 1,pgaus
               do inode = 1,pnode
                  do idime = 1,ndime
                     gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) &
                          + gpsha(DEF_VECT,inode,igaus) * elvel(DEF_VECT,idime,inode)
                  end do
               end do
            end do
         else if( kfl_advec_tem >= 2 ) then
            call tem_velfun_fast(VECTOR_DIM,pgaus,gpcod,gpvel)
    
         end if
    
         if( associated(velom) ) then
            gpmsh = 0.0_rp
            do igaus = 1,pgaus
               do inode = 1,pnode
                  do idime = 1,ndime
                     gpmsh(DEF_VECT,idime,igaus) = gpmsh(DEF_VECT,idime,igaus)&
                          + gpsha(DEF_VECT,inode,igaus) * elmsh(DEF_VECT,idime,inode)
                  end do
               end do
               gpvel(DEF_VECT,1:ndime,igaus) = gpvel(DEF_VECT,1:ndime,igaus) - gpmsh(DEF_VECT,1:ndime,igaus)
            end do        
         end if
      end if
      !
      ! Source term: GPRHS
      !
      if( kfl_sourc_tem == SOURCE_TERM_MATERIAL ) then
         !
         ! Material
         !
         if( lsour_material_tem(pmate) == 2 ) then
            do igaus = 1,pgaus
               gprhs(DEF_VECT,igaus) = xsour_material_tem(1,pmate) * sourc_tem
            end do
         end if
         
      else if( kfl_sourc_tem == SOURCE_TERM_SPACE_TIME ) then
         !
         ! Space time
         !
         do ivect=1,VECTOR_DIM
            mxdim = min(2_ip,ndime)
            do igaus = 1,pgaus
               call ker_space_time_function(&
                    kfl_sonum_tem,gpcod(ivect,1,igaus),gpcod(ivect,mxdim,igaus),gpcod(ivect,ndime,igaus),cutim,gprhs(ivect,igaus))
               gprhs(ivect,igaus) = gprhs(ivect,igaus) * sourc_tem
            end do
         end do
         
      else if( kfl_sourc_tem == SOURCE_TERM_FIELD ) then
         !
         ! Field
         !
         do ivect=1,VECTOR_DIM
            ielem=list_elements(ivect) 
            do igaus = 1,pgaus
               gprhs(ivect,igaus) = xfiel(kfl_sonum_tem) % a(1,ielem,1) * sourc_tem
            end do
         end do     
      else
         gprhs(:,:) = 0.0_rp     
      end if
      !
      ! Low-Mach: dp0/dt
      ! prefactor alpha * T ~ 1 
      !
      if( kfl_regim_tem >= 3 ) then
         do igaus = 1,pgaus
            gprhs(DEF_VECT,igaus) = gprhs(DEF_VECT,igaus) + dpthe 
         end do
    
      end if
    
    end subroutine tem_elmpre_fast

    subroutine tem_elmlen_fast(VECTOR_DIM,ndime,pnode,deriv,tragl,elcod,hnatu,hleng)
      use def_kintyp, only     : ip,rp
      implicit none
      integer(ip), intent(in)  :: VECTOR_DIM                                       !< Number of nodes
      integer(ip), intent(in)  :: ndime                                            !< 
      integer(ip), intent(in)  :: pnode                                            !< Number of nodes
      real(rp),    intent(out) :: tragl(VECTOR_DIM,ndime,ndime)
      real(rp),    intent(in)  :: hnatu,elcod(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(in)  :: deriv(ndime,pnode)
      real(rp),    intent(out) :: hleng(VECTOR_DIM,ndime)


!  integer(ip), intent(in)  :: ndime,pnode
!  real(rp),    intent(out) :: tragl(VECTOR_DIM,ndime,ndime)
!  real(rp),    intent(in)  :: hnatu,elcod(VECTOR_DIM,ndime,pnode)
!  real(rp),    intent(in)  :: deriv(ndime,pnode)
!  real(rp),    intent(out) :: hleng(VECTOR_DIM,ndime)

  integer(ip)              :: k,ivect
  real(rp)                 :: enor0(VECTOR_DIM),h_tem,gpdet(VECTOR_DIM),denom(VECTOR_DIM)
  real(rp)                 :: xjacm(VECTOR_DIM,ndime,ndime),t1(VECTOR_DIM),t2(VECTOR_DIM),t3(VECTOR_DIM)

!        call elmlen(ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
!             hnatu(pelty),hleng)


 
      if( ndime == 2 ) then
    
         xjacm(DEF_VECT,1,1) = 0.0_rp
         xjacm(DEF_VECT,1,2) = 0.0_rp
         xjacm(DEF_VECT,2,1) = 0.0_rp
         xjacm(DEF_VECT,2,2) = 0.0_rp
         do k = 1,pnode
            xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k)
            xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * deriv(2,k)
            xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * deriv(1,k)
            xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * deriv(2,k)
         end do
    
         gpdet(DEF_VECT)      =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
         denom(DEF_VECT)      =  1.0_rp/gpdet(DEF_VECT)
         tragl(DEF_VECT,1,1) =  xjacm(DEF_VECT,2,2) * denom(DEF_VECT)
         tragl(DEF_VECT,2,2) =  xjacm(DEF_VECT,1,1) * denom(DEF_VECT)
         tragl(DEF_VECT,2,1) = -xjacm(DEF_VECT,2,1) * denom(DEF_VECT)
         tragl(DEF_VECT,1,2) = -xjacm(DEF_VECT,1,2) * denom(DEF_VECT)
    
         enor0(DEF_VECT)      =  tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1) + tragl(DEF_VECT,1,2) * tragl(DEF_VECT,1,2)
         hleng(DEF_VECT,1)   =  hnatu/sqrt(enor0(DEF_VECT))
         enor0(DEF_VECT)      =  tragl(DEF_VECT,2,1) * tragl(DEF_VECT,2,1) + tragl(DEF_VECT,2,2) * tragl(DEF_VECT,2,2)
         hleng(DEF_VECT,2)   =  hnatu/sqrt(enor0(DEF_VECT))
        
         do ivect = 1,VECTOR_DIM
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem    = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
         end do
        else if( ndime == 3 ) then
     
          xjacm(DEF_VECT,1,1) = 0.0_rp ! xjacm = elcod * deriv^t
          xjacm(DEF_VECT,1,2) = 0.0_rp ! tragl = xjacm^-1
          xjacm(DEF_VECT,1,3) = 0.0_rp
          xjacm(DEF_VECT,2,1) = 0.0_rp
          xjacm(DEF_VECT,2,2) = 0.0_rp
          xjacm(DEF_VECT,2,3) = 0.0_rp
          xjacm(DEF_VECT,3,1) = 0.0_rp
          xjacm(DEF_VECT,3,2) = 0.0_rp
          xjacm(DEF_VECT,3,3) = 0.0_rp
          do k = 1,pnode
             xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k)
             xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * deriv(2,k)
             xjacm(DEF_VECT,1,3) = xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * deriv(3,k)
             xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * deriv(1,k)
             xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * deriv(2,k)
             xjacm(DEF_VECT,2,3) = xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * deriv(3,k)
             xjacm(DEF_VECT,3,1) = xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * deriv(1,k)
             xjacm(DEF_VECT,3,2) = xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * deriv(2,k)
             xjacm(DEF_VECT,3,3) = xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * deriv(3,k)
          end do
     
          t1(DEF_VECT)         =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
          t2(DEF_VECT)         = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
          t3(DEF_VECT)         =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)
    
          gpdet(DEF_VECT)      =  xjacm(DEF_VECT,1,1) * t1(DEF_VECT) + xjacm(DEF_VECT,1,2) * t2(DEF_VECT) + xjacm(DEF_VECT,1,3) * t3(DEF_VECT)
          denom(DEF_VECT)      =  1.0_rp / gpdet(DEF_VECT)
          tragl(DEF_VECT,1,1) =  t1(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,2,1) =  t2(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,3,1) =  t3(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,2,2) =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,3,2) =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom(DEF_VECT)
          tragl(DEF_VECT,3,3) =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom(DEF_VECT)
          tragl(DEF_VECT,1,2) =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,1,3) =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,2,3) =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          !
          ! Element length HLENG
          !
          enor0(DEF_VECT)    = tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1) &
                &     + tragl(DEF_VECT,1,2) * tragl(DEF_VECT,1,2) &
                &     + tragl(DEF_VECT,1,3) * tragl(DEF_VECT,1,3)
          hleng(DEF_VECT,1)  = hnatu/sqrt(enor0(DEF_VECT))
          enor0(DEF_VECT)    = tragl(DEF_VECT,2,1) * tragl(DEF_VECT,2,1) &
                &     + tragl(DEF_VECT,2,2) * tragl(DEF_VECT,2,2) &
                &     + tragl(DEF_VECT,2,3) * tragl(DEF_VECT,2,3)
          hleng(DEF_VECT,2)  = hnatu/sqrt(enor0(DEF_VECT))
          enor0(DEF_VECT)    = tragl(DEF_VECT,3,1) * tragl(DEF_VECT,3,1) &
                &     + tragl(DEF_VECT,3,2) * tragl(DEF_VECT,3,2) &
                &     + tragl(DEF_VECT,3,3) * tragl(DEF_VECT,3,3)
          hleng(DEF_VECT,3)  = hnatu/sqrt(enor0(DEF_VECT))
          !
          ! Sort hleng: hleng(1)=max; hleng(ndime)=min
          !
         do ivect = 1,VECTOR_DIM
     
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem    = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
          if( hleng(ivect,3) > hleng(ivect,1) ) then
             h_tem    = hleng(ivect,3)
             hleng(ivect,3) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
          if( hleng(ivect,3) > hleng(ivect,2) ) then
             h_tem    = hleng(ivect,3)
             hleng(ivect,3) = hleng(ivect,2)
             hleng(ivect,2) = h_tem
          end if
        end do 

      end if


    end subroutine tem_elmlen_fast

    subroutine tem_element_assembly(&
            VECTOR_DIM,pnode,pgaus,list_elements_p,&
            elcod,gpsha,gpcar,gpder,gphes,gpvol,chale,&
            cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gptem,eltem,elmat,elrhs,&
            messa, elvel)
    use def_kintyp, only     : ip,rp
    use def_temper, only     :  ADR_tem
    use def_domain, only     : mnode,ndime,ntens
   !
    ! Element dimensions
    !
    integer(ip), intent(in)            :: VECTOR_DIM                         !< Current element
    integer(ip), intent(in)            :: pnode                         !< # nodes
    integer(ip), intent(in)            :: pgaus                         !< # Gauss points
    integer(ip), intent(in)            :: list_elements_p(VECTOR_DIM)                        !< List of elements

   !
    ! Element characteristics at Gauss point
    !
    real(rp),    intent(in)            :: elcod(VECTOR_DIM,ndime,pnode)            !< Element node coordinates
    real(rp),    intent(in)            :: gpsha(VECTOR_DIM,pnode,pgaus)            !< Shape function Nk
    real(rp),    intent(in)            :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)      !< Shape function Cartesian derivatives dNk/dxi
    real(rp),    intent(in)            :: gpder(VECTOR_DIM,ndime,pnode,pgaus)      !< Shape function derivatives DNk/dsi
    real(rp),    intent(in)            :: gphes(VECTOR_DIM,ntens,mnode,pgaus)      !< Hessian dNk/dxidxj
    real(rp),    intent(in)            :: gpvol(VECTOR_DIM,pgaus)                  !< Element Jacobian
    real(rp),    intent(in)            :: chale(VECTOR_DIM,2)                      !< Element characteristic length
    !
    ! Equation coefficients
    !
    real(rp),    intent(in)            :: cutim                         !< Current time
    real(rp),    intent(in)            :: gpden(VECTOR_DIM,pgaus)                  !< Density
    real(rp),    intent(in)            :: gpvel(VECTOR_DIM,ndime,pgaus)            !< Advection vector
    real(rp),    intent(in)            :: gpdif(VECTOR_DIM,pgaus)                  !< Diffusion 
    real(rp),    intent(in)            :: gpgrd(VECTOR_DIM,ndime,pgaus)            !< Diffusion gradient
    real(rp),    intent(in)            :: gprea(VECTOR_DIM,pgaus,4)                !< Reaction
    real(rp),    intent(in)            :: gprhs(VECTOR_DIM,pgaus)                  !< RHS
    real(rp),    intent(in)            :: gptem(VECTOR_DIM,pgaus,ADR_tem % ntime)                !< Unknown at Gauss point
    real(rp),    intent(in)            :: eltem(VECTOR_DIM,pnode,2)                !< Element unknown
    !
    ! Output
    !
    real(rp),    intent(out)           :: elmat(VECTOR_DIM,pnode,pnode)            !< Element matrix
    real(rp),    intent(out)           :: elrhs(VECTOR_DIM,pnode)                  !< Element RHS
    character(*),intent(in),  optional :: messa                         !< Message
    real(rp),    intent(in),  optional :: elvel(VECTOR_DIM,ndime, pnode)           !< Element velocities

    ! To assembly of the matrix

    real(rp)    :: react(VECTOR_DIM,pgaus),sreac(VECTOR_DIM,pgaus),sgs(VECTOR_DIM,pgaus)
    real(rp)    :: fact1(VECTOR_DIM), rhsit(VECTOR_DIM, pgaus), gptau(VECTOR_DIM,pgaus), gptau_time(VECTOR_DIM,pgaus)
    real(rp)    :: gpnve(VECTOR_DIM), ka(VECTOR_DIM), aa(VECTOR_DIM),sa(VECTOR_DIM), tau(VECTOR_DIM), freq1(VECTOR_DIM),freq2(VECTOR_DIM),freq3(VECTOR_DIM)
    real(rp)    :: dtinv_elem
    integer(ip) :: itime,idime
    real(rp)    :: resi1(VECTOR_DIM,pnode), resi2(VECTOR_DIM,pnode), gpad1(VECTOR_DIM,pnode)
    real(rp)    :: gpadv(VECTOR_DIM,pnode), gppe1(VECTOR_DIM,pnode), gppe2(VECTOR_DIM,pnode)
    real(rp)    :: grvgr(VECTOR_DIM),gplap(VECTOR_DIM), d(VECTOR_DIM), CD(VECTOR_DIM,pgaus), SD(VECTOR_DIM,pgaus)
    real(rp)    :: fact2(VECTOR_DIM), fact3(VECTOR_DIM),xmuit(VECTOR_DIM),xmui3(VECTOR_DIM),alpha(VECTOR_DIM),alpha1(VECTOR_DIM)
    real(rp)    :: gppr1(VECTOR_DIM,pgaus), gpdiv(VECTOR_DIM)

    integer(ip) :: ielem,igaus,inode,jnode                     ! Indices and dimensions
    integer(ip) :: ivect,ielemone
 
  
    ielemone= list_elements_p(1) !! Just to avoid repeating terms in pseudo vectorized version


    ! Initialization 
    elrhs(DEF_VECT,1:pnode)         = 0.0_rp
    elmat(DEF_VECT,1:pnode,1:pnode) = 0.0_rp
    gppr1(DEF_VECT,1:pgaus)         = 0.0_rp
    dtinv_elem = ADR_tem % dtinv
    if (ADR_tem % kfl_time_scheme == 3 .or. ADR_tem % kfl_time_scheme == 4) then
       dtinv_elem = 0.0_rp
    end if

    rhsit(DEF_VECT,1:pgaus) = gprhs(DEF_VECT,1:pgaus)

    do igaus = 1,pgaus
        fact1(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_elem
        do itime = 2,ADR_tem % ntime
           rhsit(DEF_VECT,igaus) = rhsit(DEF_VECT,igaus) - fact1(DEF_VECT) * (ADR_tem % time_parameters(itime))*gptem(DEF_VECT,igaus,itime)
        end do
    end do

    react(DEF_VECT,1:pgaus) = gprea(DEF_VECT,1:pgaus,1) + 2.0_rp * gprea(DEF_VECT,1:pgaus,2) * gptem(DEF_VECT,1:pgaus,1)
    sreac(DEF_VECT,1:pgaus) = gprea(DEF_VECT,1:pgaus,1) + gprea(DEF_VECT,1:pgaus,2) * gptem(DEF_VECT,1:pgaus,1)
    rhsit(DEF_VECT,1:pgaus) = rhsit(DEF_VECT,1:pgaus)   + gprea(DEF_VECT,1:pgaus,2) * (gptem(DEF_VECT,1:pgaus,1)) ** 2
    if( ADR_tem % kfl_stabilization == -3 ) sreac(DEF_VECT,:) = 0.0_rp

    elements44: do ivect = 1,VECTOR_DIM
     ielem=list_elements_p(ivect)
     if (ivect==1 .or. (ivect > 1 .and. ielem /= ielemone ))  then
        ! Element characteristic length
        if( chale(ivect,1) == 0.0_rp .and. chale(ivect,2) == 0.0_rp ) then
           gptau(ivect,:) = 0.0_rp
        end if
      end if
    end do elements44


    do igaus = 1,pgaus
        gpnve(DEF_VECT)=0.0_rp
       do idime=1, ndime
            gpnve(DEF_VECT) = gpnve(DEF_VECT) + gpvel(DEF_VECT,idime,igaus)*gpvel(DEF_VECT,idime,igaus)
       end do

       gpnve(DEF_VECT) = gpden(DEF_VECT,igaus) * sqrt(gpnve(DEF_VECT)+epsilon(1.0_rp))
    
       ka(DEF_VECT)     = ADR_tem % tau_parameters(1) * gpdif(DEF_VECT,igaus) ! Diffusion 
       aa(DEF_VECT)     = ADR_tem % tau_parameters(2) * gpnve(DEF_VECT)        ! Advection
       sa(DEF_VECT)     = ADR_tem % tau_parameters(3) * sreac(DEF_VECT,igaus) ! Reaction 
       !
       ! Codina
       !
       freq1(DEF_VECT) = 4.0_rp*ka(DEF_VECT)/(chale(DEF_VECT,2)*chale(DEF_VECT,2))
       freq2(DEF_VECT) = 2.0_rp*aa(DEF_VECT)/chale(DEF_VECT,1)
       freq3(DEF_VECT) = abs(sa(DEF_VECT))
       tau(DEF_VECT)   = freq1(DEF_VECT)+freq2(DEF_VECT)+freq3(DEF_VECT)

     elements45: do ivect = 1,VECTOR_DIM
        ielem=list_elements_p(ivect)
        if (ivect==1 .or. (ivect > 1 .and. ielem /= ielemone ))  then
 
          if(tau(ivect)/=0.0_rp) tau(ivect)=1.0_rp/tau(ivect)
        end if
       end do elements45

       gptau(DEF_VECT,igaus) = tau(DEF_VECT)
    end do

    CD(DEF_VECT,:) = 0.0_rp
    SD(DEF_VECT,:) = 0.0_rp
    gptau_time(DEF_VECT,:) = gptau(DEF_VECT,:)
    sgs(DEF_VECT,1:pgaus)  = 0.0_rp 
 
    do igaus = 1,pgaus
       !
       ! Calculus of residual resid and perturbation function Pi=gppre
       !
       ! RESI1 = r1(u) =  rho/dt*Nj + rho*a.grad(Nj) + r1*Nj
       ! RESI2 = r2(u) =  -grad(k).grad(Nj) - k*Lap(Nj)
       ! GPPE1 = p1(v) =  [  Ni*(1-tau*s) + tau*rho*a.grad(Ni) ] 
       ! GPPE2 = p2(v) =  [    -Ni*tau*s  + tau*rho*a.grad(Ni) ] = p1(v) - v * |dv|
       !
       tau(DEF_VECT)   = gptau_time(DEF_VECT,igaus)
       fact1(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_elem
       d(DEF_VECT)     = fact1(DEF_VECT) * sgs(DEF_VECT,igaus)  ! rho*u'n/dt
   
       do inode = 1,pnode
          resi1(DEF_VECT,inode) = fact1(DEF_VECT) * ADR_tem % time_parameters(1) * gpsha(DEF_VECT,inode,igaus)
          gpad1(DEF_VECT,inode) = 0.0_rp
          grvgr(DEF_VECT)        = 0.0_rp
          gplap(DEF_VECT)        = 0.0_rp
          do idime = 1,ndime
             gpad1(DEF_VECT,inode)  = gpad1(DEF_VECT,inode)  + gpvel(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)           
             grvgr(DEF_VECT)        = grvgr(DEF_VECT)        + gpgrd(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)          
             gplap(DEF_VECT)        = gplap(DEF_VECT)        + gphes(DEF_VECT,idime,inode,igaus)
          end do
          gpadv(DEF_VECT,inode) = gpad1(DEF_VECT,inode) * gpden(DEF_VECT,igaus)
          resi1(DEF_VECT,inode) = resi1(DEF_VECT,inode) + gpadv(DEF_VECT,inode) + react(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
          resi2(DEF_VECT,inode) = - grvgr(DEF_VECT) - gplap(DEF_VECT) * gpdif(DEF_VECT,igaus)
          gppe1(DEF_VECT,inode) = (  gpsha(DEF_VECT,inode,igaus)*(1.0_rp-tau(DEF_VECT)*sreac(DEF_VECT,igaus)) + tau(DEF_VECT) * gpadv(DEF_VECT,inode) ) * gpvol(DEF_VECT,igaus)   
          !gppe2(inode) = ( -gpsha(inode,igaus)**sreac(igaus) + gpadv(inode) ) * tau * gpvol(igaus)   
          gppe2(DEF_VECT,inode) = gppe1(DEF_VECT,inode) - gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
       end do

       ! 
       ! Diffusion term
       !
       SD(DEF_VECT,:) =0.0
       fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * ( gpdif(DEF_VECT,igaus) + CD(DEF_VECT,igaus) )  ! Diffusion
       fact3(DEF_VECT) = SD(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)                     ! Streamline negative diffusion
       do inode = 1,pnode
          do jnode = 1,inode-1
             xmuit(DEF_VECT) = 0.0_rp
             do idime = 1,ndime
                xmuit(DEF_VECT) = xmuit(DEF_VECT) + gpcar(DEF_VECT,idime,jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
             xmuit(DEF_VECT)              = xmuit(DEF_VECT) * fact2(DEF_VECT)
             xmui3(DEF_VECT)              = gpad1(DEF_VECT,inode) * gpad1(DEF_VECT,jnode) * fact3(DEF_VECT)
             elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) + xmuit(DEF_VECT) + xmui3(DEF_VECT)
             elmat(DEF_VECT,jnode,inode) = elmat(DEF_VECT,jnode,inode) + xmuit(DEF_VECT) + xmui3(DEF_VECT)
          end do
          xmuit(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             xmuit(DEF_VECT) = xmuit(DEF_VECT) + gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          elmat(DEF_VECT,inode,inode) = elmat(DEF_VECT,inode,inode) + xmuit(DEF_VECT) * fact2(DEF_VECT) + gpad1(DEF_VECT,inode) * gpad1(DEF_VECT,inode) * fact3(DEF_VECT)
       end do
       !
       ! bemol
       !          
       do inode = 1,pnode
          fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * ADR_tem % bemol * gpvol(DEF_VECT,igaus)*gpden(DEF_VECT,igaus)
          do jnode = 1,pnode
             fact2(DEF_VECT) = gpadv(DEF_VECT,jnode) * fact1(DEF_VECT)
             elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) - fact2(DEF_VECT)
             elmat(DEF_VECT,jnode,inode) = elmat(DEF_VECT,jnode,inode) - fact2(DEF_VECT)
          end do
       end do
       !
       ! Assembly of the matrix and rhs
       !       
       ! GPPR1 <=> f - sum_i ri*u^i
       ! GPPR1 <=> - [ rho * a - grad(k) ] . grad(u) 
       ! GPPR1 <=> f - rho*(u - u^n)/dt - sum_i ri*u^i - [rho*u - grad(k)].grad(u) - k*lapl(u)
       !
       do inode = 1,pnode
          do jnode = 1,pnode
             elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) &
                  &             + resi1(DEF_VECT,jnode) * gppe1(DEF_VECT,inode) &        
                  &             + resi2(DEF_VECT,jnode) * gppe2(DEF_VECT,inode)
          end do
          elrhs(:,inode)         = elrhs(DEF_VECT,inode) &
               &                + rhsit(DEF_VECT,igaus) * gppe1(DEF_VECT,inode) &
               &                + ( d(DEF_VECT) - gppr1(DEF_VECT,igaus) ) * gppe2(DEF_VECT,inode)
       end do
       !
       ! Skew symmetric advection
       !
       if (ADR_tem % kfl_skewsymm ==1 ) then ! adds + 0.5*rho div(u)*scalar
          gpdiv(DEF_VECT) = 0.0_rp
          ! Calculates velocity divergence
          do inode=1, pnode
             do idime=1, ndime
                gpdiv(DEF_VECT) = gpdiv(DEF_VECT) + gpcar(DEF_VECT,idime, inode, igaus)*elvel(DEF_VECT,idime, inode)
             end do
          end do
          ! Adds the term + 0.5*rho div(u)*scalar
          do inode=1, pnode
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode, igaus)*gpvol(DEF_VECT,igaus)*0.5_rp*gpden(DEF_VECT,igaus)*gpdiv(DEF_VECT)
             do jnode =1, pnode
                elmat(DEF_VECT,inode, jnode) = elmat(DEF_VECT,inode, jnode) + gpsha(DEF_VECT,jnode, igaus)*fact1(DEF_VECT)
             end do
          end do
          ! Another way for skewsymm would be to use bemol =0.5, needing for boundary terms
       end if
       !
       ! Time tracking of subgrid scale
       !
       if( ADR_tem % kfl_time_sgs /= 0 ) then
          xmuit(DEF_VECT)  = gptau_time(DEF_VECT,igaus) / gptau(DEF_VECT,igaus)
          alpha(DEF_VECT)  = xmuit(DEF_VECT) * gpvol(DEF_VECT,igaus)               ! tau'/tau
          alpha1(DEF_VECT) = ( xmuit(DEF_VECT) - 1.0_rp ) * gpvol(DEF_VECT,igaus)  ! tau'/tau - 1
          do inode = 1,pnode
             do jnode = 1,pnode
                elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) &
                     + alpha1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * ( resi1(DEF_VECT,jnode) + resi2(DEF_VECT,jnode) )
             end do
             elrhs(DEF_VECT,inode) = elrhs(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) &
                  * ( alpha1(DEF_VECT) * ( rhsit(DEF_VECT,igaus) - gppr1(DEF_VECT,igaus) ) + alpha(DEF_VECT) * d(DEF_VECT) )
          end do
       end if
       ! 
       ! Lumped mass evolution matrix
       !
       if( ADR_tem % kfl_time_lumped == 1 ) then
          do inode = 1,pnode
             fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_elem
             elmat(DEF_VECT,inode,inode) = elmat(DEF_VECT,inode,inode) + fact1(DEF_VECT)
             elrhs(DEF_VECT,inode) = elrhs(DEF_VECT,inode) - fact1(DEF_VECT) * gptem(DEF_VECT,igaus,3)
             elrhs(DEF_VECT,inode) = elrhs(DEF_VECT,inode) + fact1(DEF_VECT) * eltem(DEF_VECT,inode,2)
             do jnode =1, pnode
                elmat(DEF_VECT,inode, jnode) = elmat(DEF_VECT,inode, jnode) - fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
             end do
          end do
       end if


    end do      

  do igaus=1,pgaus
       ! 
       ! Lumped mass evolution matrix
       !
       if( ADR_tem % kfl_time_lumped == 1 ) then
          do inode = 1,pnode
             fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_elem
             elmat(DEF_VECT,inode,inode) = elmat(DEF_VECT,inode,inode) + fact1(DEF_VECT)
             elrhs(DEF_VECT,inode) = elrhs(DEF_VECT,inode) - fact1(DEF_VECT) * gptem(DEF_VECT,igaus,3)
             elrhs(DEF_VECT,inode) = elrhs(DEF_VECT,inode) + fact1(DEF_VECT) * eltem(DEF_VECT,inode,2)
             do jnode =1, pnode
                elmat(DEF_VECT,inode, jnode) = elmat(DEF_VECT,inode, jnode) - fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
             end do
          end do
       end if

   end do



    end subroutine tem_element_assembly

    subroutine mod_tem_element_operations_fast(VECTOR_DIM,pnode,pgaus,list_elements,order)
      !------------------------------------------------------------------------
      !****f* Temper/tem_elmop2
      ! NAME 
      !    tem_elmop2
      ! DESCRIPTION
      !    ORDER=1:
      !      Temperature equation, elemental operations:
      !      1. Compute elemental matrix and RHS 
      !      2. Impose Dirichlet boundary conditions
      !      3. Assemble them
      !    ORDER=4:
      !      Update the subgrid scale
      ! USES
      ! USED BY
      !    tem_matrix
      !------------------------------------------------------------------------
      use def_parame
      use def_elmtyp
      use def_master
      use def_domain
      use def_kermod
      use mod_ker_proper 
      use def_temper
      use mod_tem_entropy, only : tem_entropy_viscosity
      use def_master, only : ittim
      use mod_ADR,    only : ADR_element_assembly
      use mod_ADR,    only : ADR_bubble_assembly
      use mod_ADR,    only : ADR_projections_and_sgs_assembly
      use mod_ADR,    only : ADR_add_sgs_or_bubble
      use mod_ADR,    only : ELEMENT_ASSEMBLY             ! 1
      use mod_ADR,    only : PROJECTIONS_AND_SGS_ASSEMBLY ! 4
      use mod_ADR,    only : BUBBLE_ASSEMBLY              ! 5
      use mod_ADR,    only : mreac_adr
      use mod_solver, only : solver_assemble_element_matrix
      use mod_matrix, only : matrix_assemble_element_RHS
      use mod_matrix, only : matrix_assemble_element_matrix_to_CSR
      use mod_matrix, only : matrix_assexp
      use mod_element_integration, only : element_shape_function_derivatives_jacobian
!@      use mod_tem_fast
      implicit none
      integer(ip), intent(in) :: VECTOR_DIM                                       !< Number of nodes
      integer(ip), intent(in) :: pnode                                            !< Number of nodes
      integer(ip), intent(in) :: pgaus                                            !< Number of Gauss points
      integer(ip), intent(in) :: list_elements(VECTOR_DIM)                        !< List of elements
      integer(ip), intent(in) :: order                                            ! =2: compute SGS only
    
      real(rp)    :: elmat(VECTOR_DIM,pnode,pnode),elrhs(VECTOR_DIM,pnode)
      integer(ip) :: ielem,igaus,ipoin,inode,jnode                     ! Indices and dimensions
      integer(ip) :: ivect,ielemone,irea
      integer(ip) :: pelty,pmate
      integer(ip) :: plapl,porde,ptopo
      integer(ip) :: kfl_advec_old
    
      real(rp)    :: eltem(VECTOR_DIM,pnode,ADR_tem % ntime)          ! Gather 
      real(rp)    :: elcod(VECTOR_DIM,ndime,pnode)
      real(rp)    :: elvel(VECTOR_DIM,ndime,pnode)
      real(rp)    :: eledd(VECTOR_DIM,pnode)
      real(rp)    :: elmsh(VECTOR_DIM,ndime,pnode)
    
      real(rp)    :: tragl(VECTOR_DIM,ndime,ndime),chave(VECTOR_DIM,ndime,2)    ! Stabilization
      real(rp)    :: chale(VECTOR_DIM,2),hleng(VECTOR_DIM,ndime)
      real(rp)    :: gpvol(VECTOR_DIM,pgaus)                                    ! |J|*w
      real(rp)    :: gprea(VECTOR_DIM,pgaus,mreac_adr)                          ! r
      real(rp)    :: gpvel(VECTOR_DIM,ndime,pgaus)                              ! a
      real(rp)    :: gpcon(VECTOR_DIM,pgaus),gpcod(VECTOR_DIM,ndime,pgaus)      ! k
      real(rp)    :: gpdif(VECTOR_DIM,pgaus),gpgrd(VECTOR_DIM,ndime,pgaus)      ! k+kt, grad(k+kt)
    !@ real(rp)    :: gpdiv(VECTOR_DIM,pgaus)                                    ! Divergence of convection
      real(rp)    :: gprhs(VECTOR_DIM,pgaus)                                    ! f (all terms)
      real(rp)    :: gpden(VECTOR_DIM,pgaus)                                    ! rho and then rho*cp
      real(rp)    :: gpsph(VECTOR_DIM,pgaus)                                    ! cp
      real(rp)    :: gptem(VECTOR_DIM,pgaus,ADR_tem % ntime)          ! T
      real(rp)    :: gpgrt(VECTOR_DIM,ndime,pgaus)                    ! grad(T)
      real(rp)    :: gpsgv(VECTOR_DIM,ndime,pgaus)                    ! u'
      real(rp)    :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
      real(rp)    :: gpder(VECTOR_DIM,ndime,mnode,pgaus)              ! dNk/dsj
      real(rp)    :: gpcar(VECTOR_DIM,ndime,mnode,mgaus)              ! dNk/dxj
      real(rp)    :: gphes(VECTOR_DIM,ntens,mnode,mgaus)              ! dNk/dxidxj
      real(rp)    :: gptur(VECTOR_DIM,pgaus)                          ! Turbulent viscosity
      integer(ip) :: dummi
      real(rp)    :: dtmin
    
      
      integer(ip) :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)
      integer(ip) :: lnods_loc(VECTOR_DIM,pnode)
    
      integer(ip) :: jpoin , izsol ,jcolu 
    
#ifdef EVENT
      call mpitrace_user_function(1)
#endif


      !
      ! Initialization
      !
      gpdif = 0.0_rp
      gpsph = 0.0_rp
      gpden = 0.0_rp
      gprea = 0.0_rp
      gprhs = 0.0_rp
      gptur = 0.0_rp
      gptem = 0.0_rp
      gpgrd = 0.0_rp
      gpvel = 0.0_rp
      gpgrt = 0.0_rp  
    
      ielem = list_elements(1)                            ! Select first element
      pelty = ltype(ielem)
      plapl = llapl(pelty)
      porde = lorde(pelty)
      ptopo = ltopo(pelty)
      pmate = lmate(ielem)
    
      ielemone= ielem !! Just to avoid repeating terms in pseudo vectorized version
    
      do ivect = 1,VECTOR_DIM                      
        ielem = abs(list_elements(ivect))
        if( ielem /= 0 ) then
           list_elements_p(ivect)   = list_elements(ivect)
        else
           list_elements_p(ivect)   = list_elements(1)
        end if
      end do
    
      !
      ! Loop over elements
      !  
      dtmin = 1.0e6_rp
        
      !
      ! Element dimensions
      !
      kfl_advec_old = kfl_advec_tem !Esto preguntar a Guillaume
    
      do ivect = 1, VECTOR_DIM
        ielem = abs(list_elements(ivect))
        if( ielem /= 0 ) then
           lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
           ielem                    = list_elements(ivect)
        else
           lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
           ielem                    = list_elements(1)
        end if
      end do
    
    
      !
      ! Gather operations 
      !
      call tem_elmgat_fast(VECTOR_DIM,pnode,lnods_loc,eltem,elvel,elcod,eledd,elmsh)
    
      !
      ! hleng and tragl at center of gravity
      !
    
      call tem_elmlen_fast(VECTOR_DIM,ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
    
    
      call tem_elmchl_fast(VECTOR_DIM,tragl,hleng,elcod,elvel,chave,chale,pnode,&
                porde,hnatu(pelty),kfl_advec_tem,kfl_ellen_tem)
     
      !
      ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
      !
      if( plapl == 0 ) gphes = 0.0_rp
    
      call element_shape_function_derivatives_jacobian(&
           pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
           elmar(pelty) % deriv,elmar(pelty) % heslo,&
           elcod,gpvol,gpsha,gpder,gpcar,gphes,list_elements=list_elements_p)
    
      !
      ! Temperature: GPTEM
      !
      
      gptem(DEF_VECT,:,1) = 0.0_rp
    
      do inode=1,pnode
         do igaus=1,pgaus
            gptem(DEF_VECT,igaus,1)=gptem(DEF_VECT,igaus,1)&
                 +eltem(DEF_VECT,inode,1)*gpsha(DEF_VECT,inode,igaus)
         end do
      end do
    
      !
      ! Properties: GPDEN, GPDIF, GPGRD and GPREA
      !
      call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) 
      call ker_proper('CONDU','PGAUS',dummi,list_elements_p,gpcon,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) 
      call ker_proper('SPHEA','PGAUS',dummi,list_elements_p,gpsph,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) 
      call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gptur,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) 
    
      if ( kfl_entpred_tem == 1_ip ) then
         gpcon = 0.0_rp
      end if
      
      if(kfl_grdif_tem /= 0)  call ker_proper('GRCON','PGAUS',dummi,list_elements_p,gpgrd,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
    
      !
      ! Coupling with turbul
      !
      call tem_turbul_fast(VECTOR_DIM,pnode,pgaus,gpsha,gpcon,gpsph,gpdif,gpgrd,gpden,gptur)
     
      if ( kfl_entpred_tem == 1_ip ) gpdif(DEF_VECT,:) = 0.0_rp
    
      !
      ! Equation coefficients in GP
      !
      call tem_elmpre_fast(&  
           VECTOR_DIM,pnode,pgaus,pmate,gpden,gpsph,gpsgv,gpsha,gpcar,gphes,elvel,eltem,&
           elcod,elmsh,gpvel,gptem,gprhs,gpcod,gpgrt,lnods_loc,list_elements_p)
      !
      ! Coupling with CHEMIC
      !
      call tem_chemic_fast(VECTOR_DIM,pgaus,gprhs,list_elements_p)      
     
      elements3: do ivect = 1,VECTOR_DIM
          !
          ! Entropy stable viscosity
          !
          if ((kfl_entropy_tem == 1_ip)) then   !!!!.and.(kfl_diven_tem == 0)
             call tem_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,gpsha(ivect,:,:),gpcar(ivect,:,:,:),elvel(ivect,:,:),gpden(ivect,:),hleng(ivect,:),gpvol(ivect,:),gpdif(ivect,:))
             if ( kfl_entpred_tem == 1_ip ) gpdif(ivect,:) = 0.0_rp
          end if
      end do elements3
    
     
      !
      ! Scale source terms with rho (enthalpy) or rho*cp (temperature)
      !
      if ( kfl_rhs_scal_tem > 0 ) then
          do irea=1,mreac_adr
           gprea(DEF_VECT,:,irea) = gprea(DEF_VECT,:,irea) / gpden(DEF_VECT,:)
          end do
          gprhs(DEF_VECT,:)   = gprhs(DEF_VECT,:)   / gpden(DEF_VECT,:)
      end if
    
    
      !Element Assembly
             call tem_element_assembly(&
                  VECTOR_DIM,pnode,pgaus,list_elements_p, elcod,gpsha,gpcar,elmar(pelty) % deriv, &
                  gphes,gpvol,chale,&
                  cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs, & 
                  gptem,eltem,elmat,elrhs, elvel=elvel)
     
    
      !
      ! Projections of rho*cp/dt
      !
      if(kfl_explicit_tem == 1) then
         call tem_rhocpdt_fast(VECTOR_DIM,pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden, &
            list_elements_p,dt_rho_cp_tem)
      end if
     
      !
      ! Prescribe Dirichlet boundary conditions !Hablarlo con Guillaume
      !
      if( solve(1) % kfl_iffix == 0 ) then
         call tem_elmdir_fast(VECTOR_DIM,pnode,lnods_loc,elmat,elrhs,list_elements_p)
      end if
    
      if(kfl_explicit_tem ==1 ) then
      elements5: do ivect = 1,VECTOR_DIM
         ielem=list_elements_p(ivect)
         if (ivect==1 .or. (ivect > 1 .and. ielem /= ielemone ))  then
     
            do inode = 1,pnode
               ipoin = lnods_loc(ivect,inode)
               rhsid(ipoin) = rhsid(ipoin) + elrhs(ivect,inode)
              do jnode = 1,pnode
#ifdef NO_COLORING
         !$OMP ATOMIC
#endif
                 rhsid(ipoin) = rhsid(ipoin) - elmat(ivect,inode,jnode) * eltem(ivect,jnode,1)
              end do
            end do
         end if
      end do elements5
    else
    
      elements6: do ivect = 1,VECTOR_DIM
         ielem=list_elements_p(ivect)
         if (ivect==1 .or. (ivect > 1 .and. ielem /= ielemone ))  then
           
            !
            ! Matrix assemble element RHS
            !
            do inode = 1,pnode
               ipoin = lnods_loc(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
               rhsid(ipoin) = rhsid(ipoin) + elrhs(ivect,inode)
            end do
            !
            ! Matrix assemble element matrix to CSR
            !
            if( kfl_element_to_csr ==1 ) then
               do inode = 1,pnode
                  do jnode = 1,pnode
                     izsol = lezdo(inode,jnode,ielem)
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                     amatr(izsol) = amatr(izsol) + elmat(ivect,inode,jnode)
                  end do
               end do
            else
               do inode = 1,pnode
                  ipoin = lnods_loc(ivect,inode)
                  do jnode = 1,pnode
                     jpoin = lnods_loc(ivect,jnode)
                     izsol = r_dom(ipoin)
                     jcolu = c_dom(izsol)
                     do while( jcolu /= jpoin .and. izsol < r_dom(ipoin+1)-1)
                        izsol = izsol + 1
                        jcolu = c_dom(izsol)
                     end do
      
                     if( jcolu == jpoin ) then
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                        amatr(izsol) = amatr(izsol) + elmat(ivect,inode,jnode)
                     end if
                  end do
               end do
            end if
        end if
       end do elements6
    
    end if
    
#ifdef EVENT
    call mpitrace_user_function(0)
#endif

end subroutine mod_tem_element_operations_fast

end module mod_tem_fast


