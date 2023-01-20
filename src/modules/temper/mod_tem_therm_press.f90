!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_tem_therm_press
  use def_temper
  use def_kintyp,         only : ip,rp
  use mod_messages,       only : messages_live 
  use def_kermod,         only : gasco
  use def_master,         only : zeror,INOTEMPTY,imass,kfl_prthe,prthe
  use def_master,         only : tempe_gp,wmean_gp,tempe,veloc,mass_sink 
  use def_master,         only : modul,mem_modul,dpthe,npoi1
  use def_domain,         only : vodom,elmar,mgaus,ndime,nelem,ltype,nnode
  use def_domain,         only : coord,mnode,ngaus,ntens,lnods,ndimb,nboun 
  use def_domain,         only : lelbo,ltypb,lnodb,mnodb,npoin,mgaub
  use mod_communications, only : PAR_SUM
  use mod_ADR,            only : TRAPEZOIDAL,BDF,ADAMS_BASHFORD,RUNGE_KUTTA
  use mod_memory,         only : memory_alloca,memory_deallo 
  use mod_ker_proper,     only : ker_proper 
  use def_master,         only : IMASTER, kfl_paral
  use mod_parall,         only : commd
  use mod_bouder

  implicit none

  real(rp)                    :: mass_ini=0.0_rp 
  real(rp)                    :: mass_curr=0.0_rp 
  real(rp)                    :: vintwt_ini=0.0_rp 
  real(rp)                    :: vintwt_curr=0.0_rp 
  real(rp)                    :: dmass=0.0_rp 

  integer(ip)                 :: pelty,pnode,pgaus
  integer(ip)                 :: pblty,pnodb,pgaub
  integer(ip)                 :: idime,ielem,igaus,inode,ipoin
  integer(ip)                 :: iboun,igaub,inodb

  public therm_press_update
  private
contains
  subroutine therm_press_update(itask)
     integer(ip), intent(in) :: itask
     
     integer(ip)             :: itime,jtask,nt_save
     real(rp)                :: mass_rati
     !-----------------------------------------------------------------------
     !****f* Temper/mod_therm_press/therm_press_update
     ! NAME 
     !    therm_press_update
     ! DESCRIPTION
     !    This subroutine updates the thermodynamic pressure
     !    of the low Mach number model following the ideal gass law:
     !                                                      
     !       P W                                                
     !     ------- = rho                                        
     !      Ru  T                                               
     !                                                          
     !             /   W   \
     !     P int_V| ------- | = int_V(rho) = M
     !             \ Ru T  /                                    
     !     
     !               M         / W_0 \   /      / W \
     !     P = P_0 ----- int_V| ----- | / int_V| --- |                
     !              M_0        \ T_0 / /        \ T /                 
     !
     !     kfl_prthe:
     !       1: Closed domain: M/M_0 = 1
     !       2: Domain with prescribed in and outlets or spray M/M_0 /= 1
     !
     ! USES
     ! USED BY
     !    tem_endite
     !    tem_endste
     !***
     !-----------------------------------------------------------------------

     !
     ! If thermodynamic pressure is constant, then return 
     !
     if (kfl_prthe == 0_ip) return

     jtask = abs(itask)
     select case(jtask)
        !
        ! Endite
        !
        case(1)
           !-------------------------!
           !                  / W \  !
           ! Calculate  int_V| --- | !
           !                  \ T /  !
           !-------------------------!
           call therm_press_volume_integral( vintwt_curr ) 
           
           !------------!
           ! Initialize !
           !------------!
           if (itask < 0) then
              if( imass > zeror )   then  
                 !
                 ! Set from restart or user input
                 !
                 vintwt_ini      = imass * gasco / prthe(4) * vodom
                 mass_ini        = imass * vodom              
              else 
                 !
                 ! Use vinv_tem(1) as initial 
                 !
                 call messages_live('CALCULATING INITIAL MASS IN DOMAIN')
                 vintwt_ini      = vintwt_curr
                 mass_ini        = vintwt_curr * prthe(4) / gasco  
              end if
              
              mass_curr = mass_ini
           end if

           !-----------------!
           ! Calculate mass: !
           !                 !
           ! M = M + dM      !
           !-----------------!
           mass_rati = 1.0_rp
           select case(kfl_prthe)
              case(2_ip)
                 call therm_press_mass_change(dmass)
                 mass_rati = (mass_curr + dmass)/mass_ini
           end select


           !-------------------------------------------------------!
           ! Update thermodynamic pressure and its time derivative !
           !-------------------------------------------------------!
           prthe(1) = prthe(4) * mass_rati *  vintwt_ini/vintwt_curr
           
           dpthe = 0.0_rp
           if( ADR_tem % kfl_time_scheme == TRAPEZOIDAL .or. ADR_tem % kfl_time_scheme == BDF) then
              do itime =1,ADR_tem % ntime
                 dpthe = dpthe + ADR_tem % time_parameters(itime) * prthe(itime)
              end do
              dpthe     = dpthe * ADR_tem % dtinv
           else
              !
              ! For explicit cases: 
              ! Motheau, E., & Abraham, J. (2016). A high-order numerical algorithm
              ! for DNS of low-Mach-number reactive flows with detailed chemistry and
              ! quasi-spectral accuracy. Journal of Computational Physics, 313,
              ! 430-454.
              !
              dpthe = (prthe(1) - prthe(2))*ADR_tem % dtinv

           end if
           !
           ! Update average density
           !
           xmass_tem = mass_curr      
        
        !
        ! Endste
        !
        case(2) 
           !
           ! Crank-Nicolson
           !
           if( kfl_prthe == 0 ) then
              if( ADR_tem % kfl_time_scheme == TRAPEZOIDAL .and. ADR_tem % kfl_time_order == 2 ) then
                 prthe(1) = 2.0_rp * prthe(1) - prthe(2)
              end if
           end if

           !
           ! Save previous time steps
           !
           nt_save = 3
           if( ADR_tem % kfl_time_scheme == TRAPEZOIDAL .or. ADR_tem % kfl_time_scheme == BDF) then
              nt_save = ADR_tem % ntime
           end if

           do itime = nt_save,2,-1
              prthe(itime) = prthe(itime-1)
           end do

           mass_curr = mass_curr + dmass

     end select
     
  end subroutine therm_press_update


  subroutine therm_press_volume_integral(vint)
     real(rp), intent(out) :: vint

     real(rp)              :: elcod(ndime,mnode)       
     real(rp)              :: eltem(mnode)
     real(rp)              :: gpvol(mgaus) 
     real(rp)              :: gpcar(ndime,mnode,mgaus) 
     real(rp)              :: gphes(ntens,mnode) 
     real(rp)              :: gptem(mgaus)
     real(rp)              :: gpwme(mgaus)
     logical(lg)           :: kfl_tesgs

     !
     ! Initialize
     !
     vint = 0.0_rp
     if (associated(ADR_tem % sgs)) then 
        kfl_tesgs = .true.
     else
        kfl_tesgs = .false.
     end if

     elements: do ielem = 1,nelem
        pelty = ltype(ielem) 
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)

        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           eltem(inode)=tempe(ipoin,1)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do
        call elmcar(&
             pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
             gphes,ielem)

        if (kfl_lookg_tem == 0_ip) then 
           !
           ! Gather properties
           !
           do igaus = 1,pgaus
              gptem(igaus) = 0.0_rp 
              gpwme(igaus) = 1.0_rp 
              do inode = 1,pnode
                 gptem(igaus) = gptem(igaus) + elmar(pelty) % shape(inode,igaus)*eltem(inode)
              end do
           enddo
           
        else
           !
           ! Properties already on Gauss points
           ! 
           do igaus = 1,pgaus
              gptem(igaus) = tempe_gp(ielem) % a(igaus,1,1)
              gpwme(igaus) = wmean_gp(ielem) % a(igaus,1,1)
           enddo
        endif

        !
        ! Add SGS contribution to temperature
        !
        do igaus = 1,pgaus
           if( kfl_tesgs ) gptem(igaus) = gptem(igaus) + ADR_tem % sgs(ielem)%a(1,igaus,1)
        enddo

        !
        ! Integrate properties
        !
        do igaus = 1,pgaus
           vint = vint +  gpvol(igaus) * gpwme(igaus) / gptem(igaus)
        end do
 
     end do elements

     !
     ! Sum over MPIs 
     !
     call PAR_SUM(vint,'IN MY ZONE')
               
  end subroutine therm_press_volume_integral


  subroutine therm_press_mass_change(dm)
     real(rp), intent(out) :: dm
     real(rp)              :: gbvel(ndime)
     real(rp)              :: baloc(ndime,ndime)
     real(rp)              :: elcod(ndime,mnode)       
     real(rp)              :: bocod(ndime,mnodb)       
     real(rp)              :: bovel(ndime,mnodb)       
     real(rp)              :: gbden(mgaub) 

     !real(rp)              :: elmassk(mnode)
     !real(rp)              :: gpmassk(mgaus)
     !real(rp)              :: gpvol(mgaus) 
     !real(rp)              :: gpcar(ndime,mnode,mgaus) 
     !real(rp)              :: gphes(ntens,mnode) 
     real(rp)              :: multiplicity_loc
     

     real(rp)              :: shap,eucta,gbsur
     integer(ip)           :: dummi,jj 

     !
     ! Initialize
     !
     dm = 0.0_rp

     !
     ! Loop over boundaries
     !
     boundaries: do iboun=1,nboun
        ! boundary 
        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)
        pgaub = ngaus(pblty)
        ielem = lelbo(iboun)
        
        ! volume
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
           
        !
        ! Gather operations
        !
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           do idime=1,ndime
              bovel(idime,inodb) = veloc(idime,ipoin,1)
              bocod(idime,inodb) = coord(idime,ipoin)
           end do
        end do
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do
           
        call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)

        gauss_points: do igaub=1,pgaub
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                bocod,baloc,eucta)                                   ! and Jacobian
           gbsur=elmar(pblty)%weigp(igaub)*eucta 
           call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
           do idime=1,ndime
              gbvel(idime)=0.0_rp
           end do

           do inodb=1,pnodb
              shap = elmar(pblty)%shape(inodb,igaub) 
              do idime=1,ndime
                 gbvel(idime)=gbvel(idime) + shap * bovel(idime,inodb)                       
              end do
           enddo

           
           do idime=1,ndime
              dm = dm + gbvel(idime)*baloc(idime,ndime)*gbsur*gbden(igaub)
           enddo 

        enddo gauss_points 
     enddo boundaries

     !
     ! Contribution of volumetric mass source (from particles)
     !
     if (associated(mass_sink)) then
        jj    = 0 
        do ipoin = 1,npoin
           if (ipoin <= npoi1) then
              multiplicity_loc = 1.0_rp
           else
              jj = jj + 1
              multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
           endif
           dm = dm + mass_sink(ipoin)/multiplicity_loc / ADR_tem % dtinv
        enddo   
     endif
     !elements: do ielem = 1,nelem
     !   pelty = ltype(ielem) 
     !   pnode = nnode(pelty)
     !   pgaus = ngaus(pelty)

     !   elmassk(:) = 0.0_rp
     !   if( associated(mass_sink) ) then
     !      do inode = 1,pnode
     !         ipoin = lnods(inode,ielem)
     !         elmassk(inode) = mass_sink(ipoin)
     !      enddo
     !   endif
     !   do inode=1,pnode
     !      ipoin=lnods(inode,ielem)
     !      do idime=1,ndime
     !         elcod(idime,inode)=coord(idime,ipoin)
     !      end do
     !   end do
     !   call elmcar(&
     !        pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
     !        elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
     !        gphes,ielem)

     !   !
     !   ! Gather properties
     !   !
     !   do igaus = 1,pgaus
     !      gpmassk(igaus) = 0.0_rp 
     !      do inode = 1,pnode
     !         gpmassk(igaus) = gpmassk(igaus) + elmar(pelty) % shape(inode,igaus)*elmassk(inode)
     !      end do
     !   enddo

     !   !
     !   ! Integrate properties
     !   !
     !   do igaus = 1,pgaus
     !      dm = dm +  gpvol(igaus) * gpmassk(igaus) / ADR_tem % dtinv
     !   end do
 
     !end do elements

     !
     ! Sum over MPIs 
     !
     call PAR_SUM(dm,'IN MY ZONE')
               
  end subroutine therm_press_mass_change


end module mod_tem_therm_press

