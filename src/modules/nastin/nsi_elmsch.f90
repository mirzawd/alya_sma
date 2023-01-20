!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmsch(&
     pnode,pgaus,lnods,gpcar,gpvol,gpden,gpvis,&
     gppor,gpsha,elvel,chale,gpsp1,elmap,dtinv_loc)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_elmsch
  ! NAME 
  !    nsi_elmsch
  ! DESCRIPTION
  !    Compute the Schur complement preconditioner
  ! USES 
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_nastin, only     :  nodpr_nsi,&
       &                      kfl_confi_nsi,kfl_fixpr_nsi,&
       &                      staco_nsi,kfl_taush_nsi,&
       &                      kfl_predi_nsi,kfl_advec_nsi,kfl_matdi_nsi,&
       &                      pabdf_nsi, kfl_regim_nsi,kfl_confi_nsi,&
       &                      penal_nsi
  use def_domain, only     :  mnode,ndime,lpoty
  use mod_tauadr, only     :  tauadr
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus),gpvol(pgaus)
  real(rp),    intent(in)  :: gpden(pgaus),gpvis(pgaus),gppor(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),elvel(ndime,mnode)
  real(rp),    intent(in)  :: chale(ndime),gpsp1(pgaus)
  real(rp),    intent(out) :: elmap(pnode,pnode)
  real(rp),    intent(in)  :: dtinv_loc
  integer(ip)              :: inode,jnode,kdime,igaus,ipoin,ibopo
  integer(ip)              :: idime
!  integer(ip)              :: idime,ihang
  real(rp)                 :: fact1,fact2,fact3,fact4,fact0
  real(rp)                 :: tau,adv,dif,rea,gpadv(3),gpvno, penal(pgaus)

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  do jnode = 1,pnode
     do inode = 1,pnode
        elmap(inode,jnode) = 0.0_rp
     end do
  end do

  do igaus = 1,pgaus
     penal(igaus) = 0.0_rp
  end do   

!!DMM-LM #ifdef matiaslma
  if (kfl_regim_nsi==3 .and. kfl_confi_nsi == 1) then 
     do igaus =1, pgaus
        penal(igaus) = 1.0e-4_rp*gpden(igaus)/gpvis(igaus)
     end do
  end if
!!DMM-LM #endif

  !----------------------------------------------------------------------
  !
  ! Elemental matrix
  !
  !----------------------------------------------------------------------

  if( kfl_predi_nsi == 2 ) then
     !
     ! P (+ App) = ( (dt/ (rho*pabdf_nsi(1)) )*grad p , grad q ) (+ App)
     !
     if( dtinv_loc /= 0.0_rp ) then
        fact2 = 1.0_rp / (dtinv_loc*pabdf_nsi(1))
     else
        fact2 = 0.0_rp
     end if
     do igaus = 1,pgaus
        fact3 = gpvol(igaus) * ( fact2/gpden(igaus) + gpsp1(igaus) )
        do inode = 1,pnode
           fact0 =  gpvol(igaus) * penal(igaus) * gpsha(inode,igaus) 
           do jnode = inode+1,pnode
              fact1 = 0.0_rp
              do kdime = 1,ndime
                 fact1 = fact1 + gpcar(kdime,inode,igaus)&
                      &        * gpcar(kdime,jnode,igaus)
              end do
              fact1 = fact1*fact3
              elmap(inode,jnode) = elmap(inode,jnode)+fact1 + fact0*gpsha(jnode,igaus)
              elmap(jnode,inode) = elmap(jnode,inode)+fact1 + fact0*gpsha(jnode,igaus)
           end do
           fact1 = 0.0_rp
           do kdime = 1,ndime
              fact1 = fact1 + gpcar(kdime,inode,igaus)&
                   &        * gpcar(kdime,inode,igaus)
           end do
           fact1 = fact1*fact3
           elmap(inode,inode) = elmap(inode,inode) + fact1 + fact0 * gpsha(inode,igaus)
           elmap(inode,inode) = elmap(inode,inode) + penal_nsi * gpvol(igaus) * gpsha(inode,igaus)
        end do
     end do

  else if( kfl_predi_nsi == 3 ) then
     !
     ! P (+ App) = ( tau'*grad p , grad q ) (+ App)
     !     
     do igaus = 1,pgaus

        gpadv = 0.0_rp
        if( kfl_advec_nsi /= 0 ) then
           do inode = 1,pnode
              gpadv(1:ndime) = gpadv(1:ndime) + gpsha(inode,igaus) * elvel(1:ndime,inode)
           end do
        end if
        call vecnor(gpadv,ndime,gpvno,2_ip)
        !
        !             1
        ! M^-1 = -----------
        !        rho     1
        !        --- +  ---
        !         dt    tau
        !
        adv = gpvno * gpden(igaus)            ! Convective term
        dif = gpvis(igaus)                    ! Viscous term
        rea = gppor(igaus)                    ! Reaction term
        call tauadr(&
             kfl_taush_nsi,staco_nsi,adv,dif,rea,&
             chale(1),chale(2),tau)
        fact2 = 1.0_rp / ( gpden(igaus) * dtinv_loc * pabdf_nsi(1) + 1.0_rp / tau )
        !
        ! P + App <= tau + M^-1
        !
        fact3 = ( fact2 + gpsp1(igaus) ) * gpvol(igaus) !!!FER * gpden(igaus)

        do inode = 1,pnode
           fact0 = gpvol(igaus) * penal(igaus) * gpsha(inode,igaus) 
           do jnode = inode+1,pnode
              fact1 = 0.0_rp
              do kdime = 1,ndime
                 fact1 = fact1 + gpcar(kdime,inode,igaus)&
                      &        * gpcar(kdime,jnode,igaus)
              end do
              fact1 = fact1 * fact3
              elmap(inode,jnode) = elmap(inode,jnode) + fact1  + fact0 * gpsha(jnode,igaus)
              elmap(jnode,inode) = elmap(jnode,inode) + fact1  + fact0 * gpsha(jnode,igaus)
           end do
           fact1 = 0.0_rp
           do kdime = 1,ndime
              fact1 = fact1 + gpcar(kdime,inode,igaus) * gpcar(kdime,inode,igaus)
           end do
           fact1 = fact1 * fact3
           elmap(inode,inode) = elmap(inode,inode) + fact1 + fact0 * gpsha(inode,igaus)
           elmap(inode,inode) = elmap(inode,inode) + penal_nsi * gpvol(igaus) * gpsha(inode,igaus)
        end do
     end do

  else if( kfl_predi_nsi == 4 ) then
     !
     ! P (+ App) = mu*M (+ App)
     !     
     do igaus = 1,pgaus

        fact3 = gpsp1(igaus)*gpvol(igaus) 
        fact4 = gpvol(igaus)/gpvis(igaus)

        do inode=1,pnode
           fact0 = gpvol(igaus) * ( penal(igaus)*gpsha(inode,igaus) + penal_nsi )
           do jnode=inode+1,pnode
              fact1=0.0_rp
              do kdime=1,ndime
                 fact1=fact1+gpcar(kdime,inode,igaus)&
                      &     *gpcar(kdime,jnode,igaus)
              end do
              fact1=fact1*fact3+gpsha(inode,igaus)*gpsha(jnode,igaus)*fact4
              elmap(inode,jnode)=elmap(inode,jnode)+fact1+ fact0*gpsha(jnode,igaus)
              elmap(jnode,inode)=elmap(jnode,inode)+fact1+ fact0*gpsha(jnode,igaus)
           end do
           fact1=0.0_rp
           do kdime=1,ndime
              fact1=fact1+gpcar(kdime,inode,igaus)&
                   &     *gpcar(kdime,inode,igaus)
           end do
           fact1=fact1*fact3+gpsha(inode,igaus)*gpsha(inode,igaus)*fact4
           elmap(inode,inode)=elmap(inode,inode)+fact1+ fact0*gpsha(inode,igaus)
        end do
     end do

    else if( kfl_predi_nsi == 7 ) then
       !
       ! Pure Laplacian (to be scaled later on)
       !
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                do idime = 1,ndime
                   elmap(inode,jnode) = elmap(inode,jnode) &
                        + gpvol(igaus) * gpcar(idime,inode,igaus) * gpcar(idime,jnode,igaus) 
                end do
             end do
          end do
       end do


  end if

  !----------------------------------------------------------------------
  !
  ! Prescribe boundary conditions
  !
  !----------------------------------------------------------------------

  if( kfl_matdi_nsi == 0 ) then

     do inode=1,pnode
        !
        ! Neumann nodes
        !
        ipoin = lnods(inode)
        ibopo = lpoty(ipoin)
        if( ibopo /= 0 ) then
           if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
              fact1 = elmap(inode,inode)
              do jnode = 1,pnode
                 elmap(jnode,inode) = 0.0_rp    ! Beware perhaps this line should be commented out - at least that is what 
                                                ! I do for Takuji. This leads to a non sym matrix but at least values are 
                                                ! not wiped out without sending them tp RHS - THIS is something we must anlyse
                                                ! further with guillaume
                 elmap(inode,jnode) = 0.0_rp              
              end do
              elmap(inode,inode) = fact1
           end if
        end if
     end do

     if(kfl_confi_nsi==1) then
        !
        ! Confined flow and periodicity
        !
        do inode=1,pnode
           ipoin=lnods(inode)
           if(ipoin==nodpr_nsi) then
              fact1=elmap(inode,inode)
              do jnode=1,pnode
                 elmap(jnode,inode)=0.0_rp
                 elmap(inode,jnode)=0.0_rp 
              end do
              elmap(inode,inode)=fact1
           end if
        end do

     end if

     !if(nhang<0) then
     !   !
     !   ! Hanging nodes
     !   !
     !   do inode=1,pnode
     !      ipoin=lnods(inode)
     !      ihang=0
     !      do while(ihang<nhang)
     !         ihang=ihang+1
     !         if(lhang(1,ihang)==ipoin) then
     !            ihang=nhang
     !            fact1=elmap(inode,inode)
     !            do jnode=1,pnode
     !               elmap(jnode,inode)=0.0_rp
     !               elmap(inode,jnode)=0.0_rp 
     !            end do
     !            elmap(inode,inode)=fact1              
     !         end if
     !      end do
     !   end do
     !end if

  end if

end subroutine nsi_elmsch

