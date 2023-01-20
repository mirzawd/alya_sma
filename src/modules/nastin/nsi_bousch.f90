!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bousch(&
     pnode,pnodb,pgaub,pgaus,lnods,lboel,gpcar,gbcar,shaga,&
     baloc,gbsur,gbden,gbvis,gbpor,gbsha,bovel,chale,elmap)

  !----------------------------------------------------------------------
  !****f* Nastin/nsi_bousch
  ! NAME 
  !    nsi_bousch
  ! DESCRIPTION
  !    Compute the Schur complement preconditioner
  ! USES 
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_nastin, only     :  nodpr_nsi,lperp_nsi,kfl_perip_nsi,&
       &                      kfl_confi_nsi,kfl_fixpr_nsi,&
       &                      staco_nsi,dtinv_nsi,kfl_taush_nsi,&
       &                      kfl_predi_nsi,kfl_advec_nsi,kfl_matdi_nsi,&
       &                      pabdf_nsi
  use def_domain, only     :  mnode,nperi,ndime,lpoty
  use mod_tauadr, only     :  tauadr
  implicit none
  integer(ip), intent(in)  :: pnode,pnodb,pgaub,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: gbcar(ndime,pnodb,pgaub)
  real(rp),    intent(in)  :: shaga(pgaus,pnode)
  real(rp),    intent(in)  :: baloc(ndime)
  real(rp),    intent(in)  :: gbsur(pgaub)
  real(rp),    intent(in)  :: gbden(pgaub)
  real(rp),    intent(in)  :: gbvis(pgaub)
  real(rp),    intent(in)  :: gbpor(pgaub)
  real(rp),    intent(in)  :: gbsha(pnodb,pgaub)
  real(rp),    intent(in)  :: bovel(ndime,pnodb)
  real(rp),    intent(in)  :: chale(ndime)
  real(rp),    intent(out) :: elmap(pnode,pnode)
  integer(ip)              :: inode,jnode,kdime,igaub,ipoin,ibopo  
  integer(ip)              :: ipres,iperi,idime,inodb,jnodb
!  integer(ip)              :: ihang
  integer(ip)              :: igaus
  real(rp)                 :: fact1,fact2,fact3,fact4
  real(rp)                 :: tau,adv,dif,rea,gpadv(3),gpvno

  !----------------------------------------------------------------------
  !
  ! Cartesian derivates at boundary Gauss points 
  !            
  !----------------------------------------------------------------------

  do igaub = 1,pgaub
     do inode = 1,pnode                                       
        do idime = 1,ndime                                    
           gbcar(idime,inode,igaub) = 0.0_rp                        
           do inodb = 1,pnodb                                  
              do igaus = 1,pgaus
                 gbcar(idime,inode,igaub) = gbcar(idime,inode,igaub)   &
                      + shaga(igaus,lboel(inodb)) * gbsha(inodb,igaub) &
                      * gpcar(idime,inode,igaus)
              end do
           end do
        end do
     end do
  end do

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

  !----------------------------------------------------------------------
  !
  ! Elemental matrix
  !
  !----------------------------------------------------------------------

  if(kfl_predi_nsi==2) then
     !
     ! P (+ App) = ( (dt/rho)*grad p , grad q ) (+ App)
     !     
     fact2 = 1.0_rp/(dtinv_nsi*pabdf_nsi(1))
     do igaub = 1,pgaub
        fact3 = gbsur(igaub)*( fact2/gbden(igaub) )
        do inodb = 1,pnodb
           inode = lboel(inodb)
           do jnodb = 1,pnodb
              inode = lboel(inodb)
              fact1 = 0.0_rp
              do kdime = 1,ndime
                 fact1 = fact1 + gbcar(kdime,inode,igaub)&
                      &        * gbcar(kdime,jnode,igaub)
              end do
              fact1 = fact1*fact3
              elmap(inode,jnode) = elmap(inode,jnode)+fact1
              elmap(jnode,inode) = elmap(jnode,inode)+fact1
           end do
           fact1 = 0.0_rp
           do kdime = 1,ndime
              fact1 = fact1 + gbcar(kdime,inode,igaub)&
                   &        * gbcar(kdime,inode,igaub)
           end do
           fact1 = fact1*fact3
           elmap(inode,inode) = elmap(inode,inode)+fact1
        end do
     end do

  else if(kfl_predi_nsi==3) then
     !
     ! P (+ App) = ( tau'*grad p , grad q ) (+ App)
     !     
     do igaub = 1,pgaub

        if( kfl_advec_nsi == 1 ) then
           gpadv(1) = 0.0_rp
           gpadv(2) = 0.0_rp
           gpadv(3) = 0.0_rp
           do inodb = 1,pnodb
              do idime = 1,ndime
                 gpadv(idime) = gpadv(idime) &
                      + gbsha(inodb,igaub) * bovel(idime,inodb)
              end do
           end do
           call vecnor(gpadv,ndime,gpvno,2_ip)
        else
           gpvno = 0.0_rp
        end if

        adv = gpvno*gbden(igaub)              ! Convective term
        dif = gbvis(igaub)                    ! Viscous term
        rea = gbpor(igaub)                    ! Reaction term
        call tauadr(&
             kfl_taush_nsi,staco_nsi,adv,dif,rea,&
             chale(1),chale(2),tau)
        tau = 1.0_rp / ( gbden(igaub)*( dtinv_nsi * pabdf_nsi(1)) + 1.0_rp/tau )

        fact3 = tau * gbsur(igaub)
        do inodb = 1,pnodb
           inode = lboel(inodb)
           do jnode = 1,pnode
              fact1 = 0.0_rp
              do kdime = 1,ndime
                 fact1 = fact1 + gbcar(kdime,jnode,igaub) * baloc(kdime)
              end do
              elmap(inode,jnode) = elmap(inode,jnode) &
                   - fact1 * fact3 * gbsha(inodb,igaub)
           end do
        end do
     end do

  else if(kfl_predi_nsi==4) then
     !
     ! P (+ App) = mu*M (+ App)
     !     
     do igaub=1,pgaub

        fact4 = gbsur(igaub)/gbvis(igaub)

        do inode=1,pnode
           do jnode=inode+1,pnode
              fact1=0.0_rp
              do kdime=1,ndime
                 fact1=fact1+gbcar(kdime,inode,igaub)&
                      &     *gbcar(kdime,jnode,igaub)
              end do
              fact1=gbsha(inode,igaub)*gbsha(jnode,igaub)*fact4
              elmap(inode,jnode)=elmap(inode,jnode)+fact1
              elmap(jnode,inode)=elmap(jnode,inode)+fact1
           end do
           fact1=0.0_rp
           do kdime=1,ndime
              fact1=fact1+gbcar(kdime,inode,igaub)&
                   &     *gbcar(kdime,inode,igaub)
           end do
           fact1=gbsha(inode,igaub)*gbsha(inode,igaub)*fact4
           elmap(inode,inode)=elmap(inode,inode)+fact1
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Prescribe boundary conditions
  !
  !----------------------------------------------------------------------

return
  if( kfl_matdi_nsi == 0 ) then

     do inode=1,pnode
        !
        ! Neumann nodes
        !
        ipoin=lnods(inode)
        ibopo=lpoty(ipoin)
        if(ibopo/=0) then
           if ( kfl_fixpr_nsi(1,ipoin) > 0 ) then
              fact1=elmap(inode,inode)
              do jnode=1,pnode
                 elmap(jnode,inode)=0.0_rp
                 elmap(inode,jnode)=0.0_rp              
              end do
              elmap(inode,inode)=fact1
           end if
        end if
     end do

     if(kfl_confi_nsi==1) then
        !
        ! Confined flow and periodicity
        !
        if(nperi/=0) then

           do inode=1,pnode
              ipoin=lnods(inode)
              ipres=0
              if(kfl_perip_nsi==1) then
                 iperi=1
                 loop_iperi: do while(lperp_nsi(iperi)/=0)
                    if(lperp_nsi(iperi)==ipoin) then
                       ipres=1
                       exit loop_iperi
                    end if
                    iperi=iperi+1
                 end do loop_iperi
              else
                 if(ipoin==nodpr_nsi) ipres=1
              end if
              if(ipres==1) then
                 fact1=elmap(inode,inode)
                 do jnode=1,pnode
                    elmap(jnode,inode)=0.0_rp
                    elmap(inode,jnode)=0.0_rp              
                 end do
                 elmap(inode,inode)=fact1
              end if
           end do

        else

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

end subroutine nsi_bousch

