!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updthe(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_updthe
  ! NAME 
  !    nsi_updthe
  ! DESCRIPTION
  !    This subroutine updates the thermodynamic pressure:
  !      V   p^{n+theta)-p^n                    gamma-1
  !    ----- --------------- = -p^(n+theta) M + ------- Q
  !    gamma     theta*dt                        gamma
  !    +-                  -+
  !    |        V           |                     V              gamma-1
  !    | -------------- + M | p^{n+theta) = -------------- p^n + ------- Q
  !    | gamma*theta*dt    -+               gamma*theta*dt        gamma
  !    +- 
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_nastin
  use def_domain
  use def_solver
  use def_kermod, only       :  gasco 
  use mod_communications_global, only : PAR_SUM
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: eltem(mnode),gptem,gpvol(mgaus)
  real(rp)                :: gphes(ntens,mnode)
  real(rp)                :: gpcar(ndime,mnode,mgaus), fact0, bmflx
  integer(ip)             :: ielem,inode,ipoin,idime, itime, jtime
  integer(ip)             :: pnode,pgaus,igaus
  integer(ip)             :: pelty,pmate,pevat
  real(rp)                :: avtem_xxx,volum_xxx

  if(kfl_regim_nsi==3) then
     if(abs(itask)==1) then
        !
        ! Integrate thermodynamic pressure in time, used for closed systems with inflow/outflow
        !
        if(kfl_prthe_nsi==2) then
 
           !
           ! inverse of temp and surface flux integration 
           !
           bmflx=0.0_rp
           vinvt_nsi(1)=0.0_rp
           if(INOTMASTER) then
              call nsi_smflow(bmflx)
              call nsi_vinvte(vinvt_nsi(1))
           end if        
      
           !
           ! Parall: Sum over all subdomains
           !
           call PAR_SUM(bmflx)
           call PAR_SUM(vinvt_nsi(1))
           
           if(itask==-1) then ! if initialization task
              ! if initial mass was set, then
              if (imass.gt.zeror)   then  ! set by archive
                 vinvt_nsi(2) = imass*gasco/prthe(4)*vodom           
              else 
!!$
!!$             initial mass was not fixed,                  
!!$
!!$                 call runend('NSI_UPDTHE:PLEASE, SET THE INITIAL MASS')      
                 vinvt_nsi(2)=vinvt_nsi(1) 
              end if              
              pabdf_nsi(1)=1.0_rp
              pabdf_nsi(2)=-1.0_rp
              nbdfp_nsi=2
           end if
           !
           ! Update thermodynamic pressure
           !
           fact0 = 0.0_rp
           do itime =2, nbdfp_nsi
              fact0 = fact0 -  pabdf_nsi(itime)*prthe(itime)*vinvt_nsi(itime)
           end do
           prthe(1)= fact0/(vinvt_nsi(1)*pabdf_nsi(1)+bmflx/dtinv_nsi)
           !
           ! Update thermodynamic pressure time derivative
           !
           dpthe = 0.0_rp
           do itime =1, nbdfp_nsi
              dpthe = dpthe + pabdf_nsi(itime)*prthe(itime)
           end do
           dpthe =  dpthe*dtinv_nsi
           xmass_nsi = prthe(1)*vinvt_nsi(1)/gasco



        else if(kfl_prthe_nsi==1) then   ! Thermodynamic pressure from mass conservation, 
                                         !used in closed systems without flow
           !
           ! inverse of temperature volume integration
           !
           vinvt_nsi(1)=0.0_rp           
           if(INOTMASTER) call nsi_vinvte(vinvt_nsi(1))
           call PAR_SUM(vinvt_nsi(1))
           
           if(itask==-1) then  ! initialization task
              if (imass.gt.zeror)   then   ! if initial mass was set in datafile
                 vinvt_nsi(2) = imass*gasco/prthe(4)*vodom 
              else
                 vinvt_nsi(2)=vinvt_nsi(1)
              end if
              do itime =1, nbdfp_nsi
                 pabdf_nsi(itime) = 0.0_rp ! to set initial dpthe =0 
              end do
           end if
           !
           ! Update thermodynamic pressure and its time derivative
           !
           prthe(1)  = prthe(4)*vinvt_nsi(2)/vinvt_nsi(1)
           
           dpthe = 0.0_rp
           do itime =1, nbdfp_nsi
              dpthe = dpthe + pabdf_nsi(itime)*prthe(itime)
           end do
           dpthe = dpthe*dtinv_nsi
           ! total  mass prhte/R* integr(1/T dvolu)
           xmass_nsi = prthe(1)*vinvt_nsi(1)/gasco


        end if

     else if(itask==2) then

        !
        ! Crank-Nicolson
        !
        do jtime=2, kfl_tiaor_nsi +1
           itime= kfl_tiaor_nsi +3 -jtime
           if (kfl_prthe_nsi==2) vinvt_nsi(itime) = vinvt_nsi(itime-1)
        end do
       
     end if


     ! **********************************************************************************
     ! **********************************************************************************
     !
     ! Compute reference thermodynamic pressure and temperature for hydrodynamic pressure
     ! This should run only on the first time step
     !
     if (lowtr_nsi == 0.0_rp) then
        avtem_xxx=0.0_rp
        volum_xxx=0.0_rp
        if(INOTMASTER) then
            do ielem=1,nelem
              pelty = ltype(ielem)                  ! Element properties and dimensions
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              pevat = solve(ivari_nsi)%ndofn*pnode
              pmate = 1
              if(nmate>1) then ! Check if element is a solid
                 pmate=lmate(ielem)
              end if
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 eltem(inode)=tempe(ipoin,1)
                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do
              if(pmate/=-1) then
                 call elmcar(&
                      pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                      elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
                      gphes,ielem)
                 do igaus=1,pgaus
                    gptem=0.0_rp
                    do inode=1,pnode
                       gptem=gptem+elmar(pelty)%shape(inode,igaus)*eltem(inode)
                    end do
                    avtem_xxx=avtem_xxx+gpvol(igaus)*gptem
                    volum_xxx=volum_xxx+gpvol(igaus)
                 end do
              end if
           end do 
        end if
        !
        ! Parall: Sum over all subdomains
        !
        call PAR_SUM(avtem_xxx)
        call PAR_SUM(volum_xxx)
        !
        ! Store average temperature
        !
        if (volum_xxx  /= 0.0_rp ) then 
           lowtr_nsi = avtem_xxx / volum_xxx
        else
           print *,'ZERO VOLUME !!!'
        endif
     endif 
     
  end if

end subroutine nsi_updthe


subroutine nsi_vinvte(vinvt)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_vinvt
  ! NAME 
  !    nsi_vinvt
  ! DESCRIPTION
  !    This subroutine updates the thermodynamic pressure:
  !      V   p^{n+theta)-p^n                    gamma-1
  !    ----- --------------- = -p^(n+theta) M + ------- Q
  !    gamma     theta*dt                        gamma
  !    +-                  -+
  !    |        V           |                     V              gamma-1
  !    | -------------- + M | p^{n+theta) = -------------- p^n + ------- Q
  !    | gamma*theta*dt    -+               gamma*theta*dt        gamma
  !    +- 
  ! USES
  ! USED BY
  !    nsi_updthe
  !***
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_nastin
  use def_domain
  use def_solver
  implicit none
  real(rp), intent(out)   :: vinvt
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: eltem(mnode),gptem,gpvol(mgaus)
  real(rp)                :: gphes(ntens,mnode)
  real(rp)                :: gpcar(ndime,mnode,mgaus)
  integer(ip)             :: ielem,inode,ipoin,idime
  integer(ip)             :: pnode,pgaus,igaus
  integer(ip)             :: pelty,pmate,pevat
  
  !
  ! Initializations
  !
  vinvt=0.0_rp
  !
  ! Loop on elements
  !


     elements: do ielem = 1,nelem
        pelty = ltype(ielem)                  ! Element properties and dimensions
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        pevat = solve(ivari_nsi)%ndofn*pnode
        pmate = 1
        if(nmate>1) then ! Check if element is a solid
           pmate=lmate(ielem)
        end if
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           eltem(inode)=tempe(ipoin,1)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
        end do
        if(pmate/=-1) then
           call elmcar(&
                pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
                gphes,ielem)
           do igaus=1,pgaus
              gptem=0.0_rp
              do inode=1,pnode
                 gptem=gptem+elmar(pelty)%shape(inode,igaus)*eltem(inode)
              end do
              vinvt=vinvt+gpvol(igaus)/gptem
           end do
        end if
     end do elements

end subroutine nsi_vinvte


subroutine nsi_smflow(bmflx)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_vinvt
  ! NAME 
  !    nsi_vinvt
  ! DESCRIPTION
  !    This subroutine updates the thermodynamic pressure:
  !      V   p^{n+theta)-p^n                    gamma-1
  !    ----- --------------- = -p^(n+theta) M + ------- Q
  !    gamma     theta*dt                        gamma
  !    +-                  -+
  !    |        V           |                     V              gamma-1
  !    | -------------- + M | p^{n+theta) = -------------- p^n + ------- Q
  !    | gamma*theta*dt    -+               gamma*theta*dt        gamma
  !    +- 
  ! USES
  ! USED BY
  !    nsi_updthe
  !***
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_nastin
  use def_domain
  use def_solver
  use mod_bouder
  implicit none
  real(rp), intent(out)   :: bmflx
  real(rp)                :: baloc(ndime,ndime),bovel(ndime,mnodb) , botem(mnodb)
  real(rp)                :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                :: gbsur,eucta, gbvel(ndime), gbtem
  real(rp)                :: fact0
  integer(ip)             :: ielem,inode,ipoin,idime
  integer(ip)             :: pnode,pgaus,iboun,igaub,inodb
  integer(ip)             :: pelty,pmate,pblty,pnodb,pgaub,pevat
  !
  ! Initializations
  !
  bmflx=0.0_rp
  !
  ! Loop on elements
  ! 

     !
     ! Mass M= int_S u.n ds
     !
     do iboun=1,nboun
        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaub = ngaus(pblty)
        pgaus = ngaus(pelty)
        pevat = solve(ivari_nsi)%ndofn*pnode
        pmate = 1
        if(nmate>1) then
           pmate=lmate(ielem)
        end if
        if(pmate/=-1) then
           do inodb=1,pnodb
              ipoin=lnodb(inodb,iboun)
              do idime=1,ndime
                 botem(inodb)       = tempe(ipoin,1) 
                 bovel(idime,inodb) = veloc(idime,ipoin,1)
                 bocod(idime,inodb) = coord(idime,ipoin)
              end do
           end do
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           do igaub=1,pgaub
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                   bocod,baloc,eucta)                                   ! and Jacobian
              gbsur=elmar(pblty)%weigp(igaub)*eucta 
              call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
              do idime=1,ndime
                 gbvel(idime)=0.0_rp
              end do
              gbtem =0.0_rp
              do inodb=1,pnodb
                 fact0 = elmar(pblty)%shape(inodb,igaub) 
                 do idime=1,ndime
                    gbvel(idime)=gbvel(idime)&
                         +bovel(idime,inodb)*fact0                           
                 end do
                 gbtem = gbtem + botem(inodb)*fact0
              end do
              do idime=1,ndime
                 bmflx = bmflx +gbvel(idime)*baloc(idime,ndime)*gbsur/gbtem
              end do
           end do
        end if
     end do


end subroutine nsi_smflow

