!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updsgs(&
     pgaus,pnode,ndofn,ielem,chale,elvel,gpadv,gpvis,&
     gpden,rmom1,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,&
     gpsgs,gpsgi,gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,  &
     dtsgs,resis_nsi,itsta_nsi,rmsgs_nsi,resgs_nsi,  &
     gppor)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_updsgs
  ! NAME 
  !    nsi_updsgs
  ! DESCRIPTION
  !    This subroutine updates the subgrid scale
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_nastin, only       :  staco_nsi,kfl_sgsti_nsi,&
       &                        kfl_sgsco_nsi,misgs_nsi,tosgs_nsi,&
       &                        kfl_taust_nsi,kfl_stabi_nsi,corio_nsi,&
       &                        kfl_sgsli_nsi,kfl_rmom2_nsi,relsg_nsi,&
       &                        kfl_stabi_nsi,mmsgs_nsi
  use mod_tauadr, only       :  tauadr
  implicit none
  integer(ip), intent(in)    :: pgaus,pnode,ndofn,ielem
  real(rp),    intent(in)    :: chale(2)
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gprhs(ndofn,pgaus)
  real(rp),    intent(in)    :: gprhs_sgs(ndofn,pgaus)
  real(rp),    intent(in)    :: gpgpr(ndime,pgaus)
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(inout) :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: rmom1(pnode,pgaus)
  real(rp),    intent(in)    :: rmom2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(out)   :: gpsp1(pgaus)
  real(rp),    intent(inout) :: gpsgs(ndime,pgaus,*)
  real(rp),    intent(out)   :: gpsgi(ndime,pgaus)
  real(rp),    intent(out)   :: gpgve(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpvep(ndime,pgaus)
  real(rp),    intent(in)    :: gpgrp(ndime,pgaus)
  real(rp),    intent(out)   :: gpst1(pgaus)
  real(rp),    intent(in)    :: dtsgs
  real(rp),    intent(inout) :: resis_nsi(2,*)
  integer(ip), intent(inout) :: itsta_nsi(*)
  real(rp),    intent(out)   :: rmsgs_nsi
  real(rp),    intent(out)   :: resgs_nsi(*)
  real(rp),    intent(in)    :: gppor(pgaus)
  integer(ip)                :: idime,igaus,inode,itsgs,jdime,jtsgs
  real(rp)                   :: resgs,gpnew(3),gpnor,gpnum,gpdnm,rels1
  real(rp)                   :: adv,dif,rea,h2,gpvno,dummr
  real(rp)                   :: dj_fi(ndime,ndime)       ! d(f_i(u'))/du'_j   ! for NEWTON-R (kfl_sgsli_nsi==2)
  real(rp)                   :: djinv(ndime,ndime)       ! its inverse
  real(rp)                   :: funcu(ndime)             ! f(u')
  real(rp)                   :: deltu(ndime)             ! delta u'
  real(rp)                   :: auxi1,auxi2

  if( kfl_stabi_nsi == -1 ) return
  if( kfl_sgsti_nsi == 0 .and. kfl_sgsco_nsi == 0 ) return

  h2    = chale(2) * chale(2)
  rels1 = 1.0_rp - relsg_nsi

  !----------------------------------------------------------------------
  !
  ! GPGVE: Velocity gradients
  !
  !----------------------------------------------------------------------

  do igaus = 1,pgaus
     do jdime = 1,ndime
        do idime = 1,ndime
           gpgve(idime,jdime,igaus) = 0.0_rp
           do inode = 1,pnode
              gpgve(idime,jdime,igaus) = gpgve(idime,jdime,igaus)&
                   + gpcar(idime,inode,igaus) * elvel(jdime,inode)
           end do
        end do
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! GPSGI: Never changing residual (all without PICARD advection)
  !        Ri(u) = f - rho*u/dt - sig*u - grad(p) + div[2*mu*eps(u)]
  !        where f = (body forces) + rho*u^n/dt + rho*u'^n/dt
  !        For Split OSS: add time term of the SGS
  !
  !----------------------------------------------------------------------

  do igaus = 1,pgaus
     do idime = 1,ndime
        gpsgi(idime,igaus) = gprhs(idime,igaus) + gprhs_sgs(idime, igaus) - gpgpr(idime,igaus) 

        do inode=1,pnode
           gpsgi(idime,igaus) = gpsgi(idime,igaus)&
                - rmom1(inode,igaus) * elvel(idime,inode)
        end do
     end do
  end do
  if( kfl_rmom2_nsi /= 0 ) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           do inode = 1,pnode
              do jdime = 1,ndime
                 gpsgi(idime,igaus) = gpsgi(idime,igaus)&
                      - rmom2(idime,jdime,inode,igaus) * elvel(jdime,inode)
              end do
           end do
        end do
     end do
  end if

  if( kfl_stabi_nsi == 1 ) then
     ! 
     ! OSS: rho * u'^n/dt was not included in nsi_elmre3
     ! u ' <= u ' + u'^n / dt
     !
!!$     if( kfl_sgsti_nsi /= 0 ) then
!!$        do igaus = 1,pgaus
!!$           dummr = gpden(igaus) * dtsgs
!!$           do idime = 1,ndime
!!$              gpsgi(idime,igaus) = gpsgi(idime,igaus) + dummr * gpsgs(idime,igaus,2)
!!$           end do
!!$        end do
!!$     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Solve subgrid scale equation usig PICARD (1) or NEWTON RAPHSON (2)
  !
  !----------------------------------------------------------------------

!!$  if( kfl_sgsli_nsi == 0 ) then
!!$
!!$     do igaus = 1,pgaus
!!$        !
!!$        ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*u/h + rho*|w| + sig )
!!$        !
!!$        call vecnor(gpvel(1,igaus),ndime,gpvno,2_ip)
!!$        adv = gpden(igaus)*gpvno                           ! Convective term: rho*|u|
!!$        dif = gpvis(igaus)                                 ! Viscous term:    mu
!!$        rea = gpden(igaus)*corio_nsi + abs(gppor(igaus))   ! Coriolis: w + Porosity: sig
!!$        call tauadr(&
!!$             kfl_taust_nsi,staco_nsi,adv,dif,rea,&
!!$             chale(1),chale(2),gpst1(igaus))
!!$
!!$        do idime =1, ndime
!!$           do jdime =1, ndime
!!$              gpsgi(idime, igaus) = gpsgi(idime, igaus) - gpden(igaus)&
!!$                   * gpvel(jdime,igaus) * gpgve(jdime,idime,igaus)
!!$           end do
!!$        end do
!!$        
!!$        if( kfl_sgsti_nsi == 1 ) then
!!$           !
!!$           ! tau1' = 1 / ( rho/dt + 1/tau1 )
!!$           !
!!$           gpsp1(igaus) = 1.0_rp / ( gpden(igaus) * dtsgs + 1.0_rp / gpst1(igaus) )
!!$           gpsgs(1:ndime,igaus,1) = gpsp1(igaus) * gpsgi(1:ndime,igaus) 
!!$        else
!!$           gpsgs(1:ndime,igaus,1) = gpst1(igaus) * gpsgi(1:ndime,igaus) 
!!$        end if
!!$
!!$        gpadv(1:ndime,igaus) = gpvel(1:ndime,igaus) + gpsgs(1:ndime,igaus,1)
!!$
!!$     end do
!!$
!!$  else 
     if( kfl_sgsli_nsi == 1 ) then 

     do igaus = 1,pgaus
        itsgs = 0
        resgs = 1.e06_rp

        do while( itsgs < misgs_nsi .and. resgs > tosgs_nsi )
           itsgs = itsgs + 1
           jtsgs = min(mmsgs_nsi,itsgs)
           !
           ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
           ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]
           !
           call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
           adv = gpden(igaus)*gpvno                           ! Convective term: rho*|u+u'|
           dif = gpvis(igaus)                                 ! Viscous term:    mu
           rea = gpden(igaus)*corio_nsi + abs(gppor(igaus))   ! Coriolis: w + Porosity: sig
           call tauadr(&
                kfl_taust_nsi,staco_nsi,adv,dif,rea,&
                chale(1),chale(2),gpst1(igaus))
           !
           ! tau1' = 1 / ( rho/dt + 1/tau1 )
           !
           gpsp1(igaus) = 1.0_rp / ( gpden(igaus) * dtsgs + 1.0_rp / gpst1(igaus) )
           !
           ! GPNEW: new subgrid scale 
           ! 0. ASGS ........ u' = tau1' * [ Ri(u) - rho*(uc.grad)u ]
           ! 1. Full OSS .... u' = tau1' * [ Ri(u) - rho*(uc.grad)u ] + P
           ! 2. Split OSS ... u' = tau1' * [ - rho*(uc.grad)u  + rho*u'/dt ] + P
           !
           do idime = 1,ndime
              gpnew(idime) = gpsgi(idime,igaus)
              do jdime = 1,ndime
                 gpnew(idime) = gpnew(idime) - gpden(igaus)&
                      * gpadv(jdime,igaus) * gpgve(jdime,idime,igaus)
              end do
              gpnew(idime) = gpsp1(igaus) * gpnew(idime) 
           end do
           if( kfl_stabi_nsi == 2 ) then
              call runend('NSI_UPDSGS: NOT CODED')
!!$           else if( kfl_stabi_nsi /= 0 ) then
!!$              do idime = 1,ndime
!!$                 gpnew(idime) = gpnew(idime) + gpsp1(igaus) * gpvep(idime,igaus) / gpst1(igaus)
!!$              end do
           end if
           !
           ! RESGS: Residual and update subgrid scale
           ! GPSGS: New subgrid scale=GPNEW
           !
           gpnum = 0.0_rp
           gpdnm = 0.0_rp
           do idime = 1,ndime
              gpnew(idime)         = relsg_nsi * gpnew(idime) + rels1 * gpsgs(idime,igaus,1)
              dummr                = gpsgs(idime,igaus,1) - gpnew(idime)
              gpnum                = gpnum + dummr * dummr
              gpdnm                = gpdnm + gpnew(idime) * gpnew(idime)
              gpsgs(idime,igaus,1) = gpnew(idime)
           end do
           if( gpdnm /= 0.0_rp ) then
              resgs = sqrt( gpnum/gpdnm )
           else
              resgs = sqrt( gpnum )
           end if

           if( kfl_sgsco_nsi >= 1 ) then           
              !
              ! GPADV: Update advection uc = u + u'
              !
              do idime = 1,ndime
                 gpadv(idime,igaus) = gpvel(idime,igaus) + gpsgs(idime,igaus,1)
              end do
              !
              ! RESIS_NSI: Inner residual
              !
              call vecnor(gpnew,ndime,gpnor,2_ip)
              resis_nsi(1,jtsgs) = resis_nsi(1,jtsgs) + gpnum * gpnum
              resis_nsi(2,jtsgs) = resis_nsi(2,jtsgs) + gpnor * gpnor
           end if

        end do

        rmsgs_nsi = max(rmsgs_nsi,resgs)
        !
        ! ITSTA_NSI: Subgrid scale statistics
        !
        if( kfl_sgsco_nsi >= 1 ) itsta_nsi(jtsgs) = itsta_nsi(jtsgs) + 1
        
     end do

  else ! Newton Raphson

     do igaus = 1,pgaus
        itsgs = 0
        resgs = 1.e06_rp

        do while( itsgs < misgs_nsi .and. resgs > tosgs_nsi )
           itsgs = itsgs + 1
           jtsgs = min(mmsgs_nsi,itsgs)
           !
           ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
           ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]
           !
           call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
           adv = gpden(igaus)*gpvno                           ! Convective term: rho*|u+u'|
           dif = gpvis(igaus)                                 ! Viscous term:    mu
           rea = abs(gppor(igaus))                            ! Porosity term:   sig
           call tauadr(&
                kfl_taust_nsi,staco_nsi,adv,dif,rea,&
                chale(1),chale(2),gpst1(igaus))
           !
           ! tau1' = 1 / ( rho/dt + 1/tau1 )
           !
           gpsp1(igaus) = 1.0_rp / ( gpden(igaus) * dtsgs + 1.0_rp / gpst1(igaus) )
           !
           ! obtain d(f_i(u'))/du'_j and its inverse
           !
           if ( gpvno>1.0d-8 ) then 
              auxi1 = 2.0_rp * gpden(igaus) / ( chale(1) * gpvno )     ! C2 * rho / ( h * | u + u'| ) 
           else
              auxi1 =  0.0_rp
           end if
           do idime=1,ndime
              do jdime=1,ndime
                 dj_fi(idime,jdime) = &
                      + auxi1 * gpadv(jdime,igaus) * gpsgs(idime,igaus,1) &  ! auxi1 * ( u_jd + u'_jd ) * u'_id 
                      + gpden(igaus) * gpgve(jdime,idime,igaus)              ! rho * du_id / dx_jd
              end do
           end do

           auxi1 = 1.0_rp / gpsp1(igaus) 
           do idime=1,ndime
              dj_fi(idime,idime) = dj_fi(idime,idime) + auxi1
           end do
           call invmtx(dj_fi,djinv,auxi2,ndime)
           !
           ! obtain f(u')
           !
           do idime=1,ndime
              funcu(idime) = - gpsgi(idime,igaus)  &  ! - non changing part of the residual - rho u'_n /(tita*dt) -- See before
                   + auxi1 * gpsgs(idime,igaus,1)     ! + ( inv(tau1') * u'
              do jdime=1,ndime
                 funcu(idime) = funcu(idime) + gpden(igaus) * gpadv(jdime,igaus) * gpgve(jdime,idime,igaus)    ! rho * ( u_jd + u'_jd ) * dv_i / dx_j (changing part of residual)
              end do
           end do
           if( kfl_stabi_nsi == 2 ) then
              call runend('NSI_UPDSGS: NOT CODED')
!!$           else if( kfl_stabi_nsi /= 0 ) then
!!$              do idime = 1,ndime
!!$                 funcu(idime) = funcu(idime) - gpvep(idime,igaus) / gpst1(igaus)
!!$              end do
           end if
           !
           ! from A_ij * deltau'_j =  - f(u')_i ; (where A_ij = d(f_i(u'))/du'_j )
           ! obtain deltau'_j = - invA_ji *  f(u')_i 
           !
           do jdime=1,ndime
              deltu(jdime) = 0.0_rp
              do idime=1,ndime
                 deltu(jdime) = deltu(jdime) - djinv(jdime,idime) * funcu(idime)
              end do
           end do
           !
           ! RESGS: Residual and update subgrid scale
           ! GPSGS: New subgrid scale=GPNEW
           !
           gpnum = 0.0_rp
           gpdnm = 0.0_rp
           do idime = 1,ndime
              gpnew(idime)         = gpsgs(idime,igaus,1) + relsg_nsi * deltu(idime)
              dummr                = deltu(idime)
              gpnum                = gpnum + dummr * dummr
              gpdnm                = gpdnm + gpnew(idime) * gpnew(idime)
              gpsgs(idime,igaus,1) = gpnew(idime)
           end do
           if( gpdnm /= 0.0_rp ) then
              resgs = sqrt( gpnum/gpdnm )
           else
              resgs = sqrt( gpnum )
           end if

           if( kfl_sgsco_nsi >= 1 ) then           
              !
              ! GPADV: Update advection uc = u + u'
              !
              do idime = 1,ndime
                 gpadv(idime,igaus) = gpvel(idime,igaus) + gpsgs(idime,igaus,1)
              end do
              !
              ! RESIS_NSI: Inner residual
              !
              call vecnor(gpnew,ndime,gpnor,2_ip)
              resis_nsi(1,jtsgs) = resis_nsi(1,jtsgs) + gpnum * gpnum
              resis_nsi(2,jtsgs) = resis_nsi(2,jtsgs) + gpnor * gpnor
           end if

        end do

        rmsgs_nsi = max(rmsgs_nsi,resgs)
        !
        ! ITSTA_NSI: Subgrid scale statistics
        !
        if( kfl_sgsco_nsi >= 1 ) itsta_nsi(jtsgs) = itsta_nsi(jtsgs) + 1
        
     end do

  end if

  !-------------------------------------------------------------------
  !
  ! Subgrid scale residual and update VESGS
  !
  !-------------------------------------------------------------------

  !call nsi_sgsope(3_ip,ielem,pgaus,gpsgs,dummr,resgs_nsi)

end subroutine nsi_updsgs
