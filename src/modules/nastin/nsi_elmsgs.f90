!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmsgs(&
     pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,&
     gpst2,gpsp1,gpsp2,gptt1,gptt2,rmom1,gppor,dtsgs,dtinv_loc,&
     tamin_loc,tamax_loc, elvel, gprhs, gpgpr, rmom2,gpgve)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_elmsgs
  ! NAME 
  !    nsi_elmsgs
  ! DESCRIPTION
  !    1. Compute stability parameters
  !       - tau1 [ L/(rho*U) ]
  !       - tau2 [ rho*L*U ]
  !    3. Update momentum residual
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode 
  use def_nastin, only       :  staco_nsi,kfl_sgsti_nsi,&
       &                        kfl_taust_nsi,kfl_stabi_nsi,&
       &                        kfl_ellen_nsi,corio_nsi, &
       &                        kfl_sgsli_nsi, kfl_rmom2_nsi, &
       &                        kfl_stabi_nsi
  use mod_tauadr, only       :  tauadr
  implicit none
  integer(ip), intent(in)    :: pgaus,pnode
  real(rp),    intent(in)    :: chale(2),hleng(ndime)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(inout) :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gppor(pgaus)
  real(rp),    intent(out)   :: gpst1(pgaus),gpst2(pgaus)
  real(rp),    intent(out)   :: gpsp1(pgaus),gpsp2(pgaus)
  real(rp),    intent(out)   :: gptt1(pgaus),gptt2(pgaus)
  real(rp),    intent(out)   :: rmom1(pnode,pgaus)
  real(rp),    intent(in)    :: dtsgs
  real(rp),    intent(in)    :: dtinv_loc
  real(rp),    intent(inout) :: tamin_loc
  real(rp),    intent(inout) :: tamax_loc
  real(rp),    intent(in)    :: elvel(ndime, pnode)
  real(rp),    intent(in)    :: gprhs(ndime, pgaus)
  real(rp),    intent(in)    :: gpgpr(ndime, pgaus)
  real(rp),    intent(in)    :: rmom2(ndime, ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gpgve(ndime, ndime, pnode)
  integer(ip)                :: idime,igaus,inode,jdime
  real(rp)                   :: adv,dif,rea,h2,gpvno, subgs(ndime), gpres(ndime)

  if( kfl_stabi_nsi == -1 ) then ! no stabilization 

     gpst1 = 0.0_rp
     gpst2 = 0.0_rp
     gpsp1 = 0.0_rp
     gpsp2 = 0.0_rp
     gptt1 = 0.0_rp
     gptt2 = 0.0_rp

  else

     h2 = chale(2) * chale(2)

     !----------------------------------------------------------------------
     !
     ! TAU1 and TAU2
     !
     !----------------------------------------------------------------------

     do igaus = 1,pgaus
        !
        ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
        ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]
        !
        call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
        adv   = gpden(igaus)*gpvno                                 ! Convective term: rho*|u+u'|
        dif   = gpvis(igaus)                                       ! Viscous term:    mu
        rea   = 2.0_rp*gpden(igaus)*corio_nsi + abs(gppor(igaus))  ! Coriolis: w + Porosity: sig
        call tauadr(&
             kfl_taust_nsi,staco_nsi,adv,dif,rea,&
             chale(1),chale(2),gpst1(igaus),gpden(igaus)*dtinv_loc)
        if( gpst1(igaus) /= 0.0_rp ) then
           if( kfl_ellen_nsi == 6 ) then   ! MIXLE
              gpst2(igaus) = staco_nsi(4)*(4.0_rp*dif + 2.0_rp*adv*hleng(1))    ! Use maximum elem length only for tau2
           else
              gpst2(igaus) = staco_nsi(4)*h2/gpst1(igaus)
           end if
        else
           gpst2(igaus) = 0.0_rp
        end if
        !
        ! (tau1',tau2') and (tau1'/tau1,tau2'/tau2)
        !     
        if( kfl_sgsti_nsi == 1 ) then
           gpsp1(igaus) = 1.0_rp / ( gpden(igaus) * dtsgs + 1.0_rp / gpst1(igaus) )
           gpsp2(igaus) = gpst2(igaus)
           gptt1(igaus) = gpsp1(igaus) / gpst1(igaus)
           gptt2(igaus) = 1.0_rp
        else
           gpsp1(igaus) = gpst1(igaus)
           gpsp2(igaus) = gpst2(igaus)
           gptt1(igaus) = 1.0_rp
           gptt2(igaus) = 1.0_rp
        end if
        if (kfl_stabi_nsi==1) gptt1(igaus)=1.0_rp

     end do

     !----------------------------------------------------------------------
     !
     ! RMOMU: Add advection to residual; split OSS needs only advection
     !
     !----------------------------------------------------------------------

     if (kfl_sgsli_nsi==0) then ! linear subscales model and convection tracking
        ! calculus of residual
        do igaus = 1,pgaus
           gpres(1:ndime) = gprhs(1:ndime,igaus) - gpgpr(1:ndime,igaus)
           do inode=1,pnode
              gpres(1:ndime) = gpres(1:ndime)&
                   - rmom1(inode,igaus) * elvel(1:ndime,inode)
           end do
           do jdime =1, ndime
              gpres(1:ndime) = gpres(1:ndime) - gpden(igaus)&
                   * gpadv(jdime,igaus) * gpgve(jdime,1:ndime,igaus)
           end do
           if( kfl_rmom2_nsi /= 0 ) then ! for coriolis problem
              do inode = 1,pnode
                 do jdime = 1,ndime
                    gpres(1:ndime) = gpres(1:ndime)&
                         - rmom2(1:ndime,jdime,inode,igaus) * elvel(jdime,inode)
                 end do
              end do
           end if
           ! linear subgrid scale model usgs = tau * res
           subgs(1:ndime)        = gpst1(igaus)*gpres(1:ndime)
           ! add subgrid scale to advective velocity (in residual and in test function)
           gpadv(1:ndime, igaus) = gpadv(1:ndime, igaus) + subgs(1:ndime)
        end do
     end if
     do igaus = 1,pgaus
        do inode = 1,pnode
           do idime = 1,ndime
              rmom1(inode,igaus) = rmom1(inode,igaus)& 
                   + gpden(igaus) * gpadv(idime,igaus)&
                   * gpcar(idime,inode,igaus)
           end do
        end do
     end do

     !----------------------------------------------------------------------
     !
     ! TAMIN_NSI, TAMAX_NSI: Minimum and maximum tau
     !
     !----------------------------------------------------------------------

     do igaus = 1,pgaus
        tamax_loc = max(tamax_loc,gpsp1(igaus))   
        tamin_loc = min(tamin_loc,gpsp1(igaus))   
     end do

  end if

end subroutine nsi_elmsgs
