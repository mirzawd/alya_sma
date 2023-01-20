!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmma4(&
     pnode,pgaus,pevat,gpden,gpvis,gppor,gpsp1,gpsp2,  &
     gpvol,gpsha,gpcar,gpadv,gpvep,gpprp,gpgrp,gprhs,  &
     gprhc,gpvel,gpsgs,wgrgr,agrau,elvel,elpre,elbub,  &
     elauu,elaup,elapp,elapu,elrbu,elrbp,dtinv_loc,    &
     dtsgs,pbubl,gpsha_bub,gpcar_bub,elauq,elapq,elaqu,&
     elaqp,elaqq,elrbq, densi)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmma4
  ! NAME 
  !    nsi_elmma4
  ! DESCRIPTION
  !    Compute element matrix and RHS
  ! USES
  ! USED BY
  !    nsi_elmop3
  !***
  !----------------------------------------------------------------------- 

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_master, only       :  kfl_lumped
  use def_nastin, only       :  kfl_stabi_nsi,kfl_regim_nsi,&
       &                        kfl_sgsti_nsi,kfl_limit_nsi,&
       &                        fvins_nsi,nbdfp_nsi,pabdf_nsi,&
       &                        kfl_nota1_nsi,NSI_FRACTIONAL_STEP,&
       &                        fcons_nsi,penal_nsi,&
       &                        kfl_press_nsi,NSI_GALERKIN, corio_nsi,&
       &                        fvela_nsi 
  use mod_maths,  only       :  maths_schur_complement
  implicit none

  integer(ip), intent(in)    :: pnode,pgaus,pevat
  real(rp),    intent(in)    :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)    :: gppor(pgaus)
  real(rp),    intent(in)    :: gpsp1(pgaus)
  real(rp),    intent(in)    :: gpsp2(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(inout) :: gpvep(ndime,pgaus)
  real(rp),    intent(inout) :: gpprp(pgaus)      
  real(rp),    intent(inout) :: gpgrp(ndime,pgaus)
  real(rp),    intent(inout) :: gprhs(ndime,pgaus)
  real(rp),    intent(inout) :: gprhc(pgaus)
  real(rp),    intent(in)    :: gpvel(ndime,pgaus,*)
  real(rp),    intent(in)    :: gpsgs(ndime,pgaus,*)
  real(rp),    intent(out)   :: wgrgr(pnode,pnode,pgaus)
  real(rp),    intent(out)   :: agrau(pnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime,pnode,*)
  real(rp),    intent(in)    :: elpre(pnode,*)
  real(rp),    intent(in)    :: elbub
  
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elrbu(ndime,pnode)
  real(rp),    intent(out)   :: elrbp(pnode)  
  real(rp),    intent(in)    :: dtinv_loc
  real(rp),    intent(in)    :: dtsgs
  integer(ip), intent(in)    :: pbubl
  real(rp),    intent(in)    :: gpsha_bub(pgaus)                    ! Ne
  real(rp),    intent(in)    :: gpcar_bub(ndime,pgaus)              ! dNe/dxi
  ! Enrichement Element matrices
  real(rp),    intent(out)   :: elauq(pnode*ndime,1)
  real(rp),    intent(out)   :: elapq(pnode,1)
  real(rp),    intent(out)   :: elaqu(1,pnode*ndime)
  real(rp),    intent(out)   :: elaqp(1,pnode)
  real(rp),    intent(out)   :: elaqq(1,1)
  real(rp),    intent(out)   :: elrbq(1)
  real(rp),    intent(in)    :: densi(pgaus,nbdfp_nsi)

  integer(ip)                :: inode,jnode,idofv,jdof2,jdof3
  integer(ip)                :: idof1,idof3,idof2,igaus,idime,jdof1,jdofv,itime
  real(rp)                   :: fact0,fact1,fact2,fact3,fact4,fact5,fact6
  real(rp)                   :: fact7,fact8,c1,c2,c3,c4,alpha,beta
  real(rp)                   :: gpveo(3),fact1_p
  real(rp)                   :: gpsp1_p(pgaus)
  real(rp)                   :: gpsp1_v(pgaus)
  real(rp)                   :: gpsp2_v(pgaus)
  real(rp)                   :: dtinv_mod

  if( NSI_FRACTIONAL_STEP ) then
     dtinv_mod = 0.0_rp
  else
     dtinv_mod = dtinv_loc
  end if

  !----------------------------------------------------------------------
  !
  ! possibility of using only pressure stabilization - not ready with limiter - nor with shock capturing
  !
  !----------------------------------------------------------------------

  gpsp1_p = gpsp1
  gpsp1_v = gpsp1
  gpsp2_v = gpsp2

  if( kfl_nota1_nsi == 1 ) gpsp1_v = 0.0_rp 

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  elrbp = 0.0_rp
  elrbu = 0.0_rp
  elauu = 0.0_rp
  elaup = 0.0_rp
  elapu = 0.0_rp
  elapp = 0.0_rp

  !----------------------------------------------------------------------
  !
  ! Test functions
  !
  !----------------------------------------------------------------------

  !
  ! AGRAU = rho * (a.grad) Ni
  ! WGRGR = grad(Ni) . grad(Nj)
  !
  if( ndime == 2 ) then

     do igaus = 1,pgaus
        do inode = 1,pnode
           agrau(inode,igaus) =  gpden(igaus) * (                    &
                &                gpadv(1,igaus)*gpcar(1,inode,igaus) &
                &              + gpadv(2,igaus)*gpcar(2,inode,igaus) )

           do jnode = 1,pnode
              wgrgr(inode,jnode,igaus) = &
                   &             gpcar(1,inode,igaus)*gpcar(1,jnode,igaus) &
                   &           + gpcar(2,inode,igaus)*gpcar(2,jnode,igaus) 
           end do
        end do
     end do

  else

     do igaus = 1,pgaus
        do inode = 1,pnode
           agrau(inode,igaus) =  gpden(igaus) * (                    &
                &                gpadv(1,igaus)*gpcar(1,inode,igaus) &
                &              + gpadv(2,igaus)*gpcar(2,inode,igaus) &
                &              + gpadv(3,igaus)*gpcar(3,inode,igaus) )

           do jnode = 1,pnode
              wgrgr(inode,jnode,igaus) = &
                   &             gpcar(1,inode,igaus)*gpcar(1,jnode,igaus) &
                   &           + gpcar(2,inode,igaus)*gpcar(2,jnode,igaus) & 
                   &           + gpcar(3,inode,igaus)*gpcar(3,jnode,igaus) 
           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Auu
  !
  !----------------------------------------------------------------------

  !
  ! Galerkin + ( tau2 * div(u) , div(v) ) + ( tau1 * rho*a.grad(u), rho*a.grad(v) )
  !
  if( ndime == 2 ) then

     do igaus = 1,pgaus

        fact0 = gpsp2_v(igaus) * gpvol(igaus)
        fact6 = gpvis(igaus)   * gpvol(igaus)
        fact7 = gpsp1_v(igaus) * gpvol(igaus) 
        fact8 = pabdf_nsi(1)   * gpden(igaus) * dtinv_mod + gppor(igaus)

        do inode = 1,pnode

           idof1 = 2*inode-1
           idof2 = idof1+1

           fact1 = fact0 * gpcar(1,inode,igaus)      ! div(u) * tau2' * dv/dx
           fact2 = fact0 * gpcar(2,inode,igaus)      ! div(u) * tau2' * dv/dy
           fact4 = gpsha(inode,igaus) * gpvol(igaus)

           do jnode = 1,pnode    

              jdof1 = 2*jnode-1
              jdof2 = jdof1+1

              fact5 = fact4 * ( agrau(jnode,igaus) + fact8 * gpsha(jnode,igaus) ) &          ! ( rho/dt N_j + s Nj + rho*(a.grad)Nj ) Ni
                   +  fact6 * wgrgr(inode,jnode,igaus) &                                     ! mu * grad(Ni) . grad(Nj)
                   +  fact7 * agrau(jnode,igaus) * agrau(inode,igaus)                        ! tau1 * rho*(a.grad)Nj * rho*(a.grad)Ni

              elauu(idof1,jdof1) = elauu(idof1,jdof1) + fact1 * gpcar(1,jnode,igaus) + fact5
              elauu(idof2,jdof1) = elauu(idof2,jdof1) + fact2 * gpcar(1,jnode,igaus)
              elauu(idof1,jdof2) = elauu(idof1,jdof2) + fact1 * gpcar(2,jnode,igaus) 
              elauu(idof2,jdof2) = elauu(idof2,jdof2) + fact2 * gpcar(2,jnode,igaus) + fact5

           end do

        end do
     end do
     !
     ! Skew symmetric convective operator
     ! 
     if(fcons_nsi > 0.1_rp) then 
        do igaus = 1,pgaus
           fact0 = fcons_nsi * gpden(igaus) * gpvol(igaus)
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode              
                 fact1 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) * gpsha(jnode,igaus) * fact0
                 fact2 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,jnode,igaus)) * gpsha(inode,igaus) * fact0
                 jdof1 = 2*jnode-1
                 jdof2 = jdof1+1

                 elauu(idof1,jdof1) = elauu(idof1,jdof1) - fact1 - fact2                     
                 elauu(idof2,jdof2) = elauu(idof2,jdof2) - fact1 - fact2                     
              end do
           end do
        end do
        !
        ! Low-Mach with skew symmetric (adding temporal term: d(rho*u)/dt)
        !
        if (kfl_regim_nsi==3) then 
            do itime =2, nbdfp_nsi         ! only rhs with Galerkin needs to be modified
               do igaus =1, pgaus
                  gpveo(1:ndime) = 0.0_rp  ! velocity at itime 
                  do inode =1, pnode
                     gpveo(1:ndime) = gpveo(1:ndime) + elvel(1:ndime, inode, itime)* gpsha(inode, igaus)
                  end do
                  fact0 = 0.5_rp*(gpden(igaus) - densi(igaus, itime))*pabdf_nsi(itime)*dtinv_loc*gpvol(igaus)
                  do inode =1, pnode
                     elrbu(1:ndime, inode) = elrbu(1:ndime, inode) &
                                           + fact0*gpsha(inode,igaus)*gpveo(1:ndime)
                  end do
               end do
            end do
         end if 
     end if
  
     ! EMAC
     if(fcons_nsi < -0.1_rp) then 
      do igaus = 1,pgaus
           fact0 =  gpden(igaus) * gpvol(igaus)
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode              
                 fact1 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) * gpsha(jnode,igaus) * fact0
                 fact2 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,jnode,igaus)) * gpsha(inode,igaus) * fact0
                 fact3 = gpden(igaus)*gpsha(inode,igaus)*gpvol(igaus)
                 jdof1 = 2*jnode-1
                 jdof2 = jdof1+1

                 elauu(idof1,jdof1) = elauu(idof1,jdof1) - fact1 - fact2 + (gpadv(1,igaus)*gpcar(1,jnode,igaus))*fact3 
                 elauu(idof2,jdof1) = elauu(idof2,jdof1)                 + (gpadv(1,igaus)*gpcar(2,jnode,igaus))*fact3 

                 elauu(idof2,jdof2) = elauu(idof2,jdof2) - fact1 - fact2 + (gpadv(2,igaus)*gpcar(2,jnode,igaus))*fact3 
                 elauu(idof1,jdof2) = elauu(idof1,jdof2)                 + (gpadv(2,igaus)*gpcar(1,jnode,igaus))*fact3 
              end do
           end do
        end do
     end if
  else  !ndime = 3

     do igaus = 1,pgaus

        fact0 = gpsp2_v(igaus) * gpvol(igaus)
        fact6 = gpvis(igaus)   * gpvol(igaus)
        fact7 = gpsp1_v(igaus) * gpvol(igaus)
        fact8 = pabdf_nsi(1)   * gpden(igaus) * dtinv_mod + gppor(igaus)

        do inode = 1,pnode

           idof1 = 3 * inode - 2
           idof2 = idof1     + 1
           idof3 = idof2     + 1

           fact1 = fact0 * gpcar(1,inode,igaus)      ! div(u) * tau2' * dv/dx
           fact2 = fact0 * gpcar(2,inode,igaus)      ! div(u) * tau2' * dv/dy
           fact3 = fact0 * gpcar(3,inode,igaus)      ! div(u) * tau2' * dv/dz
           fact4 = gpsha(inode,igaus) * gpvol(igaus)

           do jnode = 1,pnode    

              jdof1 = 3 * jnode - 2
              jdof2 = jdof1     + 1
              jdof3 = jdof2     + 1

              fact5 = fact4 * ( agrau(jnode,igaus) + fact8 * gpsha(jnode,igaus) ) &          ! ( rho/dt N_j + s Nj + rho*(a.grad)Nj ) Ni
                   +  fact6 * wgrgr(inode,jnode,igaus) &                                     ! mu * grad(Ni) . grad(Nj)
                   +  fact7 * agrau(jnode,igaus) * agrau(inode,igaus)                        ! t1 * rho*(a.grad)Nj * rho*(a.grad)Ni

              elauu(idof1,jdof1) = elauu(idof1,jdof1) + fact1 * gpcar(1,jnode,igaus) + fact5
              elauu(idof2,jdof1) = elauu(idof2,jdof1) + fact2 * gpcar(1,jnode,igaus)
              elauu(idof3,jdof1) = elauu(idof3,jdof1) + fact3 * gpcar(1,jnode,igaus)

              elauu(idof1,jdof2) = elauu(idof1,jdof2) + fact1 * gpcar(2,jnode,igaus) 
              elauu(idof2,jdof2) = elauu(idof2,jdof2) + fact2 * gpcar(2,jnode,igaus) + fact5
              elauu(idof3,jdof2) = elauu(idof3,jdof2) + fact3 * gpcar(2,jnode,igaus) 

              elauu(idof1,jdof3) = elauu(idof1,jdof3) + fact1 * gpcar(3,jnode,igaus) 
              elauu(idof2,jdof3) = elauu(idof2,jdof3) + fact2 * gpcar(3,jnode,igaus)
              elauu(idof3,jdof3) = elauu(idof3,jdof3) + fact3 * gpcar(3,jnode,igaus) + fact5

           end do

        end do
     end do
     !
     ! skew symmetric convective operator
     !
     if(fcons_nsi > 0.1_rp) then 
        do igaus = 1,pgaus
           fact0 = fcons_nsi * gpden(igaus) * gpvol(igaus)
           do inode = 1,pnode
              idof1 = 3 * inode - 2
              idof2 = 3 * inode - 1
              idof3 = 3 * inode
              do jnode = 1,pnode
                 jdof1 = 3 * jnode - 2
                 jdof2 = 3 * jnode - 1
                 jdof3 = 3 * jnode 
                 fact1 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) * gpsha(jnode,igaus) * fact0
                 fact2 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,jnode,igaus)) * gpsha(inode,igaus) * fact0
                 
                 elauu(idof1,jdof1) = elauu(idof1,jdof1) - fact1 - fact2   
                 elauu(idof2,jdof2) = elauu(idof2,jdof2) - fact1 - fact2   
                 elauu(idof3,jdof3) = elauu(idof3,jdof3) - fact1 - fact2                     
              end do
           end do
        end do
        !
        ! Low-Mach with skew symmetric (adding temporal term: d(rho*u)/dt)
        !
        if (kfl_regim_nsi==3) then 
           do itime =2, nbdfp_nsi ! only rhs with Galerkin needs to be modified
              do igaus =1, pgaus
                  gpveo(1:ndime) = 0.0_rp  ! velocity at itime 
                  do inode =1, pnode
                     gpveo(1:ndime) = gpveo(1:ndime) + elvel(1:ndime, inode, itime)* gpsha(inode, igaus)
                  end do
                  fact0 = 0.5_rp*(gpden(igaus) - densi(igaus, itime))*pabdf_nsi(itime)*dtinv_loc
                  do inode =1, pnode
                     elrbu(1:ndime, inode) = elrbu(1:ndime, inode) &
                          + fact0*gpsha(inode,igaus)*gpveo(1:ndime)
                  end do
               end do
           end do
        end if !end lowmach
     end if
     ! EMAC
     if(fcons_nsi < -0.1_rp) then 
        do igaus = 1,pgaus
           fact0 =  gpden(igaus) * gpvol(igaus)
           do inode = 1,pnode
              idof1 = 3 * inode - 2
              idof2 = 3 * inode - 1
              idof3 = 3 * inode
              do jnode = 1,pnode
                 jdof1 = 3 * jnode - 2
                 jdof2 = 3 * jnode - 1
                 jdof3 = 3 * jnode 
                 fact1 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) * gpsha(jnode,igaus) * fact0
                 fact2 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,jnode,igaus)) * gpsha(inode,igaus) * fact0
                 fact3 = gpden(igaus)*gpsha(inode,igaus)*gpvol(igaus)

                 elauu(idof1,jdof1) = elauu(idof1,jdof1) - fact1 - fact2 + (gpadv(1,igaus)*gpcar(1,jnode,igaus))*fact3
                 elauu(idof2,jdof1) = elauu(idof2,jdof1)                 + (gpadv(1,igaus)*gpcar(2,jnode,igaus))*fact3
                 elauu(idof3,jdof1) = elauu(idof3,jdof1)                 + (gpadv(1,igaus)*gpcar(3,jnode,igaus))*fact3

                 elauu(idof1,jdof2) = elauu(idof1,jdof2)                 + (gpadv(2,igaus)*gpcar(1,jnode,igaus))*fact3
                 elauu(idof2,jdof2) = elauu(idof2,jdof2) - fact1 - fact2 + (gpadv(2,igaus)*gpcar(2,jnode,igaus))*fact3
                 elauu(idof3,jdof2) = elauu(idof3,jdof2)                 + (gpadv(2,igaus)*gpcar(3,jnode,igaus))*fact3

                 elauu(idof1,jdof3) = elauu(idof1,jdof3)                 + (gpadv(3,igaus)*gpcar(1,jnode,igaus))*fact3                 
                 elauu(idof2,jdof3) = elauu(idof2,jdof3)                 + (gpadv(3,igaus)*gpcar(2,jnode,igaus))*fact3                 
                 elauu(idof3,jdof3) = elauu(idof3,jdof3) - fact1 - fact2 + (gpadv(3,igaus)*gpcar(3,jnode,igaus))*fact3                 

              end do
           end do
        end do
     end if
     
     ! Coriolis operator
     
     if( corio_nsi > 1.0e-8_rp ) then
        do igaus = 1,pgaus   
           fact0 = 2.0_rp * gpden(igaus) *gpvol(igaus)
!!$       fact1 = fact0  * fvela_nsi(1)
!!$       fact2 = fact0  * fvela_nsi(2)
           fact3 = fact0  * fvela_nsi(3)
!!$        do inode = 1,pnode
!!$           rmom2(1,2,inode,igaus) = - fact3 * gpsha(inode,igaus)  ! -wz*uy ! Fx  
!!$           rmom2(1,3,inode,igaus) =   fact2 * gpsha(inode,igaus)  !  wy*uz  
!!$           rmom2(2,1,inode,igaus) =   fact3 * gpsha(inode,igaus)  !  wz*ux ! Fy 
!!$           rmom2(2,3,inode,igaus) = - fact1 * gpsha(inode,igaus)  ! -wx*uz  
!!$           rmom2(3,1,inode,igaus) = - fact2 * gpsha(inode,igaus)  ! -wy*ux  
!!$           rmom2(3,2,inode,igaus) =   fact1 * gpsha(inode,igaus)  !  wx*uy                  
!!$        end do
           do inode =1, pnode
              idof1 = 3 * inode - 2
              idof2 = 3 * inode - 1
  !           idof3 = 3 * inode
              do jnode = 1,pnode
                 jdof1 = 3 * jnode - 2
                 jdof2 = 3 * jnode - 1
  !              jdof3 = 3 * jnode 
                 elauu(idof1,jdof2) = elauu(idof1,jdof2)  -  fact3 * gpsha(jnode,igaus)*gpsha(inode, igaus)
                 elauu(idof2,jdof1) = elauu(idof2,jdof1)  +  fact3 * gpsha(jnode,igaus)*gpsha(inode, igaus)  
              end do
           end do
        end do
     end if  ! corio
  end if ! ndime =3

  if( fvins_nsi > 0.9_rp ) then
     !
     ! ( mu*duj/dxi , dv/dxj ) (only div form)
     !
     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              idofv = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idofv+1
                 do jnode = 1,pnode
                    jdofv = (jnode-1)*ndime
                    fact1 = gpvis(igaus) * gpvol(igaus) * gpcar(idime,jnode,igaus)     
                    jdofv = jdofv + 1
                    elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(1,inode,igaus)
                    jdofv = jdofv + 1
                    elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(2,inode,igaus)
                 end do
                 if( fvins_nsi == 2.0_rp ) then
                    fact1 = -2.0_rp/3.0_rp * gpvis(igaus) * gpvol(igaus) * gpcar(idime,inode,igaus)
                    do jnode = 1,pnode
                       jdofv = (jnode-1)*ndime   
                       jdofv = jdofv+1
                       elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(1,jnode,igaus)
                       jdofv = jdofv+1
                       elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(2,jnode,igaus)
                    end do
                 end if
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              idofv = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idofv + 1
                 do jnode = 1,pnode
                    jdofv = (jnode-1)*ndime
                    fact1 = gpvis(igaus) * gpvol(igaus) * gpcar(idime,jnode,igaus)     
                    jdofv = jdofv + 1
                    elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(1,inode,igaus)
                    jdofv = jdofv + 1
                    elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(2,inode,igaus)
                    jdofv = jdofv + 1
                    elauu(idofv,jdofv)=  elauu(idofv,jdofv)+  fact1  *gpcar(3,inode,igaus)
                 end do
                 if( fvins_nsi == 2.0_rp ) then
                    fact1 = -2.0_rp/3.0_rp*gpvis(igaus) * gpvol(igaus) * gpcar(idime,inode,igaus)
                    do jnode = 1,pnode
                       jdofv = (jnode-1)*ndime   
                       jdofv = jdofv + 1
                       elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(1,jnode,igaus)
                       jdofv = jdofv + 1
                       elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(2,jnode,igaus)
                       jdofv = jdofv + 1
                       elauu(idofv,jdofv) = elauu(idofv,jdofv) + fact1 * gpcar(3,jnode,igaus)
                    end do
                 end if
              end do
           end do
        end do
     end if
  end if
  !
  ! Lumped evolution matrix (only backward euler)
  !
  if( kfl_lumped == 1 ) then 
     !
     ! Remove Galerkin term and add lumped term 
     ! 
     if( ndime == 2 ) then
        do igaus =1, pgaus
           do inode =1, pnode
              idof1 = (inode-1)*2+1
              idof2 = (inode-1)*2+2
              fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(1)
              elauu(idof1,idof1) =  elauu(idof1,idof1) + fact0
              elauu(idof2,idof2) =  elauu(idof2,idof2) + fact0
              do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                 jdof1 = (jnode-1)*2+1
                 jdof2 = (jnode-1)*2+2
                 elauu(idof1,jdof1) =  elauu(idof1,jdof1) - fact0*gpsha(jnode, igaus)
                 elauu(idof2,jdof2) =  elauu(idof2,jdof2) - fact0*gpsha(jnode, igaus)
              end do
           end do
        end do
        do itime =2, nbdfp_nsi ! RHS
           do igaus =1, pgaus
              gpveo(1:2) = 0.0_rp
              do inode =1, pnode
                 gpveo(1:2) = gpveo(1:2) + elvel(1:2, inode, itime)* gpsha(inode, igaus)
              end do
              do inode =1, pnode
                 fact0 = - gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(itime)
                 elrbu(1:2,inode)   =  elrbu(1:2,inode)   - fact0*gpveo(1:2)
                 elrbu(1:2,inode)   =  elrbu(1:2,inode)   + fact0*elvel(1:2, inode, itime)
              end do
           end do
        end do
     else
        do igaus =1, pgaus
           do inode =1, pnode
              idof1 = 3 * inode - 2
              idof2 = idof1 + 1
              idof3 = idof2 + 1
              fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(1)
              elauu(idof1,idof1) =  elauu(idof1,idof1) + fact0
              elauu(idof2,idof2) =  elauu(idof2,idof2) + fact0
              elauu(idof3,idof3) =  elauu(idof3,idof3) + fact0
              do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                 jdof1 = 3 * jnode - 2
                 jdof2 = jdof1 + 1
                 jdof3 = jdof2 + 1
                 elauu(idof1,jdof1) =  elauu(idof1,jdof1) - fact0*gpsha(jnode, igaus)
                 elauu(idof2,jdof2) =  elauu(idof2,jdof2) - fact0*gpsha(jnode, igaus)
                 elauu(idof3,jdof3) =  elauu(idof3,jdof3) - fact0*gpsha(jnode, igaus)
              end do
           end do
        end do
        do itime =2, nbdfp_nsi ! RHS
           do igaus =1, pgaus
              gpveo(1:3) = 0.0_rp
              do inode =1, pnode
                 gpveo(1:3) = gpveo(1:3) + elvel(1:3, inode, itime)* gpsha(inode, igaus)
              end do
              do inode =1, pnode
                 fact0 = - gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(itime)
                 elrbu(1:3,inode)   =  elrbu(1:3,inode)   - fact0*gpveo(1:3)
                 elrbu(1:3,inode)   =  elrbu(1:3,inode)   + fact0*elvel(1:3, inode, itime)
              end do
           end do
        end do
     end if

  else if( kfl_lumped == 2 ) then 
     !
     ! No time term have been added up to now: add Galerkin term
     !
     do igaus = 1,pgaus
        fact0 = gpvol(igaus) * gpden(igaus) * dtinv_mod
        do inode = 1, pnode
           fact1 = fact0 * gpsha(inode,igaus)
           do idime = 1,ndime
              idof1              = (inode-1) * ndime + idime
              elauu(idof1,idof1) = elauu(idof1,idof1) + fact1
              elrbu(idime,inode) = elrbu(idime,inode) + fact1 * elvel(idime,inode,2)
           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Apu and Aup
  !
  !----------------------------------------------------------------------
  !
  ! ( div(u) , q ) and - ( p , div(v) ) 
  !
  if( ndime == 2 ) then
     do igaus = 1,pgaus
        do inode = 1,pnode
           idof1 = 2 * inode - 1
           idof2 = idof1 + 1
           do jnode = 1,pnode
              fact0              = gpvol(igaus) * gpsha(jnode,igaus) 
              fact1              = fact0 * gpcar(1,inode,igaus)
              fact2              = fact0 * gpcar(2,inode,igaus)
              elapu(jnode,idof1) = elapu(jnode,idof1) + fact1
              elapu(jnode,idof2) = elapu(jnode,idof2) + fact2
              elaup(idof1,jnode) = elaup(idof1,jnode) - fact1
              elaup(idof2,jnode) = elaup(idof2,jnode) - fact2
           end do
        end do
     end do
  else
     do igaus = 1,pgaus
        do inode = 1,pnode
           idof1 = 3 * inode - 2
           idof2 = idof1 + 1
           idof3 = idof2 + 1
           do jnode = 1,pnode
              fact0              = gpvol(igaus) * gpsha(jnode,igaus) 
              fact1              = fact0 * gpcar(1,inode,igaus)
              fact2              = fact0 * gpcar(2,inode,igaus)
              fact3              = fact0 * gpcar(3,inode,igaus)
              elapu(jnode,idof1) = elapu(jnode,idof1) + fact1
              elapu(jnode,idof2) = elapu(jnode,idof2) + fact2
              elapu(jnode,idof3) = elapu(jnode,idof3) + fact3
              elaup(idof1,jnode) = elaup(idof1,jnode) - fact1
              elaup(idof2,jnode) = elaup(idof2,jnode) - fact2
              elaup(idof3,jnode) = elaup(idof3,jnode) - fact3
           end do
        end do
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! App
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: ( tau1' * grad(p) , grad(q) )
  ! 
  if( kfl_stabi_nsi /= -1 ) then
     do igaus = 1,pgaus
        do inode = 1,pnode
           do jnode = inode+1,pnode
              fact1 = gpsp1_p(igaus) * wgrgr(jnode,inode,igaus) * gpvol(igaus)
              elapp(jnode,inode) = elapp(jnode,inode) + fact1
              elapp(inode,jnode) = elapp(inode,jnode) + fact1
           end do
           fact1 = gpsp1_p(igaus) * wgrgr(inode,inode,igaus) * gpvol(igaus)
           elapp(inode,inode) = elapp(inode,inode) + fact1
        end do
     end do
  end if
  !
  ! Penalization
  !
  !
  ! Penalization
  !
  do igaus = 1,pgaus
     fact1 = penal_nsi * gpvol(igaus)
     do inode = 1,pnode
        elapp(inode,inode) = elapp(inode,inode) + fact1 * gpsha(inode,igaus)
        elrbp(inode)       = elrbp(inode)       + fact1 * gpsha(inode,igaus) * elpre(inode,1) 
     end do
  end do
 
  !----------------------------------------------------------------------
  !
  ! bu and bp
  !
  ! P1  = P [ tau1' * rho * a.grad(u) ]
  ! P1' = P1 + tau1' * rho * u'n / dt
  !
  ! P2  = P [ tau1' * ( grad(p) - f ) ]
  ! P2' = P2 + tau1' * rho * u'n / dt + tau1' * f 
  !
  !----------------------------------------------------------------------
  !
  ! Limiter
  !
  if( kfl_limit_nsi == -1 ) then

     do igaus = 1,pgaus
        do idime = 1,ndime
           gpvep(idime,igaus) = 0.0_rp
        end do
     end do

  else if( kfl_limit_nsi > 0 ) then

     do igaus = 1,pgaus
        c1 = 0.0_rp
        c2 = 0.0_rp
        c3 = 0.0_rp
        do idime = 1,ndime
           c4 = 0.0_rp
           do inode = 1,pnode
              c4 = c4 + agrau(inode,igaus) * elvel(idime,inode,1)
           end do
           c4 = gpsp1(igaus) * c4
           c1 = c1 + ( gpvep(idime,igaus) - c4 )**2
           c3 = c3 + gpvep(idime,igaus) * gpvep(idime,igaus)
           c2 = c2 + c4 * c4
        end do
        c3 = sqrt( c2 ) + sqrt( c3 )
        c1 = sqrt( c1 )
        if( c3 /= 0.0_rp ) then
           beta  = c1 / c3
        else
           beta  = 0.0_rp
        end if
        if( kfl_limit_nsi == 1 ) then
           alpha = min(1.0_rp,2.0_rp*(1.0_rp-beta))
        else if( kfl_limit_nsi == 2 ) then
           alpha = 0.5_rp*(tanh(20.0_rp*(beta-0.8_rp))+1.0_rp)
        end if
        do idime = 1,ndime
           gpvep(idime,igaus) = alpha * gpvep(idime,igaus)
        end do
     end do

  end if
  !
  ! P2 <= P2 + tau1' * f
  !
  if( kfl_stabi_nsi == -1 ) then
     gpgrp = 0.0_rp
     gpvep = 0.0_rp
  else
     if( ndime == 2 ) then
        do igaus = 1,pgaus
           gpgrp(1,igaus) = gpgrp(1,igaus) + gpsp1_p(igaus) * gprhs(1,igaus) 
           gpgrp(2,igaus) = gpgrp(2,igaus) + gpsp1_p(igaus) * gprhs(2,igaus) 
        end do
     else
        do igaus = 1,pgaus
           gpgrp(1,igaus) = gpgrp(1,igaus) + gpsp1_p(igaus) * gprhs(1,igaus)
           gpgrp(2,igaus) = gpgrp(2,igaus) + gpsp1_p(igaus) * gprhs(2,igaus)
           gpgrp(3,igaus) = gpgrp(3,igaus) + gpsp1_p(igaus) * gprhs(3,igaus)
        end do
     end if
     !
     ! P1 <= P1 + tau1' * rho * u'n / dt
     ! P2 <= P2 + tau1' * rho * u'n / dt
     !
     if( kfl_sgsti_nsi == 1 ) then
        if( ndime == 2 ) then
           do igaus = 1,pgaus
              fact1          = gpden(igaus) * dtsgs * gpsp1_v(igaus)
              fact1_p        = gpden(igaus) * dtsgs * gpsp1_p(igaus)
              gpvep(1,igaus) = gpvep(1,igaus) + fact1 * gpsgs(1,igaus,2)
              gpvep(2,igaus) = gpvep(2,igaus) + fact1 * gpsgs(2,igaus,2)
              gpgrp(1,igaus) = gpgrp(1,igaus) + fact1_p * gpsgs(1,igaus,2)
              gpgrp(2,igaus) = gpgrp(2,igaus) + fact1_p * gpsgs(2,igaus,2)
           end do
        else
           do igaus = 1,pgaus 
              fact1 = gpden(igaus) * dtsgs * gpsp1_v(igaus)
              fact1_p = gpden(igaus) * dtsgs * gpsp1_p(igaus)
              gpvep(1,igaus) = gpvep(1,igaus) + fact1 * gpsgs(1,igaus,2)
              gpvep(2,igaus) = gpvep(2,igaus) + fact1 * gpsgs(2,igaus,2)
              gpvep(3,igaus) = gpvep(3,igaus) + fact1 * gpsgs(3,igaus,2)
              gpgrp(1,igaus) = gpgrp(1,igaus) + fact1_p * gpsgs(1,igaus,2)
              gpgrp(2,igaus) = gpgrp(2,igaus) + fact1_p * gpsgs(2,igaus,2)
              gpgrp(3,igaus) = gpgrp(3,igaus) + fact1_p * gpsgs(3,igaus,2)
           end do
        end if
     end if
  end if
  !
  ! bu = ( f + rho*u^n/dt , v ) + ( rho * a.grad(v) , tau1' * rho u'^n/dt + P1 ) 
  !    = ( f + rho*u^n/dt , v ) + ( rho * a.grad(v) , P1' ) 
  !
  ! bp = ( f + rho*u'^n/dt , tau1' grad(q) ) + ( P2 , grad(q) )
  !    = ( P2' , grad(q) ) 
  !
  if( ndime == 2 ) then
     do igaus = 1,pgaus
        fact4 = gpden(igaus) * dtinv_mod
        do itime = 2,nbdfp_nsi
           gprhs(1,igaus) = gprhs(1,igaus) - pabdf_nsi(itime) * fact4 * gpvel(1,igaus,itime)  
           gprhs(2,igaus) = gprhs(2,igaus) - pabdf_nsi(itime) * fact4 * gpvel(2,igaus,itime)
        end do
        do inode = 1,pnode
           fact1          = gpvol(igaus) * gpsha(inode,igaus)                                ! ( f + rho*u^n/dt , v )
           fact3          = gpvol(igaus) * agrau(inode,igaus)                                ! ( rho * a.grad(v) , P1' ) 
           elrbu(1,inode) = elrbu(1,inode) + fact1 * gprhs(1,igaus) + fact3 * gpvep(1,igaus) 
           elrbu(2,inode) = elrbu(2,inode) + fact1 * gprhs(2,igaus) + fact3 * gpvep(2,igaus) 
           elrbp(inode)   = elrbp(inode)   + gpvol(igaus) * ( &                              ! ( P2' , grad(q) ) 
                &    gpsha(inode,igaus)   * gprhc(igaus)    &
                &  + gpcar(1,inode,igaus) * gpgrp(1,igaus)  &
                &  + gpcar(2,inode,igaus) * gpgrp(2,igaus)  )
        end do
     end do
  else
     do igaus = 1,pgaus
        fact4 = gpden(igaus) * dtinv_mod
        do itime = 2,nbdfp_nsi
           gprhs(1,igaus) = gprhs(1,igaus) - pabdf_nsi(itime) * fact4 * gpvel(1,igaus,itime)  
           gprhs(2,igaus) = gprhs(2,igaus) - pabdf_nsi(itime) * fact4 * gpvel(2,igaus,itime)
           gprhs(3,igaus) = gprhs(3,igaus) - pabdf_nsi(itime) * fact4 * gpvel(3,igaus,itime)
        end do
        do inode = 1,pnode
           fact1          = gpvol(igaus) * gpsha(inode,igaus)
           fact3          = gpvol(igaus) * agrau(inode,igaus)
           elrbu(1,inode) = elrbu(1,inode) + fact1 * gprhs(1,igaus) + fact3 * gpvep(1,igaus) 
           elrbu(2,inode) = elrbu(2,inode) + fact1 * gprhs(2,igaus) + fact3 * gpvep(2,igaus) 
           elrbu(3,inode) = elrbu(3,inode) + fact1 * gprhs(3,igaus) + fact3 * gpvep(3,igaus) 
           elrbp(inode)   = elrbp(inode)   + gpvol(igaus) * ( &
                &    gpsha(inode,igaus)   * gprhc(igaus)   &
                &  + gpcar(1,inode,igaus) * gpgrp(1,igaus) &
                &  + gpcar(2,inode,igaus) * gpgrp(2,igaus) &
                &  + gpcar(3,inode,igaus) * gpgrp(3,igaus) )
        end do
     end do
  end if

  if( pbubl == 1 .and. kfl_stabi_nsi == NSI_GALERKIN ) then
     !
     ! Initialization
     !
     elauq = 0.0_rp
     elapq = 0.0_rp
     elaqu = 0.0_rp
     elaqp = 0.0_rp
     elaqq = 0.0_rp
     elrbq = 0.0_rp
     !
     ! Momentum equations
     !
     if( kfl_press_nsi == 1 ) then
        do igaus = 1,pgaus
           fact1 = gpvol(igaus) * gpsha_bub(igaus)
           do inode = 1,pnode
              idof1 = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idof1 + idime
                 elauq(idofv,1) = elauq(idofv,1) - fact1 * gpcar(idime,inode,igaus) 
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idof1 + idime
                 elauq(idofv,1) = elauq(idofv,1) + gpsha(inode,igaus) * gpcar_bub(idime,igaus)
              end do
           end do
        end do
     end if
     !
     ! Bubble equation
     !
     do igaus = 1,pgaus
        fact1 = gpvol(igaus) * gpsha_bub(igaus)
        do inode = 1,pnode
           do idime = 1,ndime
              idof1 = (inode-1)*ndime+idime
              elaqu(1,idof1) = elaqu(1,idof1) + gpcar(idime,inode,igaus) * fact1
           end do
        end do
     end do
     !
     ! Also pressure equation...
     !
     do igaus = 1,pgaus
        elrbq(1) = elrbq(1) + gpvol(igaus) * gprhc(igaus) * gpsha_bub(igaus)
     end do
     !
     ! Penalization
     !
     do igaus = 1,pgaus
        elaqq(1,1) = elaqq(1,1) + penal_nsi * gpvol(igaus) * gpsha_bub(igaus) 
        elrbq(1)   = elrbq(1)   + penal_nsi * gpvol(igaus) * gpsha_bub(igaus) * elbub 
     end do
     !
     ! ( rhs , q_bub )
     !
     !do igaus = 1,pgaus
     !    elrbq(1) = elrbq(1) + gprhc(igaus) * gpvol(igaus) * gpsha_bub(igaus)
     ! end do
     !
     ! Eliminate bubble
     !
     ! Auu <= Auu - Auq Aqq^-1 Aqu
     ! Aup <= Aup - Auq Aqq^-1 Aqp
     ! bu  <= bu  - Auq Aqq^-1 bq
     !
     ! Apu <= Apu - Apq Aqq^-1 Aqu
     ! App <= App - Apq Aqq^-1 Aqp
     ! bp  <= bp  - Apq Aqq^-1 bq
     ! 
     call maths_schur_complement(pevat,pevat, 1_ip,elauu,elauq,elaqq,elaqu)
     call maths_schur_complement(pevat,pnode, 1_ip,elaup,elauq,elaqq,elaqp)
     call maths_schur_complement(pevat, 1_ip, 1_ip,elrbu,elauq,elaqq,elrbq)
     
     call maths_schur_complement(pnode,pevat, 1_ip,elapu,elapq,elaqq,elaqu)
     call maths_schur_complement(pnode,pnode, 1_ip,elapp,elapq,elaqq,elaqp)
     call maths_schur_complement(pnode, 1_ip, 1_ip,elrbp,elapq,elaqq,elrbq)
        
  end if
!!$
end subroutine nsi_elmma4

