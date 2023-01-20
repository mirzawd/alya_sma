!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_bouope(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_bouope
  ! NAME 
  !    tem_bouope
  ! DESCRIPTION
  !    ORDER=1:
  !      Temperature equation, boundary operations
  !       ITASK = 0 ... Default behavior
  !             = 1 ... Calculates variational heat flux and mass matrix on nodes for postprocessing
  !             = 2 ... Calculates heat flx and surface over boundaries
  ! USES
  ! USED BY
  !    tem_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_ker_proper 
  use def_domain
  use def_temper
  use mod_ADR,    only : ADR_add_sgs_or_bubble
  use mod_solver, only : solver_assemble_element_matrix_scalar
  use mod_matrix, only : matrix_assemble_element_RHS
  use mod_ker_functions,  only : ker_functions
  use mod_ker_tendencies, only : kfl_tendencies_ker
  use mod_ker_regularization, only : regul_k, kfl_regularization
  use mod_bouder
  use mod_tem_turbul
  implicit none
  
  integer(ip), intent(in) :: itask
  
  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: elvel(ndime,mnode),gptem(mgaus)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: bocod(ndime,mnodb),botem(mnodb)
  real(rp)    :: bovel(ndime,mnodb)
  real(rp)    :: bogrc(ndime,mnodb)
  real(rp)    :: elmas(mnodb)
  integer(ip) :: ielem,inodb,ipoin
  integer(ip) :: igaus,igaub,iboun,pblty
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  integer(ip) :: dummi
  real(rp)    :: eucta,gbsur,gpdet,adotn
  real(rp)    :: gbsph(mgaub),gbden(mgaub),gbcon(mgaub),gbvis(mgaub),gbtem,gbvel(3)
  real(rp)    :: gpsph(mgaus),gpden(mgaus),gpdif(mgaus)
  real(rp)    :: gprea(mgaus)
  real(rp)    :: gptur(mgaus)                                 ! Turbulent viscosity
  real(rp)    :: arobi,trobi,qrobi,twall
  real(rp)    :: para1,para2,para3,para4
  real(rp)    :: xmrhs,xmmat
  real(rp)    :: gpcar(ndime,mnode,mgaus)
  real(rp)    :: gpcon(mgaus),dummr(ndime*mnode),gpcod, gpvis(mgaus)
  real(rp)    :: xjaci(9),xjacm(9), kinen, roughness
  real(rp)    :: temex,velex(3),veave(3)
  integer(ip)             :: kk, ifiel
  real(rp)                :: xx, temwa,rhocp, beta, tolbemol, ustar


  tolbemol = 1.0e-4
  ! Initialize nodal traction and boundary mass
  !
  if( itask > 0_ip ) then
     if( associated(heatf_tem) ) heatf_tem = 0.0_rp
     if( associated(massb_tem) ) massb_tem = 0.0_rp
  end if
  
  !
  ! Loop over elements  
  !
  boundaries: do iboun=1,nboun

     Neumann_or_Robin: if(                 &
          &  kfl_fixbo_tem(iboun) == 2 .or. &
          &  kfl_fixbo_tem(iboun) == 3 .or. &
          &  kfl_fixbo_tem(iboun) == 4 .or. &
          &  kfl_fixbo_tem(iboun) == 5 .or. &
          &  kfl_fixbo_tem(iboun) == 20 .or. & ! stable outflow condition when inflow
          &  bemol_tem > tolbemol ) then

        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        pelty = ltype(ielem)

        if( pelty > 0 ) then

           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           pmate = 1
           if( nmate > 1 ) pmate = lmate(ielem)
           !
           ! Inititalize
           !
           elmas(1:pnodb) = 0.0_rp
           elmat = 0.0_rp
           elrhs = 0.0_rp
           !
           ! Gather operations
           !
           bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
           botem(1:pnodb)         = therm(lnodb(1:pnodb,iboun),1)
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
           if( kfl_fixbo_tem(iboun) == 5 ) then
              bogrc(1:ndime,1:pnodb) = gradc_tem(1:ndime,lnodb(1:pnodb,iboun))
           end if
           !
           ! GPTEM+GPSGS: Temperature at Gauss point
           !
           call gather(&
                1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
                elmar(pelty)%shape,therm,gptem)
           call ADR_add_sgs_or_bubble(&
                ielem,pgaus,elmar(pelty) % shape_bub,ADR_tem,gptem) 
           !
           ! Cartesian derivatives
           !
           do igaus=1,pgaus
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)        ! and Jacobian
           end do
           !
           ! Properties 
           !
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('CONDU','PGAUS',dummi,ielem,gpcon,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar) 

           call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
           call ker_proper('CONDU','PGAUB',dummi,iboun,gbcon)
           call ker_proper('SPHEA','PGAUB',dummi,iboun,gbsph)
           call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis)

           if (kfl_fixbo_tem(iboun) == 3) &
                call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,elmar(pelty)%shape,gpcar)
           ! 
           ! Coupling with turbul
           !
           call tem_turbul(&
                pnode,pgaus,1_ip,pgaus,gpcon,gpsph,gpdif,dummr,gpden,gptur)

           !
           ! Reaction term
           !
           call tem_elmrea( &
                1_ip,pnode,pgaus,1_ip,pgaus,elvel,gpden,gpcar,&
                gprea)

           !
           ! Loop over Gauss points
           !
           gauss_points: do igaub=1,ngaus(pblty)
              !
              ! Jacobian EUCTA
              !
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
                   bocod,baloc,eucta)
              gbsur = elmar(pblty)%weigp(igaub)*eucta
              if( itask == 1_ip ) then
                 do inodb = 1,pnodb
                    elmas(inodb) = elmas(inodb) + gbsur * elmar(pblty) % shape(inodb,igaub)
                 end do
              end if
              call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
              !
              ! Cylindrical coordinates
              !
              if( kfl_naxis == 1 ) then
                 gpcod = 0.0_rp
                 do inodb = 1,pnodb
                    gpcod = gpcod + bocod(1,inodb) * elmar(pblty) % shape(inodb,igaub)
                 end do
                 gbsur = gbsur * gpcod * twopi
              end if
              !
              ! Robin: a.n
              !
              
              if( kfl_advec_tem == 1 ) then
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun)
                    bovel(1:ndime,inodb)=advec(1:ndime,ipoin,1)
                 end do
              else if(kfl_advec_tem>=2) then
                 call tem_velfun(pnodb,bocod,bovel)
              end if

              gbvel=0.0_rp
              do inodb=1,pnodb
                 gbvel(1:ndime)=gbvel(1:ndime)&
                      +bovel(1:ndime,inodb)*elmar(pblty)%shape(inodb,igaub)
              end do
              adotn = dot_product(gbvel(1:ndime), baloc(1:ndime, ndime))

              gbtem = 0.0_rp
              do inodb = 1,pnodb                  
                 gbtem = gbtem + botem(inodb) *elmar(pblty) % shape(inodb,igaub)
              end do
              qrobi = 0.0_rp
              arobi = 0.0_rp
              trobi = 0.0_rp 

              if( kfl_fixbo_tem(iboun) == 2 ) then
                 !
                 ! Robin condition
                 ! k*grad(T).n = qr+ar*(T-Tr)
                 !
                 arobi = bvnat_tem(1,iboun,1)                        ! ar
                 trobi = bvnat_tem(2,iboun,1)                        ! Tr   
                 qrobi = bvnat_tem(3,iboun,1)                        ! qr    
                 
              else if( kfl_fixbo_tem(iboun) == 4 ) then
                 !
                 ! Augmented Robin condition
                 ! k*grad(T).n = P3+P1*T+P2*(exp(T/P4)-1.0)
                 !
                 para1 = bvnat_tem(1,iboun,1)                        ! P1
                 para2 = bvnat_tem(2,iboun,1)                        ! P2   
                 para3 = bvnat_tem(3,iboun,1)                        ! P3    
                 para4 = bvnat_tem(4,iboun,1)                        ! P4 

                 qrobi = -para3-para2*(exp(gbtem/para4)-1.0_rp)
                 arobi = -para1
                 trobi = 0.0_rp 

              else if( kfl_fixbo_tem(iboun) == 3 ) then
                 !
                 ! Law of the wall
                 !              
                 ! k*grad(T).n = -(rho*cp*u/T+)*(T-Tw)
                 !                  

                 ! CRITERIA
                 ! CODES BOUNDARIES
                 !   BOUNDCODE   CODE(3)   VALUE    IFIEL
                 ! END_CODES
                 ! BVNAT(1)= VALUE  
                 ! BVNAT(2) =IFIEL (from which field)


                 ifiel = nint(bvnat_tem(2, iboun,1))
                 if (ifiel==0) then ! constant value for wall temper
                    twall = bvnat_tem(1, iboun,1)
                 else  ! time dependent wall temperature
                    twall = 0.0_rp
                    if (kfl_tendencies_ker.or..true.) then ! reads wall temper from a  discrete funcion
                       call ker_functions(1_ip, ifiel, FUNCTION_DISCRETE, 1.0_rp, twall)
!                       if (iboun==1.and.igaub==1) print *, 'twall', cutim, twall
                    else ! reads from a field
                       kk = k_tran_fiel(ifiel) !indicates to which interval the current time belongs.
                       xx = x_tran_fiel(ifiel) !indicates the position between the begining and end of the interval. 
                       do inodb =1, pnodb
                          ipoin = lnodb(inodb,iboun)
                          temwa = xfiel(ifiel) % a(1,ipoin,kk) * xx + xfiel(ifiel) % a(1,ipoin,kk+1) * (1.0_rp-xx)
                          twall= twall + temwa * elmar(pblty)%shape(inodb,igaub)
                       end do
                    end if
                 end if
#ifdef GABLS1
                 twall =  265.0 - 6.9444e-05*cutim ! gabls1
#endif
                 if( kfl_rough > 0 )then
                    roughness = 0.0_rp
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)                
                       roughness=roughness+ rough(ipoin)*&
                            elmar(pblty)%shape(inodb,igaub)
                    end do
                 else
                    roughness = rough_dom
                 end if

                 kinen = 0.0_rp
                 if( kfl_ustar == 2 ) then 
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       kinen = kinen + untur(1,ipoin,1) * elmar(pblty)%shape(inodb,igaub)
                    end do
                    if (kfl_regularization) & ! then
                         kinen = regul_k(kinen)
                 end if

                 if ( kfl_waexl_ker == 1_ip ) then !if exchange location for wall law
                    temex = temel_ker(lexlo_ker(igaub,iboun))
                    velex(1:ndime) = velel_ker(1:ndime,lexlo_ker(igaub,iboun))
                 else
                    temex = 0.0_rp
                    velex(1:ndime) = 0.0_rp
                 end if

                 if ( kfl_wlaav_ker == 1_ip ) then ! Time-averaged velocity for wall law
                    veave(1:ndime) = velav_ker(1:ndime,igaub,iboun)
                 else
                    veave(1:ndime) = 0.0_rp
                 end if
                 call tem_bouwal(&
                      lboel(1,iboun),iboun, elmar(pblty)%shape(1,igaub),pnodb,&
                      pnode,gbden(igaub),gbvis(igaub),&
                      gbsph(igaub),gbcon(igaub),prtur_tem,twall,roughness,&
                      kinen,baloc,velex,veave,arobi,trobi,qrobi, gbvel, ustar)
!                 write(666,'(5(f15.5))')  cutim, twall, -arobi, therm(1,1), ustar

              else if( kfl_fixbo_tem(iboun) == 5 ) then
                 !
                 ! Robin condition for water vapor model
                 ! 
                 call tem_watvap(&
                      pnodb,ndime,baloc(1,ndime),bogrc,&
                      elmar(pblty) % shape(1,igaub),&
                      gbtem,qrobi,arobi,trobi)

              elseif( kfl_fixbo_tem(iboun) == 20 ) then
                 !
                 ! Stable outflow condition (implicit)
                 ! beta*rho * {u · n}_{-} (gptem - trobi)
                 !
                 ! boundary velocity
                 beta = bvnat_tem(1,iboun,1)  
                 trobi = bvnat_tem(2,iboun,1)
                 !
                 ! implicit boundary term to matrix
                 !                 
                 arobi =  min(adotn, 0.0_rp)*gbden(igaub)*gbsph(igaub)*beta
              
              end if
              rhocp = gbsph(igaub)*gbden(igaub)
              if (  kfl_fixbo_tem(iboun) == 3.and.&
                   kfl_waexl_imp_ker == 0 .and.kfl_waexl_ker==1) then
                 ! IF EXPLICIT EXCHANGE LOCATION all to rhs (arobi*temex to rhs)
                 xmrhs =  qrobi - arobi*trobi                         !  qr - ar * Tr
                 xmrhs =  xmrhs + (arobi - bemol_tem*adotn*rhocp)*temex     ! -ar + bemol(u.n)
                 xmmat =  0.0_rp 
              else
                 xmrhs =  qrobi - arobi*trobi             !  qr - ar * Tr
                 xmmat =- arobi + bemol_tem*adotn*rhocp   ! -ar + bemol(u.n)
              end if
              call tem_boumat(&
                   pnode,pnodb,igaub, iboun,lboel(1,iboun),lelbo(iboun),xmmat,xmrhs,&
                   elmar(pblty)%shape(1,igaub),gbsur,elmat,elrhs)

              if (itask==2) then ! boundary heat flux and surface
                 if (kfl_waexl_ker==0) then
                    temwa = gbtem
                 else
                    temwa = temex
                 end if
                 heatf_tem(iboun) = heatf_tem(iboun) + (xmrhs - xmmat*temwa)*gbsur ! boundary heat flux
                 massb_tem(iboun) = massb_tem(iboun) + gbsur                       ! boundary surface
              end if
              
           end do gauss_points
           !
           ! Prescribe Dirichlet boundary conditions
           !
           if( solve(1) % kfl_iffix == 0 ) &
                call tem_elmdir(&
                pnode,lnods(1,ielem),elmat,elrhs,ielem)
           !
           ! Assembly
           !
           !call assrhs(solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)
           call solver_assemble_element_matrix_scalar(&
                solve_sol(1),1_ip,pnode,pnode,ielem,lnods(:,ielem),elmat,amatr)

           !call assmat(&
           !     solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
           !     ielem,lnods(1,ielem),elmat,amatr)

        end if

     end if Neumann_or_Robin

  end do boundaries

end subroutine tem_bouope
