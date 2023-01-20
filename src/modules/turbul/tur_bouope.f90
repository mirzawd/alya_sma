!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_bouope()
  !------------------------------------------------------------------------
  !****f* Turbul/tur_bouope
  ! NAME 
  !    tem_bouope
  ! DESCRIPTION
  !    ORDER=1:
  !      Turbulence equation, boundary operations
  !      used only for stable outflow
  ! USES
  ! USED BY
  !    tur_matrix 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_ker_proper 
  use def_domain
  use def_turbul, only : kfl_fixbo_tur, iunkn_tur, bvnat_tur, &
       bvess_tur, kfl_fixno_tur
!  use def_turbul, only : kfl_advec_tur
  
  use mod_solver, only : solver_assemble_element_matrix_scalar
  use mod_matrix, only : matrix_assemble_element_RHS
  use mod_memory,         only : memory_alloca, memory_deallo
  use mod_bouder
  implicit none
  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: elcod(ndime,mnode), bocod(ndime,mnodb)
!  real(rp)    :: bovel(ndime,mnodb)
  integer(ip) :: ielem,ipoin
  integer(ip) :: igaub,iboun,inodb,pblty
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  integer(ip) :: dummi
  real(rp)    :: eucta,gbsur,udotn
  real(rp)    :: gbden(mgaub),gbvel(3)
  real(rp)    :: xmrhs,xmmat
  real(rp)    :: gpcod, beta, gbess, gbes2, shapb
  integer(ip),    pointer    :: fixbc(:)
  !
  ! Loop over elements  
  nullify(fixbc)
  call memory_alloca(mem_modul(1:2,modul),'FIXBC','tur_bouope',fixbc ,npoin)

  boundaries: do iboun=1,nboun
     Neumann_or_Robin: if(kfl_fixbo_tur(iboun) == 20_ip) then  ! 
  
        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        ielem = lelbo(iboun) !lboel(pnodb+1,iboun)
        pelty = ltype(ielem)

        if( pelty > 0 ) then

           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           pmate = 1
           if( nmate > 1 ) pmate = lmate(ielem)
           !
           ! Inititalize
           !
           elmat = 0.0_rp
           elrhs = 0.0_rp
          
           bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
           !
           ! Properties 
           !
           call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
!!$           if( kfl_advec_tur == 1 ) then
!!$              do inodb=1,pnodb
!!$                 ipoin=lnodb(inodb,iboun)                 
!!$                 bovel(1:ndime,inodb)=advec(1:ndime,ipoin,1)
!!$              end do
!!$           end if
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
               
              
              xmmat = 0.0_rp  ! implicit boundary term  xmmat*gptur
              xmrhs = 0.0_rp  ! explicit boundary term  xmrhs
              if( kfl_fixbo_tur(iboun) == 20_ip ) then
                 !
                 ! Stable outflow condition
                 ! beta*rho * {u · n}_{-} gptur 
                 !
                 ! boundary velocity
                 !
                 gbvel = 0.0_rp
                 gbess = 0.0_rp
                 gbes2 = 0.0_rp
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun) 
                    shapb = elmar(pblty)%shape(inodb,igaub)
                    
                    gbvel(1:ndime)=gbvel(1:ndime)&  
                         +advec(1:ndime, ipoin,1)*shapb
                    gbess= gbess &
                         + bvess_tur(1,ipoin,iunkn_tur)*shapb  
!                    gbes2 = gbes2 &
!                         + bvcod(1+iunkn_tur)%a(1, ipoin)*shapb  
                 end do
                 !
                 ! a.n
                 !
                 udotn = dot_product(gbvel(1:ndime), baloc(1:ndime, ndime))
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun) 
!!$                    if(udotn.lt.0.0.and. kfl_fixno_tur(1,ipoin, iunkn_tur) ==0) then ! incoming flow and not dirichlet
!!$                       fixbc(ipoin) =1
!!$                       kfl_fixno_tur(1,ipoin, iunkn_tur) =1
!!$                    else if (fixbc(ipoin) == 0.and.  kfl_fixno_tur(1,ipoin, iunkn_tur) ==1.and.udotn.gt.0.0) then ! leave free only 
!!$                       ! if outcoming flow for all boundaries
!!$                       kfl_fixno_tur(1,ipoin, iunkn_tur) =0
!!$                    end if
                    if (associated(kfl_geono)) then
                       if (kfl_geono(lpoty(ipoin))==20_ip) then ! outflow bc
                          if(udotn.lt.0.0_rp) then ! incoming flow
                             fixbc(ipoin) =1
                             kfl_fixno_tur(1,ipoin, iunkn_tur) =1_ip
                          else if (fixbc(ipoin) == 0 ) then ! leave free only 
                             ! if outcoming flow for all boundaries
                             kfl_fixno_tur(1,ipoin, iunkn_tur) =0_ip
                          end if
                       end if
                    else   ! not geometrical condition
                       if(udotn.lt.0.0_rp) then ! incoming flow
                          fixbc(ipoin) =1
                          if (abs(kfl_fixno_tur(1,ipoin,iunkn_tur)) == 7_ip) then 
                             kfl_fixno_tur(1,ipoin, iunkn_tur) =7_ip 
                             bvess_tur(1, ipoin, iunkn_tur) = max(1.0e-4_rp, bvess_tur(1, ipoin, iunkn_tur))
                          else
                             kfl_fixno_tur(1,ipoin, iunkn_tur) =1
                             bvess_tur(1, ipoin, iunkn_tur) = max(1.0e-4_rp, bvess_tur(1, ipoin, iunkn_tur))
                          end if
                       else if (fixbc(ipoin) == 0_ip ) then ! leave free only 
                          ! if outcoming flow for all boundaries
                          kfl_fixno_tur(1,ipoin, iunkn_tur) =0
                          if (abs(kfl_fixno_tur(1,ipoin,iunkn_tur)) == 7_ip) then 
                             kfl_fixno_tur(1,ipoin, iunkn_tur) =-7_ip 
                          else
                             kfl_fixno_tur(1,ipoin, iunkn_tur) =0_ip
                          end if
                       end if
                    end if
                 end do
                 !
                 beta = bvnat_tur(1,iboun,iunkn_tur)             
                 !
                 ! implicit boundary term to matrix
                 !                 
!                 xmmat = -  min(udotn, 0.0_rp)*gbden(igaub)*beta
!                 xmrhs = xmmat*gbes2
                 
!                 if (abs(1.0 - gbes2/gbess).gt.0.01) then
!                    print *, gbess, gbes2
                    
!                    print*, ( bvess_tur(1,lnodb(inodb, iboun), iunkn_tur), bvcod(1+iunkn_tur) %a(1,lnodb(inodb,iboun)), inodb=1, pnodb)
 !                end if
              end if
              
              call tur_boumat(&
                   pnode,pnodb,lboel(1,iboun),xmmat,xmrhs,&
                   elmar(pblty)%shape(:,igaub),gbsur,elmat,elrhs)

           end do gauss_points
           !
           ! Prescribe Dirichlet boundary conditions
           !
           if( solve(iunkn_tur) % kfl_iffix == 0 ) &
                call tur_elmdir(&
                1_ip,pnode,pnode, lnods(:,ielem),elmat,elrhs)
           !
           ! Assembly
           !          
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)
           call solver_assemble_element_matrix_scalar(&
                solve_sol(1),1_ip,pnode,pnode,ielem,lnods(:,ielem),elmat,amatr)

        
        end if

     end if Neumann_or_Robin

  end do boundaries

  call memory_deallo(mem_modul(1:2,modul),'FIXBC','tur_bouope',fixbc)

end subroutine tur_bouope
