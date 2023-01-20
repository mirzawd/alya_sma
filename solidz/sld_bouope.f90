!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_bouope.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine computes boundary operations
!> @details This routine computes boundary operations
!> @}
subroutine sld_bouope()
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_bouope
  ! NAME
  !    sld_bouope
  ! DESCRIPTION
  !    boundary operations:
  !    1. Compute elemental matrix and RHS
  !    2. Impose Dirichlet boundary conditions
  !    3. Assemble them
  ! USES
  !    sld_bouave
  !    sld_bougat
  !    sld_elmpro
  !    bouder
  !    chenor
  !    sld_bouwal
  !    elmder
  !    cartbo
  !    sld_bouopb
  !    sld_elmdir
  !    sld_assrhs
  !    assmat
  ! USED BY
  !    sld_matrix
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solidz
  use mod_sld_cardiac_cycle, only : cavities
  use mod_bouder

#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet, only: kfl_sysnet,cavities_sysnet
#endif

  implicit none

  integer(ip) :: ielem,ipoin,idime,jdime,iboun,igaub,inodb,idofn,icavi
  integer(ip) :: pelty,pblty,pnodb,pgaub,pnode,inode,kfixno,itott
  real(rp)    :: elrhs(ndofn_sld*mnode),elfex(ndofn_sld*mnode)
  real(rp)    :: baloc(ndime,ndime),gppush(ndime,ndime),system_pressure
  real(rp)    :: bocod(ndime,mnodb),elcod(ndime,mnode),elpre(mnodb),elpush(ndime,ndime,mnodb),unr,xnr,fnr,tract_ref(ndime)
  real(rp)    :: gbsur,eucta,tract(ndime),kspring,muvisco,gppre,unormal,xnormal,fnormal,tnormal,xnormPF,xtractPF(ndime)

  
  !
  ! Loop over boundaries
  !
  if (kfl_icodb .ne. 0) then

     pressEndoIntegralNew_sld   = 0.0_rp
     pushForwardIntegralNew_sld = 0.0_rp
     
     boundaries: do iboun=1,nboun

        if( ANY( kfl_fixbo_sld(iboun)== (/2_ip,3_ip,6_ip,8_ip,9_ip,10_ip,11_ip/) ) ) then
           !
           ! fixbo=2  --> Nodal pressure (in press)
           ! fixbo=3  --> Boundary element pressure (in bvnat, as a traction on a LOCAL basis)
           ! fixbo=5  --> Nodal traction coming from a field (in bvess)
           ! fixbo=6  --> Boundary element traction (in bvnat, as a traction on a GLOBAL basis)
           ! fixbo=8  --> Damper-like boundary condition
           ! fixbo=9  --> Boundary element Windkessel pressure (in bvnat, as a traction on a LOCAL basis)
           ! fixbo=10 --> Damper-like boundary condition + Cardiac Cycle
           ! fixbo=11 --> Quarte's boundary condition for the cardiac base
           !
           ! Element properties and dimensions
           !
           pblty=ltypb(iboun)
           pnodb=nnode(pblty)
           pgaub=ngaus(pblty)
           ielem=lelbo(iboun)
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           !
           ! Initialize
           !
           elrhs    = 0.0_rp
           elfex    = 0.0_rp
           elpush   = 0.0_rp
           elpre    = 0.0_rp
           !
           ! Gather operations: BOCOD, ELCOD, ELPUSH
           !
           if( kfl_follo_sld == 0_ip ) then
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 do idime=1,ndime
                    bocod(idime,inodb) = coord(idime,ipoin)
                    elpre(inodb)       = press(ipoin,1_ip)
                    if( kfl_gdepo /= 0 ) then
                       do jdime=1,ndime
                          elpush(idime,jdime,inodb)= gdeinv(idime,jdime,ipoin)*gdedet(ipoin)
                       end do
                    end if
                 end do
              end do
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do
           else
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 bocod(1:ndime,inodb) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K)
                 elpre(inodb)         = press(ipoin,1_ip)
              end do
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K) 
              end do
           end if
           !
           ! Loop over Gauss points
           !
           gauss_points: do igaub=1,pgaub
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
                   bocod,baloc,eucta)
              gbsur=elmar(pblty)%weigp(igaub)*eucta
              call chenor(&                                           ! Check normal
                   pnode,baloc,bocod,elcod)

              !
              ! Compute push forward operator in the gauss point: F^{-T}*J
              !
              gppush= 0.0_rp
              gppre= 0.0_rp
              tract = 0.0_rp
              if( kfl_gdepo /= 0 ) then
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    gppre = gppre + elmar(pblty)%shape(inodb,igaub) * elpre(inodb)
                    do idime=1,ndime
                       do jdime=1,ndime
                          gppush(idime,jdime) = gppush(idime,jdime) + elmar(pblty)%shape(inodb,igaub) * elpush(idime,jdime,inodb)
                       end do
                    end do
                 end do
              end if

              if(kfl_fixbo_sld(iboun)==2) then
                 !
                 ! pressure force(igaus) = n_i * p(igaus) * gbsur
                 !                       = baloc(1:ndime,ndime) * p(igaus) * gbsur
                 !
                 ! baloc(1:ndime,ndime) --> exterior normal vector
                 !
                 gppre= 0.0_rp
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    gppre = gppre + elmar(pblty)%shape(inodb,igaub) * elpre(inodb)
                 end do

                 tract(1:ndime) = baloc(1:ndime,ndime) * gppre

              else if (kfl_fixbo_sld(iboun)==3) then
                 !
                 ! Pressure (Hydrostatic)
                 !
                 tract(1:ndime) = baloc(1:ndime,ndime) *bvnat_sld(1,iboun,1)

              else if (kfl_fixbo_sld(iboun)==6) then
                 !
                 ! Traction 
                 !
                 tract(1:ndime) = bvnat_sld(1:ndime,iboun,1)

              else if (kfl_fixbo_sld(iboun)==8) then
                 !
                 !  normal damper-like forces
                 !  F = -k x + mu dx/dt
                 kspring= - bvnat_sld(1,iboun,2)
                 muvisco= - bvnat_sld(2,iboun,2)
         
                 tnormal = normal_damper_like_forces(kspring, muvisco)
                 tract(1:ndime) = baloc(1:ndime,ndime) * tnormal

              else if (kfl_fixbo_sld(iboun)==9) then

                 icavi= int(bvnat_sld(1,iboun,2), kind=ip)

                 system_pressure= cavities(icavi) % pres(ITER_K)     

#ifndef PROPER_ELEM_PRIVATE_OFF
                 if (kfl_sysnet) system_pressure= cavities_sysnet(icavi) % pres(TIME_N)     
#endif

                 tract(1:ndime) = - baloc(1:ndime,ndime) * system_pressure ! minus, because normal is defined pointing outwards

                 !
                 ! Compute the pressure forces push forward contribution to the integral Int_{ref} p F^{-T}*J N dA
                 ! This is used when fixbo=11 is present, i.e. the quarte's cardiac base boundary condition
                 ! It is computed now and used in the next iteration, after being PAR_SUMmed and updated in sld_matrix
                 !
                 tract_ref = - tract
                 xtractPf  = 0.0_rp
                 do idime= 1,ndime
                    do jdime= 1,ndime
                       xtractPF(idime)= xtractPF(idime)+gppush(jdime,idime)*tract_ref(jdime)
                    end do
                    pressEndoIntegralNew_sld(idime)= pressEndoIntegralNew_sld(idime) + xtractPF(idime) 
                 end do

              else if (kfl_fixbo_sld(iboun)==10) then
                 !
                 !  normal damper-like forces
                 !  F = -k x + mu dx/dt
                 !  CODES, BOUNDARIES
                 !  code 10 cavity_code (as in fixbo 9) spring_k, mu_viscosity (as in fixbo 8)     
                 !
                  icavi= int(bvnat_sld(1,iboun,2), kind=ip)
                  kspring= - bvnat_sld(2,iboun,2)
                  muvisco= - bvnat_sld(3,iboun,2)
                  
                  tnormal = normal_damper_like_forces(kspring, muvisco)
                  tract(1:ndime) = baloc(1:ndime,ndime) * tnormal - baloc(1:ndime,ndime) * cavities(icavi) % pres(ITER_K)     ! minus, because normal is defined pointing outwards

              else if (kfl_fixbo_sld(iboun)==11) then

                 !
                 ! Compute the norm of the push forward of the exterior normal
                 !
                 xnormPF  = 0.0_rp
                 xtractPF = 0.0_rp
                 do idime= 1,ndime
                    do jdime= 1,ndime
                       xtractPF(idime)= xtractPF(idime)+gppush(jdime,idime)*baloc(idime,ndime)
                    end do
                    xnormPF = xnormPF + xtractPF(idime)*xtractPF(idime)
                 end do
                 xnormPF = sqrt(xnormPF)
                 pushForwardIntegralNew_sld = pushForwardIntegralNew_sld + xnormPF

                 tract(1:ndime) = 0.0_rp
                 if (pushForwardIntegralOld_sld > 0.0_rp) then
                    tract(1:ndime)= xnormPF * pressEndoIntegralOld_sld(1:ndime) / pushForwardIntegralOld_sld
                 end if
                 
              else if (kfl_fixbo_sld(iboun)==12) then
                  call runend("Boundary condition code 12 is not programmed")
              end if  !kfl_fixbo_sld(iboun)

              do inodb = 1,pnodb
                 idofn = (lboel(inodb,iboun)-1)*ndime
                 ipoin = lnodb(inodb,iboun)
                 do idime =1,ndime
                    idofn = idofn+1
                    if( kfl_fixno_sld(idime,ipoin) == 0 .and. kfl_fixrs_sld(ipoin) == 0 ) then
                       ! Newmann part
                       elfex(idofn) = elfex(idofn) + gbsur*tract(idime) * elmar(pblty)%shape(inodb,igaub)
                    end if
                 end do
              end do

           end do gauss_points
           !
           ! Assembly
           !
           call assrhs(&
                ndofn_sld,pnode,lnods(1,ielem),elfex,rhsid)
           call assrhs(&
                ndofn_sld,pnode,lnods(1,ielem),elfex,fexte_sld)

        else if ( kfl_fixbo_sld(iboun)>0 ) then
            call runend("Unknown solidz boundary condition code "//trim(intost(kfl_fixbo_sld(iboun)))//" on boundary "//trim(intost(iboun)))
        end if  !kfl_fixbo_sld(iboun)==2 .or. 3 .or. kfl_codbo(iboun) /= 100

     end do boundaries

  end if


contains 
   real(rp) function  normal_damper_like_forces(kspring, muvisco)
      implicit none

      real(rp), intent(in)  :: kspring, muvisco 

                 !
                 !  normal damper-like forces
                 !  F = -k x + mu dx/dt
      unormal= 0.0_rp
      xnormal= 0.0_rp
      fnormal= 0.0_rp
      kfixno = 0_ip
      do inodb = 1,pnodb
         unr= 0.0_rp
         xnr= 0.0_rp
         fnr= 0.0_rp
         ipoin = lnodb(inodb,iboun)
         itott = (ipoin             -1)*ndime
         idofn=  (lboel(inodb,iboun)-1)*ndime
         if (kfl_fixno_sld(1,ipoin) == 2) kfixno = 2
         do idime=1,ndime
            idofn=idofn+1
            itott=itott+1
            fnr= fnr + baloc(idime,ndime)*rhsid(itott) ! forces nodal normal projection
            unr= unr + baloc(idime,ndime)*veloc_sld(idime,ipoin,ITER_K)  ! velocity nodal normal projection
            xnr= xnr + baloc(idime,ndime)*displ(idime,ipoin,ITER_K)  ! displacement nodal normal projection
         end do
         unormal = unormal + elmar(pblty)%shape(inodb,igaub) * unr ! velocity normal projection at the gauss point
         xnormal = xnormal + elmar(pblty)%shape(inodb,igaub) * xnr ! displacement normal projection at the gauss point
         fnormal = fnormal + elmar(pblty)%shape(inodb,igaub) * fnr ! force normal projection at the gauss point
      end do

      ! tests:
      ! if (fixed normal displacement for at least one node) then check fnormal
      !    if (fnormal < 0.0) then
      !       free node normal displacement
      !       compute tnormal and add it to elrhs
      !    else
      !       keep node normal displacement fixed
      !       keep tnormal = 0
      !    end if
      ! else
      !    if (xnormal < 0.0) then
      !       keep node normal displacement free
      !       compute tnormal and add it to elrhs
      !    else
      !       fix node normal displacement
      !       keep tnormal = 0
      !    end if
      ! end if

      tnormal= 0.0_rp

      if (kfixno == 2) then
         !if (fnormal .le. 0.0_rp) then
            kfixno = 0
            tnormal = kspring*xnormal - muvisco * unormal !normal traction component
!            tnormal = kspring*xnormal + muvisco * abs(unormal) !normal traction component
         !end if
      else
         !if (xnormal .le. 0.0_rp) then
            kfixno = 0
            tnormal = kspring*xnormal - muvisco * unormal !normal traction component
!            tnormal = kspring*xnormal + muvisco * abs(unormal) !normal traction component
         !else
         !   kfixno = 2
         !end if
      end if

      do inodb=1,pnodb
         ipoin=lnodb(inodb,iboun)
         if(kfl_fixno_sld(1,ipoin)/=1) then
            kfl_fixno_sld(1,ipoin) = kfixno
         end if
      end do

      ! tnormal = 100000.0_rp   ! debug
      !! leave this commented: it is to check that all is ok: tnormal positive means force outwards (i.e. following normal)

      !
      ! if the node moves against the normal (i.e. inwards), xnormal is negative. then
      ! a kspring negative produces a tnormal positive.
      ! then, a force in the normal direction is generated to compensate movement
      !
      ! total traction, recall that baloc(1:ndime,ndime) is the normal versor

      normal_damper_like_forces = tnormal
   end function normal_damper_like_forces

end subroutine sld_bouope
