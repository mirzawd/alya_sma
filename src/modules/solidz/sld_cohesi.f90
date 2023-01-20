!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_cohesi(itask)
  !-----------------------------------------------------------------------
  !****f* solidz/sld_cohesi
  ! NAME
  !    sld_cohesi
  ! DESCRIPTION
  !    This routine computes the cohesi and right hand side
  ! USES
  !    sld_elmope
  !    sld_bouope
  !    sld_assloa
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_solidz
  use mod_elmgeo, only  : elmgeo_natural_coordinates, elmgeo_natural_coordinates_on_boundaries 
  use mod_bouder
  implicit none
  integer(ip) :: itask,ielem,pblty,pnodb,pgaub,igaub,ifoun,ptopo,pelty
  integer(ip) :: pnode,iboun,inode,idime,inodb,ipoin,jdime,jnode
  integer(ip) :: idofn,jdofn
  real(rp)    :: eucta,baloc(9),gbsur
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: tract(ndime),stiff(ndime,ndime),xdisp(ndime)
  real(rp)    :: lmini,lmaxi
!  real(rp)    :: toler
  real(rp)    :: coloc(3),coglo(3)
  real(rp)    :: deri1(ndime,mnode)
  real(rp)    :: shap1(mnode)
  real(rp)    :: elrhs(ndofn_sld*mnode),elmat(ndofn_sld*mnode,ndofn_sld*mnode)
  
  logical     :: debugging
  
  if( kfl_cohes_sld > 0 .and. kfl_elcoh > 0 .and. kfl_xfeme_sld == 1 .and. INOTMASTER ) then

     if( ndime == 2 ) then
        pblty = BAR02
     else
        pblty = TRI03
     end if
     pnodb = nnode(pblty)
     pgaub = ngaus(pblty)
     lmini = 0.0_rp-zeror
     lmaxi = 1.0_rp+zeror

     do ielem = 1,nelem

        pelty = ltype(ielem)
        debugging = .false. !((ielem==916 .or. ielem==935 .or. ielem==944 .or. ielem==975) .and. ittim == 30)
 
        if( lelch(ielem) == ELCUT .and. pelty > 0 .and. &
           (leenr_sld(ielem) /= 0 ) ) then
           !
           ! Element
           !
           pelty = ltype(ielem)
           pnode = lnnod(ielem)
           ptopo = ltopo(pelty)
           elrhs = 0.0_rp
           do inode = 1,lnnod(ielem)
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           !
           ! Initialize
           !
           do idofn=1,ndofn_sld*mnode
              elrhs(idofn)=0.0_rp
           end do
           
           if (itask == 2) then
              do idofn=1,ndofn_sld*mnode
                 do jdofn=1,ndofn_sld*mnode
                    elmat(idofn,jdofn)=0.0_rp
                 end do
              end do
           end if
           !
           ! Boundaries
           !
           do iboun = 1,cutel(ielem) % nboun

              do inodb = 1,pnodb
                 do idime = 1,ndime
                    bocod(idime,inodb) = cutel(ielem) % lb(iboun) % bocod(idime,inodb)
                 end do
              end do

             
              do igaub = 1,pgaub
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),&
                      bocod,baloc,eucta)    
                 gbsur = elmar(pblty) % weigp(igaub) * eucta 
                 !
                 ! COGLO: absolute coordinate of IGAUB
                 !
                 do idime = 1,ndime
                    coglo(idime) = 0.0_rp
                 end do
                 
                 do inodb = 1,pnodb
                    do idime = 1,ndime
                       coglo(idime) = coglo(idime) + bocod(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                    end do
                 end do
                 !
                 ! Evaluate test function SHAP1 at gauss point IGAUB
                 !
                 call elmgeo_natural_coordinates(      &
                      ndime,pelty,pnode,elcod,shap1,   &
                      deri1,coglo,coloc,ifoun)                

                 !call elsest_chkelm(&
                 !     ndime,ptopo,pnode,elcod,shap1,deri1,&
                 !     coglo,coloc,ifoun,lmini,lmaxi)

                 !toler = 0.01_rp
                 !ptopo     = ltopo(pblty)

                 !call elmgeo_natural_coordinates_on_boundaries(&
                 !     ndime,pblty,pnodb,bocod,         &
                 !     shap1,deri1,coglo, &
                 !     coloc,ifoun,toler)


                 ! (ifoun == 0) call runend('SLD_COHESI: ELEMENT NOT FOUND')
                 !
                 ! Compute the jump in displacement (opening) at IGAUB
                 !
                 xdisp  = 0.0_rp
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    jdime = 0_ip
                    do idime = 1,ndime
                       jdime = jdime + 1
                       xdisp(jdime)  = xdisp(jdime)  + dxfem_sld(idime,ipoin,1) * shap1(inode)
                    end do
                 end do
                 !
                 ! Get the cohesive traction and tangent stiffness from cohesive law
                 !
                 tract = 0.0_rp
                 stiff = 0.0_rp
                 call sld_cohpre(itask,ielem,iboun,igaub,xdisp,tract,stiff)
                 !
                 ! Construct ELRHS based on the cohesive traction TRACT
                 !
                 do inode = 1,pnode
                    idofn = (inode-1)*ndofn_sld + ndime
                    do idime = 1,ndime
                       idofn = idofn + 1
                       elrhs(idofn) = elrhs(idofn) - gbsur * tract(idime) * shap1(inode)
                    end do
                 end do
                 !
                 ! Construct ELMAT based on the cohesive tangent stiffness STIFF
                 !
                 if (itask == 2) then
                    do inode = 1,pnode
                       idofn = (inode-1)*ndofn_sld + ndime
                       do idime = 1,ndime
                          idofn = idofn + 1
                          do jnode = 1,pnode
                             jdofn = (jnode-1)*ndofn_sld + ndime
                             do jdime = 1,ndime
                                jdofn = jdofn + 1
                                elmat(idofn,jdofn) = elmat(idofn,jdofn) + gbsur * stiff(idime,jdime) &
                                                                          * shap1(inode) * shap1(jnode)
                             end do
                          end do
                       end do
                    end do
                 end if

              end do 
              
           end do
           !
           ! Assembly ELRHS into global RHSID
           !
           call assrhs(ndofn_sld,pnode,lnods(1,ielem),elrhs,rhsid)
           if (itask == 2) then
              call assmat(&
                   solve(1)%ndofn,pnode,ndofn_sld*pnode,solve(1)%nunkn,&
                   ielem,solve(1)%kfl_algso,lnods(1,ielem),elmat,amatr)
           end if

           if (debugging) then
              write(*,*)''
              write(*,*)'ielem = ',ielem
              write(*,*)'elmat and elrhs (in sld_cohesi)'
              do idofn = 1,ndofn_sld*pnode
                write(*,'(100(1x,e10.3))')(elmat(idofn,jdofn),jdofn=1,ndofn_sld*pnode),elrhs(idofn)
              enddo
           end if

        end if
     end do

  end if
 
end subroutine sld_cohesi
