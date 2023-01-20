!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_updhfl()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_tothfl
  ! NAME 
  !    tem_tothfl
  ! DESCRIPTION
  !    This routine computes the heat injected from the boundary
  !    for Low-Mach model:
  !             _     
  !            |               
  !    TFLUX = |  k*gradt(T).n 
  !           _|S       
  !
  ! USES
  ! USED BY
  !    tem_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper
  use mod_bouder
  use mod_tem_turbul
  implicit none
  real(rp)         :: baloc(ndime,ndime)
  real(rp)         :: gptem(mgaus)
  real(rp)         :: bocod(ndime,mnodb)
  real(rp)         :: elcod(ndime,mnode),elgrt(ndime,mnode)
  real(rp)         :: eltem(mnode)
  integer(ip)      :: ielem,inode,ipoin,idime
  integer(ip)      :: igaus,igaub,iboun,inodb,pblty
  integer(ip)      :: pnodb,pmate,pnode,pelty,pgaus,dummi
  real(rp)         :: eucta,tmatr,gbsur,gbgrt(3),gbcon
  real(rp)         :: gpsph(mgaus),gpden(mgaus),gpdif(mgaus)
  real(rp)         :: gprea(mgaus),xjaci(9),xjacm(9)
  real(rp)         :: gpcar(ndime,mnode,mgaus)
  real(rp)         :: gpcon(mgaus),dummr(ndime*mnode),gpcod,detjm
  real(rp)         :: gptur(mgaus)
  real(rp)         :: gpgrt(ndime,mgaus)
  real(rp), target :: rtflu(1)

  if(kfl_regim_tem>=3) then

     rtflu(1)=0.0_rp

     if(kfl_paral/=0) then
        !
        ! Loop over elements  
        !
        boundary_nodes: do iboun=1,nboun

           pblty = ltypb(iboun)
           pnodb = nnode(pblty)
           ielem = lelbo(iboun)
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           pmate = 1
           if(nmate>1) pmate=lmate(ielem)
           !
           ! Gather operations
           !
           do inodb=1,pnodb
              ipoin=lnodb(inodb,iboun)
              do idime=1,ndime
                 bocod(idime,inodb) = coord(idime,ipoin)
              end do
           end do
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
              eltem(inode) = therm(ipoin,1)
           end do
           !
           ! GPTEM+GPSGS: Temperature at Gauss point
           !
           call gather(&
                1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
                elmar(pelty)%shape,therm,gptem)
           !
           ! Cartesian derivatives
           !
           do igaus=1,pgaus
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci)        ! and Jacobian
           end do
           !
           ! Properties
           !
           
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('CONDU','PGAUS',dummi,ielem,gpcon,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty)%shape,gpcar)
           ! 
           ! Coupling with turbul
           !
           call tem_turbul(&
                pnode,pgaus,1_ip,pgaus,gpcon,gpsph,gpdif,dummr,gpden,gptur)
                 
           !
           ! Reaction term
           !
           call tem_elmrea( &
                1_ip,pnode,pgaus,1_ip,pgaus,dummr,gpden,gpcar,&
                gprea)

           !
           ! Loop over Gauss points
           !
           gauss_points: do igaub=1,ngaus(pblty)
              !
              ! Diffusion
              !
              gbcon=0.0_rp
              do igaus=1,pgaus
                 do inodb=1,pnodb                  
                    tmatr=elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                         *elmar(pblty)%shape(inodb,igaub)
                    gbcon=gbcon+gpdif(igaus)*tmatr
                 end do
              end do
              !
              ! Jacobian EUCTA
              !
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
                   bocod,baloc,eucta)
              gbsur=elmar(pblty)%weigp(igaub)*eucta 
              call chenor(pnode,baloc,bocod,elcod)                      ! Check normal                 
              !
              ! Cylindrical coordinates
              !
              if(kfl_naxis==1) then
                 gpcod=0.0_rp
                 do inodb=1,pnodb
                    gpcod=gpcod+bocod(1,inodb)*elmar(pblty)%shape(inodb,igaub)
                 end do
                 gbsur=gbsur*gpcod*twopi
              end if
              !
              ! GBGRT: grad(T) on boundary Gauss point
              !
              elgrt=0.0_rp
              gbgrt=0.0_rp
              gpgrt=0.0_rp
              do igaus=1,pgaus
                 do inode=1,pnode
                    do idime=1,ndime
                       gpgrt(idime,igaus)=gpgrt(idime,igaus)&
                            +gpcar(idime,inode,igaus)*eltem(inode)
                    end do
                 end do
              end do
              do igaus=1,pgaus
                 do inode=1,pnode
                    ipoin=lnods(inode,ielem)
                    do idime=1,ndime
                       elgrt(idime,inode)=elgrt(idime,inode)&
                            +gpgrt(idime,igaus)*elmar(pelty)%shaga(igaus,inode) 
                    end do
                 end do
              end do
              do inodb=1,pnodb
                 inode=lboel(inodb,iboun)
                 do idime=1,ndime
                    gbgrt(idime)=gbgrt(idime)&
                         +elgrt(idime,inode)*elmar(pblty)%shape(inodb,igaub)
                 end do
              end do
              !
              ! TFLUX: Gauss point contribution
              !
              do idime=1,ndime
                 rtflu(1)=rtflu(1)+gbsur*gbgrt(idime)*baloc(idime,ndime)*gbcon
              end do

           end do gauss_points

        end do boundary_nodes
        !
        ! int_V Q dV
        !
        if(kfl_sourc_tem /= 0) then
           call runend('TEM_UPDHFL: NOT CODED')
           !!!!rtflu(1)=rtflu(1)+vodom*sourc_tem(1)
        end if

     end if
     !
     ! Parall: Sum over all subdomains
     !
     if(kfl_paral>=0) then
        nparr =  1
        parre => rtflu
        call par_operat(3_ip)
     end if
     tflux=rtflu(1)

  end if

end subroutine tem_updhfl

