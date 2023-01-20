!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_bouset(ibsec,ibset) 
  !-----------------------------------------------------------------------
  !****f* Temper/tem_bouset
  ! NAME 
  !    tem_bouset
  ! DESCRIPTION
  !    This routine computes variables on a boundary set W.
  !    The variable are:
  !    1. setsu: set surface          =  meas(W)=int_W
  !    2. setmt: set mean temperature =  int_W T/meas(W)
  !    3. sethf: set heat flux        =  int_W k*grad(T).n
  !    4. sethn: set heat flux        =  variational
  !    The heat flux represent the flo from the solid to the fluid.
  !    Positive value means incomming positive heat.
  ! USES
  !    bouder
  !    chenor
  ! USED BY
  !    tem_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper 
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_bouder
  use mod_tem_turbul
  implicit none

  integer(ip), intent(in)  :: ibsec
  integer(ip), intent(in)  :: ibset
  real(rp),    pointer     :: setsu(:),setmt(:),sethf(:),sethn(:),setde(:)
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: gbcar(ndime,mnode)
  real(rp)                 :: xjaci(ndime,ndime) 
  real(rp)                 :: xjacm(ndime,ndime) 
  integer(ip)              :: ielem,inode,ipoin
  integer(ip)              :: igaus,idime,igaub,iboun,inodb,pblty
  integer(ip)              :: pnodb,pmate,pnode,pelty,pgaus,pgaub,dummi,nn
  real(rp)                 :: eucta,tmatr,dsurf,detjm
  real(rp)                 :: gbsph,gbden,gbdif
  real(rp)                 :: gpsph(mgaus),gpden(mgaus)
  real(rp)                 :: gpdif(mgaus),gpcon(mgaus)
  real(rp)                 :: gpcar(ndime,mnode,mgaus)
  real(rp)                 :: gptur(mgaus)
  real(rp)                 :: gprea(mgaus),gptem(mgaus)
  real(rp)                 :: gbtem,gbgrt(ndime)
  real(rp)                 :: dummr(ndime,mnode)

  gpcar = 0.0_rp


  if( INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! Initialization
     !
     !----------------------------------------------------------------------

     nn    =  postp(1) % nvabs + 1
     setsu => postp(1) % vbset(nn:nn,ibset)
     setmt => postp(1) % vbset(1:1,ibset)
     sethf => postp(1) % vbset(2:2,ibset)
     sethn => postp(1) % vbset(3:3,ibset)
     setde => postp(1) % vbset(4:4,ibset)



     setsu = 0.0_rp
     setmt = 0.0_rp
     sethf = 0.0_rp
     sethn = 0.0_rp
     setde = 0.0_rp
     !
     ! Loop over elements
     !
    ! if( postp(1) % npp_setsb(2) /= 0 ) then ! postprocess boundary heat flux (-k grad T.n)
     boundaries: do iboun = 1,nboun
        
        if( lbset(iboun) == ibsec ) then

           pblty = ltypb(iboun)
           pnodb = lnnob(iboun)
           pgaub = ngaus(pblty)
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
                 bocod(idime,inodb)=coord(idime,ipoin)
              end do
           end do
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do
           !
           ! Cartesian derivatives
           !
           do igaus=1,pgaus
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
                   elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci) 
           end do
           !
           ! GPTEM: Temperature at Gauss point
           !
           call gather(&
                1_ip,pgaus,pnode,1_ip,lnods(1,ielem),&
                elmar(pelty)%shape,tempe,gptem)
           !
           ! Properties 
           !

           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('CONDU','PGAUS',dummi,ielem,gpcon,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar)
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


           gauss_points: do igaub=1,pgaub
              !
              ! Properties
              !
              gbden=0.0_rp
              gbsph=0.0_rp
              gbdif=0.0_rp
              do igaus=1,pgaus
                 do inodb=1,pnodb                  
                    tmatr=elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                         *elmar(pblty)%shape(inodb,igaub)
                    gbden=gbden+gpden(igaus)*tmatr
                    gbsph=gbsph+gpsph(igaus)*tmatr
                    gbdif=gbdif+gpdif(igaus)*tmatr
                 end do
              end do

              ! Jacobian
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
                   bocod,baloc,eucta)
              call chenor(pnode,baloc,bocod,elcod)
              dsurf=elmar(pblty)%weigp(igaub)*eucta 
              setsu=setsu+dsurf

              if(postp(1)%npp_setsb(1)/=0) then
                 !
                 ! Mean temperature
                 !
                 gbtem=0.0_rp 
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun)
                    gbtem=gbtem+elmar(pblty)%shape(inodb,igaub)*tempe(ipoin,1)
                 end do
                 setmt=setmt+gbtem*dsurf
              end if

              if( postp(1) % npp_setsb(2) /= 0 ) then
                 !
                 ! Heat flux
                 !
                 do inode=1,pnode
                    ipoin=lnods(inode,ielem)
                    elcod(1:ndime,inode)=coord(1:ndime,ipoin)
                 end do
                 !
                 ! GBCAR: Cartesian derivatives at boundary Gauss point
                 !
                 call gpcabo(&
                      pnode,pgaus,pnodb,lboel(1,iboun),elmar(pelty)%shaga,&
                      gpcar,elmar(pblty)%shape(1,igaub),gbcar)
                 gbgrt = 0.0_rp
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem) 
                    do idime = 1,ndime
                       gbgrt(idime) = gbgrt(idime) + gbcar(idime,inode) * tempe(ipoin,1)
                    end do
                 end do
                 sethf = sethf + gbdif * dsurf * dot_product(baloc(1:ndime,ndime),gbgrt(1:ndime))
              end if

              if(postp(1)%npp_setsb(4)/=0) then
                 !
                 ! Mean density
                 !
                 setde=setde+gbden*dsurf
              end if

           end do gauss_points

        end if

     end do boundaries
     !  end if ! postprocess MEANH
     !
     ! Internal force (variational heatflx)
     !
     if( postp(1) % npp_setsb(3) /= 0 ) then
        if ( kfl_waexl_ker /= 0) then ! MATIAS temporal condition
           do iboun = 1,nboun
              if( lbset(iboun) == ibsec ) then
                 sethn = sethn + heatf_tem(iboun)  ! total heat flux
!                 setsu = setsu + massb_tem(iboun)               
              end if
           end do        
        else
           call memgen(1_ip,npoin,0_ip)
           do iboun = 1,nboun
              if( lbset(iboun) == ibsec ) then
                 do inodb = 1,lnnob(iboun)
                    ipoin = lnodb(inodb,iboun)
                    gisca(ipoin) = 1
                 end do
              end if
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY CODE')
           do ipoin = 1,npoin
              if( gisca(ipoin) /= 0 ) then
                 sethn = sethn - solve_sol(1) % reaction(1,ipoin) / real(gisca(ipoin),rp)
              end if
           end do
           call memgen(3_ip,npoin,0_ip)
        end if
     end if
 endif  

     
end subroutine tem_bouset
subroutine tem_boundary_heat_flux()
  !******************************************************************************
  ! returns variational heat flux over points
  !
  !******************************************************************************
  use def_kintyp, only         :  ip,rp
  use def_domain, only         :  npoin, lpoty
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use def_master, only         :  INOTMASTER
  use def_temper, only         :  heatf_tem, massb_tem
  implicit none
  
  integer (ip)        :: ipoin
 

  call tem_bouope(1_ip)  ! loop over boundaries
  
  if( INOTMASTER) then
     call PAR_INTERFACE_NODE_EXCHANGE(heatf_tem,'SUM','IN MY CODE')
     call PAR_INTERFACE_NODE_EXCHANGE(massb_tem,'SUM','IN MY CODE')
  end if

  ! obtain variational heat flux
  do ipoin = 1,npoin
     if( lpoty(ipoin)  >= 1 ) then
        heatf_tem(ipoin) = heatf_tem(ipoin) / massb_tem(ipoin)     
     else
        heatf_tem(ipoin) = 0
        ! CHECK
        !if (abs(heatf_tem(ipoin).gt.1.0d-7)) then
        !   print *, 'heatfl is not zero'
        !   
        !end if
     end if
  end do
!  call memory_deallo(mem_modul(1:2,modul),'HEATF_TEM','tem_variational_heat_flux',massb_tem)
 
  return
end subroutine tem_boundary_heat_flux
