!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_exaerr.f90
!> @author  Guillaume Houzeaux
!> @brief   Manufactured solution handling
!> @details This routine computes the FEM errors (referred to an analytical
!>          solution defined in exacso.f90. The errors are normalized by the 
!>          appropriate norm of the exact solution, except when this norm 
!>          is zero. 
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_exaerr(itask)
  use def_parame
  use def_domain 
  use def_master
  use def_kermod
  use def_nastin
  use def_elmtyp
  use mod_ker_proper
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: gpcar(ndime,mnode) 
  real(rp)                :: xjaci(ndime,ndime) 
  real(rp)                :: xjacm(ndime,ndime) 
  real(rp)                :: elvel(ndime,mnode)                          ! Gather 
  real(rp)                :: elpre(mnode)
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: gpvol(mgaus),detjm                          ! Values at Gauss points
  real(rp)                :: gpvis(mgaus),gpden(mgaus),gppor(mgaus)
  real(rp)                :: gpgvi(ndime,mgaus)
  real(rp)                :: gpvel(ndime),grave(ndime,ndime)
  real(rp)                :: prgau,gradp(ndime),gpcod(ndime)
  real(rp)                :: exvel(ndime),exveg(ndime,ndime)
  real(rp)                :: expre,exprg(ndime),sgsno,velno
  real(rp), save          :: xerro=0.0_rp,eru02_old=0.0_rp
  real(rp)                :: difeu,abvel,diffp,abpre,dummr(3)
  real(rp)                :: erp01(2),erp02(2),erp0i(2)                  ! Pressure errors
  real(rp)                :: erp11(2),erp12(2),erp1i(2)                  ! grad(Pressure) errors
  real(rp)                :: eru01(2),eru02(2),eru0i(2)                  ! Velocity errors
  real(rp)                :: eru11(2),eru12(2),eru1i(2)                  ! grad(Velocity) errors
  integer(ip)             :: ielem,inode,ipoin,igaus,idime,jdime         ! Indices and dimensions
  integer(ip)             :: pnode,pgaus,pelty,dummi,ibopo
  integer(ip)             :: inodb,iboun
  real(rp)                :: gpdum(max(mnode,mgaus))

  if( kfl_exacs_nsi /= 0 ) then

     if( itask == 1 .and. INOTMASTER ) then

        !----------------------------------------------------------------------
        !
        ! Force dirichlet boundary condition and impose exact value on velocity
        !
        !----------------------------------------------------------------------

        if( kfl_exfix_nsi == 1 ) then 
           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 call nsi_exacso(&
                      1_ip,coord(1,ipoin),gpdum,gpdum,&
                      gpdum,gpdum,exvel,exveg,expre,exprg,&
                      dummr,dummr,dummr,dummr)
                 do idime = 1,ndime
                    kfl_fixno_nsi(idime,ipoin) = 1
                    bvess_nsi(idime,ipoin,1) = exvel(idime)
                 end do
              end if  
           end do
           !if( nodpr_nsi > 0 .and. kfl_confi_nsi >= 0 ) then
           !   call nsi_exacso(&
           !        1_ip,coord(1,nodpr_nsi),gpdum,gpdum,&
           !        gpdum,gpdum,exvel,exveg,valpr_nsi,exprg,&
           !        dummr,dummr,dummr,dummr)
           !end if
           
           call memgen(1_ip,npoin,0_ip)
           do iboun = 1,nboun
              if( kfl_fixbo_nsi(iboun) == 2 ) then
                 do inodb = 1,nnode(abs(ltypb(iboun)))
                    ipoin = lnodb(inodb,iboun)
                    gisca(ipoin) = 1
                 end do
              end if
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
           do iboun = 1,nboun
              if( kfl_fixbo_nsi(iboun) /= 2 ) then
                 do inodb = 1,nnode(abs(ltypb(iboun)))
                    ipoin = lnodb(inodb,iboun)
                    gisca(ipoin) = 0
                 end do
              end if
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MIN')
           do ipoin = 1,npoin
              if( gisca(ipoin) /= 0 ) then
                 do idime = 1,ndime
                    kfl_fixno_nsi(idime,ipoin) = 0
                 end do
              end if
           end do
           call memgen(3_ip,npoin,0_ip)

        else if( kfl_exfix_nsi == 2 ) then 
           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 call nsi_exacso(&
                      1_ip,coord(1,ipoin),gpdum,gpdum,&
                      gpdum,gpdum,exvel,exveg,expre,exprg,&
                      dummr,dummr,dummr,dummr)
                 do idime = 1,ndime
                    if( kfl_fixno_nsi(idime,ipoin) == 1 ) &
                         bvess_nsi(idime,ipoin,1) = exvel(idime)
                 end do
              end if
           end do

        end if

     else if( itask == 5 .and. INOTMASTER ) then

        !----------------------------------------------------------------------
        !
        ! Impose exact initial condition
        !
        !----------------------------------------------------------------------

        do ipoin = 1,npoin
           call nsi_exacso(&
                1_ip,coord(1,ipoin),gpdum,gpdum,&
                gpdum,gpdum,exvel,exveg,expre,exprg,&
                dummr,dummr,dummr,dummr)
           press(ipoin,1) = expre
           do idime = 1,ndime
              veloc(idime,ipoin,1) = exvel(idime)
           end do
        end do

     else if( itask == 3 .or. itask == 4 .and. INOTMASTER) then

        !----------------------------------------------------------------------
        !
        ! Compute velocity error for postprocess
        !
        !----------------------------------------------------------------------

        do ipoin = 1,npoin        
           call nsi_exacso(&
                1_ip,coord(1,ipoin),gpdum,gpdum,&
                gpdum,gpdum,exvel,exveg,expre,exprg,&
                dummr,dummr,dummr,dummr)
           if( itask == 4 ) gesca(ipoin) = abs(expre - press(ipoin,1))
           if( itask == 3 ) then
              gesca(ipoin) = 0.0_rp
              do idime = 1,ndime
                 gesca(ipoin) = gesca(ipoin) + ( exvel(idime) - veloc(idime,ipoin,1) )**2
              end do
           end if
        end do

     else if( itask /= 1 ) then

        !----------------------------------------------------------------------
        !
        ! Compute error norms of velocity and pressure
        !
        !----------------------------------------------------------------------

        erp01 = 0.0_rp
        erp02 = 0.0_rp
        erp0i = 0.0_rp
        erp11 = 0.0_rp
        erp12 = 0.0_rp
        erp1i = 0.0_rp
        eru01 = 0.0_rp
        eru02 = 0.0_rp
        eru0i = 0.0_rp
        eru11 = 0.0_rp
        eru12 = 0.0_rp
        eru1i = 0.0_rp
        gpgvi = 0.0_rp
        gpdum = 0.0_rp
        !
        ! Loop over elements to calculate global error
        !
           elements: do ielem = 1,nelem
              !
              ! Element properties and dimensions
              !
              pelty = ltype(ielem)
              if( pelty > 0 .and. lelch(ielem) /= ELEXT ) then
                 pnode = nnode(pelty)
                 pgaus = ngaus(pelty)
                 !
                 ! Gather operations
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
                    elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                    elpre(inode)         = press(ipoin,1)
                 end do
                 !
                 ! Properties
                 !
                 call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape)
                 call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,elmar(pelty)%shape)
                 call ker_proper('POROS','PGAUS',dummi,ielem,gppor,pnode,pgaus,elmar(pelty)%shape)

                 do igaus=1,pgaus 
                    call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&   ! Cartesian derivative
                         elcod,gpcar,detjm,xjacm,xjaci)                       ! and Jacobian
                    gpvol(igaus)=elmar(pelty)%weigp(igaus)*detjm
                    !
                    ! Velocity, pressure, velocity gradient and pressure gradient
                    !
                    prgau = 0.0_rp
                    gpvel = 0.0_rp
                    gradp = 0.0_rp
                    grave = 0.0_rp
                    gpcod = 0.0_rp
                    do inode = 1,pnode
                       do idime = 1,ndime
                          gpvel(idime) = gpvel(idime) + elmar(pelty)%shape(inode,igaus)*elvel(idime,inode)
                          gpcod(idime) = gpcod(idime) + elmar(pelty)%shape(inode,igaus)*elcod(idime,inode)
                          gradp(idime) = gradp(idime) + gpcar(idime,inode)*elpre(inode)
                          do jdime = 1,ndime
                             grave(jdime,idime) = grave(jdime,idime) + gpcar(jdime,inode)*elvel(idime,inode)
                          end do
                       end do
                       prgau = prgau + elmar(pelty)%shape(inode,igaus)*elpre(inode)
                    end do
                    !
                    ! Exact solution
                    ! 
                    call nsi_exacso(&
                         1_ip,gpcod,gpden(igaus),gpvis(igaus),gppor(igaus),&
                         gpgvi(1,igaus),exvel,exveg,expre,exprg,&
                         dummr,dummr,dummr,dummr)
                    !
                    ! Errors
                    !  
                    diffp    = abs(prgau-expre)
                    abpre    = abs(expre)
                    erp01(1) = erp01(1) + diffp*gpvol(igaus)
                    erp02(1) = erp02(1) + diffp*diffp*gpvol(igaus)
                    erp0i(1) = max(erp0i(1),diffp)
                    erp01(2) = erp01(2) + expre*gpvol(igaus)
                    erp02(2) = erp02(2) + expre*expre*gpvol(igaus)
                    erp0i(2) = max(erp0i(2),expre)
                    do idime=1,ndime
                       difeu    = abs(gpvel(idime)-exvel(idime))
                       abvel    = abs(exvel(idime))
                       eru01(1) = eru01(1) + difeu*gpvol(igaus)
                       eru02(1) = eru02(1) + difeu*difeu*gpvol(igaus)
                       eru0i(1) = max(eru0i(1),difeu)
                       eru01(2) = eru01(2) + abvel*gpvol(igaus)
                       eru02(2) = eru02(2) + abvel*abvel*gpvol(igaus)
                       eru0i(2) = max(eru0i(2),abvel)
                       diffp    = abs(gradp(idime)-exprg(idime))
                       abpre    = abs(exprg(idime)) 
                       erp11(1) = erp11(1) + diffp*gpvol(igaus)
                       erp12(1) = erp12(1) + diffp*diffp*gpvol(igaus)
                       erp1i(1) = max(erp1i(1),diffp)
                       erp11(2) = erp11(2) + abpre*gpvol(igaus)
                       erp12(2) = erp12(2) + abpre*abpre*gpvol(igaus)
                       erp1i(2) = max(erp1i(2),abpre)
                       do jdime=1,ndime
                          difeu    = abs(grave(idime,jdime)-exveg(idime,jdime))
                          abvel    = abs(exveg(idime,jdime))
                          eru11(1) = eru11(1) + difeu*gpvol(igaus)
                          eru12(1) = eru12(1) + difeu*difeu*gpvol(igaus)
                          eru1i(1) = max(eru1i(1),difeu)
                          eru11(2) = eru11(2) + abvel*gpvol(igaus)
                          eru12(2) = eru12(2) + abvel*abvel*gpvol(igaus)
                          eru1i(2) = max(eru1i(2),abvel)
                       end do
                    end do
                 end do
              end if
           end do elements

        call PAR_MAX(2_ip,eru0i)
        call PAR_MAX(2_ip,erp0i)
        call PAR_MAX(2_ip,eru1i)
        call PAR_MAX(2_ip,erp1i)

        call PAR_SUM(2_ip,erp02)
        call PAR_SUM(2_ip,erp12)
        call PAR_SUM(2_ip,eru02)
        call PAR_SUM(2_ip,eru12)

        call PAR_SUM(2_ip,erp01)
        call PAR_SUM(2_ip,erp11)
        call PAR_SUM(2_ip,eru01)
        call PAR_SUM(2_ip,eru11)

        erp02(1) = sqrt(erp02(1))
        erp12(1) = sqrt(erp12(1))
        eru02(1) = sqrt(eru02(1))
        eru12(1) = sqrt(eru12(1))

        xerro     = xerro + (eru02(1)+eru02_old)*dtime*0.5_rp
        eru02_old = eru02(1)
 
        if( erp01(2) > zeror ) erp01(1) = erp01(1)/erp01(2) 
        if( erp02(2) > zeror ) erp02(1) = erp02(1)/sqrt(erp02(2))
        if( erp0i(2) > zeror ) erp0i(1) = erp0i(1)/erp0i(2)
        if( eru01(2) > zeror ) eru01(1) = eru01(1)/eru01(2) 
        if( eru02(2) > zeror ) eru02(1) = eru02(1)/sqrt(eru02(2))
        if( eru0i(2) > zeror ) eru0i(1) = eru0i(1)/eru0i(2) 
        if( erp11(2) > zeror ) erp11(1) = erp11(1)/erp11(2) 
        if( erp12(2) > zeror ) erp12(1) = erp12(1)/sqrt(erp12(2))
        if( erp1i(2) > zeror ) erp1i(1) = erp1i(1)/erp1i(2)
        if( eru11(2) > zeror ) eru11(1) = eru11(1)/eru11(2) 
        if( eru12(2) > zeror ) eru12(1) = eru12(1)/sqrt(eru12(2))
        if( eru1i(2) > zeror ) eru1i(1) = eru1i(1)/eru1i(2) 

        if(INOTSLAVE) &
             write(momod(modul)%lun_outpu,100) &
             &       ittim,itcou,                                            &
             &       eru01(1),erp01(1),eru02(1),erp02(1),eru0i(1),erp0i(1),  &
             &       eru11(1),erp11(1),eru12(1),erp12(1),eru1i(1),erp1i(1)
        !if(INOTSLAVE) then
        !   write(90,*) cutim,eru02(1),xerro
        !   flush(90)
        !end if
        !
        ! Subgrid scale norm
        !
        if( kfl_sgsco_nsi /= 0 ) then
           sgsno=0.0_rp
           dummi=0
           if(INOTMASTER)then
              do ielem=1,nelem
                 pelty=ltype(ielem)
                 if(pelty>0)then
                    pnode=nnode(pelty)
                    pgaus=ngaus(pelty)
                    do igaus=1,pgaus
                       dummi=dummi+1
                       do idime=1,ndime
                          sgsno=sgsno+vesgs(ielem)%a(idime,igaus,1)**2
                       end do
                    end do
                 end if
              end do
           end if
           call PAR_SUM(sgsno)
           call PAR_SUM(dummi)
           sgsno = sqrt(sgsno) / real(dummi,rp)

           call norm2x(ndime,veloc,velno)
           dummi = npoin
           call PAR_SUM(dummi)
           velno = velno/real(dummi,rp)
           if( INOTSLAVE ) write(momod(modul)%lun_outpu,101) ittim,itcou,sgsno/velno 
           ! velno=0.0_rp
           ! do ipoin=1,npoin
           !    do idime=1,ndime
           !       velno=velno+veloc(idime,ipoin,1)**2
           !    end do
           ! end do
           ! velno=sqrt(velno)/real(npoin)
           ! sgsno=sqrt(sgsno)/real(dummi)
           ! write(momod(modul)%lun_outpu,101) ittim,itcou,sgsno/velno 
        end if

     end if
  end if
100 format(///,10X,'FINITE ELEMENT ERRORS',                              &
       &              /,10X,'=====================',//,                           &
       &  '          TIME STEP NO.',I5,',  ITERATION NO. ',I5,/,10X,44('-'),/,    &
       &  '          NORM       VELOCITY                 PRESSURE ',/,10X,44('-'),/,  &
       &  '          W(0,1) ',es16.8e3,5X,es16.8e3 ,/, &
       &  '          W(0,2) ',es16.8e3,5X,es16.8e3 ,/, &
       &  '          W(0,i) ',es16.8e3,5X,es16.8e3 ,/, &
       &  '          W(1,1) ',es16.8e3,5X,es16.8e3 ,/, &
       &  '          W(1,2) ',es16.8e3,5X,es16.8e3 ,/, &
       &  '          W(1,i) ',es16.8e3,5X,es16.8e3 ,/,10X,44('-'))
101 format(///,10X,'SUBGRID SCALE NORM',&
       &     /,10X,'==================',//,                           &
       &  '          TIME STEP NO.',I5,',  ITERATION NO. ',I5,/,10X,44('-'),/,    &
       &  '          NORM       VELOCITY ',/,10X,44('-'),/,  &
       &  '          W(0,2) ',es16.8e3)

end subroutine nsi_exaerr
