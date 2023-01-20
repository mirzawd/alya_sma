!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmset(iesec,ieset)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmset
  ! NAME 
  !    nsi_elmset
  ! DESCRIPTION
  !    This routine computes variables on an element set W.
  !    The variable are: 
  !    1. SETVO: set surface           = meas(W)=int_W
  !    2. SETVE: set mean vel. module  = int_W u^2 ]
  !    3. SETVR: set mean vort. module = int_W w^2 ]
  !                                      where w=dv/dx-du/dy
  !    SETVE and SETVR are normailzed further on in nsi_outset
  ! USES
  !    nsi_elmgat
  !    elmder
  ! USED BY
  !    nsi_outset
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  use mod_ker_proper, only : ker_proper
  use mod_arrays,     only : arrays_number
  implicit none
  integer(ip), intent(in)  :: iesec,ieset
  real(rp),    pointer     :: setvo(:),setve(:),setvr(:), setdi(:)
  real(rp),    pointer     :: setki(:),divu2(:),setmo(:), setpr(:)
  integer(ip)              :: pnode,pgaus,pelty,nvabi,dummi
  integer(ip)              :: ielem,igaus,idime,jdime,inode,ipoin, i
  real(rp)                 :: gpcar(ndime,mnode,mgaus) 
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime) 
  real(rp)                 :: elgra(ndime,ndime), strain(ndime,ndime), strain2(ndime,ndime)
  real(rp)                 :: elvel(ndime,mnode),elcod(ndime,mnode), elpre(mnode)
  real(rp)                 :: elm_vavep(ndime,mnode), elm_vminp(ndime,mnode), elm_vmaxp(ndime,mnode), elm_pinde(mnode)
  real(rp)                 :: gp_vavep(3), gp_vminp(3), gp_vmaxp(3), gp_pinde
  real(rp), pointer        :: setvma(:), setvmi(:), setvav(:), setpin(:)
  real(rp)                 :: gpden(mgaus),divu
  real(rp)                 :: dvdx,dudy, gpvol,gpdet,gpvel(3), gppre ,dummr
  logical(lg)              :: pindex
  !
  ! Initialization
  !
  elgra=0.0_rp
  nvabi =  postp(1) % nvaes+1
  setvo => veset( nvabi:nvabi ,ieset)  
  setvo =  0.0_rp  ! Set volume
  
  if( postp(1) % npp_setse(1) /= 0 )  setve  => postp(1) % veset(1:1,ieset)
  if( postp(1) % npp_setse(2) /= 0 )  setvr  => postp(1) % veset(2:2,ieset)
  if( postp(1) % npp_setse(3) /= 0 )  setki  => postp(1) % veset(3:3,ieset)
  if( postp(1) % npp_setse(4) /= 0 )  divu2  => postp(1) % veset(4:4,ieset)
  if( postp(1) % npp_setse(5) /= 0 )  setmo  => postp(1) % veset(5:7,ieset)
  if( postp(1) % npp_setse(8) /= 0 )  setpr  => postp(1) % veset(8:8,ieset)
  if( postp(1) % npp_setse(9) /= 0 )  setdi  => postp(1) % veset(9:9,ieset)
  if( postp(1) % npp_setse(10) /= 0 ) setvma => postp(1) % veset(10:10,ieset)
  if( postp(1) % npp_setse(11) /= 0 ) setvmi => postp(1) % veset(11:11,ieset)
  if( postp(1) % npp_setse(12) /= 0 ) setvav => postp(1) % veset(12:12,ieset)
  if( postp(1) % npp_setse(13) /= 0 ) setpin => postp(1) % veset(13:13,ieset)

  if( postp(1) % npp_setse(1) /= 0 ) setve = 0.0_rp  ! Set mean velocity module
  if( postp(1) % npp_setse(2) /= 0 ) setvr = 0.0_rp  ! Set mean vorticity module
  if( postp(1) % npp_setse(3) /= 0 ) setki = 0.0_rp  ! Set kinetic turbulent energy
  if( postp(1) % npp_setse(4) /= 0 ) divu2 = 0.0_rp  ! rho * (div u) * u^2
  if( postp(1) % npp_setse(5) /= 0 ) setmo = 0.0_rp  ! rho * u 
  if( postp(1) % npp_setse(8) /= 0 ) setpr = 0.0_rp  ! Set mean pressure
  if( postp(1) % npp_setse(9) /= 0 ) setdi = 0.0_rp  ! Set distortion time
  if( postp(1) % npp_setse(10) /= 0 ) setvma = 0.0_rp  ! Set veloc max
  if( postp(1) % npp_setse(11) /= 0 ) setvmi = 0.0_rp  ! Set veloc min
  if( postp(1) % npp_setse(12) /= 0 ) setvav = 0.0_rp  ! Set veloc ave
  if( postp(1) % npp_setse(13) /= 0 ) setpin = 0.0_rp  ! Set veloc pulsatility index

  pindex=.False.
  if( postp(1) % npp_setse(10) /= 0 .or. &
      postp(1) % npp_setse(11) /= 0 .or. &
      postp(1) % npp_setse(12) /= 0 .or. &
      postp(1) % npp_setse(13) /= 0 ) then
        pindex=.True.
  endif
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem

     if( leset(ielem) == iesec ) then
        ! 
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! Gather operations
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elvel(1:ndime,inode)      = veloc(1:ndime,ipoin,1)
              elcod(1:ndime,inode)      = coord(1:ndime,ipoin)
              elpre(inode)              = press(ipoin,1)
                
              if(pindex) then
                elm_vmaxp(1:ndime,inode)  = vmaxp_nsi(1:ndime,ipoin)
                elm_vminp(1:ndime,inode)  = vminp_nsi(1:ndime,ipoin)
                elm_vavep(1:ndime,inode)  = vavep_nsi(1:ndime,ipoin)
                elm_pinde(inode)          = pinde_nsi(1,ipoin)
              endif

           end do
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           !
           ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
           !
           do igaus = 1,pgaus     
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
              gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg
              setvo = setvo + gpvol
              gpvel     = 0.0_rp
              gppre     = 0.0_rp
              gp_vmaxp  = 0.0_rp
              gp_vminp  = 0.0_rp
              gp_vavep  = 0.0_rp
              gp_pinde  = 0.0_rp
              do inode = 1,pnode
                 gpvel(1:ndime) = gpvel(1:ndime)          &
                      + elmar(pelty) % shape(inode,igaus) &
                      * elvel(1:ndime,inode)
                 gppre = gppre  &
                      + elmar(pelty) % shape(inode,igaus) &
                      * elpre(inode)
                 if(pindex) then
                      gp_vmaxp(1:ndime) = gp_vmaxp(1:ndime)   &
                          + elmar(pelty) % shape(inode,igaus) &
                          * elm_vmaxp(1:ndime,inode)
                      gp_vminp(1:ndime) = gp_vminp(1:ndime)   &
                          + elmar(pelty) % shape(inode,igaus) &
                          * elm_vminp(1:ndime,inode)
                      gp_vavep(1:ndime) = gp_vavep(1:ndime)   &
                          + elmar(pelty) % shape(inode,igaus) &
                          * elm_vavep(1:ndime,inode)
                      gp_pinde = gp_pinde   &
                          + elmar(pelty) % shape(inode,igaus) &
                          * elm_pinde(inode)
                 endif
              end do
              if( postp(1) % npp_setse(1) /= 0 ) then
                 !
                 ! Velocity average module: SETVE
                 !
                 do idime=1,ndime
                    setve=setve+gpvol*gpvel(idime)**2
                 end do
              end if
              if( postp(1) % npp_setse(2) /= 0 ) then
                 !
                 ! Vorticity average module: SETVR
                 !
                 dvdx = 0.0_rp
                 dudy = 0.0_rp
                 do inode = 1,pnode
                    dvdx = dvdx + gpcar(1,inode,igaus) * elvel(2,inode)
                    dudy = dudy + gpcar(2,inode,igaus) * elvel(1,inode)
                 end do
                 setvr = setvr + gpvol*(dvdx-dudy)**2
              end if
              if( postp(1) % npp_setse(3) /= 0 ) then
                 !
                 ! Kinetic energy: SETKI
                 !
                 dummr = 0.0_rp
                 do idime = 1,ndime
                    dummr = dummr + gpvel(idime)**2
                 end do
                 setki = setki + 0.5_rp * gpden(igaus) * gpvol * dummr
              end if
              if( postp(1) % npp_setse(4) /= 0 ) then
                 !
                 ! rho * (div u) * u^2
                 !
                 dummr = 0.0_rp
                 divu  = 0.0_rp
                 do inode = 1,pnode
                    divu = divu + dot_product(gpcar(:,inode,igaus),elvel(:,inode))
                 end do
                 do idime = 1,ndime
                    dummr = dummr + gpvel(idime)**2
                 end do
                 divu2 = divu2 + gpden(igaus) * gpvol * dummr * divu
              end if
              if( postp(1) % npp_setse(5) /= 0 ) then
                 !
                 ! Momentum
                 !
                 do idime = 1,ndime
                    setmo(idime) = setmo(idime)                        &
                                  + gpden(igaus) * gpvel(idime) * gpvol
                 end do
              end if
              if( postp(1) % npp_setse(8) /= 0 ) then
                 !
                 ! Pressure
                 !
                 setpr=setpr+gpvol*gppre
              endif
              if( postp(1) % npp_setse(9) /= 0 ) then
                 !
                 ! Distortion time average: SETDI
                 !
                 dummr = 0.0_rp
                 elgra=0.0_rp
                 strain=0.0_rp
                 strain2=0.0_rp
                 do inode = 1,pnode
                    do idime=1,ndime
                        do jdime=1,ndime
                            elgra(idime,jdime)=elgra(idime,jdime) + gpcar(idime,inode,igaus)*elvel(jdime,inode)
                        enddo
                    enddo
                 end do

                 do idime = 1,ndime
                    do jdime = 1,ndime
                        strain(idime,jdime) = 0.5_rp * (elgra(idime,jdime)+elgra(jdime,idime))
                    end do
                 end do

                 strain2=matmul(strain,strain)

                 do i=1,ndime
                    dummr=dummr+strain2(i,i)
                 enddo

                 dummr = sqrt(dummr/2.0_rp)


                 setdi = setdi + gpvol/dummr

              endif
              if( postp(1) % npp_setse(10) /= 0 .and. ittim>0  ) then
                 !
                 ! Maximum velocity: VMAXP
                 !
                 if( mod(ittim, postp(1) % npp_stepi(arrays_number('VMAXP'),0) ) .eq. 0_ip) then
                   if(ndime.eq.2_ip) then
                      setvma=setvma+gpvol*sqrt(gp_vmaxp(1)**2+gp_vmaxp(2)**2)
                   else
                      setvma=setvma+gpvol*sqrt(gp_vmaxp(1)**2+gp_vmaxp(2)**2+gp_vmaxp(3)**2)
                   endif
                 endif

              endif
              if( postp(1) % npp_setse(11) /= 0 .and. ittim>0 ) then
                 !
                 ! Minimum velocity: VMINP
                 !
                 if( mod(ittim, postp(1) % npp_stepi(arrays_number('VMINP'),0) ) .eq. 0_ip ) then
                   if(ndime.eq.2_ip) then
                      setvmi=setvmi+gpvol*sqrt(gp_vminp(1)**2+gp_vminp(2)**2)
                   else
                      setvmi=setvmi+gpvol*sqrt(gp_vminp(1)**2+gp_vminp(2)**2+gp_vminp(3)**2)
                   endif
                 endif

              endif
              if( postp(1) % npp_setse(12) /= 0 .and. ittim>0 ) then
                 !
                 ! Average velocity: VAVEP
                 !
                 if( mod(ittim, postp(1) % npp_stepi(arrays_number('VAVEP'),0) ) .eq. 0_ip ) then
                   if(ndime.eq.2_ip) then
                      setvav=setvav+gpvol*sqrt(gp_vavep(1)**2+gp_vavep(2)**2)
                   else
                      setvav=setvav+gpvol*sqrt(gp_vavep(1)**2+gp_vavep(2)**2+gp_vavep(3)**2)
                   endif
                 endif

              endif
              if( postp(1) % npp_setse(13) /= 0  .and. ittim>0 ) then
                 !
                 ! Pulsatility index: PINDE
                 !
                 if( mod(ittim, postp(1) % npp_stepi(arrays_number('PINDE'),0) ) .eq. 0_ip ) then
                    setpin=setpin+gpvol*gp_pinde
                 endif

              endif
           end do
        end if
     end if

  end do elements

end subroutine nsi_elmset
