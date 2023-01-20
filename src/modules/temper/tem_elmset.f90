!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmset(iesec,ieset)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmset
  ! NAME 
  !    tem_elmset
  ! DESCRIPTION
  !    This routine computes variables on an element set W.
  !    The variable are:
  ! USES
  !    elmder
  !    tem_elmpro
  ! USED BY
  !    tem_outset
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper 
  implicit none

  integer(ip), intent(in)  :: iesec,ieset
  real(rp),    pointer     :: setvo(:)
  real(rp),    pointer     :: setmass(:)
  real(rp),    pointer     :: seth(:)
  real(rp),    pointer     :: setie(:)
  real(rp),    pointer     :: seths(:)

  integer(ip) :: ielem,inode,ipoin,igaus,idime,dummi                     
  integer(ip) :: pelty,pnode,pgaus,nvabi,pmate

  real(rp)    :: gpcar(ndime,mnode,mgaus) 
  real(rp)    :: xjaci(ndime,ndime), xjacm(ndime,ndime) 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elthe(mnode)
  real(rp)    :: eltem(mnode)
  real(rp)    :: gph(mgaus)
  real(rp)    :: gptem(mgaus)
  real(rp)    :: gpden(mgaus)
  real(rp)    :: gpsph(mgaus)
  real(rp)    :: gpcod(ndime,mgaus)
  real(rp)    :: gpvol,gpdet

  !
  ! Initialization
  !
  nvabi =  postp(1) % nvaes+1
  setvo => veset( nvabi:nvabi ,ieset)  
  setvo =  0.0_rp  ! Set volume

  if( postp(1) % npp_setse(1) /= 0 ) setmass => postp(1) % veset(1:1,ieset)
  if( postp(1) % npp_setse(2) /= 0 ) seth    => postp(1) % veset(2:2,ieset)
  if( postp(1) % npp_setse(3) /= 0 ) setie   => postp(1) % veset(3:3,ieset)
  if( postp(1) % npp_setse(4) /= 0 ) seths   => postp(1) % veset(4:4,ieset)

  if( postp(1) % npp_setse(1) /= 0 ) setmass = 0.0_rp 
  if( postp(1) % npp_setse(2) /= 0 ) seth    = 0.0_rp 
  if( postp(1) % npp_setse(3) /= 0 ) setie   = 0.0_rp 
  if( postp(1) % npp_setse(4) /= 0 ) seths   = 0.0_rp 

  !
  ! Loop over elements
  !
  elements: do ielem=1,nelem

     if(leset(ielem)==iesec) then
        ! 
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        pmate = lmate(ielem)
        
        !
        ! Gather operations
        !
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)

           elthe(inode) = therm(ipoin,1)
           eltem(inode) = tempe(ipoin,1)
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do

        !
        ! Properties
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
        call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)

        !
        ! Gauss point values
        !
        do igaus = 1,pgaus
           gph(igaus)     = 0.0_rp
           gptem(igaus)   = 0.0_rp
           gpcod(:,igaus) = 0.0_rp
        end do
        if (kfl_regim_tem == 4) then
           !
           ! Unknown is ENTHALPY
           !
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gph(igaus) = gph(igaus)                     &
                      +elmar(pelty)%shape(inode,igaus)       &
                      *elthe(inode)
                 gptem(igaus) = gptem(igaus)                 &
                      +elmar(pelty)%shape(inode,igaus)       &
                      *eltem(inode)
                 gpcod(1:ndime,igaus) = gpcod(1:ndime,igaus) &
                      +elmar(pelty)%shape(inode,igaus)       &
                      *elcod(1:ndime,inode)
              end do
           end do
        else
           !
           ! Unknown is TEMPERATURE
           !
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gptem(igaus) = gptem(igaus)             &
                      +elmar(pelty)%shape(inode,igaus)   &
                      *elthe(inode)
              end do
              gph(igaus) = gptem(igaus) * gpsph(igaus)
           end do
        endif

        gauss_points: do igaus=1,pgaus 
           !
           ! Cartesian derivatives and Jacobian
           !
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
           gpvol=elmar(pelty)%weigp(igaus)*gpdet
           setvo = setvo + gpvol
           !
           ! Set calculations
           !
           if( postp(1) % npp_setse(1) /= 0 ) then
              !
              ! Mass
              !
              setmass = setmass + gpvol * gpden(igaus)
           end if

           if( postp(1) % npp_setse(2) /= 0 ) then
              !
              ! Enthalpy in domain
              !
              seth = seth + gpvol * gph(igaus) * gpden(igaus)
           end if

           if( postp(1) % npp_setse(3) /= 0 ) then
              !
              ! Internal energy in domain
              !
              setie = setie + gpvol * (gph(igaus) * gpden(igaus) - prthe(1)) 
           end if

           if( postp(1) % npp_setse(4) /= 0 ) then
              !
              ! Source term
              !
              call tem_source(ielem,pmate,gpcod(:,igaus),gpsph(igaus))
              seths = seths + gpvol * gpsph(igaus)
           end if
           
        end do gauss_points

     end if

  end do elements

end subroutine tem_elmset
