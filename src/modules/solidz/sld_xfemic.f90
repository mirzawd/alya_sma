!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_xfemic(itask,pnode,pgaus,gphea,gpcod,elcod,ielem,eltip)
  !------------------------------------------------------------------------
  !****f* Solidz/sld_xfemic
  ! NAME
  !    sld_xfemic
  ! DESCRIPTION
  !    Standard operations for enrichment and crack propagation in XFEM
  !
  ! Check the XFEMicity of the element and act consequently
  !
  ! USES
  ! USED BY
  !    sld_elmope
  !------------------------------------------------------------------------

  use def_master                ! general global variables
  use def_kermod
  use def_domain                ! geometry information
  use def_elmtyp                ! element types
  use def_solidz                ! general solidz module information
  use mod_cutele, only : plaseg
  use def_elmgeo, only : element_type
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,ielem,itask
  integer(ip), intent(out) :: eltip
  real(rp),    intent(out) :: gphea(pnode,mgaus)
  real(rp),    intent(in)  :: gpcod(ndime,mgaus)
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  integer(ip)              :: kcrac,ifacg,inode,idime,ipoin,ielty,iface
  integer(ip)              :: igaus
  real(rp)                 :: gaupx(ndime),digau,gaumo,nodmo,dinod
  real(rp)                 :: xx(3),nodpx(ndime)
  logical(lg)              :: debugging
  
  debugging = .false.

  if (itask == 1) then
     !
     ! Finding crack tip based on existing crack path
     !
     if( leenr_sld(ielem) == 0 ) then
        !
        ! Indentify existing crack tip 
        !
        do idime = 1,ndime
           xx(idime) = 0.0_rp
        end do
        kcrac = 0
        ielty = ltype(ielem)
        do iface = 1,element_type(ielty) % number_faces  
           ifacg = lelfa(ielem) % l(iface)
!           write(*,*)'ifacg = ',ifacg,'lcrkf = ',lcrkf_sld(ifacg)
           if( lcrkf_sld(ifacg) /= 0 ) then 
              kcrac = kcrac + 1
              do idime = 1,ndime
                 xx(idime) = xx(idime) + cockf_sld(idime,ifacg)
              end do
!              if (kcrac > 0) then
!                 write(*,*)'kcrac = ',kcrac
!                 write(*,*)'ielem =',ielem,'iface =',iface,'ifacg =',ifacg
!                 write(*,*)'ielem =',ielem,'cockf =',cockf_sld(:,ifacg)
!              end if
           end if
        end do

        if( kcrac > 0 ) then
           !
           ! Compute crack position: CRAPX from previous crack tip
           !
           do idime = 1,ndime
              crapx_sld(idime,ielem) = xx(idime) / real(kcrac,rp)
           end do
           eltip = 1_rp
           if (debugging) write(*,*) 'FROM PREVIOUS'
        else
           !
           ! Compute crack position: CRAPX defined at the center of element
           !
           do idime = 1,ndime
              crapx_sld(idime,ielem) = 0.0_rp
              do inode = 1,pnode
                 crapx_sld(idime,ielem) = crapx_sld(idime,ielem) + elcod(idime,inode)
              end do
              crapx_sld(idime,ielem) = crapx_sld(idime,ielem) / real(pnode,rp)
           end do
           eltip = 0_rp
           if (debugging) write(*,*) 'ELEMENT CENTER'
        end if
        
     end if

  else if (itask == 2) then
     !
     ! Compute the value for the shifted-heaviside function:
     ! based on position of gauss point and position of node (with respect to crack)
     !
     gphea = 0.0_rp

     if( leenr_sld(ielem) /= 0 ) then
        !
        ! loop over the gauss points to see whether the gauss point is on one side or the other
        ! this should be done based either on nodal values (through level sets previously computed) or
        ! on element values (through point and normal, crapx and cranx)
        !
        if( kfl_xfcra_sld == 1 ) then
           !
           ! A: through point and normal, crapx and cranx
           !
           do igaus = 1,pgaus
              !
              ! 1. compute position vector of gauss point: gaupx= gpcod - crapx
              !
              gaumo = 0.0_rp
              do idime = 1,ndime
                 gaupx(idime) = gpcod(idime,igaus) - crapx_sld(idime,ielem)
                 gaumo        = gaumo + gaupx(idime) * gaupx(idime)
              end do
              gaumo = sqrt(gaumo)
              !
              ! 2. project gaupx onto cranx to see on which side of the crack the gauss point is
              !
              digau = 0.0_rp
              do idime = 1,ndime
                 digau = digau + cranx_sld(idime,ielem) * gaupx(idime)/gaumo
              end do

              do inode = 1,pnode
                 !
                 ! a. compute position vector of node: nodpx= coord - crapx
                 !
                 ipoin = lnods(inode,ielem)
                 nodmo = 0.0_rp
                 do idime = 1,ndime
                    nodpx(idime) = coord(idime,ipoin) - crapx_sld(idime,ielem)
                    nodmo        = nodmo + nodpx(idime) * nodpx(idime)
                 end do
                 nodmo = sqrt(nodmo)
                 !
                 ! b. project nodpx onto cranx to see on which side of the crack the node is
                 !
                 dinod = 0.0_rp
                 do idime = 1,ndime
                    dinod = dinod + cranx_sld(idime,ielem) * nodpx(idime) / nodmo
                 end do
                 !
                 ! c. compute the shifted heaviside depending on the position of gauss point and
                 !    the node w.r.t. the cranx
                 !
                 if( digau < 0.0_rp .and. dinod < 0.0_rp ) then
                    gphea(inode,igaus) =  0.0_rp
                 else if( digau > 0.0_rp .and. dinod < 0.0_rp ) then
                    gphea(inode,igaus) =  1.0_rp
                 else if( digau < 0.0_rp .and. dinod > 0.0_rp ) then
                    gphea(inode,igaus) = -1.0_rp
                 else if( digau > 0.0_rp .and. dinod > 0.0_rp ) then
                    gphea(inode,igaus) =  0.0_rp
                 end if
              end do
           end do

        else if( kfl_xfcra_sld == 2 ) then
           !
           ! B: through level sets previously computed
           !
           ! TO BE DONE...

        end if
     end if

  end if


end subroutine sld_xfemic
