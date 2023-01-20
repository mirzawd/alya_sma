!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_cratip(nctip)
  !------------------------------------------------------------------------
  !****f* Solidz/sld_cratip
  ! NAME
  !    sld_cratip
  ! DESCRIPTION
  !    Updating the position of crack tip during propagation
  !
  ! USES
  !
  ! USED BY
  !    sld_updcra
  !------------------------------------------------------------------------

  use def_master                ! general global variables
  use def_kermod
  use def_domain                ! geometry information
  use def_elmtyp                ! element types
  use def_solidz                ! general solidz module information
  use mod_cutele, only : plaseg
  use def_elmgeo, only : element_type
  implicit none
  integer(ip), intent(out) :: nctip
  integer(ip)              :: ifacg,inode,idime,ipoin,ielty,iface
  integer(ip)              :: inodb,pnodb,ifoun,jelem,pface 
  integer(ip)              :: jnodb,nfoun
  real(rp)                 :: xx(3),facoo(ndime,mnodb)
  real(rp)                 :: plapo(3)
  logical(lg)              :: debugging
  
  debugging = .false.

  !-------------------------------------------------------------------
  !
  ! Update cracked faces and positions
  !
  !-------------------------------------------------------------------
  nctip = 0_ip
     
  do jelem = 1,nelem
     if( leenr_sld(jelem) < 0 ) then
        leenr_sld(jelem) = abs(leenr_sld(jelem))
        nctip = nctip + 1 
        
        ielty = ltype(jelem)
        pface = element_type(ielty) % number_faces     
        xx    = 0.0_rp
        !
        ! Loop over faces
        !
        do iface=1,pface
           ifacg = lelfa(jelem) % l(iface)
           if( lcrkf_sld(ifacg) == 0 ) then

              pnodb = element_type(ielty) % node_faces(iface)
              do inodb = 1,pnodb
                 ! Local face node                
                 inode = element_type(ielty) % list_faces(inodb,iface)
                 ! Global face node (into domain)                
                 ipoin  = lnods(inode,jelem)
                 do idime = 1,ndime
                    ! Coordinate of face node
                    facoo(idime,inodb) = coord(idime,ipoin)                    
                 end do
              end do

              if (ndime == 2) then
                 nfoun = 0_ip
                 call plaseg(cranx_sld(:,jelem),crapx_sld(:,jelem),facoo(:,1),facoo(:,2),plapo,ifoun)
                 if (ifoun /= 0) then
                    nfoun = 2_ip
                    do idime = 1,ndime
                       xx(idime) = plapo(idime)
                    end do
                 end if
              else
                 !
                 ! Loop over edges (3D)
                 !          
                 nfoun = 0_ip
                 do idime = 1,ndime
                    xx(idime) = 0.0_rp
                 end do
                 do inodb = 1,pnodb
                    jnodb = inodb + 1
                    if (inodb == pnodb) jnodb = 1
                    call plaseg(cranx_sld(:,jelem),crapx_sld(:,jelem),facoo(:,inodb),facoo(:,jnodb),plapo,ifoun)
                    if (ifoun /= 0) then
                       nfoun = nfoun + 1_ip
                       do idime = 1,ndime
                          xx(idime) = xx(idime) + plapo(idime)
                       end do
                    end if
                 end do
                 if (nfoun > 1) then                          
                    do idime = 1,ndime
                       xx(idime) = xx(idime)/real(nfoun,rp)
                    end do
                 end if
!                 if (xx(2) > crapx_sld(2,jelem)) then  !*******
!                    write(*,*)'jelem = ',jelem
!                    write(*,*)'xx    = ',xx
!                    write(*,*)'crapx = ',crapx_sld(2,jelem)
!                 end if
              end if
              !
              ! The face intersect with the crack plane
              !
              if (nfoun > 1) then                          
                 ifacg = lelfa(jelem) % l(iface)
                 lcrkf_sld(ifacg) = lcrkf_sld(ifacg) + 1
                 do idime = 1,ndime
                    crtip_sld(idime,ifacg) = xx(idime)
!                    cockf_sld(idime,ifacg) = xx(idime)
                 end do
              end if
           end if
        end do

     end if
  end do

end subroutine sld_cratip
