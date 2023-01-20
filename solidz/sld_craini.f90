!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_craini(itask)
  !------------------------------------------------------------------------
  !****f* Solidz/sld_craini
  ! NAME
  !    sld_craini
  ! DESCRIPTION
  !    Interpreting the initial crack geometry
  !
  !
  ! USES
  !
  ! USED BY
  !    sld_iniunk
  !------------------------------------------------------------------------

  use def_master   ! general global variables
  use def_domain   ! geometry information
  use def_solidz   ! general solidz module information
  use mod_cutele
  use mod_maths,  only                   : maths_normalize_vector
  use def_elmgeo, only                   : element_type
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem,idime,inode
  integer(ip)             :: ipoin,icrak,ifoun,nfoun,nfout
  integer(ip)             :: pface,nfofa,iface,pnodb,ninit
  integer(ip)             :: ifacg,ielty,inodb,jnodb,condi
  real(rp)                :: dista,plbox(ndime,2),fabox(ndime,2)
  real(rp)                :: facoo(ndime,mnodb),plapo(ndime),avpoi(ndime)
  real(rp)                :: crnor(ndime),crpoi(ndime),avpot(ndime)

  if (itask == 1) then
     !
     ! Identify cracked elements
     !
     do icrak = 1,ncrak_sld

        !call elsest(&
        !     2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
        !     lnods,ltype,ltopo,coord,crkpo_sld(1,icrak),relse,&
        !     ielem,shap1,deri1,coloc,dummi)
        !if( ielem < 1 .and. ISEQUEN ) then
        !   call runend('SLD_XFEMIC: COULD NOT FIND CRACKED ELEMENT')
        !else if( ielem >= 1 ) then
        !
        ! Determine the normal of the initial crack 
        !
        ninit = 0_ip
        if (ndime == 2) then
           crnor(1) = -(crkco_sld(2,2,icrak) - crkco_sld(2,1,icrak)) 
           crnor(2) =   crkco_sld(1,2,icrak) - crkco_sld(1,1,icrak)

           crpoi(1) = crkco_sld(1,1,icrak)
           crpoi(2) = crkco_sld(2,1,icrak)
           ninit    = 2 

        elseif (ndime == 3) then

           ! norma(x) = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2) 
           crnor(1) = crkco_sld(2,1,icrak)*(crkco_sld(3,2,icrak) - crkco_sld(3,3,icrak)) + &
                crkco_sld(2,2,icrak)*(crkco_sld(3,3,icrak) - crkco_sld(3,1,icrak)) + &
                crkco_sld(2,3,icrak)*(crkco_sld(3,1,icrak) - crkco_sld(3,2,icrak))
           ! norma(y) = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2) 
           crnor(2) = crkco_sld(3,1,icrak)*(crkco_sld(1,2,icrak) - crkco_sld(1,3,icrak)) + &
                crkco_sld(3,2,icrak)*(crkco_sld(1,3,icrak) - crkco_sld(1,1,icrak)) + &
                crkco_sld(3,3,icrak)*(crkco_sld(1,1,icrak) - crkco_sld(1,2,icrak))
           ! norma(z) = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2) 
           crnor(3) = crkco_sld(1,1,icrak)*(crkco_sld(2,2,icrak) - crkco_sld(2,3,icrak)) + &
                crkco_sld(1,2,icrak)*(crkco_sld(2,3,icrak) - crkco_sld(2,1,icrak)) + &
                crkco_sld(1,3,icrak)*(crkco_sld(2,1,icrak) - crkco_sld(2,2,icrak))                 
           crpoi(1) = crkco_sld(1,1,icrak)
           crpoi(2) = crkco_sld(2,1,icrak)
           crpoi(3) = crkco_sld(3,1,icrak)
           !
           ! Check if the last poin is inside the plane
           !
           dista = 0.0_rp
           do idime = 1,ndime
              dista = dista  + crnor(idime)  * ( crkco_sld(idime,4,icrak) - crpoi(idime) )
           end do
           if( dista > 1.0e-8_rp ) then   !! denny.tjahjanto >> this value of tolerance ?
              call runend('SLD_XFEMIC: THE 4TH POINT IS NOT COPLANAR')           
           end if
           ninit    = 4 
        end if
        call maths_normalize_vector(ndime,crnor)
        !
        ! Boundary box of the initial crack 
        !                 
        do idime = 1,ndime
           plbox(idime,1) =  1.0e12_rp
           plbox(idime,2) = -1.0e12_rp
           do ipoin = 1,ninit
              plbox(idime,1) = min( plbox(idime,1) , crkco_sld(idime,ipoin,icrak) )
              plbox(idime,2) = max( plbox(idime,2) , crkco_sld(idime,ipoin,icrak) )
           end do
           dista = plbox(idime,2) - plbox(idime,1)
           if (dista < 2.0e-8_rp) then   !! denny.tjahjanto >> this value of tolerance ?
              plbox(idime,1) = plbox(idime,1) - 2.0e-8_rp
              plbox(idime,2) = plbox(idime,2) + 2.0e-8_rp
           else
              plbox(idime,1) = plbox(idime,1) - abs(dista*2.0e-8_rp)
              plbox(idime,2) = plbox(idime,2) + abs(dista*2.0e-8_rp)
           end if
        end do
        !
        ! Loop over elements
        !
        do ielem = 1,nelem
           ielty = ltype(ielem)
           pface = element_type(ielty) % number_faces
           nfofa = 0_ip                       
           !
           ! Loop over faces
           !
           nfout    = 0
           do idime = 1,ndime
              avpot(idime) = 0.0_rp
           end do
           do iface = 1,pface
              pnodb = nnode(element_type(ielty) % type_faces(iface))    
              do idime = 1,ndime
                 fabox(idime,1) =  1.0e12_rp
                 fabox(idime,2) = -1.0e12_rp
              end do
              do inodb = 1,pnodb
                 ! Local face node                
                 inode = element_type(ielty) % list_faces(inodb,iface)
                 ! Global face node (into domain)                
                 ipoin  = lnods(inode,ielem)
                 do idime = 1,ndime
                    facoo(idime,inodb) = coord(idime,ipoin)
                    fabox(idime,1)     = min( fabox(idime,1) , coord(idime,ipoin) )
                    fabox(idime,2)     = max( fabox(idime,2) , coord(idime,ipoin) )                    

                    dista = fabox(idime,2) - fabox(idime,1)
                    if (dista < 2.0e-8_rp) then    !! denny.tjahjanto >> this value of tolerance ?
                       fabox(idime,1) = fabox(idime,1) - 2.0e-8_rp
                       fabox(idime,2) = fabox(idime,2) + 2.0e-8_rp
                    else
                       fabox(idime,1) = fabox(idime,1) - abs(dista*2.0e-8_rp)
                       fabox(idime,2) = fabox(idime,2) + abs(dista*2.0e-8_rp)
                    end if
                 end do
              end do
              !
              ! Check if the BB of the current face intersect with the baoundary box of the initial crack 
              !
              condi = 0
              do idime = 1,ndime
                 if ( fabox(idime,1) < plbox(idime,2) .and. fabox(idime,2) > plbox(idime,1) ) then
                    condi = condi + 1_ip
                 end if
              end do
              if( condi == ndime ) then
                 if (ndime == 2) then
                    nfoun = 0
                    call plaseg(crnor(:),crpoi(:),facoo(:,1),facoo(:,2),plapo,ifoun)
                    if (ifoun /= 0) then
                       nfout = nfout + 1
                       nfoun = 2
                       do idime = 1,ndime
                          avpoi(idime) = plapo(idime)
                          avpot(idime) = avpot(idime) + plapo(idime)
                       end do
                    end if
                 else
                    !
                    ! Loop over edges (3D)
                    !          
                    nfoun = 0_ip
                    do idime = 1,ndime
                       avpoi(idime) = 0.0_rp
                    end do
                    do inodb = 1,pnodb
                       jnodb = inodb + 1
                       if (inodb == pnodb) jnodb = 1
                       call plaseg(crnor(:),crpoi(:),facoo(:,inodb),facoo(:,jnodb),plapo,ifoun)
                       if (ifoun /= 0) then
                          nfout = nfout + 1
                          nfoun = nfoun + 1_ip
                          do idime = 1,ndime
                             avpoi(idime) = avpoi(idime) + plapo(idime)
                             avpot(idime) = avpot(idime) + plapo(idime)
                          end do
                       end if
                    end do
                    if (nfoun > 0) then                          
                       do idime = 1,ndime
                          avpoi(idime) = avpoi(idime)/real(nfoun,rp)
                       end do
                    end if
                 end if
                 !
                 ! The face intersect with the crack plane
                 !
                 if (nfoun > 1) then                          
                    nfofa = nfofa + 1
                    ifacg = lelfa(ielem) % l(iface)
                    lcrkf_sld(ifacg) = lcrkf_sld(ifacg) + 1
                    do idime = 1,ndime
                       cockf_sld(idime,ifacg) = avpoi(idime)
                    end do
                 end if
              end if
           end do
           !
           ! The element intersect with the crack plane
           !
           if (nfofa > 2) then
              do idime = 1,ndime
                 cranx_sld(idime,ielem) = crnor(idime) 
                 crapx_sld(idime,ielem) = avpot(idime)/real(nfout,rp)
              end do
              leenr_sld(ielem) = 1
              lelch(ielem)     = ELCUT
           end if
        end do
     end do
     !
     ! Parall
     !
     call parari('SLX',NFACE_TYPE,nfacg,lcrkf_sld)
     call pararr('SLX',NFACE_TYPE,ndime*nfacg,cockf_sld)
     do ifacg = 1,nfacg
        lcrkf_sld(ifacg) = min(1_ip,lcrkf_sld(ifacg))
     end do

  end if

end subroutine sld_craini
