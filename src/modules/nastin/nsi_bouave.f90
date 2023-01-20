!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouave
!-----------------------------------------------------------------------
!****f* Nastin/nsi_bouave
! NAME 
!    nsi_bouave
! DESCRIPTION
!    Compute average dynamic pressure on boundaries of type 5.
!   =-(-1/2*rho*u^2) n
! USES
!    nsi_elmdir
!    nsi_assrhs
!    nsi_assmat
! USED BY
!    nsi_bouope
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper
  use mod_bouder
  implicit none
  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: bovel(ndime,mnodb)                  ! Gather 
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: gpden(mgaus),gpvis(mgaus)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode)
  integer(ip) :: ielem,igaus,idime                   ! Indices and dimensions
  integer(ip) :: pnode,iboun,igaub,inodb
  integer(ip) :: pelty,pblty,pnodb,pgaub,pgaus
  integer(ip) :: inode,ipoin,dummi
  real(rp)    :: gbsur,eucta,gbden                   ! Values at Gauss points
  real(rp)    :: gbvel(ndime),gbven
  real(rp)    :: dynpr,sutot,tmatr
  real(rp)    :: alpha,tragl(9),hleng(3)

  if( IMASTER ) return
  !
  ! Loop over boundaries
  !
  dynpr = 0.0_rp
  sutot = 0.0_rp
  alpha = 0.5_rp

  boundaries: do iboun=1,nboun

     if(  kfl_fixbo_nsi(iboun)==5.or.&
          kfl_fixbo_nsi(iboun)==10) then 
        !
        ! Element properties and dimensions
        !
        pblty=ltypb(iboun) 
        pnodb=nnode(pblty)
        ielem=lelbo(iboun)
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        pgaub=ngaus(pblty)
        pgaus=ngaus(pelty)
        !
        ! Gather operations
        !
        call nsi_bougat(ndime,pnodb,npoin,lnodb(1,iboun),&
             bovel,bocod,veloc,coord)
        if(ndime==2) then
           do inode=1,pnode 
              ipoin=lnods(inode,ielem)
              elvel(1,inode) = veloc(1,ipoin,1)
              elvel(2,inode) = veloc(2,ipoin,1)
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
           end do
        else
           do inode=1,pnode 
              ipoin=lnods(inode,ielem)
              elvel(1,inode) = veloc(1,ipoin,1)
              elvel(2,inode) = veloc(2,ipoin,1)
              elvel(3,inode) = veloc(3,ipoin,1)
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
              elcod(3,inode) = coord(3,ipoin)
           end do
        end if
        !
        ! Element length HLENG
        !
        call elmlen(&
             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Properties: mu and rho  
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden)
        call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis)

        gauss_points: do igaub=1,pgaub
           
           call bouder(&
                pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                bocod,baloc,eucta)                                   ! and Jacobian
           gbsur=elmar(pblty)%weigp(igaub)*eucta 
           !
           ! Properties: rho  
           !
           gbden=0.0_rp
           do inodb=1,pnodb 
              do igaus=1,ngaus(pelty)
                 tmatr=elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                      *elmar(pblty)%shape(inodb,igaub)
                 gbden=gbden+gpden(igaus)*tmatr
              end do
           end do 

           do idime=1,ndime
              gbvel(idime)=0.0_rp
           end do
           do inodb=1,pnodb
              do idime=1,ndime
                 gbvel(idime)=gbvel(idime)+elmar(pblty)%shape(inodb,igaub)&
                      *bovel(idime,inodb)
              end do
           end do
           call vecnor(gbvel,ndime,gbven,2_ip)
           !if(1==1) then
           !   gbcod=0.0_rp
           !   do inodb=1,pnodb
           !      do idime=1,ndime
           !         gbcod(idime)=gbcod(idime)&
           !              +elmar(pblty)%shape(inodb,igaub)&
           !              *bocod(idime,inodb)
           !      end do
           !   end do
           !   alpha=0.5_rp+1.2e-5/(gbcod(1)**1.6_rp)
           !end if
           !alpha=0.0_rp
           dynpr=dynpr-alpha*gbden*gbven*gbven*gbsur
           sutot=sutot+gbsur

        end do gauss_points
        
    end if
     
  end do boundaries
  !
  ! Update natural boundary condition
  !
  if(sutot>zensi) then 
     dynpr=dynpr/sutot
     do iboun=1,nboun
        if(  kfl_fixbo_nsi(iboun)==5.or.&
             kfl_fixbo_nsi(iboun)==10) then
           bvnat_nsi(1,iboun,1)=relbc_nsi*dynpr+(1.0_rp-relbc_nsi)*bvnat_nsi(1,iboun,1)
        end if
     end do
  end if

end subroutine nsi_bouave
