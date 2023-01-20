!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_embedd(Auu,Aup,Apu,App,Q,bu,bp,u,p)
  !-----------------------------------------------------------------------
  !****f* mathru/assma3
  ! NAME 
  !    assma3
  ! DESCRIPTION
  !    Assembly an elemental matrix ELMAT in global matrix AMATR
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_domain
  use def_master
  use def_nastin

  use mod_gradie
  use mod_kdtree

  implicit none
  real(rp),    intent(out)   :: Auu(ndime,ndime,nzdom) 
  real(rp),    intent(out)   :: Aup(ndime,nzdom)
  real(rp),    intent(out)   :: Apu(ndime,nzdom)
  real(rp),    intent(out)   :: App(nzdom)
  real(rp),    intent(out)   :: Q(nzdom)
  real(rp),    intent(out)   :: bu(ndime,npoin)
  real(rp),    intent(out)   :: bp(npoin)
  real(rp),    intent(out)   :: u(*)
  real(rp),    intent(out)   :: p(*)

  integer(ip)                :: ipoin,jpoin,kpoin
  integer(ip)                :: idime,jdime,izdom,idofn,iimbo
  integer(ip)                :: izdod,iinde,limit,dummi
  integer(ip), pointer       :: lnode(:)

  real(rp)                   :: Auud(3),Qd,Appd,coor1(ndime)
  real(rp)                   :: x(3),v(3),dumma(ndime),dummr,propo(ndime)
  real(rp),    pointer       :: shapl(:)

  if( solve(2)%kfl_symme == 1 ) call runend('SYMMETRIC ASSEMBLY NOT POSSIBLE')
  !----------------------------------------------------------------------
  !
  ! Impose velocity of fringe nodes in interpolating way 
  !
  !----------------------------------------------------------------------
  do ipoin = 1,npoin
     if( lntib(ipoin) > 0 ) then
        !
        ! IZDOD: Diagonal
        !
        izdod = r_dom(ipoin) - 1
        jpoin = 0
        do while( jpoin /= ipoin )
           izdod = izdod + 1
           jpoin = c_dom(izdod)
        end do

        do idime = 1,ndime
           Auud(idime) = 1.0_rp
        end do
        Appd = 1.0_rp
        Qd   = 1.0_rp
        !
        ! Set line to zero
        !
        do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
           do idime = 1,ndime
              do jdime = 1,ndime
                 Auu(jdime,idime,izdom) = 0.0_rp
              end do
              Aup(idime,izdom) = 0.0_rp
              Apu(idime,izdom) = 0.0_rp
           end do
           App(izdom) = 0.0_rp
           Q(izdom)   = 0.0_rp
        end do
        !
        ! Presrcibe value to zero
        !
        do idime = 1,ndime
           Auu(idime,idime,izdod) = Auud(idime)
        end do
        App(izdod) = Appd
        Q(izdod)   = Qd

        bp(ipoin)          = 0.0_rp
        p(ipoin)           = 0.0_rp
        idofn              = (ipoin-1)*ndime
        do idime = 1,ndime
           idofn           = idofn + 1
           u(idofn)        = 0.0_rp
           bu(idime,ipoin) = 0.0_rp
        end do
     elseif( lntib(ipoin) < 0 ) then

        iimbo = abs(lntib(ipoin))
        !
        ! Node information use to interpolate the particle velocity
        !
        lnode => lnint(ipoin) % lnode
        shapl => lnint(ipoin) % shapl
        limit =  lnint(ipoin) % limit

        do idime = 1,ndime
           Auud(idime) = 1.0_rp
        end do
        !
        ! Set line to zero
        !
        do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
           jpoin = c_dom(izdom)
           do idime = 1,ndime
              do jdime = 1,ndime
                 Auu(jdime,idime,izdom) = 0.0_rp
              end do
              Aup(idime,izdom) = 0.0_rp                    
           end do
        end do
        !
        ! Set the force term to zero
        !
        do idime = 1,ndime
           bu(idime,ipoin) = 0.0_rp
        end do
        !
        ! Find the projection point on the particle surface
        !
        iimbo    = abs(lntib(ipoin))
        if (limit > 1) then
           call dpopar(ndime,coord(1:ndime,ipoin),&
                imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                dummr,dumma,propo,dummi,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                imbou(iimbo) % lnele)

        elseif (limit == 1) then
           do idime = 1,ndime
              coor1(idime) =  coord(idime,ipoin) - 2.0_rp*(coord(idime,lnint(ipoin)%lnode(1)) - coord(idime,ipoin))
           end do

           call faceli(&
                imbou(iimbo) % sabox,imbou(iimbo) % blink,imbou(iimbo) %ltyib,imbou(iimbo) %lnoib,imbou(iimbo) %cooib, & 
                coor1,coord(1,lnint(ipoin)%lnode(1)),coord(1,ipoin),mnoib,imbou(iimbo) %nboib,propo,dummi,dummr)              
        end if
        !
        ! Determine the angular velocity 
        !
        x(1)     = propo(1) - imbou(iimbo) % posil(    1,1)
        x(2)     = propo(2) - imbou(iimbo) % posil(    2,1)
        x(3)     = 0.0_rp
        if ( ndime == 3)  x(ndime) = propo(ndime) - imbou(iimbo) % posil(ndime,1)
        v(1)     = 0.0_rp
        v(2)     = 0.0_rp
        v(3)     = 0.0_rp
        call vecpro(imbou(iimbo)%veloa,x,v,3_ip)            
        v(1)     = v(1)     + imbou(iimbo) % velol(    1,1)
        v(2)     = v(2)     + imbou(iimbo) % velol(    2,1)
        if ( ndime == 3) v(3) = v(3) + imbou(iimbo) % velol(3,1)

        !v(1) =  sin(pi*propo(1)-0.7_rp)*sin(pi*propo(2)+0.2_rp)
        !v(2) =  cos(pi*propo(1)-0.7_rp)*cos(pi*propo(2)+0.2_rp)
        !v(1) =  propo(1)
        !v(2) = -propo(2)
        

        !----------------------------------------------------------------------
        !
        ! Normal interpolation
        ! Set velocity line equal to kriging coefficients
        !
        !----------------------------------------------------------------------    
        do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
           jpoin = c_dom(izdom)
           if (ipoin /= jpoin) then
              kpoin = 0
              iinde = 0
              do while ( kpoin /= jpoin .and. iinde < limit )
                 iinde = iinde + 1_ip
                 kpoin = lnode(iinde)                       
              end do
              if ( kpoin == jpoin ) then  
                 do idime = 1,ndime        
                    Auu(idime,idime,izdom) = -shapl(iinde) * Auud(idime)
                 end do
              end if
           !
           ! Set the velocity line only if the current node belong to my subdomain
           !
           elseif (ipoin <= npoi1 .or. (ipoin >= npoi2 .and. ipoin <= npoi3)) then
              do idime = 1,ndime                                  
                 Auu(idime,idime,izdom) = Auud(idime)
              end do
           end if
        end do
        !----------------------------------------------------------------------
        !
        ! Constrained interpolation.
        ! Set the right-hand side only if the current node belong to my subdomain.
        !
        !----------------------------------------------------------------------   
        if (ipoin <= npoi1 .or. (ipoin >= npoi2 .and. ipoin <= npoi3)) then 
           do idime = 1,ndime
              bu(idime,ipoin) = v(idime)*shapl(limit+1)*Auud(idime) +  &
                   (1.0_rp/massc(1,ipoin))*massc(idime+1,ipoin) * & 
                   (1.0_rp/imbou(iimbo) % mass1)*(0.0_rp-imbou(iimbo) % mass2)*Auud(idime)
           end do
        end if     
     end if
  end do

end subroutine nsi_embedd
