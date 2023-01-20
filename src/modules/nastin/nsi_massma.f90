!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_massma()

  use def_kintyp
  use def_parame
  use def_domain
  use def_master
  use def_nastin
  use mod_kdtree
  use mod_communications_global, only : PAR_SUM

  implicit none

  integer(ip)                :: iimbo,ipoin,jpoin,idime,inode,limit
  integer(ip)                :: dummi
  integer(ip), pointer       :: lnode(:)
  real(rp),    pointer       :: shapl(:)
  real(rp)                   :: x(3),v(3),dummr,dumma(ndime)
  real(rp)                   :: propo(ndime),coor1(ndime)
  
  do iimbo = 1,nimbo    
     !----------------------------------------------------------------------
     !
     ! Constrained interpolation matrices
     !
     !----------------------------------------------------------------------
     imbou(iimbo) % mass1 = 0.0_rp
     imbou(iimbo) % mass2 = 0.0_rp
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           if (lntib(ipoin) == -iimbo) then     
              lnode => lnint(ipoin) % lnode
              shapl => lnint(ipoin) % shapl
              limit =  lnint(ipoin) % limit              
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
              !if ( ndime == 3) v(3) = v(3) + imbou(iimbo) % velol(3,1)
              !v(1) =  sin(pi*propo(1)-0.7_rp)*sin(pi*propo(2)+0.2_rp)
              !v(2) =  cos(pi*propo(1)-0.7_rp)*cos(pi*propo(2)+0.2_rp)
              !v(1) =   propo(1)
              !v(2) =  -propo(2)
              !
              ! ipoin is only in my subdomain
              !
              if (ipoin <= npoi1 .or. (ipoin>= npoi2 .and. ipoin<= npoi3)) then
                 do idime=1,ndime
                    imbou(iimbo) % mass1 = imbou(iimbo) % mass1 + &
                         (1.0_rp/massc(1,ipoin)) * massc(idime+1,ipoin)*massc(idime+1,ipoin)
                    !
                    ! R_in^t(ndime*npoin)*(I_nm(npoin,nnode)*b_nj(ndime*npoin))
                    !                       
                    imbou(iimbo) % mass2 = imbou(iimbo) % mass2 + &
                         massc(idime+1,ipoin) * (shapl(limit+1)*v(idime))                    
                 end do
              end if
              do inode = 1,limit
                 jpoin = lnode(inode)
                 !
                 ! jpoin is my responsability
                 !                    
                 if (jpoin <= npoin) then
                    do idime=1,ndime                       
                       !
                       ! R_in^t(ndime*npoin)*(I_nm(npoin,nnode)*a_mj(nnode,ndime*npoin))
                       !
                       imbou(iimbo) % mass2 = imbou(iimbo) % mass2 + &
                            massc(idime+1,ipoin) * (shapl(inode)*veloc(idime,jpoin,1))
                    end do
                 end if
              end do
           end if
        end do
     end if
     call PAR_SUM(imbou(iimbo) % mass1)
     call PAR_SUM(imbou(iimbo) % mass2)
  end do

end subroutine nsi_massma
