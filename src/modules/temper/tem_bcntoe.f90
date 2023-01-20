!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_bcntoe
!-----------------------------------------------------------------------
!        
! This routine transforms the boundary conditions on nodes of 
! Neumann(Robin) type to Neumann(Robin) conditions on boundaries.
!G?
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_temper
  use      mod_memchk
  implicit none
  integer(ip) :: ipoin,pnodb,iboun,inodb,iffix
  integer(ip) :: knodb(10)            ! Number of bc appearance
  real(rp)    :: vafix(10),value

  return

  do iboun=1,nboun

     ! Count number of appearance for each fixity on iboun
     ! e.g. knodb(2)= number of Neumann nodes
     !      knodb(3)= number of Robin nodes
     pnodb=nnode(ltypb(iboun))
     knodb=0
     vafix=0.0_rp
     do inodb=1,pnodb
        ipoin=lnodb(inodb,iboun)
        iffix=kfl_fixno_tem(1,ipoin)
        if(iffix>1) then                             !G? if(iffix>0) then
           knodb(iffix)=knodb(iffix)+1
           vafix(iffix)=vafix(iffix)+bvess_tem(1,ipoin,1)
        end if
     end do
     where (knodb/=0) vafix=vafix/real(knodb,rp)

     if(knodb(2)>=1) then
        !
        ! Neumann or Robin
        !
        kfl_fixbo_tem(iboun)=2
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           if(kfl_fixno_tem(1,ipoin)==2) then
              value=bvess_tem(1,ipoin,1)
           else
              value=vafix(2)
           end if
           bvnat_tem(1,iboun,1)=value
        end do
     end if

     if(knodb(3)>=1) then
        !
        ! Robin
        !
        kfl_fixbo_tem(iboun)=3
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           if(kfl_fixno_tem(1,ipoin)==3) then
              value=bvess_tem(1,ipoin,1)
           else
              value=vafix(3)
           end if
           bvnat_tem(1,iboun,1)=scond_tem
           bvnat_tem(2,iboun,1)=0.0_rp
           bvnat_tem(3,iboun,1)=value
        end do
     end if

  end do

end subroutine tem_bcntoe
