!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine finbou(pnodb,knodb,iboun)
!-----------------------------------------------------------------------
!****f* Domain/finbou
! NAME 
!    finbou
! DESCRIPTION
!    This routine looks for iboun in the list lnodb which corresponds
!    to the list knodb. This is usefull when automatic boundaries are
!    generated. In this case, the numbering of the boundaries is
!    a priori unknown as it is calculated internally by Alya.
! USES
! USED BY
!    nsi_reabcs
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  implicit none
  integer(ip), intent(in)  :: pnodb
  integer(ip), intent(in)  :: knodb(pnodb)
  integer(ip), intent(out) :: iboun
  integer(ip)              :: jpoin,inodb,pblty,onodb,inook,ibook
  integer(ip)              :: qnodb,jnodb

  iboun=0
  ibook=0
  loop_iboun: do while(iboun<nboun)
     iboun=iboun+1
     pblty=ltypb(iboun) 
     qnodb=nnode(pblty)
     if(pnodb==qnodb) then
        onodb=0
        jnodb=0
        inook=0
        loop_jnodb: do while(jnodb<pnodb)
           jnodb=jnodb+1
           jpoin=knodb(jnodb)
           inodb=0
           loop_inodb: do while(inodb<qnodb)
              inodb=inodb+1
              if(jpoin==lnodb(inodb,iboun)) then
                 onodb=onodb+1
                 inook=1
                 exit loop_inodb
              end if
           end do loop_inodb
           if(inook/=1) exit loop_jnodb
        end do loop_jnodb
        if(onodb==pnodb) then
           ibook=1
           exit loop_iboun
        end if
     end if
  end do loop_iboun
  if(ibook/=1) iboun=0

end subroutine finbou

subroutine finbou_fast(pnodb,knodb,iboun,lnodb_tmp)
!-----------------------------------------------------------------------
!****f* Domain/finbou
! NAME 
!    finbou
! DESCRIPTION
!    This routine looks for iboun in the list lnodb which corresponds
!    to the list knodb. This is usefull when automatic boundaries are
!    generated. In this case, the numbering of the boundaries is
!    a priori unknown as it is calculated internally by Alya.
! USES
! USED BY
!    nsi_reabcs
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  implicit none
  integer(ip), intent(in)    :: pnodb
  integer(ip), intent(in)    :: knodb(pnodb)
  integer(ip), intent(out)   :: iboun
  integer(ip), intent(inout) :: lnodb_tmp(mnodb,nboun)
  integer(ip)                :: inodb
  logical(lg)                :: ifoun

  iboun = 1
  loop_iboun: do while( iboun <= nboun )
     if( lnodb_tmp(1,iboun) /= 0 ) then
        ifoun = .true.
        inodb_loop: do inodb = 1,pnodb
           if( knodb(inodb) /= lnodb_tmp(inodb,iboun) ) then
              ifoun = .false.
              exit inodb_loop 
           end if
        end do inodb_loop
        if( ifoun ) then
           lnodb_tmp(1,iboun) = 0
           return
        end if
     end if
     iboun = iboun + 1
  end do loop_iboun

end subroutine finbou_fast
