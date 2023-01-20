!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_senset()
  !------------------------------------------------------------------------
  !****f* Parall/par_senset
  ! NAME
  !    par_senset
  ! DESCRIPTION
  !    Send/receive sets and witness points
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_parall
  use def_domain
  use def_master
  use mod_memchk
  use mod_domain, only : domain_memory_allocate
  use mod_domain, only : domain_memory_deallocate
  implicit none  
  integer(ip)          :: jpoin,ipoin
  integer(ip)          :: jelem,jboun,iboun,ielem,indice0

  if( IMASTER ) then
     
     if( neset > 0 ) then

        !----------------------------------------------------------------
        !
        ! Element sets
        !
        !----------------------------------------------------------------

        call memgen(1_ip,nelem,0_ip)
        do ielem = 1,nelem
           jelem        = leper_par(ielem)
           gisca(jelem) = leset(ielem)
        end do
        do kfl_desti_par = 1,npart_par
           npari =  nelem_par(kfl_desti_par)
           parin => gisca(leind_par(kfl_desti_par):)
           strin =  'GISCA'
           call par_sendin()
        end do
        call memgen(3_ip,nelem,0_ip)
        call domain_memory_deallocate('LESET')
        
     end if

     if( nbset > 0 ) then

        !----------------------------------------------------------------
        !
        ! Boundary sets
        !
        !----------------------------------------------------------------

        call memgen(1_ip,nboun,0_ip)
        do iboun = 1,nboun
           jboun        = lbper_par(iboun)
           gisca(jboun) = lbset(iboun)
        end do 
        do kfl_desti_par = 1,npart_par
           npari =  nboun_par(kfl_desti_par)
           parin => gisca(lbind_par(kfl_desti_par):)
           strin =  'LBSET'
           call par_sendin()
        end do 
        call memgen(3_ip,nboun,0_ip)
        call domain_memory_deallocate('LBSET')

     end if

     if( nnset > 0 ) then
        
        !----------------------------------------------------------------
        !
        ! Node sets
        !
        !----------------------------------------------------------------

        call memgen(1_ip,npoin_total,0_ip)
        do ipoin = 1,npoin_total
           jpoin        = lninv_loc(ipoin)
           gisca(ipoin) = lnset(jpoin)
        end do
        indice0 = 1
        do kfl_desti_par = 1,npart_par
           npari =  npoin_par(kfl_desti_par)
           parin => gisca(indice0:)
           strin =  'LNSET'
           call par_sendin()
           indice0 = indice0 + npoin_par(kfl_desti_par)
        end do
        call memgen(3_ip,npoin,0_ip)
        call domain_memory_deallocate('LNSET')

     end if
     
  else if( ISLAVE ) then

     kfl_desti_par = 0

     if( neset > 0 ) then
        !
        ! Element sets
        !
        call domain_memory_allocate('LESET')
        npari =  nelem
        parin => leset(:)
        call par_receiv()
     end if

     if( nbset > 0 ) then
        !
        ! Boundary sets
        !
        call domain_memory_allocate('LBSET')
        npari =  nboun
        parin => lbset(:)
        call par_receiv()
     end if

     if( nnset > 0 ) then
        !
        ! Node sets 
        !
        call domain_memory_allocate('LNSET')
        npari =  npoin
        parin => lnset(:)
        call par_receiv()
     end if

  end if

end subroutine par_senset
