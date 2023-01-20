!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_scatte()
  !-----------------------------------------------------------------------
  !****f* Parall/par_scatte
  ! NAME
  !    par_scatte
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use def_master
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip)          :: ipoin
  integer(ip)          :: jpoin
  integer(ip)          :: ii,i1
!  integer(ip)          :: i2,i3
  integer(4)           :: istat
  integer(ip), pointer :: loc_pari1(:)
  real(rp),    pointer :: loc_parr1(:),loc_parr2(:,:)
!  real(rp),    pointer :: loc_parr3(:,:,:)
  !
  ! Initialize
  !
  npari=0
  nparr=0
  nparc=0
  if(ISLAVE) kfl_desti_par=0

  if(party==1) then
     !
     ! Element
     !
     if(pardi==1.and.parki==1) then

     end if

  else if(party==2) then
     !
     ! Boundary
     !
     if(pardi==1.and.parki==1) then


     else if(pardi==1.and.parki==2) then


     else if(pardi==3.and.parki==2) then


     end if

  else if(party==3) then
     !
     ! Node
     !
     if(pardi==1.and.parki==1) then

        if(IMASTER) then 
           allocate(pari1(npoin),stat=istat)
           call memchk(zero,istat,par_memor,'PARI1','par_scatte',pari1)
           ipoin=1
           do kfl_desti_par= 1, npart_par
              allocate(loc_pari1(npoin_par(kfl_desti_par)),stat=istat)
              call memchk(zero,istat,par_memor,'LOC_PARI1','par_scatte',loc_pari1)
              npari =  npoin_par(kfl_desti_par)
              parin => loc_pari1
              call par_receiv()
              do ii=1,npoin_par(kfl_desti_par)
                 jpoin=lninv_loc(ipoin)
                 pari1(jpoin)=loc_pari1(ii)
                 ipoin=ipoin+1
              end do
              call memchk(two,istat,par_memor,'LOC_PARI1','par_scatte',loc_pari1)
              deallocate(loc_pari1,stat=istat)
              if(istat/=0) call memerr(two,'LOC_PARI1','par_scatte',0_ip)
           end do
        else
           npari =  npoin
           parin => pari1
           call par_sendin()
        end if

     else if(pardi==1.and.parki==2) then

        if(IMASTER) then 
           allocate(parr1(npoin),stat=istat)
           call memchk(zero,istat,par_memor,'PARR1','par_scatte',parr1)
           ipoin=1
           do kfl_desti_par= 1, npart_par
              allocate(loc_parr1(npoin_par(kfl_desti_par)),stat=istat)
              call memchk(zero,istat,par_memor,'LOC_PARR1','par_scatte',loc_parr1)
              nparr =  npoin_par(kfl_desti_par)
              parre => loc_parr1
              call par_receiv()
              do ii=1,npoin_par(kfl_desti_par)
                 jpoin=lninv_loc(ipoin)
                 parr1(jpoin)=loc_parr1(ii)
                 ipoin=ipoin+1
              end do
              call memchk(two,istat,par_memor,'LOC_PARR1','par_scatte',loc_parr1)
              deallocate(loc_parr1,stat=istat)
              if(istat/=0) call memerr(two,'LOC_PARR1','par_scatte',0_ip)
           end do
        else
           nparr =  npoin
           parre => parr1
           call par_sendin()
        end if

     else if(pardi==2.and.parki==2) then

        if(IMASTER) then 
           allocate(parr2(ndime,npoin),stat=istat)
           call memchk(zero,istat,par_memor,'PARR2','mod_postpr',parr2)
           ipoin=1
           do kfl_desti_par= 1, npart_par
              allocate(loc_parr2(pard1,npoin_par(kfl_desti_par)),stat=istat)
              call memchk(zero,istat,par_memor,'LOC_PARR2','par_scatte',loc_parr2)
              nparr =  pard1*npoin_par(kfl_desti_par)
              call par_srreal(2_ip,nparr,loc_parr2)
              do ii=1,npoin_par(kfl_desti_par)
                 jpoin=lninv_loc(ipoin)
                 do i1=1,pard1
                    parr2(i1,jpoin)=loc_parr2(i1,ii)
                 end do
                 ipoin=ipoin+1
              end do
              call memchk(two,istat,par_memor,'LOC_PARR2','par_scatte',loc_parr2)
              deallocate(loc_parr2,stat=istat)
              if(istat/=0) call memerr(two,'LOC_PARR2','par_scatte',0_ip)
           end do
        else
           nparr =  pard1*npoin
           call par_srreal(1_ip,nparr,parr2)
        end if

     else if(pardi==2.and.parki==2) then


     end if

  else if(party==4) then
     !
     ! Boundary node: REAL ARRAY(NBOPO)
     !
     call runend('PAR_SCATTE: NO LONGER SUPPORTED')
!!$     if(pardi==3.and.parki==2) then
!!$        if(IMASTER) then
!!$           ii = 0
!!$           do kfl_desti_par= 1, npart_par
!!$              nparr = ndime*ndime*nbopo_par(kfl_desti_par)
!!$              allocate(loc_parr3(pard1,pard2,nbopo_par(kfl_desti_par)),stat=istat)
!!$              call memchk(zero,istat,par_memor,'LOC_PARR3','par_scatte',loc_parr3)
!!$              call par_receiv()
!!$              do i1= 1, nbopo_par(kfl_desti_par)
!!$                 do i2= 1, ndime
!!$                    do i3= 1, ndime
!!$                       parr3(i3,i2,i1+ii) = loc_parr3(i3,i2,i1)
!!$                    end do
!!$                 end do
!!$              end do
!!$              call memchk(two,istat,par_memor,'LOC_PARR3','par_scatte',loc_parr3)
!!$              deallocate(loc_parr3,stat=istat)
!!$              if(istat/=0) call memerr(two,'LOC_PARR3','par_scatte',0_ip)
!!$              ii = ii + nbopo_par(kfl_desti_par)
!!$           end do
!!$        else
!!$           nparr =  pard1*pard2*nbopo
!!$           call par_srreal(1_ip,nparr,parr3)
!!$        end if
!!$     end if

  end if

  npari=0
  nparr=0
  nparc=0

end subroutine par_scatte
