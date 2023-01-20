!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_scafil()
  !-----------------------------------------------------------------------
  !****f* Parall/par_scafil
  ! NAME
  !    par_scafil
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves using
  !    a filter.
  !    The filter is saved in PARI1 and should be deallocated afterwards
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip)          :: ipoin,jpoin,jfilt,pfilt,ii,kpoin,idime
  integer(4)           :: istat
  real(rp),    pointer :: loc_parr1(:),loc_parr2(:,:)
  integer(ip), pointer :: perfi(:)
  integer(ip), target  :: dummi_par(npart_par)
  !
  ! Initialize
  !
  npari = 0
  nparr = 0
  nparc = 0
  jfilt = 0
  if( ISLAVE ) kfl_desti_par = 0

  !----------------------------------------------------------------------
  !
  ! Count permutation
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) then

     do ipoin = 1,npoi1
        if( gefil(ipoin) /= 0 ) jfilt = jfilt + 1
     end do
     do ipoin = npoi2,npoi3
        if( gefil(ipoin) /= 0 ) jfilt = jfilt + 1
     end do
     dummi_par(1) =  jfilt
     npari        =  1
     parin        => dummi_par(1:)
     call par_sendin()

  else

     do kfl_desti_par = 1, npart_par
        npari =  1
        parin => dummi_par(kfl_desti_par:)
        call par_receiv()
     end do
     pfilt = 0
     do kfl_desti_par = 1, npart_par
        pfilt = pfilt + dummi_par(kfl_desti_par)
     end do

  end if

  if(party==3) then

     if(pardi==1.and.parki==2) then

        !----------------------------------------------------------------------
        !
        ! Compute and send permutation: REAL(NPOIN)
        !
        !----------------------------------------------------------------------

        if( ISLAVE ) then

           allocate(perfi(jfilt),stat=istat)
           call memchk(zero,istat,par_memor,'PERFI','par_scafil',perfi)

           allocate(loc_parr1(jfilt),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_scafil',loc_parr1)

           jfilt = 0
           do ipoin = 1,npoi1
              if( gefil(ipoin) /= 0 ) then
                 jfilt            = jfilt + 1
                 perfi(jfilt)     = ipoin
                 loc_parr1(jfilt) = parr1(ipoin)
              end if
           end do
           do ipoin = npoi2,npoi3
              if( gefil(ipoin) /= 0 ) then
                 jfilt            = jfilt + 1
                 perfi(jfilt)     = ipoin
                 loc_parr1(jfilt) = parr1(ipoin)
              end if
           end do

           npari =  jfilt
           parin => perfi
           call par_sendin()
           nparr =  jfilt
           parre => loc_parr1
           call par_sendin()

           call memchk(two,istat,par_memor,'LOC_PARR1','par_scafil',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR1','par_scafil',0_ip)
           call memchk(two,istat,par_memor,'PERFI','par_scafil',perfi)
           deallocate(perfi,stat=istat)
           if(istat/=0) call memerr(two,'PERFI','par_scafil',0_ip)

        else

           allocate(parr1(pfilt),stat=istat)
           call memchk(zero,istat,par_memor,'PARR1','par_scafil',parr1)
           allocate(pari1(pfilt),stat=istat)
           call memchk(zero,istat,par_memor,'PARI1','par_scafil',pari1)
           ipoin = 1

           kpoin = 0

           do kfl_desti_par = 1, npart_par

              jfilt = dummi_par(kfl_desti_par)
              allocate(perfi(jfilt),stat=istat)
              call memchk(zero,istat,par_memor,'PERFI','par_scafil',perfi)        
              allocate(loc_parr1(jfilt),stat=istat)
              call memchk(zero,istat,par_memor,'LOC_PARR1','par_scafil',loc_parr1)

              npari =  jfilt       ! Receive permuation
              parin => perfi
              call par_receiv()
              nparr =  jfilt       ! Receive array
              parre => loc_parr1
              call par_receiv()        

              do ii = 1,jfilt
                 jpoin        = lninv_loc(kpoin+perfi(ii))
                 pari1(ipoin) = jpoin
                 parr1(ipoin) = loc_parr1(ii)
                 ipoin        = ipoin + 1
              end do
              kpoin = kpoin + npoin_par(kfl_desti_par)

              call memchk(two,istat,par_memor,'LOC_PARR1','par_scafil',loc_parr1)
              deallocate(loc_parr1,stat=istat)
              if( istat /= 0 ) call memerr(two,'LOC_PARR1','par_scafil',0_ip)
              call memchk(two,istat,par_memor,'PERFI','par_scafil',perfi)
              deallocate(perfi,stat=istat)
              if( istat /= 0 ) call memerr(two,'PERFI','par_scafil',0_ip)

           end do

        end if

     else if(pardi==2.and.parki==2) then

        !----------------------------------------------------------------------
        !
        ! Compute and send permutation: REAL(NDIME,NPOIN)
        !
        !----------------------------------------------------------------------

        if( ISLAVE ) then

           allocate(perfi(jfilt),stat=istat)
           call memchk(zero,istat,par_memor,'PERFI','par_scafil',perfi)

           allocate(loc_parr2(ndime,jfilt),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR2','par_scafil',loc_parr2)

           jfilt = 0
           do ipoin = 1,npoi1
              if( gefil(ipoin) /= 0 ) then
                 jfilt            = jfilt + 1
                 perfi(jfilt)     = ipoin
                 do idime = 1,ndime
                    loc_parr2(idime,jfilt) = parr2(idime,ipoin)
                 end do
              end if
           end do
           do ipoin = npoi2,npoi3
              if( gefil(ipoin) /= 0 ) then
                 jfilt            = jfilt + 1
                 perfi(jfilt)     = ipoin
                 do idime = 1,ndime
                    loc_parr2(idime,jfilt) = parr2(idime,ipoin)
                 end do
              end if
           end do

           npari =  jfilt
           parin => perfi
           call par_sendin()
           nparr =  jfilt*ndime
           call par_srreal(1_ip,nparr,loc_parr2)

           call memchk(two,istat,par_memor,'LOC_PARR2','par_scafil',loc_parr2)
           deallocate(loc_parr2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR2','par_scafil',0_ip)
           call memchk(two,istat,par_memor,'PERFI','par_scafil',perfi)
           deallocate(perfi,stat=istat)
           if(istat/=0) call memerr(two,'PERFI','par_scafil',0_ip)

        else

           allocate(parr2(ndime,pfilt),stat=istat)
           call memchk(zero,istat,par_memor,'PARR2','par_scafil',parr2)
           allocate(pari1(pfilt),stat=istat)
           call memchk(zero,istat,par_memor,'PARI1','par_scafil',pari1)
           ipoin = 1

           kpoin = 0

           do kfl_desti_par = 1, npart_par

              jfilt = dummi_par(kfl_desti_par)
              allocate(perfi(jfilt),stat=istat)
              call memchk(zero,istat,par_memor,'PERFI','par_scafil',perfi)        
              allocate(loc_parr2(ndime,jfilt),stat=istat)
              call memchk(zero,istat,par_memor,'LOC_PARR2','par_scafil',loc_parr2)

              npari =  jfilt             ! Receive permuation
              parin => perfi
              call par_receiv()
              nparr =  jfilt*ndime       ! Receive array
              call par_srreal(2_ip,nparr,loc_parr2)      

              do ii = 1,jfilt
                 jpoin        = lninv_loc(kpoin+perfi(ii))
                 pari1(ipoin) = jpoin
                 do idime = 1,ndime
                    parr2(idime,ipoin) = loc_parr2(idime,ii)
                 end do
                 ipoin        = ipoin + 1
              end do
              kpoin = kpoin + npoin_par(kfl_desti_par)

              call memchk(two,istat,par_memor,'LOC_PARR2','par_scafil',loc_parr2)
              deallocate(loc_parr2,stat=istat)
              if( istat /= 0 ) call memerr(two,'LOC_PARR2','par_scafil',0_ip)
              call memchk(two,istat,par_memor,'PERFI','par_scafil',perfi)
              deallocate(perfi,stat=istat)
              if( istat /= 0 ) call memerr(two,'PERFI','par_scafil',0_ip)

           end do

        end if

     end if

  end if

  npari=0
  nparr=0
  nparc=0

end subroutine par_scafil
