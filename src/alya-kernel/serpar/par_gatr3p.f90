!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_gatr3p(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_gather
  ! NAME
  !    par_gather
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
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip)                :: ielem,idime,igaus,itime,ii
  integer(4)                 :: istat
  integer(ip), save, target  :: isize(1)
  integer(ip)                :: ksize(npart_par)
  real(rp),    pointer       :: loc_parr1(:)
!  integer(ip), pointer       :: loc_pari1(:)
  !
  ! Initialize
  !
  npari = 0
  nparr = 0
  nparc = 0
  if( ISLAVE ) kfl_desti_par = 0

  !----------------------------------------------------------------------
  !
  ! Sizes
  !
  !----------------------------------------------------------------------

  if( IMASTER ) then
     !
     ! Master receives sizes
     !
     parii = 0
     do kfl_desti_par = 1,npart_par
        npari =  1
        parin => isize
        call par_receiv()
        ksize(kfl_desti_par) = isize(1)
        parii =  parii + isize(1)
     end do

  else if( ISLAVE ) then
     !
     ! Slaves send sizes
     !
     isize(1) =  0
     do ielem = 1,nelem
        isize(1) = isize(1) + ngaus(ltype(ielem))
     end do
     isize(1) = isize(1) * ndime * 2
     npari =  1
     parin => isize
     call par_sendin()

  end if

  if( itask == 3 ) return ! Master needs only size

  !----------------------------------------------------------------------
  !
  ! Slaves send to master
  !
  !----------------------------------------------------------------------

  if( IMASTER .and. itask == 1 ) then
     !
     ! Master receives array
     !
     ii = 1
     do kfl_desti_par = 1,npart_par
        nparr =  ksize(kfl_desti_par)
        parre => parr1(ii:)
        call par_receiv()
        ii    =  ii + ksize(kfl_desti_par)
     end do

  else if( ISLAVE .and. itask == 1 ) then
     !
     ! Slave send array
     !
     allocate(loc_parr1(isize(1)),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_PARR1','par_gather',loc_parr1)
     nparr = 0
     do ielem = 1,nelem
        do idime = 1,ndime
           do igaus = 1,ngaus(ltype(ielem))
              do itime = 1,2
                 nparr = nparr + 1
                 loc_parr1(nparr) = par3p(ielem)%a(idime,igaus,itime)
              end do
           end do
        end do
     end do
     parre => loc_parr1
     call par_sendin()
     call memchk(two,istat,par_memor,'LOC_PARR1','par_gather',loc_parr1)
     deallocate(loc_parr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_PARR1','par_gather',0_ip) 

  end if

  !----------------------------------------------------------------------
  !
  ! Master sends to slaves
  !
  !----------------------------------------------------------------------

  if( IMASTER .and. itask == 2 ) then
     !
     ! Master reorders
     !
     !allocate(loc_pari1(parii),stat=istat)
     !call memchk(zero,istat,par_memor,'LOC_PARI1','par_gather',loc_pari1)
     !call memgen(1_ip,nelem,0_ip)
     !gisca(1) = 1
     !do ielem = 2,nelem
     !   gisca(ielem) = gisca(ielem-1) + pari1(ielem-1) * ndime * 2
     !end do
     !nparr = 0
     !do ielem = 1,nelem
     !   jelem = leinv_par(ielem)
     !   ndofj = gisca(jelem) - 1
     !   do idime = 1,ndime
     !      do igaus = 1,ngaus(ltype(ielem))
     !         do itime = 1,2
     !            nparr = nparr + 1
     !            ndofj = ndofj + 1
     !            loc_parr1(nparr) = parr1(ndofj)
     !         end do
     !      end do
     !   end do
     !end do
     !
     ! Master sends
     !
     ii = 1
     do kfl_desti_par = 1,npart_par
        nparr =  ksize(kfl_desti_par)
        parre => parr1(ii:)
        !parre => loc_parr1(ii:)
        call par_sendin()
        ii    =  ii + ksize(kfl_desti_par)
     end do

     !call memgen(3_ip,nelem,0_ip)
     !call memchk(two,istat,par_memor,'LOC_PARI1','par_gather',loc_pari1)
     !deallocate(loc_pari1,stat=istat)
     !if(istat/=0) call memerr(two,'par_gather','LOC_PARI1',0_ip)           

  else if( ISLAVE .and. itask == 2 ) then
     !
     ! Slaves receive
     !
     allocate(loc_parr1(isize(1)),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_PARR1','par_gather',loc_parr1)
     nparr =  isize(1)
     parre => loc_parr1
     call par_receiv()
     nparr =  0
     do ielem = 1,nelem
        do idime = 1,ndime
           do igaus = 1,ngaus(ltype(ielem))
              do itime = 1,2
                 nparr = nparr + 1
                 par3p(ielem)%a(idime,igaus,itime) = loc_parr1(nparr)
              end do
           end do
        end do
     end do
     call memchk(two,istat,par_memor,'LOC_PARR1','par_gather',loc_parr1)
     deallocate(loc_parr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_PARR1','par_gather',0_ip) 

  end if


end subroutine par_gatr3p
