!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine outset(itask,lun_setsa,nvars,npp_setsa,woase,vaset)
  !-----------------------------------------------------------------------
  !****f* outrut/outset
  ! NAME 
  !    nsi_outset
  ! DESCRIPTION
  !    Output sets values
  !    ITASK = 1 ... Element sets
  !            2 ... Boundary sets
  !            3 ... Node sets
  !            4 ... IB sets
  !            5 ... Witness
  !    In the case of node sets, this subroutine must be called with
  !    nvars-1 as the dimension of vnset is vnset(nvars,nnset).
  ! USES
  ! USED BY
  !    ***_outset
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  ittim,itcou,itinn,cutim,ioutp,coutp,&
       &                      modul,routp,INOTSLAVE,kfl_rstar,&
       &                      ITASK_ENDRUN
  use def_kermod, only     :  nwitn,nwitg
  use def_domain, only     :  neset,nbset,nnset,&
       &                      lesec,lbsec,lnsec
  use mod_outfor, only     :  outfor
  use mod_iofile, only     :  iofile_flush_unit  
  implicit none
  integer(ip),  intent(in) :: itask,lun_setsa,nvars,npp_setsa(*)
  real(rp),     intent(in) :: vaset(nvars+1,*)
  character(5), intent(in) :: woase(*)
  integer(ip)              :: ipase,ipasb,ipasn,ivars,ieset,jbset,inset
  integer(ip)              :: iwitn
 
  if( INOTSLAVE ) then

     if(itask==-1.and.kfl_rstar/=2.and.maxval(npp_setsa(1:nvars))>0) then
        !
        ! Element sets: Header
        !
        ipase=0
        coutp(1)='Element'
        call outfor(34_ip,lun_setsa,' ')
        ipase=ipase+1
        ioutp(1)=ipase
        coutp(1)='ISET '
        call outfor(35_ip,lun_setsa,' ')
        ipase=ipase+1
        ioutp(1)=ipase
        coutp(1)='SIZES'
        call outfor(35_ip,lun_setsa,' ')
        do ivars=1,nvars
           if(npp_setsa(ivars)/=0) then
              ipase=ipase+1
              ioutp(1)=ipase
              coutp(1)=woase(ivars)
              call outfor(35_ip,lun_setsa,' ')
           end if
        end do
        ioutp(1)=ipase
        ioutp(2)=neset
        call outfor(33_ip,lun_setsa,' ')

     else if(itask==1) then
        !
        ! Element sets: Iteration header
        !
        ioutp(1)=ittim
        ioutp(2)=itcou
        ioutp(3)=itinn(modul)
        routp(1)=cutim
        call outfor(32_ip,lun_setsa,' ')

     else if(itask==10) then
        !
        ! Element sets: Write results
        !
        do ieset=1,neset
           write(lun_setsa,5,advance='no') lesec(ieset)
           write(lun_setsa,1,advance='no') vaset(nvars+1,ieset)
           do ivars=1,nvars
              if(npp_setsa(ivars)/=0)& 
                   write(lun_setsa,1,advance='no') vaset(ivars,ieset)
           end do
           write(lun_setsa,*)
        end do
        call iofile_flush_unit(lun_setsa)

     else if(itask==-2.and.kfl_rstar/=2.and.maxval(npp_setsa(1:nvars))>0) then
        !
        ! Boundary sets: File header
        !
        ipasb=0
        coutp(1)='Boundary'
        call outfor(34_ip,lun_setsa,' ')
        ipasb=ipasb+1
        ioutp(1)=ipasb
        coutp(1)='ISET '
        call outfor(35_ip,lun_setsa,' ')
        ipasb=ipasb+1
        ioutp(1)=ipasb
        coutp(1)='SIZES'
        call outfor(35_ip,lun_setsa,' ')
        do ivars=1,nvars
           if(npp_setsa(ivars)/=0) then
              ipasb=ipasb+1
              ioutp(1)=ipasb
              coutp(1)=woase(ivars)
              call outfor(35_ip,lun_setsa,' ')
           end if
        end do
        ioutp(1)=ipasb
        ioutp(2)=nbset
        call outfor(33_ip,lun_setsa,' ')

     else if(itask==2) then
        !
        ! Boundary sets: Iteration header
        !
        ioutp(1)=ittim
        ioutp(2)=itcou
        ioutp(3)=itinn(modul)
        routp(1)=cutim
        call outfor(32_ip,lun_setsa,' ')

     else if(itask==20) then
        !
        ! Boundary sets: Write results
        !
        do jbset=1,nbset

           write(lun_setsa,5,advance='no') lbsec(jbset)
           write(lun_setsa,1,advance='no') vaset(nvars+1,jbset)
           do ivars=1,nvars   
              if(npp_setsa(ivars)/=0)& 
                   write(lun_setsa,1,advance='no') vaset(ivars,jbset)
           end do           
           write(lun_setsa,*)
        end do
        call iofile_flush_unit(lun_setsa)

     else if(itask==-3.and.kfl_rstar/=2.and.maxval(npp_setsa(1:nvars+1))>0) then
        !
        ! Node sets: File header
        !
        ipasn=0
        coutp(1)='Node'
        call outfor(34_ip,lun_setsa,' ')
        ipasn=ipasn+1
        ioutp(1)=ipasn
        coutp(1)='ISET '
        call outfor(35_ip,lun_setsa,' ')
        do ivars=1,nvars+1
           if(npp_setsa(ivars)/=0) then
              ipasn=ipasn+1
              ioutp(1)=ipasn
              coutp(1)=woase(ivars)
              call outfor(35_ip,lun_setsa,' ')
           end if
        end do
        ioutp(1)=ipasn
        ioutp(2)=nnset
        call outfor(33_ip,lun_setsa,' ')

     else if(itask==3) then
        !
        ! Node sets: Iteration header
        !
        ioutp(1)=ittim
        ioutp(2)=itcou
        ioutp(3)=itinn(modul)
        routp(1)=cutim
        call outfor(32_ip,lun_setsa,' ')

     else if(itask==30) then
        !
        ! Node sets: Write results
        !
        do inset=1,nnset
           write(lun_setsa,5,advance='no') lnsec(inset)
           do ivars=1,nvars+1
              if(npp_setsa(ivars)/=0)& 
                   write(lun_setsa,1,advance='no') vaset(ivars,inset)
           end do
           write(lun_setsa,*)
        end do
        call iofile_flush_unit(lun_setsa)

     else if(itask==-5.and.kfl_rstar/=2.and.maxval(npp_setsa(1:nvars+1))>0) then
        !
        ! Wintess nodes: File header
        !
        ipasn=0
        coutp(1)='Witness'
        call outfor(34_ip,lun_setsa,' ')
        ipasn=ipasn+1
        ioutp(1)=ipasn
        coutp(1)='ISET '
        call outfor(35_ip,lun_setsa,' ')
        do ivars=1,nvars+1
           if(npp_setsa(ivars)/=0) then
              ipasn=ipasn+1
              ioutp(1)=ipasn
              coutp(1)=woase(ivars)
              call outfor(35_ip,lun_setsa,' ')
           end if
        end do
        ioutp(1)=ipasn
        ioutp(2)=nwitn
        call outfor(33_ip,lun_setsa,' ')

     else if(itask==-6.and.kfl_rstar/=2.and.maxval(npp_setsa(1:nvars+1))>0) then
        !
        ! Wintess nodes: File header
        !
        ipasn=0
        coutp(1)='Witness geometry'
        call outfor(34_ip,lun_setsa,' ')
        ipasn=ipasn+1
        ioutp(1)=ipasn
        coutp(1)='ISET '
        call outfor(35_ip,lun_setsa,' ')
        do ivars=1,nvars+1
           if(npp_setsa(ivars)/=0) then
              ipasn=ipasn+1
              ioutp(1)=ipasn
              coutp(1)=woase(ivars)
              call outfor(35_ip,lun_setsa,' ')
           end if
        end do
        ioutp(1)=ipasn
        ioutp(2)=nwitg
        call outfor(33_ip,lun_setsa,' ')

     else if(itask==5) then
        !
        ! Witness pointes: Iteration header
        !
        ioutp(1)=ittim
        ioutp(2)=itcou
        ioutp(3)=itinn(modul)
        routp(1)=cutim
        call outfor(32_ip,lun_setsa,' ')

     else if(itask==6) then
        !
        ! Witness points: Iteration header
        !
        ioutp(1)=ittim
        ioutp(2)=itcou
        ioutp(3)=itinn(modul)
        routp(1)=cutim
        call outfor(32_ip,lun_setsa,' ')

     else if(itask==50) then
        !
        ! Witness point: Write results
        !
        do iwitn=1,nwitn
           write(lun_setsa,5,advance='no') iwitn
           do ivars=1,nvars+1
              if(npp_setsa(ivars)/=0)& 
                   write(lun_setsa,1,advance='no') vaset(ivars,iwitn)
           end do
           write(lun_setsa,*)
        end do
        call iofile_flush_unit(lun_setsa)

     else if(itask==60) then
        !
        ! Witness geometry: Write results
        !
        do iwitn=1,nwitg
           write(lun_setsa,5,advance='no') iwitn
           do ivars=1,nvars+1
              if(npp_setsa(ivars)/=0) then
                 write(lun_setsa,1,advance='no') vaset(ivars,iwitn)
              end if
           end do
           write(lun_setsa,*)
        end do
        call iofile_flush_unit(lun_setsa)

     end if

  end if

1 format(1x,e16.8E3)
5 format(1x,i9)
100 format (10(e16.8E3,','))
end subroutine outset
