!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_dync01(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_dync01
  ! NAME 
  !    tem_dync01
  ! DESCRIPTION
  !    Marek: Cooling system + typical room
  ! USES
  ! USED BY
  !    tem_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_temper
  use def_domain
  use mod_messages, only : livinf
  use mod_messages, only : messages_live
  use mod_ecoute,   only : ecoute_set_read_unit
  use mod_ecoute,   only : ecoute_set_write_unit
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  integer(ip)             :: ipoin,icode,ibset,ioerr,ichan,lchan(16)
  real(rp)                :: surf,vold=0.0_rp
  real(rp),    target     :: rpoin(2)
  real(rp),    save       :: tinf,vinf,tout,vint
  character(200)          :: filna
  character(5)            :: winf

  select case(itask)

  case(1_ip)
     !
     ! Before a time step: reads input data
     !
     if( ISEQUEN .or. IMASTER ) then

        call livinf(59_ip,'OPEN  INPUT FILE DYNAMIC COUPLING',1_ip)
        icode=0
        ioerr=1
        icode0:do while(icode==0) 
           do while(ioerr/=0)
              open(unit=lun_dynin_tem,file=trim(fil_dynin_tem),status='old',iostat=ioerr)
           end do
           read(lun_dynin_tem,*,err=2,end=2) tinf
           read(lun_dynin_tem,*,err=2,end=2) vinf
           read(lun_dynin_tem,*,err=2,end=2) icode
           if(icode==1) exit icode0
           icode=0
           ioerr=1
2          close(lun_dynin_tem)
        end do icode0
        close(lun_dynin_tem,status='delete')

        write(lun_dynlo_tem,'(a)') '-------------------------------------------' 
        write(lun_dynlo_tem,*) 
        call livinf(59_ip,'CLOSE INPUT FILE DYNAMIC COUPLING',1_ip)
        write(lun_dynlo_tem,10) cutim
        write(lun_dynlo_tem,12) tinf
        write(lun_dynlo_tem,13) vinf 
        !
        ! Get file name according to vinf 
        !
        if(kfl_inter_tem==1 ) then
           call tem_dync01_finam(vinf,winf,vint)
           write(lun_dynlo_tem,20) vint
        end if

     end if
     if( IPARALL .and. kfl_inter_tem==1 ) then
        nparr    =  1
        rpoin(1) =  vint
        parre    => rpoin
        call par_broadc()
        vint     =  parre(1)
     end if
     !
     ! Interpolate velocity and turbulent viscosity
     !
     if(kfl_inter_tem==1 ) then
        if( vold/=vint ) then
           if(vint==0.0_rp) then
              kfl_advec_tem = 0
              if( ISEQUEN .or. IMASTER ) write(lun_dynlo_tem,15) 
           else
              kfl_advec_tem = 1
              if( ISEQUEN .or. IMASTER ) then
                 filna='results/'//trim(winf)//'.post.'
                 write(lun_dynlo_tem,31) trim(filna)//'msh'
                 write(lun_dynlo_tem,32) trim(filna)//'res'
                 write(lun_dynlo_tem,*) 
                 open(unit=90,file='interpolation.tem.dat',status='unknown')
                 write(90,'(a)') 'PHYSICAL_PROBLEM'  
                 write(90,'(a)') '  INTERPOLATE: FROM_2D_FLOW'  
                 write(90,'(a)') '    INCLUDE '//trim(filna)//'msh'
                 write(90,'(a)') '    INCLUDE '//trim(filna)//'res'
                 write(90,'(a)') '    CORNERS'
                 write(90,'(a)') '       1  -0.01 0.01 0.23'
                 write(90,'(a)') '       2  -0.01 0.05 0.23'
                 write(90,'(a)') '       3  -0.01 0.09 0.23' 
                 write(90,'(a)') '       4  -0.01 0.13 0.23'
                 write(90,'(a)') '       5  -0.01 0.01 0.16'  
                 write(90,'(a)') '       6  -0.01 0.05 0.16'
                 write(90,'(a)') '       7  -0.01 0.09 0.16'  
                 write(90,'(a)') '       8  -0.01 0.13 0.16'   
                 write(90,'(a)') '       9  -0.01 0.01 0.09'      
                 write(90,'(a)') '      10  -0.01 0.05 0.09'    
                 write(90,'(a)') '      11  -0.01 0.09 0.09'     
                 write(90,'(a)') '      12  -0.01 0.13 0.09'  
                 write(90,'(a)') '      13  -0.01 0.01 0.02'
                 write(90,'(a)') '      14  -0.01 0.05 0.02'
                 write(90,'(a)') '      15  -0.01 0.09 0.02'
                 write(90,'(a)') '      16  -0.01 0.13 0.02'
                 write(90,'(a)') '    END_CORNERS'
                 write(90,'(a)') '  END_INTERPOLATE'
                 write(90,'(a)') 'END_PHYSICAL_PROBLEM'
                 close(unit=90)
                 open(unit=90,file='interpolation.tem.dat',status='old')
                 
                 call ecoute_set_read_unit (90_ip)                    ! Reading file
                 call ecoute_set_write_unit(momod(modul) % lun_outpu) ! Writing file
              end if
              call livinf(59_ip,'INTERPOLATE FLOW VARIABLES FOR DYNAMIC COUPLING',1_ip)
              call runend('TEM_DYNC01: OBSOLETE OPTION')
              if( ISEQUEN .or. IMASTER ) close(unit=90)
           end if
        end if
        vold=vint
     end if

     if( ISEQUEN .or. IMASTER ) then
        write(lun_dynlo_tem,*) 
        flush(lun_dynlo_tem)
     end if
     !
     ! Impose temperature
     !
     if( IPARALL ) then
        nparr    =  1
        rpoin(1) =  tinf
        parre    => rpoin
        call par_broadc()
        tinf     =  parre(1)
     end if

     if( INOTMASTER ) then
        do ipoin=1,npoin
           if(kfl_fixno_tem(1,ipoin)==1) then
              bvess_tem(1,ipoin,1) = tinf
              tempe(ipoin,1)       = tinf
              tempe(ipoin,2)       = tinf
           end if
        end do
     end if

  case(2_ip)
     !
     ! After a time step: writes output data: average temperature
     !
     if( IMASTER .or. ISEQUEN ) then

        call livinf(59_ip,'OPEN  OUTPUT FILE DYNAMIC COUPLING',1_ip)
        open(unit=lun_dynou_tem,file=trim(fil_dynou_tem),status='unknown',iostat=ioerr)
        if(ioerr/=0) then
           call runend('DYNAMICAL COUPLING 1: COULD NOT OPEN OUTPUT FILE')
        end if

        if(postp(1) % npp_setsb(1)==0) &
             call runend('DYNAMICAL COUPLING 1: OUTPUT BOUNDARY SET AVERAGE TEMPERATURE')

        lchan( 1) = 305
        lchan( 2) = 289
        lchan( 3) = 291
        lchan( 4) = 249
        lchan( 5) = 283
        lchan( 6) = 274
        lchan( 7) = 243
        lchan( 8) = 268
        lchan( 9) = 263
        lchan(10) = 253
        lchan(11) = 300
        lchan(12) = 302
        lchan(13) = 276
        lchan(14) = 256
        lchan(15) = 278
        lchan(16) = 316
        tout      = 0.0_rp
        surf      = 0.0_rp

        do ibset=1,nbset
           ichan=0
           do while(ichan<16)
              ichan=ichan+1
              if(lbsec(ibset)==lchan(ichan)) ichan=17
           end do
           if(ichan==17) then
              surf=surf+postp(1)%vbset(postp(1)%nvabs+1,ibset)
              tout=tout+postp(1)%vbset(1,ibset)*postp(1)%vbset(postp(1)%nvabs+1,ibset)
           end if
        end do
        tout=tout/surf
        write(lun_dynou_tem,*) tout      
        write(lun_dynou_tem,*) '1' 
        close(lun_dynou_tem)

        write(lun_dynlo_tem,40) tout
        write(lun_dynlo_tem,*) 
        flush(lun_dynlo_tem)

        if(ipass==0) then
           ipass=1
           write(lun_dynre_tem,100)
        end if
        write(lun_dynre_tem,50) cutim,tinf,tout,vinf,vint
        flush(lun_dynre_tem)

        call livinf(59_ip,'CLOSE OUTPUT FILE DYNAMIC COUPLING',1_ip)

     end if

  end select

1 format(i1)
10 format('TIME=                                  ',e13.6)
12 format('INFLOW TEMPERATURE   (READ FROM FILE)= ',e13.6)
13 format('MEAN INFLOW VELOCITY (READ FROM FILE)= ',e13.6) 
15 format('INTERPOLATION NOT NEEDED: MEAN VELOCITY HAS NOT CHANGED') 
20 format('CLOSEST MEAN INFLOW VELOCITY=          ',e13.6)  
30 format('RESULTS INTERPOLATED FROM ANOTHER MESH:') 
31 format('MESH FILE:                             ',a)  
32 format('RESULT FILE:                           ',a)  
40 format('OUTPUT TEMPERATURE   (WRITE TO FILE)=  ',e13.6)
50 format(10(1x,e13.6))
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Temper Inflow      3. Temper. outflow  '  ,/,&
       &   '# --| 4. Velocity inflow   5. Nearest velocity   ' ,/,&
       &   '# ','          1','            2','            3','            4',&
       &   '            5')

end subroutine tem_dync01

subroutine tem_dync01_finam(vinf,winf,vout)
  use def_kintyp, only : ip,rp  
  use mod_messages, only : livinf
  use mod_messages, only : messages_live
  implicit none
  real(rp),     intent(in)  :: vinf
  character(5), intent(out) :: winf
  real(rp),     intent(out) :: vout
  integer(ip),  parameter   :: nvelo=50
  integer(ip)               :: ivelo
  real(rp)                  :: vint(2,0:nvelo),dist1,dist2
  character(5)              :: wint(0:nvelo)

  ! File name        Pressure gradient    Mean velocity
  !------------------------------------------------------------------
  wint( 0)='0';       vint(1, 1)= 0.0_rp;     vint(2, 0)= 0.0_rp 
                                                                            
  wint( 1)='005';     vint(1, 1)= 0.05_rp;    vint(2, 1)= 0.82417129E-001_rp
  wint( 2)='00625';   vint(1, 2)= 0.0625_rp;  vint(2, 2)= 0.10302141E+000_rp
  wint( 3)='0075';    vint(1, 3)= 0.075_rp;   vint(2, 3)= 0.12362570E+000_rp
  wint( 4)='01';      vint(1, 4)= 0.1_rp;     vint(2, 4)= 0.16472611E+000_rp
  wint( 5)='011';     vint(1, 5)= 0.11_rp;    vint(2, 5)= 0.18103794E+000_rp
  wint( 6)='0125';    vint(1, 6)= 0.125_rp;   vint(2, 6)= 0.20543896E+000_rp
  wint( 7)='015';     vint(1, 7)= 0.15_rp;    vint(2, 7)= 0.24581825E+000_rp
  wint( 8)='0175';    vint(1, 8)= 0.175_rp;   vint(2, 8)= 0.28573969E+000_rp
  wint( 9)='019';     vint(1, 9)= 0.19_rp;    vint(2, 9)= 0.30942981E+000_rp
  wint(10)='02';      vint(1,10)= 0.2_rp;     vint(2,10)= 0.32510479E+000_rp
  wint(11)='0225';    vint(1,11)= 0.225_rp;   vint(2,11)= 0.36382267E+000_rp
  wint(12)='025';     vint(1,12)= 0.25_rp;    vint(2,12)= 0.40182928E+000_rp
  wint(13)='03';      vint(1,13)= 0.3_rp;     vint(2,13)= 0.47548647E+000_rp
  wint(14)='0325';    vint(1,14)= 0.325_rp;   vint(2,14)= 0.51106377E+000_rp
  wint(15)='035';     vint(1,15)= 0.35_rp;    vint(2,15)= 0.54577862E+000_rp
  wint(16)='0375';    vint(1,16)= 0.375_rp;   vint(2,16)= 0.57962351E+000_rp
  wint(17)='04';      vint(1,17)= 0.4_rp;     vint(2,17)= 0.61261557E+000_rp
  wint(18)='05';      vint(1,18)= 0.5_rp;     vint(2,18)= 0.73608711E+000_rp
  wint(19)='06';      vint(1,19)= 0.6_rp;     vint(2,19)= 0.84726147E+000_rp
  wint(20)='08';      vint(1,20)= 0.8_rp;     vint(2,20)= 0.10397879E+001_rp
  wint(21)='10';      vint(1,21)= 1.0_rp;     vint(2,21)= 0.12032424E+001_rp
  wint(22)='12';      vint(1,22)= 1.2_rp;     vint(2,22)= 0.13467132E+001_rp
  wint(23)='14';      vint(1,23)= 1.4_rp;     vint(2,23)= 0.14759501E+001_rp
  wint(24)='16';      vint(1,24)= 1.6_rp;     vint(2,24)= 0.15946090E+001_rp
  wint(25)='18';      vint(1,25)= 1.8_rp;     vint(2,25)= 0.17050904E+001_rp
  wint(26)='20';      vint(1,26)= 2.0_rp;     vint(2,26)= 0.18090326E+001_rp
  wint(27)='22';      vint(1,27)= 2.2_rp;     vint(2,27)= 0.19075974E+001_rp
  wint(28)='24';      vint(1,28)= 2.4_rp;     vint(2,28)= 0.20016382E+001_rp
  wint(29)='26';      vint(1,29)= 2.6_rp;     vint(2,29)= 0.20918000E+001_rp
  wint(30)='28';      vint(1,30)= 2.8_rp;     vint(2,30)= 0.21785861E+001_rp
  wint(31)='30';      vint(1,31)= 3.0_rp;     vint(2,31)= 0.22623954E+001_rp
  wint(32)='32';      vint(1,32)= 3.2_rp;     vint(2,32)= 0.23435524E+001_rp
  wint(33)='34';      vint(1,33)= 3.4_rp;     vint(2,33)= 0.24223240E+001_rp
  wint(34)='36';      vint(1,34)= 3.6_rp;     vint(2,34)= 0.24989339E+001_rp
  wint(35)='38';      vint(1,35)= 3.8_rp;     vint(2,35)= 0.25735718E+001_rp
  wint(36)='40';      vint(1,36)= 4.0_rp;     vint(2,36)= 0.26463996E+001_rp
  wint(37)='42';      vint(1,37)= 4.2_rp;     vint(2,37)= 0.27175590E+001_rp
  wint(38)='44';      vint(1,38)= 4.4_rp;     vint(2,38)= 0.27871713E+001_rp
  wint(39)='46';      vint(1,39)= 4.6_rp;     vint(2,39)= 0.28553448E+001_rp
  wint(40)='48';      vint(1,40)= 4.8_rp;     vint(2,40)= 0.29221749E+001_rp
  wint(41)='50';      vint(1,41)= 5.0_rp;     vint(2,41)= 0.29877448E+001_rp
  wint(42)='52';      vint(1,42)= 5.2_rp;     vint(2,42)= 0.30521319E+001_rp
  wint(43)='54';      vint(1,43)= 5.4_rp;     vint(2,43)= 0.31154031E+001_rp
  wint(44)='56';      vint(1,44)= 5.6_rp;     vint(2,44)= 0.31776203E+001_rp
  wint(45)='58';      vint(1,45)= 5.8_rp;     vint(2,45)= 0.32388382E+001_rp
  wint(46)='60';      vint(1,46)= 6.0_rp;     vint(2,46)= 0.32991081E+001_rp
  wint(47)='62';      vint(1,47)= 6.2_rp;     vint(2,47)= 0.33584756E+001_rp
  wint(48)='64';      vint(1,48)= 6.4_rp;     vint(2,48)= 0.34169832E+001_rp 
  wint(49)='66';      vint(1,49)= 6.6_rp;     vint(2,49)= 0.34746692E+001_rp
  wint(50)='68';      vint(1,50)= 6.8_rp;     vint(2,50)= 0.35315695E+001_rp

  ivelo=-1
  do while(ivelo<nvelo)
     ivelo=ivelo+1
     if(vinf<=vint(2,ivelo)) then
        dist1=abs(vinf-vint(2,ivelo-1))
        dist2=abs(vinf-vint(2,ivelo))
        if(dist1<dist2) then
           winf=wint(ivelo-1)
           vout=vint(2,ivelo-1)
           return
        else
           winf=wint(ivelo)
           vout=vint(2,ivelo)
           return
        end if
     end if
  end do
  call messages_live('TEM_DYNCO01: MASS FLOW RATE IS TOO HIGH')
  winf=wint(nvelo)
  vout=vint(2,nvelo)

end subroutine tem_dync01_finam
