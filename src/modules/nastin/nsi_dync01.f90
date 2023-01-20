!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_dync01(itask)
  !------------------------------------------------------------------------
  !****f* Temper/nsi_dync01
  ! NAME 
  !    nsi_dync01
  ! DESCRIPTION
  !    Marek: Cooling system + typical room
  ! USES
  ! USED BY
  !    nsi_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,icode,ioerr
  real(rp),    target     :: rpoin(2)
  real(rp),    save       :: vinf

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
              open(unit=lun_dynin_nsi,file=trim(fil_dynin_nsi),status='old',iostat=ioerr)
           end do
           read(lun_dynin_nsi,*,err=2,end=2) vinf
           read(lun_dynin_nsi,*,err=2,end=2) icode
           if(icode==1) exit icode0
           icode=0
           ioerr=1
2          close(lun_dynin_nsi)
        end do icode0
        close(lun_dynin_nsi,status='delete')
        call livinf(59_ip,'CLOSE INPUT FILE DYNAMIC COUPLING',1_ip)

        !write(lun_dynlo_nsi,'(a)') '-------------------------------------------' 
        !write(lun_dynlo_nsi,*) 
        !call livinf(59_ip,'CLOSE INPUT FILE DYNAMIC COUPLING',1_ip)
        !write(lun_dynlo_nsi,10) cutim
        !write(lun_dynlo_nsi,13) vinf 

     end if
     if( IPARALL ) then
        nparr    =  1
        rpoin(1) =  vinf
        parre    => rpoin
        call par_broadc()
        vinf     =  parre(1)
     end if

     if( ISEQUEN .or. IMASTER ) then
        !write(lun_dynlo_nsi,*) 
        !flush(lun_dynlo_nsi)
     end if

     if( INOTMASTER ) then
        do ipoin=1,npoin
           if(kfl_fixno_nsi(1,ipoin)==8) then
              bvess_nsi(1,ipoin,1) = vinf
              veloc(1,ipoin,1)     = vinf
              veloc(1,ipoin,2)     = vinf
           end if
        end do
     end if

  end select

1 format(i1)
10 format('TIME=                                  ',e12.6)
12 format('INFLOW TEMPERATURE   (READ FROM FILE)= ',e12.6)
13 format('MEAN INFLOW VELOCITY (READ FROM FILE)= ',e12.6) 
15 format('INTERPOLATION NOT NEEDED: MEAN VELOCITY HAS NOT CHANGED') 
20 format('CLOSEST MEAN INFLOW VELOCITY=          ',e12.6)  
30 format('RESULTS INTERPOLATED FROM ANOTHER MESH:') 
31 format('MESH FILE:                             ',a)  
32 format('RESULT FILE:                           ',a)  
40 format('OUTPUT TEMPERATURE   (WRITE TO FILE)=  ',e12.6)
50 format(10(1x,e12.6))
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Temper Inflow      3. Temper. outflow  '  ,/,&
       &   '# --| 4. Velocity inflow   5. Nearest velocity   ' ,//,&
       &   '$ ','          1','            2','            3','            4',&
       &   '            5')

end subroutine nsi_dync01
