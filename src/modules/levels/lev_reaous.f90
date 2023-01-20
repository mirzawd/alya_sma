!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_reaous

  !-----------------------------------------------------------------------
  !
  ! This routine reads the output strategy for the level set advection
  ! equation.
  !
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_levels
  use      def_domain
  use mod_ecoute, only : ecoute
  use mod_output_postprocess, only : output_postprocess_read
  implicit none
  integer(ip) :: ngaug,idime

  if( INOTSLAVE ) then
     !
     ! Initializations.
     !
     npp_gauge_lev = 0
     npp_nbgau_lev = 0
     cogau_lev     = 0.0_rp
     ngaug         = 0
     npp_inter_lev = 0
     !
     ! Reach the section.
     !
     call ecoute('lev_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('lev_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('lev_reaous')

        call output_postprocess_read()

        if(words(1)=='INTER') then
           npp_inter_lev  = getint('STEPS',1_ip,'Frequency to output interface mesh')
        else if(words(1)=='GAUGE') then
           npp_gauge_lev = 1
           if(words(2)=='NUMBE')&
                npp_nbgau_lev=getint('NUMBE',0_ip,' number of interface gauge to postprocess')
           call ecoute('lev_reaous')
           if(words(1)/='ENDGA') then
              do while(words(1)/='ENDGA')
                 if(ngaug>=ngaug_lev) then
                    call runend('LEV_REAOUS: TOO MANY GAUGES ')
                 else
                    ngaug=ngaug+1
                    typga_lev(ngaug) = int(param(1))
                    do idime=2,ndime+1
                       cogau_lev(idime-1,ngaug)=param(idime)
                    end do
                 end if
                 call ecoute('lev_reaous')
              enddo
           endif


        end if
     end do
  end if

end subroutine lev_reaous
    
