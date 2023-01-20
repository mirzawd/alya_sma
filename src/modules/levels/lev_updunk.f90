!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Wavequ/lev_updunk
  ! NAME 
  !    lev_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates for the level set
  !    function
  ! USED BY
  !    lev_begste (itask=1)
  !    lev_begite (itask=2)
  !    lev_endite (itask=3, inner loop) 
  !    lev_endite (itask=4, outer loop) 
  !    lev_endste (itask=5)
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_levels
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,ipoin
  integer(ip)             :: iredi

  select case (itask)

  case(1_ip)
     !
     ! Assign u(n,0,*) <-- u(n-1,*,*), initial guess for outer iterations
     !
     if( INOTMASTER ) then
        if(kfl_conbc_lev==0)then
           do ipoin=1,npoin
              if(kfl_funno_lev(ipoin)==0)fleve(ipoin,2) = fleve(ipoin,min(3_ip,ncomp_lev))
           end do

        else
           do ipoin=1,npoin
              fleve(ipoin,2) = fleve(ipoin,min(3_ip,ncomp_lev))
           end do
        endif
     end if
  case(2_ip)
     !
     ! Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
     !
     if( INOTMASTER ) then
        if(kfl_conbc_lev==0)then
           do ipoin=1,npoin
              if(kfl_funno_lev(ipoin)==0)fleve(ipoin,1) = fleve(ipoin,2)
           end do
           do ipoin=1,npoin
              if(kfl_funno_lev(ipoin)==0)unkno(ipoin) = fleve(ipoin,2)
           end do
        else
           do ipoin=1,npoin
              fleve(ipoin,1) = fleve(ipoin,2)
           end do
           do ipoin=1,npoin
              unkno(ipoin) = fleve(ipoin,2)
           end do

        endif
     end if

  case(3_ip)
     !
     ! Assign u(n,i,j-1) <-- u(n,i,j), update of the wave amplitude
     !
     if( INOTMASTER ) then
        if(kfl_conbc_lev==0)then
           do ipoin=1,npoin
              if(kfl_funno_lev(ipoin)==0)fleve(ipoin,1) = unkno(ipoin)
           end do
        else
           do ipoin=1,npoin
              fleve(ipoin,1) = unkno(ipoin)
           end do

        endif
     end if

  case(4_ip)
     !
     ! Assign u(n,i-1,*) <-- u(n,i,*)
     !        
     if( INOTMASTER ) then
        do ipoin=1,npoin
           fleve(ipoin,2) = fleve(ipoin,1)
        end do
     endif

  case(5_ip)
     !
     ! Update previous time steps
     !         
     if(kfl_timet_lev==1) then ! Explicit treatment case 

        iredi=0
        if(mod(ittim,nfred_lev)==0.and.iredi==0_ip) then
           call lev_redist
           if(tyred_lev==1) then
              call lev_redist()
           else if( tyred_lev >1 ) then
              call lev_redieq()
           end if
           iredi=1
        endif


        if( INOTMASTER ) then
           do ipoin=1,npoin
              fleve(ipoin,4)=fleve(ipoin,3) ! u^{n-1}
              fleve(ipoin,3)=fleve(ipoin,1) ! u^{n}
           end do
        endif

     else if(kfl_timet_lev==2) then ! Implicit treatment case 

        if(kfl_tiacc_lev==1) then

           if(kfl_corvo_lev==1) then
              call lev_corvol()
           else if(kfl_corvo_lev==2) then
              call lev_corvo2()
           else
              call lev_calvol()
           endif

           if( INOTSLAVE ) then   
              write(lun_volum_lev,*) cutim,' ',volit_lev
              flush(lun_volum_lev)
           endif

           iredi=0

           if(mod(ittim,nfred_lev)==0.and.iredi==0_ip) then
              if(tyred_lev==1) then
                 call lev_redist()
              else if(tyred_lev==-1) then
                 call lev_redist_generalized_distance()
              else if(tyred_lev==-2) then
                 call lev_redist_geometrical_distance()
              else if( tyred_lev >1 ) then
                 call lev_redieq()
              end if
              iredi=1
           endif


           if( INOTMASTER ) then
              do ipoin=1,npoin
                 if((kfl_tiaor_lev==2).and.(kfl_tisch_lev==2)) then 
                    fleve(ipoin,4)=fleve(ipoin,3) ! u^{n-2}
                 end if
                 fleve(ipoin,3)=fleve(ipoin,1) ! u^{n-1}
              end do
           endif

        else if(kfl_tiacc_lev==2) then 

           if(kfl_tisch_lev==1) then
              !
              ! Crank-Nicholson method 
              !        
              if( INOTMASTER ) then
                 do ipoin=1,npoin
                    if(kfl_funno_lev(ipoin)==0)fleve(ipoin,1) = 2.0_rp*fleve(ipoin,1)-fleve(ipoin,3)
                 end do
              endif

              if(kfl_corvo_lev==1) then
                 call lev_corvol()
              else if(kfl_corvo_lev==2) then
                 call lev_corvo2()
              else
                 call lev_calvol()
              endif

              if( INOTSLAVE ) then   
                 write(lun_volum_lev,*) cutim,' ',volit_lev
                 flush(lun_volum_lev)
                 !print*,' vol it ',volit_lev
              endif

              iredi=0
              if(mod(ittim,nfred_lev)==0.and.iredi==0_ip) then
                 if(tyred_lev==1) then
                    call lev_redist()
                 else if(tyred_lev==-1) then
                    call lev_redist_generalized_distance()
                 else if(tyred_lev==-2) then
                    call lev_redist_geometrical_distance()
                 else if( tyred_lev >1 ) then
                    call lev_redieq()
                 end if
                 iredi=1
              endif


              if( INOTMASTER ) then
                 do ipoin=1,npoin
                    fleve(ipoin,3)=fleve(ipoin,1) ! u^{n-1}
                 end do
              endif

           else if(kfl_tisch_lev==2) then 

              !
              ! BDF2 scheme
              !

              if(kfl_corvo_lev==1) then
                 call lev_corvol()
              else if(kfl_corvo_lev==2) then
                 call lev_corvo2()
              else
                 call lev_calvol()
              endif

              if( INOTSLAVE ) then   
                 write(lun_volum_lev,*) cutim,' ',volit_lev
                 flush(lun_volum_lev)
              endif

              iredi=0
              if(mod(ittim,nfred_lev)==0.and.iredi==0_ip) then
                 if(tyred_lev==1) then
                    call lev_redist()
                 else if(tyred_lev==-1) then
                    call lev_redist_generalized_distance()
                 else if(tyred_lev==-2) then
                    call lev_redist_geometrical_distance()
                 else if( tyred_lev >1 ) then
                    call lev_redieq()
                 end if
                 iredi=1
              endif


              if( INOTMASTER ) then
                 do ipoin=1,npoin
                    fleve(ipoin,4)=fleve(ipoin,3) ! u^{n-2}
                    fleve(ipoin,3)=fleve(ipoin,1) ! u^{n-1}
                 end do
              end if

           endif


        end if


     end if


  case(6_ip) 
     !
     ! Assign u(n,1,*) <-- u(n-1,*,*), when reading from restart file
     ! 
     icomp=min(3_ip,ncomp_lev) 
     do ipoin=1,npoin
        if(kfl_funno_lev(ipoin)==0)fleve(ipoin,1) = fleve(ipoin,icomp)
     end do

  case(7_ip) 
     !
     ! UNKNO <= FLEVE
     ! 
     if( INOTMASTER ) then
        do ipoin=1,npoin
           unkno(ipoin) = fleve(ipoin,1) 
        end do
     end if

  case(11_ip)
     !
     ! Actualize initial condition after redistantiation
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           fleve(ipoin,3)         = fleve(ipoin,1)
           fleve(ipoin,2)         = fleve(ipoin,1)
           fleve(ipoin,ncomp_lev) = fleve(ipoin,1)
           flsex_lev(ipoin)       = fleve(ipoin,1)
           if(kfl_timet_lev==1) then
              !
              ! Explicit 
              !
              fleve(ipoin,4) = fleve(ipoin,1)
           end if
        end do
     end if

  end select

end subroutine lev_updunk

