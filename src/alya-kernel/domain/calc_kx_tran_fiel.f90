!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine calc_kx_tran_fiel()
  !-----------------------------------------------------------------------
  !****f* Domain/calc_kx_tran_fiel
  ! NAME
  !    calc_kx_tran_fiel
  ! DESCRIPTION
  !    Obtains k_tran_fiel and x_tran_fiel for field with more than 1 step
  !    k_tran_fiel(ifiel) indicates to which interval the current time belongs.
  !    x_tran_fiel(ifiel) indicates the position between the begining and end of the interval.
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_mpio,            only : mpio_memor
  use mod_mpio_files,      only : fil_field,xfiel_name
  use mod_mpio_generic_io
  use mod_opfpos,          only : postpr_intto8
  use mod_messages,        only : livinf
  use mod_memory,          only : memory_deallo


  implicit none

  integer(ip)        :: ifiel,kk,istep
  real(rp)           :: mintime
  real(rp)           :: maxtime
  real(rp)           :: cortime                  ! corrected time, within one period
  real(rp), pointer  :: buf(:,:)                 ! temporary storage for the fields to copy to xfiel().a
  integer(ip)        :: kdime, size_p
  real(rp)           :: delay, new_period


  nullify(buf)

  do ifiel = 1,nfiel

     if ( kfl_field(4,ifiel) > 1 ) then !number of timesteps

        kexist_tran_fiel = 1
        mintime        = time_field(ifiel) % a(1)
        maxtime        = time_field(ifiel) % a(kfl_field(4,ifiel))
        !
        ! Periodic field, special treament
        !
        if( kfl_field(7,ifiel) == 1 ) then
           cortime = mod(cutim,maxtime)
           if( kfl_field(3,ifiel)>0 .and. cutim > maxtime )then
               delay   = time_field(ifiel) % a(kfl_field(3,ifiel)) - mintime
               new_period = maxtime - time_field(ifiel) % a(kfl_field(3,ifiel)) 
               cortime = mod(cutim - delay, new_period) + delay 
           end if
        else
           cortime = cutim
        end if

        if(      cortime <= mintime ) then
           !
           ! Below minimum time
           !
           if ( kfl_field(3,ifiel) > 1 ) then 
               if( cutim >= time_field(ifiel) % a(kfl_field(3,ifiel)))then
                  k_tran_fiel(ifiel) = kfl_field(3,ifiel)
               else
                k_tran_fiel(ifiel) = 1
               end if
           else
               k_tran_fiel(ifiel) = 1
           end if
           x_tran_fiel(ifiel) = 1.0_rp
        else if( cortime >= maxtime ) then
           !
           ! Above maximum time
           !
           k_tran_fiel(ifiel) = kfl_field(4,ifiel)-1 
           x_tran_fiel(ifiel) = 0.0_rp
        else
           !
           ! Right in the time range
           !
           kk = kfl_field(4,ifiel)
            steps: do istep = 1,kfl_field(4,ifiel)
              if( cortime <= time_field(ifiel) % a(istep) ) then
                 kk = istep - 1
                 exit steps 
              end if
            end do steps
            kk = max(kk,1_ip)
            kk = min(kk,kfl_field(4,ifiel)-1)
            k_tran_fiel(ifiel) = kk
            x_tran_fiel(ifiel) = &
                ( time_field(ifiel) % a(kk+1) - cortime                   ) / &
                ( time_field(ifiel) % a(kk+1) - time_field(ifiel) % a(kk) )
          
        end if


        ! ON DEMAND field loading. 
        ! Check if field ifiel are available for real kk and kk+1
        ! We need:
        ! - xfiel(ifiel) % a(:,:,kk)
        ! - xfiel(ifiel) % a(:,:,kk+1)
        ! Save them to 
        ! - xfiel(ifiel) % a(:,:,1)
        ! - xfiel(ifiel) % a(:,:,2)
        ! reset kk to 1

        if (kfl_field(6,ifiel) == 1 ) then !on demand
           !print *, "0 = ", 0_rp
           !print *, "1 = ", 1_rp
           !print *, "cutim = ",cutim
           !print *, "k_tran_fiel_real(ifiel) = ",k_tran_fiel_real(ifiel)
           !print *, "k_tran_fiel(ifiel) ) =",k_tran_fiel(ifiel) 
           if ( k_tran_fiel_real(ifiel) /= k_tran_fiel(ifiel) )  then
              !this runs only if the kk changed from the previous time this block was executed
              !save the real k_tran_field
              !set the k_tran_field to 1 because we will read only fields 1 and 2  

              if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                 size_p=nelem
              else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
                 size_p=npoin
              else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
                 size_p=nboun
              else
                 call runend('UNDEFINED FIELD TYPE')
              end if
              kdime = kfl_field(6,ifiel)

              !print *,'old k_tran_fiel_real(ifiel) =',k_tran_fiel_real(ifiel)

              k_tran_fiel_real(ifiel) = k_tran_fiel(ifiel)
              k_tran_fiel(ifiel) = 1
              !print *,'k_tran_fiel_real(ifiel) =',k_tran_fiel_real(ifiel)
              !print *,'k_tran_fiel(ifiel) =',k_tran_fiel(ifiel)

              !
              ! Read fields kk and kk+1 and put them at timesteps 1 and 2
              !
              ! Field kk
              !print *,'Number of points: ', npoin

              if(INOTMASTER) then
                 fil_field = xfiel_name(ifiel, k_tran_fiel_real(ifiel))
                 call livinf(0_ip,'READING FIELD 1 FROM '//fil_field,0_ip)

                 !                print *,'before: size_p= ', size_p, ', kdime=',kdime

                 !if ( kdime == 1 ) then
                 ! This fails for some reason
                 !call MPIO_READ(buf, fil_field, size_p )
                 !else
                 call MPIO_READ(buf, fil_field, size_p, kdime )
                 !end if            
                 !print *,'buf size:',size(buf,1,8),'x',size(buf,2,8)

                 !                print *,'after: size_p= ', size_p, ', kdime=',kdime
                 !                print *,'buf=', buf
                 !                do istep=1,size(buf, 2, ip)
                 !                    print *, buf(:,istep)
                 !                end do


                 !print *,"Assigning buf to k_tran_fiel(ifiel) = ",k_tran_fiel(ifiel)
                 xfiel(ifiel) % a(:,:,k_tran_fiel(ifiel)) = buf
                 call memory_deallo(mpio_memor,'buf', "calc_kx_tran_fiel_MPIO_READ", buf)
                 nullify(buf)


                 !                  
                 ! Field kk+1
                 !
                 fil_field = xfiel_name(ifiel, k_tran_fiel_real(ifiel)+1)
                 call livinf(0_ip,'READING FIELD 2 FROM '//fil_field,0_ip)
                 !if ( kdime == 1 ) then
                 !this fails for some reason
                 !call MPIO_READ( buf, fil_field,  size_p )
                 !else    
                 call MPIO_READ( buf, fil_field,  size_p, kdime )
                 !end if


                 !print *,"Assigning buf to k_tran_fiel(ifiel)+1 = ",k_tran_fiel(ifiel)+1
                 xfiel(ifiel) % a(:,:,k_tran_fiel(ifiel)+1) = buf
                 call memory_deallo(mpio_memor,'buf', "calc_kx_tran_fiel_MPIO_READ", buf)
                 nullify(buf)

                 !print *,'calc_kx_tran_fiel'
                 !print *,'-------------------------------------------------------------'
                 !print *,'T=',cutim
                 !print *,'read veloc 1='
                 !do istep=1,size(xfiel(ifiel) % a, 2, ip)
                 !    print *, xfiel(ifiel) % a(:,istep,1)
                 !end do
                 !print *,'read veloc 2='
                 !do istep=1,size(xfiel(ifiel) % a, 2, ip)
                 !    print *, xfiel(ifiel) % a(:,istep,2)
                 !end do
                 !print *,'-------------------------------------------------------------'

              end if !inotmaster
           else
              !timestep did not change
              k_tran_fiel(ifiel) = 1
           end if !if timestep has changed

        else !if ondemand fields required

           !assign the value in any case, in case anyone decides to use this variable
           k_tran_fiel_real(ifiel) = k_tran_fiel(ifiel)

        end if

     end if ! kfl_field(4,ifiel) > 1
  end do ! ifiel = 1,nfiel


end subroutine calc_kx_tran_fiel
