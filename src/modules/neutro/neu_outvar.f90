!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_outvar.f90
!> @date    01/04/2016
!> @author  Guillaume Houzeaux
!> @brief   Postprocess
!> @details Postprocess
!> @} 
!-----------------------------------------------------------------------

subroutine neu_outvar(ivari,iter)

  use def_master
  use def_domain
  use def_neutro
  use mod_postpr
  use mod_outvar
  use def_parame, only     :  in4pi
  implicit none
  integer(ip), intent(in) :: ivari,iter
  integer(ip)             :: idire,idime,ipoin,iener, ielem,iboun,ipoin2,pnodb
  integer(ip)             :: unit_aux=5678
  character(5)            :: wopos(3)
  character(20)           :: wdire,zener
  character(30)    :: file_aux
  character(6)    :: ini='ITERA'
  character(4)            :: fin='.dat'
!   character(3)            :: it='  '
  
  character(1)            :: coma=','
  real(8)             :: modul_mio
!   real(8)    ::   facG=1.0

  external :: memgen

  ! select case (arrays_name(ivari))
  select case ( ivari )  

  case ( 0_ip )
      return
!   case ( 10_ip )
!      !
!      ! output in iteration N
!      !
!        if( INOTMASTER ) then
!           write(file_aux,'(i4)') iter

!          file_aux=ini//TRIM(file_aux)//fin
!          open(unit=unit_aux,file=file_aux)
         
!          if(ndime==2) then
!              write(unit_aux,*) 'X , Y, Cx, Cy, ModC '
!          elseif(ndime==3) then
!              write(unit_aux,*) 'X , Y, Z, Cx, Cy, Cz, ModC '
!          endif 


!         call memgen(0_ip,ndime,npoin)!> Allocate memory for a vector array of length ndime and element dimension npoin which 
!                                      !>  will be allocated in gevec (for a scalar it would be gesca)
!       !   gevec(1:ndime,1:npoin)=0.0_rp
!         do ipoin = 1,npoin  !> Iterating over the points of the mesh
!            do iener = 1,num_energies_neu !> Iterating over the energies in each point
!             do idire = 1,num_directions_neu !> Iterating over the directions in each point
!               do idime = 1,ndime !> Iterating over the spatial dimensions
!                  !> The current is the vector sum of the flux*direction*weight
!                  gevec(idime,ipoin) = gevec(idime,ipoin) + neutr(iener,idire,ipoin,1) * &
!                                           weigd_neu(idire) * direc_neu(idime,idire) !*ener_weigd_neu(iener)
!                !   gevec(idime,ipoin) = gevec(idime,ipoin) + neutr(iener,idire,ipoin,1) * &
!                !                            weigd_neu(idire) * direc_neu(idime,idire)*in4pi !*ener_weigd_neu(iener)
!               end do
!             end do
!            end do
!         end do
     
!         do ipoin = 1,npoin 
!              modul_mio = 0.0_rp
!              do idime=1,ndime
!                 modul_mio = modul_mio  + gevec(idime,ipoin)*gevec(idime,ipoin)
!               enddo
!               modul_mio = sqrt(modul_mio)
!             if(ndime==2) then
!                write(unit_aux,'(E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5)')  &
!                   coord(1,ipoin),coma,coord(2,ipoin),coma,gevec(1,ipoin),coma,gevec(2,ipoin),coma,modul_mio
!             elseif(ndime==3) then
!                 write(unit_aux,'(E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5)')  &
!                   coord(1,ipoin),coma,coord(2,ipoin),coma,coord(3,ipoin),coma,gevec(1,ipoin),coma,&
!                   gevec(2,ipoin),coma,gevec(3,ipoin),coma,modul_mio
!             endif
           
!         enddo 
!         close(unit_aux)
!      end if


!   case ( 'CURRE' )
  case ( 1_ip )
     !
     ! CURRE: neutron current
     !

   
  ! current 

     if( INOTMASTER ) then

         ! open(unit=unit_curr_neu,file=fil_outCUR_neu)
         
         ! if(ndime==2) then
         !     write(unit_curr_neu,*) 'X , Y, Cx, Cy, ModC '
         ! elseif(ndime==3) then
         !     write(unit_curr_neu,*) 'X , Y, Z, Cx, Cy, Cz, ModC '
         ! endif 


        call memgen(0_ip,ndime,npoin)!> Allocate memory for a vector array of length ndime and element dimension npoin which 
                                     !>  will be allocated in gevec (for a scalar it would be gesca)
      !   gevec(1:ndime,1:npoin)=0.0_rp
        do ipoin = 1,npoin  !> Iterating over the points of the mesh
           do iener = 1,num_energies_neu !> Iterating over the energies in each point
             do idire = 1,num_directions_neu !> Iterating over the directions in each point
               do idime = 1,ndime !> Iterating over the spatial dimensions
                  !> The current is the vector sum of the flux*direction*weight
                  gevec(idime,ipoin) = gevec(idime,ipoin) + neutr(iener,idire,ipoin,1) * &
                                          direc_neu(idime,idire) * weigd_neu(idire) !* ener_weigd_neu(iener)
                  ! gevec(idime,ipoin) = gevec(idime,ipoin) + neutr(iener,idire,ipoin,1) * &
                  !                         direc_neu(idime,idire) * weigd_neu(idire)*in4pi !* ener_weigd_neu(iener) 
               end do
             end do
           end do
        end do
     
      !   do ipoin = 1,npoin 
      !        modul_mio =0.0
      !        do idime=1,ndime
      !           modul_mio = modul_mio  + gevec(idime,ipoin)*gevec(idime,ipoin)
      !        enddo
      !        modul_mio = sqrt(modul_mio)
      !        if(ndime==2) then
      !             write(unit_curr_neu,'(E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5)')  coord(1,ipoin),coma,coord(2,ipoin),coma,&
      !                                                                       gevec(1,ipoin),coma,gevec(2,ipoin),coma,modul_mio
      !        elseif(ndime==3) then
      !             write(unit_curr_neu,'(E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5,a1,E15.5)')  coord(1,ipoin),coma,&
      !                                                                       coord(2,ipoin),coma,coord(3,ipoin),coma,gevec(1,ipoin),&
      !                                                                       coma,gevec(2,ipoin),coma,gevec(3,ipoin),coma,modul_mio
      !        endif
      !   enddo 
      !   close(unit_curr_neu)
     end if

!   case ( 'FLUX' )
  case ( 2_ip )
     !
     ! FLUX: neutron flux
     !
     if( INOTMASTER ) then
        
      !   open(unit=unit_flux_neu,file=fil_outFLX_neu)
      !   if(ndime==2) then
      !       write(unit_flux_neu,*) 'X , Y, Flux '
      !   elseif(ndime==3) then
      !       write(unit_flux_neu,*) 'X , Y, Z, Flux '
      !   endif
        call memgen(0_ip,npoin,0_ip) !> Allocate memory for an npoin scalar array in gesca
      !   gesca(1:npoin)=0.0_rp
        do ipoin = 1,npoin  
          do iener = 1,num_energies_neu !> Iterating over the energies in each point
            do idire = 1,num_directions_neu !> Iterating over the directions in each point
               !> We add the contributions for each direction (scalarly) multiplying by the weight of that direction
               gesca(ipoin) = gesca(ipoin) + neutr(iener,idire,ipoin,1) * weigd_neu(idire) !*ener_weigd_neu(iener)
               ! gesca(ipoin) = gesca(ipoin) + neutr(iener,idire,ipoin,1) * weigd_neu(idire)*in4pi !*ener_weigd_neu(iener)
            end do
          enddo 
        end do

      !    do ipoin = 1,npoin 
      !       if(ndime==2) then
      !          write(unit_flux_neu,*)  coord(1,ipoin),coma,coord(2,ipoin),coma,gesca(ipoin)
      !       elseif(ndime==3) then
      !          write(unit_flux_neu,'(E15.5,a1,E15.5,a1,E15.5,a1,E15.5)')  &
      !              coord(1,ipoin),coma,coord(2,ipoin),coma,coord(3,ipoin),coma,gesca(ipoin)
      !       endif
           
      !   enddo 
      !   close(unit_flux_neu)

     end if

!   case ( 'RADIA' )
  case ( 3_ip )
     !
     ! NEUTRONS: Neutrons (what is this??) - per energy group (angle integrated)
     !
         ! Mis salidas para chequear en Studio
!  open(unit=unit_radd_neu,file=fil_outRAD_neu)


     wopos(2:3) = postp(1) % wopos(2:3,ivari)

     if( INOTMASTER ) call memgen(0_ip,npoin,0_ip)
     do iener = 1,num_energies_neu !> Iterating over the energies in each point
    ! do idire = 1,num_directions_neu
       
        zener =adjustl(intost(iener))
    !    zener =adjustl(intost(idire))
        if( iener < 10 ) then
           wopos(1) = postp(1) % wopos(1,ivari)(1:2) // '00' // trim(zener)
        elseif( iener >= 10 .and. iener < 100 ) then
          wopos(1) = postp(1) % wopos(1,ivari)(1:2)// '0' // trim(zener)
        elseif( iener >= 100 .and. iener < 1000 ) then
           wopos(1) = postp(1) % wopos(1,ivari)(1:2)// trim(zener)
        end if
!        if( idire < 10 ) then
!           wopos(1) = postp(1) % wopos(1,3)(1:2) // '00' // trim(zener)
!        elseif( idire >= 10 .and. idire < 100 ) then
!           wopos(1) = postp(1) % wopos(1,3)(1:2)// '0' // trim(zener)
!        elseif( idire >= 100 .and. idire < 1000 ) then
!          wopos(1) = postp(1) % wopos(1,3)(1:2)// trim(zener)
!        end if
        if( INOTMASTER ) then
           gesca(1:npoin)=0.0_rp
           do idire = 1,num_directions_neu
              !        do iener = 1,num_energies_neu !> Iterating over the energies in each point
              gesca(1:npoin) = gesca(1:npoin) + neutr(iener,idire,1:npoin,1)* weigd_neu(idire)
            !   gesca(1:npoin) = gesca(1:npoin) + neutr(iener,idire,1:npoin,1)* weigd_neu(idire)*in4pi
           enddo
        end if

        call postpr(gesca,wopos,ittim,cutim)
     
     end do


!     do iener = 1,num_energies_neu !> Iterating over the energies in each point
!       do idire = 1,num_directions_neu
!         
!        wdire = adjustl(intost(idire))
!        zener =adjustl(intost(iener))
!        if( idire < 10 ) then
!           wopos(1) = postp(1) % wopos(1,3)(1:3) // '0' // trim(wdire)
!        else
!           wopos(1) = postp(1) % wopos(1,3)(1:3) // trim(wdire)
!        end if
!        if( iener < 10 ) then
!           wopos(1) = wopos(1) // '0' // trim(zener)
!        else
!           wopos(1) = wopos(1) // trim(zener)
!        end if
!
!        if( INOTMASTER ) gesca(1:npoin) = neutr(iener,idire,1:npoin,1)
!        call postpr(gesca,wopos,ittim,cutim)
!       enddo
!     end do
!

     if( INOTMASTER ) call memgen(2_ip,npoin,0_ip)
     return

   ! case ( 'DIREC' )
   case ( 4_ip )
      !
      ! NEUTRONS: Neutrons (what is this??) - per angle (energy integrated)
      !
      wopos(2:3) = postp(1) % wopos(2:3,ivari)
 
      ! if( INOTMASTER ) then
          if( INOTMASTER ) call memgen(0_ip,npoin,0_ip)
      
          do idire = 1,num_directions_neu !> Iterating over the directions in each point
        
            wdire =adjustl(intost(idire)) 
            ! zener =adjustl(intost(iener))
      !    zener =adjustl(intost(idire))
            if( idire < 10 ) then
               wopos(1) = postp(1) % wopos(1,ivari)(1:2) // '00' // trim(wdire)
            elseif( idire >= 10 .and. idire < 100 ) then
               wopos(1) = postp(1) % wopos(1,ivari)(1:2)// '0' // trim(wdire)
            elseif( idire >= 100 .and. idire < 1000 ) then
               wopos(1) = postp(1) % wopos(1,ivari)(1:2)// trim(wdire)
            end if

            if( INOTMASTER ) then
               gesca(1:npoin)=0.0_rp
               do iener = 1,num_energies_neu !> Iterating over the energies in each point
               !   gesca(1:npoin) = gesca(1:npoin) + neutr(iener,idire,1:npoin,1)* weigd_neu(idire)
                  gesca(1:npoin) = gesca(1:npoin) + neutr(iener,idire,1:npoin,1) ! * weigd_neu(idire)*in4pi
               enddo
            end if
   
            call postpr(gesca,wopos,ittim,cutim)
         
         end do
         if( INOTMASTER ) call memgen(2_ip,npoin,0_ip)
      ! end if
      ! if( INOTMASTER ) call memgen(2_ip,npoin,0_ip)
      return

   ! case ( 'MATER' )
   case (5_ip)
      !DEFINE ELEMENT MATERIAL TO VISUALIZE IN OUTPUT
      if (INOTMASTER) then
         ! call memgen(1_ip,nelem,0_ip) !> Allocate memory for an nelem integer scalar array in gisca
         call memgen(0_ip,nelem,0_ip) !> Allocate memory for an nelem integer scalar array in gesca
         ! gesca(1:nelem)=0.0_rp
         do ielem=1, nelem
            ! gisca(ielem) = lmate(ielem)
            gesca(ielem) = real(lmate(ielem),rp)
         end do
      end if

   ! case ( 'HEAT' )
   case (6_ip)
      
      ! HEAT: flux * kerma (macroscopic)

      if( INOTMASTER ) then

         call memgen(0_ip,npoin,0_ip) !> Allocate memory for an npoin scalar array in gesca
         ! gesca(1:npoin)=0.0_rp

         if (associated(heat_sink)) then ! EG: is this neccesary? In NEUTRO, heat_sink should always be associated
            do ipoin = 1,npoin  
               gesca(ipoin) = heat_sink(ipoin)
            end do
         else
            do ipoin=1,npoin
               gesca(ipoin) = 0.0_rp 
            end do
         end if

      end if

   ! case ( 'SOURC' )
   case (7_ip)
      if( INOTMASTER ) then
         call memgen(0_ip,npoin,0_ip) !> Allocate memory for an npoin scalar array in gesca
         do iener = 1,num_energies_neu
            do idire = 1,num_directions_neu
               do iboun = 1,nboun
                  if ( kfl_fixbo_neu(iener,idire) % l(iboun) == 4 ) then
                     pnodb = nnode(ltypb(iboun))
                     do ipoin2=1,pnodb
                        ipoin=lnodb(ipoin2,iboun)
                        gesca(ipoin) = real(kfl_funbo_neu(current_energy_neu,current_direction_neu) % l(iboun),rp)
                     end do
                  endif
               enddo
            end do
         end do
      end if

   ! case ( 'VACUU' )
   case (8_ip)
      if( INOTMASTER ) then
         call memgen(0_ip,npoin,0_ip) !> Allocate memory for an npoin scalar array in gesca
         do iener = 1,num_energies_neu
            do idire = 1,num_directions_neu
               do iboun = 1,nboun
                  if ( kfl_fixbo_neu(iener,idire) % l(iboun) == 1 ) then
                     pnodb = nnode(ltypb(iboun))
                     do ipoin2=1,pnodb
                        ipoin=lnodb(ipoin2,iboun)
                        gesca(ipoin) = 1.0_rp
                     end do
                  endif
               enddo
            end do
         end do
      end if

   ! case ( 'REFLE' )
   case (9_ip)
      if( INOTMASTER ) then
         call memgen(0_ip,npoin,0_ip) !> Allocate memory for an npoin scalar array in gesca
         do iener = 1,num_energies_neu
            do idire = 1,num_directions_neu
               do iboun = 1,nboun
                  if ( kfl_fixbo_neu(iener,idire) % l(iboun) == 2 .OR. &
                       kfl_fixbo_neu(iener,idire) % l(iboun) == 3) then
                     pnodb = nnode(ltypb(iboun))
                     do ipoin2=1,pnodb
                        ipoin=lnodb(ipoin2,iboun)
                        gesca(ipoin) = real(kfl_fixbo_neu(iener,idire) % l(iboun),rp)
                     end do
                  endif
               enddo
            end do
         end do
      end if

   case ( 10_ip )
      !
      ! HEAT BY ENERGY GROUP
      !
      
      wopos(2:3) = postp(1) % wopos(2:3,ivari)
 
      if( INOTMASTER ) call memgen(0_ip,npoin,0_ip)
      do iener = 1,num_energies_neu !> Iterating over the energies in each point
     ! do idire = 1,num_directions_neu
        
         zener =adjustl(intost(iener))
     !    zener =adjustl(intost(idire))
         if( iener < 10 ) then
            wopos(1) = postp(1) % wopos(1,ivari)(1:2) // '00' // trim(zener)
         elseif( iener >= 10 .and. iener < 100 ) then
           wopos(1) = postp(1) % wopos(1,ivari)(1:2)// '0' // trim(zener)
         elseif( iener >= 100 .and. iener < 1000 ) then
            wopos(1) = postp(1) % wopos(1,ivari)(1:2)// trim(zener)
         end if

         if( INOTMASTER ) then
            gesca(1:npoin)=0.0_rp
            do idire = 1,num_directions_neu
               !        do iener = 1,num_energies_neu !> Iterating over the energies in each point
               gesca(1:npoin) = gesca(1:npoin) + neutr(iener,idire,1:npoin,1) * weigd_neu(idire) * kerma_poin_neu(1:ipoin, iener)
             !   gesca(1:npoin) = gesca(1:npoin) + neutr(iener,idire,1:npoin,1)* weigd_neu(idire)*in4pi
            enddo
         end if
 
         call postpr(gesca,wopos,ittim,cutim)
      
      end do
 
 
 !     do iener = 1,num_energies_neu !> Iterating over the energies in each point
 !       do idire = 1,num_directions_neu
 !         
 !        wdire = adjustl(intost(idire))
 !        zener =adjustl(intost(iener))
 !        if( idire < 10 ) then
 !           wopos(1) = postp(1) % wopos(1,3)(1:3) // '0' // trim(wdire)
 !        else
 !           wopos(1) = postp(1) % wopos(1,3)(1:3) // trim(wdire)
 !        end if
 !        if( iener < 10 ) then
 !           wopos(1) = wopos(1) // '0' // trim(zener)
 !        else
 !           wopos(1) = wopos(1) // trim(zener)
 !        end if
 !
 !        if( INOTMASTER ) gesca(1:npoin) = neutr(iener,idire,1:npoin,1)
 !        call postpr(gesca,wopos,ittim,cutim)
 !       enddo
 !     end do
 !
 
      if( INOTMASTER ) call memgen(2_ip,npoin,0_ip)
      return

  end select
  !
  ! Output GESCA for a scalar or GEVEC for a vector
  !
  if(ivari/=0_ip) then
     call outvar(ivari,ittim,cutim,postp(1) % wopos(1,ivari))
  endif

end subroutine neu_outvar
