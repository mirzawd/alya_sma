!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_isochr(icomp)
   !-----------------------------------------------------------------------
   !****f* exmedi/exm_isochr
   ! NAME 
   !    exm_isochr
   ! DESCRIPTION
   ! USES
   ! USED BY
   !    exm_isochr
   !***
   !-----------------------------------------------------------------------

   use      def_parame
   use      def_master
   use      def_domain
   use      def_postpr
   use      def_exmedi
   use      mod_communications, only : PAR_MAX
   use      mod_memory,         only : memory_alloca, memory_deallo


   implicit none
   integer(ip), intent(in)  :: icomp
   integer(ip)              :: ipoin
   real(rp)                 :: xauxi

   if( INOTMASTER) then
      if (fisoc_exm(1)>0.5_rp) then    ! isochrones trigger    
         
         !the deafult behaviour is that it saves only the first upstroke
         do ipoin= 1,npoin

            xauxi = elmag(ipoin,icomp) 
            
            if (fisoc(1,ipoin)<0.0_rp) then !isochrones not saved yet
               if (xauxi > fisoc_exm(2)) then ! isochrones threshold
                  fisoc(1,ipoin) = cutim
               end if
            end if

            if ( (fisoc(2,ipoin) < 0.0_rp) .and. (fisoc(1,ipoin)>epsilon(0.0_rp)) ) then !downstroke not measured yet
               if (elmag_minmax(2,ipoin)>0.0_rp) then !max voltage should cross 0
                  if (xauxi < elmag_minmax(1,ipoin) + 0.1*(elmag_minmax(2,ipoin)-elmag_minmax(1,ipoin))) then ! isochrones threshold
                     fisoc(2,ipoin) = cutim
                  end if
               end if
            end if
         end do !ipoin
      end if  ! isochrones trigger   
      
   end if !INOTMASTER


end subroutine exm_isochr
