!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name    Writing table of particle parameters to a file
!> @file    mod_result_io.f90
!> @author  Constantine Butakoff
!> @date    09/03/2021
!> @brief   Output of particles
!> @details Writing table of particle parameters to a ASCII or binary file
!> @}
!------------------------------------------------------------------------

module mod_result_io
   use def_master
   use def_partis, only : lun_resul_pts, fil_resul_pts

   implicit none
   private

   integer(ip), public, parameter   ::      & ! Format for the pts.res file
      PTS_PTSRES_ASCII         =  0,&
      PTS_PTSRES_BIN           =  1

   integer(ip), parameter   ::      & ! One pts.res for all time steps or one per time step
      PTS_PTSRES_LUMPED        =  0,&
      PTS_PTSRES_SEPARATE      =  1




   type, public :: PTS_RESULT_FILE
      private
      integer(ip)                      :: &
         ptsres_format, &  ! One of the PTS_PTSRES_* values
         ptsres_splitting    ! One of the PTS_PTSRES_* values
      contains
         procedure, public :: open_file          
         procedure, public :: create_file        
         procedure, public :: write              
         procedure, public :: set_format         
         procedure, public :: set_format_ascii   
         procedure, public :: set_format_bin     
         procedure, public :: get_format         
         procedure, public :: set_splitting_off  
         procedure, public :: set_splitting_on   
         procedure, public :: is_splitting_on    
         procedure, public :: flush

         procedure, private :: write_header
   end type PTS_RESULT_FILE

   type(PTS_RESULT_FILE), public :: pts_result_io = PTS_RESULT_FILE( ptsres_format=PTS_PTSRES_ASCII, ptsres_splitting=PTS_PTSRES_LUMPED)

   public :: result_io_exchange_add ! calls exchange_add for everything relevant

   contains 

   subroutine flush(self)
      use mod_iofile, only : iofile_flush_unit
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self

      call iofile_flush_unit(lun_resul_pts)
   end subroutine flush

   subroutine result_io_exchange_add()
      use mod_exchange, only : exchange_add

      implicit none
      call exchange_add( pts_result_io % ptsres_format )
      call exchange_add( pts_result_io % ptsres_splitting )
      
   end subroutine result_io_exchange_add
      

   subroutine write(self, ftime, ilagr, itype, iexist, ppvars)
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self
      real(rp),               intent(in) :: ftime  ! TIME columnn
      integer(ip),            intent(in) :: ilagr  ! ilagr columnn
      integer(ip),            intent(in) :: itype  ! itype columnn
      integer(ip),            intent(in) :: iexist ! exist columnn
      real(rp), dimension(:), intent(in) :: ppvars ! other variables to store

      if ( self % ptsres_format == PTS_PTSRES_BIN ) then !save binary
         write(lun_resul_pts)     ftime, ilagr, itype, iexist, ppvars
      else if (self % ptsres_format == PTS_PTSRES_ASCII ) then !save text
         write(lun_resul_pts,100) ftime, ilagr, itype, iexist, ppvars
      end if 

   100 format(1x,es16.8e3,3(1x,i16),400(1x,es16.8e3))

   end subroutine write


   subroutine open_file(self)
      use def_parame, only : zero
      use mod_iofile, only : iofile
      use mod_outfor, only : outfor

      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self

      !No idea what is kfl_reawr and why it was needed here
      !kfl_reawr = 1    
      if( kfl_rstar == 2_ip ) then !continue

         if ( self % ptsres_format == PTS_PTSRES_ASCII ) then
            call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION','old','formatted','append')
         else if ( self % ptsres_format == PTS_PTSRES_BIN ) then
            open(lun_resul_pts, file=fil_resul_pts, access='stream', position='append')
         end if

      else if (kfl_rstar == 1_ip) then !init

         if ( self % ptsres_format == PTS_PTSRES_ASCII ) then
            call iofile(4_ip,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION','old','unformatted')
         else if ( self % ptsres_format == PTS_PTSRES_BIN ) then
            open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
         end if

         call self % write_header()

      else !clean start

         if ( self % ptsres_format == PTS_PTSRES_ASCII ) then
            call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION')
         else if ( self % ptsres_format == PTS_PTSRES_BIN ) then
            open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
         end if

         call self % write_header()

      end if
      !kfl_reawr = 0

   end subroutine open_file

   subroutine create_file(self)
      use def_parame, only : zero
      use mod_opfpos, only : postpr_intto8
      use mod_iofile, only : iofile

      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self

      if( self % is_splitting_on() ) then

         !
         ! self task is for splitted writing of pts.res
         ! open pts.res.000000<time_step_number> for writing
         ! ignore restart flag since there is no appending
         ! 
         fil_resul_pts = adjustl(trim(namda))//'.'//exmod(modul)//'.'//postpr_intto8(ittim)//'.res'

         close(lun_resul_pts)

         if ( self % ptsres_format == PTS_PTSRES_ASCII ) then
            call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION')
         else if ( self % ptsres_format==PTS_PTSRES_BIN ) then
            open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
         end if

         call self % write_header()
      end if

   end subroutine create_file


   subroutine set_format( self, file_format )
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self
      integer(ip), intent(in) :: file_format ! one of the PTS_PTSRES

      if (file_format .ne. PTS_PTSRES_ASCII .and. &
          file_format .ne. PTS_PTSRES_BIN ) then
            call runend("Unknown file format for particles result file :", file_format)
      end if 

      self % ptsres_format = file_format
   end subroutine set_format


   subroutine set_format_ascii( self )
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self

      self % ptsres_format = PTS_PTSRES_ASCII
   end subroutine set_format_ascii


   subroutine set_format_bin( self )
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self

      self % ptsres_format = PTS_PTSRES_BIN
   end subroutine set_format_bin

   integer(ip) function get_format(self)
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self
      get_format = self % ptsres_format
   end function get_format

   subroutine set_splitting_off(self)
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self
      self % ptsres_splitting = PTS_PTSRES_LUMPED
   end subroutine set_splitting_off

   subroutine set_splitting_on(self)
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self      
      self % ptsres_splitting = PTS_PTSRES_SEPARATE
   end subroutine set_splitting_on

   logical(lg) function is_splitting_on(self)
      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self
      is_splitting_on = ( self % ptsres_splitting == PTS_PTSRES_SEPARATE )
   end function is_splitting_on


   subroutine write_header(self)
      use def_partis
      use mod_outfor, only : outfor

      implicit none
      class(PTS_RESULT_FILE), intent(inout) :: self
      integer(ip) :: ii,ivarp
      character(1000) :: binary_header=""
    
      if( nvarp_pts > 0 ) then
         if ( self % ptsres_format==PTS_PTSRES_ASCII ) then
    
            ii = 0
            call outfor(53_ip,lun_resul_pts,' ')
    
            do ivarp = 1,nvarp_pts
               if( postprocess_list_pts(ivarp) /= 0 ) then
                  ii       = ii+1
                  ioutp(1) = ii
                  coutp(1) = postprocess_name_pts(postprocess_list_pts(ivarp))
                  call outfor(35_ip,lun_resul_pts,' ')
               end if
            end do
            write(lun_resul_pts,1) (adjustr(trim(postprocess_name_pts(postprocess_list_pts(ivarp)))),ivarp=1,nvarp_pts)
    
         else if( self % ptsres_format == PTS_PTSRES_BIN ) then
    
            write(binary_header,2) ip, (adjustr(trim(postprocess_name_pts(postprocess_list_pts(ivarp)))),ivarp=1,nvarp_pts)
            write(lun_resul_pts) binary_header
    
         end if 
      end if 
    
    1 format('#           ',a5,100('            ',a5))
    2 format("PTSRES_I",i1,100(' ',a5))
        
    end subroutine write_header

end module mod_result_io