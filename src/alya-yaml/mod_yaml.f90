!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module mod_yaml

#ifdef ALYA_YAML
   use YAMLInterface
   use YAMLRead
   use def_kintyp, only: ip, rp, lg

   implicit none

   private

   public :: AlyaYAMLFile
   public :: AlyaYAMLMap

   type AlyaYAMLFile
      type(YAMLHandler) :: domain
   contains
      procedure :: open => file_open
      procedure :: close => file_close
      procedure :: has_map => file_has_map
      procedure :: get_map => file_get_map
   end type AlyaYAMLFile

   type AlyaYAMLMap
      type(YAMLMap) :: map
      integer(ip) :: ilabel = 1_ip
   contains
      procedure :: has_label => map_has_label
      procedure :: get_map => map_get_map
      procedure :: get_next_label => map_get_next_label
      procedure :: set_logical => map_set_logical
      procedure :: set_integer => map_set_integer
      procedure :: set_real => map_set_real
      procedure :: set_string => map_set_string
      procedure :: destroy => map_destroy
   end type AlyaYAMLMap

   interface AlyaYAMLMap
      procedure :: new_AlyaYAMLMap
   end interface

contains

   subroutine array_to_string(array, string)
      implicit none
      character(*), intent(in) :: array(:)
      character(:), allocatable, intent(out) :: string
      integer(ip) :: i
      i = size(array) + 1
      allocate (character(len=i) :: string)
      write (string, *) array
      string = trim(adjustl(string))
   end subroutine array_to_string

   subroutine file_open(this, filename)
      implicit none
      class(AlyaYAMLFile), intent(inout) :: this

      character(*), intent(in) :: filename
      this%domain = yaml_open_file(filename)
   end subroutine file_open

   subroutine file_close(this)
      implicit none
      class(AlyaYAMLFile), intent(inout) :: this

      call yaml_close_file(this%domain)
   end subroutine file_close

   function file_has_map(this, label) result(ret)
      implicit none
      class(AlyaYAMLFile), intent(inout) :: this

      character(*), intent(in) :: label
      type(YAMLMap)            :: map
      logical(lg)              :: ret

      map = yaml_start_from_map(this%domain, label)
      ret = .false.
      if (int(size(map%labels), ip) .gt. 0_ip) then
         ret = .true.
      end if
      call map%destroy()
   end function file_has_map

   function file_get_map(this, label) result(alyamap)
      implicit none
      class(AlyaYAMLFile), intent(inout) :: this

      character(*), intent(in) :: label
      type(YAMLMap)            :: map
      type(AlyaYAMLMap)        :: alyamap

      map = yaml_start_from_map(this%domain, label)
      alyamap = AlyaYAMLMap(map)
   end function file_get_map

   type(AlyaYAMLMap) function new_AlyaYAMLMap(map)
      implicit none
      type(YAMLMap), intent(in) :: map

      new_AlyaYAMLMap%map = map
   end function new_AlyaYAMLMap

   function map_has_label(this, label) result(ret)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      character(*), intent(in)  :: label
      logical(lg)               :: ret
      integer(ip)               :: i
      character(:), allocatable :: buffer

      ret = .false.
      do i = 1, int(size(this%map%labels), ip)
         call array_to_string(this%map%labels(i)%str, buffer)
         if (buffer .eq. label) then
            ret = .true.
            deallocate (buffer)
            exit
         else
            deallocate (buffer)
         end if
      end do
   end function map_has_label

   function map_get_next_label(this) result(label)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this
      character(256) :: label
      character(:), allocatable :: alabel

      if (this%ilabel .le. int(size(this%map%labels), ip)) then
         call array_to_string(this%map%labels(this%ilabel)%str, alabel)
         this%ilabel = this%ilabel + 1_ip
         label = trim(adjustl(alabel))
         deallocate(alabel)
      else
         label = ""
      end if
   end function map_get_next_label

   function map_get_map(this, label) result(alyamap)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      character(*), intent(in)  :: label
      type(YAMLMap)             :: map
      type(AlyaYAMLMap)         :: alyamap

      map = this%map%value_map(label)
      alyamap = AlyaYAMLMap(map)
   end function map_get_map

   subroutine map_set_logical(this, label, option)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      character(*), intent(in)   :: label
      logical(lg), intent(inout) :: option
      integer                    :: code
      integer(ip)                :: i
      character(:), allocatable  :: buffer

      do i = 1, int(size(this%map%labels), ip)
         call array_to_string(this%map%labels(i)%str, buffer)
         if (buffer .eq. label) then
            deallocate (buffer)
            call array_to_string(this%map%value_str(label, code), buffer)
            if (buffer .eq. "true") then
               option = .true.
            else
               option = .false.
            end if
            deallocate (buffer)
            exit
         else
            deallocate (buffer)
         end if
      end do
   end subroutine map_set_logical

   subroutine map_set_integer(this, label, option)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      character(*), intent(in)   :: label
      integer(ip), intent(inout) :: option
      integer                    :: code
      integer(ip)                :: i
      character(:), allocatable  :: buffer

      do i = 1, int(size(this%map%labels), ip)
         call array_to_string(this%map%labels(i)%str, buffer)
         if (buffer .eq. label) then
            deallocate (buffer)
            option = int(this%map%value_int(label, code), ip)
            exit
         else
            deallocate (buffer)
         end if
      end do
   end subroutine map_set_integer

   subroutine map_set_real(this, label, option)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      character(*), intent(in)   :: label
      real(rp), intent(inout)    :: option
      integer                    :: code
      integer(ip)                :: i
      character(:), allocatable  :: buffer

      do i = 1, int(size(this%map%labels), ip)
         call array_to_string(this%map%labels(i)%str, buffer)
         if (buffer .eq. label) then
            deallocate (buffer)
            option = real(this%map%value_double(label, code), rp)
            exit
         else
            deallocate (buffer)
         end if
      end do
   end subroutine map_set_real

   subroutine map_set_string(this, label, option)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      character(*), intent(in)    :: label
      character(*), intent(inout) :: option
      integer                     :: code
      integer(ip)                 :: i
      character(:), allocatable   :: buffer

      do i = 1, int(size(this%map%labels), ip)
         call array_to_string(this%map%labels(i)%str, buffer)
         if (buffer .eq. label) then
            deallocate (buffer)
            call array_to_string(this%map%value_str(label, code), buffer)
            option = trim(adjustl(buffer))
            deallocate (buffer)
            exit
         else
            deallocate (buffer)
         end if
      end do
   end subroutine map_set_string

   subroutine map_destroy(this)
      implicit none
      class(AlyaYAMLMap), intent(inout) :: this

      call this%map%destroy()
   end subroutine map_destroy

#endif

end module mod_yaml
