!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    mod_postpr_mpio.f90
!> @author  houzeaux
!> @date    2020-05-15
!> @brief   Postprocess with MPIO
!> @details Postprocess with MPIO. Bridge to mpio modules
!-----------------------------------------------------------------------

module mod_postpr_on_basic_mesh
  
  use def_kintyp_basic,        only : ip,rp,lg
  use def_kintyp_mesh,         only : mesh_type_basic
  use def_master,              only : lun_mesh
  use def_postpr,              only : fsize_pos  
  use def_postpr,              only : postpr_intto8
  use def_postpr,              only : postpr_mesh_variable
  use def_master,              only : namda
  use def_mpio,                only : mpio_ext
  use def_mpio,                only : header_ip
  use def_mpio,                only : mpio_header
  use def_mpio,                only : header_size
  use mod_mpio_seq_io,         only : FILL_HEADER
  use mod_mpio_par_io,         only : PAR_FILE_WRITE_HEADER
  use mod_mpio_par_mpiwrapper, only : PAR_FILE_OPEN_WRITE
  use mod_mpio_par_mpiwrapper, only : PAR_FILE_CLOSE
  use mod_mpio_par_mpiwrapper, only : PAR_FILE_WRITE
  use mod_mpio_par_mpiwrapper, only : PAR_FILE_SET_VIEW
  use mod_mpio_par_mpiwrapper, only : PAR_FILE_SEEK_SET
  use mod_mpio_par_mpiwrapper, only : PAR_FILE_SET_SIZE
  use mod_communications,      only : PAR_SUM
  use mod_communications,      only : PAR_MAX
  use mod_communications,      only : PAR_SEND
  use mod_communications,      only : PAR_RECEIVE
  use mod_postpr_tools,        only : postpr_file_name
  use mod_postpr_tools,        only : postpr_is_mpio
  use mod_iofile,              only : iofile_available_unit
  use mod_iofile,              only : iofile_open_unit
  use mod_iofile,              only : iofile_close_unit
  use mod_iofile,              only : iofile_flush_unit
  use mod_iofile,              only : iofile_opened
  use mod_iofile,              only : iofile_delete_file
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_optional_argument,   only : optional_argument
  use def_mpi
#include "def_mpi.inc"
  implicit none

  private
  
  interface postpr_mesh_mpio
     module procedure postpr_mesh_mpio_IP1,&
          &           postpr_mesh_mpio_IP2,&
          &           postpr_mesh_mpio_RP1,&
          &           postpr_mesh_mpio_RP2
  end interface postpr_mesh_mpio
  
  interface postpr_mesh_alya
     module procedure postpr_mesh_alya_IP1,&
          &           postpr_mesh_alya_IP2,&
          &           postpr_mesh_alya_RP1,&
          &           postpr_mesh_alya_RP2
  end interface postpr_mesh_alya

  interface postpr_mesh_permutation
     module procedure postpr_mesh_permutation_IP1,&
          &           postpr_mesh_permutation_IP2,&
          &           postpr_mesh_permutation_RP1,&
          &           postpr_mesh_permutation_RP2
  end interface postpr_mesh_permutation

  integer(8)   :: memor_loc(2)=0_8
  character(5) :: wopos_loc(5)
  
  public :: postpr_mesh_mpio
  public :: postpr_mesh_alya
  public :: postpr_mesh_file_list
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Postprocess
  !> @details Postprocess variables associated ot mesh MESH
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_mesh_mpio_IP1(mesh,xx,wopos,itste,ttime,TAG1,TAG2)

    class(mesh_type_basic),           intent(in) :: mesh
    character(5),                     intent(in) :: wopos(*)
    integer(ip),            pointer,  intent(in) :: xx(:)
    integer(ip),                      intent(in) :: itste
    real(rp),                         intent(in) :: ttime
    integer(ip),            optional, intent(in) :: TAG1
    integer(ip),            optional, intent(in) :: TAG2
    integer(4)                                   :: comm_size4,my_rank4
    integer(ip)                                  :: nsize
    MY_MPI_COMM                                  :: PAR_COMM
    MY_MPI_FILE                                  :: fh
    character(200)                               :: filename
    type(mpio_header)                            :: header
    integer(8)                                   :: offset
    integer(ip),            pointer              :: yy(:)

    if( mesh % comm % rank4 >= 0_4 ) then
       nullify(yy)
       fh                 = MPI_FILE_NULL
       my_rank4           = mesh % comm % RANK4
       comm_size4         = mesh % comm % SIZE4
       PAR_COMM           = mesh % comm % PAR_COMM_WORLD
       filename           = mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       wopos_loc(1:3)     = wopos(1:3)
       wopos_loc(5)       = 'INT'
       header             = postpr_mesh_header_mpio(mesh,wopos_loc,1_ip)
       header % columns   = int(1_ip,header_ip)
       
       call postpr_mesh_permutation(mesh,xx,wopos,yy) 
       if(      wopos(3) == 'NELEM' ) then
          nsize  = mesh % nelem
          offset = int(mesh % comm % offset_nelem,8)
       else if( wopos(3) == 'NPOIN' ) then
          nsize  = mesh % npoin
          offset = int(mesh % comm % offset_npoin,8)
       end if

       call FILL_HEADER          (header,0_ip,0_ip,0_ip,itste,ttime)
       call PAR_FILE_OPEN_WRITE  (fh, filename, PAR_COMM= PAR_COMM)
       call PAR_FILE_SET_SIZE    (fh, header % file_size)
       call PAR_FILE_WRITE_HEADER(fh, header,RANK4=my_rank4)
       call PAR_FILE_SET_VIEW    (fh, yy,int(header_size,8))
       call PAR_FILE_SEEK_SET    (fh, offset)
       call PAR_FILE_WRITE       (fh, yy,nsize)
       call PAR_FILE_CLOSE       (fh)
       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
       
    end if
    
  end subroutine postpr_mesh_mpio_IP1

  subroutine postpr_mesh_mpio_IP2(mesh,xx,wopos,itste,ttime,TAG1,TAG2)

    class(mesh_type_basic),           intent(in) :: mesh
    character(5),                     intent(in) :: wopos(*)
    integer(ip),            pointer,  intent(in) :: xx(:,:)
    integer(ip),                      intent(in) :: itste
    real(rp),                         intent(in) :: ttime
    integer(ip),            optional, intent(in) :: TAG1
    integer(ip),            optional, intent(in) :: TAG2
    integer(4)                                   :: comm_size4,my_rank4
    integer(ip)                                  :: nsize
    MY_MPI_COMM                                  :: PAR_COMM
    MY_MPI_FILE                                  :: fh
    character(200)                               :: filename
    type(mpio_header)                            :: header
    integer(8)                                   :: offset
    integer(ip),            pointer              :: yy(:,:)

    if( mesh % comm % rank4 >= 0_4 ) then
       nullify(yy)
       fh                 =  MPI_FILE_NULL
       my_rank4           =  mesh % comm % RANK4
       comm_size4         =  mesh % comm % SIZE4
       PAR_COMM           =  mesh % comm % PAR_COMM_WORLD
       filename           =  mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       wopos_loc(1:3)     =  wopos(1:3)
       wopos_loc(5)       = 'INT'
       header             =  postpr_mesh_header_mpio(mesh,wopos_loc,1_ip)       
       header % columns   =  size(xx,1,KIND=header_ip)

       call PAR_MAX(header % columns,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)

       call postpr_mesh_permutation(mesh,xx,wopos,yy) 
       if(      wopos(3) == 'NELEM' ) then
          nsize  = mesh % nelem * size(xx,1)
          offset = int(mesh % comm % offset_nelem * size(xx,1),8)
       else if( wopos(3) == 'NPOIN' ) then
          nsize  = mesh % npoin * size(xx,1)
          offset = int(mesh % comm % offset_npoin * size(xx,1),8)
       end if
              
       call FILL_HEADER          (header,0_ip,0_ip,0_ip,itste,ttime)
       call PAR_FILE_OPEN_WRITE  (fh, filename, PAR_COMM= PAR_COMM)
       call PAR_FILE_SET_SIZE    (fh, header % file_size)
       call PAR_FILE_WRITE_HEADER(fh, header,RANK4=my_rank4)
       call PAR_FILE_SET_VIEW    (fh, yy,int(header_size,8))
       call PAR_FILE_SEEK_SET    (fh, offset)
       call PAR_FILE_WRITE       (fh, yy,nsize)
       call PAR_FILE_CLOSE       (fh)
       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
       
    end if
    
  end subroutine postpr_mesh_mpio_IP2

  subroutine postpr_mesh_mpio_RP1(mesh,xx,wopos,itste,ttime,TAG1,TAG2)

    class(mesh_type_basic),           intent(in) :: mesh
    character(5),                     intent(in) :: wopos(*)
    real(rp),               pointer,  intent(in) :: xx(:)
    integer(ip),                      intent(in) :: itste
    real(rp),                         intent(in) :: ttime
    integer(ip),            optional, intent(in) :: TAG1
    integer(ip),            optional, intent(in) :: TAG2
    integer(4)                                   :: comm_size4,my_rank4
    integer(ip)                                  :: nsize
    MY_MPI_COMM                                  :: PAR_COMM
    MY_MPI_FILE                                  :: fh
    character(200)                               :: filename
    type(mpio_header)                            :: header
    integer(8)                                   :: offset
    real(rp),               pointer              :: yy(:)

    if( mesh % comm % rank4 >= 0_4 ) then
       nullify(yy)
       fh                 =  MPI_FILE_NULL
       my_rank4           =  mesh % comm % RANK4
       comm_size4         =  mesh % comm % SIZE4
       PAR_COMM           =  mesh % comm % PAR_COMM_WORLD
       filename           =  mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       wopos_loc(1:3)     =  wopos(1:3)
       wopos_loc(5)       = 'REAL'
       header             =  postpr_mesh_header_mpio(mesh,wopos_loc,1_ip)       
       header % columns   =  int(1,KIND=header_ip)
       
       call postpr_mesh_permutation(mesh,xx,wopos,yy) 
       if(      wopos(3) == 'NELEM' ) then
          nsize  = mesh % nelem 
          offset = int(mesh % comm % offset_nelem,8)
       else if( wopos(3) == 'NPOIN' ) then
          nsize  = mesh % npoin 
          offset = int(mesh % comm % offset_npoin,8)
       end if

       
       call FILL_HEADER          (header,0_ip,0_ip,0_ip,itste,ttime)
       call PAR_FILE_OPEN_WRITE  (fh, filename, PAR_COMM= PAR_COMM)
       call PAR_FILE_SET_SIZE    (fh, header % file_size)
       call PAR_FILE_WRITE_HEADER(fh, header,RANK4=my_rank4)
       call PAR_FILE_SET_VIEW    (fh, yy,int(header_size,8))
       call PAR_FILE_SEEK_SET    (fh, offset)
       call PAR_FILE_WRITE       (fh, yy,nsize)
       call PAR_FILE_CLOSE       (fh)
       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
       
    end if
    
  end subroutine postpr_mesh_mpio_RP1

  subroutine postpr_mesh_mpio_RP2(mesh,xx,wopos,itste,ttime,TAG1,TAG2)

    class(mesh_type_basic),           intent(in) :: mesh
    character(5),                     intent(in) :: wopos(*)
    real(rp),               pointer,  intent(in) :: xx(:,:)
    integer(ip),                      intent(in) :: itste
    real(rp),                         intent(in) :: ttime
    integer(ip),            optional, intent(in) :: TAG1
    integer(ip),            optional, intent(in) :: TAG2
    integer(4)                                   :: comm_size4,my_rank4
    integer(ip)                                  :: nsize
    MY_MPI_COMM                                  :: PAR_COMM
    MY_MPI_FILE                                  :: fh
    character(200)                               :: filename
    type(mpio_header)                            :: header
    integer(8)                                   :: offset
    real(rp),               pointer              :: yy(:,:)
    
    if( mesh % comm % rank4 >= 0_4 ) then
       nullify(yy)
       fh                 =  MPI_FILE_NULL
       my_rank4           =  mesh % comm % RANK4
       comm_size4         =  mesh % comm % SIZE4
       PAR_COMM           =  mesh % comm % PAR_COMM_WORLD

       filename           =  mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       wopos_loc(1:3)     =  wopos(1:3)
       wopos_loc(5)       = 'REAL'
       header             =  postpr_mesh_header_mpio(mesh,wopos_loc,1_ip)       
       header % columns   =  size(xx,1,KIND=header_ip)

       call PAR_MAX(header % columns,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       
       call postpr_mesh_permutation(mesh,xx,wopos,yy) 
       if(      wopos(3) == 'NELEM' ) then
          nsize  = mesh % nelem * size(xx,1)
          offset = int(mesh % comm % offset_nelem * size(xx,1),8)
       else if( wopos(3) == 'NPOIN' ) then
          nsize  = mesh % npoin * size(xx,1)
          offset = int(mesh % comm % offset_npoin * size(xx,1),8)
       end if
       
       call FILL_HEADER          (header,0_ip,0_ip,0_ip,itste,ttime)
       call PAR_FILE_OPEN_WRITE  (fh, filename, PAR_COMM= PAR_COMM)
       call PAR_FILE_SET_SIZE    (fh, header % file_size)
       call PAR_FILE_WRITE_HEADER(fh, header,RANK4=my_rank4)
       call PAR_FILE_SET_VIEW    (fh, yy,int(header_size,8))
       call PAR_FILE_SEEK_SET    (fh, offset)
       call PAR_FILE_WRITE       (fh, yy,nsize)
       call PAR_FILE_CLOSE       (fh)

       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
       
    end if
    
  end subroutine postpr_mesh_mpio_RP2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   File name
  !> @details Compose a file name for postprocess with MPIO
  !> 
  !-----------------------------------------------------------------------

  function postpr_mesh_header_mpio(mesh,wopos,columns) result(header)

    class(mesh_type_basic),           intent(in) :: mesh
    character(5),                     intent(in) :: wopos(5)
    integer(ip),            optional, intent(in) :: columns
    type(mpio_header)                            :: header
    integer(ip)                                  :: nlines

    header % nsubd     = max(1_header_ip,int(mesh % comm % SIZE4,KIND=header_ip))
    header % object    = wopos(1)//'00'//char(0)
    header % resultson = wopos(3)//'00'//char(0)
    
    if(      wopos(3) == 'NPOIN') then
       nlines = mesh % npoin
    else if( wopos(3) == 'NELEM') then
       nlines = mesh % nelem
    else
       call runend('HEADER: CANNOT POSTPROCESS THIS')
    end if
    call PAR_SUM(nlines,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)

    header % lines = int(nlines,KIND=header_ip)
    if( wopos(5) == 'REAL' ) then
       header % size = '8BYTE00'//char(0)
       header % type = 'REAL000'//char(0)       
    else if( wopos(5) == 'INT' ) then
       header % type='INTEG00'//char(0)
       if( ip == 4 ) then
          header % size = '4BYTE00'//char(0)
       else
          header % size = '8BYTE00'//char(0)       
       end if
    end if

    if( present(columns) ) then
       header % columns   = int(columns,header_ip)
    else
       header % columns   = 1_header_ip
    end if

  end function postpr_mesh_header_mpio

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Postprocess
  !> @details Postprocess variables associated ot mesh MESH
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_mesh_alya_IP1(mesh,xx,wopos,itste,ttime,TAG1,TAG2,READ_MODE)

    class(mesh_type_basic),           intent(in)    :: mesh
    character(5),                     intent(in)    :: wopos(*)
    integer(ip),            pointer,  intent(inout) :: xx(:)
    integer(ip),                      intent(in)    :: itste
    real(rp),                         intent(in)    :: ttime
    integer(ip),            optional, intent(in)    :: TAG1
    integer(ip),            optional, intent(in)    :: TAG2
    logical(lg),            optional, intent(in)    :: READ_MODE
    integer(4)                                      :: comm_size4,my_rank4
    MY_MPI_COMM                                     :: PAR_COMM4
    integer(ip)                                     :: pleng,nenti
    integer(ip)                                     :: nunit,ipart,ii,nenti_tot
    character(200)                                  :: filename
    integer(ip),            pointer                 :: yy(:)

    if( mesh % comm % RANK4 >= 0 ) then

       nullify(yy)
       my_rank4       = mesh % comm % RANK4
       comm_size4     = mesh % comm % SIZE4
       PAR_COMM4      = mesh % comm % PAR_COMM_WORLD
       filename       = mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       nunit          = lun_mesh
       wopos_loc(1:3) = wopos(1:3)
       wopos_loc(5)   = 'INT'

       call postpr_mesh_dim(mesh,wopos,nenti)
       nenti_tot = nenti
       call PAR_SUM(nenti_tot,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)

       if( mesh % comm % RANK4 == 0 ) then
          
          call iofile_open_unit(nunit,filename,formo='unformatted')       
          call postpr_mesh_header_alya(mesh,wopos_loc,nunit,itste,ttime,nenti_tot,1_ip,TAG1,TAG2,READ_MODE)

          pleng = nenti

          if( optional_argument(.false.,READ_MODE) ) then
             !
             ! Mater read mode
             !
             do ipart = 1,int(comm_size4,ip)-1
                read(nunit) pleng
                allocate(yy(pleng))
                if( associated(yy) ) read(nunit) ( yy(ii), ii = 1,pleng )
                call PAR_SEND(pleng,   DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                call PAR_SEND(pleng,yy,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                deallocate(yy) 
             end do
             
          else
             !
             ! Master write mode
             !
             call postpr_mesh_permutation(mesh,xx,wopos,yy)           
             write(nunit) pleng
             if( associated(yy) ) write(nunit) ( yy(ii), ii = 1,pleng )
             if( .not. postpr_mesh_variable(wopos(1)) ) call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
             
             do ipart = 1,int(comm_size4,ip)-1
                call PAR_RECEIVE(pleng,   DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                allocate(yy(pleng))
                call PAR_RECEIVE(pleng,yy,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)             
                write(nunit) pleng
                if( pleng > 0 ) write(nunit) ( yy(ii), ii = 1,pleng )
                deallocate(yy)
             end do
             
          end if
          
       else if( mesh % comm % RANK4 > 0 ) then
          
          if( optional_argument(.false.,READ_MODE) ) then
             !
             ! Slave read mode
             !
             call PAR_RECEIVE(nenti,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
             call PAR_RECEIVE(nenti,xx,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
          else
             !
             ! Slave write mode
             !
             call postpr_mesh_permutation(mesh,xx,wopos,yy) 
             call PAR_SEND(nenti,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
             call PAR_SEND(nenti,yy,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
          end if
          
       end if

       call iofile_close_unit(nunit)       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
      
    end if 
   
  end subroutine postpr_mesh_alya_IP1
  
  subroutine postpr_mesh_alya_IP2(mesh,xx,wopos,itste,ttime,TAG1,TAG2,READ_MODE)
    class(mesh_type_basic),           intent(in)    :: mesh
    character(5),                     intent(in)    :: wopos(*)
    integer(ip),            pointer,  intent(inout) :: xx(:,:)
    integer(ip),                      intent(in)    :: itste
    real(rp),                         intent(in)    :: ttime
    integer(ip),            optional, intent(in)    :: TAG1
    integer(ip),            optional, intent(in)    :: TAG2
    logical(lg),            optional, intent(in)    :: READ_MODE
    integer(4)                                      :: comm_size4,my_rank4
    MY_MPI_COMM                                     :: PAR_COMM4
    integer(ip)                                     :: pleng,pdime,nenti
    integer(ip)                                     :: nunit,ipart,ii,kk,nenti_tot
    character(200)                                  :: filename
    integer(ip),            pointer                 :: yy(:,:)

    if( mesh % comm % RANK4 >= 0 ) then

       nullify(yy)
       my_rank4       = mesh % comm % RANK4
       comm_size4     = mesh % comm % SIZE4
       PAR_COMM4      = mesh % comm % PAR_COMM_WORLD
       filename       = mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       nunit          = lun_mesh
       wopos_loc(1:3) =  wopos(1:3)
       wopos_loc(5)   = 'INT'

       call postpr_mesh_dim(mesh,wopos,nenti)
       pdime = size(xx,1)
       nenti_tot = nenti
       call PAR_MAX(pdime,PAR_COMM4)
       call PAR_SUM(nenti_tot,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)

       if( mesh % comm % RANK4 == 0 ) then

          call iofile_open_unit(nunit,filename,formo='unformatted')       
          call postpr_mesh_header_alya(mesh,wopos_loc,nunit,itste,ttime,nenti_tot,pdime,TAG1,TAG2,READ_MODE)

          if( optional_argument(.false.,READ_MODE) ) then
             !
             ! Mater read mode
             !
             read(nunit) pleng
             if( associated(yy) ) read(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
             if( .not. postpr_mesh_variable(wopos(1)) ) call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)                    

             do ipart = 1,int(comm_size4,ip)-1
                read(nunit) pleng
                allocate(yy(pdime,pleng))
                call PAR_SEND(pleng,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                if( pleng > 0 ) read(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
                call PAR_SEND(yy,   DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                deallocate(yy)
             end do
             
          else
             !
             ! Master write mode
             !
             call postpr_mesh_permutation(mesh,xx,wopos,yy)
             pleng = nenti
             write(nunit) pleng
             if( associated(yy) ) write(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
             if( .not. postpr_mesh_variable(wopos(1)) ) call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)                    

             do ipart = 1,int(comm_size4,ip)-1
                call PAR_RECEIVE(pleng,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                allocate(yy(pdime,pleng))
                call PAR_RECEIVE(yy,   DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                write(nunit) pleng
                if( pleng > 0 ) write(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
                deallocate(yy)
             end do

          end if

       else if( mesh % comm % RANK4 > 0 ) then

         if( optional_argument(.false.,READ_MODE) ) then
            !
            ! Slave read mode
            !
            call PAR_RECEIVE(nenti,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
            call PAR_RECEIVE(xx,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)            
         else
            !
            ! Slave write mode
            !
            call postpr_mesh_permutation(mesh,xx,wopos,yy) 
            call PAR_SEND(nenti,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
            call PAR_SEND(yy,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
         end if
         
       end if

       call iofile_close_unit(nunit)       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if

    end if

  end subroutine postpr_mesh_alya_IP2
  
  subroutine postpr_mesh_alya_RP1(mesh,xx,wopos,itste,ttime,TAG1,TAG2,READ_MODE)
    class(mesh_type_basic),           intent(in)    :: mesh
    character(5),                     intent(in)    :: wopos(*)
    real(rp),               pointer,  intent(inout) :: xx(:)
    integer(ip),                      intent(in)    :: itste
    real(rp),                         intent(in)    :: ttime
    integer(ip),            optional, intent(in)    :: TAG1
    integer(ip),            optional, intent(in)    :: TAG2
    logical(lg),            optional, intent(in)    :: READ_MODE
    integer(4)                                      :: comm_size4,my_rank4
    MY_MPI_COMM                                     :: PAR_COMM4
    integer(ip)                                     :: pleng,nenti
    integer(ip)                                     :: nunit,ipart,ii,nenti_tot
    character(200)                                  :: filename
    real(rp),               pointer                 :: yy(:)
    
    if( mesh % comm % RANK4 >= 0 ) then
       
       nullify(yy)
       my_rank4   = mesh % comm % RANK4
       comm_size4 = mesh % comm % SIZE4
       PAR_COMM4  = mesh % comm % PAR_COMM_WORLD
       filename   = mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       nunit      = lun_mesh
       wopos_loc(1:3) =  wopos(1:3)
       wopos_loc(5)   = 'REAL'

       call postpr_mesh_dim(mesh,wopos,nenti)
       nenti_tot = nenti
       call PAR_SUM(nenti_tot,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)

       if( mesh % comm % RANK4 == 0 ) then

          call iofile_open_unit(nunit,filename,formo='unformatted')       
          call postpr_mesh_header_alya(mesh,wopos_loc,nunit,itste,ttime,nenti_tot,1_ip,TAG1,TAG2,READ_MODE)
          
          if( optional_argument(.false.,READ_MODE) ) then
             !
             ! Mater read mode
             !
             do ipart = 1,int(comm_size4,ip)-1
                read(nunit) pleng
                allocate(yy(pleng))
                if( associated(yy) ) read(nunit) ( yy(ii), ii = 1,pleng )
                call PAR_SEND(pleng,   DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                call PAR_SEND(pleng,yy,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                deallocate(yy) 
             end do

          else
             !
             ! Master write mode
             !
             call postpr_mesh_permutation(mesh,xx,wopos,yy)
             pleng = nenti
             write(nunit) pleng
             if( associated(yy) ) write(nunit) ( yy(ii), ii=1,pleng)
             if( .not. postpr_mesh_variable(wopos(1)) ) call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
             
             do ipart = 1,int(comm_size4,ip)-1
                call PAR_RECEIVE(pleng,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                allocate(yy(pleng))
                call PAR_RECEIVE(yy,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)             
                write(nunit) pleng
                if( pleng > 0 ) write(nunit) ( yy(ii), ii=1,pleng)
                deallocate(yy)
             end do
             
          end if
          
       else if( mesh % comm % RANK4 > 0 ) then
          
          if( optional_argument(.false.,READ_MODE) ) then
             !
             ! Slave read mode
             !
             call PAR_RECEIVE(nenti,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
             call PAR_RECEIVE(nenti,xx,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
          else
             call postpr_mesh_permutation(mesh,xx,wopos,yy) 
             call PAR_SEND(nenti,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
             call PAR_SEND(yy,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
          end if
          
       end if

       call iofile_close_unit(nunit)       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
      
    end if 

  end subroutine postpr_mesh_alya_RP1
  
  subroutine postpr_mesh_alya_RP2(mesh,xx,wopos,itste,ttime,TAG1,TAG2,READ_MODE)
    class(mesh_type_basic),           intent(in)    :: mesh
    character(5),                     intent(in)    :: wopos(*)
    real(rp),               pointer,  intent(inout) :: xx(:,:)
    integer(ip),                      intent(in)    :: itste
    real(rp),                         intent(in)    :: ttime
    integer(ip),            optional, intent(in)    :: TAG1
    integer(ip),            optional, intent(in)    :: TAG2
    logical(lg),            optional, intent(in)    :: READ_MODE
    integer(4)                                      :: comm_size4,my_rank4
    MY_MPI_COMM                                     :: PAR_COMM4
    integer(ip)                                     :: pleng,pdime,nenti
    integer(ip)                                     :: nunit,ipart,ii,kk,nenti_tot
    character(200)                                  :: filename
    real(rp),               pointer                 :: yy(:,:)

    if( mesh % comm % RANK4 >= 0 ) then

       nullify(yy)
       my_rank4       = mesh % comm % RANK4
       comm_size4     = mesh % comm % SIZE4
       PAR_COMM4      = mesh % comm % PAR_COMM_WORLD
       filename       = mesh_file_name(mesh,wopos(1),itste,TAG1,TAG2)
       nunit          = lun_mesh
       wopos_loc(1:3) =  wopos(1:3)
       wopos_loc(5)   = 'REAL'
       
       call postpr_mesh_dim(mesh,wopos,nenti)
       pdime = size(xx,1)
       nenti_tot = nenti
       call PAR_MAX(pdime,PAR_COMM4)
       call PAR_SUM(nenti_tot,mesh % comm % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)

       if( mesh % comm % RANK4 == 0 ) then

          call iofile_open_unit(nunit,filename,formo='unformatted')       
          call postpr_mesh_header_alya(mesh,wopos_loc,nunit,itste,ttime,nenti_tot,pdime,TAG1,TAG2,READ_MODE)

         if( optional_argument(.false.,READ_MODE) ) then
             !
             ! Mater read mode
             !
             read(nunit) pleng
             if( associated(yy) ) read(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
             if( .not. postpr_mesh_variable(wopos(1)) ) call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)                    

             do ipart = 1,int(comm_size4,ip)-1
                read(nunit) pleng
                allocate(yy(pdime,pleng))
                call PAR_SEND(pleng,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                if( pleng > 0 ) read(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
                call PAR_SEND(yy,   DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                deallocate(yy)
             end do
          else 
             !
             ! Master write mode
             !
             call postpr_mesh_permutation(mesh,xx,wopos,yy)          
             pleng = nenti
             write(nunit) pleng
             if( associated(yy) ) write(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
             if( .not. postpr_mesh_variable(wopos(1)) ) call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
             
             do ipart = 1,int(comm_size4,ip)-1
                call PAR_RECEIVE(pleng,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)
                allocate(yy(pdime,pleng))
                call PAR_RECEIVE(yy,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)             
                write(nunit) pleng
                if( pleng > 0 ) write(nunit) ( (yy(kk,ii),kk=1,pdime), ii=1,pleng)
                deallocate(yy)
             end do
          end if
          
       else if( mesh % comm % RANK4 > 0 ) then
          
         if( optional_argument(.false.,READ_MODE) ) then
            !
            ! Slave read mode
            !
            call PAR_RECEIVE(nenti,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
            call PAR_RECEIVE(yy,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
         else
            call postpr_mesh_permutation(mesh,xx,wopos,yy) 
            call PAR_SEND(nenti,DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
            call PAR_SEND(yy,   DOM_I=0_ip,PAR_COMM_IN=PAR_COMM4)
         end if
         
       end if

       call iofile_close_unit(nunit)       
       if( .not. postpr_mesh_variable(wopos(1)) ) then
          call memory_deallo(memor_loc,'YY','postpr_mesh_permutation',yy)          
       end if
      
    end if 
          
  end subroutine postpr_mesh_alya_RP2
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   File name
  !> @details Compose a file name for postprocess with MPIO
  !>          case-SUBMESH-LTYPE.post.mpio.bin
  !> 
  !-----------------------------------------------------------------------

   subroutine postpr_mesh_dim(mesh,wopos,nenti,pdime) 
     
     class(mesh_type_basic),           intent(in)  :: mesh
     character(5),                     intent(in)  :: wopos(5)
     integer(ip),                      intent(out) :: nenti
     integer(ip),            optional, intent(out) :: pdime

    if(       wopos(3) == 'NPOIN' ) then
        nenti = mesh % npoin
     else if( wopos(3) == 'NELEM' ) then
        nenti = mesh % nelem
     end if

     if( present(pdime) ) then
        if(      wopos(2) == 'SCALA' ) then
           pdime = 1
        else if( wopos(2) == 'VECTO' ) then
           pdime = mesh % ndime
        end if
     end if
     
   end subroutine postpr_mesh_dim

   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2020-05-14
   !> @brief   File name
   !> @details Compose a file name for postprocess with MPIO
   !>          case-SUBMESH-LTYPE.post.mpio.bin
   !> 
   !-----------------------------------------------------------------------
   
   subroutine postpr_mesh_file_list(mesh) 
     
    class(mesh_type_basic), intent(in) :: mesh
    character(200)                     :: dir
    character(200)                     :: prb
    character(500)                     :: fil

    if( mesh % comm % RANK4 == 0 ) then       
       dir = trim(mesh % name)//'/'
       prb = trim(namda)//'-'//trim(mesh % name)
       fil = trim(dir)//trim(prb)//'.post.alyafil'
       call iofile_delete_file(trim(fil),lun_mesh)
       call iofile_open_unit (lun_mesh,trim(fil),'POSTPROCESS LIST OF FILES')
       call iofile_close_unit(lun_mesh,trim(fil),'POSTPROCESS LIST OF FILES')
    end if
  end subroutine postpr_mesh_file_list
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   File name
  !> @details Compose a file name for postprocess with MPIO
  !>          case-SUBMESH-LTYPE.post.mpio.bin
  !> 
  !-----------------------------------------------------------------------

  function mesh_file_name(mesh,wopos,itste,TAG1,TAG2) result(filename)

    class(mesh_type_basic),             intent(in) :: mesh
    character(5),                       intent(in) :: wopos
    integer(ip),                        intent(in) :: itste
    integer(ip),            optional,   intent(in) :: TAG1
    integer(ip),            optional,   intent(in) :: TAG2
    character(len=:),       allocatable            :: filename
    character(len=fsize_pos)                       :: wtags
    character(20)                                  :: ext
    character(20)                                  :: tim
    character(200)                                 :: dir
    character(200)                                 :: prb
    character(200)                                 :: tag
    !
    ! Deal with tags
    !
    wtags = ''
    if( present(TAG1) ) then
       if( TAG1 /= -1 ) then 
          wtags = trim(wtags)//'.'//postpr_intto8(TAG1)
       end if
    end if
    if( present(TAG2) ) then
       if( TAG2 /= -1 ) then
          wtags = trim(wtags)//'.'//postpr_intto8(TAG2)
       end if
    end if
    !
    ! Extension
    !
    if( postpr_is_mpio() ) then
       ext = '.post'//trim(mpio_ext)
    else
       ext = '.post.alyabin'
    end if
    !
    ! Mesh file
    !
    if( postpr_mesh_variable(wopos) ) then
       tim = ''
    else
       tim = '-'//postpr_intto8(itste)
    end if
    !
    ! File name
    !
    dir      = trim(mesh % name)//'/'
    prb      = trim(namda)//'-'//trim(mesh % name)
    tag      = '-'//trim(wopos)//trim(wtags)//adjustl(trim(tim))//trim(ext)
    filename = trim(dir)//trim(prb)//trim(tag)
    !
    ! Add file to list
    !
    if( .not. postpr_mesh_variable(wopos) .and. mesh % comm % RANK4 == 0 ) then
       call iofile_open_unit(lun_mesh,trim(dir)//trim(prb)//'.post.alyafil','POSTPROCESS LIST OF FILES',POSIO='append')
       write(lun_mesh,'(a)') trim(prb)//trim(tag)
       call iofile_close_unit(lun_mesh)
    end if
    
  end function mesh_file_name

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Header
  !> @details Alya header
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_mesh_header_alya(mesh,wopos,nunit,itste,ttime,pleng,pdime,TAG1,TAG2,READ_MODE)

    class(mesh_type_basic),           intent(in) :: mesh
    character(5),                     intent(in) :: wopos(5)
    integer(ip),                      intent(in) :: nunit
    integer(ip),                      intent(in) :: itste
    real(rp),                         intent(in) :: ttime
    integer(ip),                      intent(in) :: pleng
    integer(ip),            optional, intent(in) :: pdime
    integer(ip),            optional, intent(in) :: TAG1
    integer(ip),            optional, intent(in) :: TAG2
    logical(lg),            optional, intent(in) :: READ_MODE
    integer(4)                                   :: nunit4
    integer(4)                                   :: ihead4
    integer(ip),            parameter            :: nwwww=9
    character(5)                                 :: wwwww(nwwww)
    character(8)                                 :: wwww8(nwwww)
    integer(4)                                   :: iiiii(10)=0_4
    real(8)                                      :: rrrrr(10)=0.0_8
    integer(ip)                                  :: tag1_opt
    integer(ip)                                  :: tag2_opt
    integer(ip)                                  :: ii

    if( mesh % comm % RANK4 == 0 ) then

       nunit4   = int(nunit,4)
       ihead4   = 1234_4
       tag1_opt = optional_argument(0_ip,TAG1)
       tag2_opt = optional_argument(0_ip,TAG2)
       !
       ! Type of postprocess
       !
       wwwww(1) = 'ALYA '
       wwwww(2) = 'V0003'
       wwwww(3) = wopos(1)
       wwwww(4) = wopos(2)

       iiiii(5) = int(tag1_opt,4)
       iiiii(6) = int(tag2_opt,4)

       if( trim(wopos(5)) == 'REAL' ) then
          wwwww(6) = 'REAL '
          wwwww(7) = '8BYTE'
       else
          wwwww(6) = 'INTEG'
          if( ip == 8 ) then
             wwwww(7) = '4BYTE'
          else
             wwwww(7) = '8BYTE'             
          end if
       end if

       if( wopos(2) == 'R3PVE' ) then 
          wwwww(5) = 'R3P  '
          iiiii(1) = int(mesh % ndime,4)
       else
          wwwww(5) = wopos(3)
          if( present(pdime) ) then
             iiiii(1) = int(pdime,4)
          else
             if( wopos(2) == 'SCALA') then
                iiiii(1) = 1_4
             else
                iiiii(1) = int(mesh % ndime,4)
             end if
          end if
       end if
       iiiii(2) = int(pleng,4)
       iiiii(4) = int(itste,4)
       rrrrr(1) = ttime
       wwwww(8) = 'PARAL'
       iiiii(3) = max(1_ip,int(mesh % comm % SIZE4,ip))

       wwww8(1) = 'ALYAPOST'
       do ii = 2,nwwww
          wwww8(ii)=wwwww(ii)//'   '
       end do

       if( optional_argument(.false.,READ_MODE) ) then
          read (nunit4) ihead4
          read (nunit4) wwww8(1) ! ALYA
          read (nunit4) wwww8(2) ! Version
          read (nunit4) wwww8(3) ! Variable name
          read (nunit4) wwww8(4) ! Variable type: SCALAR, VECTOR
          read (nunit4) wwww8(5) ! Variable entity: NPOIN, NELEM
          read (nunit4) wwww8(6) ! Real or integer
          read (nunit4) wwww8(7) ! Byte size
          read (nunit4) wwww8(8) 
          read (nunit4) wwww8(9)
          read (nunit4) iiiii(1)
          read (nunit4) iiiii(2)
          read (nunit4) iiiii(3)
          read (nunit4) iiiii(4)
          read (nunit4) iiiii(5)
          read (nunit4) iiiii(6)
          read (nunit4) rrrrr(1)
       else
          write(nunit4) ihead4
          write(nunit4) wwww8(1) ! ALYA
          write(nunit4) wwww8(2) ! Version
          write(nunit4) wwww8(3) ! Variable name
          write(nunit4) wwww8(4) ! Variable type: SCALAR, VECTOR
          write(nunit4) wwww8(5) ! Variable entity: NPOIN, NELEM
          write(nunit4) wwww8(6) ! Real or integer
          write(nunit4) wwww8(7) ! Byte size
          write(nunit4) wwww8(8) 
          write(nunit4) wwww8(9)
          write(nunit4) iiiii(1)
          write(nunit4) iiiii(2)
          write(nunit4) iiiii(3)
          write(nunit4) iiiii(4)
          write(nunit4) iiiii(5)
          write(nunit4) iiiii(6)
          write(nunit4) rrrrr(1)
       end if

    end if

  end subroutine postpr_mesh_header_alya

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Header
  !> @details Alya header
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_mesh_permutation_IP1(mesh,xx,wopos,yy)

    class(mesh_type_basic),           intent(in)    :: mesh
    integer(ip),            pointer,  intent(in)    :: xx(:)
    character(5),                     intent(in)    :: wopos(5)
    integer(ip),            pointer,  intent(inout) :: yy(:)
    integer(ip)                                     :: ii

    if( associated(xx) ) then
       if( postpr_mesh_variable(wopos(1)) ) then
          yy => xx
       else
          if( wopos(3) == 'NPOIN' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,mesh % npoin)
             if( associated(mesh % permn) ) then
                do ii = 1,mesh % npoin
                   yy(ii) = xx(mesh % permn(ii))
                end do
             else
                do ii = 1,mesh % npoin
                   yy(ii) = xx(ii)
                end do
             end if
          else if( wopos(3) == 'NELEM' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,mesh % nelem)
             if( associated(mesh % perme) ) then
                do ii = 1,mesh % nelem
                   yy(ii) = xx(mesh % perme(ii))
                end do
             else
                do ii = 1,mesh % nelem
                   yy(ii) = xx(ii)
                end do
             end if
          end if
       end if
    end if

  end subroutine postpr_mesh_permutation_IP1

  subroutine postpr_mesh_permutation_IP2(mesh,xx,wopos,yy)

    class(mesh_type_basic),           intent(in)    :: mesh
    integer(ip),            pointer,  intent(in)    :: xx(:,:)
    character(5),                     intent(in)    :: wopos(5)
    integer(ip),            pointer,  intent(inout) :: yy(:,:)
    integer(ip)                                     :: ii,pdime

    if( associated(xx) ) then
       if( postpr_mesh_variable(wopos(1)) ) then
          yy => xx
       else
          pdime = size(xx,DIM=1)
          if( wopos(3) == 'NPOIN' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,pdime,mesh % npoin)
             if( associated(mesh % permn) ) then
                do ii = 1,mesh % npoin
                   yy(:,ii) = xx(:,mesh % permn(ii))
                end do
             else
                do ii = 1,mesh % npoin
                   yy(:,ii) = xx(:,ii)
                end do
             end if
          else if( wopos(3) == 'NELEM' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,pdime,mesh % nelem)
             if( associated(mesh % perme) ) then
                do ii = 1,mesh % nelem
                   yy(:,ii) = xx(:,mesh % perme(ii))
                end do
             else
                do ii = 1,mesh % nelem
                   yy(:,ii) = xx(:,ii)
                end do
             end if
          end if
       end if
    end if

  end subroutine postpr_mesh_permutation_IP2

  subroutine postpr_mesh_permutation_RP1(mesh,xx,wopos,yy)

    class(mesh_type_basic),           intent(in)    :: mesh
    real(rp),               pointer,  intent(in)    :: xx(:)
    character(5),                     intent(in)    :: wopos(5)
    real(rp),               pointer,  intent(inout) :: yy(:)
    integer(ip)                                     :: ii

    if( associated(xx) ) then
       if( postpr_mesh_variable(wopos(1)) ) then
          yy => xx
       else
          if( wopos(3) == 'NPOIN' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,mesh % npoin)
             if( associated(mesh % permn) ) then
                do ii = 1,mesh % npoin
                   yy(ii) = xx(mesh % permn(ii))
                end do
             else
                do ii = 1,mesh % npoin
                   yy(ii) = xx(ii)
                end do
             end if
          else if( wopos(3) == 'NELEM' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,mesh % nelem)
             if( associated(mesh % perme) ) then
                do ii = 1,mesh % nelem
                   yy(ii) = xx(mesh % perme(ii))
                end do
             else
                do ii = 1,mesh % nelem
                   yy(ii) = xx(ii)
                end do
             end if
          end if
       end if
    end if

  end subroutine postpr_mesh_permutation_RP1

  subroutine postpr_mesh_permutation_RP2(mesh,xx,wopos,yy)

    class(mesh_type_basic),           intent(in)    :: mesh
    real(rp),               pointer,  intent(in)    :: xx(:,:)
    character(5),                     intent(in)    :: wopos(5)
    real(rp),               pointer,  intent(inout) :: yy(:,:)
    integer(ip)                                     :: ii,pdime

    if( associated(xx) ) then

       if( postpr_mesh_variable(wopos(1)) ) then
          yy => xx
       else
          pdime = size(xx,DIM=1)
          if( wopos(3) == 'NPOIN' ) then
             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,pdime,mesh % npoin)
             if( associated(mesh % permn) ) then
                do ii = 1,mesh % npoin
                   yy(:,ii) = xx(:,mesh % permn(ii))
                end do
             else
                do ii = 1,mesh % npoin
                   yy(:,ii) = xx(:,ii)
                end do
             end if
          else if( wopos(3) == 'NELEM' ) then

             call memory_alloca(memor_loc,'YY','postpr_mesh_permutation',yy,pdime,mesh % nelem)
             if( associated(mesh % perme) ) then
                do ii = 1,mesh % nelem
                   yy(:,ii) = xx(:,mesh % perme(ii))
                end do
             else
                do ii = 1,mesh % nelem
                   yy(:,ii) = xx(:,ii)
                end do
             end if
          end if
       end if
    end if

  end subroutine postpr_mesh_permutation_RP2

end module mod_postpr_on_basic_mesh
!> @}
