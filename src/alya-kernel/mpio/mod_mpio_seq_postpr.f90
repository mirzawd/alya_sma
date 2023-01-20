!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_seq_postpr.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO post process (parallel)
!> @details This module writes/reads restart and post process files in sequential using the MPI-IO format
!> @}
!-----------------------------------------------------------------------

 module mod_mpio_seq_postpr

    use def_kintyp,             only : ip,rp,lg
    use def_parame
    use def_master
    use def_domain
    use def_elmtyp
    use def_postpr
    use def_mpio
    use mod_mpio_config,        only : mpio_config
    use def_parall
    use mod_mpio_seq_io
    use mod_permut
    use mod_memory,             only : memory_alloca, memory_deallo
    use mod_messages, only : livinf

    implicit none

    private

    character(8)                            :: resultson, object
    integer(ip)                             :: dim
    integer(ip)                             :: imesh
    logical(lg)                             :: kfl_perm


    integer(ip), pointer                    :: isca(:)
    integer(ip), pointer                    :: ivec(:,:)
    real(rp),    pointer                    :: rsca(:)
    real(rp),    pointer                    :: rvec(:,:)

    integer(ip), pointer                    :: permuter(:)

    public                                  :: seq_posmpio_int_v, seq_posmpio_int_m, seq_posmpio_real_v, &
                                               seq_posmpio_real_m, mpio_postpr_init, mpio_postpr_permu_end, &
                                               mpio_postpr_permu_int_v, mpio_postpr_permu_int_m, mpio_postpr_permu_real_v, &
                                               mpio_postpr_permu_real_m, mpio_postpr_header
    public                                  :: resultson, object, dim, imesh, kfl_perm, isca, ivec, rsca, rvec, permuter
    public                                  :: mod_mpio_seq_postpr_initialization
    
  contains

    subroutine mod_mpio_seq_postpr_initialization()

      nullify(isca)
      nullify(ivec)
      nullify(rsca) 
      nullify(rvec)
      nullify(permuter)
      nullify(gescar_mpio)
      nullify(gevecr_mpio)
      nullify(gescai_mpio)
      nullify(geveci_mpio)

    end subroutine mod_mpio_seq_postpr_initialization
      
    subroutine seq_posmpio_int_v()
        type(mpio_header)                   :: header
        call mpio_postpr_init()
        call mpio_postpr_permu_int_v()
        if (kfl_reawr/=1) then
            call livinf(0_ip,'WRITE '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_WRITE_ALL(isca, fil_postp, object, resultson, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=1_ip)
            call mpio_postpr_permu_end()
        else if (kfl_reawr==1) then
            call livinf(0_ip,'READ '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_READ_HEADER(fil_postp, header)
            call mpio_postpr_header(header)
            call SEQ_FILE_READ_ALL(gescai_mpio, fil_postp, dim, header=header)
        end if
    end subroutine

    subroutine seq_posmpio_int_m()
        type(mpio_header)                   :: header
        call mpio_postpr_init()
        call mpio_postpr_permu_int_m()
        if (kfl_reawr/=1) then
            call livinf(0_ip,'WRITE '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_WRITE_ALL(ivec, fil_postp, object, resultson, pdime_mpio, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=1_ip)
            call mpio_postpr_permu_end()
        else if (kfl_reawr==1) then
            call livinf(0_ip,'READ '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_READ_HEADER(fil_postp, header)
            call mpio_postpr_header(header)
            call SEQ_FILE_READ_ALL(geveci_mpio, fil_postp, pdime_mpio, dim, header=header)
        end if
    end subroutine

    subroutine seq_posmpio_real_v()
        type(mpio_header)                   :: header
        call mpio_postpr_init()
        call mpio_postpr_permu_real_v()
        if (kfl_reawr/=1) then
            call livinf(0_ip,'WRITE '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_WRITE_ALL(rsca, fil_postp, object, resultson, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=1_ip)
            call mpio_postpr_permu_end()
        else if (kfl_reawr==1) then
            call livinf(0_ip,'READ '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_READ_HEADER(fil_postp, header)
            call mpio_postpr_header(header)
            call SEQ_FILE_READ_ALL(gescar_mpio, fil_postp, dim, header=header)
        end if
    end subroutine

    subroutine seq_posmpio_real_m()
        type(mpio_header)                   :: header
        call mpio_postpr_init()
        call mpio_postpr_permu_real_m()
        if (kfl_reawr/=1) then
            call livinf(0_ip,'WRITE '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_WRITE_ALL(rvec, fil_postp, object, resultson, pdime_mpio, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=1_ip)
            call mpio_postpr_permu_end()
        else if (kfl_reawr==1) then
            call livinf(0_ip,'READ '//object(1:5)//' (POST/RST - SEQUENTIAL)',0_ip)
            call SEQ_FILE_READ_HEADER(fil_postp, header)
            call mpio_postpr_header(header)
            call SEQ_FILE_READ_ALL(gevecr_mpio, fil_postp, pdime_mpio, dim, header=header)
        end if
    end subroutine


!TO BE SHARED

    subroutine mpio_postpr_dimensions()
        object=wopos_mpio(1)//c5_f
        resultson=wopos_mpio(3)//c5_f
        if( resultson(1:5) == 'NPOIN' ) then
            dim=meshe(imesh) % npoin
        else if( resultson(1:5) == 'NELEM' ) then
            dim=meshe(imesh) % nelem
        else if( resultson(1:5) == 'NBOUN' ) then
            dim=meshe(imesh) % nboun
        else if( resultson(1:5) == 'WHATE' ) then
            dim=pleng_mpio
        else
            call runend("The result type "//resultson(1:5)//" can not be processed by the MPI-IO module...")
        end if
    end subroutine

    subroutine mpio_postpr_init()
        call mpio_postpr_wopos()
        call mpio_postpr_permu_init()
        call mpio_opfpos(tag1_mpio, tag2_mpio)
        call mpio_postpr_dimensions()
   end subroutine

    subroutine mpio_postpr_wopos()
        wopos_pos(1) = wopos_mpio(1)
        wopos_pos(2) = wopos_mpio(2)
        wopos_pos(3) = wopos_mpio(3)
    end subroutine

    subroutine mpio_postpr_permu_init()
        if( kfl_reawr == 0 ) then
           imesh = kfl_posdi
        else
           imesh = ndivi
        end if
        kfl_perm = .false.
        if( kfl_reawr == 0 .and. kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos_pos(3) == 'NPOIN' ) then
            permuter => lpmsh
            kfl_perm = .true.
        else if( kfl_reawr == 0 .and. kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos_pos(3) == 'NELEM' ) then
            permuter => lemsh
            kfl_perm = .true.
        else if( kfl_reawr == 0 .and. kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos_pos(3) == 'NBOUN' ) then
            permuter => lbmsh
            kfl_perm = .true.
        end if
        if (kfl_perm .and. mpio_config%output%merge) then
          call runend("CANNOT USE POST MERGE AND PERMUTATION TOGETHER!")
        end if
    end subroutine

    subroutine mpio_postpr_permu_int_v()
      character(100), PARAMETER :: vacal = "mpio_postpr_permu_int_v"
      integer(ip) :: ii
      if (kfl_perm) then
         nullify(isca)
         call memory_alloca(mpio_memor,'isca', vacal, isca, dim)
         if( dim > 0 ) call permut(dim,permuter,gescai_mpio,isca)
      else if (mpio_config%output%merge) then
         call memory_alloca(mpio_memor,'isca', vacal, isca, dim)
         if (dim>0) then
            do ii = 1,dim
               isca(ii) = gescai_mpio(ii)
            end do
         end if
      else
         isca => gescai_mpio
      end if
    end subroutine mpio_postpr_permu_int_v

    subroutine mpio_postpr_permu_int_m()
      character(100), PARAMETER :: vacal = "mpio_postpr_permu_int_m"
      integer(ip) :: ii
        if (kfl_perm) then
            nullify(ivec)
            call memory_alloca(mpio_memor,'ivec', vacal, ivec, pdime_mpio, dim)
            if( dim > 0 ) call permut(pdime_mpio,dim,permuter,geveci_mpio,ivec)
        else if (mpio_config%output%merge) then
            call memory_alloca(mpio_memor,'ivec', vacal, ivec, pdime_mpio, dim)
            if (dim>0) then
               do ii = 1,dim
                  ivec(1:pdime_mpio,ii) = geveci_mpio(1:pdime_mpio,ii)
               end do
            end if
        else
            ivec => geveci_mpio
        end if
    end subroutine

    subroutine mpio_postpr_permu_real_v()
      character(100), PARAMETER :: vacal = "mpio_postpr_permu_real_v"
      integer(ip) :: ii      
        if (kfl_perm) then
            nullify(rsca)
            call memory_alloca(mpio_memor,'rsca', vacal, rsca, dim)
            if( dim > 0 ) call permut(dim,permuter,gescar_mpio,rsca)
        else if (mpio_config%output%merge) then
            call memory_alloca(mpio_memor,'rsca', vacal, rsca, dim)
            if (dim>0) then
               do ii = 1,dim
                  rsca(ii) = gescar_mpio(ii)
               end do
             end if
        else
            rsca => gescar_mpio
        end if
    end subroutine

    subroutine mpio_postpr_permu_real_m()
      character(100), PARAMETER :: vacal = "mpio_postpr_permu_real_m"
      integer(ip) :: ii
        if (kfl_perm) then
            nullify(rvec)
            call memory_alloca(mpio_memor,'rvec', vacal, rvec, pdime_mpio, dim)
            if( dim > 0 ) call permut(pdime_mpio,dim,permuter,gevecr_mpio,rvec)
        else if (mpio_config%output%merge) then
            call memory_alloca(mpio_memor,'rvec', vacal, rvec, pdime_mpio, dim)
            if (dim>0) then
               do ii = 1,dim
                  rvec(1:pdime_mpio,ii) = gevecr_mpio(1:pdime_mpio,ii)
               end do
            end if
         else
            rvec => gevecr_mpio
        end if
    end subroutine

    subroutine mpio_postpr_permu_end(VARIABLE_NAME)
      
      character(len=*), intent(in), optional :: VARIABLE_NAME
      character(100),   PARAMETER            :: vacal = "mpio_postpr_permu_end"

      if (kfl_perm .or. mpio_config%output%merge) then
         if( present(VARIABLE_NAME) ) then
            if (associated(isca)) call memory_deallo(mpio_memor,trim(VARIABLE_NAME), vacal, isca)
            if (associated(ivec)) call memory_deallo(mpio_memor,trim(VARIABLE_NAME), vacal, ivec)
            if (associated(rsca)) call memory_deallo(mpio_memor,trim(VARIABLE_NAME), vacal, rsca)
            if (associated(rvec)) call memory_deallo(mpio_memor,trim(VARIABLE_NAME), vacal, rvec)
         else
            if (associated(isca)) call memory_deallo(mpio_memor,'isca', vacal, isca)
            if (associated(ivec)) call memory_deallo(mpio_memor,'ivec', vacal, ivec)
            if (associated(rsca)) call memory_deallo(mpio_memor,'rsca', vacal, rsca)
            if (associated(rvec)) call memory_deallo(mpio_memor,'rvec', vacal, rvec)
         end if
       else
          nullify(isca,ivec,rsca,rvec)
        end if
    end subroutine

    subroutine mpio_postpr_header(header)
        type(mpio_header), intent(in)                   :: header
        if (ISEQUEN .or. IMASTER) then
            if( ncoun_pos == 0 ) then
                ncoun_pos = ncoun_pos + 1
                write( lun_pos02,'(a,1x,i9,1x,e12.5)') 'START',header%nsubd,header%time
            end if
            write( lun_pos02,'(a)') header%object(1:5)
            flush(lun_pos02)
        end if
    end subroutine

    subroutine mpio_opfpos(TAG1, TAG2)
  !-----------------------------------------------------------------------
  !****f* Domain/opfpos
  ! NAME
  !    opfpos
  ! DESCRIPTION
  !    This subroutine gets ALL the file names to be used by Alya in two
  !    possible ways and them open them:
  !
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  !
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked".
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
        use def_parame
        use def_master
        use def_domain
        use def_inpout
        use mod_iofile
        use def_postpr
        use mod_opfpos
        implicit none

        integer(ip),  intent(in)            :: TAG1
        integer(ip),  intent(in)            :: TAG2
        logical                             :: mesh

        call opfpos_name(wopos_pos(1), mpio_ext, mesh, TAG1, TAG2)
        if (.not.mesh) then
            if( kfl_reawr == 0 ) then
              !
              ! Open postprocess file name
              !
              if( INOTSLAVE ) then
                 write(lun_pos01,'(a)') trim(fil_postp)
                 flush(lun_pos01)
              end if

          end if

       end if
    end subroutine

end module
