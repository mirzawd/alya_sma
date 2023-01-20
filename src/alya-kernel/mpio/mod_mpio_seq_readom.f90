!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_seq_readom.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   Mesh sequential reading
!> @details This module reads mesh files in sequential using MPI-IO format
!> @}
!-----------------------------------------------------------------------


module mod_mpio_seq_readom

  use def_master,             only : IIOSLAVE, IMASTER, INOTMASTER, ISEQUEN
  use def_master,             only : namda
  use def_master,             only : intost, retost, kfl_ptask
  use def_kintyp,             only : ip, rp, lg
  use def_domain,             only : coord, ltype, lnods, nnode, ndime, nelem, npoin, nelty, mnode,&
       lnodb, lelbo, ltypb, mnodb, nboun, ndimb, kfl_autbo, kfl_bouel, utype, lmate,&
       leset, lnset, lbset, lesub, kfl_codbo, kfl_codno, kfl_icodb, kfl_icodn, neset, nbset, nnset, mcono,&
       lnoch, lelch, lboch, lmast, nsteps_fiel_ondemand
  use def_domain,             only : nnset_origi,neset_origi,nbset_origi
  use def_domain,             only : memor_dom
  use mod_memory,             only : memory_alloca, memory_deallo
  use mod_messages,           only : livinf
  use mod_domain,             only : domain_memory_allocate
  use mod_messages,           only : messages_live
  use def_mpio
  use mod_mpio_files
  use mod_mpio_config,        only : mpio_config

  use mod_mpio_seq_iowrapper
  use mod_mpio_seq_io


  implicit none

  private

  !-----------------------------------------------------------------------------------------------------------------------------
  !                              Private Variables
  !-----------------------------------------------------------------------------------------------------------------------------


  logical                                 :: skip_ltype=.FALSE.
  logical                                 :: present_ltypb=.FALSE.
  logical                                 :: skip_boundaries=.FALSE.

  integer(ip)                             :: i,j



  public                                  :: reageo_seq_mpio
  public                                  :: reaset_seq_mpio
  public                                  :: reabcs_seq_mpio
  public                                  :: reafie_seq_mpio


contains

  !-----------------------------------------------------------------------
  !
  !> @author    Damien Dosimont
  !> @name      Mesh parallel reading subroutine
  !> @brief     This routine is the entry point to the parallel reading
  !> @details   This routine is the entry point to the parallel reading
  !>
  !-----------------------------------------------------------------------

  subroutine reageo_seq_mpio()
    if (IMASTER .or. ISEQUEN) then
       call messages_live('GEOMETRY READING (SEQUENTIAL)','START SECTION')
       call domain_memory_allocate('GEOMETRY')
       call livinf(0_ip,'READ COORD FILE',0_ip)
       call seq_parse_coord()
       call livinf(0_ip,'READ LNOCH FILE',0_ip)
       call seq_parse_lnoch()
       call livinf(0_ip,'READ LTYPE FILE',0_ip)
       call seq_parse_ltype()
       call livinf(0_ip,'READ LNODS FILE',0_ip)
       call seq_parse_lnods()
       call livinf(0_ip,'READ LMATE FILE',0_ip)
       call seq_parse_lmate()
       call livinf(0_ip,'READ LMAST FILE',0_ip)
       call seq_parse_lmast()
       call livinf(0_ip,'READ LESUB FILE',0_ip)
       call seq_parse_lesub()
       call livinf(0_ip,'READ LELCH FILE',0_ip)
       call seq_parse_lelch()
       if (nboun==0) then
          skip_boundaries=.TRUE.
       end if
       if ((.not.skip_boundaries).and.kfl_autbo /= 1) then
          call livinf(0_ip,'READ LTYPB FILE',0_ip)
          call seq_parse_ltypb()
          call livinf(0_ip,'READ LNODB FILE',0_ip)
          call seq_parse_lnodb()
          call livinf(0_ip,'READ LELBO FILE',0_ip)
          call seq_parse_lelbo()
          call livinf(0_ip,'READ LBOCH FILE',0_ip)
          call seq_parse_lboch()
       end if
       call messages_live('GEOMETRY READING (SEQUENTIAL)','END SECTION')
       !call livinf(-4_ip,'CHECK ARRAYS',0_ip)
       !call livinf(-5_ip,'END CHECK ARRAYS',0_ip)
    end if

  end subroutine reageo_seq_mpio

    subroutine reaset_seq_mpio()
        character(100), PARAMETER :: vacal = "reaset_seq_mpio"
        neset     = 0                   ! There is no element set
        nbset     = 0                   ! There is no boundary set
        nnset     = 0                   ! There is no node set
        if (IMASTER .or. ISEQUEN) then
           call messages_live('SET READING (SEQUENTIAL)','START SECTION')
            if (FILE_EXIST(fil_leset)) then
                call livinf(0_ip,'READ LESET FILE',0_ip)
                call memory_alloca(memor_dom,'leset',vacal,leset,nelem)
                call SEQ_FILE_READ_ALL(leset, fil_leset, nelem)
                neset = 1_ip
            end if
            if (FILE_EXIST(fil_lbset)) then
                call livinf(0_ip,'READ LBSET FILE',0_ip)
                call memory_alloca(memor_dom,'lbset',vacal,lbset,nboun)
                call SEQ_FILE_READ_ALL(lbset, fil_lbset, nboun)
                nbset = 1_ip
            end if
            if (FILE_EXIST(fil_lnset)) then
                call livinf(0_ip,'READ LNSET FILE',0_ip)
                call memory_alloca(memor_dom,'lnset',vacal,lnset,npoin)
                call SEQ_FILE_READ_ALL(lnset, fil_lnset, npoin)
                nnset = 1_ip
            end if
           call messages_live('SET READING (SEQUENTIAL)','END SECTION')
         end if
        neset_origi = neset  
        nbset_origi = nbset   
        nnset_origi = nnset    
         
    end subroutine

    subroutine reabcs_seq_mpio()
        character(100), PARAMETER :: vacal = "reabcs_seq_mpio"
        kfl_icodn = 0
        kfl_icodb = 0
        if (IMASTER .or. ISEQUEN) then
           call messages_live('BCS READING (SEQUENTIAL)','START SECTION')
            if (FILE_EXIST(fil_codno)) then
                call livinf(0_ip,'READ CODNO FILE',0_ip)
                call memory_alloca(memor_dom,'KFL_CODNO',vacal,kfl_codno, mcono, npoin)
                call SEQ_FILE_READ_ALL(kfl_codno, fil_codno, mcono, npoin)
                kfl_icodn = 1_ip
            end if
            if (FILE_EXIST(fil_codbo)) then
                call livinf(0_ip,'READ CODBO FILE',0_ip)
                call memory_alloca(memor_dom,'codbo',vacal,kfl_codbo,nboun)
                call SEQ_FILE_READ_ALL(kfl_codbo, fil_codbo, nboun)
                kfl_icodb = 1_ip
            end if
           call messages_live('BCS READING (SEQUENTIAL)','END SECTION')
        endif
    end subroutine

    subroutine reafie_seq_mpio()
      use def_master
      use def_domain, only: kfl_field, xfiel, nfiel, time_field
      use mod_opfpos, only: postpr_intto8
      character(100), PARAMETER :: vacal = "reafie_seq_mpio"
      integer(ip)                              :: kdime, size_p, ii,jj
      real(rp), pointer                        :: a1(:)
      real(rp), pointer                        :: a2(:,:)
      integer(ip)                              :: ifiel, istep
      type(mpio_header)                        :: header
      integer(ip)                              :: nsteps

      nullify(a1,a2)
      call messages_live('FIELD READING (SEQUENTIAL)','START SECTION')
      if (IMASTER .or. ISEQUEN) then
         do ifiel = 1,nfiel
            call domain_memory_allocate('XFIEL % A',NUMBER1=ifiel)
            if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
               size_p=nelem
            else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
               size_p=npoin
            else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
               size_p=nboun
            else
               call runend('UNDEFINED FIELD TYPE')
            end if
            kdime = kfl_field(1,ifiel)


            if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_config%output%post_process%export_only) ) then 
               nsteps = kfl_field(4,ifiel)
            else
               nsteps = nsteps_fiel_ondemand
            end if

            !read timesteps for all the files, needed to decide later which wiles to load if ONDEMAND is specified
            do istep = 1,kfl_field(4,ifiel)
               fil_field = adjustl(trim(namda))//trim(field_ext)//"."//postpr_intto8(ifiel)//"."//postpr_intto8(istep)//trim(mpio_ext)
               call SEQ_FILE_READ_HEADER(fil_field, header)
               time_field(ifiel) % a(istep)=header%time
            end do
            do istep = 1,nsteps
               call livinf(0_ip,'READ XFIEL '//postpr_intto8(ifiel)//" "//postpr_intto8(istep),0_ip)
               fil_field = adjustl(trim(namda))//trim(field_ext)//"."//postpr_intto8(ifiel)//"."//postpr_intto8(istep)//trim(mpio_ext)

               if( kdime == 1 ) then
                  call memory_alloca(memor_dom,'a1',vacal,a1,size_p)
                  call SEQ_FILE_READ_ALL(a1, fil_field, size_p)
                  xfiel(ifiel) % a(1,:,istep)=a1(:)
                  call memory_deallo(memor_dom,'a1',vacal,a1)
               else
                  call memory_alloca(memor_dom,'a2',vacal,a2,kdime,size_p)
                  call SEQ_FILE_READ_ALL(a2, fil_field, kdime, size_p)
                  do ii=1,kdime
                     do jj=1,size_p
                        xfiel(ifiel) % a(ii,jj,istep)=a2(ii,jj)
                     end do
                  end do
                  call memory_deallo(memor_dom,'a2',vacal,a2)
               end if
            end do
         end do
      call messages_live('FIELD READING (SEQUENTIAL)','END SECTION')
      end if
    end subroutine reafie_seq_mpio

    subroutine seq_parse_coord()
        call SEQ_FILE_READ_ALL(coord, fil_coord, ndime, npoin)
    end subroutine

    subroutine seq_parse_lnoch()
        if (FILE_EXIST(fil_lnoch)) then
            call SEQ_FILE_READ_ALL(lnoch, fil_lnoch, npoin)
        end if
    end subroutine

    subroutine seq_parse_lelch()
        if (FILE_EXIST(fil_lelch)) then
            call SEQ_FILE_READ_ALL(lelch, fil_lelch, nelem)
        end if
    end subroutine

    subroutine seq_parse_lboch()
        if (FILE_EXIST(fil_lboch)) then
            call SEQ_FILE_READ_ALL(lboch, fil_lboch, nboun)
        end if
    end subroutine

    subroutine seq_parse_ltype()
        call seq_eval_types()
        if (.NOT.skip_ltype) then
            call SEQ_FILE_READ_ALL(ltype, fil_ltype, nelem)
        else
            ltype(1:nelem)=utype
        end if
    end subroutine

    subroutine seq_eval_types()
        if (utype/=-1) then
            skip_ltype=.TRUE.
        end if
    end subroutine

    subroutine seq_parse_lnods()
        call SEQ_FILE_READ_ALL(lnods, fil_lnods, mnode, nelem)
    end subroutine

    subroutine seq_parse_lmate()
        if (FILE_EXIST(fil_lmate)) then
            call SEQ_FILE_READ_ALL(lmate, fil_lmate, nelem)
        end if
    end subroutine

    subroutine seq_parse_lmast()
        if (FILE_EXIST(fil_lmast)) then
            call SEQ_FILE_READ_ALL(lmast, fil_lmast, npoin)
        end if
    end subroutine

    subroutine seq_parse_lesub()
        if (FILE_EXIST(fil_lesub)) then
            call SEQ_FILE_READ_ALL(lesub, fil_lesub, nelem)
        end if
    end subroutine

    subroutine seq_parse_ltypb()
        if (FILE_EXIST(fil_ltypb)) then
            call SEQ_FILE_READ_ALL(ltypb, fil_ltypb, nboun)
            present_ltypb=.true.
        end if
    end subroutine

    subroutine seq_parse_lnodb()
        integer(ip)                                  ::              knode, iblty
        call SEQ_FILE_READ_ALL(lnodb, fil_lnodb, mnodb, nboun)
        if (.not.present_ltypb) then
            do i=1, nboun
                do j=1, mnodb
                    if (lnodb(j,i) /=0 ) then
                        knode=j
                    else
                        exit
                    end if
                end do
                call fintyp(ndimb,knode,iblty)
                ltypb(i) = iblty
            end do
        end if
    end subroutine

    subroutine seq_parse_lelbo()
        if (FILE_EXIST(fil_lelbo)) then
            kfl_bouel=1
            call SEQ_FILE_READ_ALL(lelbo, fil_lelbo, nboun)
        else
            kfl_bouel=0
        end if
    end subroutine

end module


