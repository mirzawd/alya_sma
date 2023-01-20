!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_reageo.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   Mesh parallel reading
!> @details This module reads mesh files in parallel using MPI-IO format
!> @}
!-----------------------------------------------------------------------

module mod_mpio_par_readom

    use def_master,                 only : IIOSLAVE, IMASTER, INOTMASTER
    use def_master,                 only : namda
    use def_master,                 only : intost, retost, kfl_ptask
    use def_kintyp,                 only : ip, rp, lg
    use def_domain,                 only : coord, ltype, lnods, nnode, ndime, nelem, npoin, nelty, mnode,&
                                           lnodb, lelbo, ltypb, mnodb, nboun, ndimb, kfl_autbo, kfl_bouel,&
                                           utype, kfl_codbo, kfl_codno, lesub, lboch, lnoch, lelch, lmast,&
                                           lmate, kfl_icodb, kfl_icodn, lbset, lnset, leset,&
                                           mcono, lmast, nsteps_fiel_ondemand
    use def_domain,                 only : nnset,neset,nbset
    use def_domain,                 only : nnset_origi,neset_origi,nbset_origi
    use def_domain,                 only : memor_dom
    use mod_messages,               only : livinf
    use def_mpio
    use mod_mpio_files
    use mod_mpio_config,            only : mpio_config


    use mod_communications,         only : PAR_COMM_SPLIT, PAR_GATHERV, PAR_GATHER, PAR_SEND_RECEIVE
    use mod_communications,         only : PAR_BROADCAST, PAR_BARRIER, PAR_SUM
    use mod_memory,                 only : memory_alloca, memory_deallo
    use mod_memory,                 only : memory_size
    use mod_parall,                 only : PAR_COMM_WORLD
    use mod_parall,                 only : PAR_COMM_MPIO, PAR_COMM_MPIO_WM, PAR_COMM_MPIO_RANK_WM, PAR_COMM_MPIO_WM_SIZE
    use mod_parall,                 only : PAR_MY_CODE_RANK
    use mod_mpio_par_mpiwrapper,    only : PAR_FILE_OPEN_READ, PAR_FILE_CLOSE, PAR_FILE_SEEK_SET, PAR_FILE_READ
    use mod_mpio_par_io,            only : PAR_FILE_READ_ALL, PAR_FILE_WRITE_ALL, PAR_FILE_READ_HEADER, PAR_FILE_OPTION_ENABLED
    use mod_domain,                 only : domain_memory_deallocate
    use mod_domain,                 only : domain_memory_allocate
    use mod_messages,               only : messages_live
    implicit none

    private

!-----------------------------------------------------------------------------------------------------------------------------
!                              Private Variables
!-----------------------------------------------------------------------------------------------------------------------------

    !
    ! Indices
    !

    integer                                 :: i, j                                 !

    integer(ip)                             :: nelem_aver     ! Number of elements of all partitions except last one
    integer(ip)                             :: npoin_aver     ! Number of points of all partitions except last one
    integer(ip)                             :: nboun_aver     ! Number of points of all partitions except last one

    integer(ip)                             :: nelem_all      ! Number of elements of all partitions except last one
    integer(ip)                             :: npoin_all      ! Number of points of all partitions except last one
    integer(ip)                             :: nboun_all      ! Number of points of all partitions except last one

    integer(ip)                             :: nelem_blok_p   ! Number of elements of current partition
    integer(ip)                             :: npoin_blok_p   ! Number of points of current partition
    integer(ip)                             :: nboun_blok_p   ! Number of points of current partition

    integer(ip)                             :: nelem_part_p   ! Number of elements of current partition
    integer(ip)                             :: npoin_part_p   ! Number of points of current partition
    integer(ip)                             :: nboun_part_p   ! Number of points of current partition

    integer(4), pointer                     :: nelem_vec(:)
    integer(4), pointer                     :: npoin_vec(:)
    integer(4), pointer                     :: nboun_vec(:)

    integer(ip)                             :: nbou1,kpoin,kelem,kboun,ipart

    !
    ! Where
    !

    character(150)                          ::  wherein_code="IN MY CODE"                       !
    character(150)                          ::  wherein_p="IN MPIO WITH MASTER"     !

    logical                                 ::  skip_ltype=.FALSE.
    logical                                 ::  present_ltypb=.FALSE.
    logical                                 ::  skip_boundaries=.FALSE.

    interface mpio_gather
      module procedure                          mpio_gather_int_sca,  &
                                                mpio_gather_int_vec,  &
                                                mpio_gather_real_sca, &
                                                mpio_gather_real_vec
    end interface


    public                                  :: nelem, nelem_aver, npoin, npoin_aver, nboun, nboun_aver
    public                                  :: nelem_vec, npoin_vec, nboun_vec
    public                                  :: coord, ltype, lnods
    public                                  :: par_deallocate_vecs
    public                                  :: nelem_blok_p,npoin_blok_p,nboun_blok_p
    public                                  :: nelem_part_p,npoin_part_p,nboun_part_p
    public                                  :: mpio_gather,par_vec

    public                                  :: reageo_par_mpio, reaset_par_mpio, reabcs_par_mpio, reafie_par_mpio, par_readom_finalize
    public                                  :: mpio_initialization,par_init_dimensions,par_init_share

    contains

      subroutine mpio_initialization()

#ifndef MPI_OFF
        use def_mpi, only : MPI_INFO_NULL
#endif
        
        nullify(nelem_vec)
        nullify(npoin_vec)
        nullify(nboun_vec)
#ifndef MPI_OFF
        info = MPI_INFO_NULL
        ierr = 0
#endif
        
      end subroutine mpio_initialization
      
!-----------------------------------------------------------------------
!
!> @author    Damien Dosimont
!> @name      Mesh parallel reading subroutine
!> @brief     This routine is the entry point to the parallel reading
!> @details   This routine is the entry point to the parallel reading
!>
!-----------------------------------------------------------------------

    subroutine reageo_par_mpio()
#ifndef MPI_OFF
        call messages_live('GEOMETRY READING (PARALLEL)','START SECTION')
        call livinf(0_ip,'INIT DIMENSIONS',0_ip)
        call par_init_share()
        call par_init_dimensions()
        call par_allocate_arrays()
        !Points
        call livinf(0_ip,'READ COORD FILE',0_ip)
        call par_parse_coord()
        call livinf(0_ip,'READ LNOCH FILE',0_ip)
        call par_parse_lnoch()
        !Elements
        call livinf(0_ip,'READ LTYPE FILE',0_ip)
        call par_parse_ltype()
        call livinf(0_ip,'READ LNODS FILE',0_ip)
        call par_parse_lnods()
        call livinf(0_ip,'READ LMATE FILE',0_ip)
        call par_parse_lmate()
        call livinf(0_ip,'READ LMAST FILE',0_ip)
        call par_parse_lmast()
        call livinf(0_ip,'READ LESUB FILE',0_ip)
        call par_parse_lesub()
        call livinf(0_ip,'READ LELCH FILE',0_ip)
        call par_parse_lelch()
        !Boundaries
        if ((.not.skip_boundaries).and.kfl_autbo /= 1) then
            call livinf(0_ip,'READ LTYPB FILE',0_ip)
            call par_parse_ltypb()
            call livinf(0_ip,'READ LNODB FILE',0_ip)
            call par_parse_lnodb()
            call livinf(0_ip,'READ LELBO FILE',0_ip)
            call par_parse_lelbo()
            call livinf(0_ip,'READ LBOCH FILE',0_ip)
            call par_parse_lboch()
        end if
        call par_gather_missing()
        call PAR_BARRIER(wherein_p)
        call messages_live('GEOMETRY READING (PARALLEL)','END SECTION')
#endif
        return
    end subroutine

    subroutine reaset_par_mpio()
#ifndef MPI_OFF
      call messages_live('SET READING (PARALLEL)','START SECTION')

        neset       = 0                   ! There is no element set
        nbset       = 0                   ! There is no boundary set
        nnset       = 0                   ! There is no node set
        neset_origi = 0                   ! There is no element set
        nbset_origi = 0                   ! There is no boundary set
        nnset_origi = 0                   ! There is no node set

        if (FILE_EXIST(fil_leset)) then
            call livinf(0_ip,'READ LESET FILE',0_ip)
            call par_parse_int_sca(fil_leset, leset, nelem, VARIABLE_NAME='LESET')
            if (.not. mpio_config%input%parallel_partition) then
               call mpio_gather(leset, leset, nelem_vec, VARIABLE_NAME='LESET')
            end if
            neset = 1_ip
         end if

        if (FILE_EXIST(fil_lbset)) then
            call livinf(0_ip,'READ LBSET FILE',0_ip)
            call par_parse_int_sca(fil_lbset, lbset, nboun, VARIABLE_NAME='LBSET')
            if (.not. mpio_config%input%parallel_partition) then
                call mpio_gather(lbset, lbset, nboun_vec, VARIABLE_NAME='LBSET')
            end if
            nbset = 1_ip
         end if

        if (FILE_EXIST(fil_lnset)) then
            call livinf(0_ip,'READ LNSET FILE',0_ip)
            call par_parse_int_sca(fil_lnset, lnset, npoin, VARIABLE_NAME='LNSET')
            if (.not. mpio_config%input%parallel_partition) then
               call mpio_gather(lnset, lnset, npoin_vec, VARIABLE_NAME='LNSET')
            end if
            nnset = 1_ip
         end if

         neset_origi = neset 
         nbset_origi = nbset 
         nnset_origi = nnset

      call messages_live('SET READING (PARALLEL)','END SECTION')

#endif
    end subroutine

    subroutine reabcs_par_mpio()
#ifndef MPI_OFF
      kfl_icodn = 0
      kfl_icodb = 0
      call messages_live('BCS READING (PARALLEL)','START SECTION')
      if (FILE_EXIST(fil_codno)) then
         call livinf(0_ip,'READ CODNO FILE',0_ip)
         call par_parse_int_vec(fil_codno, kfl_codno, mcono, npoin, VARIABLE_NAME='KFL_CODNO')
         if (.not. mpio_config%input%parallel_partition) then
            call mpio_gather(kfl_codno, kfl_codno, npoin_vec, mcono, VARIABLE_NAME='KFL_CODNO')
         end if
         kfl_icodn = 1_ip
      end if
      if (FILE_EXIST(fil_codbo)) then
         call livinf(0_ip,'READ CODBO FILE',0_ip)
         call par_parse_int_sca(fil_codbo, kfl_codbo, nboun, VARIABLE_NAME='KFL_CODBO')
         if (.not. mpio_config%input%parallel_partition) then
            call mpio_gather(kfl_codbo, kfl_codbo, nboun_vec, VARIABLE_NAME='KFL_CODBO')
         end if
         kfl_icodb = 1_ip
      end if
      call messages_live('BCS READING (PARALLEL)','END SECTION')
#endif
    end subroutine reabcs_par_mpio

    subroutine reafie_par_mpio()
#ifndef MPI_OFF
      use def_master
      use def_domain, only: kfl_field, xfiel, nfiel, time_field
      use mod_opfpos, only: postpr_intto8
      integer(4), pointer                      :: var_vec(:)
      integer(ip)                              :: kdime, size_p,ii,jj
      real(rp), pointer                        :: a1_tmp(:)
      real(rp), pointer                        :: a2_tmp(:,:)
      real(rp), pointer                        :: a1(:)
      real(rp), pointer                        :: a2(:,:)
      integer(ip)                              :: ifiel, istep
      integer(ip)                              :: nsteps
      type(mpio_header)                        :: header


      nullify(a1_tmp,a2_tmp)

      call messages_live('FIELD READING (PARALLEL)','START SECTION')
      do ifiel = 1,nfiel
         if (IMASTER) then
            if(  .not. mpio_config%input%parallel_partition ) then
               nelem=nelem_all
               npoin=npoin_all
               nboun=nboun_all
            end if
         else
         end if
         call domain_memory_allocate('XFIEL % A',NUMBER1=ifiel)
         if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
            var_vec => nelem_vec
            size_p  = nelem
         else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
            var_vec => npoin_vec
            size_p  =  npoin
         else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
            var_vec => nboun_vec
            size_p  =  nboun
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
           fil_field = xfiel_name(ifiel, istep)
           call PAR_FILE_READ_HEADER(fil_field, header, wherein=wherein_p)
           time_field(ifiel) % a(istep)=header%time
         end do



         do istep = 1,nsteps
            call livinf(0_ip,'READ XFIEL '//postpr_intto8(ifiel)//" "//postpr_intto8(istep),0_ip)
            fil_field = xfiel_name(ifiel, istep)
            if( kdime == 1 ) then
               call par_parse_real_sca(fil_field, a1_tmp, size_p,VARIABLE_NAME='XFIEL % A') !, time_field(ifiel) % a(istep))
               if (.not. mpio_config%input%parallel_partition) then
                  call mpio_gather(a1, a1_tmp, var_vec,VARIABLE_NAME='XFIEL % A')
                  if (IMASTER.and.size_p>0) then
                     do jj=1, size_p
                        xfiel(ifiel) % a(1,jj,istep)=a1(jj)
                     end do
                  end if
               else
                  if (ISLAVE.and.size_p>0) then
                     do jj=1, size_p
                        xfiel(ifiel) % a(1,jj,istep)=a1_tmp(jj)
                     end do
                  end if
               end if
            else
               call par_parse_real_vec(fil_field, a2_tmp, kdime, size_p,VARIABLE_NAME='XFIEL % A') !, time_field(ifiel) % a(istep))
                if( size_p > 0 ) then
                  if (.not. mpio_config%input%parallel_partition) then
                     call mpio_gather(a2, a2_tmp, var_vec, kdime,VARIABLE_NAME='XFIEL % A')
                     if (IMASTER.and.size_p>0) then
                        do jj=1,size_p
                           do ii=1,kdime
                              xfiel(ifiel) % a(ii,jj,istep)=a2(ii,jj)
                           end do
                        end do
                     end if
                  else
                     if (ISLAVE.and.size_p>0) then
                        do jj=1,size_p
                           do ii=1,kdime
                              xfiel(ifiel) % a(ii,jj,istep)=a2_tmp(ii,jj)
                           end do
                        end do
                     end if
                  end if
               end if
            end if
         end do
         if (ISLAVE .and. .not. mpio_config%input%parallel_partition) then
            call domain_memory_deallocate('XFIEL % A',NUMBER1=ifiel)
         end if
      end do
      call messages_live('FIELD READING (PARALLEL)','END SECTION')
#endif
    end subroutine reafie_par_mpio

    subroutine par_readom_finalize()
#ifndef MPI_OFF
      if( IMASTER ) then
         if(  .not. mpio_config%input%parallel_partition ) then
            nelem=nelem_all
            npoin=npoin_all
            nboun=nboun_all
         else
            nelem=0
            npoin=0
            nboun=0
         end if
      else
         nelem=0
         npoin=0
         nboun=0
      end if
      call par_deallocate_vecs
#endif
    end subroutine par_readom_finalize


    subroutine par_init_share()
        !
        ! Send sizes and stuff the slaves don't have
        !
        call PAR_BROADCAST(nelem,wherein_code)
        call PAR_BROADCAST(npoin,wherein_code)
        call PAR_BROADCAST(nboun,wherein_code)
#ifndef PNODE_VALUE
        call PAR_BROADCAST(mnode,wherein_code)
#endif
        call PAR_BROADCAST(mnodb,wherein_code)
        call PAR_BROADCAST(ndimb,wherein_code)
        call PAR_BROADCAST(kfl_autbo, wherein_code)
        nelem_all=nelem
        npoin_all=npoin
        nboun_all=nboun
    end subroutine

    subroutine par_init_dimensions()

      if (nboun_all==0) then
         skip_boundaries=.TRUE.
      end if
      nbou1 = max(1_ip,nboun_all)

      if (IIOSLAVE) then

         nelem_aver = nelem_all / PAR_COMM_MPIO_WM_SIZE
         npoin_aver = npoin_all / PAR_COMM_MPIO_WM_SIZE
         nboun_aver = nbou1 / PAR_COMM_MPIO_WM_SIZE
         nelem = nelem_aver
         npoin = npoin_aver
         nboun = nboun_aver

         if( 1 == 1 ) then

            nelem = nelem_aver + 1
            npoin = npoin_aver + 1
            nboun = nboun_aver + 1

            nelem_blok_p = nelem
            npoin_blok_p = npoin
            nboun_blok_p = nboun

            ipart = PAR_COMM_MPIO_RANK_WM
            kelem = nelem_all-nelem*ipart
            if( kelem < nelem ) nelem = max(kelem,0_ip)

            kpoin = npoin_all-npoin*ipart
            if( kpoin < npoin ) npoin = max(kpoin,0_ip)

            kboun = nboun_all-nboun*ipart
            if( kboun < nboun ) nboun = max(kboun,0_ip)

         else

            if (mod(nelem, PAR_COMM_MPIO_WM_SIZE)>PAR_COMM_MPIO_RANK_WM) then
               nelem = nelem+1
            end if
            if (mod(npoin, PAR_COMM_MPIO_WM_SIZE)>PAR_COMM_MPIO_RANK_WM) then
               npoin = npoin+1
            end if
            if (mod(nboun, PAR_COMM_MPIO_WM_SIZE)>PAR_COMM_MPIO_RANK_WM) then
               nboun = nboun+1
            end if
         end if

      end if
      if (IMASTER .or. .not. IIOSLAVE) then
         nelem = 0
         npoin = 0
         nboun = 0
      end if
      nelem_part_p=nelem
      npoin_part_p=npoin
      nboun_part_p=nboun
      call par_vec()
    end subroutine par_init_dimensions

    subroutine par_force_dimensions()
        call par_vec()
    end subroutine

    subroutine par_vec()
        character(100), PARAMETER :: vacal = "par_vec"
        nullify(nelem_vec)
        nullify(npoin_vec)
        nullify(nboun_vec)
        call memory_alloca(memor_dom,'nelem_vec',vacal,nelem_vec, PAR_COMM_MPIO_WM_SIZE+1)
        call memory_alloca(memor_dom,'npoin_vec',vacal,npoin_vec, PAR_COMM_MPIO_WM_SIZE+1)
        call memory_alloca(memor_dom,'nboun_vec',vacal,nboun_vec, PAR_COMM_MPIO_WM_SIZE+1)
        call PAR_GATHER(nelem,nelem_vec,wherein_p)
        call PAR_GATHER(npoin,npoin_vec,wherein_p)
        call PAR_GATHER(nboun,nboun_vec,wherein_p)
    end subroutine

    subroutine par_allocate_arrays()
      character(100), PARAMETER :: vacal = "par_allocate_arrays_reageo_par_mpio"
      
      if (IIOSLAVE) then
         call domain_memory_allocate('GEOMETRY')
      end if
    end subroutine par_allocate_arrays

    subroutine par_deallocate_vecs()
        character(100), PARAMETER :: vacal = "par_deallocate_vecs_reageo_par_mpio"
        !if(associated(nelem_vec)) call memory_deallo(memor_dom,'nelem_vec',vacal,nelem_vec)
        !if(associated(npoin_vec)) call memory_deallo(memor_dom,'npoin_vec',vacal,npoin_vec)
        !if(associated(nboun_vec)) call memory_deallo(memor_dom,'nboun_vec',vacal,nboun_vec)
    end subroutine

    subroutine par_parse_coord()
        call PAR_FILE_READ_ALL(coord, fil_coord, ndime, npoin, wherein=wherein_p)
        if (.not. mpio_config%input%parallel_partition) then
            call mpio_gather(coord, coord, npoin_vec, ndime,VARIABLE_NAME='COORD')
        end if
    end subroutine

    subroutine par_parse_lnoch()
      if (FILE_EXIST(fil_lnoch)) then
          call PAR_FILE_READ_ALL(lnoch, fil_lnoch, npoin, wherein=wherein_p)
      end if
      if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lnoch, lnoch, npoin_vec,VARIABLE_NAME='LNOCH')
      end if
    end subroutine

    subroutine par_parse_lelch()
      if (FILE_EXIST(fil_lelch)) then
          call PAR_FILE_READ_ALL(lelch, fil_lelch, nelem, wherein=wherein_p)
      end if
      if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lelch, lelch, nelem_vec,VARIABLE_NAME='LELCH')
      end if
    end subroutine

    subroutine par_parse_lboch()
      if (FILE_EXIST(fil_lboch)) then
          call PAR_FILE_READ_ALL(lboch, fil_lboch, nboun, wherein=wherein_p)
      end if
      if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lboch, lboch, nboun_vec,VARIABLE_NAME='LBOCH')
      end if
    end subroutine

    subroutine par_parse_ltype()
      integer(ip) :: ielem
        call par_eval_types()
        if (.NOT.skip_ltype) then
            call PAR_FILE_READ_ALL(ltype, fil_ltype, nelem, wherein=wherein_p)
         else
            if(IIOSLAVE) then
               do ielem = 1,nelem
                  ltype(ielem)=utype
               end do
            end if
        end if
        if (.not. mpio_config%input%parallel_partition) then
            call mpio_gather(ltype, ltype, nelem_vec,VARIABLE_NAME='LTYPE')
        end if
    end subroutine

    subroutine par_eval_types()
        if (IMASTER) then
            if (utype/=-1) then
                skip_ltype=.TRUE.
            end if
        end if
        call PAR_BROADCAST(skip_ltype,wherein_code)
        call PAR_BROADCAST(utype,wherein_code)
    end subroutine

    subroutine par_parse_lnods()
        call PAR_FILE_READ_ALL(lnods, fil_lnods, mnode, nelem, wherein=wherein_p)
        if (.not. mpio_config%input%parallel_partition) then
            call mpio_gather(lnods, lnods, nelem_vec, mnode,VARIABLE_NAME='LNODS')
        end if
    end subroutine

    subroutine par_parse_lmate()
      if (FILE_EXIST(fil_lmate)) then
          call PAR_FILE_READ_ALL(lmate, fil_lmate, nelem, wherein=wherein_p)
      end if
      if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lmate, lmate, nelem_vec,VARIABLE_NAME='LMATE')
      end if
    end subroutine

    subroutine par_parse_lmast()
      if (FILE_EXIST(fil_lmast)) then
          call PAR_FILE_READ_ALL(lmast, fil_lmast, npoin, wherein=wherein_p)
      end if
      if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lmast, lmast, npoin_vec,VARIABLE_NAME='LMAST')
      end if
    end subroutine

    subroutine par_parse_lesub()
      if (FILE_EXIST(fil_lesub)) then
          call PAR_FILE_READ_ALL(lesub, fil_lesub, nelem, wherein=wherein_p)
      end if
      if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lesub, lesub, nelem_vec,VARIABLE_NAME='LESUB')
      end if
    end subroutine

    subroutine par_parse_ltypb()
      if (FILE_EXIST (fil_ltypb)) then
          call PAR_FILE_READ_ALL(ltypb, fil_ltypb, nboun, wherein=wherein_p)
          present_ltypb=.true.
      end if
    end subroutine

    subroutine par_parse_lnodb()
      integer(ip)                                  ::              knode, iblty
        call PAR_FILE_READ_ALL(lnodb, fil_lnodb, mnodb, nboun, wherein=wherein_p)
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
        if (.not. mpio_config%input%parallel_partition) then
          call mpio_gather(lnodb, lnodb, nboun_vec, mnodb,VARIABLE_NAME='LNODB')
          call mpio_gather(ltypb, ltypb, nboun_vec,VARIABLE_NAME='LTYPB')
        end if
    end subroutine

    subroutine par_parse_lelbo()
        if (FILE_EXIST (fil_lelbo)) then
            kfl_bouel=1
            call PAR_FILE_READ_ALL(lelbo, fil_lelbo, nboun, wherein=wherein_p)
        else
            kfl_bouel=0
         end if
        if (.not. mpio_config%input%parallel_partition) then
            call mpio_gather(lelbo, lelbo, nboun_vec,VARIABLE_NAME='LELBO')
         end if         
    end subroutine

    subroutine par_gather_missing()
        if (.not. mpio_config%input%parallel_partition) then
            if (.not.((.not.skip_boundaries).and.kfl_autbo /= 1)) then
              ! ugly as fuck
              if (IMASTER) then
                call memory_alloca(memor_dom,'LNODB'    ,'memgeo' , lnodb     , mnodb   , nbou1 )
                call memory_alloca(memor_dom,'LTYPB'    ,'memgeo' , ltypb     , nbou1   )
                call memory_alloca(memor_dom,'LBOCH'    ,'memgeo' , lboch     , nbou1   )
                call memory_alloca(memor_dom,'LELBO'    ,'memgeo' , lelbo     , nbou1   )
              end if
            end if
        end if
    end subroutine

    subroutine par_parse_int_sca(fil, var_p, size_p,VARIABLE_NAME)
      character(*),           intent(in)           ::    fil
      integer(ip), pointer,   intent(inout)        ::    var_p(:)
      integer(ip),            intent(in)           ::    size_p
      character(100), PARAMETER                    :: vacal = "par_parse_int_sca"
      character(len=*),       intent(in), optional :: VARIABLE_NAME
      nullify(var_p)

      if (IIOSLAVE) then
         if( memory_size(var_p) /= size_p ) then
            if( present(VARIABLE_NAME) ) then
               call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
               call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_p,size_p)
            else
               call memory_deallo(memor_dom,'var_p',vacal,var_p)
               call memory_alloca(memor_dom,'var_p',vacal,var_p,size_p)
            end if
         end if
      end if

      call PAR_FILE_READ_ALL(var_p, fil, size_p, wherein=wherein_p, force_subd=.true.)

    end subroutine par_parse_int_sca

    subroutine par_parse_int_vec(fil, var_p, size1_p, size2_p, time,VARIABLE_NAME)
        character(*),           intent(in)            ::    fil
        integer(ip), pointer,   intent(inout)         ::    var_p(:,:)
        integer(ip),            intent(in)            ::    size1_p, size2_p
        real(rp),               intent(out), optional ::    time
        character(len=*),       intent(in),  optional :: VARIABLE_NAME
        type(mpio_header)                             :: header
        character(100), PARAMETER                     :: vacal = "par_parse_real_vec"
        
        nullify(var_p)
        if (IIOSLAVE) then
           if( memory_size(var_p) /= size1_p*size2_p ) then
              if( present(VARIABLE_NAME) ) then
                 call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
                 call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_p,size1_p,size2_p)
              else
                 call memory_deallo(memor_dom,'var_p',vacal,var_p)
                 call memory_alloca(memor_dom,'var_p',vacal,var_p,size1_p,size2_p)
              end if
           end if
        end if
        if (present(time)) then
            call PAR_FILE_READ_HEADER(fil, header, wherein_p)
            time=header%time
        end if
        call PAR_FILE_READ_ALL(var_p, fil, size1_p, size2_p, wherein=wherein_p, force_subd=.true.)
    end subroutine

    subroutine par_parse_real_sca(fil, var_p, size_p, time,VARIABLE_NAME)
        character(*),           intent(in)            ::    fil
        real(rp), pointer,      intent(inout)         ::    var_p(:)
        integer(ip),            intent(in)            ::    size_p
        real(rp),               intent(out), optional ::    time
        character(len=*),       intent(in),  optional :: VARIABLE_NAME
        type(mpio_header)                             :: header
        character(100), PARAMETER                     :: vacal = "par_parse_real_sca"
        
        nullify(var_p)
        if (IIOSLAVE) then
           if( memory_size(var_p) /= size_p ) then
              if( present(VARIABLE_NAME) ) then
                 call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
                 call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_p,size_p)
              else
                 call memory_deallo(memor_dom,'var_p',vacal,var_p)
                 call memory_alloca(memor_dom,'var_p',vacal,var_p,size_p)
              end if
           end if
        end if
        if (present(time)) then
            call PAR_FILE_READ_HEADER(fil, header, wherein_p)
            time=header%time
        end if
        call PAR_FILE_READ_ALL(var_p, fil, size_p, wherein=wherein_p, force_subd=.true.)
    end subroutine

    subroutine par_parse_real_vec(fil, var_p, size1_p, size2_p, time,VARIABLE_NAME)
        character(*),           intent(in)            ::    fil
        real(rp), pointer,      intent(inout)         ::    var_p(:,:)
        integer(ip),            intent(in)            ::    size1_p, size2_p
        real(rp),               intent(out), optional ::    time
        character(len=*),       intent(in),  optional :: VARIABLE_NAME
        type(mpio_header)                             :: header
        character(100), PARAMETER                     :: vacal = "par_parse_real_vec"
        nullify(var_p)
        if (IIOSLAVE) then
           !if( memory_size(var_p) /= size1_p*size2_p ) then
              if( present(VARIABLE_NAME) ) then
                 call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
                 call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_p,size1_p,size2_p)
              else
                 call memory_deallo(memor_dom,'var_p',vacal,var_p)
                 call memory_alloca(memor_dom,'var_p',vacal,var_p,size1_p,size2_p)
              end if
           !end if
        end if
        if (present(time)) then
            call PAR_FILE_READ_HEADER(fil, header, wherein_p)
            time=header%time
        end if
        call PAR_FILE_READ_ALL(var_p, fil, size1_p, size2_p, wherein=wherein_p, force_subd=.true.)
    end subroutine

    subroutine mpio_gather_int_sca(var, var_p, var_vec,VARIABLE_NAME)
    
        integer(ip), pointer,  intent(inout)          ::    var(:)
        integer(ip), pointer                          ::    var_tmp(:)
        integer(ip), pointer,  intent(inout)          ::    var_p(:)
        integer(4),  pointer,  intent(in)             ::    var_vec(:)
        character(len=*),      intent(in),   optional :: VARIABLE_NAME
        integer(ip)                                   ::    size_p
        character(100), PARAMETER                     :: vacal = "par_gather_int_sca"

        if (IMASTER) then
            nullify(var)
        end if
        nullify(var_tmp)
        size_p=sum(var_vec)
        if (IMASTER) then
           if( present(VARIABLE_NAME) ) then
              call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_tmp,size_p)
           else
              call memory_alloca(memor_dom,'var_tmp',vacal,var_tmp,size_p)
           end if
        end if
        call PAR_GATHERV(var_p, var_tmp, var_vec,wherein_p)
        if( present(VARIABLE_NAME) ) then
           call memory_deallo(memor_dom,'var_p',vacal,var_p)
        else
           call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)           
        end if
        if (IMASTER) then
            var => var_tmp
         end if

    end subroutine

    subroutine mpio_gather_real_sca(var, var_p, var_vec,VARIABLE_NAME)
        real(rp),   pointer,   intent(inout)          ::    var(:)
        real(rp),   pointer                           ::    var_tmp(:)
        real(rp),   pointer,   intent(inout)          ::    var_p(:)
        integer(4), pointer,   intent(in)             ::    var_vec(:)
        character(len=*),      intent(in),   optional :: VARIABLE_NAME
        integer(ip)                                   ::    size_p
        character(100), PARAMETER                     :: vacal = "par_gather_int_sca"
        if (IMASTER) then
            nullify(var)
        end if
        nullify(var_tmp)
        size_p=sum(var_vec)
        if (IMASTER) then
           if( present(VARIABLE_NAME) ) then
              call memory_alloca(memor_dom,'var_tmp',vacal,var_tmp,size_p)
           else
              call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_tmp,size_p)
           end if
        end if
        call PAR_GATHERV(var_p, var_tmp, var_vec, wherein_p)
        if( present(VARIABLE_NAME) ) then
           call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
        else
           call memory_deallo(memor_dom,'var_p',vacal,var_p)
        end if
        if (IMASTER) then
            var => var_tmp
        end if
    end subroutine

    subroutine mpio_gather_int_vec(var, var_p, var_vec, size1_p,VARIABLE_NAME)
      integer(ip), pointer,   intent(inout) ::    var(:,:)
      integer(ip), pointer                  ::    var_tmp(:,:)
      integer(ip), pointer,   intent(inout) ::    var_p(:,:)
      integer(4), pointer, intent(in)       ::    var_vec(:)
      integer(4), pointer                   ::    var_vec2(:)
      integer(ip),         intent(in)       ::    size1_p
        character(len=*),      intent(in),   optional :: VARIABLE_NAME
      integer(ip)                           ::    size2_p
      character(100), PARAMETER  :: vacal = "par_gather_real_vec"

      if (IMASTER) then
         nullify(var)
      end if
      nullify(var_tmp)
      nullify(var_vec2)
      size2_p=sum(var_vec)
      if( present(VARIABLE_NAME) ) then
         call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_vec2,PAR_COMM_MPIO_WM_SIZE+1)
      else
         call memory_alloca(memor_dom,'var_vec2',vacal,var_vec2,PAR_COMM_MPIO_WM_SIZE+1)
      end if
      var_vec2=var_vec*size1_p
      if (IMASTER) then
         if( present(VARIABLE_NAME) ) then
            call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_tmp,size1_p,size2_p)
         else
            call memory_alloca(memor_dom,'var',vacal,var_tmp,size1_p,size2_p)
         end if
      end if
      call PAR_GATHERV(var_p, var_tmp, var_vec2, wherein_p)
      if( present(VARIABLE_NAME) ) then
         call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
         call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_vec2)
      else
         call memory_deallo(memor_dom,'var_p',vacal,var_p)
         call memory_deallo(memor_dom,'var_vec2',vacal,var_vec2)
      end if
      if (IMASTER) then
         var => var_tmp
      end if

    end subroutine mpio_gather_int_vec

    subroutine mpio_gather_real_vec(var, var_p, var_vec, size1_p,VARIABLE_NAME)
      real(rp), pointer,   intent(inout)    ::    var(:,:)
      real(rp), pointer                     ::    var_tmp(:,:)
      real(rp), pointer,   intent(inout)    ::    var_p(:,:)
      integer(4), pointer, intent(in)       ::    var_vec(:)
      integer(4), pointer                   ::    var_vec2(:)
      integer(ip),         intent(in)       ::    size1_p
        character(len=*),      intent(in),   optional :: VARIABLE_NAME
      integer(ip)                           ::    size2_p
      character(100), PARAMETER  :: vacal = "par_gather_real_vec"
      if (IMASTER) then
         nullify(var)
      end if
      nullify(var_tmp)
      nullify(var_vec2)
      size2_p=sum(var_vec)
      if( present(VARIABLE_NAME) ) then
         call memory_alloca(memor_dom,trim(VARIABLE_NAME),vacal,var_vec2,PAR_COMM_MPIO_WM_SIZE+1)
      else
         call memory_alloca(memor_dom,'var_vec2',vacal,var_vec2,PAR_COMM_MPIO_WM_SIZE+1)
      end if
      var_vec2=var_vec*size1_p
      if (IMASTER) then
         call memory_alloca(memor_dom,'var',vacal,var_tmp,size1_p,size2_p)
      end if
      
      call PAR_GATHERV(var_p, var_tmp, var_vec2, wherein_p)

      if( present(VARIABLE_NAME) ) then
         call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_p)
         call memory_deallo(memor_dom,trim(VARIABLE_NAME),vacal,var_vec2)
      else
         call memory_deallo(memor_dom,'var_p',vacal,var_p)
         call memory_deallo(memor_dom,'var_vec2',vacal,var_vec2)
      end if
      if (IMASTER) then
         var => var_tmp
      end if
    end subroutine mpio_gather_real_vec

  end module mod_mpio_par_readom


