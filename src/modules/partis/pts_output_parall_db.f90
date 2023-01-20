!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Starts an iteration
!! @file    pts_output.f90
!> @author  Guillaume Houzeaux
!> @date    28/01/2013
!! @brief   Output of particles
!! @details Postprocess on mesh (particle density) and output
!> @}
!------------------------------------------------------------------------

subroutine pts_output_parall_db()
  use def_parame
  use def_master
  use def_kermod
  use def_partis
  use mod_memory
  use mod_ker_timeline,   only : ker_timeline
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_BARRIER
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess_now

#ifdef DBPARTICLES
  use, intrinsic :: iso_c_binding
#endif
  implicit none
  
#ifdef DBPARTICLES
  interface
    subroutine sendparticles (recvcounts, nvar2, kfl_posla_pts, ilimi, cutim, sendbuf_rp) bind ( c )
      use iso_c_binding
      integer ( c_int ), VALUE :: recvcounts
      integer ( c_int ), VALUE :: nvar2
      integer ( c_int ), VALUE :: kfl_posla_pts
      integer ( c_int ), VALUE :: ilimi
      real ( c_double ), VALUE :: cutim
      real ( c_double ) :: sendbuf_rp(*)
    end subroutine sendparticles

    subroutine sendsurfacedepos (cutim, subdomain, sum_parts) bind ( c )
          use iso_c_binding
          real ( c_double ), VALUE :: cutim
          integer ( c_int ), VALUE :: subdomain
          integer ( c_int ), VALUE :: sum_parts
    end subroutine sendsurfacedepos

    subroutine writedeposition (pos_tim, partid, parttype, x, y, z, kfl_exist, pos_set, stk_1, stk_2) bind ( c )
          use iso_c_binding
          real ( c_double ), VALUE :: pos_tim
          integer ( c_int ), VALUE :: partid
          integer ( c_int ), VALUE :: parttype
          real ( c_double ), VALUE :: x
          real ( c_double ), VALUE :: y
          real ( c_double ), VALUE :: z
          integer ( c_int ), VALUE :: kfl_exist
          integer ( c_int ), VALUE :: pos_set
          real ( c_double ), VALUE :: stk_1
          real ( c_double ), VALUE :: stk_2
    end subroutine writedeposition

  end interface
  integer ( c_int ), parameter :: ll = 3
#endif  
 
  integer(ip)           :: ivarp,nvar2,ilagr,ipars,sendsize
#ifdef DBPARTICLES
  integer(ip)           :: ipart
#endif
  integer(ip)           :: kfl_ifpos,kfl_ifdep,ilimi,pos_set
  integer(ip)           :: pos_tim
  integer(4)            :: nvar2_4
  integer(ip), save     :: ittim_last=-1
  real(rp)              :: dt1,cutim_pts
  real(rp), allocatable :: depos_surface(:)


  !----------------------------------------------------------------------
  !
  ! Mesh dependent postprocess
  !
  !----------------------------------------------------------------------

  do ivarp = 1,nvarp
     if( output_postprocess_check_variable_postprocess_now(ivarp) ) &
          call pts_outvar(ivarp)
  end do

  !----------------------------------------------------------------------
  !
  ! Particles output
  !
  !----------------------------------------------------------------------
  !
  ! NVAR2: Number of variables for postprocess
  !
  if(      kfl_posla_pts == 0 ) then
     nvar2 = 17
  else if( kfl_posla_pts == 1 ) then
     nvar2 = 19
  else if( kfl_posla_pts == 2 ) then
     nvar2 = 22
  else if( kfl_posla_pts == 3 ) then
     nvar2 = 31
  end if
  nvar2_4 = int(nvar2,4)
  pos_set = nvar2-1         ! Position of boundary set for deposition
  pos_tim = nvar2
  !
  ! Check what kind of postprocess
  ! KFL_ISPOS = 1: All particles info is required
  ! KFL_IFDEP = 1: Deposited particles info is required
  ! Initial solution if this is a restart is not written
  !
  kfl_ifpos = 0
  kfl_ifdep = 0
  if( mod(ittim, kfl_dbfre_pts) == 0 )       kfl_ifpos = 1
  if( kfl_injec == 1 )                       kfl_ifpos = 1
  if( nlagr_going_out_pts > 0 .and. kfl_oudep_pts /= 0 ) kfl_ifdep = 1

  if( ittyp == ITASK_INITIA .and. kfl_rstar == 2 ) then
     kfl_ifpos = 0
     kfl_ifdep = 0
  end if

  if(      kfl_injec == 1 ) then
     ilimi = -5                                                ! Only just deposited particles information is gathered
  else if( kfl_injec == 0 ) then
     if( kfl_ifpos == 0 .and. kfl_ifdep == 1 ) then
        ilimi = -2                                             ! Only deposited and vanishing particles information is gathered
     else if ( kfl_ifpos == 0 .and. nlagr_going_out_pts > 0 ) then
        ilimi = -7
     else
        ilimi = -1                                             ! Existing, deposited and vanishing particles information is gathered
     end if
  end if
  if( ittim == ittim_last ) kfl_ifpos = 0                      ! Last time was already postprocessed
  !
  ! Current time
  !
  if( kfl_injec == 0 ) then
     cutim_pts = cutim
  else
     cutim_pts = cutim-dtime
  end if

  if( kfl_ifpos == 1 .or. kfl_ifdep == 1 ) then
     !
     ! MPI_GATHER size of particles info
     !
     call ker_timeline('INI_OUTPUT')
     if( kfl_injec == 1 ) then
        ittim_last = -1
     else
        ittim_last = ittim
     end if

     if( IPARALL ) then

        nullify(sendbuf_rp)

         if( ISLAVE ) then
            sendsize = 0
            do ilagr = 1,mlagr
                if( lagrtyp(ilagr) % kfl_exist <= ilimi ) then
                    sendsize = sendsize + nvar2_4
                end if
            end do

           call memory_alloca(mem_modul(1:2,modul),'SENDBUF_RP','pts_output',sendbuf_rp,sendsize,'DO_NOT_INITIALIZE')
           ipars = 0
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist <= ilimi ) then
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = TRANSFER(abs(int(lagrtyp(ilagr) % ilagr,rp)),dt1)     ! 10
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % coord(1)           !  1
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % coord(2)           !  2
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % coord(3)           !  3
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % veloc(1)           !  4
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % veloc(2)           !  5
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % veloc(3)           !  6
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % accel(1)           !  7
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % accel(2)           !  8
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % accel(3)           !  9

                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = real(lagrtyp(ilagr) % itype,rp)     ! 11
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = real(lagrtyp(ilagr) % ielem,rp)     ! 12
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % dt_k               ! 13
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = real(lagrtyp(ilagr) % kfl_exist,rp) ! 14
                 if( kfl_posla_pts >= 1 ) then
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % Cd              ! 15
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % Stk(1)          ! 16
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % Stk(2)          ! 17
                 end if
                 if( kfl_posla_pts >= 2 ) then
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % v_fluid_k(1)    ! 18
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % v_fluid_k(2)    ! 19
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % v_fluid_k(3)    ! 20
                 end if
                 if( kfl_posla_pts >= 3 ) then
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % acced(1)        ! 21
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % acced(2)        ! 22
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % acced(3)        ! 23
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % accee(1)        ! 24
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % accee(2)        ! 25
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % accee(3)        ! 26
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % acceg(1)        ! 27
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % acceg(2)        ! 28
                    ipars             = ipars + 1
                    sendbuf_rp(ipars) = lagrtyp(ilagr) % acceg(3)        ! 29
                 end if
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = real(lagrtyp(ilagr) % boundary_set,rp) ! 30
                 ipars             = ipars + 1
                 sendbuf_rp(ipars) = lagrtyp(ilagr) % t                     ! 31
              end if
           end do
#ifdef DBPARTICLES
           if( kfl_ifpos == 1 ) then
            call sendparticles (ipars, nvar2_4, kfl_posla_pts, ilimi, cutim_pts, sendbuf_rp)
           end if
        else
            !MASTER
            !TODO Do I have the injected particles? Isend is 1 sometimes -> no
            if( kfl_ifpos == 1 ) then
              call sendparticles (0, nvar2_4, kfl_posla_pts, ilimi, cutim_pts, sendbuf_rp)
            end if    
#endif
        end if

        !
        ! MPI_BARRIER
        !
        call PAR_BARRIER()
     end if

     !----------------------------------------------------------------------
     !
     ! Postprocess
     !
     !----------------------------------------------------------------------

#ifdef EVENT_POINTS
     call mpitrace_eventandcounters(500,4)
#endif

  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess deposition
  !
  !----------------------------------------------------------------------
  if( kfl_ifdep == 1 ) then
     if(IPARALL .and. ISLAVE .and. kfl_oudep_pts /= 0) then
        ipars = 0
        !Write to lun_oudep_pts
        do ilagr = 1,mlagr
            !/home/bsc31/bsc31226/Alya_hadrien/Sources/modules/partis/pts_solite.f90 : L1770
           if(int(lagrtyp(ilagr) % ilagr,rp) < 0 .and. (lagrtyp(ilagr) % kfl_exist) <= ilimi) then

#ifdef DBPARTICLES
              call writedeposition(lagrtyp(ilagr) % t, abs(lagrtyp(ilagr) % ilagr), lagrtyp(ilagr) % itype, &
                   lagrtyp(ilagr) % coord(1), lagrtyp(ilagr) % coord(2), lagrtyp(ilagr) % coord(3),         &
                   lagrtyp(ilagr) % kfl_exist, lagrtyp(ilagr) % boundary_set,                               &
                   lagrtyp(ilagr) % Stk(1), lagrtyp(ilagr) % Stk(2))
#endif
           end if
           ipars = ipars + nvar2
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess surface deposition
  !
  !----------------------------------------------------------------------
  if( kfl_ifdep == 1 .AND. kfl_depos_surface_pts ==1 ) then

     allocate( depos_surface(ntyla_pts) )
     call pts_deposition_surface(depos_surface)
#ifdef DBPARTICLES
     !TODO verify when should we call this method (ISLAVE;ISEQ...)
     !TODO ipart should be the slave id
     !Write to lun_depsu_pts
     ipart = 1
     if( IMASTER ) then
        call sendsurfacedepos(cutim_pts,ipart,int(sum(depos_surface)))
     else if( ISEQUEN ) then
        call sendsurfacedepos(cutim_pts,ipart,int(sum(depos_surface)))
     end if
#endif

     deallocate( depos_surface )

  end if
#ifdef EVENT_POINTS
  call mpitrace_eventandcounters(500,0)
#endif
  !
  ! Deallocate memory
  !
  if( IPARALL .and. ( kfl_ifdep == 1 .or. kfl_ifpos == 1 ) ) then
     call memory_deallo(mem_modul(1:2,modul),'SENDBUF_RP','pts_output',sendbuf_rp)
  end if
  !
  ! Formats
  !
  !100 format(e12.6,2x,i7,      6(2x,e12.6),2(2x,i7),40(2x,e12.6))
  !102 format(e12.6,',',i7,',',i7,3(',',e12.6),',',i2,',',i4)
100 format(es16.8e3,2x,i7,      6(2x,es16.8e3),2(2x,i7),40(2x,es16.8e3))
!102 format(es16.8e3,',',i7,',',i7,3(',',es16.8e3),',',i2,',',i4)
102 format(es16.8e3,',',i7,',',i7,3(',',es16.8e3),',',i2,',',i4,3(',',es16.8e3))
103 format(es16.8e3,',',es16.8e3)
end subroutine pts_output_parall_db
