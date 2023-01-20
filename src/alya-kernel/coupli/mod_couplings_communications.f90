!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> Generate transmission matrices and associated communication schemes
!> @{
!> @file    mod_couplings_communications.f90
!> @author  borrell
!> @date    2018-03-27
!> @brief   Coupling transmission matrices
!> @details generates transmission matrices and associated communication schemes
!>          
!-----------------------------------------------------------------------
module mod_couplings_communications

  use def_kintyp_basic,   only : ip,rp,lg,r1p,i1p
  use def_kintyp_comm,    only : comm_data_par  
  use def_spmat,          only : spmat
  use def_master,         only : INOTMASTER,ISEQUEN
  use def_domain,         only : lnods
  use def_domain,         only : lnnod
  use def_domain,         only : npoin
  use def_domain,         only : nelem
  use def_domain,         only : nboun
  use mod_parall,         only : color_target
  use mod_parall,         only : color_source
  use mod_parall,         only : PAR_COMM_CURRENT
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_parall,         only : PAR_COMM_COLOR_ARRAY
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_memory,         only : memory_size
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_interpolation,  only : COU_GENERATE_PROJECTION_MATRIX
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_BARRIER
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_START_NON_BLOCKING_COMM
  use mod_communications, only : PAR_END_NON_BLOCKING_COMM
  use mod_communications, only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE      
  use def_coupli,         only : coupling_type
  use def_coupli,         only : typ_color_coupling
  use def_coupli,         only : ELEMENT_INTERPOLATION
  use def_coupli,         only : BOUNDARY_INTERPOLATION
  use def_coupli,         only : NEAREST_BOUNDARY_NODE
  use def_coupli,         only : NEAREST_ELEMENT_NODE
  use def_coupli,         only : GLOBAL_NUMBERING
  use def_coupli,         only : BOUNDARY_VECTOR_PROJECTION
  use def_coupli,         only : TRANSPOSE_MIRROR
  use def_coupli,         only : STRESS_PROJECTION
  use def_coupli,         only : PROJECTION
  use def_coupli,         only : VALUES_ON_ELEMENTS
  use def_coupli,         only : VALUES_ON_NODES
  use def_coupli,         only : VALUES_ON_BOUNDARIES
  use def_coupli,         only : memor_cou
  use def_coupli,         only : lnodb_cou
  use def_coupli,         only : lnnob_cou
  use mod_matrix,         only : matrix_COO_spgemm
  use mod_matrix,         only : matrix_COO_aggregate
  use def_mpi
#include "def_mpi.inc"
  
  implicit none 
  private

  public :: COU_GENERATE_LOCAL_TRANSMISSION_MATRICES               
  public :: COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES               
  public :: COU_PARALLELIZE_TRANSMISSION_MATRICES               

contains

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    27/03/2018
  !> @brief   Generate local transmission matrices  
  !> @details Generate the local transission matrices for a coupling
  !>
  !----------------------------------------------------------------------
  subroutine COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling,all_interp_)

    implicit none
    type(typ_color_coupling), intent(inout)        :: coupling
    logical(lg),              intent(in), optional :: all_interp_
    logical(lg)                                    :: all_interp
    real(rp)                                       :: time1, time2, val
    integer(ip)                                    :: dom_i,ineig,jelem
    integer(ip)                                    :: inise,finse
    integer(ip)                                    :: kelem,ielem
    integer(ip)                                    :: pnode, inode, pnodb
    integer(ip)                                    :: ipoin, jpoin, kpoin
    integer(ip)                                    :: iboun, jboun, inodb, kboun
    integer(ip)                                    :: kgaub, itype
    integer(ip)                                    :: icont, jcont,kcont
    logical(lg)                                    :: I_AM_SOURCE
    logical(lg)                                    :: I_AM_TARGET
    type(r1p),   pointer                           :: spmat_send(:)
    type(r1p),   pointer                           :: spmat_recv(:)
    integer(ip), pointer                           :: nentr_send(:)
    integer(ip), pointer                           :: nentr_recv(:)
    integer(ip), pointer                           :: npoin_send(:)
    integer(ip), pointer                           :: npoin_recv(:)
    MY_MPI_COMM                                    :: PAR_COMM_SAVE
    integer(ip)                                    :: PAR_CURRENT_SIZE
    integer(ip)                                    :: PAR_CURRENT_RANK
    real(rp)                                       :: dummr(1)
    integer(ip)                                    :: nsend
    integer(ip)                                    :: nrecv
    integer(ip)                                    :: nentr_ineig
    integer(ip), pointer                           :: permm(:)
    type(i1p), pointer                             :: perms(:)
    integer(ip)                                    :: nneig
    integer(ip), pointer                           :: inv_status(:)
    integer(ip)                                    :: cont_recv
    integer(ip)                                    :: cont_wetp
    type(spmat)                                    :: M_proj
    type(spmat)                                    :: M_aux
    integer(ip)                                    :: nenti

    character(100), PARAMETER :: vacal = "COU_GENERATE_LOCAL_TRANSMISSION_MATRICES"
    real(rp),PARAMETER :: tolval = epsilon(1.0_rp) ! tolval = 1.0e-10_rp 
    !--------------------------------------------------------------------
    !
    ! Initializations
    !
    !--------------------------------------------------------------------
    all_interp = .false.
    if(present(all_interp_)) all_interp=all_interp_
    !
    ! Timer
    !
    call cputim(time1)
    !
    ! Target and source colors
    !
    color_target = coupling % color_target
    color_source = coupling % color_source
    I_AM_SOURCE = .false.
    I_AM_TARGET = .false.
    if(I_AM_IN_COLOR(color_source)) then 
       I_AM_SOURCE = .true.
    end if
    if(I_AM_IN_COLOR(color_target)) then
       I_AM_TARGET = .true.
    end if
    !
    ! Communicator
    !
    PAR_COMM_SAVE    = PAR_COMM_CURRENT 
    PAR_COMM_CURRENT = PAR_COMM_COLOR(color_target,color_source) 
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    !
    ! Allocate arrays
    !
    nullify(spmat_send)
    nullify(spmat_recv)
    nullify(nentr_send)
    nullify(nentr_recv)
    nullify(npoin_send)
    nullify(npoin_recv)
    nullify(permm)
    nullify(perms)
    nneig = coupling % commd % nneig


    if(      coupling % kfl_source_value == VALUES_ON_NODES      ) then
       nenti = npoin
    else if( coupling % kfl_source_value == VALUES_ON_ELEMENTS   ) then
       nenti = nelem
    else if( coupling % kfl_source_value == VALUES_ON_BOUNDARIES ) then
       nenti = nboun
    else
       call runend("GENERATE_LOCAL_TRANSMISSION_MATRICES: incorrect source value")
    endif
        
    itype = coupling % itype
    call memory_alloca(memor_cou,'SPMAT_SEND',vacal,spmat_send,nneig)
    call memory_alloca(memor_cou,'SPMAT_RECV',vacal,spmat_recv,nneig)
    call memory_alloca(memor_cou,'PERMS'     ,vacal,perms,nneig)
    call memory_alloca(memor_cou,'NENTR_SEND',vacal,nentr_send,nneig)
    call memory_alloca(memor_cou,'NPOIN_SEND',vacal,npoin_send,nneig)
    call memory_alloca(memor_cou,'NENTR_RECV',vacal,nentr_recv,nneig)
    call memory_alloca(memor_cou,'NPOIN_RECV',vacal,npoin_recv,nneig)
    call memory_alloca(memor_cou,'PERMM'     ,vacal,permm,nenti)

    if(nneig > 0) then
       npoin_send(1:nneig) = 0
       npoin_recv(1:nneig) = 0
       nentr_send(1:nneig) = 0
       nentr_recv(1:nneig) = 0
    endif
    !--------------------------------------------------------------------
    !
    ! Evaluate Transmission matrices
    !
    !--------------------------------------------------------------------

    if(I_AM_SOURCE) then

       if( itype == ELEMENT_INTERPOLATION ) then
          !
          ! Element interpolation
          !
          if( coupling % geome % nelem_source > 0 ) then
             kelem = 0
             do ineig = 1,nneig
                dom_i = coupling % commd % neights(ineig)          
                inise = coupling % commd % lsend_size(ineig) 
                finse = coupling % commd % lsend_size(ineig+1) - 1
                nsend = finse - inise + 1

                icont=0
                jcont=0
                kcont=0

                permm(1:nenti)=0
                if( coupling % kfl_source_value == VALUES_ON_NODES) then
                   do jelem = 1,nsend
                      kelem = kelem + 1
                      jcont = jcont + 1 
                      ielem = coupling % geome % lelem_source(kelem)
                      pnode = lnnod(ielem)
                      do inode = 1,pnode
                         val = coupling % geome % shapf(inode,kelem)
                         if(abs(val) > tolval) then
                            icont = icont +1
                            ipoin = lnods(inode,ielem)
                            if(.not. permm(ipoin) ==1) then
                               permm(ipoin) = 1
                               kcont = kcont +1
                            end if
                         end if
                      end do
                   end do
                else
                   do jelem = 1,nsend
                      kelem = kelem + 1
                      jcont = jcont + 1
                      icont = icont + 1 
                      ielem = coupling % geome % lelem_source(kelem)
                      if( .not. permm(ielem) == 1 ) then
                         permm(ielem) = 1
                         kcont = kcont +1
                      end if
                   end do
                endif
                kelem = kelem - jcont

                call memory_alloca(memor_cou,'SPMAT_SEND % A',vacal,spmat_send(ineig) % a,icont*3)
                call memory_alloca(memor_cou,'PERMS % L'      ,vacal,perms(ineig) % l,kcont)
                nentr_send(ineig) = icont*3
                npoin_send(ineig) = kcont

                icont=1
                do jcont = 1, nenti
                   if(permm(jcont)>0) then 
                      permm(jcont) = icont
                      perms(ineig) % l(icont) = jcont
                      icont = icont + 1
                   end if
                end do

                icont=0
                if( coupling % kfl_source_value == VALUES_ON_NODES ) then
                   do jelem = 1,nsend
                      kelem = kelem + 1 
                      ielem = coupling % geome % lelem_source(kelem)
                      pnode = lnnod(ielem)
                      do inode = 1,pnode
                         val = coupling % geome % shapf(inode,kelem)
                         if(abs(val) > tolval) then
                            icont                                 = icont + 1
                            ipoin                                 = lnods(inode,ielem)
                            spmat_send(ineig) % a ((icont-1)*3+1) = real(jelem,rp)
                            spmat_send(ineig) % a ((icont-1)*3+2) = real(permm(ipoin) ,rp)
                            spmat_send(ineig) % a ((icont)*3)     = val
                         endif
                      end do
                   end do
                else
                   do jelem = 1,nsend
                      kelem                                 = kelem + 1 
                      ielem                                 = coupling % geome % lelem_source(kelem)
                      icont                                 = icont + 1
                      spmat_send(ineig) % a ((icont-1)*3+1) = real(jelem,rp)
                      spmat_send(ineig) % a ((icont-1)*3+2) = real(permm(ielem) ,rp)
                      spmat_send(ineig) % a ((icont)*3)     = 1.0_rp
                   end do
                end if
                   
             end do
          endif

       else if( itype == NEAREST_BOUNDARY_NODE .or. &
            &   itype == NEAREST_ELEMENT_NODE  .or. &
            &   itype == GLOBAL_NUMBERING      ) then

          !
          ! Nearest boundary/element node
          !
          if( coupling % geome % npoin_source > 0 ) then
             kpoin = 0
             do ineig = 1,nneig
                dom_i = coupling % commd % neights(ineig)          
                inise = coupling % commd % lsend_size(ineig) 
                finse = coupling % commd % lsend_size(ineig+1) - 1
                nsend = finse - inise + 1

                permm(1:npoin)=0
                kcont=0
                jcont=0
                do jpoin = 1,nsend
                   kpoin = kpoin + 1
                   jcont = jcont + 1
                   ipoin = coupling % geome % lpoin_source(kpoin)
                   if(.not. permm(ipoin) ==1) then
                      permm(ipoin) = 1
                      kcont = kcont +1
                   end if
                end do
                kpoin = kpoin - jcont

                call memory_alloca(memor_cou,'SPMAT_SEND % A',vacal,spmat_send(ineig) % a,nsend*3)
                call memory_alloca(memor_cou,'PERMS % L     ',vacal,perms(ineig) % l,kcont)
                nentr_send(ineig) = nsend*3

                icont=1
                do jcont = 1, npoin
                   if(permm(jcont)>0) then 
                      permm(jcont) = icont
                      perms(ineig) % l(icont) = jcont
                      icont = icont + 1
                   end if
                end do
                npoin_send(ineig) = icont-1

                do jpoin = 1,nsend
                   kpoin = kpoin + 1
                   ipoin = coupling % geome % lpoin_source(kpoin)
                   spmat_send(ineig) % a ((jpoin-1)*3+1) = real(jpoin,rp)
                   spmat_send(ineig) % a ((jpoin-1)*3+2) = real(permm(ipoin),rp)
                   spmat_send(ineig) % a (jpoin*3) = 1
                end do
             end do
          end if

       else if( itype == BOUNDARY_INTERPOLATION .or. itype == BOUNDARY_VECTOR_PROJECTION ) then

          !
          ! Boundary interpolation
          !
          if( coupling % geome % nboun_source > 0 ) then
             kboun = 0
             do ineig = 1,nneig
                dom_i = coupling % commd % neights(ineig)          
                inise = coupling % commd % lsend_size(ineig) 
                finse = coupling % commd % lsend_size(ineig+1) - 1
                nsend = finse - inise + 1
                icont = 0
                jcont = 0
                kcont = 0 

                permm(1:nenti)=0
                if( coupling % kfl_source_value == VALUES_ON_NODES) then
                   do jboun =1, nsend
                      kboun = kboun + 1
                      jcont = jcont + 1
                      iboun = coupling % geome % lboun_source(kboun)
                      pnodb = lnnob_cou(iboun)
                      do inodb = 1,pnodb
                         val = coupling % geome % shapf(inodb,kboun)
                         if(abs(val) > tolval) then
                            icont = icont + 1
                            ipoin = lnodb_cou(inodb,iboun)
                            if(.not. permm(ipoin) ==1) then
                               permm(ipoin) = 1
                               kcont = kcont +1
                            end if
                         end if
                      end do
                   end do
                else
                   do jboun =1, nsend
                      kboun = kboun + 1
                      jcont = jcont + 1
                      icont = icont + 1
                      iboun = coupling % geome % lboun_source(kboun)
                      if( .not. permm(iboun) == 1 ) then
                         permm(iboun) = 1
                         kcont = kcont +1
                      end if
                   end do
                endif
                kboun = kboun - jcont

                call memory_alloca(memor_cou,'SPMAT_SEND % A',vacal,spmat_send(ineig) % a,icont*3)
                call memory_alloca(memor_cou,'PERMS % L'     ,vacal,perms(ineig) % l,kcont)
                nentr_send(ineig) = icont*3

                icont=1
                do jcont = 1, nenti
                   if(permm(jcont)>0) then 
                      permm(jcont) = icont
                      perms(ineig) % l(icont) = jcont
                      icont = icont + 1
                   end if
                end do
                npoin_send(ineig) = icont-1

                icont = 0
                if( coupling % kfl_source_value == VALUES_ON_NODES) then
                   do jboun = 1,nsend
                      kboun = kboun + 1
                      iboun = coupling % geome % lboun_source(kboun)
                      pnodb = lnnob_cou(iboun)
                      do inodb = 1,pnodb
                         val = coupling % geome % shapf(inodb,kboun)
                         if(abs(val) > tolval) then
                            icont                                 = icont + 1
                            ipoin                                 = lnodb_cou(inodb,iboun)
                            spmat_send(ineig) % a ((icont-1)*3+1) = real(jboun,rp)
                            spmat_send(ineig) % a ((icont-1)*3+2) = real(permm(ipoin),rp)
                            spmat_send(ineig) % a (icont*3)       = val
                         endif
                      end do
                   end do
                else
                   do jboun = 1,nsend
                      kboun                                 = kboun + 1
                      iboun                                 = coupling % geome % lboun_source(kboun)
                      icont                                 = icont + 1
                      spmat_send(ineig) % a ((icont-1)*3+1) = real(jboun,rp)
                      spmat_send(ineig) % a ((icont-1)*3+2) = real(permm(iboun),rp)
                      spmat_send(ineig) % a (icont*3)       = 1.0_rp
                   end do
                endif
             end do

          end if

       else if( itype == STRESS_PROJECTION ) then

          !
          ! Projection
          !
          !
          ! XX_IN                   XX_IN
          !  ||                       ||
          !  \/      XX_SEND          \/
          !  o---------x--------------o     <= IBOUN
          !         
          !  <------------------------>
          !          SIZE_IBOUN
          !
          !           PNODB
          !           +---
          !           \
          !  XX_SEND= /    XX_IN(INODB) * N_INODB(x) / MASS(IPOIN)
          !           +---
          !           INODB=1
          !
          !  N_INODB(x) is the shape function at the wet point position in boundary IBOUN
          !  N_INODB(x) = COUPLING % GEOME % SHAPF(INODB,KBOUN)
          !  KBOUN      = Local nummbering of global boundary IBOUN
          !  SIZE_IBOUN = Surface of IBOUN
          !
          if( coupling % geome % nboun_source > 0 ) then
             kboun = 0
             kgaub = 0
             do ineig = 1,nneig
                dom_i = coupling % commd % neights(ineig)          
                inise = coupling % commd % lsend_size(ineig) 
                finse = coupling % commd % lsend_size(ineig+1) - 1
                nsend = finse - inise + 1

                icont = 0
                jcont = 0
                kcont = 0

                permm(1:npoin)=0
                do jboun =1, nsend
                   kboun = kboun + 1
                   jcont = jcont + 1
                   iboun = coupling % geome % lboun_source(kboun)
                   pnodb = lnnob_cou(iboun)
                   do inodb = 1,pnodb
                      ipoin = lnodb_cou(inodb,iboun)
                      val = coupling % geome % shapf(inodb,kboun)/coupling % geome % vmasb(ipoin)
                      if(abs(val) > tolval) then
                         icont = icont + 1
                         ipoin = lnodb_cou(inodb,iboun)
                         if(.not. permm(ipoin) ==1) then
                            permm(ipoin) = 1
                            kcont = kcont +1
                         end if
                      endif
                   end do
                end do
                kboun = kboun - jcont

                call memory_alloca(memor_cou,'SPMAT_SEND % A',vacal,spmat_send(ineig) % a,icont*3)
                call memory_alloca(memor_cou,'PERMS % L'     ,vacal,perms(ineig) % l,kcont)
                nentr_send(ineig) = icont*3

                icont=1
                do jcont = 1, npoin
                   if(permm(jcont)>0) then 
                      permm(jcont) = icont
                      perms(ineig) % l(icont) = jcont
                      icont = icont + 1
                   end if
                end do
                npoin_send(ineig) = icont-1

                icont = 0
                do jboun = 1,nsend
                   kboun = kboun + 1
                   kgaub = kgaub + 1
                   iboun = coupling % geome % lboun_source(kboun)
                   pnodb = lnnob_cou(iboun)
                   do inodb = 1,pnodb
                      ipoin = lnodb_cou(inodb,iboun)
                      val = coupling % geome % shapf(inodb,kboun)/coupling % geome % vmasb(ipoin)
                      if(abs(val) > tolval) then
                         icont = icont + 1
                         spmat_send(ineig) % a ((icont-1)*3+1) = real(jboun,rp)
                         spmat_send(ineig) % a ((icont-1)*3+2) = real(permm(ipoin),rp)
                         spmat_send(ineig) % a (icont*3) = val
                      endif

                   end do
                end do
             end do
          end if

       else if( itype == PROJECTION ) then

          !
          ! Projection
          !
          if( coupling % geome % nboun_source > 0 ) then
             kboun = 0
             kgaub = 0
             do ineig = 1,nneig
                dom_i = coupling % commd % neights(ineig)          
                inise = coupling % commd % lsend_size(ineig) 
                finse = coupling % commd % lsend_size(ineig+1) - 1
                nsend = finse - inise + 1

                icont = 0
                jcont = 0
                kcont = 0

                permm(1:npoin)=0
                do jboun =1, nsend
                   kboun = kboun + 1
                   jcont = jcont + 1
                   iboun = coupling % geome % lboun_source(kboun)
                   pnodb = lnnob_cou(iboun)
                   do inodb = 1,pnodb
                      val = coupling % geome % shapf(inodb,kboun)
                      if(abs(val) > tolval) then 
                         icont = icont + 1
                         ipoin = lnodb_cou(inodb,iboun)
                         if(.not. permm(ipoin) ==1) then
                            permm(ipoin) = 1
                            kcont = kcont +1
                         end if
                      endif
                   end do
                end do
                kboun = kboun - jcont

                call memory_alloca(memor_cou,'SPMAT_SEND % A',vacal,spmat_send(ineig) % a,icont*3)
                call memory_alloca(memor_cou,'PERMS % L'     ,vacal,perms(ineig) % l,kcont)
                nentr_send(ineig) = icont*3

                icont=1
                do jcont = 1, npoin
                   if(permm(jcont)>0) then 
                      permm(jcont) = icont
                      perms(ineig) % l(icont) = jcont
                      icont = icont + 1
                   end if
                end do
                npoin_send(ineig) = icont-1

                icont = 0
                do jboun = 1,nsend
                   kboun = kboun + 1
                   kgaub = kgaub + 1
                   iboun = coupling % geome % lboun_source(kboun)
                   pnodb = lnnob_cou(iboun)
                   do inodb = 1,pnodb
                      val = coupling % geome % shapf(inodb,kboun)
                      if(abs(val) > tolval) then
                         icont = icont + 1
                         ipoin = lnodb_cou(inodb,iboun)
                         spmat_send(ineig) % a ((icont-1)*3+1) = real(jboun,rp)
                         spmat_send(ineig) % a ((icont-1)*3+2) = real(permm(ipoin),rp)
                         spmat_send(ineig) % a (icont*3) = val
                      endif
                   end do
                end do
             end do
          end if

       end if
    end if

    !--------------------------------------------------------------------
    !
    ! Send transmission matrix to target  
    !
    !--------------------------------------------------------------------
    !
    ! Send number of entries of matrix
    !
    call cputim(time1)
    call PAR_START_NON_BLOCKING_COMM(1_ip,nneig)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,nneig
       dom_i = coupling % commd % neights(ineig)          
       call PAR_SEND_RECEIVE(1_ip,1_ip,nentr_send(ineig:ineig),nentr_recv(ineig:ineig),'IN CURRENT COUPLING',dom_i,'NON BLOCKING')
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Send number of points required 
    !
    call PAR_START_NON_BLOCKING_COMM(1_ip,nneig)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,nneig
       dom_i = coupling % commd % neights(ineig)          
       call PAR_SEND_RECEIVE(1_ip,1_ip,npoin_send(ineig:ineig),npoin_recv(ineig:ineig),'IN CURRENT COUPLING',dom_i,'NON BLOCKING')
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    ! note: two previous comunications can be joined in one
    !
    ! Send matrix entries
    !
    call PAR_START_NON_BLOCKING_COMM(1_ip,nneig)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,nneig
       dom_i = coupling % commd % neights(ineig)
       nsend = nentr_send(ineig)
       nrecv = nentr_recv(ineig)
       if( nsend > 0 .and. nrecv > 0 ) then
          call memory_alloca(memor_cou,'SPMAT_RECV % A',vacal, spmat_recv(ineig) % a,nrecv)
          call PAR_SEND_RECEIVE(nsend,nrecv,spmat_send(ineig) % a(1:nsend),spmat_recv(ineig) % a(1:nrecv),'IN CURRENT COUPLING',dom_i,'NON BLOCKING')
       else
          if(nsend>0) then
             call PAR_SEND_RECEIVE(nsend,0_ip,spmat_send(ineig) % a(1:nsend),dummr(1:1),'IN CURRENT COUPLING',dom_i,'NON BLOCKING')
          endif
          if(nrecv>0) then
             call memory_alloca(memor_cou,'SPMAT_RECV % A',vacal, spmat_recv(ineig) % a,nrecv)
             call PAR_SEND_RECEIVE(0_ip,nrecv,dummr(1:1),spmat_recv(ineig) % a(1:nrecv),'IN CURRENT COUPLING',dom_i,'NON BLOCKING')
          end if
       end if
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)

    !--------------------------------------------------------------------
    !
    ! Store matrices in coupling % ltransmat_target
    !
    !--------------------------------------------------------------------

    if(I_AM_TARGET) then
       !
       ! Eval permutation of row ordinates
       !
       nullify(inv_status)
       call memory_alloca(memor_cou,'INV_STATUS',vacal,inv_status,coupling % commd % lrecv_dim)
       if( itype /= PROJECTION .and. itype /= STRESS_PROJECTION) then
          do icont = 1, coupling % wet   % number_wet_points
             if( coupling % geome % status(icont) /= 0 ) then
                inv_status(coupling % geome % status(icont)) = icont
             end if
          end do
       else
          do icont = 1, coupling % wet % number_wet_points
             inv_status(icont) = icont
          end do
       endif
       !
       ! Introduce entries to spmat structure
       !
       call memory_alloca(memor_cou,'LTRANSMAT_TARGET',vacal, coupling % ltransmat_target,nneig) ! SPMAT
       cont_wetp = 0
       cont_recv = 0
       icont = 1
       do ineig = 1,nneig
          dom_i = coupling % commd % neights(ineig)   
          nentr_ineig = int(nentr_recv(ineig)/3,ip)
          cont_wetp = coupling % commd % lrecv_size(ineig)-1

          !call coupling % ltransmat_target(ineig) % alloca(1_ip,1_ip,nentr_ineig,memor_cou)
          call memory_alloca(memor_cou,'LTRANSMAT_TARGET',vacal,coupling % ltransmat_target(ineig),1_ip,nentr_ineig) ! SPMAT
      
          do jcont = 1, nentr_ineig
             coupling % ltransmat_target(ineig) % iA (jcont) = inv_status(int(spmat_recv(ineig) % a((jcont-1)*3+1), ip)+cont_wetp)
             coupling % ltransmat_target(ineig) % jA (jcont) = int(spmat_recv(ineig) % a((jcont-1)*3+2), ip) + cont_recv
             coupling % ltransmat_target(ineig) % vA (1,1,jcont) = spmat_recv(ineig) % a((jcont)*3)
             icont = icont + 1
          end do
          cont_recv = cont_recv + npoin_recv(ineig)
       end do

       do ineig = 1,nneig
          coupling % ltransmat_target(ineig) % nrows = coupling % wet % number_wet_points
          coupling % ltransmat_target(ineig) % ncols = cont_recv
          coupling % ltransmat_target(ineig) % ndof1 = 1
          coupling % ltransmat_target(ineig) % ndof2 = 1
       end do
       !
       ! For projection schemes, integrate gaussian points: multiply received matrix by projection matrix
       !
       if( itype == PROJECTION .or. itype == STRESS_PROJECTION) then
          if(coupling % wet % npoin_wet > 0) then 
             call COU_GENERATE_PROJECTION_MATRIX(coupling,M_proj)
             do ineig = 1,nneig
                if(associated(coupling % ltransmat_target(ineig) % iA)) then
                   call matrix_COO_SPGEMM(M_proj,coupling % ltransmat_target(ineig),M_aux)                   
                   !call coupling % ltransmat_target(ineig) % deallo()
                   !call coupling % ltransmat_target(ineig) % alloca(1_ip,1_ip,int(size(M_aux % iA),ip))
                   call memory_deallo(memor_cou,'LTRANSMAT_TARGET',vacal,coupling % ltransmat_target(ineig)) ! SPMAT
                   call memory_alloca(memor_cou,'LTRANSMAT_TARGET',vacal,coupling % ltransmat_target(ineig), 1_ip, int(size(M_aux % iA),ip)) ! SPMAT
                   coupling % ltransmat_target(ineig) % ndof1 = M_aux % ndof1
                   coupling % ltransmat_target(ineig) % ndof2 = M_aux % ndof2
                   coupling % ltransmat_target(ineig) % nrows = M_aux % nrows
                   coupling % ltransmat_target(ineig) % ncols = M_aux % ncols
                   coupling % ltransmat_target(ineig) % iA(:) = M_aux % iA(:)
                   coupling % ltransmat_target(ineig) % jA(:) = M_aux % jA(:)
                   coupling % ltransmat_target(ineig) % vA(1,1,:) = M_aux % vA(1,1,:)
                   !call coupling % M_aux % deallo()                   
                   call memory_deallo(memor_cou,'M_AUX',vacal,M_aux) 
                endif
             end do
             call memory_deallo(memor_cou,'M_PROJ',vacal,M_proj)
          endif
       endif
       !
       ! Agregatte (join) splited entries in transmission matrix
       !
       call cputim(time1)
       do ineig = 1,nneig
          if(associated(coupling % ltransmat_target(ineig) % iA)) then
             call matrix_COO_aggregate(coupling % ltransmat_target(ineig))
          endif
       enddo
       if(associated(inv_status)) &
            call memory_deallo(memor_cou,'INV_STATUS',vacal,inv_status)
    endif

    !--------------------------------------------------------------------
    !
    ! Generate new communication scheme   
    !
    !--------------------------------------------------------------------
    call cputim(time1)

    if( nneig > 0) then
       if(I_AM_SOURCE) then 
          coupling % commd % lsend_size(1) = 1 
          do ineig = 2,nneig + 1
             coupling % commd % lsend_size(ineig) = coupling % commd % lsend_size(ineig-1) + npoin_send(ineig-1)
          end do
          coupling % commd % lsend_dim = coupling % commd % lsend_size(nneig + 1) - 1

          call memory_alloca(memor_cou,'coupling % commd % lsend_perm',vacal,coupling % commd % lsend_perm,&
               coupling % commd % lsend_dim)
          do ineig = 1, nneig
             do icont =1, npoin_send(ineig)
                coupling % commd % lsend_perm(coupling % commd % lsend_size(ineig)+ icont -1) = perms(ineig) % l(icont)
             end do
          end do
       endif
       if(I_AM_TARGET) then
          coupling % commd % lrecv_size(1) = 1
          do ineig = 2,nneig + 1
             coupling % commd % lrecv_size(ineig) = coupling % commd % lrecv_size(ineig-1) + npoin_recv(ineig-1)
          end do
          coupling % commd % lrecv_dim = coupling % commd % lrecv_size(nneig + 1) - 1

       endif
    endif

    !--------------------------------------------------------------------
    !
    ! Finish 
    !
    !--------------------------------------------------------------------
    !
    ! Deallocate pointers
    !
    call memory_deallo(memor_cou,'nentr_send',vacal,nentr_send)
    call memory_deallo(memor_cou,'nentr_recv',vacal,nentr_recv)
    call memory_deallo(memor_cou,'npoin_send',vacal,npoin_send)
    call memory_deallo(memor_cou,'npoin_recv',vacal,npoin_recv)
    call memory_deallo(memor_cou,'spmat_send',vacal,spmat_send)
    call memory_deallo(memor_cou,'spmat_recv',vacal,spmat_recv)
    call memory_deallo(memor_cou,'perms'     ,vacal,perms)
    call memory_deallo(memor_cou,'permm'     ,vacal,permm)
    !
    ! Recover communicator
    !

    PAR_COMM_CURRENT = PAR_COMM_SAVE
    !
    ! Timer
    !
    call cputim(time2) 
    coupling % cputim(7) = coupling % cputim(7) + time2 - time1

  end subroutine COU_GENERATE_LOCAL_TRANSMISSION_MATRICES

 !----------------------------------------------------------------------
 !
 !> @author  Ricard Borrell
 !> @date    27/03/2018
 !> @brief   Generate transposed transmission matrices  
 !> @details For the TRANSPOSE_MIRROR coupling obtain the transmat from its mirror
 !>
 !----------------------------------------------------------------------
  subroutine COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES(coupling)

    implicit none
    type(typ_color_coupling), intent(inout)        :: coupling
    real(rp)                                       :: time1, time2
    MY_MPI_COMM                                    :: PAR_COMM_SAVE
    integer(ip)                                    :: PAR_CURRENT_SIZE
    integer(ip)                                    :: PAR_CURRENT_RANK
    logical(lg)                                    :: I_AM_SOURCE
    logical(lg)                                    :: I_AM_TARGET
    integer(ip)                                    :: micou
    integer(ip)                                    :: nneig, ineig, ipart, ipoin 
    integer(ip)                                    :: inentr, ientr, isend
    integer(ip)                                    :: nsend, nrecv, npose
    integer(ip)                                    :: cont_recv
    integer(ip)                                    :: npoin_wet, iwetp
    integer(ip)                                    :: isize
    type(r1p),   pointer                           :: spmat_send(:)
    type(r1p),   pointer                           :: spmat_recv(:)
    integer(ip), pointer                           :: nentr_send(:)
    integer(ip), pointer                           :: nentr_recv(:)
    integer(ip), pointer                           :: npoin_send(:)
    integer(ip), pointer                           :: npoin_recv(:)
    integer(ip), pointer                           :: permm(:)
    type(i1p),   pointer                           :: perms(:)
    integer(ip), pointer                           :: inv_lpwet(:)
    real(rp),    pointer                           :: multip_source(:)
    real(rp)                                       :: dummr(1)

    character(100), PARAMETER :: vacal = "COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES"
    !
    ! Nullify
    !
    nullify(spmat_send)
    nullify(spmat_recv)
    nullify(nentr_send)
    nullify(nentr_recv)
    nullify(npoin_send)
    nullify(npoin_recv)
    nullify(permm)   
    nullify(perms)   
    nullify(inv_lpwet)
    nullify(multip_source)
    !
    ! Target and source colors
    !
    color_target = coupling % color_target
    color_source = coupling % color_source
    I_AM_SOURCE = .false.
    I_AM_TARGET = .false.
    if(I_AM_IN_COLOR(color_source)) then 
       I_AM_SOURCE = .true.
    end if
    if(I_AM_IN_COLOR(color_target)) then
       I_AM_TARGET = .true.
    end if
    !
    ! Communicator
    !
    PAR_COMM_SAVE    = PAR_COMM_CURRENT 
    PAR_COMM_CURRENT = PAR_COMM_COLOR(color_target,color_source) 
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    !
    ! Transpose and send_recv
    !  
    if( coupling % itype == TRANSPOSE_MIRROR) then

       micou = coupling % mirror_coupling
       coupling % commd % nneig = coupling_type(micou) % commd % nneig
       nneig = coupling % commd % nneig

       call cputim(time1)
       nullify(spmat_send,spmat_recv)
       nullify(nentr_recv,nentr_send)
       nullify(npoin_send,npoin_recv)
       call memory_alloca(memor_cou,'coupling % commd % neights'   ,vacal,coupling % commd % neights   ,nneig)
       call memory_alloca(memor_cou,'spmat_send'                   ,vacal,spmat_send                   ,nneig)
       call memory_alloca(memor_cou,'spmat_recv'                   ,vacal,spmat_recv                   ,nneig)
       call memory_alloca(memor_cou,'nentr_recv'                   ,vacal,nentr_recv                   ,nneig)
       call memory_alloca(memor_cou,'nentr_send'                   ,vacal,nentr_send                   ,nneig)
       call memory_alloca(memor_cou,'npoin_send'                   ,vacal,npoin_send                   ,nneig)
       call memory_alloca(memor_cou,'npoin_recv'                   ,vacal,npoin_recv                   ,nneig)
       call memory_alloca(memor_cou,'coupling % commd % lsend_size',vacal,coupling % commd % lsend_size,nneig+1)
       call memory_alloca(memor_cou,'coupling % commd % lrecv_size',vacal,coupling % commd % lrecv_size,nneig+1)
       if( nneig > 0 ) then
          coupling % commd % neights               = coupling_type(micou) % commd % neights
          nentr_recv(1:nneig)                      = 0
          nentr_send(1:nneig)                      = 0
          npoin_send(1:nneig)                      = 0
          npoin_recv(1:nneig)                      = 0
          coupling % commd % lsend_size(1:nneig+1) = 0
          coupling % commd % lrecv_size(1:nneig+1) = 0
          coupling % commd % lsend_dim             = 0
          coupling % commd % lrecv_dim             = 0
       end if

       !--------------------------------------------------------------------
       !
       ! 1.SOURCE: Transpose and eval lsend_size, lsend_dim and lsend_perm 
       !
       !--------------------------------------------------------------------
       
       if( I_AM_SOURCE .and. INOTMASTER) then

          coupling % geome % npoin_source = coupling_type(micou) % wet % npoin_wet
          call memory_alloca(memor_cou,'coupling % geome % lpoin_source',vacal,coupling % geome % lpoin_source,coupling % geome % npoin_source)  
          if(coupling % geome % npoin_source>0) coupling % geome % lpoin_source = coupling_type(micou) % wet % lpoin_wet

          call memory_alloca(memor_cou,'perms',vacal,perms,nneig)

          isize=max(npoin,1_ip)
          call memory_alloca(memor_cou,'multip_source',vacal,multip_source,isize)
          multip_source(1:isize) = 1.0_rp


          if( coupling_type(micou) % itype /= PROJECTION .and. &
              coupling_type(micou) % itype /= STRESS_PROJECTION ) then
          
               call PAR_INTERFACE_NODE_EXCHANGE(multip_source,'SUM','IN CURRENT SOURCE COLOR')
               multip_source = 1.0_rp / multip_source

         endif

          do ineig = 1, nneig
             if(associated(coupling_type(micou) % ltransmat_target(ineig) % iA)) then
                inentr = size(coupling_type(micou) % ltransmat_target(ineig) % iA)
                nentr_send( ineig ) = inentr
                npoin_wet = coupling_type(micou) % wet % npoin_wet 
                call memory_alloca(memor_cou,'spmat_send % a',vacal,spmat_send(ineig) % a,inentr*3_ip)
                call memory_alloca(memor_cou,'permm',vacal,permm,npoin_wet)
                permm(1:npoin_wet) = 0_ip
                npose = 0_ip
                do ientr = 1,inentr
                   if( .not. permm(coupling_type(micou) % ltransmat_target(ineig) % iA(ientr)) == 1_ip) then
                      permm(coupling_type(micou) % ltransmat_target(ineig) % iA(ientr)) = 1_ip   
                      npose = npose + 1_ip
                   end if
                end do
                call memory_alloca(memor_cou,'perms % l',vacal,perms(ineig) % l,npose)
                do iwetp = 1,npoin_wet
                   if(permm(iwetp) == 1) then
                      npoin_send(ineig) = npoin_send(ineig) + 1_ip
                      permm(iwetp) = npoin_send(ineig)
                      perms(ineig) % l(npoin_send(ineig)) = coupling_type(micou) % wet % lpoin_wet(iwetp) 
                   endif
                end do
                do ientr = 1,inentr
                   ipoin = coupling_type(micou) % wet % lpoin_wet(coupling_type(micou) % ltransmat_target(ineig) % iA (ientr))                      
                   spmat_send(ineig) % a ((ientr-1)*3+1) = real(coupling_type(micou) % ltransmat_target(ineig) % jA (ientr) - &
                        coupling_type(micou) % commd % lrecv_size(ineig) + 1_ip,rp)
                   spmat_send(ineig) % a ((ientr-1)*3+2) = real(permm(coupling_type(micou) % ltransmat_target(ineig) % iA (ientr)),rp)
                   spmat_send(ineig) % a (ientr*3)       = coupling_type(micou) % ltransmat_target(ineig) % vA (1,1,ientr) * multip_source(ipoin)
                end do
                call memory_deallo(memor_cou,'permm',vacal,permm) 
             endif
          end do

          call memory_deallo(memor_cou,'multip_source',vacal,multip_source)

          if( nneig > 0) then
             coupling % commd % lsend_size(1) = 1 
             do ineig = 2,nneig + 1
                coupling % commd % lsend_size(ineig) = coupling % commd % lsend_size(ineig-1) + npoin_send(ineig-1)
             end do
             coupling % commd % lsend_dim = coupling % commd % lsend_size(nneig + 1) - 1

             call memory_alloca(memor_cou,'coupling % commd % lsend_perm',vacal,coupling % commd % lsend_perm,&
                  coupling % commd % lsend_dim)
             do ineig = 1, nneig
                do isend =1, npoin_send(ineig)
                   coupling % commd % lsend_perm(coupling % commd % lsend_size(ineig)-1+isend ) = perms(ineig) % l(isend)
                end do
             end do
          endif
          call memory_deallo(memor_cou,'perms',vacal,perms)
       end if
       !--------------------------------------------------------------------
       !
       ! 2.BOTH: Send/Recv
       !
       !--------------------------------------------------------------------
       !
       ! Send_recv #entries
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,nneig)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,nneig
          ipart = coupling % commd % neights(ineig)          
          call PAR_SEND_RECEIVE(1_ip,1_ip,nentr_send(ineig:ineig),nentr_recv(ineig:ineig),&
               'IN CURRENT COUPLING',ipart,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !
       ! Send_recv #(points to be transmited) 
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,nneig)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,nneig
          ipart = coupling % commd % neights(ineig)          
          call PAR_SEND_RECEIVE(1_ip,1_ip,npoin_send(ineig:ineig),npoin_recv(ineig:ineig),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       ! note: two previous comunications can be joined in one
       !
       ! Send_recv matrix entries
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,nneig)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,nneig
          ipart = coupling % commd % neights(ineig)
          nsend = nentr_send(ineig)*3
          nrecv = nentr_recv(ineig)*3
          !if( nsend > 0 .and. nrecv > 0 ) then
          !   call memory_alloca(memor_cou,'spmat_recv',vacal, spmat_recv(ineig) % a,nrecv)
          !   call PAR_SEND_RECEIVE(nsend,nrecv,spmat_send(ineig) % a(1:nsend),spmat_recv(ineig) % a(1:nrecv),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
          !else
          if(nsend>0) then
             call PAR_SEND_RECEIVE(nsend,0_ip,spmat_send(ineig) % a(1:nsend),dummr(1:1),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
          endif
          if(nrecv>0) then
             call memory_alloca(memor_cou,'spmat_recv % a',vacal, spmat_recv(ineig) % a,nrecv)
             call PAR_SEND_RECEIVE(0_ip,nrecv,dummr(1:1),spmat_recv(ineig) % a(1:nrecv),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
          end if
          !end if
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !--------------------------------------------------------------------
       !
       ! 3.TARGET: Store matrix and eval lsend_size and lsend dim
       !
       !--------------------------------------------------------------------
       if(I_AM_TARGET .and. INOTMASTER) then

          cont_recv = 0
          call memory_alloca(memor_cou,'ltransmat_target',vacal,coupling % ltransmat_target,nneig)
          call memory_alloca(memor_cou,'inv_lpwet'       ,vacal,inv_lpwet                  ,npoin)
          npoin_wet = coupling % wet % npoin_wet 

          inv_lpwet(1:npoin) = 0
          do iwetp=1,npoin_wet
             inv_lpwet(coupling % wet % lpoin_wet(iwetp)) = iwetp
          end do

          do ineig = 1,nneig
             ipart = coupling % commd % neights(ineig)   
             inentr = int(nentr_recv(ineig),ip)
             if(inentr > 0) then
                call memory_alloca(memor_cou,'ltransmat_target',vacal,coupling % ltransmat_target(ineig),1_ip,inentr)
                do ientr = 1, inentr
                   coupling % ltransmat_target(ineig) % iA (ientr) = inv_lpwet( coupling_type(micou) % commd % &
                        lsend_perm(int(spmat_recv(ineig) % a((ientr-1)*3+1)) + coupling_type(micou) % commd % lsend_size(ineig) -1 )) 
                   coupling % ltransmat_target(ineig) % jA (ientr) = int(spmat_recv(ineig) % a((ientr-1)*3+2), ip) + cont_recv
                   coupling % ltransmat_target(ineig) % vA (1,1,ientr) = spmat_recv(ineig) % a((ientr)*3)
                end do
                cont_recv = cont_recv + npoin_recv(ineig)
             endif
          end do
          call memory_deallo(memor_cou,'inv_lpwet',vacal,inv_lpwet)

          do ineig = 1,nneig
             inentr = int(nentr_recv(ineig),ip)
             if(inentr > 0) then
                coupling % ltransmat_target(ineig) % nrows = coupling % wet % npoin_wet
                coupling % ltransmat_target(ineig) % ncols = cont_recv
                coupling % ltransmat_target(ineig) % ndof1 = 1
                coupling % ltransmat_target(ineig) % ndof2 = 1
             end if
          end do

          if( nneig > 0) then
             coupling % commd % lrecv_size(1) = 1
             do ineig = 2,nneig + 1
                coupling % commd % lrecv_size(ineig) = coupling % commd % lrecv_size(ineig-1) + npoin_recv(ineig-1)
             end do
             coupling % commd % lrecv_dim = coupling % commd % lrecv_size(nneig + 1) - 1
          endif
       end if
       !
       ! Deallocate arrays
       !
       call memory_deallo(memor_cou,'nentr_send',vacal,nentr_send)
       call memory_deallo(memor_cou,'nentr_recv',vacal,nentr_recv)
       call memory_deallo(memor_cou,'npoin_send',vacal,npoin_send)
       call memory_deallo(memor_cou,'npoin_recv',vacal,npoin_recv)
       call memory_deallo(memor_cou,'spmat_send',vacal,spmat_send)
       call memory_deallo(memor_cou,'spmat_recv',vacal,spmat_recv)
       call cputim(time2) 
       coupling % cputim(7) = coupling % cputim(7) + time2 - time1
    end if
    !
    ! Recover communicator
    !
    PAR_COMM_CURRENT = PAR_COMM_SAVE

  end subroutine COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    27/03/2018
  !> @brief   Parallelize transmission matrices 
  !> @details Complete the rows of the interface nodes in order to avoid 
  !>          one communication episode in the coupling process
  !>
  !----------------------------------------------------------------------
  subroutine COU_PARALLELIZE_TRANSMISSION_MATRICES(coupling)

    implicit none
    type(typ_color_coupling), intent(inout)        :: coupling

    logical(lg)                                    :: I_AM_SOURCE,I_AM_TARGET,auxlg
    real(rp)                                       :: time1, time2
    real(rp)                                       :: dummr(1), weight
    integer(ip)                                    :: dummi(1) 
    integer(ip)                                    :: PAR_COUP_RANK, PAR_COUP_SIZE
    MY_MPI_COMM                                    :: PAR_COMM_COUP,PAR_COMM_TARG
    integer(ip)                                    :: PAR_TARG_RANK, PAR_TARG_SIZE
    integer(ip)                                    :: aux,aux1(1), permval
    integer(ip)                                    :: ipart,jpart,icol,irow,iwetp,ineig,ientr,jcol
    integer(ip)                                    :: icont,jcont,kcont,lcont,ibuff,ipoin
    integer(ip)                                    :: npoin_wet,nneig_tar,nsend,nrecv,nentr,nentr_prev
    integer(ip)                                    :: nneig_prev,nneig_new_src,nneig_new_tar
    integer(ip)                                    :: dim_prev_tar,dim_prev_src,dim_new_src,dim_new_tar 
    integer(ip)                                    :: bound_dim,incr 
    integer(ip)                                    :: max_recv_size,max_send_size
    integer(ip)                                    :: sum_entr, sum_send
    integer(ip)                                    :: sum_rep, max_rep
    real(rp), pointer                              :: multip(:)
    integer(ip), pointer                           :: rank_targ2coup(:)
    integer(ip), pointer                           :: lpoin_wet(:)
    integer(ip), pointer                           :: bound_perm(:),bound_size(:)
    integer(ip), pointer                           :: lnent(:)
    integer(ip), pointer                           :: lnent_send(:), lnent_recv(:)
    integer(ip), pointer                           :: lncol_demand_tar(:),lncol_send_src(:)
    integer(ip), pointer                           :: neights_new_tar(:),neights_new_src(:)
    integer(ip), pointer                           :: lsend_size_new(:),lrecv_size_new(:)
    integer(ip), pointer                           :: perm_new_src(:)
    integer(ip), pointer                           :: order_neights(:),llcol_recv_ord(:,:),llcol_send_ord(:,:)
    integer(ip), pointer                           :: lnent_add(:),lncol_add(:),lnsend_add(:)
    integer(ip), pointer                           :: inv_lpwet(:)
    integer(ip), pointer                           :: norep_a(:,:)                                   
    integer(ip), pointer                           :: norep_b(:,:)                                   
    integer(ip), pointer                           :: norep_c(:)                                   
    type(i1p), pointer                             :: llcol_demand_tar(:),llcol_send_src(:)
    type(i1p), pointer                             :: llsend_add(:)
    type(r1p), pointer                             :: llent_send(:),llent_recv(:),llent(:),llent_add(:)
    type(spmat), pointer                           :: M(:),M_new(:)
    type(comm_data_par), pointer                   :: comm_target

    character(100), PARAMETER :: vacal = "COU_PARALLELIZE_TRANSMISSION_MATRICES"

    !
    ! Timer
    !
    call cputim(time1)
    !
    ! Target and source colors
    !
    color_target = coupling % color_target
    color_source = coupling % color_source
    I_AM_SOURCE = .false.
    I_AM_TARGET = .false.
    if(I_AM_IN_COLOR(color_source))  I_AM_SOURCE = .true.
    if(I_AM_IN_COLOR(color_target))  I_AM_TARGET = .true.
    !
    ! Communicators
    !
    PAR_COMM_COUP = PAR_COMM_COLOR(color_target,color_source) 
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_COUP,PAR_COUP_RANK,PAR_COUP_SIZE)
    if(I_AM_TARGET) then
       comm_target   => PAR_COMM_COLOR_ARRAY(color_target)
       PAR_COMM_TARG =  comm_target % PAR_COMM_WORLD
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TARG,PAR_TARG_RANK,PAR_TARG_SIZE)
    endif
    !
    ! Global pointers
    !
    nullify(lncol_demand_tar,lncol_send_src,order_neights,llcol_demand_tar,M_new,rank_targ2coup)
    call memory_alloca(memor_cou,'order_neights',vacal,order_neights,PAR_COUP_SIZE)
    call memory_alloca(memor_cou,'lncol_send_src',vacal,lncol_send_src,PAR_COUP_SIZE)
    call memory_alloca(memor_cou,'lncol_demand_tar',vacal,lncol_demand_tar,PAR_COUP_SIZE)
    call memory_alloca(memor_cou,'llcol_demand_tar',vacal,llcol_demand_tar,PAR_COUP_SIZE)
    lncol_demand_tar(:) = 0 
    lncol_send_src(:) = 0
    nneig_prev = coupling % commd % nneig
    nneig_new_tar = 0_ip
    nneig_new_src = 0_ip
    !
    ! Target operations
    !
    if(I_AM_TARGET) then
       !
       ! 1. 
       ! a) Fill-in rank_targ2coup (correspondence between rank of the parallel
       !    process in the target communicator, PAR_TARG_RANK, and in the coupling 
       !    communicator, PAR_COUP_RANK) 
       ! b) Eval max_recv_size (for all the target processes in the coupling)
       ! c) Eval multiplicity of nodes in target subgomain
       !
       nullify(multip)
       call memory_alloca(memor_cou,'rank_targ2coup',vacal,rank_targ2coup,PAR_TARG_SIZE)
       call PAR_ALLGATHER(PAR_COUP_RANK,rank_targ2coup,1_4,'IN CURRENT TARGET COLOR')
       call memory_alloca(memor_cou,'multip',vacal,multip,npoin)
       aux1(1) = 0
       do ineig = 1, coupling % commd % nneig
          aux = coupling % commd % lrecv_size(ineig+1) - coupling % commd % lrecv_size(ineig)
          if(aux > aux1(1)) aux1(1) = aux 
       end do
       call PAR_MAX(1_ip,aux1,'IN CURRENT TARGET COLOR')
       max_recv_size = aux1(1)
       if(npoin>0) multip(1:npoin) = 1.0_rp
       call PAR_INTERFACE_NODE_EXCHANGE(multip,'SUM','IN CURRENT TARGET COLOR')
    endif

    if(I_AM_TARGET .and. INOTMASTER) then
       !
       ! 2. Initialize reference pointers
       !
       lpoin_wet  => coupling % wet % lpoin_wet
       npoin_wet  =  coupling % wet % npoin_wet
       M          => coupling % ltransmat_target
       bound_size => comm_target % bound_size
       bound_perm => comm_target % bound_perm
       bound_dim  =  comm_target % bound_dim
       nneig_tar  =  comm_target % nneig
       !
       ! 3. Send/Recv number of entries to share 
       !
       ! 3.1 Fill in:
       !     a) lnent: number of entres of the transmission matrix corresponding 
       !               to each point (only wet points can have lnent > 0)
       !     b) llent: list of entries from the transmission matrix correspoding
       !               to each point
       !
       nullify(lnent,llent)
       call memory_alloca(memor_cou,'LNENT',vacal,lnent,npoin)
       call memory_alloca(memor_cou,'LLENT',vacal,llent,npoin)
       lnent(1:npoin) = 0

       if(associated(M)) then
          do ineig = 1, coupling % commd % nneig
             if(associated(M(ineig) % iA)) then
                do ientr = 1,size(M(ineig) % iA)
                   ipoin = lpoin_wet(M(ineig) % iA(ientr))
                   lnent(ipoin) = lnent(ipoin)+1
                end do
             end if
          end do
       end if
       do ipoin =1,npoin
          if(lnent(ipoin)>0) then
             call memory_alloca(memor_cou,'LLENT % A',vacal,llent(ipoin) % a,lnent(ipoin)*3)
             lnent(ipoin)=0
          endif
       end do
       if(associated(M)) then
          do ineig = 1, coupling % commd % nneig
             if(associated(M(ineig) % iA)) then
                do ientr = 1,size(M(ineig) % iA,kind=ip)
                   ipoin = lpoin_wet(M(ineig) % iA(ientr))
                   llent(ipoin) % a((lnent(ipoin))*3+1) = real(M(ineig) % jA(ientr) - coupling % commd % lrecv_size(ineig) + 1,rp) ! column
                   llent(ipoin) % a((lnent(ipoin))*3+2) = M(ineig) % vA(1,1,ientr)                                                 ! value
                   llent(ipoin) % a((lnent(ipoin))*3+3) = real(coupling % commd % neights(ineig),rp)                               ! proc
                   lnent(ipoin) = lnent(ipoin)+1
                end do
             end if
          end do
       end if
       !
       ! 3.2 Fill in:
       !     a) lnent_send: of dimention nneig_tar stores the number of entries
       !     to be sent to each neighbour process within the target subdomain
       !
       !     b) lnent_recv: of dimention nneig_tar stores the number of entries
       !     to be sent to each neighbour process within the target subdomain
       !
       nullify(norep_a,norep_b,norep_c,lnent_recv,lnent_send)
       call memory_alloca(memor_cou,'NOREP_A',vacal,norep_a,PAR_COUP_SIZE,PAR_COUP_SIZE)
       call memory_alloca(memor_cou,'NOREP_B',vacal,norep_b,PAR_COUP_SIZE,PAR_COUP_SIZE)
       norep_a(:,:)=0
       norep_b(:,:)=0
       sum_rep = 0

       !    calculate 
       !    a) norep_a: this table show the coupling between ranks in
       !    the coupling that each process will generate by assembling the
       !    interfaces of the transmission matrices. This happens in the wet
       !    points that are in the interior interfaces of the target subdomain:
       !    norep_a(jrank,irank) > 0 means irank will be coupled with jrank at
       !    extending the transmission matrices. #norep_a(jrank,irank) are the
       !    number of entries of ltransmat_target producing this new coupling.
       !
       !    b) sum_rep: are the non null entries in norep_a, so the new
       !    couplings that will be generated

       do ineig = 1, nneig_tar
          do icont = bound_size(ineig),bound_size(ineig+1)-1
             ipoin = bound_perm(icont)
             jpart = rank_targ2coup(comm_target % neights(ineig)+1) + 1
             do jcont = 1,lnent(ipoin)
                ipart = int(llent(ipoin) % a((jcont-1)*3+3),ip) + 1
                if( norep_a(jpart,ipart) == 0) sum_rep = sum_rep + 1
                norep_a(jpart,ipart) = norep_a(jpart,ipart) + 1
             enddo
          end do
       end do
       !
       ! norep_c is allocated, this buffer will be used for sotore columns that
       ! are sent discarding repited values
       ! 
       max_rep=maxval(norep_a)
       call memory_alloca(memor_cou,'NOREP_C',vacal,norep_c,max_rep*sum_rep)
       !
       ! An ordering is generated for the new complings
       !
       kcont = 1
       do icont   = 1, PAR_COUP_SIZE
          do jcont = 1, PAR_COUP_SIZE
             if(norep_a(jcont,icont) > 0)then
                norep_a(jcont,icont) = kcont
                kcont = kcont + 1
             endif
          enddo
       enddo

       call memory_alloca(memor_cou,'LNENT_RECV',vacal,lnent_recv,nneig_tar)
       call memory_alloca(memor_cou,'LNENT_SEND',vacal,lnent_send,nneig_tar)
       do ineig = 1, nneig_tar
          do icont = bound_size(ineig),bound_size(ineig+1)-1
             ipoin = bound_perm(icont)
             lnent_send(ineig) = lnent_send(ineig) + lnent(ipoin)
             !
             ! Apart than evaluating lnent_send, here bellow: 
             ! a) in norep_c we store the columns to be sent to each neighbour of
             ! the target subdomain (two entries with same column value only
             ! produce one entre in norep_c)
             ! b) lncol_demand_tar: number of coordinate values that will demand
             ! the curent target process from each coupling process to be
             ! transfered to the neighbour processes in the target subdomain
             !
             jpart = rank_targ2coup(comm_target % neights(ineig)+1) + 1
             do jcont = 1,lnent(ipoin)
                ipart = int(llent(ipoin) % a((jcont-1)*3+3),ip) + 1
                icol  = int(llent(ipoin) % a((jcont-1)*3+1),ip)
                auxlg = .true.
                do kcont = 1,norep_b(jpart,ipart)
                   if(norep_c((norep_a(jpart,ipart)-1) * max_rep + kcont) == icol) auxlg = .false.
                enddo
                if(auxlg) then
                   lncol_demand_tar(ipart) = lncol_demand_tar(ipart) + 1
                   norep_b(jpart,ipart) = norep_b(jpart,ipart) + 1
                   norep_c((norep_a(jpart,ipart)-1) * max_rep + norep_b(jpart,ipart))=icol
                endif
             enddo
          end do
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !
       ! 3.3 Send_recv #entries
       ! OPTIMIZATION: this can be done with an Alltoall!!
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,nneig_tar)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,nneig_tar
          ipart = comm_target % neights(ineig)          
          call PAR_SEND_RECEIVE(1_ip,1_ip,lnent_send(ineig:ineig),lnent_recv(ineig:ineig),&
               'IN CURRENT TARGET COLOR',ipart,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       !
       ! 4. Send/Recv entries 
       !
       ! 4.1 Fill buffers
       norep_a(:,:) = 0
       norep_b(:,:) = 0
       nullify(llent_send,llent_recv)
       call memory_alloca(memor_cou,'LLENT_SEND',vacal,llent_send,nneig_tar)
       call memory_alloca(memor_cou,'LLENT_RECV',vacal,llent_recv,nneig_tar)
       do ipart = 1, PAR_COUP_SIZE
          if(lncol_demand_tar(ipart) > 0 ) then
             call memory_alloca(memor_cou,'LLCOL_DEMAND_TAR % L',vacal,llcol_demand_tar(ipart) % l,&
                lncol_demand_tar(ipart)*2)
             lncol_demand_tar(ipart) = 0
          endif
       end do

       !
       ! OPTIMIZATION: Why do we reevaluate norep_a, and norep_c?????
       !
       do ineig = 1, nneig_tar
          jpart = rank_targ2coup(comm_target % neights(ineig)+1) + 1 
          do icont = bound_size(ineig),bound_size(ineig+1)-1
             ipoin = bound_perm(icont)
             do jcont = 1,lnent(ipoin)
                ipart = int(llent(ipoin) % a((jcont-1)*3+3),ip) + 1 
                if( norep_a(jpart,ipart) == 0) sum_rep = sum_rep + 1
                norep_a(jpart,ipart) = norep_a(jpart,ipart) + 1
             end do
          end do
       end do
       max_rep=maxval(norep_a)
       call memory_deallo(memor_cou,'NOREP_C',vacal,norep_c)
       nullify(norep_c)
       call memory_alloca(memor_cou,'NOREP_C',vacal,norep_c,max_rep*sum_rep)
       kcont = 1
       do icont = 1, PAR_COUP_SIZE
          do jcont = 1, PAR_COUP_SIZE
             if(norep_a(jcont,icont) > 0)then
                norep_a(jcont,icont) = kcont
                kcont = kcont + 1
             endif
          enddo
       enddo

       !
       ! In this loop are evaluated:
       ! a) llent_send: list of entries sent to each process, four values are
       ! sent for each entry [row,col,val,proc]
       ! b) llent_recv: buffer is just allocated according to lnent_recv
       ! c) lncol_demand_tar: number of coordinate values that will demand
       ! the curent target process from each coupling process to be transfered
       ! to the target processe in the target subdomain
       ! d) llcol_demand_tar: list of entries that link current target process others
       ! we send column and destination process
       ! 
       ! OPTIMIZATION: norep_a, norep_c are reevaluated, what for...mmm to be
       ! able to evaluate llcol_demand_tar
       ! 
       do ineig = 1, nneig_tar
          call memory_alloca(memor_cou,'LLENT_SEND % A',vacal,llent_send(ineig) % a,4*lnent_send(ineig))
          call memory_alloca(memor_cou,'LLENT_RECV % A',vacal,llent_recv(ineig) % a,4*lnent_recv(ineig))
          jpart = rank_targ2coup(comm_target % neights(ineig)+1) 
          kcont = 1
          do icont = bound_size(ineig),bound_size(ineig+1)-1
             ipoin = bound_perm(icont)
             do jcont = 1,lnent(ipoin)
                llent_send(ineig) % a(kcont) = real(icont - bound_size(ineig) + 1,rp) !row
                kcont = kcont + 1
                llent_send(ineig) % a(kcont:kcont+2) = llent(ipoin) % a((jcont-1)*3+1:jcont*3) !col,val,proc
                kcont = kcont + 3
                ipart = int(llent(ipoin) % a((jcont-1)*3+3),ip) + 1 
                icol  = int(llent(ipoin) % a((jcont-1)*3+1),ip)
                auxlg = .true.
                do lcont = 1,norep_b(jpart + 1,ipart)
                   if(norep_c((norep_a(jpart + 1,ipart)-1) * max_rep + lcont) == icol) auxlg = .false.
                enddo
                if(auxlg) then
                   llcol_demand_tar(ipart) % l(lncol_demand_tar(ipart)*2+1) = icol
                   llcol_demand_tar(ipart) % l(lncol_demand_tar(ipart)*2+2) = jpart 
                   lncol_demand_tar(ipart) = lncol_demand_tar(ipart) + 1
                   norep_b(jpart+1,ipart) = norep_b(jpart+1,ipart) + 1
                   norep_c((norep_a(jpart + 1,ipart)-1) * max_rep + norep_b(jpart+1,ipart))=icol
                endif
             end do
          end do
       end do
       call memory_deallo(memor_cou,'NOREP_A',vacal,norep_a)
       call memory_deallo(memor_cou,'NOREP_B',vacal,norep_b)
       call memory_deallo(memor_cou,'NOREP_C',vacal,norep_c)
       !
       ! 4.2 Send/Recv lists of entries
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,nneig_tar)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,nneig_tar
          ipart = comm_target % neights(ineig)
          nsend = lnent_send(ineig)*4
          nrecv = lnent_recv(ineig)*4
          if( nsend > 0 .and. nrecv > 0 ) then
             call PAR_SEND_RECEIVE(nsend,nrecv,llent_send(ineig) % a(1:nsend),llent_recv(ineig) % a(1:nrecv),'IN CURRENT TARGET COLOR',ipart,'NON BLOCKING')
          else
             if(nsend>0) then
                call PAR_SEND_RECEIVE(nsend,0_ip,llent_send(ineig) % a(1:nsend),dummr(1:1),'IN CURRENT TARGET COLOR',ipart,'NON BLOCKING')
             endif
             if(nrecv>0) then
                call PAR_SEND_RECEIVE(0_ip,nrecv,dummr(1:1),llent_recv(ineig) % a(1:nrecv),'IN CURRENT TARGET COLOR',ipart,'NON BLOCKING')
             end if
          end if
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !
       ! 5. Evaluate new column indices for the received entries (keep them on buffers)
       !
       ! Note: the algorithm will not work if there are repited entries in COO format
       ! this is avoided by calling before  matrix_COO_aggregate
       nullify(llcol_recv_ord,lncol_add,lnent_add,llent_add)
       call memory_alloca(memor_cou,'LNCOL_ADD',     vacal,lncol_add,PAR_COUP_SIZE)
       call memory_alloca(memor_cou,'LNENT_ADD',     vacal,lnent_add,PAR_COUP_SIZE)
       call memory_alloca(memor_cou,'LLENT_ADD',     vacal,llent_add,PAR_COUP_SIZE)
       call memory_alloca(memor_cou,'LLCOL_RECV_ORD',vacal,llcol_recv_ord,max_recv_size,PAR_COUP_SIZE)
       order_neights(:) = 0
       lnent_add(:) = 0
       lncol_add(:)  = 0
       do ineig = 1, nneig_tar
          order_neights(rank_targ2coup(comm_target % neights(ineig)+1)+1) = ineig
       end do
       do icont = 1,PAR_COUP_SIZE 
          ineig = order_neights(icont) !Is necessary to follow this order...
          if(ineig > 0) then
             if(lnent_recv(ineig)>0) then
                llcol_recv_ord(:,:) = 0
                do ibuff = 1,lnent_recv(ineig)  
                   ipart = int(llent_recv(ineig) % a((ibuff-1)*4+4),ip)
                   icol  = int(llent_recv(ineig) % a((ibuff-1)*4+2),ip)
                   lnent_add(ipart+1) = lnent_add(ipart+1) + 1
                   llcol_recv_ord(icol,ipart+1) = 1 
                end do
                do jcont = 1,PAR_COUP_SIZE
                   do kcont = 1,max_recv_size
                      if(llcol_recv_ord(kcont,jcont)>0) then
                         lncol_add(jcont) = lncol_add(jcont) + 1
                         llcol_recv_ord(kcont,jcont) = lncol_add(jcont);
                      end if
                   end do
                end do
                do ibuff = 1,lnent_recv(ineig)   
                   ipart = int(llent_recv(ineig) % a((ibuff-1)*4+4),ip) + 1
                   icol  = int(llent_recv(ineig) % a((ibuff-1)*4+2),ip)
                   llent_recv(ineig) % a((ibuff-1)*4+2) = real(llcol_recv_ord(icol,ipart),rp)
                end do
             end if
          end if
       end do
       do ipart = 1,PAR_COUP_SIZE
         if(lnent_add(ipart) > 0) then
            call memory_alloca(memor_cou,'LLENT_ADD % A',vacal,llent_add(ipart) % a,lnent_add(ipart)*4)
            lnent_add(ipart) = 0
         endif
       end do
       do ineig =1,nneig_tar
          do ibuff=1,lnent_recv(ineig)
             ipart = int(llent_recv(ineig) % a((ibuff-1)*4+4),ip)
             llent_add(ipart + 1) % a(lnent_add(ipart + 1)*4 + 1) = llent_recv(ineig) % a((ibuff-1)*4+1) !row
             llent_add(ipart + 1) % a(lnent_add(ipart + 1)*4 + 2) = llent_recv(ineig) % a((ibuff-1)*4+2) !col
             llent_add(ipart + 1) % a(lnent_add(ipart + 1)*4 + 3) = llent_recv(ineig) % a((ibuff-1)*4+3) !val
             llent_add(ipart + 1) % a(lnent_add(ipart + 1)*4 + 4) = real(ineig,rp)
             lnent_add(ipart + 1) = lnent_add(ipart + 1) + 1
          enddo
       enddo
       !
       ! 6. Evaluate new comm_data_par for coupling (target)
       !
       dim_prev_tar   = coupling % commd % lrecv_dim
       nneig_new_tar  = nneig_prev
       dim_new_tar    = dim_prev_tar
       order_neights(:) = 0
       do ineig = 1, coupling % commd % nneig
          order_neights(coupling % commd % neights(ineig)+1) = ineig
       end do
       do ipart = 1,PAR_COUP_SIZE
          if(lncol_add(ipart) > 0 ) then
             if(order_neights(ipart) == 0) then
                nneig_new_tar = nneig_new_tar + 1_ip
                order_neights(ipart) = nneig_new_tar
             endif
             dim_new_tar = dim_new_tar + lncol_add(ipart)
          endif
       end do
       if(dim_prev_tar /= dim_new_tar ) then
          nullify(lrecv_size_new,neights_new_tar);
          call memory_alloca(memor_cou,'LRECV_SIZE_NEW' ,vacal,lrecv_size_new,nneig_new_tar+1_ip)
          call memory_alloca(memor_cou,'NEIGHTS_NEW_TAR',vacal,neights_new_tar,nneig_new_tar)
          coupling % commd % nneig = nneig_new_tar
          coupling % commd % lrecv_dim = dim_new_tar
          neights_new_tar(:)=0
          do ipart = 1,PAR_COUP_SIZE
             ineig = order_neights(ipart)
             if(ineig /= 0) then
                neights_new_tar(ineig) = ipart-1
             endif
          end do
          sum_entr =0
          lrecv_size_new(1) = 1
          do ineig = 1, coupling % commd % nneig
             ipart = neights_new_tar(ineig)
             sum_entr = sum_entr + lncol_add(ipart + 1)
             if(ineig < nneig_prev+1_ip) then
                lrecv_size_new(ineig+1) = coupling % commd % lrecv_size(ineig+1) + sum_entr
             else
                lrecv_size_new(ineig+1) = lrecv_size_new(ineig) + lncol_add(ipart + 1)
             endif
          end do
          if(lrecv_size_new(nneig_new_tar+1_ip)-1_ip /= dim_new_tar)&
               call runend("Something wrong at the target comm_data_par")
       endif
       !
       ! 7. Introduce new entries to the transmission matrices
       ! 
       if(dim_prev_tar /= dim_new_tar) then
          !7.1 Invert lpoint_wet
          nullify(inv_lpwet)
          call memory_alloca(memor_cou,'INV_LPWET',vacal,inv_lpwet,npoin)
          inv_lpwet(1:npoin) = 0
          do iwetp=1,npoin_wet
             inv_lpwet(coupling % wet % lpoin_wet(iwetp)) = iwetp
          end do
          !7.3 Allocate and fill new transmission matrices
          call memory_alloca(memor_cou,'M_NEW',vacal,M_new,nneig_new_tar)
          do ineig = 1,nneig_new_tar
             nentr_prev = 0
             incr = 0
             if(ineig < nneig_prev + 1_ip) then
                if (associated(M(ineig) % iA )) nentr_prev = size(M(ineig) % iA)
                incr = lrecv_size_new(ineig) - coupling % commd % lrecv_size(ineig)
                if(incr < 0) call runend("PARALLELIZE TRANSMAT: incr<0 not expected situation")
             endif
             nentr = nentr_prev + lnent_add(neights_new_tar(ineig)+1)
             call memory_alloca(memor_cou,'M_NEW(INEIG)',vacal,M_new(ineig),1_ip,nentr)
             
             M_new(ineig) % nrows = npoin_wet
             M_new(ineig) % ncols = coupling % commd % lrecv_dim
             M_new(ineig) % ndof1 = 1
             M_new(ineig) % ndof2 = 1
             do ientr = 1,nentr_prev
                weight = 1.0_rp
                if( coupling % itype == ELEMENT_INTERPOLATION  .or. &
                     coupling % itype == NEAREST_BOUNDARY_NODE  .or. &
                     coupling % itype == BOUNDARY_INTERPOLATION .or. &
                     coupling % itype == NEAREST_ELEMENT_NODE   .or. &
                     coupling % itype == GLOBAL_NUMBERING       ) then
                   aux = coupling % wet % lpoin_wet(M(ineig) % iA(ientr))
                   weight = 1.0_rp / multip(aux)
                endif
                M_new(ineig) % iA(ientr) = M(ineig) % iA(ientr)
                M_new(ineig) % jA(ientr) = M(ineig) % jA(ientr) + incr
                M_new(ineig) % vA(1,1,ientr) = M(ineig) % vA(1,1,ientr) * weight
             end do

             ipart = neights_new_tar(ineig) + 1
             do ientr = nentr_prev+1,nentr
                icont = ientr - nentr_prev
                if(ineig < nneig_prev + 1_ip) then
                   incr = lrecv_size_new(ineig) + coupling % commd % lrecv_size(ineig+1)&
                        - coupling % commd % lrecv_size(ineig) - 1
                else
                   incr = lrecv_size_new(ineig) - 1
                endif
                weight = 1.0_rp
                aux = bound_perm(bound_size(int(llent_add(ipart) % a((icont-1)*4+4),ip)) + int(llent_add(ipart) % a((icont-1)*4+1),ip) - 1)
                if( coupling % itype == ELEMENT_INTERPOLATION  .or. &
                     coupling % itype == NEAREST_BOUNDARY_NODE  .or. &
                     coupling % itype == BOUNDARY_INTERPOLATION .or. &
                     coupling % itype == NEAREST_ELEMENT_NODE   .or. &
                     coupling % itype == GLOBAL_NUMBERING       ) then
                   weight = 1.0_rp / multip(aux)
                endif
                irow = inv_lpwet(aux)
                M_new(ineig) % iA(ientr) = irow               
                M_new(ineig) % jA(ientr) = int(llent_add(ipart) % a((icont-1)*4+2),ip) + incr
                M_new(ineig) % vA(1,1,ientr) = llent_add(ipart) % a((icont-1)*4+3) * weight
             end do
          end do
          call memory_deallo(memor_cou,'INV_LPWET',vacal,inv_lpwet)
       end if
       !
       ! 8. Deallocate arrays at target
       !
       call memory_deallo(memor_cou,'LNENT'         ,vacal,lnent)
       call memory_deallo(memor_cou,'LLENT'         ,vacal,llent)
       call memory_deallo(memor_cou,'LLENT_ADD'     ,vacal,llent_add)
       call memory_deallo(memor_cou,'LLENT_SEND'    ,vacal,llent_send)
       call memory_deallo(memor_cou,'LLENT_RECV'    ,vacal,llent_recv)
       call memory_deallo(memor_cou,'LNENT_SEND'    ,vacal,lnent_send)
       call memory_deallo(memor_cou,'LNENT_RECV'    ,vacal,lnent_recv)
       call memory_deallo(memor_cou,'LLCOL_RECV_ORD',vacal,llcol_recv_ord)
       call memory_deallo(memor_cou,'LNCOL_ADD'     ,vacal,lncol_add)
       call memory_deallo(memor_cou,'LNENT_ADD'     ,vacal,lnent_add)
    endif
    if(associated(rank_targ2coup)) call memory_deallo(memor_cou,'RANK_TARG2COUP',vacal,rank_targ2coup)
    !
    ! 9. Evaluate new comm_data_par for coupling (target)
    !
    ! 9.1 Sorce send num. of new elements it must send
    call PAR_BARRIER('IN CURRENT COUPLING')
    call PAR_START_NON_BLOCKING_COMM(1_ip, nneig_prev)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,nneig_prev
       ipart = coupling % commd % neights(ineig)
       aux = ipart + 1      
       call PAR_SEND_RECEIVE(1_ip,1_ip,lncol_demand_tar(aux:aux),lncol_send_src(aux:aux),&
            'IN CURRENT COUPLING',ipart,'NON BLOCKING')
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    ! 9.2 Send/Recv the elements
    nullify(llcol_send_src)
    call memory_alloca(memor_cou,'llcol_send_src',vacal,llcol_send_src,PAR_COUP_SIZE)
    call PAR_START_NON_BLOCKING_COMM(1_ip,nneig_prev)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,nneig_prev
       ipart = coupling % commd % neights(ineig)
       nsend = lncol_demand_tar(ipart+1)*2
       nrecv = lncol_send_src(ipart+1)*2
       if( nsend > 0 .and. nrecv > 0 ) then
          call memory_alloca(memor_cou,'llcol_send_src(ipart)',vacal,llcol_send_src(ipart+1) % l,nrecv)
          call PAR_SEND_RECEIVE(nsend,nrecv,llcol_demand_tar(ipart+1) % l(1:nsend),llcol_send_src(ipart+1) % l(1:nrecv),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
       else
          if(nsend>0) then
             call PAR_SEND_RECEIVE(nsend,0_ip,llcol_demand_tar(ipart+1) % l(1:nsend),dummi(1:1),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
          endif
          if(nrecv>0) then
             call memory_alloca(memor_cou,'llcol_send_src(ipart)',vacal,llcol_send_src(ipart+1) % l,nrecv)
             call PAR_SEND_RECEIVE(0_ip,nrecv,dummi(1:1),llcol_send_src(ipart+1) % l(1:nrecv),'IN CURRENT COUPLING',ipart,'NON BLOCKING')
          end if
       end if
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    if(I_AM_SOURCE) then
       !
       ! 10. Fill-in max_send_size
       !
       aux1(1) = 1000 !seciroty threshold
       if( INOTMASTER ) then 
          do ineig = 1, int(memory_size(coupling % commd % lsend_size),ip)-1_ip
             aux = coupling % commd % lsend_size(ineig+1) - coupling % commd % lsend_size(ineig)
             if(aux > aux1(1)) aux1(1) = aux 
          end do
       endif
       call PAR_MAX(1_ip,aux1,'IN CURRENT SOURCE COLOR')
       max_send_size = aux1(1)
    endif
    if(I_AM_SOURCE .and. INOTMASTER) then
       !
       ! 11. Continue with comm_data_par
       !
       ! Note: the algorithm will not work if there are repited entries in COO format
       !order elements of llcol_send_src
       nullify(llcol_send_ord)
       call memory_alloca(memor_cou,'llcol_send_ord',vacal,llcol_send_ord,max_send_size,PAR_COUP_SIZE)
       do ipart=1,PAR_COUP_SIZE
          if(lncol_send_src(ipart)>0) then
             llcol_send_ord(:,:)=0
             do icont = 1,lncol_send_src(ipart)
                jcol = llcol_send_src(ipart) % l((icont-1)*2+1)
                jpart =llcol_send_src(ipart) % l((icont-1)*2+2) + 1
                llcol_send_ord(jcol,jpart) = 1
             enddo
             kcont = 0
             do icont = 1,PAR_COUP_SIZE
                do jcont = 1,max_send_size
                   if(llcol_send_ord(jcont,icont) > 0) then
                      llcol_send_src(ipart) % l(kcont*2+1) = jcont
                      llcol_send_src(ipart) % l(kcont*2+2) = icont - 1
                      kcont = kcont + 1
                   endif
                enddo
             enddo
          endif
       enddo

       nullify(lnsend_add,llsend_add) 
       call memory_alloca(memor_cou,'lnsend_add',vacal,lnsend_add,PAR_COUP_SIZE)
       call memory_alloca(memor_cou,'llsend_add',vacal,llsend_add,PAR_COUP_SIZE)
       lnsend_add(:)=0
       do ineig = 1,nneig_prev
          ipart = coupling % commd % neights(ineig)
          do icont = 1,lncol_send_src(ipart+1)
             jpart = llcol_send_src(ipart+1) % l((icont-1)*2+2) + 1
             lnsend_add(jpart) = lnsend_add(jpart) + 1
          end do
       end do
       do ipart = 1,PAR_COUP_SIZE
          if(lnsend_add(ipart) > 0) then
             call memory_alloca(memor_cou,'LNSEND_ADD % L',vacal,llsend_add(ipart) % l,lnsend_add(ipart))
             lnsend_add(ipart) = 0
          end if
       end do
       do ineig = 1,nneig_prev
          ipart = coupling % commd % neights(ineig)
          do icont = 1,lncol_send_src(ipart+1)
             jpart   = llcol_send_src(ipart+1) % l((icont-1)*2+2) + 1
             permval = coupling % commd % lsend_perm(coupling % commd % lsend_size(ineig) + llcol_send_src(ipart+1) % &
                  l((icont-1)*2+1)-1)
             lnsend_add(jpart) = lnsend_add(jpart) + 1
             llsend_add(jpart) % l(lnsend_add(jpart)) = permval
          end do
       end do
       dim_prev_src   = coupling % commd % lsend_dim
       nneig_new_src  = nneig_prev
       dim_new_src = dim_prev_src
       order_neights = 0
       do ineig = 1, nneig_prev
          order_neights(coupling % commd % neights(ineig)+1) = ineig
       end do
       do ipart = 1,PAR_COUP_SIZE
          if(lnsend_add(ipart) > 0) then
             if(order_neights(ipart) == 0) then
                nneig_new_src = nneig_new_src + 1_ip
                order_neights(ipart) = nneig_new_src
             endif
             dim_new_src = dim_new_src + lnsend_add(ipart)
          endif
       end do
       if(dim_prev_src /= dim_new_src ) then
          nullify(lsend_size_new,neights_new_src,perm_new_src)
          call memory_alloca(memor_cou,'neights_new_src',vacal,neights_new_src,nneig_new_src)
          call memory_alloca(memor_cou,'lsend_size_new' ,vacal,lsend_size_new,nneig_new_src+1_ip)
          call memory_alloca(memor_cou,'perm_new_src'   ,vacal,perm_new_src,dim_new_src)
          coupling % commd % nneig = nneig_new_src
          coupling % commd % lsend_dim = dim_new_src
          do ipart = 1,PAR_COUP_SIZE
             ineig = order_neights(ipart)
             if(ineig /= 0) then
                neights_new_src(ineig) = ipart-1
             endif
          end do
          sum_entr =0
          lsend_size_new(1) = 1
          do ineig = 1, coupling % commd % nneig
             ipart = neights_new_src(ineig)
             sum_entr = sum_entr + lnsend_add(ipart + 1)
             if(ineig < nneig_prev+1_ip) then
                lsend_size_new(ineig+1) = coupling % commd % lsend_size(ineig+1) + sum_entr
             else
                lsend_size_new(ineig+1) = lsend_size_new(ineig) + lnsend_add(ipart + 1)
             endif
          end do
          if(lsend_size_new(nneig_new_src+1_ip)-1_ip /= dim_new_src)&
               call runend("Something wrong at the target comm_data_par 2")
          ! Fill-in new lsend_perm
          do ineig = 1, nneig_prev
             sum_send = coupling % commd % lsend_size(ineig+1)&
                  - coupling % commd % lsend_size(ineig)
             jcont = lsend_size_new(ineig)
             kcont = coupling % commd % lsend_size(ineig)
             do icont = 1,sum_send
                perm_new_src(jcont) = coupling % commd % lsend_perm(kcont)
                jcont = jcont + 1 
                kcont = kcont + 1
             end do
          end do
          do ineig = 1, coupling % commd % nneig
             ipart = neights_new_src(ineig) + 1
             if(lnsend_add(ipart) > 0) then
                if(ineig <= nneig_prev) then
                   jcont = lsend_size_new(ineig) + coupling % commd % lsend_size(ineig+1)&
                        - coupling % commd % lsend_size(ineig)
                else
                   jcont = lsend_size_new(ineig)
                endif
                do icont = 1,lnsend_add(ipart)
                   perm_new_src(jcont) = llsend_add(ipart) % l(icont)
                   jcont = jcont + 1
                end do
             end if
          end do
       endif
       !
       ! 12. Deallocate pointers at source
       !
       call memory_deallo(memor_cou,'llcol_send_ord',vacal,llcol_send_ord)
       call memory_deallo(memor_cou,'lnsend_add'    ,vacal,lnsend_add)
       call memory_deallo(memor_cou,'llsend_add'    ,vacal,llsend_add)
    endif
    if(INOTMASTER .and. I_AM_TARGET) then
       if(dim_prev_tar /= dim_new_tar) then
          nullify(M)
          if(associated(coupling % commd % neights)) &
               call memory_deallo(memor_cou,'coupling % commd % neights',vacal, coupling % commd % neights)
          if(associated(coupling % commd % lrecv_size)) &
             call memory_deallo(memor_cou,'coupling % commd % lrecv_size',vacal, coupling % commd % lrecv_size)
          if(associated(coupling % commd % lsend_size) .and. (.not. I_AM_SOURCE)) &
             call memory_deallo(memor_cou,'coupling % commd % lsend_size',vacal, coupling % commd % lsend_size)
          if(associated(coupling % ltransmat_target)) &
               call memory_deallo(memor_cou,'ltransmat_target',vacal,coupling % ltransmat_target)
          call memory_alloca(memor_cou,'coupling % commd % neights',vacal, coupling % commd % neights,nneig_new_tar)
          if(.not. I_AM_SOURCE) call memory_alloca(memor_cou,'coupling % commd % lsend_size',vacal, &
               coupling % commd % lsend_size,nneig_new_tar+1_ip)
          call memory_alloca(memor_cou,'coupling % commd % lrecv_size',vacal, coupling % commd % lrecv_size,nneig_new_tar+1_ip)
          call memory_alloca(memor_cou,'ltransmat_target',vacal,coupling % ltransmat_target,nneig_new_tar)
          M => coupling % ltransmat_target
          coupling % commd % lrecv_size(:) = lrecv_size_new(:)
          coupling % commd % neights(:) = neights_new_tar(:)
          do ineig = 1,nneig_new_tar
             if(associated(M_new(ineig) % iA)) then
                nentr = size(M_new(ineig) % iA)
                call memory_alloca(memor_cou,'M(INEIG)',vacal,M(ineig),1_ip,nentr)
                M(ineig) % nrows = M_new(ineig) % nrows
                M(ineig) % ncols = M_new(ineig) % ncols
                M(ineig) % ndof1 = M_new(ineig) % ndof1 
                M(ineig) % ndof2 = M_new(ineig) % ndof2 
                do ientr = 1,nentr
                   M(ineig) % iA(ientr) = M_new(ineig) % iA(ientr)
                   M(ineig) % jA(ientr) = M_new(ineig) % jA(ientr)
                   M(ineig) % vA(1,1,ientr) = M_new(ineig) % vA(1,1,ientr)
                end do
             endif
          end do
          call memory_deallo(memor_cou,'neights_new_tar',vacal,neights_new_tar)
          call memory_deallo(memor_cou,'lrecv_size_new' ,vacal,lrecv_size_new)
       endif
    endif
    if(INOTMASTER .and. I_AM_SOURCE) then
       if(dim_prev_src /= dim_new_src ) then
          if(associated(coupling % commd % neights)) &
               call memory_deallo(memor_cou,'coupling % commd % neights',vacal, coupling % commd % neights)
          if(associated(coupling % commd % lsend_size)) &
               call memory_deallo(memor_cou,'coupling % commd % lsend_size',vacal, coupling % commd % lsend_size)
          if(associated(coupling % commd % lrecv_size) .and. (.not. I_AM_TARGET)) &
               call memory_deallo(memor_cou,'coupling % commd % lsend_size',vacal, coupling % commd % lrecv_size)
          if(associated(coupling % commd % lsend_perm)) &
               call memory_deallo(memor_cou,'coupling % commd % lsend_perm',vacal, coupling % commd % lsend_perm)
          call memory_alloca(memor_cou,'coupling % commd % neights',vacal, coupling % commd % neights,nneig_new_src)
          if(.not. I_AM_TARGET) call memory_alloca(memor_cou,'coupling % commd % lrecv_size',vacal,&
               coupling % commd % lrecv_size,nneig_new_src+1_ip)
          call memory_alloca(memor_cou,'coupling % commd % lsend_size',vacal, coupling % commd % lsend_size,nneig_new_src+1_ip)
          call memory_alloca(memor_cou,'coupling % commd % lsend_perm',vacal, coupling % commd % lsend_perm,dim_new_src)
          coupling % commd % lsend_size(:) = lsend_size_new(:)
          coupling % commd % neights(:) = neights_new_src(:)
          coupling % commd % lsend_perm(:) = perm_new_src(:)
          call memory_deallo(memor_cou,'neights_new_src',vacal,neights_new_src)
          call memory_deallo(memor_cou,'lsend_size_new',vacal,lsend_size_new)
          call memory_deallo(memor_cou,'perm_new_src',vacal,perm_new_src)
       endif
    endif

    if(I_AM_SOURCE .and. I_AM_TARGET .and. nneig_new_tar /= nneig_new_src) then
       print *,"HOLA: ",PAR_MY_CODE_RANK," ",coupling % commd % nneig," ",I_AM_SOURCE,&
            " ",nneig_new_src," ",I_AM_TARGET," ",nneig_new_tar
       call runend("There is a problem to be fixed, ask Ricard")
    endif
    !
    ! 11. Deallocate global pointers
    !
    if(associated(M_new)) call memory_deallo(memor_cou,'M_new',vacal,M_new)
    call memory_deallo(memor_cou,'order_neights'   ,vacal,order_neights)
    call memory_deallo(memor_cou,'lncol_send_src'  ,vacal,lncol_send_src)
    call memory_deallo(memor_cou,'lncol_demand_tar',vacal,lncol_demand_tar)
    call memory_deallo(memor_cou,'llcol_demand_tar',vacal,llcol_demand_tar)
    call memory_deallo(memor_cou,'llcol_send_src'  ,vacal,llcol_send_src)
    if(I_AM_TARGET) call memory_deallo(memor_cou,'multip',vacal,multip)
    !
    ! Timer
    !
    call cputim(time2)
    coupling % cputim(7) = coupling % cputim(7) + time2 - time1

  end subroutine COU_PARALLELIZE_TRANSMISSION_MATRICES

end module mod_couplings_communications
!-----------------------------------------------------------------------
