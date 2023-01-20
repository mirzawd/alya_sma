!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_redistribution

  use def_kintyp_basic,    only : ip, rp, r1p, lg, i1p
  use def_kintyp_comm,     only : comm_data_par
  use mod_parall,          only : par_memor,PAR_COMM_MY_CODE
  use mod_parall,          only : comm_redistribution
  use def_parall,          only : method_redistribution_par
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_deallo
  use mod_memory,          only : memory_size
  use mod_communications,  only : PAR_MAX
  use mod_communications,  only : PAR_SEND_RECEIVE 
  use mod_communications,  only : PAR_ALLTOALL
  use mod_communications,  only : PAR_ALLTOALLV
  use mod_communications,  only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications,  only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,  only : PAR_START_NON_BLOCKING_COMM
  use mod_communications,  only : PAR_END_NON_BLOCKING_COMM
  use mod_communications,  only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use def_domain,          only : nnode,memor_dom
  use def_domain
  use def_master
  use mod_messages
  use def_mpi
#include "def_mpi.inc"
  
  use mod_communications

  implicit none

  private

  interface redistribution_array
     module procedure &
          &           redistribution_array_IP_1 , &
          &           redistribution_array_IP_2 , &
          &           redistribution_array_RP_1 , &
          &           redistribution_array_RP_2 , &
          &           redistribution_array_RP_3 , &
          &           redistribution_array_RP_4 , &
          &           redistribution_array_R3P_1, &
          &           redistribution_array_R2P_1, &
          &           redistribution_array_R1P_1           
  end interface redistribution_array

  MY_MPI_COMM         :: PAR_COMM_RED
  MY_MPI_COMM         :: PAR_COMM_RED4
  integer(ip)         :: PAR_COMM_RED_RANK
  integer(ip)         :: PAR_COMM_RED_SIZE
  type(comm_data_par) :: commd_nelem      
  type(comm_data_par) :: commd_npoin
  type(comm_data_par) :: commd_nboun     
  type(comm_data_par) :: commd_nbopo     
  character(50)       :: my_variable_name

  type :: hash_chunkdist
     integer(ip)             :: chunk
     integer(ip)             :: nchunk
     integer(ip),  pointer   :: offset(:)
     type(i1p),    pointer   :: perm(:)
  end type hash_chunkdist

  public :: commd_npoin
  public :: commd_nelem
  public :: commd_nboun
  public :: commd_nbopo
  
  public :: redistribution_renumbering_npoin_communicator
  public :: redistribution_renumbering_nelem_communicator
  public :: redistribution_renumbering_nboun_communicator

  character(100), PARAMETER :: vacal = "mod_redistribution"
  interface generate_comm_chunkdist
     module procedure generate_comm_chunkdist_a,&
          &           generate_comm_chunkdist_b
  end interface generate_comm_chunkdist
  !
  ! Public functions
  !
  public :: generate_comm_chunkdist
  public :: generate_comm_sizes
  public :: redistribution_domain
  public :: redistribution_array
  public :: redistribution_array_RP_3
  public :: redistribution_comm_initialization
  public :: redistribution_lid_chunkdist
  public :: redistribution_generate_hash_chunkdist
  public :: redistribution_deallocate_hash
  public :: hash_chunkdist
  public :: redistribution_copy_communicator
  public :: redistribution_initialization

contains

  subroutine redistribution_initialization()

      call commd_npoin % init(COMM_NAME='COMMD_NPOIN')
      call commd_nelem % init(COMM_NAME='COMMD_NELEM')
      call commd_nboun % init(COMM_NAME='COMMD_NBOUN')
      call commd_nbopo % init(COMM_NAME='COMMD_NBOPO')

    end subroutine redistribution_initialization
    
  subroutine redistribution_comm_initialization(PAR_COMM)

    MY_MPI_COMM   , intent(in) :: PAR_COMM

    PAR_COMM_RED  = PAR_COMM
    PAR_COMM_RED4 = PAR_COMM
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_RED,PAR_COMM_RED_RANK,PAR_COMM_RED_SIZE)
    
  end subroutine redistribution_comm_initialization

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell 
  !> @date    28/05/2018 
  !> @brief   Generate communicators
  !> @details Generate elements, nodes and boundary nodes communicators
  !           for mesh redistribution
  !
  !----------------------------------------------------------------------

  subroutine redistribution_copy_communicator(comm)

    implicit none

    type(comm_redistribution),  intent(inout) :: comm

    call comm % commd_npoin % copy(commd_npoin)
    call comm % commd_nelem % copy(commd_nelem)
    call comm % commd_nboun % copy(commd_nboun)
    call comm % commd_nbopo % copy(commd_nbopo)    

  end subroutine redistribution_copy_communicator
  
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell 
  !> @date    28/05/2018 
  !> @brief   Generate communicators
  !> @details Generate elements, nodes and boundary nodes communicators
  !           for mesh redistribution
  !
  !----------------------------------------------------------------------

  subroutine generate_comms(new_part,comm_elem,comm_poin,comm_boun,comm_bopo)

    implicit none

    integer(ip), pointer,           intent(in)    :: new_part(:)    
    type(comm_data_par),            intent(inout) :: comm_elem      
    type(comm_data_par),            intent(inout) :: comm_poin 
    type(comm_data_par),            intent(inout) :: comm_boun      
    type(comm_data_par),  optional, intent(inout) :: comm_bopo      

    call generate_comm_elem(new_part,comm_elem)
    call generate_comm_poin(new_part,comm_poin)
    call generate_comm_boun(new_part,comm_boun)
    if( present(comm_bopo) ) call generate_comm_boun(new_part,comm_bopo)

  end subroutine generate_comms

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate elements communicator
  !> @details
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_elem(new_part,comm_elem, nenti_,COMM_NAME)

    use def_domain,  only : nelem

    implicit none

    integer(ip), pointer, intent(in)    :: new_part(:)   
    type(comm_data_par),  intent(inout) :: comm_elem      
    integer(ip), optional,intent(in)    :: nenti_
    character(*),optional, intent(in)   :: COMM_NAME
    integer(ip)           :: ielem, ipart
    integer(ip)           :: isend, nsendp, ibuff, iperm
    integer(ip), pointer  :: all2send(:)
    integer(ip), pointer  :: nperm(:)
    integer(4),  pointer  :: lnsend(:)
    type(i1p),   pointer  :: send_perm(:)
    integer(ip)           :: nenti
    integer(4)            :: PAR_COMM_RED_SIZE4
    character(20)         :: my_comm_name

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM_ELEM'
    end if

    call comm_elem % deallo(par_memor,COMM_NAME=trim(my_comm_name))
    call comm_elem % init  (COMM_NAME='COM_ELEM')
    !
    ! Nullify
    !
    nullify(all2send)
    nullify(nperm)
    nullify(lnsend)
    nullify(send_perm)      
    !
    ! Check arguments and allocate pointers
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti) call runend("mod_redistribution: wrong size new_part")           
    !
    ! Evaluate lnsend
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'LNSEND',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)
    lnsend = 0_4
    do ielem = 1,nenti
       if( new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip ) then
          print*,'error=', new_part(ielem),PAR_COMM_RED_SIZE
          call runend("mod_redistribution: wrong values partition ranks in new_part")
       else
          lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4
       endif
    enddo
    !
    ! Eval comm_elem sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_elem,COMM_NAME=comm_name)      
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'ALL2SEND', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'NPERM',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'SEND_PERM',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0_ip,PAR_COMM_RED_SIZE-1
       if(lnsend(ipart) > 0_4) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'SEND_PERM % L',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do ielem = 1,nenti
       isend = all2send(new_part(ielem))
       nperm(isend) = nperm(isend) + 1_ip
       send_perm(isend) % l( nperm(isend) ) = ielem
    enddo

    iperm=1_ip
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_elem % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,'LNSEND'   ,vacal,lnsend)
    call memory_deallo(par_memor,'SEND_PERM',vacal,send_perm)
    call memory_deallo(par_memor,'ALL2SEND' , vacal,all2send)
    call memory_deallo(par_memor,'NPERM'    , vacal,nperm)

  end subroutine generate_comm_elem

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate nodes communicator
  !> @details nenti_ refers to number of elements to consider (outer loop is
  !           through elements).
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_bopo(new_part,comm_bopo,nenti_,COMM_NAME)

    use def_domain,  only : npoin, nnode
    use def_domain,  only : lnods, nelem, ltype, lpoty

    implicit none

    integer(ip),         pointer,  intent(in)    :: new_part(:)   
    type(comm_data_par),           intent(inout) :: comm_bopo      
    integer(ip),         optional, intent(in)    :: nenti_
    character(*),        optional, intent(in)    :: COMM_NAME

    integer(ip)                                  :: ipoin, ielem, inode
    integer(ip)                                  :: isend, ipart, nsendp
    integer(ip)                                  :: melem, ibuff, iperm
    integer(ip)                                  :: nenti, isin,  ibopo
    integer(ip),         pointer                 :: lelems(:,:)
    integer(ip),         pointer                 :: lnelem(:)
    integer(ip),         pointer                 :: all2send(:)
    integer(ip),         pointer                 :: nperm(:)
    integer(ip),         pointer                 :: lpoty_inv(:)
    integer(4),          pointer                 :: lnsend(:)
    type(i1p),           pointer                 :: send_perm(:)
    integer(4)                                   :: PAR_COMM_RED_SIZE4
    character(20)                                :: my_comm_name

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM_BOPO'
    end if
 
    call comm_bopo % deallo(par_memor,COMM_NAME=trim(my_comm_name))
    call comm_bopo % init  (COMM_NAME='COMM_BOPO')
    !
    ! Nullify
    !
    nullify(lelems)
    nullify(lnelem)
    nullify(all2send)
    nullify(nperm)
    nullify(lpoty_inv)
    nullify(lnsend)
    nullify(send_perm)
    !
    ! Check arguments
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti)   call runend("mod_redistribution: wrong size new_part")
    if(nenti /= 0_ip .and. nenti /= nelem) call runend("mod_redistribution: option nenti not encoded")
    !
    ! LPOTY_INV
    !
    call memory_alloca(par_memor,'lpoty_inv',vacal,lpoty_inv,max(1_ip,npoin))
    do ibopo = 1,nbopo
       ipoin = lpoty(ibopo)
       lpoty_inv(ipoin) = 1
    end do
    !
    ! Generate lnelem and lelems 
    !
    call memory_alloca(par_memor,'lnelem',vacal,lnelem,max(1_ip,npoin))

    lnelem = 0_ip
    melem = 0_ip
    do ielem = 1,nenti
       do ipoin = 1,nnode(ltype(ielem))
          inode = lnods(ipoin,ielem)
          if( lpoty_inv(inode) > 0 ) then
             lnelem(inode) = lnelem(inode) + 1_ip 
             if(lnelem(inode) > melem ) melem = lnelem(inode)
          end if
       enddo
    enddo

    lnelem = 0_ip
    call memory_alloca(par_memor,'lelems',vacal,lelems,melem,max(1_ip,npoin))
    do ielem = 1,nenti
       do ipoin = 1,nnode(ltype(ielem))
          inode = lnods(ipoin,ielem)
          if( lpoty_inv(inode) > 0 ) then
             isin=0_ip;
             do ibuff=1_ip,lnelem(inode)
                if(new_part(lelems(ibuff,inode))==new_part(ielem)) isin = isin + 1_ip;
             enddo
             if(isin == 0_ip) then
                lnelem(inode) = lnelem(inode) + 1_ip
                lelems(lnelem(inode),inode) = ielem 
             else if(isin > 1_ip) then
                call runend("mod_redistribution: something broken isin > 1")
             endif
          end if
       enddo
    enddo
    !
    ! Eval lnsend
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)
    do ipoin = 1,npoin
       if( lpoty_inv(ipoin) > 0 ) then
          do ibuff = 1,lnelem(ipoin)
             ielem = lelems(ibuff,ipoin)
             if(new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip) then
                call runend("mod_redistribution: wrong values partition ranks in new_part")
             else
                lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4 
             endif
          enddo
       end if
    enddo
    !
    ! Eval comm_bopo sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_bopo,COMM_NAME=comm_name)
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'ALL2SEND', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'nperm',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'send_perm',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0,PAR_COMM_RED_SIZE-1
       if( lnsend(ipart) > 0_4 ) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'send_perm % l',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do ipoin = 1,npoin
       if( lpoty_inv(ipoin) > 0 ) then
          do ibuff = 1,lnelem(ipoin)
             ielem = lelems(ibuff,ipoin)
             isend = all2send(new_part(ielem))
             nperm(isend) = nperm(isend) + 1_ip
             send_perm(isend) % l( nperm(isend) ) = ipoin
          enddo
       end if
    enddo

    iperm=1
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_bopo % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,'lnelem'   ,vacal,lnelem)
    call memory_deallo(par_memor,'lelems'   ,vacal,lelems)
    call memory_deallo(par_memor,'lnsend'   ,vacal,lnsend)
    call memory_deallo(par_memor,'send_perm',vacal,send_perm)
    call memory_deallo(par_memor,'all2send' ,vacal,all2send)
    call memory_deallo(par_memor,'nperm'    ,vacal,nperm)
    call memory_deallo(par_memor,'lpoty_inv' ,vacal,lpoty_inv)

  end subroutine generate_comm_bopo
  
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate nodes communicator
  !> @details nenti_ refers to number of elements to consider (outer loop is
  !           through elements).
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_poin(new_part,comm_poin,nenti_,COMM_NAME)

    use def_domain,  only : npoin, nnode
    use def_domain,  only : lnods, nelem, ltype

    implicit none

    integer(ip),         pointer,  intent(in)    :: new_part(:)   
    type(comm_data_par),           intent(inout) :: comm_poin      
    integer(ip),         optional, intent(in)    :: nenti_
    character(*),        optional, intent(in)    :: COMM_NAME

    integer(ip)                                  :: ipoin, ielem, inode
    integer(ip)                                  :: isend, ipart, nsendp
    integer(ip)                                  :: melem, ibuff, iperm
    integer(ip)                                  :: nenti, isin
    integer(ip),         pointer                 :: lelems(:,:)
    integer(ip),         pointer                 :: lnelem(:)
    integer(ip),         pointer                 :: all2send(:)
    integer(ip),         pointer                 :: nperm(:)
    integer(4),          pointer                 :: lnsend(:)
    type(i1p),           pointer                 :: send_perm(:)
    integer(4)                                   :: PAR_COMM_RED_SIZE4
    character(20)                                :: my_comm_name

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM_POIN'
    end if

    call comm_poin % deallo(par_memor,COMM_NAME=trim(my_comm_name))
    call comm_poin % init  (COMM_NAME='COMM_POINT')
    !
    ! Nullify
    !
    nullify(lelems)
    nullify(lnelem)
    nullify(all2send)
    nullify(nperm)
    nullify(lnsend)
    nullify(send_perm)
    !
    ! Check arguments
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti)   call runend("mod_redistribution: wrong size new_part")
    if(nenti /= 0_ip .and. nenti /= nelem) call runend("mod_redistribution: option nenti not encoded")
    !
    ! Generate lnelem and lelems 
    !
    call memory_alloca(par_memor,'lnelem',vacal,lnelem,max(1_ip,npoin))

    lnelem = 0_ip
    melem = 0_ip
    do ielem = 1,nenti
       do ipoin = 1,nnode(ltype(ielem))
          inode = lnods(ipoin,ielem)
          lnelem(inode) = lnelem(inode) + 1_ip 
          if(lnelem(inode) > melem ) melem = lnelem(inode)
       enddo
    enddo

    lnelem = 0_ip
    call memory_alloca(par_memor,'lelems',vacal,lelems,melem,max(1_ip,npoin))
    do ielem = 1,nenti
       do ipoin = 1,nnode(ltype(ielem))
          inode = lnods(ipoin,ielem)
          isin=0_ip;
          do ibuff=1_ip,lnelem(inode)
             if(new_part(lelems(ibuff,inode))==new_part(ielem)) isin = isin + 1_ip;
          enddo
          if(isin == 0_ip) then
             lnelem(inode) = lnelem(inode) + 1_ip
             lelems(lnelem(inode),inode) = ielem 
          else if(isin > 1_ip) then
             call runend("mod_redistribution: something broken isin > 1")
          endif
       enddo
    enddo
    !
    ! Eval lnsend
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)
    do ipoin = 1,npoin
       do ibuff = 1,lnelem(ipoin)
          ielem = lelems(ibuff,ipoin)
          if(new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip) then
             call runend("mod_redistribution: wrong values partition ranks in new_part")
          else
             lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4 
          endif
       enddo
    enddo
    !
    ! Eval comm_poin sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_poin,COMM_NAME=comm_name)
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'ALL2SEND', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'nperm',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'send_perm',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0,PAR_COMM_RED_SIZE-1
       if( lnsend(ipart) > 0_4 ) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'send_perm % l',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do ipoin = 1,npoin
       do ibuff = 1,lnelem(ipoin)
          ielem = lelems(ibuff,ipoin)
          isend = all2send(new_part(ielem))
          nperm(isend) = nperm(isend) + 1_ip
          send_perm(isend) % l( nperm(isend) ) = ipoin
       enddo
    enddo

    iperm=1
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_poin % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,'lnelem'   ,vacal,lnelem)
    call memory_deallo(par_memor,'lelems'   ,vacal,lelems)
    call memory_deallo(par_memor,'lnsend'   ,vacal,lnsend)
    call memory_deallo(par_memor,'send_perm',vacal,send_perm)
    call memory_deallo(par_memor,'all2send' ,vacal,all2send)
    call memory_deallo(par_memor,'nperm'    ,vacal,nperm)

  end subroutine generate_comm_poin

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate elements communicator
  !> @details
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_boun(new_part,comm_boun,nenti_,COMM_NAME)

    use def_domain,  only : lelbo, nboun
    use def_domain,  only : nelem
    use mod_memory      
    implicit none

    integer(ip), pointer,          intent(in)    :: new_part(:)  
    type(comm_data_par),           intent(out)   :: comm_boun
    integer(ip),         optional, intent(in)    :: nenti_
    character(*),        optional, intent(in)    :: COMM_NAME

    integer(ip)                                  :: ielem
    integer(ip)                                  :: isend, ipart, nsendp
    integer(ip)                                  :: ibuff, iperm
    integer(ip)                                  :: iboun, nenti
    integer(ip),         pointer                 :: all2send(:)
    integer(ip),         pointer                 :: nperm(:)
    integer(4),          pointer                 :: lnsend(:)
    type(i1p),           pointer                 :: send_perm(:)
    integer(4)                                   :: PAR_COMM_RED_SIZE4
    character(20)                                :: my_comm_name

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM_BOUN'
    end if
    
    call comm_boun % deallo(par_memor,COMM_NAME=trim(my_comm_name))
    call comm_boun % init(COMM_NAME='COMM_BOUN')
    !
    ! Nullify pointers and communicator 
    !
    nullify(all2send)
    nullify(nperm)
    nullify(lnsend)
    nullify(send_perm)
    !
    ! Check arguments
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti)   call runend("mod_redistribution: wrong size new_part")
    if(nenti /= 0_ip .and. nenti /= nelem) call runend("mod_redistribution: option nenti not encoded")
    !
    ! Generate lnelem and lelems
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)

    do iboun = 1,nboun
       ielem = lelbo(iboun)
       if(new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip) then
          call runend("mod_redistribution: wrong values partition ranks in new_part")
       else
          lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4 
       endif
    end do
    !
    ! Eval comm_poin sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_boun,COMM_NAME=comm_name)      
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'ALL2SEND', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'NPERM',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'SEND_PERM',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0,PAR_COMM_RED_SIZE-1
       if( lnsend(ipart) > 0_4 ) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'SEND_PERM % L',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do iboun = 1,nboun
       ielem        = lelbo(iboun)
       isend        = all2send(new_part(ielem))
       nperm(isend) = nperm(isend) + 1_ip
       send_perm(isend) % l( nperm(isend) ) = iboun
    enddo

    iperm=1_ip
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_boun % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,'LNSEND',   vacal,lnsend)
    call memory_deallo(par_memor,'SEND_PERM',vacal,send_perm)
    call memory_deallo(par_memor,'ALL2SEND', vacal,all2send)
    call memory_deallo(par_memor,'NPERM',    vacal,nperm)

  end subroutine generate_comm_boun
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    29/05/2018 
  !> @brief   Generate communicator sizes (private subroutine)
  !> @details Generate communicator size structures from the send 
  !           the send requirements   
  ! 
  !----------------------------------------------------------------------

  subroutine generate_comm_sizes(lnsend,nsendp,comm,PAR_COMM_IN,COMM_NAME)

    implicit none

    integer(4), pointer,  intent(in)            :: lnsend(:)    
    integer(ip),          intent(out)           :: nsendp 
    type(comm_data_par),  intent(out)           :: comm 
    MY_MPI_COMM   ,       intent(in),  optional :: PAR_COMM_IN
    character(len=*),     intent(in),  optional :: COMM_NAME
    
    integer(ip)          :: ipart, ineig
    integer(4), pointer  :: lnrecv(:)
    integer(4)           :: istat
    integer(4)           :: PAR_MYCOMM_SIZE4
    integer(4)           :: PAR_MYCOMM_RANK4    
    MY_MPI_COMM          :: PAR_MYCOMM4
    character(20)        :: my_comm_name

    if( present(COMM_NAME) ) then
        my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if
    !
    ! Allocate pointers
    !
    nullify(lnrecv)
    if(.not. present(PAR_COMM_IN)) then
       PAR_MYCOMM_SIZE4 = int(PAR_COMM_RED_SIZE,4)
       PAR_MYCOMM4      = PAR_COMM_RED4
    else
       PAR_MYCOMM4 = PAR_COMM_IN
       call PAR_COMM_RANK_AND_SIZE(PAR_MYCOMM4,PAR_MYCOMM_RANK4,PAR_MYCOMM_SIZE4)
    endif
    call memory_alloca(par_memor,'LNRECV',vacal,lnrecv,PAR_MYCOMM_SIZE4,lboun=0_4)
    lnrecv = 0_ip
    nsendp = 0_ip 
    ! 
    ! Communication of communication requirements
    !
#ifndef MPI_OFF      
    call MPI_Alltoall(lnsend,1_4,MPI_INTEGER4,lnrecv,1_4,MPI_INTEGER4,PAR_MYCOMM4,istat)
#endif
    ! 
    ! Allocate and fill send structures (size, dime, neights)
    ! 
    comm  % PAR_COMM_WORLD = PAR_MYCOMM4
    comm  % nneig = 0_ip
    do ipart = 0_ip,PAR_MYCOMM_SIZE4-1
       if(lnrecv(ipart) > 0_ip .or. lnsend(ipart) > 0_ip) then
          comm  % nneig = comm  % nneig + 1_ip
       endif
       if( lnsend(ipart) > 0_ip ) nsendp = nsendp + 1_ip
    end do
    nullify(comm % neights,comm % lsend_size,comm % lrecv_size)
    call memory_alloca(par_memor,trim(comm_name)//' % NEIGHTS'   ,vacal,comm % neights,   max(1_ip, comm % nneig   ))
    call memory_alloca(par_memor,trim(comm_name)//' % LSEND_SIZE',vacal,comm % lsend_size,max(1_ip, comm % nneig+1 ))
    call memory_alloca(par_memor,trim(comm_name)//' % LRECV_SIZE',vacal,comm % lrecv_size,max(1_ip, comm % nneig+1 ))

    comm % lsend_dim = 0_ip
    comm % lrecv_dim = 0_ip
    comm % lsend_size(1_ip) = 1_ip
    comm % lrecv_size(1_ip) = 1_ip

    ineig = 1_ip
    do ipart = 0_ip,PAR_MYCOMM_SIZE4-1
       if( lnrecv(ipart) > 0_ip .or. lnsend(ipart) > 0_ip ) then
          comm % neights(ineig)      = ipart
          comm % lsend_dim           = comm % lsend_dim + lnsend(ipart)
          comm % lsend_size(ineig+1) = comm % lsend_size(ineig) + lnsend(ipart)
          comm % lrecv_dim           = comm % lrecv_dim + lnrecv(ipart)
          comm % lrecv_size(ineig+1) = comm % lrecv_size(ineig) + lnrecv(ipart)
          ineig                      = ineig + 1_ip
       endif
    end do

    nullify( comm % lsend_perm )
    call memory_alloca(par_memor,trim(comm_name)//' % LSEND_PERM',vacal,comm % lsend_perm, comm % lsend_dim)
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,'LNRECV',vacal,lnrecv)

  end subroutine generate_comm_sizes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-30
  !> @brief   ???
  !> @details NENTI = NELEM, NPOIN, NBOUN
  !>          XX(NENTI)
  !>          XX(:,NENTI), XX(NENTI,:)
  !>          XX(:,:,NENTI), XX(:,NENTI,:), XX(NENTI,:,:)
  !>          XX(:,:,:,NENTI), XX(:,:,NENTI,:), XX(:,NENTI,:,:), XX(NENTI,:,:,:)
  !> 
  !-----------------------------------------------------------------------

!!$   subroutine redistribution_array_ip_2(xx,WENTI,posit_opt)
!!$  
!!$     integer(ip),  intent(inout)        :: xx(:,:)
!!$     character(*), intent(in)           :: WENTI
!!$     integer(ip),  intent(in), optional :: posit_opt
!!$     integer(ip)                        :: ndim1,ndim2,posit,ndofn
!!$     integer(ip)                        :: ndim,ldim(10),ii
!!$     integer(ip)                        :: ldim_perm(10)
!!$     !
!!$     ! Find dimensions of XX
!!$     !
!!$     if( associated(xx) ) then
!!$        ndim = size(shape(xx))
!!$        ldim = shape(xx)
!!$     else
!!$        ndim = 0
!!$        ldim = 0
!!$     end if
!!$     call PAR_MAX(ndim)
!!$     call PAR_MAX(ndim,ldim)
!!$     !
!!$     ! Position of entity
!!$     !
!!$     if( present(posit_opt) ) then
!!$        posit = posit_opt
!!$     else
!!$        posit = ndim
!!$     end if
!!$     !
!!$     ! Check dimensions
!!$     !
!!$     if( associated(xx) ) then
!!$        do ii = 1,ndim
!!$           if( ii /= posit .and. size(xx,ii) /= ldim(ii) ) call runend('WRONG DIMENSIONS')
!!$        end do
!!$     end if
!!$
!!$     call redistribution_array(xx,WENTI,posit,ndim,ldim)
!!$     
!!$   end subroutine redistribution_array_ip_2
!!$
!!$   subroutine redistribution_array(xx,WENTI,posit,ndim,ldim)
!!$
!!$     integer(ip),  intent(inout)           :: xx(*)
!!$     character(*), intent(in)              :: WENTI
!!$     integer(ip)                           :: posit
!!$     integer(ip)                           :: ndim
!!$     integer(ip)                           :: ldim(10)
!!$     type(comm_data_par),         pointer  :: COMMD 
!!$     integer(ip)                           :: ndim1,ndim2,posit,ndofn
!!$     integer(ip)                           :: ldim_perm(10),ii
!!$     integer(ip),                 pointer  :: bsend(:,:),brecv(:,:)
!!$     !
!!$     ! => Independent of XX
!!$     !
!!$     nullify(bsend,brecv)
!!$     !
!!$     ! Dimensions
!!$     ! LDIM      = n1,n2,nelem,n3
!!$     ! LDIM_PERM = n1,n2,n3,nelem
!!$     ! NODFN     = n1*n2*n3
!!$     !
!!$     ndofn = 1
!!$     jj    = 0
!!$     do ii = 1,ndim
!!$        if( ii /= posit ) then
!!$           ndofn         = ndofn * ldim(ii)
!!$           jj            = jj + 1
!!$           ldim_perm(jj) = ldim(ii)
!!$        end if
!!$     end do
!!$     ldim_perm(ndim) = ndim
!!$     !
!!$     !
!!$     !
!!$     if(      WENTI == 'NELEM' ) then
!!$        COMMD => commd_nelem
!!$        nenti =  nelem
!!$     else if( WENTI == 'NPOIN' ) then
!!$        COMMD => commd_npoin
!!$        nenti =  npoin
!!$     else if( WENTI == 'NBOUN' ) then
!!$        COMMD => commd_nboun
!!$        nenti =  nboun
!!$     end if
!!$     !
!!$     ! Prepare buffers
!!$     !
!!$     call memory_alloca(par_memor,vacal,'bsend',bsend,ndofn,COMMD % lsend_dim)
!!$     call memory_alloca(par_memor,vacal,'brecv',brecv,ndofn,COMMD % lrecv_dim)
!!$     !
!!$     ! Reshape
!!$     !
!!$     if( posit == ndim ) then
!!$        
!!$     else
!!$     do ii = 1_ip,nenti
!!$        idofn = 0
!!$        do idim = 1,ndim-1
!!$           jdim = ldim_perm(idim) ! Dimension in XX
!!$           do idofi = 1,ldim(idim)
!!$              idofn = idofn + 1       ! Location in BSEND
!!$              !jdofn = 
!!$              bsend(idim,ii) = xx(jj)
!!$           end do
!!$        end do
!!$     end do
!!$     !
!!$     ! Perform redistribution
!!$     !
!!$     call PAR_SEND_RECEIVE_TO_ALL(ndofn,bsend,brecv,COMMD,'SYNCHRONOUS',COMMD % lsend_perm)
!!$     !
!!$     ! Reallocate depends on XX
!!$     !        
!!$     !call memory_deallo(par_memor,vacal,'xx',xx)
!!$     !if(      posit == ndim ) then
!!$     !   call memory_alloca(par_memor,vacal,'xx',xx,ldim(1),COMMD % lrecv_dim)
!!$     !else if( posit == 1    ) then           
!!$     !   call memory_alloca(par_memor,vacal,'xx',xx,COMMD % lrecv_dim,ldim(2))
!!$     !end if
!!$     
!!$   end subroutine redistribution_array_ip
!!$

  subroutine redistribution_domain(lpart,PAR_COMM_IN,VERBOSE)
    use mod_parall,           only : PAR_COMM_MY_CODE
    use mod_htable,           only  : HtableMaxPrimeNumber
    use mod_htable,           only  : hash_t
    use mod_htable,           only  : htaini
    use mod_htable,           only  : htaadd
    use mod_htable,           only  : htades

    integer(ip),        intent(in),           pointer     :: lpart(:)
    MY_MPI_COMM   ,     intent(in), optional              :: PAR_COMM_IN
    logical(lg),        intent(in), optional              :: VERBOSE
    integer(ip)                                           :: ipoin
    integer(ip)                                           :: inode,ii,icont,ielem,ineig
    integer(ip)                                           :: nenti,kpoin,offset,iboun
    integer(ip)                                           :: dom_i,jneig,pnodb,kelem
    integer(ip)                                           :: ifiel,knode,knodb
    MY_MPI_COMM                                           :: PAR_COMM
    logical(lg)                                           :: lfoun
    integer(ip),                              pointer     :: bsend(:,:)
    integer(ip),                              pointer     :: irecv(:,:)
    integer(ip),                              pointer     :: irecv1(:)
    real(rp),                                 pointer     :: rrecv(:,:)
    integer(ip),                              allocatable :: permn(:)
    integer(ip),                              allocatable :: perme(:)
    integer(ip),                              pointer     :: lninv_tmp(:)
    logical(lg)                                           :: if_verbose
    type(hash_t)                                          :: ht
    integer(ip)                                           :: lid
    
    if_verbose = .true.
    if( present(VERBOSE) ) if_verbose = VERBOSE

    if( if_verbose ) call messages_live('REDISTRIBUTION OF DOMAIN ARRAYS','START SECTION')
    !
    ! Compute communicator
    !
    if( present(PAR_COMM_IN) ) then
       PAR_COMM = PAR_COMM_IN
    else
       PAR_COMM = PAR_COMM_MY_CODE
    end if
    call redistribution_comm_initialization(PAR_COMM)
    !
    ! Nullify
    !
    nullify(bsend)
    nullify(irecv)
    nullify(irecv1)
    nullify(rrecv)
    nullify(lninv_tmp)
    !
    ! Domain dimensions
    !
    knode = memory_size(lnods,1_ip)
    knodb = memory_size(lnodb,1_ip)
    call PAR_MAX(knode,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    call PAR_MAX(knodb,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    !
    ! Generate communicator
    !
    nenti = memory_size(lpart) 

    if( nenti /= nelem .and. nenti /= 0 ) call runend('WRONG LPART')

    if( if_verbose ) call messages_live('GENERATE ELEMENT COMMUNICATOR')
    call generate_comm_elem(lpart,commd_nelem,nenti,COMM_NAME='COMMD_NELEM')
    if( if_verbose ) call messages_live('GENERATE NODE COMMUNICATOR')
    call generate_comm_poin(lpart,commd_npoin,nenti,COMM_NAME='COMMD_NPOIN')
    if( if_verbose ) call messages_live('GENERATE BOUNDARY COMMUNICATOR')
    call generate_comm_boun(lpart,commd_nboun,nenti,COMM_NAME='COMMD_NBOUN')      
    !if( if_verbose ) call messages_live('GENERATE BOUNDARY NODE COMMUNICATOR')
    !call generate_comm_bopo(lpart,commd_nbopo,nenti)     
    !
    ! Renumber LNODS/LNODB/LBOEL/LELBO in fragmented local numbering
    !
    if( if_verbose ) call messages_live('RENUMBER ELEMENT/BOUNDARY CONNECTIVITY')
    allocate( permn(max(1_ip,npoin)) )
    permn = -1_ip
    do ineig = 1,commd_npoin % nneig
       icont = 0_ip
       do ii = commd_npoin % lsend_size(ineig),commd_npoin % lsend_size(ineig+1)-1
          ipoin        = commd_npoin % lsend_perm(ii)
          icont        = icont + 1_ip
          permn(ipoin) = icont
       end do
       jneig = 1
       lfoun = .false.
       do while( jneig <= commd_nelem % nneig )
          if( commd_npoin % neights(ineig) == commd_nelem % neights(jneig) ) then
             lfoun = .true.
             exit
          end if
          jneig = jneig + 1
       end do
       if( lfoun ) then
          do ii = commd_nelem  % lsend_size(jneig),commd_nelem % lsend_size(jneig+1)-1
             ielem = commd_nelem % lsend_perm(ii)
             do inode = 1,nnode(abs(ltype(ielem)))
                ipoin              = lnods(inode,ielem)
                kpoin              = permn(ipoin)
                if( kpoin == -1 ) call runend('WE ARE IN BIG TROUBLES 1')
                lnods(inode,ielem) = kpoin
             end do
          end do
       end if
       if( commd_nboun % nneig /= 0 ) then
          jneig = 1
          lfoun = .false.
          do while( jneig <= commd_nboun % nneig )
             if( commd_npoin % neights(ineig) == commd_nboun % neights(jneig) ) then
                lfoun = .true.
                exit
             end if
             jneig = jneig + 1
          end do
          if( lfoun ) then
             do ii = commd_nboun  % lsend_size(jneig),commd_nboun % lsend_size(jneig+1)-1
                iboun = commd_nboun % lsend_perm(ii)
                pnodb = nnode(abs(ltypb(iboun)))
                do inode = 1,pnodb
                   ipoin              = lnodb(inode,iboun)
                   kpoin              = permn(ipoin)
                   if( kpoin == -1 ) call runend('WE ARE IN BIG TROUBLES 2')
                   lnodb(inode,iboun) = kpoin
                end do
             end do
          end if
       end if
       permn = -1_ip
    end do
    !
    ! Renumber LELBO
    !
    allocate(perme(nelem))
    do ineig = 1,commd_nelem % nneig        
       kelem = 0
       do ii = commd_nelem  % lsend_size(ineig),commd_nelem % lsend_size(ineig+1)-1
          ielem = commd_nelem % lsend_perm(ii)
          kelem = kelem + 1
          perme(ielem) = kelem
       end do
    end do
    do iboun = 1,nboun
       ielem = lelbo(iboun)
       kelem = perme(ielem)
       lelbo(iboun) = kelem
    end do
    deallocate( perme )
    deallocate( permn )
    !
    ! New dimensions
    !
    nelem = commd_nelem % lrecv_dim
    nboun = commd_nboun % lrecv_dim
    !
    ! LNODS
    !
    call redistribution_array(lnods,'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LNODS')
    !
    ! LNODB
    !
    call redistribution_array(lnodb,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LNODB')
    !
    ! LBOEL and LELBO
    !
    call redistribution_array(lboel,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBOEL') 
    call redistribution_array(lelbo,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LELBO')
    !
    ! LNINV_TMP
    !
    call messages_live('REDISTRIBUTE GLOBAL NUMBERING')
    call memory_alloca(par_memor,'IRECV1',vacal,irecv1,commd_npoin % lrecv_dim)
    call PAR_SEND_RECEIVE_TO_ALL(lninv_loc,irecv1,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    !
    ! NPOIN and node permutation: unique node numbering from fragmented geometries
    !
    call memory_alloca(par_memor,'COMMD_NPOIN % LRECV_PERM',vacal,commd_npoin % lrecv_perm,commd_npoin % lrecv_dim,REALLOCATE=.true.)         
    npoin = 0_ip
    call htaini( ht, commd_npoin % lrecv_dim+1, lidson=.true., AUTOMATIC_SIZE=.true. )
    do ipoin = 1,commd_npoin % lrecv_dim
       call htaadd( ht, irecv1(ipoin), lid) 
       commd_npoin % lrecv_perm(ipoin) = lid
    end do
    npoin = ht % nelem 
    call memory_deallo(par_memor,'IRECV1',vacal,irecv1)
    !
    ! NELEM and NBOUN permutation: identity! This is required in case of renumbering
    !
    call memory_alloca(par_memor,'COMMD_NELEM % LRECV_PERM',vacal,commd_nelem % lrecv_perm,commd_nelem % lrecv_dim,REALLOCATE=.true.)         
    do ielem = 1,commd_nelem % lrecv_dim
       commd_nelem % lrecv_perm(ielem) = ielem
    end do
    call memory_alloca(par_memor,'COMMD_NBOUN % LRECV_PERM',vacal,commd_nboun % lrecv_perm,commd_nboun % lrecv_dim,REALLOCATE=.true.)         
    do iboun = 1,commd_nboun % lrecv_dim
       commd_nboun % lrecv_perm(iboun) = iboun
    end do
    !
    ! Renumber LNODS
    !
    ielem = 0
    do ineig = 1,commd_nelem % nneig
       offset = commd_npoin % lrecv_size(ineig)-1

       do ii = commd_nelem % lrecv_size(ineig),commd_nelem % lrecv_size(ineig+1)-1
          ielem = ielem + 1
          inode = 0
          loop_inode: do while( inode < knode )
             inode              = inode + 1
             ipoin              = lnods(inode,ielem)
             if( ipoin == 0 ) exit loop_inode
             lnods(inode,ielem) = commd_npoin % lrecv_perm(ipoin+offset)
          end do loop_inode
       end do

    end do
    !
    ! Renumber LNODB
    !
    iboun = 0
    do ineig = 1,commd_nboun % nneig

       dom_i = commd_nboun % neights(ineig)
       jneig = 1
       do while( commd_npoin % neights(jneig) /= dom_i )
          jneig = jneig + 1
       end do
       offset = commd_npoin % lrecv_size(jneig)-1

       do ii = commd_nboun % lrecv_size(ineig),commd_nboun % lrecv_size(ineig+1)-1
          iboun = iboun + 1
          inode = 0
          loop_inodb: do while( inode < knodb )
             inode              = inode + 1
             ipoin              = lnodb(inode,iboun)
             if( ipoin == 0 ) exit loop_inodb
             lnodb(inode,iboun) = commd_npoin % lrecv_perm(ipoin+offset)
          end do loop_inodb
       end do

    end do
    !
    ! Deallocate hash table
    !
    call htades( ht )

    !-------------------------------------------------------------------
    !
    ! NELEM arrays
    !
    !-------------------------------------------------------------------

    call redistribution_array(leinv_loc,'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LEINV_LOC')
    call redistribution_array(ltype,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LTYPE')
    call redistribution_array(lelch,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LELCH')
    call redistribution_array(lesub,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LESUB')
    call redistribution_array(lmate,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LMATE')
    call redistribution_array(leset,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LESET')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NELEM_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call redistribution_array_RP_3(xfiel(ifiel) % a,'NELEM',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A')
          else
             call redistribution_array_RP_3(xfiel(ifiel) % a,'NELEM',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A',EXPLODE_DIM=3_ip)
          end if
       end if
    end do

    !-------------------------------------------------------------------
    !
    ! NPOIN arrays
    !
    !-------------------------------------------------------------------

    call redistribution_array(coord,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='COORD')
    call redistribution_array(lninv_loc,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNINV_LOC')
    call redistribution_array(lnoch,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNOCH')
    call redistribution_array(lmast,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LMAST')
    call redistribution_array(kfl_codno,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='KFL_CODNO')
    call redistribution_array(lnset,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNSET')
    call redistribution_array(lgrou_dom,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LGROU_DOM')

    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call redistribution_array_RP_3(xfiel(ifiel) % a,'NPOIN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A')
          else
             call redistribution_array_RP_3(xfiel(ifiel) % a,'NPOIN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A',EXPLODE_DIM=3_ip)             
          end if
       end if
    end do

    !-------------------------------------------------------------------
    !
    ! NBOUN arrays
    !
    !-------------------------------------------------------------------

    call redistribution_array(ltypb,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LTYPB')
    call redistribution_array(lboch,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBOCH')
    call redistribution_array(lbset,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBSET')
    call redistribution_array(kfl_codbo,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='KFL_CODBO')
    call redistribution_array(lbinv_loc,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBINV_LOC')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
          call redistribution_array_RP_3(xfiel(ifiel) % a,'NBOUN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A')
       end if
    end do
    !
    ! Renumber LELBO
    !
    iboun = 0
    do ineig = 1,commd_nboun % nneig
       offset = 0_ip
       jneig_loop: do jneig = 1,commd_nelem % nneig
          if( commd_nboun % neights(ineig) == commd_nelem % neights(jneig) ) then
             offset = commd_nelem % lrecv_size(jneig)-1
             exit jneig_loop
          endif
       end do jneig_loop
       do ii = commd_nboun % lrecv_size(ineig),commd_nboun % lrecv_size(ineig+1)-1
          iboun = iboun + 1
          kelem = lelbo(iboun)
          lelbo(iboun) = kelem + offset
       end do
    end do
    if(iboun/=nboun) call runend("Something wrong here")

    if( if_verbose ) call messages_live('REDISTRIBUTION OF DOMAIN ARRAYS','END SECTION')

  end subroutine redistribution_domain

  subroutine redistribution_array_IP_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY,VERBOSE)

    integer(ip),                   pointer,  intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    integer(ip),         optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:)
    logical(lg),         optional,           intent(in)    :: VERBOSE
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: nsize
    integer(ip)                                            :: ienti,jenti,nenti_send
    integer(ip)                                            :: ndime_send,ndime_recv
    integer(ip)                                            :: nenti_recv
    logical(lg)                                            :: if_permute_recv
    logical(lg)                                            :: if_verbose
    integer(ip),                   pointer                 :: xsend(:)
    integer(ip),                   pointer                 :: xrecv(:)
    integer(ip),                   pointer                 :: xx_out(:)
    integer(ip),                   pointer                 :: lrecv_perm(:)
    
    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
    
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    
   if( present(VERBOSE) ) then
       if_verbose = VERBOSE
    else
       if_verbose = .true.
    end if
    
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(xrecv)
    nullify(xsend)
    nullify(lrecv_perm)
    !
    ! Check dimensions
    !
    nenti_send = memory_size(xx)
    !
    ! Should the variable be exchanged
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
    else
       call PAR_MAX(nsize,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Redistribute and permute if required
    !
    if( present(variable_name) .and. if_verbose ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if
    !
    ! Dimensions
    !
    if( present(COMM) ) then
       ndime_recv =  COMM % lrecv_dim
       ndime_send =  COMM % lsend_dim
       lrecv_perm => COMM % lrecv_perm
       nenti_recv =  COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv =  commd_nelem % lrecv_dim
       ndime_send =  commd_nelem % lsend_dim
       lrecv_perm => commd_nelem % lrecv_perm
       nenti_recv =  nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv =  commd_npoin % lrecv_dim
       ndime_send =  commd_npoin % lsend_dim
       lrecv_perm => commd_npoin % lrecv_perm
       nenti_recv =  npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv =  commd_nboun % lrecv_dim
       ndime_send =  commd_nboun % lsend_dim
       lrecv_perm => commd_nboun % lrecv_perm
       nenti_recv =  nboun
    else if( trim(wtype) == 'NBOPO' ) then
       ndime_recv =  commd_nbopo % lrecv_dim
       ndime_send =  commd_nbopo % lsend_dim
       lrecv_perm => commd_nbopo % lrecv_perm
       nenti_recv =  nbopo
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN ENTITY TYPE')
    end if

    call memory_alloca(memor_loc,'XSEND',vacal,xsend,nenti_send)
    call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndime_recv)
    !
    ! Put entity at last position
    !
    if( ndime_send > 0 ) then       
       do ienti = 1,nenti_send
          xsend(ienti) = xx(ienti)
       end do       
    end if
    !
    ! Send receive
    !
    if( present(COMM) ) then 
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NELEM' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NPOIN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOUN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOPO' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nbopo,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    end if
    !
    ! Allocate
    !
    if( present(OUTPUT_ARRAY) ) then
       call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
       call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti_recv)
       xx_out => OUTPUT_ARRAY
    else
       call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
       call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv)
       xx_out => xx
    end if
    !
    ! Substitute
    !
    if( associated(lrecv_perm) .and. if_permute_recv ) then
       do ienti = 1,ndime_recv
          jenti = lrecv_perm(ienti)
          xx_out(jenti) = xrecv(ienti)
       end do
    else
       do ienti = 1,ndime_recv
          xx_out(ienti) = xrecv(ienti)
       end do
    end if
    
    call memory_deallo(memor_loc,'XRECV',vacal,xrecv)     
    call memory_deallo(memor_loc,'XSEND',vacal,xsend)
    
    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_IP_1

  subroutine redistribution_array_IP_2(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    integer(ip),                   pointer,  intent(inout) :: xx(:,:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    integer(ip),         optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ndim1,nsize
    integer(ip)                                            :: ienti,jenti,nenti_send,my_posit
    integer(ip)                                            :: ndime_send,ndime_recv
    integer(ip)                                            :: nenti_recv
    logical(lg)                                            :: if_permute_recv
    integer(ip),                   pointer                 :: xsend(:,:)
    integer(ip),                   pointer                 :: xrecv(:,:)
    integer(ip),                   pointer                 :: xx_out(:,:)
    integer(ip),                   pointer                 :: lrecv_perm(:)

    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
     
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(xrecv)
    nullify(xsend)
    nullify(lrecv_perm)
    !
    ! Check dimensions
    !
    ndim1 = 0
    if( present(posit) ) then
       my_posit = posit
    else
       my_posit = 2
    end if
   
    if( my_posit == 1 ) then
       ndim1      = memory_size(xx,2_ip)
       nenti_send = memory_size(xx,1_ip)
    else if( my_posit == 2 ) then
       ndim1      = memory_size(xx,1_ip)
       nenti_send = memory_size(xx,2_ip)
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN POSITION')
    end if
    !
    ! Should the variable be exchanged
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
    else
       call PAR_MAX(nsize,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Redistribute and permute if required
    !
    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if
    !
    ! Dimensions
    !
    if( present(COMM) ) then
       ndime_recv =  COMM % lrecv_dim
       ndime_send =  COMM % lsend_dim
       lrecv_perm => COMM % lrecv_perm
       nenti_recv =  COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv =  commd_nelem % lrecv_dim
       ndime_send =  commd_nelem % lsend_dim
       lrecv_perm => commd_nelem % lrecv_perm
       nenti_recv =  nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv =  commd_npoin % lrecv_dim
       ndime_send =  commd_npoin % lsend_dim
       lrecv_perm => commd_npoin % lrecv_perm
       nenti_recv =  npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv =  commd_nboun % lrecv_dim
       ndime_send =  commd_nboun % lsend_dim
       lrecv_perm => commd_nboun % lrecv_perm
       nenti_recv =  nboun
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN ENTITY TYPE')
    end if

    call memory_alloca(memor_loc,'XSEND',vacal,xsend,ndim1,nenti_send)
    call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndim1,ndime_recv)
    !
    ! Put entity at last position
    !
    if( ndime_send > 0 ) then
       if(      my_posit == 1 ) then
          do ienti = 1,nenti_send
             xsend(1:ndim1,ienti) = xx(ienti,1:ndim1)
          end do
       else if( my_posit == 2 ) then 
          do ienti = 1,nenti_send
             xsend(1:ndim1,ienti) = xx(1:ndim1,ienti)
          end do
       end if
    end if
    !
    ! Send receive
    !
    if( present(COMM) ) then 
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NELEM' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NPOIN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOUN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    end if
    !
    ! Allocate
    !
    if( present(OUTPUT_ARRAY) ) then
       call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
       if(      my_posit == 1 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti_recv,ndim1)
       else if( my_posit == 2 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti_recv)
       end if
       xx_out => OUTPUT_ARRAY
    else
       call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
       if(      my_posit == 1 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv,ndim1)
       else if( my_posit == 2 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,nenti_recv)
       end if
       xx_out => xx
    end if
    !
    ! Substitute
    !
    if( associated(lrecv_perm) .and. if_permute_recv ) then
       if( my_posit == 1 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(jenti,1:ndim1) = xrecv(1:ndim1,ienti)
          end do
       else if( my_posit == 2 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(1:ndim1,jenti) = xrecv(1:ndim1,ienti)
          end do
       end if
    else
       if( my_posit == 1 ) then
          do ienti = 1,ndime_recv
             xx_out(ienti,1:ndim1) = xrecv(1:ndim1,ienti)
          end do
       else if( my_posit == 2 ) then
          do ienti = 1,ndime_recv
             xx_out(1:ndim1,ienti) = xrecv(1:ndim1,ienti)
          end do
       end if
    end if
    
    call memory_deallo(memor_loc,'XRECV',vacal,xrecv)     
    call memory_deallo(memor_loc,'XSEND',vacal,xsend)
    
    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_IP_2

  subroutine redistribution_array_RP_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY,VERBOSE)

    real(rp),                      pointer,  intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    real(rp),            optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:)
    logical(lg),         optional,           intent(in)    :: VERBOSE
    logical(lg)                                            :: if_verbose
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: nsize
    integer(ip)                                            :: ienti,jenti,nenti_send
    integer(ip)                                            :: ndime_send,ndime_recv
    integer(ip)                                            :: nenti_recv
    logical(lg)                                            :: if_permute_recv
    real(rp),                      pointer                 :: xsend(:)
    real(rp),                      pointer                 :: xrecv(:)
    real(rp),                      pointer                 :: xx_out(:)
    integer(ip),                   pointer                 :: lrecv_perm(:)

    nullify(xrecv)
    nullify(xsend)
    nullify(xx_out)
    nullify(lrecv_perm)
    !
    ! Optional arguments
    !
    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
    
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    
    if( present(VERBOSE) ) then
       if_verbose = VERBOSE
    else
       if_verbose = .true.
    end if
    
    if( present(PERMUTE_RECV) ) then
       if_permute_recv = PERMUTE_RECV
    else
       if_permute_recv = .true.
    end if
    !
    ! Check dimensions
    !
    nenti_send = memory_size(xx)
    !
    ! Should the variable be exchanged
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
    else
       call PAR_MAX(nsize,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Redistribute and permute if required
    !
    if( present(variable_name) .and. if_verbose ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if
    !
    ! Dimensions
    !
    if( present(COMM) ) then
       ndime_recv =  COMM % lrecv_dim
       ndime_send =  COMM % lsend_dim
       lrecv_perm => COMM % lrecv_perm
       nenti_recv =  COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv =  commd_nelem % lrecv_dim
       ndime_send =  commd_nelem % lsend_dim
       lrecv_perm => commd_nelem % lrecv_perm
       nenti_recv =  nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv =  commd_npoin % lrecv_dim
       ndime_send =  commd_npoin % lsend_dim
       lrecv_perm => commd_npoin % lrecv_perm
       nenti_recv =  npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv =  commd_nboun % lrecv_dim
       ndime_send =  commd_nboun % lsend_dim
       lrecv_perm => commd_nboun % lrecv_perm
       nenti_recv =  nboun
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN ENTITY TYPE')
    end if

    call memory_alloca(memor_loc,'XSEND',vacal,xsend,nenti_send)
    call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndime_recv)
    !
    ! Put entity at last position
    !
    if( ndime_send > 0 ) then       
       do ienti = 1,nenti_send
          xsend(ienti) = xx(ienti)
       end do       
    end if
    !
    ! Send receive
    !
    if( present(COMM) ) then 
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NELEM' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NPOIN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOUN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    end if
    !
    ! Allocate
    !
    if( present(OUTPUT_ARRAY) ) then
       call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
       call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti_recv)
       xx_out => OUTPUT_ARRAY
    else
       call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
       call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv)
       xx_out => xx
    end if
    !
    ! Substitute
    !
    if( associated(lrecv_perm) .and. if_permute_recv ) then
       do ienti = 1,ndime_recv
          jenti = lrecv_perm(ienti)
          xx_out(jenti) = xrecv(ienti)
       end do
    else
       do ienti = 1,ndime_recv
          xx_out(ienti) = xrecv(ienti)
       end do
    end if
    
    call memory_deallo(memor_loc,'XRECV',vacal,xrecv)     
    call memory_deallo(memor_loc,'XSEND',vacal,xsend)
    
    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_RP_1

  subroutine redistribution_array_RP_2(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY,VARIABLE_TAG,VERBOSE)

    real(rp),                      pointer,  intent(inout) :: xx(:,:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    real(rp),            optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:)
    character(*),        optional,           intent(in)    :: VARIABLE_TAG
    logical(lg),         optional,           intent(in)    :: VERBOSE
    logical(lg)                                            :: if_verbose
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ndim1,nsize
    integer(ip)                                            :: ienti,jenti,nenti_send,my_posit
    integer(ip)                                            :: ndime_send,ndime_recv
    integer(ip)                                            :: nenti_recv
    logical(lg)                                            :: if_permute_recv
    real(rp),                      pointer                 :: xsend(:,:)
    real(rp),                      pointer                 :: xrecv(:,:)
    real(rp),                      pointer                 :: xx_out(:,:)
    integer(ip),                   pointer                 :: lrecv_perm(:)
    
    nullify(xrecv)
    nullify(xsend)
    nullify(lrecv_perm)
    nullify(xx_out)
    !
    ! Optional arguments
    !
    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
    
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    
    if( present(VERBOSE) ) then
       if_verbose = VERBOSE
    else
       if_verbose = .true.
    end if

    if( present(PERMUTE_RECV) ) then
       if_permute_recv = PERMUTE_RECV
    else
       if_permute_recv = .true.   
    end if
    
    if( present(posit) ) then
       my_posit = posit
    else
       my_posit = 2
    end if
    !
    ! Check dimensions
    !
    ndim1 = 0
   
    if( my_posit == 1 ) then
       ndim1      = memory_size(xx,2_ip)
       nenti_send = memory_size(xx,1_ip)
    else if( my_posit == 2 ) then
       ndim1      = memory_size(xx,1_ip)
       nenti_send = memory_size(xx,2_ip)
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN POSITION')
    end if
    !
    ! Should the variable be exchanged
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
    else
       call PAR_MAX(nsize,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Redistribute and permute if required
    !
    if( present(variable_name) .and. if_verbose ) then
       if( present(VARIABLE_TAG) ) then
          call messages_live('REDISTRIBUTE '//trim(variable_name)//'_'//trim(VARIABLE_TAG))
       else
          call messages_live('REDISTRIBUTE '//trim(variable_name))
       end if
    end if
    !
    ! Dimensions
    !
    if( present(COMM) ) then
       ndime_recv =  COMM % lrecv_dim
       ndime_send =  COMM % lsend_dim
       lrecv_perm => COMM % lrecv_perm
       nenti_recv =  COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv =  commd_nelem % lrecv_dim
       ndime_send =  commd_nelem % lsend_dim
       lrecv_perm => commd_nelem % lrecv_perm
       nenti_recv =  nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv =  commd_npoin % lrecv_dim
       ndime_send =  commd_npoin % lsend_dim
       lrecv_perm => commd_npoin % lrecv_perm
       nenti_recv =  npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv =  commd_nboun % lrecv_dim
       ndime_send =  commd_nboun % lsend_dim
       lrecv_perm => commd_nboun % lrecv_perm
       nenti_recv =  nboun
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN ENTITY TYPE')
    end if

    call memory_alloca(memor_loc,'XSEND',vacal,xsend,ndim1,nenti_send)
    call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndim1,ndime_recv)
    !
    ! Put entity at last position
    !
    if( ndime_send > 0 ) then
       if(      my_posit == 1 ) then
          do ienti = 1,nenti_send
             xsend(1:ndim1,ienti) = xx(ienti,1:ndim1)
          end do
       else if( my_posit == 2 ) then 
          do ienti = 1,nenti_send
             xsend(1:ndim1,ienti) = xx(1:ndim1,ienti)
          end do
       end if
    end if
    !
    ! Send receive
    !
    if( present(COMM) ) then 
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NELEM' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NPOIN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOUN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    end if
    !
    ! Allocate
    !
    if( present(OUTPUT_ARRAY) ) then
       call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
       if(      my_posit == 1 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti_recv,ndim1)
       else if( my_posit == 2 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti_recv)
       end if
       xx_out => OUTPUT_ARRAY
    else
       call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
       if(      my_posit == 1 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv,ndim1)
       else if( my_posit == 2 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,nenti_recv)
       end if
       xx_out => xx
    end if
    !
    ! Substitute
    !
    if( associated(lrecv_perm) .and. if_permute_recv ) then
       if( my_posit == 1 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(jenti,1:ndim1) = xrecv(1:ndim1,ienti)
          end do
       else if( my_posit == 2 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(1:ndim1,jenti) = xrecv(1:ndim1,ienti)
          end do
       end if
    else
       if( my_posit == 1 ) then
          do ienti = 1,ndime_recv
             xx_out(ienti,1:ndim1) = xrecv(1:ndim1,ienti)
          end do
       else if( my_posit == 2 ) then
          do ienti = 1,ndime_recv
             xx_out(1:ndim1,ienti) = xrecv(1:ndim1,ienti)
          end do
       end if
    end if
    
    call memory_deallo(memor_loc,'XRECV',vacal,xrecv)     
    call memory_deallo(memor_loc,'XSEND',vacal,xsend)
    
    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_RP_2

  subroutine redistribution_array_RP_3(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY,EXPLODE_DIM,LAST_COMPONENT)

    real(rp),                      pointer,  intent(inout) :: xx(:,:,:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    real(rp),            optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:,:)
    integer(ip),         optional,           intent(in)    :: EXPLODE_DIM
    integer(ip),         optional,           intent(in)    :: LAST_COMPONENT
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ndim1,ndim2,ndim3,ipass
    integer(ip)                                            :: idim1,idim2,idim3,jdim3
    integer(ip)                                            :: ienti,jenti,nenti_send,my_posit
    integer(ip)                                            :: ndime_send,ndime_recv
    integer(ip)                                            :: nenti_recv,nsize
    logical(lg)                                            :: if_permute_recv
    real(rp),                      pointer                 :: xsend(:,:,:)
    real(rp),                      pointer                 :: xrecv(:,:,:)
    real(rp),                      pointer                 :: xx_out(:,:,:)
    real(rp),                      pointer                 :: xx2(:,:)
    real(rp),                      pointer                 :: xx_cpy(:,:,:)
    integer(ip),                   pointer                 :: lrecv_perm(:)
    logical(lg)                                            :: if_explode
    character(100)                                         :: my_variable_name
    integer(ip)                                            :: my_explode_dim
    
    nullify(xrecv)
    nullify(xsend)
    nullify(xx2)
    nullify(xx_cpy)
    nullify(lrecv_perm)
    !
    ! Optional arguments
    !    
    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if

    if( present(EXPLODE_DIM) ) then
       if_explode     = .true.
       my_explode_dim = EXPLODE_DIM
    else
       if_explode     = .false.
    end if

    if( present(PERMUTE_RECV) ) then
       if_permute_recv = PERMUTE_RECV
    else
       if_permute_recv = .true.       
    end if

    if( present(LAST_COMPONENT) ) then
       jdim3 = LAST_COMPONENT
       if_explode = .true.
       my_explode_dim = 3
    else
       jdim3 = 1
    end if    
    !
    ! Check dimensions
    !
    if( present(posit) ) then
       my_posit = posit
    else
       my_posit = 3
    end if

    if(      my_posit == 1 ) then
       ndim1      = memory_size(xx,2_ip)
       ndim2      = memory_size(xx,3_ip)
       nenti_send = memory_size(xx,1_ip)
    else if( my_posit == 2 ) then
       ndim1      = memory_size(xx,1_ip)
       ndim2      = memory_size(xx,3_ip)
       nenti_send = memory_size(xx,2_ip)
    else if( my_posit == 3 ) then
       ndim1      = memory_size(xx,1_ip)
       ndim2      = memory_size(xx,2_ip)
       nenti_send = memory_size(xx,3_ip)
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN POSITION')
    end if
    !
    ! Should the variable be exchanged
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim2,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
    else
       call PAR_MAX(nsize,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim2,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Dimensions
    !
    if( present(COMM) ) then
       ndime_recv =  COMM % lrecv_dim
       ndime_send =  COMM % lsend_dim
       lrecv_perm => COMM % lrecv_perm
       nenti_recv =  COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv =  commd_nelem % lrecv_dim
       ndime_send =  commd_nelem % lsend_dim
       lrecv_perm => commd_nelem % lrecv_perm
       nenti_recv =  nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv =  commd_npoin % lrecv_dim
       ndime_send =  commd_npoin % lsend_dim
       lrecv_perm => commd_npoin % lrecv_perm
       nenti_recv =  npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv =  commd_nboun % lrecv_dim
       ndime_send =  commd_nboun % lsend_dim
       lrecv_perm => commd_nboun % lrecv_perm
       nenti_recv =  nboun
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN ENTITY TYPE')
    end if
    !
    ! Explode one the dimensions, basically to consume less memory
    !
    if( if_explode ) then
       if( my_explode_dim == 3 .and. my_posit /= 3 ) then
          !
          !    NDIM1  NENTI_SEND NDIM2
          ! XX( DIM , ENTITY ,  STEPS  ) => FOR EACH NDIM3, REDISTRIBUTE XX(NDIM,ENTITY)
          !
          ndim1 = ndim1
          ndim3 = ndim2
          ndim2 = nenti_send
          call memory_alloca(memor_loc,'XX_CPY',vacal,xx_cpy,ndim1,ndim2,ndim3)
          do idim3 = jdim3,ndim3
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   xx_cpy(idim1,idim2,idim3) = xx(idim1,idim2,idim3)
                end do
             end do
          end do
          ipass = 0
          do idim3 = jdim3,ndim3
             ipass = ipass + 1
             call memory_alloca(memor_loc,trim(my_variable_name),vacal,xx2,ndim1,ndim2)
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   xx2(idim1,idim2) = xx_cpy(idim1,idim2,idim3)
                end do
             end do
             call redistribution_array(xx2,wtype,POSIT=2_ip,MEMOR=memor_loc,VARIABLE_NAME=trim(my_variable_name),COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)
             if( ipass == 1 ) then
                call memory_deallo(memor_loc,trim(my_variable_name),vacal,xx)
                call memory_alloca(memor_loc,trim(my_variable_name),vacal,xx,ndim1,nenti_recv,ndim3)
             end if
             do idim2 = 1,nenti_recv
                do idim1 = 1,ndim1
                   xx(idim1,idim2,idim3) = xx2(idim1,idim2) 
                end do
             end do
             call memory_deallo(memor_loc,trim(my_variable_name),vacal,xx2)
          end do
          call memory_deallo(memor_loc,'XX_CPY',vacal,xx_cpy)
       else
          call runend('REDISTRIBUTION_ARRAY_RP_3: YOU SHOULD CODE THIS TO USE IT ;o)')
       end if
       
    else
       !
       ! Redistribute and permute if required
       !
       if( present(variable_name) ) then
          call messages_live('REDISTRIBUTE '//trim(variable_name))
       end if

       call memory_alloca(memor_loc,'XSEND',vacal,xsend,ndim1,ndim2,nenti_send)
       call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndim1,ndim2,ndime_recv)
       !
       ! Put entity at last position
       !
       if( ndime_send > 0 ) then
          if(      my_posit == 1 ) then
             do ienti = 1,nenti_send
                xsend(1:ndim1,1:ndim2,ienti) = xx(ienti,1:ndim1,1:ndim2)
             end do
          else if( my_posit == 2 ) then 
             do ienti = 1,nenti_send
                xsend(1:ndim1,1:ndim2,ienti) = xx(1:ndim1,ienti,1:ndim2)
             end do
          else if( my_posit == 3 ) then 
             do ienti = 1,nenti_send
                xsend(1:ndim1,1:ndim2,ienti) = xx(1:ndim1,1:ndim2,ienti)
             end do
          end if
       end if
       !
       ! Send receive
       !
       if( present(COMM) ) then 
          call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
       else if( trim(wtype) == 'NELEM' ) then
          call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       else if( trim(wtype) == 'NPOIN' ) then
          call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       else if( trim(wtype) == 'NBOUN' ) then
          call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       end if
       !
       ! Allocate
       !
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          if(      my_posit == 1 ) then
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti_recv,ndim1,ndim2)
          else if( my_posit == 2 ) then
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti_recv,ndim2)
          else if( my_posit == 3 ) then
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,ndim2,nenti_recv)
          end if
          xx_out => OUTPUT_ARRAY
       else
          call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
          if(      my_posit == 1 ) then
             call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv,ndim1,ndim2)
          else if( my_posit == 2 ) then
             call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,nenti_recv,ndim2)
          else if( my_posit == 3 ) then
             call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,ndim2,nenti_recv)
          end if
          xx_out => xx
       end if
       !
       ! Substitute
       !
       if( associated(lrecv_perm) .and. if_permute_recv ) then
          if( my_posit == 1 ) then
             do ienti = 1,ndime_recv
                jenti = lrecv_perm(ienti)
                xx_out(jenti,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
             end do
          else if( my_posit == 2 ) then
             do ienti = 1,ndime_recv
                jenti = lrecv_perm(ienti)
                xx_out(1:ndim1,jenti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
             end do
          else if( my_posit == 3 ) then
             do ienti = 1,ndime_recv
                jenti = lrecv_perm(ienti)
                xx_out(1:ndim1,1:ndim2,jenti) = xrecv(1:ndim1,1:ndim2,ienti)
             end do
          end if
       else
          if( my_posit == 1 ) then
             do ienti = 1,ndime_recv
                xx_out(ienti,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
             end do
          else if( my_posit == 2 ) then
             do ienti = 1,ndime_recv
                xx_out(1:ndim1,ienti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
             end do
          else if( my_posit == 3 ) then
             do ienti = 1,ndime_recv
                xx_out(1:ndim1,1:ndim2,ienti) = xrecv(1:ndim1,1:ndim2,ienti)
             end do
          end if
       end if

       call memory_deallo(memor_loc,'XRECV',vacal,xrecv)     
       call memory_deallo(memor_loc,'XSEND',vacal,xsend)

    end if

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_RP_3

  subroutine redistribution_array_RP_4(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    real(rp),                      pointer,  intent(inout) :: xx(:,:,:,:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    real(rp),            optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:,:,:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ndim1,nsize
    integer(ip)                                            :: ienti,jenti,nenti_send,my_posit
    integer(ip)                                            :: ndime_send,ndime_recv
    integer(ip)                                            :: ndim2,ndim3,nenti_recv
    logical(lg)                                            :: if_permute_recv
    real(rp),                      pointer                 :: xsend(:,:,:,:)
    real(rp),                      pointer                 :: xrecv(:,:,:,:)
    real(rp),                      pointer                 :: xx_out(:,:,:,:)
    integer(ip),                   pointer                 :: lrecv_perm(:)

    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
    
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(xrecv)
    nullify(xsend)
    nullify(lrecv_perm)
    !
    ! Check dimensions
    !
    ndim1 = 0
    if( present(posit) ) then
       my_posit = posit
    else
       my_posit = 4
    end if
   
    if( my_posit == 1 ) then
       ndim1      = memory_size(xx,2_ip)
       ndim2      = memory_size(xx,3_ip)
       ndim3      = memory_size(xx,4_ip)
       nenti_send = memory_size(xx,1_ip)
    else if( my_posit == 2 ) then
       ndim1      = memory_size(xx,1_ip)
       ndim2      = memory_size(xx,3_ip)
       ndim3      = memory_size(xx,4_ip)
       nenti_send = memory_size(xx,2_ip)
    else if( my_posit == 3 ) then
       ndim1      = memory_size(xx,1_ip)
       ndim2      = memory_size(xx,2_ip)
       ndim3      = memory_size(xx,4_ip)
       nenti_send = memory_size(xx,3_ip)
    else if( my_posit == 4 ) then
       ndim1      = memory_size(xx,1_ip)
       ndim2      = memory_size(xx,2_ip)
       ndim3      = memory_size(xx,3_ip)
       nenti_send = memory_size(xx,4_ip)
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN POSITION')
    end if
    !
    ! Should the variable be exchanged
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim2,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim3,COMM % PAR_COMM_WORLD,INCLUDE_ROOT=.true.)
    else
       call PAR_MAX(nsize,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim1,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim2,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndim3,PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Redistribute and permute if required
    !
    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if
    !
    ! Dimensions
    !
    if( present(COMM) ) then
       ndime_recv =  COMM % lrecv_dim
       ndime_send =  COMM % lsend_dim
       lrecv_perm => COMM % lrecv_perm
       nenti_recv =  COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv =  commd_nelem % lrecv_dim
       ndime_send =  commd_nelem % lsend_dim
       lrecv_perm => commd_nelem % lrecv_perm
       nenti_recv =  nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv =  commd_npoin % lrecv_dim
       ndime_send =  commd_npoin % lsend_dim
       lrecv_perm => commd_npoin % lrecv_perm
       nenti_recv =  npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv =  commd_nboun % lrecv_dim
       ndime_send =  commd_nboun % lsend_dim
       lrecv_perm => commd_nboun % lrecv_perm
       nenti_recv =  nboun
    else
       call runend('MOD_REDISTRIBUTION: UNKNOWN ENTITY TYPE')
    end if

    call memory_alloca(memor_loc,'XSEND',vacal,xsend,ndim1,ndim2,ndim3,nenti_send)
    call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndim1,ndim2,ndim3,ndime_recv)
    !
    ! Put entity at last position
    !
    if( ndime_send > 0 ) then
       if(      my_posit == 1 ) then
          do ienti = 1,nenti_send
             xsend(1:ndim1,1:ndim2,1:ndim3,ienti) = xx(ienti,1:ndim1,1:ndim2,1:ndim3)
          end do
       else if( my_posit == 2 ) then 
          do ienti = 1,nenti_send
             xsend(1:ndim1,1:ndim2,1:ndim3,ienti) = xx(1:ndim1,ienti,1:ndim2,1:ndim3)
          end do
       else if( my_posit == 3 ) then 
          do ienti = 1,nenti_send
             xsend(1:ndim1,1:ndim2,1:ndim3,ienti) = xx(1:ndim1,1:ndim2,ienti,1:ndim3)
          end do
       else if( my_posit == 4 ) then 
          do ienti = 1,nenti_send
             xsend(1:ndim1,1:ndim2,1:ndim3,ienti) = xx(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       end if
    end if
    !
    ! Send receive
    !
    if( present(COMM) ) then 
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NELEM' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NPOIN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOUN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    end if
    !
    ! Allocate
    !
    if( present(OUTPUT_ARRAY) ) then
       call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
       if(      my_posit == 1 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti_recv,ndim1,ndim2,ndim3)
       else if( my_posit == 2 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti_recv,ndim2,ndim3)
       else if( my_posit == 3 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,ndim2,nenti_recv,ndim3)
       else if( my_posit == 4 ) then
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,ndim2,ndim3,nenti_recv)
       end if
       xx_out => OUTPUT_ARRAY
    else
       call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
       if(      my_posit == 1 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv,ndim1,ndim2,ndim3)
       else if( my_posit == 2 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,nenti_recv,ndim2,ndim3)
       else if( my_posit == 3 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,ndim2,nenti_recv,ndim3)
       else if( my_posit == 4 ) then
          call memory_alloca(memor_loc,trim(variable_name),vacal,xx,ndim1,ndim2,ndim3,nenti_recv)
       end if
       xx_out => xx
    end if
    !
    ! Substitute
    !
    if( associated(lrecv_perm) .and. if_permute_recv ) then
       if( my_posit == 1 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(jenti,1:ndim1,1:ndim2,1:ndim3) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       else if( my_posit == 2 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(1:ndim1,jenti,1:ndim2,1:ndim3) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       else if( my_posit == 3 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(1:ndim1,1:ndim2,jenti,1:ndim3) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       else if( my_posit == 4 ) then
          do ienti = 1,ndime_recv
             jenti = lrecv_perm(ienti)
             xx_out(1:ndim1,1:ndim2,1:ndim3,jenti) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       end if
    else
       if( my_posit == 1 ) then
          do ienti = 1,ndime_recv
             xx_out(ienti,1:ndim1,1:ndim2,1:ndim3) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       else if( my_posit == 2 ) then
          do ienti = 1,ndime_recv
             xx_out(1:ndim1,ienti,1:ndim2,1:ndim3) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       else if( my_posit == 3 ) then
          do ienti = 1,ndime_recv
             xx_out(1:ndim1,1:ndim2,ienti,1:ndim3) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       else if( my_posit == 4 ) then
          do ienti = 1,ndime_recv
             xx_out(1:ndim1,1:ndim2,1:ndim3,ienti) = xrecv(1:ndim1,1:ndim2,1:ndim3,ienti)
          end do
       end if
    end if
    
    call memory_deallo(memor_loc,'XRECV',vacal,xrecv)     
    call memory_deallo(memor_loc,'XSEND',vacal,xsend)
    
    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_RP_4

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    8/11/2018 
  !> @brief   Generate comm_data_par to transfer requests from/to an 
  !>          array distributed in chunks
  !> @details
  !>
  !>        Inputs:
  !>            - reqs:     global indices of the requested components
  !>            - chunk:    chunk size
  !>            - dist_dom: SOURCE (defalult) / TARGET
  !>
  !>        Note: when the distributed domain (dist_dom) is the  TARGET,    
  !>              the requests (reqs) are in fact the inputs that will 
  !>              be tranfered to the target domain partitioned into 
  !>              chunks of size "chunk"
  !>
  !>        Outputs:
  !>             - comm:   comm_data_par filled in
  !>
  !----------------------------------------------------------------------------

  subroutine generate_comm_chunkdist_a(reqs,chunk,PAR_COMM4,comm,dist_dom)

    integer(ip),         pointer,  intent(in)    :: reqs(:)    !< Nu,bering (size not given)
    integer(ip),                   intent(in)    :: chunk      !< Chunk size
    MY_MPI_COMM   ,                intent(in)    :: PAR_COMM4  !< Communicator
    type(comm_data_par),           intent(inout) :: comm       !< Communication type
    character(*),        optional, intent(in)    :: dist_dom
    integer(ip)                                  :: nreqs

    nreqs = memory_size(reqs)
    call generate_comm_chunkdist_b(reqs,nreqs,chunk,par_comm4,comm,dist_dom)
    
  end subroutine generate_comm_chunkdist_a

  subroutine generate_comm_chunkdist_b(reqs,nreqs_in,chunk,par_comm4,comm,dist_dom)

    integer(ip),         pointer,  intent(in)    :: reqs(:)   
    integer(ip),                   intent(in)    :: nreqs_in 
    integer(ip),                   intent(in)    :: chunk
    MY_MPI_COMM   ,                intent(in)    :: PAR_COMM4
    type(comm_data_par),           intent(inout) :: comm
    character(*),        optional, intent(in)    :: dist_dom

    logical(lg)                                  :: dist_src
    integer(ip)                                  :: nrank, rank, nneig    
    integer(ip)                                  :: nsend, nrecv,ineig  
    integer(ip)                                  :: ireq,irank,idloc,nreqs
    integer(ip)                                  :: ibuff, iperm
    integer(ip),                   pointer       :: dummi(:)
    integer(ip),                   pointer       :: lnreq_send(:)
    integer(ip),                   pointer       :: lnreq_recv(:)
    type(i1p),                     pointer       :: buff_send(:)
    type(i1p),                     pointer       :: buff_recv(:)
    integer(ip),                   pointer       :: libuf(:)
    integer(ip),                   pointer       :: rank2neig(:)
    integer(ip),                   pointer       :: lsend_perm(:)
    integer(ip),                   pointer       :: lrecv_perm(:)
    integer(ip),                   pointer       :: lrecv_size(:)
    integer(ip)                                  :: mpi_sumsend
    integer(ip)                                  :: mpi_sumrecv
    integer(ip),                   pointer       :: mpi_sendbuf(:)
    integer(ip),                   pointer       :: mpi_recvbuf(:)
    integer(4),                    pointer       :: mpi_sendcounts(:)
    integer(4),                    pointer       :: mpi_recvcounts(:)

    character(100), PARAMETER :: vacal = "generate_comm_chunkdist"
    !
    ! Process inputs
    !
    if( .not. associated(reqs) ) then
       nreqs = 0
    else
       nreqs = nreqs_in
    end if
    
    dist_src = .true.
    if(present(dist_dom))then
       if(dist_dom == 'TARGET') then
          dist_src = .false.
       elseif(dist_dom /= 'SOURCE') then
          call runend('generate_comm_chunkdist: wrong dist_dom argument')
       endif
    endif

    !
    ! Eval rank and nrank and allocate memory
    !
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM4,rank,nrank)
    call comm % init(COMM_NAME='COMM')

    nullify(lnreq_send,lnreq_recv,libuf,buff_send,buff_recv,rank2neig,dummi)
    call memory_alloca(par_memor,'lnreq_send',vacal,lnreq_send,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'lnreq_recv',vacal,lnreq_recv,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'buff_send' ,vacal,buff_send ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'buff_recv' ,vacal,buff_recv ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'libuf'     ,vacal,libuf     ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'rank2neig' ,vacal,rank2neig ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'dummi'     ,vacal,dummi,2_ip,lboun=0_ip)

    !
    !  Generate lists with #requests to send/recv
    !
    do ireq = 1_ip,nreqs
       irank             = min((reqs(ireq)-1_ip)/chunk,nrank-1_ip) 
       lnreq_send(irank) = lnreq_send(irank) + 1_ip
    enddo

    !
    ! Send/recv #requests   
    !
    call  PAR_ALLTOALL(lnreq_send,lnreq_recv,PAR_COMM_IN = PAR_COMM4)

    !
    ! Send/recv requests
    !

    !1) Allocate buffers to receive and send
    do irank = 0,nrank-1
       call memory_alloca(par_memor,"BUFF_SEND % L",vacal,buff_send(irank) % l,max(lnreq_send(irank),1_ip))
       call memory_alloca(par_memor,"BUFF_RECV % L",vacal,buff_recv(irank) % l,max(lnreq_recv(irank),1_ip))
    enddo

    !2) Fill buffers with to send (with receiver local order)
    do ireq = 1_ip,nreqs
       irank                              = min((reqs(ireq)-1_ip)/chunk,nrank-1_ip) 
       idloc                              = mod((reqs(ireq)-1_ip),chunk) + 1_ip
       libuf(irank)                       = libuf(irank) + 1_ip
       buff_send(irank) % l(libuf(irank)) = idloc 
    end do

    !3) Send and recv
    ! isize = 2_ip*nrank
    ! call PAR_START_NON_BLOCKING_COMM(1_ip,isize)
    ! call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    nneig = 0_ip
    do irank = 0,nrank-1
       nsend = lnreq_send(irank)
       nrecv = lnreq_recv(irank)       
       !    if(nrecv > 0_ip) then 
       !       call PAR_SEND_RECEIVE(0_ip,nrecv, dummi,buff_recv(irank) % l,&
       !          & dom_i=irank, wsynch='NON BLOCKING', PAR_COMM_IN=PAR_COMM4)
       !    endif
       if(nsend > 0_ip .or. nrecv > 0_ip) then 
          nneig = nneig+1_ip
          rank2neig(irank) = nneig
       endif
       ! enddo
       ! do irank = 0,nrank-1
       !    nsend = lnreq_send(irank)
       !    if(nsend > 0_ip) then 
       !       call PAR_SEND_RECEIVE(nsend,0_ip, buff_send(irank) % l,dummi,&
       !          & dom_i=irank, wsynch='NON BLOCKING', PAR_COMM_IN=PAR_COMM4)
       !    endif
    end do
    ! call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Fill mpi arrays
    !
    nullify(mpi_sendbuf,mpi_recvbuf,mpi_sendcounts,mpi_recvcounts)
    call memory_alloca(par_memor,'mpi_sendcounts',vacal,mpi_sendcounts,int(nrank,4),lboun=0_4)
    call memory_alloca(par_memor,'mpi_recvcounts',vacal,mpi_recvcounts,int(nrank,4),lboun=0_4)

    mpi_sumsend = 0
    mpi_sumrecv = 0
    do irank = 0,nrank-1
       mpi_sendcounts(irank) = int(lnreq_send(irank),4)
       mpi_recvcounts(irank) = int(lnreq_recv(irank),4)
       mpi_sumsend           = mpi_sumsend + lnreq_send(irank)
       mpi_sumrecv           = mpi_sumrecv + lnreq_recv(irank)
    end do
    call memory_alloca(par_memor,'mpi_sendbuf',vacal,mpi_sendbuf,max(mpi_sumsend,1_ip),lboun=0_ip)
    call memory_alloca(par_memor,'mpi_recvbuf',vacal,mpi_recvbuf,max(mpi_sumrecv,1_ip),lboun=0_ip)
    mpi_sumsend = 0_4
    do irank = 0,nrank-1
       do ibuff = 1,mpi_sendcounts(irank) 
          mpi_sendbuf(mpi_sumsend + ibuff -1) = buff_send(irank) % l(ibuff)
       enddo
       mpi_sumsend = mpi_sumsend + lnreq_send(irank)
    enddo
    !
    ! Perform mpi communications
    !
    call PAR_ALLTOALLV(mpi_sendbuf,mpi_recvbuf,mpi_sendcounts,mpi_recvcounts,PAR_COMM_IN=PAR_COMM4)     
    ! 
    ! Deallocate buffers
    !
    mpi_sumrecv = 0_4
    do irank = 0,nrank-1
       do ibuff = 1,  mpi_recvcounts(irank) 
          buff_recv(irank) % l(ibuff) = mpi_recvbuf(mpi_sumrecv + ibuff -1)
       end do
       mpi_sumrecv = mpi_sumrecv + lnreq_recv(irank)
    end do
    call memory_deallo(par_memor,'mpi_sendcounts',vacal,mpi_sendcounts)
    call memory_deallo(par_memor,'mpi_recvcounts',vacal,mpi_recvcounts)
    call memory_deallo(par_memor,'mpi_sendbuf'   ,vacal,mpi_sendbuf)
    call memory_deallo(par_memor,'mpi_recvbuf'   ,vacal,mpi_recvbuf)
    !
    ! Generate comm_data_par 
    !
    call memory_alloca(par_memor,'COMM % neights'   ,vacal,comm % neights,   max(1_ip, nneig   ))
    call memory_alloca(par_memor,'COMM % lsend_size',vacal,comm % lsend_size,max(1_ip, nneig+1 ))
    call memory_alloca(par_memor,'COMM % lrecv_size',vacal,comm % lrecv_size,max(1_ip, nneig+1 ))

    comm % PAR_COMM_WORLD           = PAR_COMM4
    comm % RANK4                    = int(rank,4)
    comm % lsend_dim                = 0_ip
    comm % lrecv_dim                = 0_ip
    comm % lsend_size(1_ip)         = 1_ip
    comm % lrecv_size(1_ip)         = 1_ip
    comm % nneig                    = nneig

    ineig = 1_ip
    do irank = 0_ip,nrank-1
       if(dist_src) then
          nrecv = lnreq_send(irank) ! you receive the requests sent
          nsend = lnreq_recv(irank) ! you send the requests received
       else
          nsend = lnreq_send(irank) ! inputs sent
          nrecv = lnreq_recv(irank) ! inputs received
       endif

       if( nsend > 0_ip .or. nrecv > 0_ip ) then
          comm % neights(ineig)      = irank
          comm % lsend_dim           = comm % lsend_dim + nsend
          comm % lsend_size(ineig+1) = comm % lsend_size(ineig) + nsend
          comm % lrecv_dim           = comm % lrecv_dim + nrecv
          comm % lrecv_size(ineig+1) = comm % lrecv_size(ineig) + nrecv
          ineig                      = ineig + 1_ip
       endif
    end do

    call memory_alloca(par_memor,'COMM % LSEND_PERM',vacal,comm % lsend_perm, comm % lsend_dim)
    call memory_alloca(par_memor,'COMM % LRECV_PERM',vacal,comm % lrecv_perm, comm % lrecv_dim)

    nullify(lsend_perm,lrecv_perm,lrecv_size)

    if(dist_src) then
       lsend_perm => comm % lsend_perm
       lrecv_perm => comm % lrecv_perm
       lrecv_size => comm % lrecv_size
    else
       lsend_perm => comm % lrecv_perm
       lrecv_perm => comm % lsend_perm
       lrecv_size => comm % lsend_size
    endif

    !lsend_perm
    iperm = 1_ip
    do irank = 0_ip,nrank-1
       do ibuff = 1_ip, lnreq_recv(irank)
          lsend_perm(iperm)=buff_recv(irank) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    end do

    !lrecv_perm
    libuf(:)=0_ip
    do ireq = 1_ip,nreqs
       irank = min((reqs(ireq)-1_ip)/chunk,nrank-1_ip) 
       ineig = rank2neig(irank)
       lrecv_perm(lrecv_size(ineig)+libuf(irank)) = ireq
       libuf(irank) = libuf(irank) + 1_ip 
    enddo

    nullify(lsend_perm,lrecv_perm,lrecv_size) 
    !
    ! Deallocate memory
    !
    call memory_deallo(par_memor,'lnreq_send',vacal,lnreq_send)     
    call memory_deallo(par_memor,'lnreq_recv',vacal,lnreq_recv)     
    call memory_deallo(par_memor,'buff_send', vacal,buff_send)
    call memory_deallo(par_memor,'buff_recv', vacal,buff_recv)
    call memory_deallo(par_memor,'libuf',     vacal,libuf)
    call memory_deallo(par_memor,'rank2neig', vacal,rank2neig)
    call memory_deallo(par_memor,'dummi'    , vacal,dummi)

  end subroutine generate_comm_chunkdist_b

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    13/11/2018 
  !> @brief   Generate a hash_chunkdist in order to perform permuations from
  !>          gloval to local id of the nodes. The global id (gid) is the 
  !>          global identification of the node, the local id (lid) is its
  !>          position in the list of requests (reqs) used in the subroutione: 
  !>          generate_comm_chunkdist and sent here also as an argument.
  !>  
  !> @details
  !>
  !>        Inputs:
  !>            - reqs:     global indices of the requested components
  !>            - chunk:    chunk size
  !>            - comm:     comm_data_par evaluated in generate_comm_chunkdist
  !>
  !>        Outputs:
  !>             - htable: hash_chunkdist for permutation
  !>
  !----------------------------------------------------------------------------
  subroutine redistribution_generate_hash_chunkdist(reqs,chunk,comm,htable)

    implicit none

    integer(ip),          pointer, intent(in)    :: reqs(:)   
    integer(ip),                   intent(in)    :: chunk
    type(comm_data_par),           intent(inout) :: comm
    type(hash_chunkdist),          intent(inout) :: htable

    integer(ip)                                  :: nrank, rank
    integer(ip)                                  :: ineig,neigh
    integer(ip)                                  :: irecv,nrecv
    integer(ip)                                  :: iperm,ilid
    integer(ip)                                  :: ireq
    integer(ip)                                  :: mingid,maxgid

    character(100), PARAMETER :: vacal = "redistribution_generate_hash_chunkdist"

    call PAR_COMM_RANK_AND_SIZE( comm % PAR_COMM_WORLD,rank,nrank)
    nullify(htable % offset,htable % perm)
    call memory_alloca(par_memor,'htable % offset',vacal,htable % offset, nrank,lboun=0_ip)
    call memory_alloca(par_memor,'htable % perm'  ,vacal,htable % perm,   nrank,lboun=0_ip)

    htable % chunk  = chunk
    htable % nchunk = nrank

    ilid = 1_ip
    do ineig = 1_ip, comm % nneig

       neigh  = comm % neights(ineig)

       mingid = huge(1_ip)
       maxgid = -1_ip
       nrecv  =  comm % lrecv_size(ineig+1_ip)-comm % lrecv_size(ineig)  

       if(nrecv > 0_ip) then

          !
          ! Eval range within each block
          !
          do irecv = comm % lrecv_size(ineig) , comm % lrecv_size(ineig+1_ip) - 1_ip
             ireq = reqs(comm % lrecv_perm(irecv))
             if(ireq < mingid) mingid = ireq
             if(ireq > maxgid) maxgid = ireq
          enddo

          !
          ! Store offset and allocte local permutation array
          !
          htable % offset(neigh) = mingid - 1_ip
          call memory_alloca(par_memor,'hatable % perm % l',vacal,htable % perm(neigh) % l,maxgid-mingid+1_ip)

          !
          ! Eval local permutaiton array
          !
          do irecv = comm % lrecv_size(ineig) , comm % lrecv_size(ineig+1_ip) - 1_ip
             ireq = reqs(comm % lrecv_perm(irecv))
             htable % perm(neigh) % l(ireq-mingid+1) = 1_ip
          enddo
          do iperm = 1_ip,maxgid-mingid+1_ip
             if(htable % perm(neigh) % l(iperm) == 1_ip) then
                htable % perm(neigh) % l(iperm) = ilid
                ilid = ilid + 1_ip 
             endif
          enddo
       endif

    enddo

  end subroutine redistribution_generate_hash_chunkdist

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    13/11/2018 
  !> @brief   Deallocate a hash_chunkdist (htable)
  !> @details Deallocate pointers within the type hash_chunkdist
  !>
  !----------------------------------------------------------------------------
  subroutine redistribution_deallocate_hash(htable)

    type(hash_chunkdist), intent(inout) :: htable

    character(100), PARAMETER :: vacal = "redistribution_deallocate_hash"

    call memory_deallo(par_memor,'htable % offset',vacal,htable % offset)
    call memory_deallo(par_memor,'htable % perm'  ,vacal,htable % perm)

  end subroutine redistribution_deallocate_hash

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    13/11/2018 
  !> @brief   Permute from global id (gid) to local id (lid) using a  
  !>          hash_chunkdist
  !> @details
  !>
  !>        Inputs:
  !>            - gid: global id of the node
  !>
  !>        Outputs:
  !>             - local id (position of the gid in reqs)
  !>
  !----------------------------------------------------------------------------
  integer(ip) function redistribution_lid_chunkdist(gid,htable)

    integer(ip),                   intent(in)   :: gid 
    type(hash_chunkdist),          intent(in)   :: htable

    integer(ip)                                 :: ichunk

    ichunk = min((gid-1_ip)/htable % chunk,htable % nchunk-1_ip) 
    redistribution_lid_chunkdist = htable % perm(ichunk) % l(gid-htable % offset(ichunk))

  end function redistribution_lid_chunkdist

  subroutine redistribution_array_R1P_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    type(r1p),           pointer,            intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    type(r1p),           optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ienti,nenti_send
    integer(ip)                                            :: nenti_recv
    integer(ip)                                            :: idim1
    integer(ip)                                            :: kdime_max(1),kdime_min(1)
    integer(ip)                                            :: all_equal
    integer(ip),                   pointer                 :: kdime(:,:)
    real(rp),                      pointer                 :: xx_tmp(:,:)
 
    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
    
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    nullify(kdime)
    nullify(xx_tmp)
    !
    ! Check if all dimensions are equal... if this is the case, no need to exchange them
    ! so that we can save time
    !
    nenti_send = memory_size(xx)
    call memory_alloca(memor_loc,'KDIME',vacal,kdime,1_ip,nenti_send)

    all_equal =  1
    kdime_min =  huge(1_ip)
    kdime_max = -huge(1_ip)
    do ienti = 1,nenti_send
       if( associated(xx(ienti) % a) ) then
          kdime(1,ienti) = size(xx(ienti) % a,DIM=1,KIND=ip)
          kdime_min(1)   = min(kdime_min(1),kdime(1,ienti))
          kdime_max(1)   = max(kdime_max(1),kdime(1,ienti))
       end if
    end do
    call PAR_MAX(1_ip,kdime_max)
    if( kdime_max(1) <= 0 ) goto 10

    call PAR_MIN(1_ip,kdime_min)
    
    if( kdime_max(1) /= kdime_min(1) ) all_equal = 0

    if( present(OUTPUT_ARRAY) ) call runend('THIS OPTION DOES NOT EXIST WITH OUTPUT ARRAY')
    !
    ! Copy XX to temporary array XX_TMP
    !
    call memory_alloca(memor_loc,trim(variable_name),vacal,xx_tmp,kdime_max(1),nenti_send)       
    do ienti = 1,nenti_send
       do idim1 = 1,kdime(1,ienti)
          xx_tmp(idim1,ienti) = xx(ienti) % a(idim1)
       end do
    end do
    !
    ! If dimensions are not the same, then exchange them!
    !
    if( all_equal == 0 ) then
       call redistribution_array(kdime,WTYPE=wtype,POSIT=2_ip,MEMOR=memor_loc,VARIABLE_NAME='KDIME',COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)
    end if
    !
    ! Redistribute XX_TMP
    !
    call redistribution_array(xx_tmp,WTYPE=wtype,POSIT=2_ip,MEMOR=memor_loc,VARIABLE_NAME=variable_name,COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)
    nenti_recv = memory_size(xx_tmp,2_ip)
    !
    ! Copy XX_TMP into XX
    !
    call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
    call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv)

    if( all_equal == 1 ) then
       do ienti = 1,nenti_recv
          call memory_alloca(memor_loc,trim(variable_name)//' % A',vacal,xx(ienti)%a,kdime_max(1))  
          do idim1 = 1,kdime_max(1)
             xx(ienti) % a(idim1) = xx_tmp(idim1,ienti) 
          end do
       end do
    else
       do ienti = 1,nenti_recv
          call memory_alloca(memor_loc,trim(variable_name)//' % A',vacal,xx(ienti)%a,kdime(1,ienti))      
          do idim1 = 1,kdime(1,ienti)
             xx(ienti) % a(idim1) = xx_tmp(idim1,ienti) 
          end do
       end do
    end if
    !
    ! Deallocate memory
    !
10  continue
    call memory_deallo(memor_loc,trim(variable_name),vacal,xx_tmp)
    call memory_deallo(memor_loc,'KDIME',vacal,kdime)

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_R1P_1

  subroutine redistribution_array_R3P_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    type(r3p),           pointer,            intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    type(r3p),           optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:,:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ienti,nenti_send
    integer(ip)                                            :: nenti_recv
    integer(ip)                                            :: idim1,idim2,idim3
    integer(ip)                                            :: kdime_max(3),kdime_min(3)
    integer(ip)                                            :: all_equal
    integer(ip),                   pointer                 :: kdime(:,:)
    real(rp),                      pointer                 :: xx_tmp(:,:,:,:)

    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    nullify(kdime)
    nullify(xx_tmp)
    !
    ! Check if all dimensions are equal... if this is the case, no need to exchange them
    ! so that we can save time
    !
    nenti_send = memory_size(xx)
    call memory_alloca(memor_loc,'KDIME',vacal,kdime,3_ip,nenti_send)

    all_equal =  1
    kdime_min =  huge(1_ip)
    kdime_max = -huge(1_ip)
    do ienti = 1,nenti_send
       if( associated(xx(ienti) % a) ) then
          kdime(1,ienti) = size(xx(ienti) % a,DIM=1,KIND=ip)
          kdime(2,ienti) = size(xx(ienti) % a,DIM=2,KIND=ip)
          kdime(3,ienti) = size(xx(ienti) % a,DIM=3,KIND=ip)
          kdime_min(1)   = min(kdime_min(1),kdime(1,ienti))
          kdime_min(2)   = min(kdime_min(2),kdime(2,ienti))
          kdime_min(3)   = min(kdime_min(3),kdime(3,ienti))
          kdime_max(1)   = max(kdime_max(1),kdime(1,ienti))
          kdime_max(2)   = max(kdime_max(2),kdime(2,ienti))
          kdime_max(3)   = max(kdime_max(3),kdime(3,ienti))
       end if
    end do
    call PAR_MAX(3_ip,kdime_max)
    if( kdime_max(1) <= 0 .or. kdime_max(2) <= 0 .or. kdime_max(3) <= 0 ) goto 10

    call PAR_MIN(3_ip,kdime_min)
    
    if( kdime_max(1) /= kdime_min(1) ) all_equal = 0
    if( kdime_max(2) /= kdime_min(2) ) all_equal = 0
    if( kdime_max(3) /= kdime_min(3) ) all_equal = 0

    if( present(OUTPUT_ARRAY) ) call runend('THIS OPTION DOES NOT EXIST WITH OUTPUT ARRAY')
    !
    ! Copy XX to temporary array XX_TMP
    !
    call memory_alloca(memor_loc,trim(variable_name),vacal,xx_tmp,kdime_max(1),kdime_max(2),kdime_max(3),nenti_send)       
    do ienti = 1,nenti_send
       do idim3 = 1,kdime(3,ienti)
          do idim2 = 1,kdime(2,ienti)
             do idim1 = 1,kdime(1,ienti)
                xx_tmp(idim1,idim2,idim3,ienti) = xx(ienti) % a(idim1,idim2,idim3)
             end do
          end do
       end do
    end do
    !
    ! If dimensions are not the same, then exchange them!
    !
    if( all_equal == 0 ) then
       call redistribution_array(kdime,WTYPE=wtype,POSIT=2_ip,MEMOR=memor_loc,VARIABLE_NAME='KDIME',COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)
    end if
    !
    ! Redistribute XX_TMP
    !
    call redistribution_array(xx_tmp,WTYPE=wtype,POSIT=4_ip,MEMOR=memor_loc,VARIABLE_NAME=variable_name,COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)
    nenti_recv = memory_size(xx_tmp,4_ip)
    !
    ! Copy XX_TMP into XX
    !
    call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
    call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv)

    if( all_equal == 1 ) then
       do ienti = 1,nenti_recv
          call memory_alloca(memor_loc,trim(variable_name)//' % A',vacal,xx(ienti)%a,kdime_max(1),kdime_max(2),kdime_max(3))      
          do idim3 = 1,kdime_max(3)
             do idim2 = 1,kdime_max(2)
                do idim1 = 1,kdime_max(1)
                   xx(ienti) % a(idim1,idim2,idim3) = xx_tmp(idim1,idim2,idim3,ienti) 
                end do
             end do
          end do
       end do
    else
       do ienti = 1,nenti_recv
          call memory_alloca(memor_loc,trim(variable_name)//' % A',vacal,xx(ienti)%a,kdime(1,ienti),kdime(2,ienti),kdime(3,ienti))      
          do idim3 = 1,kdime(3,ienti)
             do idim2 = 1,kdime(2,ienti)
                do idim1 = 1,kdime(1,ienti)
                   xx(ienti) % a(idim1,idim2,idim3) = xx_tmp(idim1,idim2,idim3,ienti) 
                end do
             end do
          end do
       end do
    end if
    !
    ! Deallocate memory
    !
10  continue
    call memory_deallo(memor_loc,trim(variable_name),vacal,xx_TMP)
    call memory_deallo(memor_loc,'KDIME',vacal,kdime)

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_R3P_1

  subroutine redistribution_array_R2P_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    type(r2p),           pointer,            intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    type(r2p),           optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ienti,nenti_send
    integer(ip)                                            :: nenti_recv
    integer(ip)                                            :: idim1,idim2
    integer(ip)                                            :: kdime_max(2),kdime_min(2)
    integer(ip)                                            :: all_equal
    integer(ip),                   pointer                 :: kdime(:,:)
    real(rp),                      pointer                 :: xx_tmp(:,:,:)

    if( present(variable_name) ) then
       my_variable_name = trim(variable_name)
    else
       my_variable_name = 'XX'
    end if
    
    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    nullify(kdime)
    nullify(xx_tmp)
    !
    ! Check if all dimensions are equal... if this is the case, no need to exchange them
    ! so that we can save time
    !
    nenti_send = memory_size(xx)
    call memory_alloca(memor_loc,'KDIME',vacal,kdime,2_ip,nenti_send)

    all_equal =  1
    kdime_min =  huge(1_ip)
    kdime_max = -huge(1_ip)
    do ienti = 1,nenti_send
       if( associated(xx(ienti) % a) ) then
          kdime(1,ienti) = size(xx(ienti) % a,DIM=1,KIND=ip)
          kdime(2,ienti) = size(xx(ienti) % a,DIM=2,KIND=ip)
          kdime_min(1)   = min(kdime_min(1),kdime(1,ienti))
          kdime_min(2)   = min(kdime_min(2),kdime(2,ienti))
          kdime_max(1)   = max(kdime_max(1),kdime(1,ienti))
          kdime_max(2)   = max(kdime_max(2),kdime(2,ienti))
       end if
    end do
    call PAR_MAX(2_ip,kdime_max)
    if( kdime_max(1) <= 0 .or. kdime_max(2) <= 0 ) goto 10

    call PAR_MIN(2_ip,kdime_min)
    
    if( kdime_max(1) /= kdime_min(1) ) all_equal = 0
    if( kdime_max(2) /= kdime_min(2) ) all_equal = 0

    if( present(OUTPUT_ARRAY) ) call runend('THIS OPTION DOES NOT EXIST WITH OUTPUT ARRAY')
    !
    ! Copy XX to temporary array XX_TMP
    !
    call memory_alloca(memor_loc,trim(variable_name),vacal,xx_tmp,kdime_max(1),kdime_max(2),nenti_send)       
    do ienti = 1,nenti_send
       do idim2 = 1,kdime(2,ienti)
          do idim1 = 1,kdime(1,ienti)
             xx_tmp(idim1,idim2,ienti) = xx(ienti) % a(idim1,idim2)
          end do
       end do
    end do
    !
    ! If dimensions are not the same, then exchange them!
    !
    if( all_equal == 0 ) then
       call redistribution_array(kdime,WTYPE=wtype,POSIT=2_ip,MEMOR=memor_loc,VARIABLE_NAME='KDIME',COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)!,OUTPUT_ARRAY)
    end if
    !
    ! Redistribute XX_TMP
    !
    call redistribution_array(xx_tmp,WTYPE=wtype,POSIT=3_ip,MEMOR=memor_loc,VARIABLE_NAME=variable_name,COMM=COMM,PERMUTE_RECV=PERMUTE_RECV)!,OUTPUT_ARRAY)
    nenti_recv = memory_size(xx_tmp,3_ip)
    !
    ! Copy XX_TMP into XX
    !
    call memory_deallo(memor_loc,trim(variable_name),vacal,xx)
    call memory_alloca(memor_loc,trim(variable_name),vacal,xx,nenti_recv)

    if( all_equal == 1 ) then
       do ienti = 1,nenti_recv
          call memory_alloca(memor_loc,trim(variable_name)//' % A',vacal,xx(ienti)%a,kdime_max(1),kdime_max(2))      
          do idim2 = 1,kdime_max(2)
             do idim1 = 1,kdime_max(1)
                xx(ienti) % a(idim1,idim2) = xx_tmp(idim1,idim2,ienti) 
             end do
          end do
       end do
    else
       do ienti = 1,nenti_recv
          call memory_alloca(memor_loc,trim(variable_name)//' % A',vacal,xx(ienti)%a,kdime(1,ienti),kdime(2,ienti))      
          do idim2 = 1,kdime(2,ienti)
             do idim1 = 1,kdime(1,ienti)
                xx(ienti) % a(idim1,idim2) = xx_tmp(idim1,idim2,ienti) 
             end do
          end do
       end do
    end if
    !
    ! Deallocate memory
    !
10  continue
    call memory_deallo(memor_loc,trim(variable_name),vacal,xx_TMP)
    call memory_deallo(memor_loc,'KDIME',vacal,kdime)

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_R2P_1
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-07
  !> @brief   Renumber node communicator for redistribution
  !> @details Renumber node communicator for redistribution
  !> 
  !-----------------------------------------------------------------------

  subroutine redistribution_renumbering_npoin_communicator(permr)
    
    integer(ip), intent(in), pointer :: permr(:)
    integer(ip)                      :: ii,ipoin,jpoin

    do ii = 1,commd_npoin % lrecv_dim
       ipoin = commd_npoin % lrecv_perm(ii)
       jpoin = permr(ipoin)
       commd_npoin % lrecv_perm(ii) = jpoin
    end do

  end subroutine redistribution_renumbering_npoin_communicator

  subroutine redistribution_renumbering_nboun_communicator(permr)
    
    integer(ip), intent(in), pointer :: permr(:)
    integer(ip)                      :: ii,iboun,jboun

    do ii = 1,commd_nboun % lrecv_dim
       iboun = commd_nboun % lrecv_perm(ii)
       jboun = permr(iboun)
       commd_nboun % lrecv_perm(ii) = jboun
    end do

  end subroutine redistribution_renumbering_nboun_communicator

  subroutine redistribution_renumbering_nelem_communicator(permr)
    
    integer(ip), intent(in), pointer :: permr(:)
    integer(ip)                      :: ii,ielem,jelem

    do ii = 1,commd_nelem % lrecv_dim
       ielem = commd_nelem % lrecv_perm(ii)
       jelem = permr(ielem)
       commd_nelem % lrecv_perm(ii) = jelem
    end do

  end subroutine redistribution_renumbering_nelem_communicator

end module mod_redistribution
!> @}
