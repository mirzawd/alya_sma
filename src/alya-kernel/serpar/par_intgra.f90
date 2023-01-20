!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_intgra(kpoin,permR,ia)
  !-----------------------------------------------------------------------
  !****f* domain/par_intgra
  ! NAME
  !    domain
  ! DESCRIPTION
  !    Compute interior nodes graph
  ! OUTPUT
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_domain
  use def_master
  use def_parall
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_parall, only : par_memor
  implicit none 
  integer(ip), intent(in)  :: kpoin
  integer(ip), intent(in)  :: permR(*)
  integer(ip), intent(out) :: ia(*)
  integer(ip), pointer     :: nepoi_tmp(:)
  integer(ip), pointer     :: lelpo_tmp(:)
  integer(ip), pointer     :: pelpo_tmp(:)
  integer(ip), pointer     :: lista(:)
  integer(ip)              :: ielem,inode,ipoin,jelem,mpopo_tmp,lsize
  integer(ip)              :: pnode,jpoin,nzdom_tmp

  nullify(nepoi_tmp)
  nullify(lelpo_tmp)
  nullify(pelpo_tmp)
  nullify(lista)

  !----------------------------------------------------------------------
  !
  ! NELPO_tmp, PELPO_tmp
  !
  !----------------------------------------------------------------------
  !
  ! Allocate memory for NEPOI and compute it
  ! 
  call memory_alloca(par_memor,'NEPOI_TMP','par_intgra',nepoi_tmp,npoin)

  do ielem = 1,nelem
     do inode = 1,nnode(ltype(ielem)) ! LNNOD STUFF
     !do inode = 1,lnnod(ielem) 
        ipoin = lnods(inode,ielem)
        nepoi_tmp(ipoin) = nepoi_tmp(ipoin) + 1
     end do
  end do
  !
  ! Allocate memory for PELPO_tmp and compute it
  !
  call memory_alloca(par_memor,'PELPO_TMP','par_intgra',pelpo_tmp,npoin+1_ip)
  pelpo_tmp(1) = 1
  do ipoin = 1,npoin
     pelpo_tmp(ipoin+1) = pelpo_tmp(ipoin) + nepoi_tmp(ipoin)
  end do
  !
  ! Allocate memory for LELPO_tmp and construct the list
  !
  call memory_alloca(par_memor,'LELPO_TMP','par_intgra',lelpo_tmp,pelpo_tmp(npoin+1))
  mpopo_tmp = 0
  do ielem = 1,nelem
     pnode = nnode(ltype(ielem)) ! LNNOD STUFF
     !pnode = lnnod(ielem)
     mpopo_tmp = mpopo_tmp + pnode * pnode
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        lelpo_tmp(pelpo_tmp(ipoin)) = ielem
        pelpo_tmp(ipoin) = pelpo_tmp(ipoin)+1
     end do
  end do
  !
  ! Recompute PELPO_tmp 
  !
  pelpo_tmp(1) = 1
  do ipoin = 1,npoin
     pelpo_tmp(ipoin+1) = pelpo_tmp(ipoin) + nepoi_tmp(ipoin)
  end do
  !
  ! Deallocate memory for temporary node/element connectivity
  !
  call memory_deallo(par_memor,'NEPOI_TMP','par_intgra',nepoi_tmp)

  !----------------------------------------------------------------------
  !
  ! IA, JA (in GISCA): graph does not include own node
  !
  !----------------------------------------------------------------------

  call memory_alloca(par_memor,'LISTA','par_intgra',lista,mpopo_tmp)
  !
  ! Construct the array of indexes
  !     
  ia(1) = 1
 
  do ipoin = 1,npoin
     lsize = 0
     jpoin = permR(ipoin)
     if( jpoin > 0 ) then
        do ielem = pelpo_tmp(ipoin),pelpo_tmp(ipoin+1)-1
           jelem = lelpo_tmp(ielem)
           pnode = nnode(ltype(jelem)) ! LNNOD STUFF
           !pnode = lnnod(jelem)
           call par_mergli( lista(ia(jpoin):), lsize, lnods(:,jelem), &
                pnode, jpoin , permR )
        end do
        ia(jpoin+1) = ia(jpoin) + lsize
     end if
  end do
 
  nzdom_tmp = ia(kpoin+1)-1
  call memgen(1_ip,max(1_ip,nzdom_tmp),0_ip) ! A minimu size is required
  do ipoin = 1,ia(kpoin+1)-1
     gisca(ipoin) = lista(ipoin)
  end do
  !
  ! Deallocate memory
  !
  call memory_deallo(par_memor,'LISTA'    ,'par_intgra',lista)
  call memory_deallo(par_memor,'PELPO_TMP','par_intgra',pelpo_tmp)
  call memory_deallo(par_memor,'LELPO_TMP','par_intgra',lelpo_tmp)

end subroutine par_intgra

subroutine par_mergli( lista, lsize, lnods, pnode, me , permR )
  !-----------------------------------------------------------------------
  !****f* Domain/mergli
  ! NAME
  !    mergli
  ! DESCRIPTION
  !    This routine merges to list of lnods
  ! OUTPUT
  !    LISTA
  !    LSIZE
  ! USED BY
  !    par_intgra
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,lg
  implicit none
  integer(ip), intent(inout) :: lsize,lista(*)
  integer(ip), intent(in)    :: pnode,me
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip), intent(in)    :: permR(*)
  integer(ip)                :: inode, jnode, n1, n2
  logical(lg)                :: noEncontrado

  do inode = 1,pnode
     n1 = permR(lnods(inode))
     if( n1 > 0 .and. n1 /= me ) then
        jnode = 1
        noEncontrado = .true.
        do while( jnode <= lsize .and. noEncontrado )
           n2 = lista(jnode)
           if( n1 == n2 .or. n2 <= 0 ) then
              noEncontrado = .false.
           end if
           jnode = jnode + 1
        end do
        if( noEncontrado ) then
           lsize = lsize + 1
           lista(lsize) = n1
        end if
     end if
  end do

end subroutine par_mergli
