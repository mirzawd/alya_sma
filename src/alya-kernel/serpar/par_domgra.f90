!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_domgra.f90
!> @author  Guillaume Houzeaux
!> @brief   Parallel subdomain arrays
!> @details Compute some parallel subdomain arrays.
!!          \verbatim
!!          NEIGHDOM(I,J)... I =J : total # boundary nodes (including 
!!                                  repeated communications)
!!                           I/=J : # boundary node shared between I and J
!!          LNPAR_PAR(I) ... Subdomain # if I=internal node
!!                       ... 0           if I=boundary node
!!          LNEIG_PAR(I) ... Number of neighbors of subdomain I
!!          XADJDOM(I) ..... Adjacency pointer for subdomain I
!!          ADJDOM ......... Adjacancy array of list of neighboring subdomains
!!                           List of neighbors of I: ADJDOM(XADJDOM(I):XADJDOM(I+1)-1)
!!          \endverbatim
!!
!!          A simple example to illustrate arrays:
!!
!!          \verbatim
!!
!!            +----+------------+
!!            |    |      2     |
!!            | 1  +------+-----+
!!            |    |      |     |
!!            +----+      |  3  |
!!            |           |     |
!!            |     4     |     |
!!            |           |     |
!!            +-----------+-----+
!!
!!            XADJDOM and ADJDOM: Subdomain adjancancy graph
!!
!!            I = 1 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 2,4
!!            I = 2 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 1,3,4
!!            I = 3 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 2,4
!!            I = 3 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 1,2,3
!!
!!          \endverbatim
!> @} 
!------------------------------------------------------------------------
subroutine par_domgra(inter,front,memor)
  use def_parame
  use def_kintyp
  use def_domain
  use def_master
  use def_parall
  use mod_memory
  use mod_parall, only : par_memor
  implicit none
  integer(8),  intent(out) :: memor(2)
  integer(ip), intent(out) :: inter
  integer(ip), intent(out) :: front
  integer(ip)              :: inode,dom2,ii,jj,kk
  integer(ip)              :: ndomi,domin,dom1
  integer(ip), pointer     :: domli(:) => null()
  !
  ! Create a half-matrix with the neighbourhood domains (including main diagonal)
  ! 
  call memory_alloca(par_memor,'NEIGHDOM','par_domgra' , neighDom , (npart_par*(npart_par+1))/2 )

  do dom1 = 1,npart_par
     do dom2 = 1,dom1
        neighDom( (dom1*(dom1-1))/2 + dom2) = 0
     end do
     npoin_par(dom1) = 0
  end do

  call memory_alloca(par_memor,'DOMLI','par_domgra' , domli , mepoi )

  front = 0
  inter = 0

  do inode = 1,npoin

     call par_domlis( pelpo, lelpo, inode, lepar_par, ndomi, domli )

     if( ndomi == 1 ) then
        !
        ! This is an internal node
        !
        inter            = inter + 1
        domin            = domli(1)
        lnpar_par(inode) = domin
        npoin_par(domin) = npoin_par(domin) + 1  
    
     else
        !
        ! This is a boundary node
        ! neighDom(i,j)
        !      if i.ne.j   boundary nodes shared by domains i and j
        !      if i.eq.j   boundary nodes of domain i
        !
        front = front + 1
        do ii = 1,ndomi
           dom1             = domli(ii)
           lnpar_par(inode) = 0
           npoin_par(dom1)  = npoin_par(dom1) + 1
           do jj = 1,ii
              dom2 = domli(jj)
              if( dom1 > dom2 ) then
                 kk = (dom1*(dom1-1))/2 + dom2
              else
                 kk = (dom2*(dom2-1))/2 + dom1
              end if
              neighDom(kk) = neighDom(kk) + 1          
           end do
        end do

     end if

  end do
  call memory_deallo(par_memor,'DOMLI','par_domgra' , domli )
  !
  ! The sum of the nodes that will have every domain
  !
  npoin_total = 0
  do domin = 1,npart_par
     npoin_total = npoin_total + npoin_par(domin)
  enddo
  nelem_total = nelem
  call par_disbou( memor )
  !
  ! Count the number of neighbours per domain 
  !
  do dom1 = 1, npart_par
     lneig_par(dom1) = 0
  enddo

  do dom1 = 1, npart_par
     kk = (dom1*(dom1-1))/2
     do dom2 = 1, dom1-1
        if( neighDom(kk + dom2) /= 0 ) then
           lneig_par(dom1) = lneig_par(dom1) + 1
           lneig_par(dom2) = lneig_par(dom2) + 1
        end if
     end do
  end do
  !
  ! Build the domain interconnection graph     
  !
  call par_memory(16_ip)

  xadjDom(1) = 1
  do dom1 = 1, npart_par
     xadjDom(dom1+1) = xadjDom(dom1) + lneig_par(dom1)
  end do

  call par_memory(17_ip)
  !
  ! iwa is a vector with indices to insert adjacencies
  !
  do dom1= 1, npart_par
     kk = (dom1*(dom1-1))/2
     do dom2= 1, dom1-1
        if( neighDom(kk + dom2) /= 0 ) then
           adjDom(xadjDom(dom1)) = dom2
           xadjDom(dom1)         = xadjDom(dom1) + 1
           adjDom(xadjDom(dom2)) = dom1
           xadjDom(dom2)         = xadjDom(dom2) + 1
        end if
     end do
  end do

  xadjDom(1) = 1
  do dom1 = 1,npart_par
     xadjDom(dom1+1) = xadjDom(dom1) + lneig_par(dom1)
  end do

end subroutine par_domgra
