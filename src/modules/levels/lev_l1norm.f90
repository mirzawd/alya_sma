!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_l1norm
  !-----------------------------------------------------------------------
  !****f* Levels/lev_l1norm
  ! NAME 
  !    lev_l1norm
  ! DESCRIPTION
  !    Compute the norm of the L1 error ||H(phi)-H(phi_ex)||_L1
  !                                    
  ! USES
  !    
  ! USED BY
  !    lev_
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_levels
  use      def_domain
  use      def_solver
  implicit none
  integer(ip)             :: ielem,idime     ! Indices and dimensions
  integer(ip)             :: pelty,pnode,pgaus
  integer(ip)             :: inode,ipoin
  integer(ip)             :: tetra(4_ip,6_ip)
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode),ellex(mnode)
  real(rp)                :: volc1, volc2

  l1nor_lev(1)=0_rp

  if( INOTMASTER ) then

     if(ndime==2_ip) then

        do ielem=1,nelem
           !
           ! Element dimensions
           !
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           pgaus=ngaus(pelty)

           !
           ! Gather operations
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode)=fleve(ipoin,1)
              ellex(inode)=flsex_lev(ipoin)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do

           volc1 = 0_ip
           volc2 = 0_ip

           if(pnode==3_ip) then
              call lev_elmv23(pnode,pgaus,pelty,elcod,ellev,volc1)
              call lev_elmv23(pnode,pgaus,pelty,elcod,ellex,volc2)
           else  if(pnode==4_ip) then
              call lev_elmv24(pnode,pgaus,pelty,elcod,ellev,volc1)
              call lev_elmv24(pnode,pgaus,pelty,elcod,ellex,volc2)
           endif   
           
           l1nor_lev(1) = l1nor_lev(1) + abs(volc1-volc2)

        enddo

     else if(ndime==3_ip) then

        if(mnode==8_ip) then

           tetra(1,1)=1
           tetra(2,1)=3
           tetra(3,1)=4
           tetra(4,1)=5

           tetra(1,2)=3
           tetra(2,2)=4
           tetra(3,2)=5
           tetra(4,2)=8

           tetra(1,3)=1
           tetra(2,3)=2
           tetra(3,3)=3
           tetra(4,3)=5

           tetra(1,4)=2
           tetra(2,4)=3
           tetra(3,4)=5
           tetra(4,4)=6

           tetra(1,5)=3
           tetra(2,5)=5
           tetra(3,5)=6
           tetra(4,5)=7

           tetra(1,6)=3
           tetra(2,6)=5
           tetra(3,6)=7
           tetra(4,6)=8
        endif

        do ielem=1,nelem

           !
           ! Element dimensions
           !

           pelty=ltype(ielem)
           pnode=nnode(pelty)
           pgaus=ngaus(pelty)

           !
           ! Gather operations
           !

           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              ellev(inode)=fleve(ipoin,1)
              ellex(inode)=flsex_lev(ipoin)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do

           volc1 = 0.0_rp
           if( pnode == 8_ip ) then
              call lev_elmv38(pnode,pgaus,pelty,elcod,ellev,tetra,volc1)
              call lev_elmv38(pnode,pgaus,pelty,elcod,ellex,tetra,volc2)
           else if( pnode == 4_ip ) then
              call lev_elmv34(pnode,elcod,ellev,volc1)
              call lev_elmv34(pnode,elcod,ellex,volc2)
           endif

           l1nor_lev(1) = l1nor_lev(1) + abs(volc1-volc2)

        enddo

     endif

  endif

  if( IPARALL ) then
     nparr =  1
     parre => l1nor_lev
     call par_operat(3_ip)
  end if

  if(kfl_paral<=0_ip) then
     print*,' l1nor ',l1nor_lev(1)
     write(lun_volum_lev,*) ' l1 norm ',l1nor_lev(1)
  endif

end subroutine lev_l1norm
