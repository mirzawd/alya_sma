!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_calvol
  !-----------------------------------------------------------------------
  !****f* Levels/lev_calvol
  ! NAME 
  !    lev_calvol
  ! DESCRIPTION
  !    Compute the volume of the phase with phi positive
  !                                    
  ! USES
  !    lev_elmv23
  !    lev_elmv24
  !    lev_elmv34
  !    lev_elmv38  
  !    lev_elmv36  
  ! USED BY
  !    lev_updunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip)             :: ielem, idime
  integer(ip)             :: pelty, pnode, pgaus
  integer(ip)             :: inode, ipoin
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode)
  integer(ip)             :: tetra(4_ip,6_ip), tetrp(4_ip,3_ip), tetpy(4_ip,2_ip)
  integer(ip)             :: ipasq, ipasp, ipapy
  real(rp)                :: voloc  
  !
  ! Loop over the elements to find cuts (2d case)!!!
  !
  volit_lev = 0.0_rp

  if( INOTMASTER ) then

     if( ndime == 2 ) then

        do ielem = 1,nelem

           pelty = ltype(ielem)

           if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              do inode = 1,pnode
                 ipoin                = lnods(inode,ielem)
                 ellev(inode)         = fleve(ipoin,1)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              end do   
              if( pnode == 3_ip ) then
                 call lev_elmv23(pnode,pgaus,pelty,elcod,ellev,voloc)
              else if( pnode == 4_ip ) then
                 call lev_elmv24(pnode,pgaus,pelty,elcod,ellev,voloc)
              else
                 call runend('LEV_CALVOL: ELEMENT NOT CODED')
              end if

              volit_lev = volit_lev + voloc

           end if

        end do

     else if( ndime == 3 ) then

        ipasq = 0_ip
        ipasp = 0_ip
        ipapy = 0_ip

        do ielem=1,nelem

           !
           ! Element dimensions
           !

           pelty=ltype(ielem)
           if( pelty > 0 ) then
              pnode=nnode(pelty)
              pgaus=ngaus(pelty)

              if((pnode==8_ip).and.(ipasq==0_ip)) then

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

                 ipasq = 1

              else if((pnode==6_ip).and.(ipasp==0_ip)) then

                 ipasp = 1

                 tetrp(1,1)=1
                 tetrp(2,1)=2
                 tetrp(3,1)=3
                 tetrp(4,1)=4

                 tetrp(1,2)=2
                 tetrp(2,2)=3
                 tetrp(3,2)=4
                 tetrp(4,2)=5

                 tetrp(1,3)=3
                 tetrp(2,3)=4
                 tetrp(3,3)=5
                 tetrp(4,3)=6

              else if((pnode==5_ip).and.(ipapy==0_ip)) then

                 ipapy = 1

                 tetpy(1,1)=1
                 tetpy(2,1)=2
                 tetpy(3,1)=3
                 tetpy(4,1)=5

                 tetpy(1,2)=1
                 tetpy(2,2)=3
                 tetpy(3,2)=4
                 tetpy(4,2)=5

              endif

              !
              ! Gather operations
              !

              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 ellev(inode)=fleve(ipoin,1)
                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do

              voloc = 0.0_rp
              ! Case Q1
              if( pnode == 8_ip ) then
                 call lev_elmv38(pnode,pgaus,pelty,elcod,ellev,tetra,voloc)
                 ! Case P1
              else if( pnode == 4_ip ) then
                 call lev_elmv34(pnode,elcod,ellev,voloc)
                 ! Case Prism
              else if( pnode == 6_ip ) then
                 call lev_elmv36(pnode,elcod,ellev,voloc,tetrp)
                 ! Case Pyramid
              else if( pnode == 5_ip) then   
                 call lev_elmv35(pnode,elcod,ellev,voloc,tetpy)
              endif
              volit_lev = volit_lev + voloc
           end if

        enddo

     endif
  endif

  call PAR_SUM(volit_lev,'IN MY CODE')

end subroutine lev_calvol
