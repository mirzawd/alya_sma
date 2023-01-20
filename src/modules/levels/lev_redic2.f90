!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_redic2()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_redisc
  ! NAME 
  !    lev_redisc
  ! DESCRIPTION
  !    Compute the level set function redistanciation geometrically only on cut nodes when tyred_lev==2
  !    Geometrically in cut nodes but searching surfaces on the element connected to the node.  
  ! USES
  !    
  ! USED BY
  !    lev_redieq
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use def_solver
  use mod_memchk

  implicit none
  integer(ip)             :: ielem,idime
  integer(ip)             :: inode,ipoin,compt
  integer(ip)             :: pelty,pnode
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: ellev(mnode)


  !----------------------------------------------------------------------
  !
  ! Initializes the distance
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     do ipoin=1,npoin
        dista_lev(ipoin)=1.0e+06
!!$        print*,ipoin,' ',fleve(ipoin,3)
     end do
     nelcr_lev=0_ip
  end if
!!$  stop

  !----------------------------------------------------------------------
  !
  ! Put a flag on points in element crossed by the interface
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then

     do ipoin = 1,npoin
        icupt_lev(ipoin) = 0_ip
     enddo

     do ielem=1,nelem
        !
        ! Element dimensions
        !
        pelty=ltype(ielem)
        pnode=nnode(pelty)
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

        compt=0_ip
        !
        ! Determine type of element
        ! type= 1 empty element (phi<0) in each node
        ! type= 2 full element (phi>0) in each node
        ! type= 3 partly full element (phi<0) in some nodes and (phi>0) in some nodes
        !
        do inode=1,pnode
           if(ellev(inode)>=0.0_rp) then
              compt=compt+1
           endif
        end do

        if((compt/=0_ip).and.(compt/=pnode)) then

           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              icupt_lev(ipoin)=1_ip
           end do
           nelcr_lev = nelcr_lev + 1_ip
           elcro_lev(nelcr_lev) = ielem

        endif

     enddo
     !
     ! Compute level set function geometrically on cut nodes
     !
     call lev_caldis()

  endif

  if( IPARALL ) then
     call vocabu(NPOIN_REAL_1DIM,1_ip,0_ip)
     parr1 => dista_lev
     call par_minxch() 
  endif

  if(INOTMASTER) then

     do ipoin = 1,npoin
        if( dista_lev(ipoin) < 1.0e+06 ) then
           icupt_lev(ipoin) = 1_ip
           if( fleve(ipoin,1) >= 0.0_rp ) then
              fleve(ipoin,1) =  dista_lev(ipoin)
           else
              fleve(ipoin,1) = -dista_lev(ipoin)
           endif
        endif
     enddo

  endif


end subroutine lev_redic2
