!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outwit()
  !------------------------------------------------------------------------
  !****f* Temper/tem_outwit
  ! NAME 
  !    tem_outwit
  ! DESCRIPTION
  !    Output results on witness points
  ! USES
  ! USED BY
  !    tem_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_temper
  use mod_memory, only : memory_alloca, memory_deallo
  implicit none
  integer(ip)          :: iwitn,ielem,inode,pnode,pelty,ipoin

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
     !
     ! Results on witness points
     !
     witne => postp(1) % witne

     if( INOTMASTER ) then

        do iwitn = 1, nwitn
           ielem = lewit(iwitn)
           if( ielem > 0 ) then
              pelty = ltype(ielem)
              pnode = nnode(pelty)

              if( postp(1) % npp_witne(1) == 1 ) then
                 !
                 ! T
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) * tempe(ipoin,1) * postp(1) % witne_dt(1)
                 end do
              end if

              if( postp(1) % npp_witne(2) == 1 ) then
                 !
                 ! dT/dx
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(2,iwitn) = witne(2,iwitn) + dewit(1,inode,iwitn) * tempe(ipoin,1) * postp(1) % witne_dt(2)
                 end do
              end if

              if( postp(1) % npp_witne(3) == 1 ) then
                 !
                 ! dT/dy
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(3,iwitn) = witne(3,iwitn) + dewit(2,inode,iwitn) * tempe(ipoin,1) * postp(1) % witne_dt(3)
                 end do
              end if
 
              if( postp(1) % npp_witne(4) == 1 ) then
                 !
                 ! dT/dz
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(4,iwitn) = witne(4,iwitn) + dewit(3,inode,iwitn) * tempe(ipoin,1) * postp(1) % witne_dt(4)
                 end do
              end if

              if( postp(1) % npp_witne(5) == 1 ) then
                 !
                 ! T or h
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(5,iwitn) = witne(5,iwitn) + shwit(inode,iwitn) * therm(ipoin,1) * postp(1) % witne_dt(5)
                 end do
              end if

           end if
        end do
        
     end if

  end if

end subroutine tem_outwit

