!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_outwit()
  !------------------------------------------------------------------------
  !****f* Turbul/tur_outwit
  ! NAME 
  !    tur_outwit
  ! DESCRIPTION
  !    Output results on witness points
  ! USES
  ! USED BY
  !    tur_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_turbul
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  
  implicit none
  integer(ip) :: iwitn,ielem,inode,pnode,pelty,ipoin,ivawi
  real(rp)    :: turb1, turb2

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then

     if( INOTMASTER ) then 

        do iwitn = 1, nwitn
           ielem = lewit(iwitn)
           if( ielem > 0 ) then

              pelty = ltype(ielem)
              pnode = nnode(pelty)
              do ivawi = 1,postp(1) % nvawi
                 witne(ivawi,iwitn) = 0.0_rp
              end do

              if( postp(1) % npp_witne(1) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    if ( kfl_regularization) then 
                       turb1 = regul_k(untur(1,ipoin,1))
                    else
                       turb1 = untur(1,ipoin,1)
                    end if 
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) * turb1 * postp(1) % witne_dt(1)
                 end do
              end if

              if( postp(1) % npp_witne(2) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    if ( kfl_regularization) then 
                       turb2 = regul_e(untur(2,ipoin,1))
                    else
                       turb2 = untur(2,ipoin,1)
                    end if
                    witne(2,iwitn) = witne(2,iwitn) + shwit(inode,iwitn) * turb2 * postp(1) % witne_dt(2)
                 end do
              end if

           end if
        end do

     end if

  end if

end subroutine tur_outwit

