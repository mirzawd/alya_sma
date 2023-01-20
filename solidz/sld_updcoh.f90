!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_updcoh(itask)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_updcoh
  ! NAME
  !    sld_updcoh
  ! DESCRIPTION
  !    This routine performs several updates for cohesive law parameters
  ! USED BY
  !    sld_begite (itask=1)
  !    sld_endite (itask=2)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solidz
  use def_elmtyp
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iboun,ielem,igaub,pgaub,pblty,kboun
  logical(lg)             :: debugging
  
  debugging = .false.

  if( kfl_xfeme_sld > 0 .and. INOTMASTER ) then
  

     select case (itask)

     case(0)
          !
          ! Initialize internal variables of cohesive laws and friction
          !
          dceff_sld(:,:,:,1) = 0.0_rp
          dslip_sld(:,:,:,1) = 0.0_rp
          dfric_sld(:,:,:,1) = 0.0_rp

     case(1)
          !
          ! Wind forward the effective opening and friction/contact
          !
          dceff_sld(:,:,:,2) = dceff_sld(:,:,:,1)
          dslip_sld(:,:,:,2) = dslip_sld(:,:,:,1)
          dfric_sld(:,:,:,2) = dfric_sld(:,:,:,1)
                  
     case(2)
          !
          ! Update the maximum opening
          !
          do ielem = 1,nelem

            if (lelch(ielem)==ELCUT) then
               if( ndime == 2 ) then
                 pblty = BAR02
                 kboun = cutel(ielem) % nboun
               else
                 pblty = TRI03
                 kboun = cutel(ielem) % nboun
               end if
               pgaub = ngaus(pblty)
            else
               if( ndime == 2 ) then
                 pblty = BAR02
                 kboun = 1
               else
                 pblty = TRI03
                 kboun = 1
               end if
               pgaub = ngaus(pblty)
            end if               
               
            do iboun = 1,kboun
              do igaub = 1,pgaub
                dcmax_sld(ielem,iboun,igaub) = max(dcmax_sld(ielem,iboun,igaub), &
                                                   dceff_sld(ielem,iboun,igaub,1))
              end do
            end do
          end do

          if (debugging) write(*,*)' '
          if (debugging) write(*,*)'dceff_old = ',dceff_sld(22,1,2,2)
          if (debugging) write(*,*)'dceff_new = ',dceff_sld(22,1,2,1)
          if (debugging) write(*,*)'dcmax_sld = ',dcmax_sld(22,1,2)
          if (debugging) write(*,*)'dfric_old = ',dfric_sld(22,1,2,2)
          if (debugging) write(*,*)'dfric_new = ',dfric_sld(22,1,2,1)

      end select

  end if

end subroutine sld_updcoh

