!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_upwmea(itask)
  !-----------------------------------------------------------------------
  ! Update
  ! NAME
  !    chm_upwmea
  ! DESCRIPTION
  !    Update each of the species specific heat and viscosity
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_master, only      :  ITASK_BEGSTE
  use def_master, only      :  ITASK_ENDSTE
  use def_master, only      :  ITASK_ENDITE
  use def_domain, only      :  npoin,ltype,nelem,ngaus
  use def_chemic, only      :  kfl_warni_chm,nspec_chm,kfl_tiaor_chm, &
                               kfl_tisch_chm,kfl_lookg_chm
  use def_master, only      :  conce,wmean,INOTMASTER,speci,wmean_gp
  implicit none
  integer(ip),intent(in)    :: itask   ! 1 = init step, 2=local update, 3=final update
  integer(ip)               :: ipoin,ispec,itime,ielem,pelty,&
                               igaus,pgaus

  external                  :: runend

  if (INOTMASTER) then
     select case(itask)

     case (ITASK_BEGSTE) ! At Begste

        if (associated(wmean)) then
           do ipoin = 1,npoin
              wmean(ipoin,2) = wmean(ipoin,1)
           enddo
        endif

     case (ITASK_ENDITE) ! At Endite(2)

        if (kfl_lookg_chm > 0 .and. associated(wmean_gp)) then
            do ielem = 1,nelem
               pelty = ltype(ielem)
               pgaus = ngaus(pelty)

               do igaus = 1,pgaus
                 wmean_gp(ielem) % a(igaus,2,1) = wmean_gp(ielem) % a(igaus,1,1)
               enddo
            end do
            if(kfl_tisch_chm==3) then
                call runend('Chemic upwmea: GP lookup is not implemented for Adams-Bashforth')
            endif
        else
           if (associated(wmean)) then
              do ipoin = 1,npoin
                 wmean(ipoin,2) = wmean(ipoin,1)
              enddo
           endif
        endif


     case (ITASK_ENDSTE) ! At Endste

        if (kfl_lookg_chm > 0 .and. associated(wmean_gp)) then
            if(kfl_tisch_chm/=3) then
                do ielem = 1,nelem
                   pelty = ltype(ielem)
                   pgaus = ngaus(pelty)

                   do igaus = 1,pgaus
                     wmean_gp(ielem) % a(igaus,3,1) = wmean_gp(ielem) % a(igaus,1,1)
                   enddo
                end do

            end if
        else
           if (associated(wmean)) then
              if(kfl_tisch_chm==2) then
                 !
                 ! BDF scheme
                 !
                 do ipoin=1,npoin
                    do itime=2+kfl_tiaor_chm,4,-1
                       wmean(ipoin,itime) = wmean(ipoin,itime-1)
                    end do
                 end do
              end if
              if(kfl_tisch_chm/=3) then
                 do ipoin = 1,npoin
                    wmean(ipoin,3) = wmean(ipoin,1)
                 enddo
              end if
           endif
        endif


     case (6_ip) ! Update mean molecular weight, usually at endite(1)

        if (associated(wmean)) then
           do ipoin = 1,npoin
              wmean(ipoin,1) = 0.0_rp
              do ispec = 1,nspec_chm
                 wmean(ipoin,1) = wmean(ipoin,1) +  conce(ipoin,ispec,1) / speci(ispec)%weigh
              enddo
              if (wmean(ipoin,1) /= 0.0_rp) then
                 wmean(ipoin,1) = 1.0_rp / wmean(ipoin,1)
              else
                 if (kfl_warni_chm==1) print *, ' WARNING: point with zero mean weight',ipoin
              endif
           enddo
        endif

     case (7_ip) ! For restart: wmean(:,3) = wmean(:,1)
        if (associated(wmean)) then
           do ipoin = 1,npoin
              wmean(ipoin,3) = wmean(ipoin,1)
           enddo
        endif

     end select
  endif

end subroutine chm_upwmea
