!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_fields
  !------------------------------------------------------------------------
  !****f* chemic/chm_fields
  ! NAME
  !    chm_fields
  ! DESCRIPTION
  !    This routine:
  !    Reads the initial values for the cocentration
  !    It is important to keep the same order of the species in the fields
  !    as used in ALYA (check the order because it may differ from *chm.dat)
  !
  ! USES
  ! USED BY
  !    chm_iniunk
  !***
  !------------------------------------------------------------------------
  use def_master,      only : therm, intost
  use def_domain,      only : xfiel, nfiel, npoin
  use def_chemic,      only : mixedEq_eqs_chm, bvess_chm, kfl_model_chm, kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm,&
                              kfl_spec_name_chm, kfl_weigh_in_eq_CMC_chm, nclas_chm, nspec_chm, nvar_CMC_chm, kfl_field_chm,&
                              kfl_fixno_chm, Field_ind_chm
  use def_kintyp,      only : ip
  use mod_chm_mixedEq, only : CHM_INI_OLDWAY
  use mod_chm_mixedEq, only : CHM_INI_FIELD
  use mod_chm_mixedEq, only : CHM_INI_CONST
  implicit none
  integer(ip)             :: ipoin,iclas,ifield,nField

  external                :: runend

  !----------------------------------------------------------------------
  !
  ! Write concentration fields
  !
  !----------------------------------------------------------------------

   if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then

      if (kfl_weigh_in_eq_CMC_chm == 1) then  ! Start from alfa field
         do ipoin=1,npoin
            bvess_chm(nvar_CMC_chm+1,ipoin) = xfiel(kfl_field_chm(3))%a(1,ipoin,1)
         end do

         if (kfl_solve_enth_CMC_chm /= 0) then
            do ipoin = 1, npoin
               bvess_chm(nvar_CMC_chm,ipoin) = xfiel(kfl_field_chm(2))%a(1,ipoin,1)
            end do
         end if

      else  ! Providing the conditional values
         call runend('CHEMIC REAPHY: providing conditional values for each field in CMC not implemented yet')
      end if


   else
      if( kfl_spec_name_chm > 1_ip ) then
         do ipoin=1,npoin
            do iclas = 1,kfl_spec_name_chm
               if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
                  bvess_chm( (Field_ind_chm( iclas ) + 1_ip),ipoin) = xfiel(iclas+kfl_field_chm(1)-1_ip)%a(1,ipoin,1)
               end if
            enddo
         enddo
      else
         do ipoin=1,npoin
            do iclas = 1,nclas_chm

               select case(mixedEq_eqs_chm(iclas) % kfl_ini_type)
                  
               case (CHM_INI_OLDWAY)
                  !
                  ! Soot model
                  !
                  nField = nspec_chm
                  if (kfl_field_chm(4) > 0_ip) &
                     nField = nclas_chm
                  !
                  ! Initialize based on FIELDS: SPECI=n
                  !
                  if (iclas <= nField ) then
                     ifield = iclas+kfl_field_chm(1)-1_ip
                  else
                     ifield = 0
                  endif

               case (CHM_INI_FIELD)
                  !
                  ! Initialize based on individual field given to equation
                  !
                  ifield = mixedEq_eqs_chm(iclas) % kfl_ini_field

               case (CHM_INI_CONST)
                  !
                  ! Initialize to constant value
                  !
                  if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
                      bvess_chm(iclas,ipoin) = mixedEq_eqs_chm(iclas) % ini_value
                  end if
               end select

               !
               ! Field based approaches
               !
               select case(mixedEq_eqs_chm(iclas) % kfl_ini_type)
               case (CHM_INI_OLDWAY, CHM_INI_FIELD)
                  if (ifield > 0) then
                     if (ifield > nfiel) then
                        !
                        ! Throw human-readable error for non-existent fields
                        !
                        call runend('chm_fields: trying to access field '//trim(intost(ifield))//' that is not defined in .dom.dat')
                     else
                        !
                        ! Apply field
                        !
                        if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
                            bvess_chm(iclas,ipoin) = xfiel(ifield)%a(1,ipoin,1)
                        end if
                     endif
                  endif
               end select
            enddo
         enddo
      end if
      !
      ! Initialize therm field if input is present
      !
      if ( kfl_field_chm(2) > 0_ip .and. associated(therm) ) then
         ifield = kfl_field_chm(2)
         if (ifield > 0 .and. ifield <= nfiel) then
            do ipoin = 1,npoin
               therm(ipoin,1) =  xfiel(ifield)%a(1,ipoin,1)
            end do
         elseif ( ifield > nfiel ) then
             call runend('chm_fields: trying to access field '//trim(intost(ifield))//' that is not defined in .dom.dat')
         endif
      endif

   end if

end subroutine chm_fields

