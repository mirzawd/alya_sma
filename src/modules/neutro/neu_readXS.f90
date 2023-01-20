!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NeutroInput
!> @{
!> @file    neu_readXC.f90
!> @author  A.S
!> @date    24/04/2019
!> @brief   Read XC data
!> @} 
!-----------------------------------------------------------------------
subroutine neu_readXS()

  use def_parame
  use def_inpout
  use def_master 
  use def_neutro
  use def_domain
  use mod_messages,          only : messages_live
  implicit none


  integer(ip) :: kk,jj,ii,hh,ll,G1,G2,cuenta,file_estado, nlineas, nlineas_tot, isour
  real(rp)    :: leo_leg(num_legendre_lee+1),graux,abs_aux,graux2,abs_aux2
  real(rp)    :: Factor_a_macroXS(num_materials_neu), group_energy_orig(num_energies_neu*2)
  character(10) :: word,word2              
  integer(ip) :: naux
  character (256) :: my_iomsg


   integer(ip) :: efectivos_neu_in_isotope(num_isotopes_max_neu,num_energies_neu), &
                  efectivos_neu_out_isotope(num_isotopes_max_neu,num_energies_neu)

   if( INOTSLAVE ) then

      if ( num_isotopes_max_neu>1_ip ) then
         ALLOCATE( grupo_energias_isotope(num_materials_neu, num_isotopes_max_neu, num_energies_neu), &
                    absor_neu_isotope(num_materials_neu, num_isotopes_max_neu, num_energies_neu), &
                    scatt_neu_isotope(num_materials_neu, num_isotopes_max_neu, num_energies_neu, &
                                       num_energies_neu, num_legendre_lee+1), &
                    kerma_isotope_neu(num_materials_neu, num_isotopes_max_neu, num_energies_neu))
         efectivos_neu_in_isotope = 0_ip
         efectivos_neu_out_isotope = 0_ip
         grupo_energias_isotope = 0.0_rp
         absor_neu_isotope = 0.0_rp
         scatt_neu_isotope = 0.0_rp
         kerma_isotope_neu = 0.0_rp
         
      endif
  
     
      do kk=1,num_materials_neu

        if(units_factor_neu==0) then
             Factor_a_macroXS(kk) = 1.0_rp
        elseif(units_factor_neu==1) then
             Factor_a_macroXS(kk) = Densidad_(kk)/(At_weight(kk)*uma_a_g)*barns_cm
        elseif(units_factor_neu==2) then
             Factor_a_macroXS(kk) = Densidad_(kk)/(At_weight(kk)*uma_a_kg)*barns_metros
        endif

     !
     ! Read Total XS
     !
         if(fil_totalXS_neu(kk)/='   ') then  
            open(unit=9876,file=fil_totalXS_neu(kk),action='READ',iostat=file_estado,iomsg=my_iomsg)
            if (file_estado/=0) then
               call runend('NEUTRO: '//trim(my_iomsg))
            end if
            if ( num_isotopes_neu(kk) == 1_ip ) then
               read(9876,*,IOSTAT=file_estado) word,word2

               ii=1
               do jj=1,num_energies_neu
                  read(9876,*,IOSTAT=file_estado) graux,abs_aux
                  read(9876,*,IOSTAT=file_estado) graux2,abs_aux2
                  group_energy_orig(jj*2-1) = graux
                  group_energy_orig(jj*2  ) = graux2
                  grupo_energias(kk,jj)=(graux+graux2)*0.5_rp
                  absor_neu(kk,jj) = (abs_aux+abs_aux2)*0.5_rp*Factor_a_macroXS(kk)
                  ! ener_weigd_neu(jj) = 1.0 ! /real(num_energies_neu) 
               enddo
            else
               do ii=1,num_isotopes_neu(kk)
                  read(9876,*,IOSTAT=file_estado) word,word2
                  do jj=1,num_energies_neu
                     read(9876,*,IOSTAT=file_estado) graux,abs_aux
                     read(9876,*,IOSTAT=file_estado) graux2,abs_aux2
                     group_energy_orig(jj*2-1) = graux
                     group_energy_orig(jj*2  ) = graux2
                     grupo_energias_isotope(kk,ii,jj) = (graux+graux2)*0.5_rp
                     absor_neu_isotope(kk,ii,jj) = (abs_aux+abs_aux2)*0.5_rp
                  enddo
               enddo
               do jj=1,num_energies_neu
                  do ii=1,num_isotopes_neu(kk)
                     grupo_energias(kk,jj) = grupo_energias(kk,jj) + grupo_energias_isotope(kk,ii,jj)*&
                                                                     atom_percentage_isotope_neu(kk,ii)
                     absor_neu(kk,jj) = absor_neu(kk,jj) + absor_neu_isotope(kk, ii, jj)*atom_percentage_isotope_neu(kk,ii)
                  enddo
                  absor_neu(kk,jj) = absor_neu(kk,jj)*Factor_a_macroXS(kk)
               enddo

            end if
            close(9876)
         else
            call runend('(NEUTRO) NO TOTAL FILE DEFINED FOR MATERIAL '//trim(intost(kk)))
         endif
     !
     ! Read scattering XS
     !
         ! ii = 0
       if(fil_scattXS_neu(kk)/='   ') then  
          open(unit=9876,file=fil_scattXS_neu(kk),action='READ',iostat=file_estado,iomsg=my_iomsg)
          if (file_estado/=0) then
            call runend('NEUTRO: '//trim(my_iomsg))
          end if
          if ( num_isotopes_neu(kk) == 1_ip ) then
            cuenta=num_energies_neu+1_ip
            read(9876,*,IOSTAT=file_estado) word,word,word
            read(9876,*,IOSTAT=file_estado) G1,G2,(leo_leg(ii),ii=1,num_legendre_lee+1)
            do while(file_estado>=0_ip)
               do ii=1,num_legendre_lee+1_ip
                  scatt_neu(kk,G2,G1,ii)=leo_leg(ii)*Factor_a_macroXS(kk)
               enddo
               if(cuenta/=G1) then
                  ! cuenta=cuenta-1_ip
                  cuenta=G1
                  efectivos_neu_out(cuenta)=G1
                  efectivos_neu_in(cuenta)=G2
                  n_efectivos_neu_sparse(kk, cuenta) = 1_ip
               else
                  efectivos_neu_in(cuenta)=G2
                  n_efectivos_neu_sparse(kk, cuenta) = n_efectivos_neu_sparse(kk, cuenta) + 1_ip
               endif
               read(9876,*,IOSTAT=file_estado) G1,G2,(leo_leg(ii),ii=1,num_legendre_lee+1)
            enddo
         else
            do jj=1,num_isotopes_neu(kk)
               nlineas = 0_ip
               naux = 0_ip
               cuenta=num_energies_neu+1_ip
               read(9876,*,IOSTAT=file_estado) word,word,word,nlineas_tot
               do while(file_estado>=0 .and. nlineas < nlineas_tot)
                  nlineas = nlineas + 1_ip
                  read(9876,*,IOSTAT=file_estado) G1,G2,(leo_leg(ii),ii=1,num_legendre_lee+1)
                  do ii=1,num_legendre_lee+1_ip
                     scatt_neu_isotope(kk,jj,G2,G1,ii)=leo_leg(ii)
                  enddo
                  if(cuenta/=G1) then
                     ! cuenta=cuenta-1_ip
                     cuenta=G1
                     efectivos_neu_out_isotope(jj, cuenta)=G1
                     efectivos_neu_in_isotope(jj, cuenta)=G2
                  else
                     efectivos_neu_in_isotope(jj, cuenta)=G2
                  endif
               enddo
            enddo
            do ii=1,num_energies_neu
               n_efectivos_neu_sparse(kk, ii) = 0_ip
               do hh=1,num_energies_neu
                  do ll=1,num_legendre_lee+1
                     do jj=1,num_isotopes_neu(kk)
                        scatt_neu(kk,ii,hh,ll) = scatt_neu(kk,ii,hh,ll) + &
                                                   scatt_neu_isotope(kk,jj,ii,hh,ll)*atom_percentage_isotope_neu(kk,jj)
                     enddo
                     scatt_neu(kk,ii,hh,ll) = scatt_neu(kk,ii,hh,ll)*Factor_a_macroXS(kk)
                  enddo
                  if (scatt_neu(kk, ii, hh, 1)/=0.0_rp) n_efectivos_neu_sparse(kk, hh) = n_efectivos_neu_sparse(kk, hh) + 1_ip
               enddo
               do jj=2,num_isotopes_neu(kk)
                  if ( efectivos_neu_out_isotope(1,ii) /= efectivos_neu_out_isotope(jj,ii) ) &
                     call runend('neu_readXS: error in final group between components of effective material')
               enddo
               efectivos_neu_out(ii) = efectivos_neu_out_isotope(1,ii)
               efectivos_neu_in(ii) = maxval(efectivos_neu_in_isotope(:,ii))
            enddo
         endif
       
         close(9876)
      else
         call runend('(NEUTRO) NO ELAST FILE DEFINED FOR MATERIAL '//trim(intost(kk)))
      endif   

      ! Guardo valores acumulados de no nulos (grupos conectados) por grupo 
      do ii=1,num_energies_neu
         n_efectivos_neu_sparse_acum(kk, ii+1) = n_efectivos_neu_sparse_acum(kk, ii)+n_efectivos_neu_sparse(kk, ii)
      enddo

        !
     ! Read fission XS
     !
      if(fil_fisionXS_neu(kk)/='   ') then  
         open(unit=9876,file=fil_fisionXS_neu(kk),action='READ',iostat=file_estado,iomsg=my_iomsg)
         if (file_estado/=0) then
            call runend('NEUTRO: '//trim(my_iomsg))
         end if
          read(9876,*,IOSTAT=file_estado) word,word,word
          do while(file_estado>=0)
              read(9876,*,IOSTAT=file_estado) G1,G2,fiss_neu(kk,G1,G2)
            fiss_neu(kk,G1,G2)=fiss_neu(kk,G1,G2)*Factor_a_macroXS(kk)
         enddo
       
         close(9876)
      else
         call messages_live('(NEUTRO) NO FISSION FILE DEFINED FOR MATERIAL '//trim(intost(kk)),'WARNING')
      endif
       
      ! Read source
!       if(fil_source_neu(kk)/='   ') then  
!          open(unit=9876,file=fil_source_neu(kk))
!          !read(9876,*) word,word2
!          read(9876,*,IOSTAT=file_estado)  word,word2
!  !        do while(.not. eof(9876))
!  !           read(9876,*) G1,source_bound_neu(G1)
!           do while(file_estado>=0)
!               read(9876,*,IOSTAT=file_estado)G1,source_bound_neu(G1)
!          enddo
       
!          close(9876)

!       else
!          call messages_live('(NEUTRO) NO SOURCE FILE FOR MATERIAL '//trim(intost(kk)),'WARNING')
!       endif

     !
     ! Read KERMA factors 
     ! (initially in [eV*barns], converted to [J*m^2] and *m^-3 = [J*m^-1])
     !
      if(fil_kerma_neu(kk)/='   ') then  
         open(unit=9876,file=fil_kerma_neu(kk),action='READ',iostat=file_estado,iomsg=my_iomsg)
         if (file_estado/=0) then
            call runend('NEUTRO: '//trim(my_iomsg))
         end if
         if ( num_isotopes_neu(kk) == 1 ) then
            read(9876,*,IOSTAT=file_estado) word,word2

            ii=1
            do jj=1,num_energies_neu
               read(9876,*,IOSTAT=file_estado) G1, kerma_neu(kk, jj)
               kerma_neu(kk, jj) = kerma_neu(kk, jj)*Factor_a_macroXS(kk)*ev_a_joule
            enddo
         else
            do ii=1,num_isotopes_neu(kk)
               read(9876,*,IOSTAT=file_estado) word,word2
               do jj=1,num_energies_neu
                  read(9876,*,IOSTAT=file_estado) G1, kerma_isotope_neu(kk, ii, jj)
               enddo
            enddo
            do jj=1,num_energies_neu
               do ii=1,num_isotopes_neu(kk)
                  kerma_neu(kk, jj) = kerma_neu(kk, jj) + kerma_isotope_neu(kk, ii, jj)*atom_percentage_isotope_neu(kk,ii)
               enddo
               kerma_neu(kk, jj) = kerma_neu(kk, jj)*Factor_a_macroXS(kk)*ev_a_joule
            enddo

         end if
         close(9876)
      else
         call messages_live('(NEUTRO) NO KERMA FILE DEFINED FOR MATERIAL '//trim(intost(kk)),'WARNING')
      endif

    enddo

    !READ SOURCES
      do isour=1,num_sources_neu
         ! if(fil_source_neu(kk)/='   ') then  
            open(unit=9876,file='SOURF'//trim(intost(isour))//'.dat',action='READ',iostat=file_estado,iomsg=my_iomsg)
            if (file_estado/=0) then
               call runend('NEUTRO: '//trim(my_iomsg))
            end if
            read(9876,*,IOSTAT=file_estado)  word,word2
             do while(file_estado>=0)
                 read(9876,*,IOSTAT=file_estado)G1,source_bound_neu(isour,G1)
            enddo
          
            close(9876)
         ! endif
      enddo

    ! Armo el vector sparse con los no-nulos (grupos finales) por grupo y material
    max_val_acum_neu = maxval(n_efectivos_neu_sparse_acum)
    allocate(efectivos_neu_sparse(num_materials_neu, max_val_acum_neu))
    efectivos_neu_sparse=0
    do kk=1,num_materials_neu
      naux = 0_ip
      do ii=1,num_energies_neu
         do jj=1,num_energies_neu
            if (scatt_neu(kk, jj, ii, 1) /= 0.0_rp) then
               naux = naux+1_ip
               efectivos_neu_sparse(kk,naux) = jj
            end if
         enddo
      enddo
    enddo

    ! Escribo valores finales de factor a macro, total XS y scatter xs por material en NEU.LOG
    do kk=1,num_materials_neu
      write(momod(modul)%lun_outpu,*) ''
      write(momod(modul)%lun_outpu,*) '       Material: ',kk,' Factor macro XS= ',Factor_a_macroXS(kk)
      write(momod(modul)%lun_outpu,*) ''

      do isour=1,num_sources_neu
         do hh=1,num_energies_neu
            write(momod(modul)%lun_outpu,'(2(A,I3),A,ES15.6)') &
               '        Source ',isour, ' energy group: ', hh, ' = ', source_bound_neu(isour,hh)
         enddo
      enddo

      write(momod(modul)%lun_outpu,*) ''
      do hh=num_energies_neu,1,-1
         write(momod(modul)%lun_outpu,'(A,I3,2ES15.6,A,ES15.6)') &
            '        Energy group: ', hh, group_energy_orig(2*hh-1),group_energy_orig(2*hh), &
            ' Total XS= ', absor_neu(kk,hh)
      enddo

      write(momod(modul)%lun_outpu,*) ''
      write(momod(modul)%lun_outpu,*) '        1. In group  2. Out group  3. Leg coeff*Factor macro'
      do hh=num_energies_neu,1,-1
         do jj=n_efectivos_neu_sparse_acum(kk, hh+1)-1, n_efectivos_neu_sparse_acum(kk,hh), -1
            ii = efectivos_neu_sparse(kk, jj)
            write(momod(modul)%lun_outpu,'(2(A,I3), '//intost(num_legendre_lee+1)//'ES15.6)') &
               '        ',hh,'   ',ii, (scatt_neu(kk,ii,hh,ll),ll=1,num_legendre_lee+1)
         enddo
      end do
   end do

    if ( num_isotopes_max_neu > 1_ip ) DEALLOCATE( grupo_energias_isotope, absor_neu_isotope, scatt_neu_isotope, &
                                                atom_percentage_isotope_neu, kerma_isotope_neu)

endif


end subroutine neu_readXS
