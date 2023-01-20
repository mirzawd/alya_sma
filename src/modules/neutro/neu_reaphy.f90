!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NeutroInput
!> @{
!> @file    neu_reaphy.f90
!> @author  Guillaume Houzeaux
!> @date    29/03/2016
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine neu_reaphy()
  use def_parame
  use def_inpout
  use def_master 
  use def_neutro
  use def_domain
  use mod_ecoute, only :  ecoute
  use mod_ADR, only : ADR_read_data
  implicit none
  
  real(rp) :: source_bound
  integer(ip) :: group_bound
  integer(ip) :: nmat = 1, nisotope=1
!   integer :: ii, jj
  character(1):: sm

!   real(rp),allocatable :: At_weight_isotope_aux(:,:), mass_percentage_isotope_neu_aux(:,:), &
!                           atom_percentage_isotope_neu_aux(:,:) 
  real(rp) :: sum_mass, sum_atom, sum_aux
  logical :: mass_percentage, atomic_percentage

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     ! MODULEX
     !
     num_energies_neu   = 1_ip      ! default Number of energy groups 
     num_directions_neu = 10_ip      ! default Number of directions
     num_materials_neu   = 1_ip      !default materials
     num_sources_neu    = 1_ip      ! default number of sources
     kfl_icosa_neu      = 0_ip
     kfl_snord_neu      = 0_ip
     mass_percentage = .FALSE.
     atomic_percentage = .FALSE.
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('neu_reaphy') 
     do while( words(1) /= 'PHYSI' )
        call ecoute('neu_reaphy')
     end do
     call ecoute('neu_reaphy')
     !--><group>
     !-->    <groupName>PHYSICAL_PROBLEM</groupName>
     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Physical properties definition
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> PHYSICAL_PROBLEM
     !
     
        !
        ! Read ADR data
        !
        call ADR_read_data(ADR_neu) 
     do while( words(1) /= 'ENDPH' )

        if(      words(1) == 'ENERG' ) then 
           !
           ! Energy groups
           !
           num_energies_neu = getint('ENERG',1_ip,'#Number of energy groups')

        elseif( words(1) == 'DIREC' ) then 
           !
           ! Direction 
           !
           num_directions_neu = getint('DIREC',1_ip,'#Number of directions')

                 
        elseif( words(1) == 'NUNMA' ) then 
           
           num_materials_neu=getint('NUNMA',1_ip,'#Number of materials')
           
           allocate(Isotope(num_materials_neu))
           allocate(At_weight(num_materials_neu))
           allocate(Densidad_(num_materials_neu))
           allocate(aniso_neu(num_materials_neu))
           allocate(absor_neu_cte(num_materials_neu))
           allocate(scatt_neu_cte(num_materials_neu))
           allocate(funsource(num_materials_neu))
           
           funsource=0.0_rp
          
           allocate(fil_totalXS_neu(num_materials_neu))
           allocate(fil_scattXS_neu(num_materials_neu))
           allocate(fil_fisionXS_neu(num_materials_neu))
           allocate(fil_source_neu(num_materials_neu))
           allocate(fil_kerma_neu(num_materials_neu))

           allocate(num_isotopes_neu(num_materials_neu))
           num_isotopes_neu=1_ip

           fil_totalXS_neu = '   '
           fil_scattXS_neu = '   '
           fil_fisionXS_neu = '   '
           fil_source_neu = '   '
           fil_kerma_neu = '   '
        
         else if( words(1)=='METHO')   then       
           !             
           ! Numerical method
           !
           if( words(2) =='DOM  ') then
              if(      words(3) == 'ICO20')  then 
                 num_directions_neu = 20
                 kfl_icosa_neu = 1  
              else if( words(3) == 'ICO80')  then 
                 num_directions_neu = 80         
                 if( ndime == 2 ) num_directions_neu = 40 
                 kfl_icosa_neu = 1
              else if( words(3) == 'SN8  ')   then 
                 num_directions_neu = 80
                 if( ndime == 2 ) num_directions_neu = 40 
                 kfl_snord_neu = 8
              else if( words(3) == 'SN6  ')   then 
                 num_directions_neu = 48
                 if( ndime == 2 ) num_directions_neu = 24        
                 kfl_snord_neu = 6
              else if( words(3) == 'SN4  ')   then 
                 num_directions_neu = 24
                 if( ndime == 2 ) num_directions_neu = 12        
                 kfl_snord_neu = 4
              else if( words(3) == 'SN10 ')   then 
                 num_directions_neu = 120
                 if( ndime == 2 ) num_directions_neu = 60 
                 kfl_snord_neu = 10
              else if( words(3) == 'SN12 ')   then 
                 num_directions_neu = 168
                 if( ndime == 2 ) num_directions_neu = 84 
                 kfl_snord_neu = 12
              else if( words(3) == 'SN16 ')   then 
                 num_directions_neu = 288
                 if( ndime == 2 ) num_directions_neu = 144 
                 kfl_snord_neu = 16
              else if( words(3) == 'SN64 ')   then 
                 num_directions_neu = 512
                 kfl_snord_neu = 64
              else 
                 call runend('neu_reaphy: error reading number of directions')                  
              end if
              
           end if
        else if( words(1)=='SOURC')   then   
            if( words(2) =='USER ') then
                 funsource=-1.0_rp
            elseif(words(2) =='MATER') then
                 funsource(1:num_materials_neu)=getrea('MATER',0.0_rp,'#sources')

            endif
         else if( words(1)=='NSOUR')   then     
          
            num_sources_neu=getint('NSOUR',1_ip,'#number of sources')  
        
         else if( words(1)=='OUTPU')   then     
          
            output_level_neu=getint('OUTPU',0_ip,'#Verbosity for NEUTRO')  

        else if( words(1)=='BOUNS')   then     
          
           source_bound=getrea('BOUNS',0.0_rp,'#source in bound')  

        else if( words(1)=='GRBOU')   then     
          
           group_bound = getint('GRBOU',1_ip,'#group in boundary')

        else if( words(1)=='ALBED')   then     
          
           albedo_neu=getrea('ALBED',1.0_rp,'#albedo in bound')  

        else if( words(1)=='NUNIT')   then       
             
           units_factor_neu=getint('NUNIT',0_ip,'#UNITS to use')  
        
        else if( words(1)=='NUMGR')   then       
             
             num_energies_neu = getint('NUMGR',1_ip,'#Number of energy groups')
        
        else if( words(1)=='REALE')   then      
        
             num_legendre_lee=getint('REALE',1_ip,'#Number of legendre leeido')

        else if( words(1)=='USELE')   then      
        
             num_legendre_neu=getint('USELE',1_ip,'#Number of legendre Usado')
        
         else if( words(1) == 'PROPE' ) then                    
            
            
            do while(words(1) /= 'ENDPR')
               
               if( words(1) == 'MATER' ) then
                  nmat = getint('MATER',1_ip,'#Number of mat')
               elseif( words(1)=='NUMIS')   then

                  num_isotopes_neu(nmat)=getint('NUMIS',1_ip,'#Number of isotopes of effective material')
      
                  if ( num_isotopes_neu(nmat)>num_isotopes_max_neu ) then
                     
                     if ( num_isotopes_neu(nmat) > 1_ip ) then
                        if ( associated(At_weight_isotope) ) then
                           allocate(At_weight_isotope_aux(num_materials_neu,num_isotopes_max_neu), &
                                    mass_percentage_isotope_neu_aux(num_materials_neu,num_isotopes_max_neu), &
                                    atom_percentage_isotope_neu_aux(num_materials_neu,num_isotopes_max_neu))
                           At_weight_isotope_aux = At_weight_isotope
                           mass_percentage_isotope_neu_aux = mass_percentage_isotope_neu
                           atom_percentage_isotope_neu_aux = atom_percentage_isotope_neu
                           deallocate(At_weight_isotope, mass_percentage_isotope_neu, atom_percentage_isotope_neu)
                        endif
                        allocate(At_weight_isotope(num_materials_neu,num_isotopes_neu(nmat)), &
                                 mass_percentage_isotope_neu(num_materials_neu,num_isotopes_neu(nmat)), &
                                 atom_percentage_isotope_neu(num_materials_neu,num_isotopes_neu(nmat)))
                        if ( allocated(At_weight_isotope_aux) ) then
                           At_weight_isotope(:, :num_isotopes_max_neu) = At_weight_isotope_aux(:, :num_isotopes_max_neu)
                           mass_percentage_isotope_neu(:, :num_isotopes_max_neu) = &
                              mass_percentage_isotope_neu_aux(:, :num_isotopes_max_neu)
                           atom_percentage_isotope_neu(:, :num_isotopes_max_neu) = &
                              atom_percentage_isotope_neu_aux(:, :num_isotopes_max_neu)
                           deallocate(At_weight_isotope_aux, mass_percentage_isotope_neu_aux, atom_percentage_isotope_neu_aux)
                        else
                           At_weight_isotope = 0.0_rp
                           mass_percentage_isotope_neu = 0.0_rp
                           atom_percentage_isotope_neu = 0.0_rp
                        end if

                     else
                        allocate(mass_percentage_isotope_neu(num_materials_neu,num_isotopes_neu(nmat)), &
                                 atom_percentage_isotope_neu(num_materials_neu,num_isotopes_neu(nmat)))
                        mass_percentage_isotope_neu = 0.0_rp
                        atom_percentage_isotope_neu = 0.0_rp
                     end if

                     num_isotopes_max_neu=num_isotopes_neu(nmat)

                  endif

               else if( words(1) == 'ISOTO' ) then
                  nisotope = getint('ISOTO',1_ip,'#Number of component')

               elseif( words(1) == 'ANISO' ) then
                  aniso_neu(nmat) =  getrea('ANISO',0.0_rp,'#anisotropy')
                  if( aniso_neu(nmat) < -1.0_rp ) call runend('neu_reaphy: anisotropy factor cannot be lower than -1.0')     
                  if( aniso_neu(nmat) >  1.0_rp ) call runend('neu_reaphy: anisotropy factor cannot be greater than 1.0')     
            
               elseif( words(1) == 'ATOMI' ) then 
                  if ( num_isotopes_neu(nmat)>1_ip ) then
                     At_weight_isotope(nmat,nisotope) = getrea('ATOMI',0.0_rp,'#At_weights')
                  else
                     At_weight(nmat) = getrea('ATOMI',0.0_rp,'#At_weights')
                  endif
               elseif( words(1) == 'DENSI' ) then 
                  Densidad_(nmat) = getrea('DENSI',0.0_rp,'#Densidad')
               ! else if( words(1) == 'ISOTO' ) then 
               !    Isotope(nmat) = getrea('ISOTO',0.0_rp,'#isotopes')
               else if( words(1)=='ABSOR')   then     
                  absor_neu_cte(nmat)=getrea('ABSOR',0.0_rp,'#absorciones')
               else if( words(1)=='SCATT')   then     
                  scatt_neu_cte(nmat)=getrea('SCATT',0.0_rp,'#scattering')
               elseif( words(1) == 'MASSP' ) then
                  if (atomic_percentage)  call runend('neu_reaphy: both mass percentage and atomic '//&
                     'percentege have been defined. Can only use one or the other.')
                  mass_percentage_isotope_neu(nmat,nisotope) = getrea('MASSP',0.0_rp,'#Mass percentege of component')
                  mass_percentage = .TRUE.
               elseif( words(1) == 'ATOMP' ) then 
                  if (mass_percentage) call runend('neu_reaphy: both mass percentage and atomic '//&
                     'percentege have been defined. Can only use one or the other.')
                  atom_percentage_isotope_neu(nmat,nisotope) = getrea('ATOMP',0.0_rp,'#Atomic percentege of component')
                  atomic_percentage = .TRUE.
               else if( words(1)=='TOTXS')   then      
      
                  fil_totalXS_neu(nmat)=getcha('TOTXS','TOTXS','#total XC file')
                  write(sm,'(i1)') nmat
                  if(fil_totalXS_neu(nmat)/='  ') then
                     fil_totalXS_neu(nmat)=TRIM(fil_totalXS_neu(nmat))//sm//'.dat'
                  endif
               else if( words(1)=='SCTXS')   then      
      
                  fil_scattXS_neu(nmat)=getcha('SCTXS','SCTXS','#Scatter XC file')
                  write(sm,'(i1)') nmat
                  if(fil_scattXS_neu(nmat)/='  ') then
                     fil_scattXS_neu(nmat)=TRIM(fil_scattXS_neu(nmat))//sm//'.dat'
                  endif

               ! else if( words(1)=='SOURF')   then      
      
               !    fil_source_neu(nmat)=getcha('SOURF','SOURF','#Source bound file')
               !    write(sm,'(i1)') nmat
               !    if(fil_source_neu(nmat)/='  ') then
               !       fil_source_neu(nmat)=TRIM(fil_source_neu(nmat))//sm//'.dat'
               !    endif

               else if( words(1)=='FISXS')   then      
      
                  fil_fisionXS_neu(nmat)=getcha('FISXS','FISXS','#Fision XC file')
                  write(sm,'(i1)') nmat
                  if(fil_fisionXS_neu(nmat)/='  ') then
                     fil_fisionXS_neu(nmat)=TRIM(fil_fisionXS_neu(nmat))//sm//'.dat'
                  endif

               else if( words(1)=='KERMA')   then      
      
                  fil_kerma_neu(nmat)=getcha('KERMA','KERMA','#Kerma factors file')
                  write(sm,'(i1)') nmat
                  if(fil_kerma_neu(nmat)/='  ') then
                     fil_kerma_neu(nmat)=TRIM(fil_kerma_neu(nmat))//sm//'.dat'
                  endif

               endif
               
               call ecoute('neu_reaphy')

            enddo

            ! If mass percentage in input, calculate atomic percentage of each isotope
            if (mass_percentage) then
               do nmat=1,num_materials_neu
                  if ( num_isotopes_neu(nmat) > 1_ip ) then
                     sum_aux = 0.0_rp
                     do nisotope=1,num_isotopes_neu(nmat)
                        atom_percentage_isotope_neu(nmat,nisotope) = &
                           mass_percentage_isotope_neu(nmat,nisotope)/At_weight_isotope(nmat,nisotope)
                        sum_aux = sum_aux + atom_percentage_isotope_neu(nmat,nisotope)
                     enddo
                     do nisotope=1,num_isotopes_neu(nmat)
                        atom_percentage_isotope_neu(nmat,nisotope) = atom_percentage_isotope_neu(nmat,nisotope)/sum_aux
                     enddo
                  endif
               enddo
            else ! viceversa 
               do nmat=1,num_materials_neu
                  if ( num_isotopes_neu(nmat) > 1_ip ) then
                     sum_aux = 0.0_rp
                     do nisotope=1,num_isotopes_neu(nmat)
                        mass_percentage_isotope_neu(nmat,nisotope) = &
                           atom_percentage_isotope_neu(nmat,nisotope)*At_weight_isotope(nmat,nisotope)
                        sum_aux = sum_aux + mass_percentage_isotope_neu(nmat,nisotope)
                     enddo
                     do nisotope=1,num_isotopes_neu(nmat)
                        mass_percentage_isotope_neu(nmat,nisotope) = mass_percentage_isotope_neu(nmat,nisotope)/sum_aux
                     enddo
                  endif
               enddo
            endif

            do nmat=1,num_materials_neu
               if ( num_isotopes_neu(nmat) > 1_ip ) then
                  sum_atom = 0.0_rp
                  sum_mass = 0.0_rp
                  At_weight(nmat) = 0.0_rp
                  ! Densidad_(nmat) = 0.0_rp
                  do nisotope=1,num_isotopes_neu(nmat)
                     sum_atom = sum_atom + atom_percentage_isotope_neu(nmat,nisotope)
                     sum_mass = sum_mass + mass_percentage_isotope_neu(nmat,nisotope)

                     At_weight(nmat) = At_weight(nmat) + mass_percentage_isotope_neu(nmat,nisotope)/At_weight_isotope(nmat,nisotope)
                     ! Densidad_(nmat) = Densidad_(nmat) + mass_percentage_isotope_neu(nmat,nisotope)/Densidad_comp(nmat,nisotope)
                  end do
                  ! Densidad_(1) = 1/Densidad_(1)
                  if (abs(sum_mass-1.0_rp)>1E-2_rp) call runend('neu_reaphy: error in sum of mass '//&
                     'percentages (not close enough to 1) in material '//trim(intost(nmat)))
                  if (abs(sum_atom-1.0_rp)>1E-2_rp) call runend('neu_reaphy: error in sum of atomic '//&
                     'percentages (not close enough to 1) in material '//trim(intost(nmat)))
                  
                  At_weight(nmat) = 1.0_rp/At_weight(nmat)

               end if
            enddo
            if ( num_isotopes_max_neu>1_ip ) deallocate(At_weight_isotope)

        end if
        
    
        call ecoute('neu_reaphy')


        
     end do
     !--></group>
     !
     ! ADOC[0]> END_PHYSICAL_PROBLEM
     !
     ! OJO esto debe masnadarse s todos los porcesos en parall!!

     allocate(absor_neu(num_materials_neu,num_energies_neu),grupo_energias(num_materials_neu,num_energies_neu))
     allocate(scatt_neu(num_materials_neu,num_energies_neu,num_energies_neu,num_legendre_lee+1))
     allocate(fiss_neu(num_materials_neu,num_energies_neu,num_energies_neu))
     allocate(source_bound_neu(num_sources_neu,num_energies_neu))
     allocate(efectivos_neu_in(num_energies_neu),efectivos_neu_out(num_energies_neu))
     allocate(kerma_neu(num_materials_neu, num_energies_neu))
     allocate(n_efectivos_neu_sparse(num_materials_neu, num_energies_neu), &
               n_efectivos_neu_sparse_acum(num_materials_neu, num_energies_neu+1))
     allocate(resid_energy_group_neu(num_energies_neu))

     absor_neu=0.0_rp
     scatt_neu=0.0_rp
     fiss_neu=0.0_rp
     grupo_energias=0.0_rp
     source_bound_neu=0.0_rp
     efectivos_neu_in=0_ip
     efectivos_neu_out=0_ip
     kerma_neu=0.0_rp

     n_efectivos_neu_sparse=0_ip
     n_efectivos_neu_sparse_acum = 0_ip
     n_efectivos_neu_sparse_acum(:,1)=1_ip
     resid_energy_group_neu=-1.0_rp

  end if
  
end subroutine neu_reaphy
