!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_drugs

    use def_master
    use def_domain 
    use def_exmedi
    use def_kermod
    
    implicit none

    private

    !The numbers IONIC_* reflect the order in which the doses are stored in DRUG_COLLECTION and read from exm.dat
    !these numbers should be sequential
    integer(ip), parameter, public ::&                    
        IONIC_CAL = 1_ip,&
        IONIC_KR  = 2_ip,&
        IONIC_NA  = 3_ip,&
        IONIC_K1  = 4_ip,&
        IONIC_NAL = 5_ip,&
        IONIC_KS  = 6_ip,&
        IONIC_ITO = 7_ip,&
        NVARS_DRUG = 7_ip  ! number of variables affected by the drug

    character(5), dimension(NVARS_DRUG), parameter :: &
        drug_var_names = ['ICaL ','KR   ','NA   ','IK1  ','INaL ','IKs  ','Ito  ']

    TYPE :: DRUG_CONDUCTANCES
        real(rp)    :: pca  = 0.0_rp
        real(rp)    :: gkr  = 0.0_rp
        real(rp)    :: gna  = 0.0_rp
        real(rp)    :: gk1  = 0.0_rp
        real(rp)    :: gnal = 0.0_rp
        real(rp)    :: gks  = 0.0_rp
        real(rp)    :: Ito  = 0.0_rp
    END TYPE DRUG_CONDUCTANCES



    type(DRUG_CONDUCTANCES) :: precalculated_drug_effect
    logical(lg)             :: drug_conductances_precalculated = .False.

    TYPE :: DRUG
        character(5) :: name
        real(rp), dimension(NVARS_DRUG) :: &
            drug_ic50,&       ! parameters to introduce a drug effect on cell models, 12xnmate
            drug_dosis,&      ! parameters to introduce a drug effect on cell models, 12xnmate
            drug_h            ! h value. As in g = g_control ( 1 + (dosis/ic50)^h )^(-1)  
    END TYPE DRUG


    TYPE, PUBLIC :: DRUG_COLLECTION
        integer(ip)                       :: ndrugs          ! number of drugs
        type(DRUG), dimension(:), pointer :: drugs => null() ! drug container

        contains
        procedure :: allocate          => DRUG_COLLECTION_allocate
        procedure :: print             => DRUG_COLLECTION_print
        procedure :: compare_to_file   => DRUG_COLLECTION_compare_to_file
        procedure :: save_to_file      => DRUG_COLLECTION_save_to_file
    END TYPE DRUG_COLLECTION

    type(DRUG_COLLECTION), dimension(:), pointer :: drugs_per_material => null() ! drug container for all the materials

    interface memory_alloca
    module procedure memory_alloca_drugcollection   ! drug collection per material
    module procedure memory_alloca_drugs            ! drugs in each colection
    module procedure memory_alloca_ohara_drugconduct
    end interface memory_alloca

    interface memory_deallo
    module procedure memory_deallo_drugcollection
    module procedure memory_deallo_drugs             
    module procedure memory_deallo_ohara_drugconduct
    end interface memory_deallo


    public :: DRUG_CONDUCTANCES

    public :: exm_drugs_exchange
    public :: exm_drugs_allocate
    public :: exm_drugs_ondrugs
    public :: exm_drugs_calculate_effect
    public :: exm_drugs_read_data
    public :: exm_drugs_get_drugdmate ! temporary, for TOR MODEL, DO NOT USE!!!
    public :: exm_drugs_write_log

    !Serialization
    public :: exm_drugs_save_to_file
    public :: exm_drugs_compare_to_file

    !Functions for the unit tests
    public :: unitt_exm_drugs_set      !set values of drugs
    public :: unitt_exm_drugs_allocate !allocate the arrays for unit tests, since arrays are private


    !memory allocation
    public :: memory_alloca
    public :: memory_deallo
    public :: memory_alloca_ohara_drugconduct
    public :: memory_deallo_ohara_drugconduct

contains

logical(lg) function exm_drugs_compare_to_file(file_handle)
    !!! Make sure integers are saved with precicion 8
    implicit none
    integer(ip), intent(in) :: file_handle
    integer(ip)             :: imate

    exm_drugs_compare_to_file = .TRUE.

    do imate=1,nmate
        exm_drugs_compare_to_file = exm_drugs_compare_to_file .and. drugs_per_material(imate) % compare_to_file( file_handle )
        if ( .not. exm_drugs_compare_to_file ) exit
    end do
end function exm_drugs_compare_to_file


subroutine exm_drugs_save_to_file(file_handle)
    !!! Make sure integers are saved with precicion 8
    implicit none
    integer(ip), intent(in) :: file_handle
    integer(ip)             :: imate

    do imate=1,nmate
        call drugs_per_material(imate) % save_to_file( file_handle )
    end do
end subroutine exm_drugs_save_to_file



logical(lg) function  DRUG_COLLECTION_compare_to_file(this, file_handle)
    !!! Make sure integers are read with precicion 8
    implicit none
    class(DRUG_COLLECTION)  :: this
    integer(ip), intent(in) :: file_handle
    integer(ip)             :: idrug
    type(DRUG)              :: drug_data
    integer(ip)             :: iostat
    DRUG_COLLECTION_compare_to_file = .TRUE.

    do idrug=1_ip, this % ndrugs    
        read(file_handle, iostat = iostat) drug_data % name 
        read(file_handle, iostat = iostat) drug_data % drug_ic50
        read(file_handle, iostat = iostat) drug_data % drug_dosis
        read(file_handle, iostat = iostat) drug_data % drug_h

        if ( iostat .ne. 0 ) then
            DRUG_COLLECTION_compare_to_file = .FALSE.
            exit
        end if

        DRUG_COLLECTION_compare_to_file = DRUG_COLLECTION_compare_to_file .and. ( this % drugs(idrug) % name == drug_data % name )
        DRUG_COLLECTION_compare_to_file = DRUG_COLLECTION_compare_to_file .and. all( abs( this % drugs(idrug) % drug_ic50  - drug_data % drug_ic50  ) < epsilon(1.0_rp) )
        DRUG_COLLECTION_compare_to_file = DRUG_COLLECTION_compare_to_file .and. all( abs( this % drugs(idrug) % drug_dosis - drug_data % drug_dosis ) < epsilon(1.0_rp) )
        DRUG_COLLECTION_compare_to_file = DRUG_COLLECTION_compare_to_file .and. all( abs( this % drugs(idrug) % drug_h     - drug_data % drug_h     ) < epsilon(1.0_rp) )
    end do
end function DRUG_COLLECTION_compare_to_file


subroutine DRUG_COLLECTION_save_to_file(this, file_handle)
    implicit none
    class(DRUG_COLLECTION)  :: this
    integer(ip), intent(in) :: file_handle
    integer(ip)             :: idrug

    do idrug=1_ip, this % ndrugs    
        write(file_handle) this % drugs(idrug) % name 
        write(file_handle) this % drugs(idrug) % drug_ic50
        write(file_handle) this % drugs(idrug) % drug_dosis
        write(file_handle) this % drugs(idrug) % drug_h
    end do
end subroutine DRUG_COLLECTION_save_to_file



subroutine DRUG_COLLECTION_print(this)
    implicit none
    class(DRUG_COLLECTION)  :: this
    integer(ip) :: idrug

    do idrug=1_ip, this % ndrugs
        print *,idrug,' DRUG NAME :',this % drugs(idrug) % name 
        print *,idrug,' IC50 :'
        print *,this % drugs(idrug) % drug_ic50
        print *,idrug,' DOSIS :'
        print *,this % drugs(idrug) % drug_dosis
        print *,idrug,' H :'
        print *,this % drugs(idrug) % drug_h
    end do
end subroutine DRUG_COLLECTION_print


subroutine DRUG_COLLECTION_allocate(this)
    implicit none
    class(DRUG_COLLECTION)  :: this
    integer(ip)             :: idrug

    call memory_alloca(mem_modul(1:2,modul),'DRUG_COLLECTION % drugs', 'mod_exm_drugs', this % drugs, this % ndrugs)

    do idrug=1_ip, this % ndrugs
        this % drugs(idrug) % name = ""
        this % drugs(idrug) % drug_ic50  = 0.0_rp
        this % drugs(idrug) % drug_dosis = 0.0_rp
        this % drugs(idrug) % drug_h     = 1.0_rp
    end do

    precalculated_drug_effect % pca  = 0.0_rp
    precalculated_drug_effect % gkr  = 0.0_rp
    precalculated_drug_effect % gna  = 0.0_rp
    precalculated_drug_effect % gk1  = 0.0_rp
    precalculated_drug_effect % gnal = 0.0_rp
    precalculated_drug_effect % gks  = 0.0_rp
    precalculated_drug_effect % Ito  = 0.0_rp
    
end subroutine DRUG_COLLECTION_allocate



subroutine exm_drugs_get_drugdmate(imate, drugdmate, nvars_per_drug, ndrugs)
    use mod_messages, only : messages_live

    ! Left here for compatibility with torord model
    !!! DO NOT USE !!!
    !Copy out first nvars of drugdmate
    implicit none
    integer(ip), intent(in) :: imate
    integer(ip), intent(in) :: nvars_per_drug
    integer(ip), intent(in) :: ndrugs
    real(rp), dimension(:), intent(inout) :: drugdmate
    integer(ip) :: i, idrug, ivar

    call messages_live("USING DEPRECATED exm_drugs_get_drugdmate, REFACTOR YOUR CODE TO USE OTHER FUNCTIONS TO ACCESS THE DRUG INFORMATION","WARNING")

    !flatten the arrays for TOR model until someone refactors TOR model
    if ( size(drugdmate,1_ip, kind=ip) < nvars_per_drug*ndrugs ) then
        call runend("exm_drugs_get_drugdmate: size(drugdmate) < nvars")
    end if

    i = 1_ip
    do idrug = 1_ip, min(ndrugs, drugs_per_material(imate) % ndrugs)
        do ivar = 1_ip, min(nvars_per_drug, NVARS_DRUG)
            drugdmate(i) = drugs_per_material(imate) % drugs(idrug) % drug_dosis(ivar) !drugdmate_exm(i,imate)
            i = i + 1_ip
            drugdmate(i) = drugs_per_material(imate) % drugs(idrug) % drug_ic50(ivar)  !drugdmate_exm(i,imate)
            i = i + 1_ip
        end do
    end do

end subroutine exm_drugs_get_drugdmate

logical(lg) function exm_drugs_ondrugs(imate)
    implicit none
    integer(ip), intent(in) :: imate

    exm_drugs_ondrugs = ( drugs_per_material(imate) % ndrugs > 0_ip )
end function exm_drugs_ondrugs


subroutine exm_drugs_allocate()
    use mod_exm_memory, only : memory_alloca
    implicit none
    integer(ip) :: imate

    call memory_alloca(mem_modul(1:2,modul), 'drugs_per_material', 'mod_exm_drugs', drugs_per_material, nmate)

    do imate = 1_ip, nmate
        drugs_per_material(imate) % ndrugs = 0_ip
        nullify( drugs_per_material(imate) % drugs )
    end do

    drug_conductances_precalculated = .False.
end subroutine exm_drugs_allocate


subroutine exm_drugs_exchange()
    use      mod_exchange,            only : exchange_add
    use      mod_exchange,            only : exchange_end
    use      mod_exchange,            only : exchange_init

    implicit none
    integer(ip) :: imate, idrug

    call exchange_init()
    do imate = 1_ip, nmate
        call exchange_add( drugs_per_material(imate) % ndrugs )
    end do

    call exchange_end() 
    call exchange_init()

    do imate = 1_ip, nmate
        if (ISLAVE) call drugs_per_material(imate) % allocate()

        do idrug = 1_ip, drugs_per_material(imate) % ndrugs
            call exchange_add( drugs_per_material(imate) % drugs(idrug) % name       )
            call exchange_add( drugs_per_material(imate) % drugs(idrug) % drug_ic50  )
            call exchange_add( drugs_per_material(imate) % drugs(idrug) % drug_dosis )
            call exchange_add( drugs_per_material(imate) % drugs(idrug) % drug_h     )
        end do
    end do

    call exchange_end()


end subroutine exm_drugs_exchange



subroutine exm_drugs_read_data(imate)
    use def_inpout,   only : param, exists, words, nnpar, nnwor, getint, getcha
    use mod_messages, only: messages_live
    use mod_ecoute,   only : ecoute

    implicit none
    integer(ip), intent(in) :: imate
    integer(ip)             :: ivalu
    integer(ip)             :: idrug

    drug_conductances_precalculated = .False.

    if(words(1)=='ONDRU') call runend("EXMEDI: Material "//trim(intost(imate))//" is using old style drug definition (Found ON_DRUGS)")
    

    if(words(1)=='DRUGS') then
        drugs_per_material(imate) % ndrugs = getint('NUMBE',0_ip,'!Number of drugs') 

        call ecoute('exm_reaphy') !skip to next line

        if ( drugs_per_material(imate) % ndrugs > 0_ip ) then
            call messages_live("           ON "//trim(intost(drugs_per_material(imate) % ndrugs))//" DRUGS")

            call drugs_per_material(imate) % allocate()

            idrug = 0_ip
            do while (words(1)/='ENDDR')
                idrug = idrug + 1_ip

                !for each drug read the table of pairs (dosis, ic50)
                if( words(1) == "DOSIS" ) then
                    drugs_per_material(imate) % drugs(idrug) % name = getcha('NAME ','     ','!Name of the drug') 
                    call ecoute('exm_reaphy') !skip to next line

                    ivalu = 0_ip
                    do while (words(1)/='ENDDO')
                        ivalu = ivalu + 1_ip

                        if ( nnpar < 2_ip .or. nnpar > 3_ip ) then 
                            call runend('exm_drugs_read_data: MATERIAL '//trim(intost(imate))//&
                            ', DRUG '//trim(intost(idrug))//&
                            ', ROW '//trim(intost(ivalu))//&
                            ': NUMBER OF VALUES IS '//trim(intost(nnpar))//' INSTEAD OF EXPECTED 2 OR 3 (dosis, ic50, h)')
                        end if  
                        if ( nnwor > 0_ip ) then 
                            call runend('exm_drugs_read_data: MATERIAL '//trim(intost(imate))//&
                            ', DRUG '//trim(intost(idrug))//&
                            ', ROW '//trim(intost(ivalu))//&
                            ': READ '//trim(intost(nnwor))//' WORDS. NO WORDS ARE EXPECTED. FIRST 3 WORDS ARE: '//words(1)//','//words(2)//','//words(3))
                        end if  

                        if ( any(param(1:nnpar)<0.0_rp) ) then 
                            call runend('EXM_REAPHY: MATERIAL '//trim(intost(imate))//', DRUG '//trim(drugs_per_material(imate) % drugs(idrug) % name)//&
                                ' CONTAINS NEGATIVE DOSIS, IC50, or H')
                        end if  

                        drugs_per_material(imate) % drugs(idrug) % drug_dosis(ivalu) = param(1)
                        drugs_per_material(imate) % drugs(idrug) % drug_ic50 (ivalu) = param(2)

                        if ( nnpar==3_ip ) then
                            drugs_per_material(imate) % drugs(idrug) % drug_h (ivalu) = param(3)
                        else 
                            drugs_per_material(imate) % drugs(idrug) % drug_h (ivalu) = 1.0_rp
                        end if

                        call ecoute('exm_reaphy')
                    end do

                    if (ivalu .ne. NVARS_DRUG) then
                        call runend('exm_drugs_read_data: MATERIAL '//trim(intost(imate))//&
                        ', DRUG '//trim(intost(idrug))//&
                        ' HAS '//trim(intost(ivalu))//' ROWS, INSTEAD OF EXPECTED '//trim(intost(NVARS_DRUG)) )
                    end if
                end if

                call ecoute('exm_reaphy')
            end do

            if ( idrug .ne. drugs_per_material(imate) % ndrugs ) then
                call runend("EXMEDI: Number of drugs read "//trim(intost(idrug))//&
                " is not equal to the number of drugs in DRUGS tag "//trim(intost(drugs_per_material(imate) % ndrugs)))
            end if
            
        else
            call ecoute('exm_reaphy')                   
        end if

        !call drugs_per_material(imate) % print()
    end if


end subroutine exm_drugs_read_data


type(DRUG_CONDUCTANCES) function exm_drugs_calculate_effect(imate)
    use mod_messages, only : messages_live
    implicit none
    integer(ip), intent(in) :: imate

    real(rp), dimension(NVARS_DRUG)    :: additive_effect
    integer(ip)                        :: ivar, idrug
    real(rp)                           :: g_drug
    real(rp)                           :: d_ic50

    if (.not. exm_drugs_ondrugs(imate)) then
        call runend("MOD_EXM_DRUGS: ATTEMPT TO CALCULATE DRUGS ON A MATERIAL "//trim(intost(imate))//" WITHOUT DRUGS")
    end if 

    if( .not. drug_conductances_precalculated ) then
        !additive effect
        additive_effect(:) = 1.0_rp

        do idrug = 1_ip, drugs_per_material(imate) % ndrugs
            do ivar = 1_ip, NVARS_DRUG
                if ( drugs_per_material(imate) % drugs(idrug) % drug_dosis(ivar)  > 0.00001_rp ) then
                    d_ic50 = drugs_per_material(imate) % drugs(idrug) % drug_dosis(ivar) / drugs_per_material(imate) % drugs(idrug) % drug_ic50(ivar)  
                    
                    !this if is to minimize the number of exponentiations, since real numer is never equal exactly 1
                    if ( abs(drugs_per_material(imate) % drugs(idrug) % drug_h(ivar) - 1.0_rp ) > 1.0e-5_rp ) then
                        g_drug = 1.0_rp / (1.0_rp + d_ic50 ** drugs_per_material(imate) % drugs(idrug) % drug_h(ivar)  )
                    else
                        g_drug = 1.0_rp / (1.0_rp + d_ic50 )
                    end if

                    additive_effect(ivar) = additive_effect(ivar) - (1.0_rp - g_drug)
                end if
            end do    
        end do

        precalculated_drug_effect % pca  = max( 0.0_rp, additive_effect(IONIC_CAL ) )
        precalculated_drug_effect % gkr  = max( 0.0_rp, additive_effect(IONIC_KR  ) )
        precalculated_drug_effect % gna  = max( 0.0_rp, additive_effect(IONIC_NA  ) )
        precalculated_drug_effect % gk1  = max( 0.0_rp, additive_effect(IONIC_K1  ) )
        precalculated_drug_effect % gnal = max( 0.0_rp, additive_effect(IONIC_NAL ) )
        precalculated_drug_effect % gks  = max( 0.0_rp, additive_effect(IONIC_KS  ) )
        precalculated_drug_effect % Ito  = max( 0.0_rp, additive_effect(IONIC_ITO ) )


        if      (precalculated_drug_effect % pca  < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: ICaL drug combination produced a negative or zero conductance for material '//trim(intost(imate)),'WARNING')
        else if (precalculated_drug_effect % gkr  < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: IKr drug combination produced a negative or zero conductance for material '//trim(intost(imate)),'WARNING')
        else if (precalculated_drug_effect % gna  < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: INa drug combination produced a negative or zero conductance for material '//trim(intost(imate)),'WARNING')
        else if (precalculated_drug_effect % gk1  < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: IK1 drug combination produced a negative or zero conductance for material '//trim(intost(imate)),'WARNING')
        else if (precalculated_drug_effect % gnal < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: INaL drug combination produced a negativeor zero conductance for material '//trim(intost(imate)),'WARNING')
        else if (precalculated_drug_effect % gks  < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: IKs drug combination produced a negative or zero conductance for material '//trim(intost(imate)),'WARNING')
        else if (precalculated_drug_effect % Ito  < epsilon(0.0_rp) ) then
            call messages_live('exm_oceohr: Ito drug combination produced a negative or zero conductance for material '//trim(intost(imate)),'WARNING')
        end if
    end if


    drug_conductances_precalculated = .True.
    exm_drugs_calculate_effect % pca  = precalculated_drug_effect % pca
    exm_drugs_calculate_effect % gkr  = precalculated_drug_effect % gkr
    exm_drugs_calculate_effect % gna  = precalculated_drug_effect % gna
    exm_drugs_calculate_effect % gk1  = precalculated_drug_effect % gk1
    exm_drugs_calculate_effect % gnal = precalculated_drug_effect % gnal
    exm_drugs_calculate_effect % gks  = precalculated_drug_effect % gks
    exm_drugs_calculate_effect % Ito  = precalculated_drug_effect % Ito
end function exm_drugs_calculate_effect






!-----------------------------------------------------------------------
!> 
!> @author  Constantine Butakoff
!> @date    2020-09-16
!> @brief   Allotate DRUG_COLLECTION
!> @details Allotate DRUG_COLLECTION
!> 
!-----------------------------------------------------------------------

subroutine memory_alloca_drugcollection(memor,vanam,vacal,varia,ndim1)
    use mod_memory,         only : lbytm, kfl_alloc
    use mod_memory_tools
    use mod_memory_basic
    implicit none

    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(DRUG_COLLECTION),  intent(inout), pointer        :: varia(:)      !< Variable
    integer(ip),  intent(in)                              :: ndim1         !< Dimension
    integer(4)                                            :: istat

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = storage_size(varia)*ndim1/8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type')

    else

       nullify(varia)

    end if

 end subroutine memory_alloca_drugcollection

 subroutine memory_deallo_drugcollection(memor,vanam,vacal,varia)
    use mod_memory,         only : lbytm
    use mod_memory_tools
    use mod_memory_basic
    implicit none

    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(DRUG_COLLECTION),  intent(inout), pointer        :: varia(:)      !< Variable
    integer(4)                                            :: istat
    
    if( associated(varia) ) then
       
       lbytm = -storage_size(varia)/8*size(varia,1,kind=ip)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type')

    else

       lbytm = 0
    
    end if

 end subroutine memory_deallo_drugcollection



!-----------------------------------------------------------------------
!> 
!> @author  Constantine Butakoff
!> @date    2020-09-16
!> @brief   Allocate array of DRUG
!> @details Allocate array of DRUG
!> 
!-----------------------------------------------------------------------

 subroutine memory_alloca_drugs(memor,vanam,vacal,varia,ndim1)
    use mod_memory,         only : lbytm, kfl_alloc
    use mod_memory_tools
    use mod_memory_basic
    implicit none

    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(DRUG),  intent(inout), pointer        :: varia(:)      !< Variable
    integer(ip),  intent(in)                              :: ndim1         !< Dimension
    integer(4)                                            :: istat

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = storage_size(varia)*ndim1/8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type')

    else

       nullify(varia)

    end if

 end subroutine memory_alloca_drugs

 subroutine memory_deallo_drugs(memor,vanam,vacal,varia)
    use mod_memory,         only : lbytm
    use mod_memory_tools
    use mod_memory_basic
    implicit none

    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(DRUG),  intent(inout), pointer        :: varia(:)      !< Variable
    integer(4)                                            :: istat
    
    if( associated(varia) ) then
       
       lbytm = -storage_size(varia)/8*size(varia,1,kind=ip)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type')

    else

       lbytm = 0
    
    end if

 end subroutine memory_deallo_drugs




 subroutine exm_drugs_write_log()
    use mod_outfor
    implicit none


    integer(ip) :: imate, idrug, ivar

        call outfor(25_ip,momod(modul) % lun_outpu,'DRUGS APPLIED TO MATERIALS')
        write(momod(modul)%lun_outpu,*) NEW_LINE('a')

        do imate=1,nmate
            write(momod(modul)%lun_outpu,*) '   Material '//trim(intost(imate))
            if ( drugs_per_material(imate) % ndrugs > 0_ip ) then
                do idrug = 1_ip, drugs_per_material(imate) % ndrugs
                    write(momod(modul)%lun_outpu,*) '      ', drugs_per_material(imate) % drugs(idrug) % name
                    write(momod(modul)%lun_outpu,*) '      DOSIS, IC50, H'
                    do ivar = 1_ip, NVARS_DRUG
                        write(momod(modul)%lun_outpu,*) '      '//drug_var_names(ivar)//                        &
                            " : "//trim(retost(drugs_per_material(imate) % drugs(idrug) % drug_dosis(ivar)))//  &
                            ", "//trim(retost(drugs_per_material(imate) % drugs(idrug) % drug_ic50(ivar)))//    &
                            ", "//trim(retost(drugs_per_material(imate) % drugs(idrug) % drug_h(ivar)))
                    end do
                    write(momod(modul)%lun_outpu,*) '      '
                end do
            else
                write(momod(modul)%lun_outpu,*) '      None'
                write(momod(modul)%lun_outpu,*) '      '
            end if
        end do

        write(momod(modul)%lun_outpu,*) NEW_LINE('a')

end subroutine exm_drugs_write_log




!--------------------------------------------------------------------
!
!   Functions for the unit tests
!
!-------------------------------------------------------------------
subroutine unitt_exm_drugs_set(imate, idrug, name, dosis, ic50, h)
    implicit none
    integer(ip), intent(in)  :: imate, idrug
    character(5), intent(in) :: name
    real(rp), dimension(NVARS_DRUG), intent(in) :: dosis, ic50, h

    drug_conductances_precalculated = .False.
    
    drugs_per_material(imate) % drugs(idrug) % name = name
    drugs_per_material(imate) % drugs(idrug) % drug_dosis(:) = dosis(:)
    drugs_per_material(imate) % drugs(idrug) % drug_ic50(:)  = ic50(:)
    drugs_per_material(imate) % drugs(idrug) % drug_h(:)     = h(:)
end subroutine unitt_exm_drugs_set

subroutine unitt_exm_drugs_allocate(ndrugs)
    !set nmate before calling
    implicit none
    integer(ip), intent(in) :: ndrugs
    integer(ip)             :: imate

    call memory_alloca(mem_modul(1:2,modul), 'drugs_per_material', 'mod_exm_drugs', drugs_per_material, nmate)

    do imate = 1_ip, nmate
        drugs_per_material(imate) % ndrugs = ndrugs
        nullify( drugs_per_material(imate) % drugs )
        call drugs_per_material(imate) % allocate()
    end do 
end subroutine unitt_exm_drugs_allocate









subroutine memory_alloca_ohara_drugconduct(memor,vanam,vacal,varia,ndim1)
    use mod_memory,         only : lbytm, kfl_alloc
    use mod_memory_tools
    use mod_memory_basic
    implicit none
    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(DRUG_CONDUCTANCES),  intent(inout), pointer      :: varia(:)      !< Variable
    integer(ip),  intent(in)                              :: ndim1         !< Dimension
    integer(4)                                            :: istat

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = storage_size(varia)*ndim1/8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type')

    else

       nullify(varia)

    end if

 end subroutine memory_alloca_ohara_drugconduct



 subroutine memory_deallo_ohara_drugconduct(memor,vanam,vacal,varia)
    use mod_memory,         only : lbytm
    use mod_memory_tools
    use mod_memory_basic
    implicit none
    character(*), intent(in)                              :: vanam         !< Variable name
    character(*), intent(in)                              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)                           :: memor(2)      !< Memory counter
    type(DRUG_CONDUCTANCES),  intent(inout), pointer      :: varia(:)      !< Variable
    integer(4)                                            :: istat
    
    if( associated(varia) ) then
       
       lbytm = -storage_size(varia)/8*size(varia,1,kind=ip)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type')

    else

       lbytm = 0
    
    end if

 end subroutine memory_deallo_ohara_drugconduct

end module mod_exm_drugs
  
