!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @addtogroup Fiber manipulation
!> @{
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @date    March, 2020
!> @}
!------------------------------------------------------------------------------
module mod_sld_fibers
    ! =========================================================================
    use def_kintyp, only                    :  ip, rp, lg
    use def_solidz, only                    :  kfl_fiber_sld
    ! ----------------------------------------
    implicit none
    ! ----------------------------------------
    type bio_fibers
        logical(lg)                         :: initialised
        integer(ip)                         :: location  
        integer(ip)                         :: field_fiber
        integer(ip)                         :: field_sheet
        integer(ip)                         :: field_norma
        real(rp),          dimension(3)     :: vector
        real(rp),          dimension(3,3)   :: basis
        real(rp), pointer, dimension(:,:)   :: fiber
        real(rp), pointer, dimension(:,:)   :: sheet
        real(rp), pointer, dimension(:,:)   :: norma
    end type 
    ! ----------------------------------------
    integer(ip),  parameter                 :: FIBER_NODAL     = 1_ip
    integer(ip),  parameter                 :: FIBER_ELEMENTAL = 2_ip
    integer(ip),  parameter                 :: FIBER_VECTOR    = 3_ip
    integer(ip),  parameter                 :: FIBER_BASIS     = 4_ip
    ! ----------------------------------------
    type(bio_fibers)                        :: fibers_sld
    ! =========================================================================
    contains 

        ! ---------------------------------------------------------------------
        !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date    
        !> @brief    Read input variables concerning the fiber model
        !> @details  
        !> @todo     Improve details explanation.
        ! ---------------------------------------------------------------------
        subroutine sld_fib_read_data( &
            imate )
            use def_domain,           only :  ndime
            use def_inpout,           only :  words, param, exists, getint, getrea, getcha
            use def_master,           only :  intost
            use mod_ecoute,           only :  ecoute
            use mod_messages,         only :  livinf, messages_live
            use def_kermod,           only :  kfl_fiber_long_fun
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)        :: imate
            ! -------------------------------
            call messages_live('   READING BIO-FIBERS MODEL')

            ! Initialize data
            fibers_sld % location    = 0_ip
            fibers_sld % field_fiber = 0_ip
            fibers_sld % field_sheet = 0_ip
            fibers_sld % field_norma = 0_ip

            ! Set key-flag indicating at least one model
            kfl_fiber_sld = 1_ip 
        
            !this function is called without going to the next line. Forcing read of the next line here
            call ecoute('sld_fib_read_data')

            ! Read the fiber model options
            do while(words(1)/='ENDFI')

                if( words(1) == 'LOCAT' )then
                    if( words(2) == 'NODAL' )then
                        fibers_sld % location = FIBER_NODAL
                        call messages_live('      FIBERS DEFINED AT NODES')

                    elseif( words(2) == 'ELEME' )then
                        fibers_sld % location = FIBER_ELEMENTAL 
                        call messages_live('      FIBERS DEFINED AT ELEMENTS')

                    endif
                    
                elseif( words(1) == 'VECTO' )then
                    fibers_sld % location = FIBER_VECTOR
                    fibers_sld % vector(1:3) = param(1:3)
                    call messages_live('      FIBERS DEFINED THROUGH A VECTOR')
               
                elseif( words(1) == 'BASIS' )then
                    fibers_sld % location = FIBER_BASIS
                    fibers_sld % basis(1:3,1) = param(1:3)
                    fibers_sld % basis(1:3,2) = param(4:6)
                    fibers_sld % basis(1:3,3) = param(7:9)
                    call messages_live('      FIBERS DEFINED THROUGH A BASIS')   

                elseif( words(1) == 'FIBER' )then
                    fibers_sld % field_fiber = int(param(1)) 
                    call messages_live('      FIBER VECTORS READ FROM FIELD '//intost(fibers_sld % field_fiber))
                    if ( fibers_sld % field_fiber .ne. -kfl_fiber_long_fun ) then
                        call messages_live("SOLIDS FIBER FIELD "//trim(intost(fibers_sld % field_fiber))//" IS DIFFERENT FROM THE KERNEL FIBER FIELD "//trim(intost(-kfl_fiber_long_fun)),'WARNING')
                    end if
                              
                elseif( words(1) == 'SHEET' )then
                    fibers_sld % field_sheet = int(param(1))
                    call messages_live('      SHEET VECTORS READ FROM FIELD '//intost(fibers_sld % field_sheet))

                elseif( words(1) == 'NORMA' )then
                    fibers_sld % field_norma = int(param(1)) 
                    call messages_live('      NORMA VECTORS READ FROM FIELD '//intost(fibers_sld % field_norma))

                end if

                call ecoute('sld_fib_read_data')
  
            enddo

            ! Set that the bio-fiber model has been initilised
            fibers_sld % initialised = .true.

            ! Perform some comprovations
            if( fibers_sld % initialised )then
                if( fibers_sld % location == 0_ip )then
                    call runend('SLD_FIB_READ_DATA: LOCATION ONLY CAN BE NODAL OR ELEMENTAL OR BASIS OR VECTOR')
                elseif( fibers_sld % location == FIBER_VECTOR )then
                    if( sum(fibers_sld % vector(:)) == 0_rp )then
                        call runend('SLD_FIB_READ_DATA: THE VECTOR MUST BE DEFINED ') 
                    endif
                elseif( fibers_sld % location == FIBER_BASIS )then
                    ! TODO : DEFINE THE BASIS CHECKING
                elseif( fibers_sld % location == FIBER_ELEMENTAL .or. fibers_sld % location == FIBER_NODAL  )then
                    if( fibers_sld % field_fiber == 0_ip )then
                        call runend('SLD_FIB_READ_DATA: AT LEAST A FIELD FOR THE FIBERS HAS TO BE DEFINED') 
                    endif
                endif
                if( ndime /= 3_ip )then
                    call runend('SLD_FIB_READ_DATA: BIO-FIBERS ONLY AVAILABLE FOR 3D ANALYSIS')
                endif
            endif            
            ! -------------------------------
        end subroutine sld_fib_read_data


        subroutine sld_fib_read_data_warning()
            use mod_messages,         only :  livinf, messages_live
            ! -------------------------------
            call messages_live("-------------------------------------------------------","WARNING")
            call messages_live(" A deprecated syntax is used to define the FIBER_MODEL:",'WARNING')
            call messages_live("    FIBER_MODEL FIELD=1                                ",'WARNING')
            call messages_live("    ORTHOTROPIC_MODEL FIELD=1,2                        ",'WARNING')
            call messages_live(" Instead of the old syntax, use the following one:     ","WARNING")
            call messages_live("     FIBER_MODEL                                       ","WARNING")
            call messages_live("         LOCATION: NODAL | ELEMENTAL                   ","WARNING")
            call messages_live("         FIBER:    1     $ (mandatory ID_FIELD)        ","WARNING")
            call messages_live("         SHEET:    2     $ (optional ID_FIELD)         ","WARNING")
            call messages_live("         NORMA:    3     $ (optional ID_FIELD)         ","WARNING")
            call messages_live("     END_FIBER_MODEL                                   ","WARNING")
            call messages_live("-------------------------------------------------------","WARNING")
            call runend('SLD_FIB_READ_DATA: READ THE PREVIOUS WARNINGS')
            ! -------------------------------
        end subroutine

        ! ---------------------------------------------------------------------
        !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date     July 2020 
        !> @brief    Send data read in the sld file from master to slaves 
        !> @details  
        !> @todo     
        ! ---------------------------------------------------------------------
        subroutine sld_fib_send_data(&
             itask )
          ! -------------------------------
          use mod_exchange, only : exchange_add
          implicit none
          ! -------------------------------
          integer(ip), intent(in)        :: itask
          ! -------------------------------
          select case(itask)

          case(1_ip)
             !
             ! Interchange flags 
             !
             !call iexcha( kfl_fiber_sld )

          case(2_ip)
             ! 
             ! Interchange readed data
             !
             if( kfl_fiber_sld == 1_ip )then
                call exchange_add( fibers_sld % initialised )
                call exchange_add( fibers_sld % location    )
                call exchange_add( fibers_sld % field_fiber )
                call exchange_add( fibers_sld % field_sheet )
                call exchange_add( fibers_sld % field_norma ) 
                call exchange_add( fibers_sld % vector      )
                call exchange_add( fibers_sld % basis       )
             endif

          end select
          ! -------------------------------
        end subroutine sld_fib_send_data

        ! ---------------------------------------------------------------------
        !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date     July 2020 
        !> @brief    Send data read in the sld file from master to slaves 
        !> @details  
        !> @todo     Delete fbts_sld and fibtn_sld     
        ! ---------------------------------------------------------------------
        subroutine sld_fib_initialise(&
            )
            use def_master,           only  : INOTMASTER, fiber
            use def_domain,           only  : xfiel  
            use def_solidz,           only  : fibts_sld, fibtn_sld
            ! -------------------------------
            implicit none
            ! -------------------------------
            if( INOTMASTER .and. kfl_fiber_sld == 1_ip ) then
       
                ! Associate fiber direction to the corresponding field
                if( fibers_sld % field_fiber > 0_ip )then
                    fibers_sld % fiber => xfiel(fibers_sld % field_fiber) % a(:,:,1)
                    call normalize(fibers_sld % fiber) 
                endif
                 
                ! Associate sheet direction to the corresponding field
                if( fibers_sld % field_sheet > 0_ip )then
                    fibers_sld % sheet => xfiel(fibers_sld % field_sheet) % a(:,:,1)
                    call normalize(fibers_sld % sheet) 
                endif
          
                ! Associate norma direction to the corresponding field
                if( fibers_sld % field_norma > 0_ip )then
                    fibers_sld % norma => xfiel(fibers_sld % field_norma) % a(:,:,1)
                    call normalize(fibers_sld % norma) 
                endif
 
            end if     
            ! -------------------------------
            ! TODO : Delete this point done for retrocompatibility
            fiber => fibers_sld % fiber
            if( fibers_sld % field_sheet > 0_ip  )then
                fibts_sld => fibers_sld % sheet
            endif
            if( fibers_sld % field_norma > 0_ip )then
                fibtn_sld => fibers_sld % norma
            endif
            ! -------------------------------
            contains 

            subroutine normalize(array)
                ! ------------------------
                implicit none
                ! ------------------------
                real(rp), intent(inout) :: array(:,:) 
                real(rp)                :: norm
                integer(ip)             :: ii, jj
                ! ------------------------
                do jj = 1, size(array,dim=2)
                    norm = 0.0_rp
                    do ii = 1, size(array,dim=1)
                        norm = norm + array(ii,jj)**2
                    enddo
                    if( abs(norm) > 0.0_rp )then
                        array(:,jj) = array(:,jj)/sqrt(norm)
                    endif
                enddo
                ! ------------------------
            end subroutine normalize            
            ! -------------------------------
        end subroutine sld_fib_initialise

        ! ---------------------------------------------------------------------
        !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
        !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
        !> @date     July 2020 
        !> @brief    Gather fibers according to the model
        !> @details  
        ! ---------------------------------------------------------------------
        subroutine sld_fib_get_normalized_fiber_directions_at_gauss_points( &
            ielem, pnode, pgaus, lnods, gp_shape, gp_fiber, gp_sheet, gp_norma)
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)        :: ielem
            integer(ip), intent(in)        :: pnode
            integer(ip), intent(in)        :: pgaus
            integer(ip), intent(in)        :: lnods(pnode)
            real(rp),    intent(in)        :: gp_shape(pnode,pgaus)
            real(rp),    intent(out)       :: gp_fiber(3,pgaus)
            real(rp),    intent(out)       :: gp_sheet(3,pgaus)
            real(rp),    intent(out)       :: gp_norma(3,pgaus)
            real(rp)                       :: el_fiber(3,pnode)
            real(rp)                       :: el_sheet(3,pnode)
            real(rp)                       :: el_norma(3,pnode)
            integer(ip)                    :: ii
            ! -------------------------------
            ! Initialize
            gp_fiber = 0.0_rp 
            gp_sheet = 0.0_rp 
            gp_norma = 0.0_rp 

            ! Gather according to the location
            if( fibers_sld % initialised  )then 
 
                if(     fibers_sld % location == FIBER_VECTOR )then 
                    do ii = 1, pgaus
                        gp_fiber(1:3,ii) = fibers_sld % vector(1:3)
                    enddo

                elseif( fibers_sld % location == FIBER_BASIS )then
                    do ii = 1, pgaus
                        gp_fiber(1:3,ii) = fibers_sld % basis(1:3,1)
                    enddo
                    do ii = 1, pgaus
                        gp_sheet(1:3,ii) = fibers_sld % basis(1:3,2)
                    enddo
                    do ii = 1, pgaus
                        gp_norma(1:3,ii) = fibers_sld % basis(1:3,3)
                    enddo

                elseif( fibers_sld % location == FIBER_ELEMENTAL )then
                    if( fibers_sld % field_fiber > 0_ip)then
                        do ii = 1, pgaus
                            gp_fiber(1:3,ii) = fibers_sld % fiber(1:3,ielem)
                        enddo
                    endif
                    if( fibers_sld % field_sheet > 0_ip)then
                        do ii = 1, pgaus
                            gp_sheet(1:3,ii) = fibers_sld % sheet(1:3,ielem)
                        enddo
                    endif
                    if( fibers_sld % field_norma > 0_ip)then
                        do ii = 1, pgaus
                            gp_norma(1:3,ii) = fibers_sld % norma(1:3,ielem)
                        enddo
                    endif
                   
                elseif( fibers_sld % location == FIBER_NODAL )then
                    if( fibers_sld % field_fiber > 0_ip )then
                        el_fiber = gather_from_alya_to_array(lnods,pnode,fibers_sld % fiber)
                        gp_fiber = interpolate_and_normalize(pnode,pgaus,gp_shape,el_fiber)
                    endif
                    if( fibers_sld % field_sheet > 0_ip )then
                        el_sheet = gather_from_alya_to_array(lnods,pnode,fibers_sld % sheet)
                        gp_sheet = interpolate_and_normalize(pnode,pgaus,gp_shape,el_sheet)
                    endif 
                    if( fibers_sld % field_norma > 0_ip )then
                        el_norma = gather_from_alya_to_array(lnods,pnode,fibers_sld % norma)
                        gp_norma = interpolate_and_normalize(pnode,pgaus,gp_shape,el_norma)
                    endif

                endif

            endif
            ! -------------------------------
            contains

                function gather_from_alya_to_array(lnods,nnode,alyavar) result(elemarr)
                    ! ------------------------
                    implicit none
                    ! ------------------------
                    integer(ip), intent(in) :: nnode
                    integer(ip), intent(in) :: lnods(:)
                    real(rp),    intent(in) :: alyavar(:,:)
                    real(rp)                :: elemarr(3,nnode)
                    integer(ip)             :: inode, ipoin
                    ! ------------------------
                    do inode = 1,nnode
                        ipoin = lnods(inode)
                        elemarr(1:3,inode) = alyavar(1:3,ipoin)
                    end do
                    ! ------------------------
                end function gather_from_alya_to_array

                function interpolate_and_normalize(nnode,ngaus,gpsha,var_el) result(var_gp)
                    ! ------------------------
                    implicit none
                    ! ------------------------
                    integer(ip), intent(in) :: nnode, ngaus
                    real(rp),    intent(in) :: gpsha(:,:)
                    real(rp),    intent(in) :: var_el(:,:)
                    real(rp)                :: var_gp(3,ngaus)
                    real(rp)                :: norm
                    integer(ip)             :: idof, inode, igaus
                    ! ------------------------
                    var_gp(:,:) = 0.0_rp
                    do igaus = 1, ngaus
                        ! Interpolate
                        do inode = 1, pnode
                            do idof = 1, 3
                                var_gp(idof,igaus) = var_gp(idof,igaus) + gpsha(inode,igaus)*var_el(idof,inode)
                            enddo
                        end do
                        ! Normalise
                        norm = 0.0_rp
                        do idof = 1, 3
                            norm = norm + var_gp(idof,igaus)**2
                        enddo
                        if( abs(norm) > 0.0_rp )then
                            var_gp(:,igaus) = var_gp(:,igaus)/sqrt(norm)
                        endif
                    enddo
                    ! ------------------------
                end function interpolate_and_normalize
            
            ! -------------------------------
        end subroutine sld_fib_get_normalized_fiber_directions_at_gauss_points

    ! =========================================================================
end module mod_sld_fibers
