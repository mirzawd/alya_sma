!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_cellmodel
    use def_master, only : rp, ip

    implicit none

    integer(ip), parameter  :: &
        EXM_CELLMODEL_CONVERGED      = 0_ip, &
        EXM_CELLMODEL_NOTCONVERGED   = 1_ip, &
        EXM_CELLMODEL_NOTINITIALIZED = 2_ip, &
        EXM_CELLMODEL_LOADED         = 3_ip

    integer(ip), parameter  ::&
        EXM_CELL_STEADY_VOLTAGE=0_ip, &                          !variables to decide on steady state in cell model
        EXM_CELL_STEADY_CALCIUM=1_ip

    integer(ip), parameter  ::&  ! for kfl_hfmodmate_exm
        EXM_MYOCYTE_NORMAL   = 0_ip, &
        EXM_MYOCYTE_MODIFIED = 1_ip

    !integer(ip), parameter  ::&  ! for kfl_hfmod_exm
    !    EXM_MYOCYTE_NONE         = 0_ip, &
    !    EXM_MYOCYTE_PIG          = 1_ip, &
    !    EXM_MYOCYTE_MALE         = 2_ip, &
    !    EXM_MYOCYTE_FEMALE       = 3_ip, &
    !    EXM_MYOCYTE_MODIFIED_PIG = 4_ip, &
    !    EXM_MYOCYTE_TABLE        = 10_ip   
        
    integer(ip), parameter  ::&
        EXM_CELLTYPE_ENDO = 1_ip, &
        EXM_CELLTYPE_MID  = 2_ip, &                                     !3-Epicardium, 1-endocardium, 2-midmyocardium
        EXM_CELLTYPE_EPI  = 3_ip 

    integer(ip), parameter ::&   !defined for Courtemanche model
        EXM_CELLTYPE_RA = 1_ip, &
        EXM_CELLTYPE_CTBBRA = 2_ip, &
        EXM_CELLTYPE_BBLA = 3_ip, & 
        EXM_CELLTYPE_TVR = 4_ip, & 
        EXM_CELLTYPE_MVR = 5_ip, & 
        EXM_CELLTYPE_RAA = 6_ip, & 
        EXM_CELLTYPE_LAA = 7_ip, & 
        EXM_CELLTYPE_LA = 8_ip, & 
        EXM_CELLTYPE_PV = 9_ip 
 
    TYPE :: CELL_MODEL_OUTPUTS
        real(rp)    :: S =           0.0_rp
        real(rp)    :: W =           0.0_rp
        real(rp)    :: CaTRPN =      0.0_rp
        real(rp)    :: B =           0.0_rp
        real(rp)    :: zeta_s =      0.0_rp
        real(rp)    :: zeta_w =      0.0_rp
        real(rp)    :: Ca50 =        0.805_rp
        real(rp)    :: Lambda =      1.0_rp     ! Lambda for solidz
        real(rp)    :: Vinit  =      0.0_rp     ! initial voltage
        integer(ip) :: nbeats =     -1_ip       ! number of taken to converge
        real(rp)    :: toler =      -1.0_rp     ! tolerance to decide if the cell model converged
        real(rp)    :: rmse =       -1.0_rp     ! rmse between the last two beats of the cell model use to decide convergence
        integer(ip) :: success =    -1_ip       ! if not 0, ohara cell model failed
        real(rp)    :: dt =         -1.0_rp     ! cell model timestep used, to save in log
    END TYPE CELL_MODEL_OUTPUTS

contains

subroutine exm_init_ncelltypes( array )
    use mod_eccoupling
    implicit none
    integer(ip), dimension(:) :: array

    if( size(array, 1, kind=ip)<EXMSLD_CELL_MAXMODELID ) call runend("exm_init_ncelltypes: argument has too few elements")

    array(EXMSLD_CELL_FITZHUGH    ) = 1_ip
    array(EXMSLD_CELL_FENTON      ) = 1_ip
    array(EXMSLD_CELL_TENTUSCHER  ) = 1_ip
    array(EXMSLD_CELL_TT2006      ) = 1_ip
    array(EXMSLD_CELL_OHARA       ) = 3_ip
    array(EXMSLD_CELL_OHARA_INAPA ) = 3_ip
    array(EXMSLD_CELL_SCATRIA     ) = 1_ip
    array(EXMSLD_CELL_SCVENTRI    ) = 1_ip        
    array(EXMSLD_CELL_TORORD      ) = 3_ip        
    array(EXMSLD_CELL_COURTE      ) = 9_ip     

end subroutine  



function exm_get_conductance_name(index) result(name)
    use def_kintyp, only : ip
    implicit none
    integer(ip), intent(in) :: index
    character(10)            :: name

    name = '     '

    select case (index) 
    case (1_ip)          
        name = 'GIto      ' 
    case (2_ip)          
        name = 'GKs       ' 
    case (3_ip)          
        name = 'GK1       ' 
    case (4_ip)          
        name = 'GKr       ' 
    case (5_ip)          
        name = 'GNa       ' 
    case (6_ip)          
        name = 'GNaL      ' 
    case (7_ip)          
        name = 'GNaCa     ' 
    case (8_ip)          
        name = 'gKb       ' 
    case (9_ip)          
        name = 'pCa       ' 
    case (10_ip)          
        name = 'pNaK      ' 
    case (11_ip)          
        name = 'Calmo     ' 
    case (12_ip)          
        name = 'Jup       ' 
    case (13_ip)          
        name = 'NaCaK     ' 
    case (14_ip)          
        name = 'Ikatp     ' 
    case (15_ip)          
        name = 'HFjlk     ' 
    case (16_ip)          
        name = 'JrNP      ' 
    case (17_ip)          
        name = 'Inal_tauHL' 
    case default          
        call runend("exm_get_conductance_name: index out of range")
    end select 

end function exm_get_conductance_name

subroutine exm_ohara_conductances_write_log()
    use def_master,     only : momod,modul
    use def_domain,     only : nmate
    use mod_outfor,     only : outfor
    use mod_eccoupling, only : eccou_get_cellmod, EXMSLD_CELL_NOMODEL
    use def_exmedi,     only : ttparmate_exm
    use mod_memory,     only : memory_size
    use def_master,     only : retost
    use def_master,     only : intost
    implicit none


    integer(ip) :: imate, ivar

        call outfor(25_ip,momod(modul) % lun_outpu,'CONDUCTANCES APPLIED TO MATERIALS (AS READ FROM .DAT, DRUG EFFECT NOT INCLUDED)')
        write(momod(modul)%lun_outpu,*) NEW_LINE('a')

        do imate=1,nmate
            write(momod(modul)%lun_outpu,*) '   Material '//trim(intost(imate))
            if ( eccou_get_cellmod(imate) .ne. EXMSLD_CELL_NOMODEL ) then
                do ivar = 1_ip, memory_size(ttparmate_exm,2_ip)
                    write(momod(modul)%lun_outpu,*) '      '//exm_get_conductance_name(ivar)//" : "//trim(retost(ttparmate_exm(1_ip,ivar, imate))) &
                                                                                            // ", "//trim(retost(ttparmate_exm(2_ip,ivar, imate))) &   
                                                                                            // ", "//trim(retost(ttparmate_exm(3_ip,ivar, imate)))
                end do
                write(momod(modul)%lun_outpu,*) '      '
            else
                write(momod(modul)%lun_outpu,*) '      None (NOMODEL)'
                write(momod(modul)%lun_outpu,*) '      '
            end if
        end do

        write(momod(modul)%lun_outpu,*) NEW_LINE('a')
end subroutine exm_ohara_conductances_write_log


end module mod_exm_cellmodel
