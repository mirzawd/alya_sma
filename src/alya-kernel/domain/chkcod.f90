!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    chkcod.f90
!> @author  houzeaux
!> @date    2018-12-04
!> @brief   Boundary conditions codes on nodes
!> @details Output the different combinationsof boundary conditions
!>          on nodes
!> @} 
!-----------------------------------------------------------------------

subroutine chkcod()

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use mod_memory
  use mod_parall
  use mod_communications
  use mod_outfor,         only : outfor
  use mod_std
  implicit none
  
  integer(ip)          :: ipoin,icono,icode,kcomb,icomb,mcod1
  integer(ip)          :: icomb_par
  logical(lg)          :: lcodt(-mcodb-1:mcodb+1)
  integer(ip)          :: ncomb_par
  integer(ip)          :: ncomb
  integer(4)           :: ncomb4
  integer(4),  pointer :: ncomb_gat4(:) 
  integer(ip), pointer :: lcomb(:,:) 
  integer(ip), pointer :: lcomb_gat(:,:) 

  if( kfl_icodn == 0 ) return

  nullify(lcomb)
  nullify(lcomb_gat)
  nullify(ncomb_gat4)

  !----------------------------------------------------------------------
  !
  ! Compute possible combinations LCOMB
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     !
     ! NCOMB: Number of possible combinations
     !
     mcod1 = mcodb + 1
     lcodt = .false.
     do ipoin = 1,npoin
        do icono = 1,mcono
           icode = kfl_codno(icono,ipoin)
           lcodt(icode) = .true.           
        end do
     end do
     
     ncomb = count(lcodt(-mcodb:mcodb),KIND=ip)

     ncomb =    ncomb &
          &  +  ncomb * ( ncomb - 1 ) / 2 &
          &  +  ncomb * ( ncomb - 1 ) * ( ncomb - 2 ) / 3
     
     call memory_alloca(memor_dom,'LCOMB','chkcod',lcomb,mcono,ncomb)
     !
     ! Fill in combination table
     !
     if( ncomb > 0 ) then
        !
        ! Fill in combination table
        !
        lcomb = mcodb + 1
        kcomb = 0

        do ipoin = 1,npoin

           if( kfl_codno(1,ipoin) /= mcod1 ) then
              icomb = 1
              icode = lcomb(1,1) 
              do while( icode /= mcodb + 1 )
                 if(    kfl_codno(1,ipoin) == lcomb(1,icomb) .and. &
                      & kfl_codno(2,ipoin) == lcomb(2,icomb) .and. &
                      & kfl_codno(3,ipoin) == lcomb(3,icomb) ) then
                    icomb = -1
                    icode = mcodb + 1
                 else
                    icomb = icomb + 1
                    icode = lcomb(1,icomb)
                 end if
              end do
              !
              ! One new combination has been found
              !
              if( icomb /= -1 ) then
                 kcomb            = kcomb + 1
                 lcomb(1:3,kcomb) = kfl_codno(1:3,ipoin) 
              end if

           end if

        end do
        !
        ! Sort combinations
        !
        ncomb = kcomb
        do icomb = 1,ncomb
           if( lcomb(2,icomb) == mcod1 ) lcomb(2,icomb) = -mcod1 
           if( lcomb(3,icomb) == mcod1 ) lcomb(3,icomb) = -mcod1 
        end do
        call memory_resize(memor_dom,'LCOMB','chkcod',lcomb,mcono,ncomb)
        if( ncomb > 0 ) call sorti3(ncomb,lcomb)

     end if
     
  end if

  !----------------------------------------------------------------------
  !
  ! Gather number of combinations NCOMB_GAT4 and list of combinations LCOMB_GAT
  !
  !----------------------------------------------------------------------

  if( IPARALL ) then
     if( IMASTER ) then
        ncomb = 0
        call memory_alloca(memor_dom,'NCOMB_GAT4','chkcod',ncomb_gat4,int(PAR_WORLD_SIZE,4),lboun=0_4)
     end if
     ncomb4 = int(ncomb,4)
     call PAR_GATHER(ncomb4,ncomb_gat4)
     if( IMASTER ) then
        ncomb_par  = int(sum(ncomb_gat4),ip)
        ncomb_gat4 = ncomb_gat4*int(mcono,4)
        call memory_alloca(memor_dom,'LCOMB_GAT','chkcod',lcomb_gat,mcono,ncomb_par)
     end if
     call PAR_GATHERV(lcomb,lcomb_gat,ncomb_gat4)
     if( IMASTER ) then
        call memory_alloca(memor_dom,'LCOMB','chkcod',lcomb,mcono,ncomb_par)
        ncomb = 0    
        do icomb_par = 1,ncomb_par
           icomb = 1
           icode = lcomb(1,icomb)
           do while( icode > 0 ) 
              if(          lcomb(1,icomb) == lcomb_gat(1,icomb_par) &
                   & .and. lcomb(2,icomb) == lcomb_gat(2,icomb_par) &
                   & .and. lcomb(3,icomb) == lcomb_gat(3,icomb_par) ) then
                 icode = -1
              else
                 icomb = icomb + 1
                 icode = lcomb(1,icomb)
              end if
           end do
           if( icode /= -1 ) then
              ncomb            = ncomb + 1
              lcomb(1:3,ncomb) = lcomb_gat(1:3,icomb_par)
           end if
        end do
     end if
  end if
  
  !----------------------------------------------------------------------
  !
  ! Output bombinations LCOMB
  !
  !----------------------------------------------------------------------

  if( INOTSLAVE .and. ncomb /= 0 ) then
     call sorti3(ncomb,lcomb)
     call outfor(42_ip,lun_outpu,' ')
     do icomb = 1,ncomb
        ioutp(1) = lcomb(1,icomb)
        ioutp(2) = lcomb(2,icomb)
        ioutp(3) = lcomb(3,icomb)
        call outfor(43_ip,lun_outpu,' ')
     end do
  end if

  call memory_deallo(memor_dom,'NCOMB_GAT4','chkcod',ncomb_gat4)
  call memory_deallo(memor_dom,'LCOMB_GAT' ,'chkcod',lcomb_gat)
  call memory_deallo(memor_dom,'LCOMB'     ,'chkcod',lcomb)

end subroutine chkcod
