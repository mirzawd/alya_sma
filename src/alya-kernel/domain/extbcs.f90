!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine extbcs()
  !-----------------------------------------------------------------------
  !****f* Domain/extbcs
  ! NAME
  !    extbcs
  ! DESCRIPTION
  !    Extrapolate cvonditions on boundaries to conditions on nodes
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_elmtyp
  use mod_communications, only : PAR_GHOST_BOUNDARY_EXCHANGE
  use mod_communications, only : PAR_MAX
  use mod_memory,         only : memory_size
  use mod_outfor,         only : outfor
  implicit none
  integer(ip)   :: ipoin,iboun,inodb,kodbo
  integer(ip)   :: kcode,icono,mcod1,isize
  character(20) :: mess1,mess2
  !
  ! Check if we have to extrapolate something
  !
  isize = memory_size(kfl_codbo)
  call PAR_MAX(isize,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
  if( isize == 0 ) return
  
  mess1 = intost(mcono)
  !
  ! Parall: Exchange KFL_CODBO on fringe nodes
  !  
  if( ISLAVE ) then
     call memgen(1_ip,nboun_2,0_ip)
     do iboun = 1,nboun
        gisca(iboun) = kfl_codbo(iboun)
     end do  
     call PAR_GHOST_BOUNDARY_EXCHANGE(gisca,'SUBSTITUTE','IN MY CODE')
  else 
     gisca => kfl_codbo
  end if
  !
  ! Extrapolate to a dummi array GIVEC
  !
  call memgen(1_ip,mcono,npoin_2)
  mcod1 = mcodb+1
  do ipoin = 1,npoin_2
     do icono = 1,mcono
        givec(icono,ipoin) = mcod1
     end do
  end do

  do iboun = 1,nboun_2
     if( lboch(iboun) /= BOEXT ) then
        kodbo = gisca(iboun)
        if( abs(kodbo) <= mcodb ) then
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              !
              ! Look for a place
              !
              kcode = 1
              mcono_loop2: do while( givec(kcode,ipoin) /= mcod1 )
                 if( kodbo == givec(kcode,ipoin) ) kcode = mcono + 1           
                 kcode = kcode + 1
                 if( kcode > mcono ) exit mcono_loop2
              end do mcono_loop2

              if( kcode <= mcono ) then
                 if( givec(kcode,ipoin) /= mcod1 ) then
                    mess2=intost(ipoin)
                    call outfor(2_ip,lun_outpu,'NODE '//trim(mess1)&
                         //' HAVE MORE THAN '//trim(mess2))
                 else
                    givec(kcode,ipoin) = kodbo
                 end if
              end if

           end do

        end if
     end if

  end do
  !
  ! Deallocate memory
  !
  if( ISLAVE ) then
     call memgen(3_ip,nboun_2,0_ip)
  else
     nullify(gisca)
  end if
  !
  ! Sort codes
  !
  do ipoin = 1,npoin
     kcode = 1
     kcode_loop2: do while( givec(kcode,ipoin) /= mcod1 )        
        kcode = kcode + 1
        if( kcode > mcono ) exit kcode_loop2
     end do kcode_loop2
     kcode = kcode - 1
     call heapsorti1(2_ip,kcode,givec(1,ipoin))
  end do
  !
  ! Write over KFL_CODNO
  ! 
  do ipoin = 1,npoin
     if( kfl_codno(1,ipoin) == mcod1 ) then
        do icono = 1,mcono
           kfl_codno(icono,ipoin) = givec(icono,ipoin)
        end do
     end if
  end do
  call memgen(3_ip,mcono,npoin_2)
  
end subroutine extbcs

