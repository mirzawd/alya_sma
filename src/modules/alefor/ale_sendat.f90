!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_sendat(itask)
  !-----------------------------------------------------------------------
  !****f* master/ale_sendat
  ! NAME
  !    ale_sendat
  ! DESCRIPTION
  !    This routines creates the IB structure
  ! for the moment only itask=0 & 1 are ready
  ! y borre algo de 7 
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,ii,jj

  if( itask == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------
     !
     ! Read in ale_reaphy
     !
     call exchange_add( kfl_rigid_ale )
     do ii = 1,3
        call exchange_add( nstli_ale(ii) )
        call exchange_add( nstro_ale(ii) )
     end do
     call exchange_add( kfl_mvext_ale )
     call exchange_add( kfl_ralei_ale )
     call exchange_add( nstra_ale )
     call exchange_add( nenra_ale )
     call exchange_add( kfl_catfo_ale )             ! Catamaran force
     call exchange_add( kfl_grafo_ale )             ! Gravity force
     call exchange_add( kfl_forca_res )             ! Residual based forces for rigid body motion
     call exchange_add( kfl_sprin_ale )             ! Spring activated for rigid body motion
     call exchange_add( kfl_topmo_ale )             ! Top model activated for rigid body motion
     call exchange_add( kfl_pertu_ale )             ! Use perturbation of generalized forces
     call exchange_add( kfl_genco_ale )             ! Generalized coordinates solver for rigid body
     call exchange_add( ralei_ale )
     !
     ! Spring constants for rigid body motion
     !
     do ii = 1_ip, 3_ip
        call exchange_add( sprin_ale(ii) )
     end do

     do iimbo = 1,nrbod
        
        call exchange_add( rbbou(iimbo) % npoib )
        call exchange_add( rbbou(iimbo) % nboib )
        call exchange_add( rbbou(iimbo) % nrbse )
        do ii = 1,10
           call exchange_add( rbbou(iimbo) % lrbse(ii) )
        end do

        call exchange_add( rbbou(iimbo) % massa )
        call exchange_add( rbbou(iimbo) % densi )
        call exchange_add( rbbou(iimbo) % volum )
        do ii = 1,6
           call exchange_add( rbbou(iimbo) % momin(ii) )
        end do
        do ii = 1,3
           call exchange_add( rbbou(iimbo) % posgr(ii) ) 
        end do
        !
        ! Linear motion
        !
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % posil(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % velol(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % accel(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % force(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % vpfor(ii,jj) )
           end do
        end do
        do ii = 1,3
           call exchange_add( rbbou(iimbo) % vforce(ii) )
        end do
        do ii = 1,3
           call exchange_add( rbbou(iimbo) % pforce(ii) )
        end do
        !
        ! Angular motion
        !
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % posia(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % veloa(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % accea(ii,jj) )
           end do
        end do
        do jj = 1,3
           do ii = 1,3
              call exchange_add( rbbou(iimbo) % rotac(ii,jj) )
           end do
        end do
        do ii = 1,4
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % quate(jj,ii) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % torqu(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,4
              call exchange_add( rbbou(iimbo) % vptor(ii,jj) )
           end do
        end do
        do ii = 1,3
           call exchange_add( rbbou(iimbo) % vtorqu(ii) )
        end do
        do ii = 1,3
           call exchange_add( rbbou(iimbo) % ptorqu(ii) )
        end do
     end do
  end if

end subroutine ale_sendat
