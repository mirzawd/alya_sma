!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_updunk.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Updates.
!> @details Updates.
!> @} 
!-----------------------------------------------------------------------
subroutine ale_updunk(itask)
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use mod_memory, only : memory_size
  use def_kermod, only : kfl_adj_prob
  implicit none
  integer(ip), intent(in) :: itask  !> where the subrutine is called
  integer(ip)             :: ipoin,idime,iimbo,itotn,icomp
  !
  ! For the adjoint case, not to enter here
  !
  if (kfl_adj_prob == 1_ip) return

  !----------------------------------------------------------------------
  !
  ! Alefor main variables
  !
  !----------------------------------------------------------------------

  select case ( itask )

  case ( ITASK_INIUNK )
     !
     ! (:,all) <= (:,1): Initial solution
     !
     do icomp = 2,memory_size(dispm,3_ip)
        do ipoin = 1,npoin
           dispm(:,ipoin,icomp) = dispm(:,ipoin,1)
        end do
     end do
     do icomp = 2,memory_size(coord_ale,3_ip)
        do ipoin = 1,npoin
           coord_ale(:,ipoin,icomp) = coord_ale(:,ipoin,1)
        end do
     end do

  case ( ITASK_BEGSTE )
     !
     ! Initial guess for dispm: d(n,0,*) <-- d(n-1,*,*)
     !
     do ipoin = 1, npoin
        dispm    (:,ipoin,2) = dispm    (:,ipoin,1)
        coord_ale(:,ipoin,2) = coord_ale(:,ipoin,1)
     end do

  case ( ITASK_BEGITE )
     !
     ! Initial guess for dispm to be used in the inner iterations 
     !
     do ipoin = 1, npoin
        do idime = 1, ndime
           unkno( (ipoin-1)*ndime+idime) = 0.0_rp
        end do
     end do

  case ( ITASK_ENDINN )
     ! 
     ! Update displacement, mesh velocity and new mesh coordinate
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !
     do ipoin = 1,npoin        
        do idime = 1,ndime          
           itotn = (ipoin-1_ip) * ndime + idime
           
           if( kfl_fixno_ale(idime,ipoin) == -1 .or. kfl_fixno_ale(idime,ipoin) == 3 ) then
              !
              ! FMALE type
              !
              dispm(idime,ipoin,1)     =  unkno(itotn)
              coord_ale(idime,ipoin,1) =  coord_ale(idime,ipoin,1) + dispm(idime,ipoin,1)
              velom(idime,ipoin)       = -dtinv * ( coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,3) )
              bvess_ale(idime,ipoin,1) =  0.0_rp                 
           else
              !
              ! Normal type
              !
              dispm(idime,ipoin,1)     =  unkno(itotn)
              coord_ale(idime,ipoin,1) =  coord_ale(idime,ipoin,1) + dispm(idime,ipoin,1)
              velom(idime,ipoin)       =  dtinv * ( coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,3) )
              coord(idime,ipoin)       =  coord_ale(idime,ipoin,1)
           end if

        end do
     end do

  case ( ITASK_BEGZON )   
     !
     ! Should check that
     !
     do ipoin = 1, npoin
        coord_ale(:,ipoin,1) = coord_ale(:,ipoin,3) ! CHECK THIS I DO NOT UNDERTSAND
        coord(:,ipoin)       = coord_ale(:,ipoin,1)
     end do

  case ( ITASK_ENDITE ) 
     !
     ! (:,2) <= (:,1): End of inner iteration
     !        
     do ipoin = 1, npoin
        dispm    (:,ipoin,2) = dispm    (:,ipoin,1)
        coord_ale(:,ipoin,2) = coord_ale(:,ipoin,1)
     end do

  case ( ITASK_ENDSTE )
     !
     ! (:,3) <= (:,1): End of time step
     ! (:,4) <= (:,3)
     ! (:,5) <= (:,4)
     ! ...
     !    
     do ipoin = 1,npoin
        dispm    (:,ipoin,3) = dispm    (:,ipoin,1)
        coord_ale(:,ipoin,3) = coord_ale(:,ipoin,1)
     end do

  end select

  if( kfl_rigid_ale == 1 ) then

     !----------------------------------------------------------------------
     !
     ! Rigid body
     !
     !----------------------------------------------------------------------

     select case ( itask )

     case ( ITASK_INIUNK )

        do iimbo = 1,nrbod
           do icomp = 2,4
              rbbou(iimbo) % quate(1,icomp) = rbbou(iimbo) % quate(1,1)
              rbbou(iimbo) % q_dot(1,icomp) = rbbou(iimbo) % q_dot(1,1)
              rbbou(iimbo) % posil(:,icomp) = rbbou(iimbo) % posil(:,1)
              rbbou(iimbo) % velol(:,icomp) = rbbou(iimbo) % velol(:,1)
              rbbou(iimbo) % accel(:,icomp) = rbbou(iimbo) % accel(:,1)
              rbbou(iimbo) % posia(:,icomp) = rbbou(iimbo) % posia(:,1)
              rbbou(iimbo) % veloa(:,icomp) = rbbou(iimbo) % veloa(:,1)
              rbbou(iimbo) % accea(:,icomp) = rbbou(iimbo) % accea(:,1)
              rbbou(iimbo) % quate(:,icomp) = rbbou(iimbo) % quate(:,1)
              rbbou(iimbo) % q_dot(:,icomp) = rbbou(iimbo) % q_dot(:,1)
           end do
        end do

     case ( ITASK_BEGSTE )  
        !
        ! Assign a(n,0,*) <-- a(n-1,*,*),  RB initial guess for outer iterations
        ! For the Force and Moment this might be a good place to extrapolate, force( ,nprev_ale) and torqu( ,nprev_ale)
        ! must have been set in nsi after solving if in the previous step - This supposes ale is solved before
        ! nsi if this is not the case some other strategy must be thought
        !
        if ( kfl_genco_ale == 1_ip ) then
           !
           ! For the moment, only angular coordinates are being used as generalized coordinates
           !
           do iimbo = 1,nrbod
              rbbou(iimbo) % posia(:,2_ip) = rbbou(iimbo) % posia(:,1)
              rbbou(iimbo) % veloa(:,2_ip) = rbbou(iimbo) % veloa(:,1)
              rbbou(iimbo) % vpfor(:,2_ip) = rbbou(iimbo) % vpfor(:,1)
           end do

        else

           do iimbo = 1,nrbod
              rbbou(iimbo) % quate(1,2) = rbbou(iimbo) % quate(1,1)
              rbbou(iimbo) % q_dot(1,2) = rbbou(iimbo) % q_dot(1,1)
              rbbou(iimbo) % accel(:,2) = rbbou(iimbo) % accel(:,1) 
              rbbou(iimbo) % velol(:,2) = rbbou(iimbo) % velol(:,1) 
              rbbou(iimbo) % posil(:,2) = rbbou(iimbo) % posil(:,1)
              rbbou(iimbo) % accea(:,2) = rbbou(iimbo) % accea(:,1) 
              rbbou(iimbo) % veloa(:,2) = rbbou(iimbo) % veloa(:,1) 
              rbbou(iimbo) % posia(:,2) = rbbou(iimbo) % posia(:,1)                 
              rbbou(iimbo) % quate(:,2) = rbbou(iimbo) % quate(:,1)
              rbbou(iimbo) % q_dot(:,2) = rbbou(iimbo) % q_dot(:,1)
              
              if ( ( kfl_foexo_ale == 1_ip ) .or. ( kfl_crist_ale == 1_ip ) ) then ! Force ( & Torque) extrapolation order
                 rbbou(iimbo) % vpfor(:,2) = 1.0_rp * rbbou(iimbo) % vpfor(:,1)
                 rbbou(iimbo) % vptor(:,2) = 1.0_rp * rbbou(iimbo) % vptor(:,1)
              else
                 rbbou(iimbo) % vpfor(:,2) = 1.5_rp * rbbou(iimbo) % vpfor(:,1) - 0.5_rp * rbbou(iimbo) % vpfor(:,nprev_ale+1_ip)
                 rbbou(iimbo) % vptor(:,2) = 1.5_rp * rbbou(iimbo) % vptor(:,1) - 0.5_rp * rbbou(iimbo) % vptor(:,nprev_ale+1_ip)
              end if
           end do

        end if  ! kfl_genco_ale

     case ( ITASK_ENDITE )
        !
        ! keep the last dispm for the next coupling iteration
        !
        if ( kfl_genco_ale == 1_ip ) then
           !
           ! For the moment, only angular coordinates are being used as generalized coordinates
           !
           do iimbo = 1,nrbod
              rbbou(iimbo) % posia(:,2_ip) = rbbou(iimbo) % posia(:,1_ip)
              rbbou(iimbo) % veloa(:,2_ip) = rbbou(iimbo) % veloa(:,1_ip)
              rbbou(iimbo) % force(:,2_ip) = rbbou(iimbo) % force(:,1_ip)
           end do

        else
           
           do iimbo = 1,nrbod
              rbbou(iimbo) % quate(1,2) = rbbou(iimbo) % quate(1,1)
              rbbou(iimbo) % q_dot(1,2) = rbbou(iimbo) % q_dot(1,1)
              rbbou(iimbo) % accel(:,2) = rbbou(iimbo) % accel(:,1) 
              rbbou(iimbo) % velol(:,2) = rbbou(iimbo) % velol(:,1) 
              rbbou(iimbo) % posil(:,2) = rbbou(iimbo) % posil(:,1)
              rbbou(iimbo) % accea(:,2) = rbbou(iimbo) % accea(:,1) 
              rbbou(iimbo) % veloa(:,2) = rbbou(iimbo) % veloa(:,1) 
              rbbou(iimbo) % posia(:,2) = rbbou(iimbo) % posia(:,1)                 
              rbbou(iimbo) % quate(:,2) = rbbou(iimbo) % quate(:,1)
              rbbou(iimbo) % q_dot(:,2) = rbbou(iimbo) % q_dot(:,1)
              rbbou(iimbo) % force(:,2) = rbbou(iimbo) % force(:,1) 
              rbbou(iimbo) % torqu(:,2) = rbbou(iimbo) % torqu(:,1) 
           end do

        end if ! kfl_genco_ale

     case( ITASK_ENDSTE )
        !
        ! Assign a(n,i-1,*) <-- a(n,i,*)
        !
        if ( kfl_genco_ale == 1_ip ) then
           !
           ! For the moment, only angular coordinates are being used as generalized coordinates
           !
           do iimbo = 1,nrbod
              rbbou(iimbo) % posia(:,4_ip) = rbbou(iimbo) % posia(:,3_ip)
              rbbou(iimbo) % posia(:,3_ip) = rbbou(iimbo) % posia(:,1_ip)
              rbbou(iimbo) % veloa(:,4_ip) = rbbou(iimbo) % veloa(:,3_ip)
              rbbou(iimbo) % veloa(:,3_ip) = rbbou(iimbo) % veloa(:,1_ip)
              rbbou(iimbo) % force(:,4_ip) = rbbou(iimbo) % force(:,3_ip)
              rbbou(iimbo) % force(:,3_ip) = rbbou(iimbo) % force(:,1_ip)
           end do

        else

           do iimbo = 1,nrbod
              rbbou(iimbo) % quate(1,4) = rbbou(iimbo) % quate(1,3)
              rbbou(iimbo) % q_dot(1,4) = rbbou(iimbo) % q_dot(1,3)
              rbbou(iimbo) % accel(:,4) = rbbou(iimbo) % accel(:,3) 
              rbbou(iimbo) % velol(:,4) = rbbou(iimbo) % velol(:,3) 
              rbbou(iimbo) % posil(:,4) = rbbou(iimbo) % posil(:,3)
              rbbou(iimbo) % accea(:,4) = rbbou(iimbo) % accea(:,3) 
              rbbou(iimbo) % veloa(:,4) = rbbou(iimbo) % veloa(:,3) 
              rbbou(iimbo) % posia(:,4) = rbbou(iimbo) % posia(:,3)                 
              rbbou(iimbo) % quate(:,4) = rbbou(iimbo) % quate(:,3)
              rbbou(iimbo) % q_dot(:,4) = rbbou(iimbo) % q_dot(:,3)
              rbbou(iimbo) % force(:,4) = rbbou(iimbo) % force(:,3) 
              rbbou(iimbo) % torqu(:,4) = rbbou(iimbo) % torqu(:,3) 
              !
              rbbou(iimbo) % quate(1,3) = rbbou(iimbo) % quate(1,1)
              rbbou(iimbo) % q_dot(1,3) = rbbou(iimbo) % q_dot(1,1)
              rbbou(iimbo) % accel(:,3) = rbbou(iimbo) % accel(:,1) 
              rbbou(iimbo) % velol(:,3) = rbbou(iimbo) % velol(:,1) 
              rbbou(iimbo) % posil(:,3) = rbbou(iimbo) % posil(:,1)
              rbbou(iimbo) % accea(:,3) = rbbou(iimbo) % accea(:,1) 
              rbbou(iimbo) % veloa(:,3) = rbbou(iimbo) % veloa(:,1) 
              rbbou(iimbo) % posia(:,3) = rbbou(iimbo) % posia(:,1)                 
              rbbou(iimbo) % quate(:,3) = rbbou(iimbo) % quate(:,1)
              rbbou(iimbo) % q_dot(:,3) = rbbou(iimbo) % q_dot(:,1)
              
              if ( kfl_crist_ale /=1 ) rbbou(iimbo) % force(:,3) = rbbou(iimbo) % force(:,1)   ! for cristobal's case this is not needed
              if ( kfl_crist_ale /=1 ) rbbou(iimbo) % torqu(:,3) = rbbou(iimbo) % torqu(:,1) 
           end do

        end if  ! kfl_genco_ale

     end select

  end if

end subroutine ale_updunk

