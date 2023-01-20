!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine opebcs(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/opebcs
  ! NAME
  !    opebcs
  ! DESCRIPTION
  !    Operation on boundaries:
  !    - Geometrical boundary conditions
  !    - Extrapolaiton from boundaries to nodes
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use mod_messages,      only : livinf
  use mod_domain,        only : domain_memory_allocate
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iboun,icond,ipara,ncond,ivalu
  integer(ip)             :: ipoin

  select case( itask )

  case ( 1_ip )

     if( INOTMASTER ) then
        !
        ! Extrapolate from boundary to node
        ! 
        if( kfl_extra == 1 ) then
           if( kfl_icodn == 0 ) then
              kfl_icodn = 1
              call domain_memory_allocate('KFL_CODNO')
              do ipoin = 1,npoin
                 kfl_codno(1,ipoin) = mcodb+1
              end do
           end if
           call extbcs()
        end if

     else

        if( kfl_extra == 1 ) then
           call livinf(0_ip,'EXTRAPOLATE BOUNDARY CODES TO NODE CODES',0_ip)
           kfl_icodn = 1
        end if

     end if

  case ( 2_ip )

     if( INOTMASTER ) then
        !
        ! Initializations and allocations
        !
        if( kfl_geome == 1 ) then 

           !-------------------------------------------------------------
           !
           ! Geometrical boundaries
           !
           ! Define boundary codes: 
           ! Prescribed ........ 10
           ! Freestream ........ 20
           ! Wall law .......... 30
           ! Symmetry .......... 40
           ! No slip / Wall .... 50
           ! Free / Outflow .... 60
           !
           !-------------------------------------------------------------
           !
           ! Allocate memory for KFL_GEOBO: includes fringe boundaries
           !
           call domain_memory_allocate('KFL_GEOBO')
           if( kfl_icodb <= 0 ) then
              call runend('NO BOUNDARY CODE HAS BEEN IMPOSED: CANNOT GENERATE GEOMETRICAL NORMALS')
           end if
           ncond = size(npbcs,KIND=ip)
           do iboun = 1,nboun
              do icond = 1,ncond
                 ivalu = 10*icond
                 do ipara = 1,npbcs(icond)
                    if( kfl_codbo(iboun) == lsbcs(ipara,icond) ) kfl_geobo(iboun) = ivalu
                 end do
              end do
           end do
           call geonor(1_ip)
        end if

     else

        if( kfl_geome == 1 ) then
           call geonor(1_ip)
        end if

     end if

  end select

end subroutine opebcs
