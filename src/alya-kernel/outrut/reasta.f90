!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine reasta(&
     kfl_ellen,kfl_stabi,kfl_taust,kfl_sgsno,kfl_sgsti,&
     kfl_shock,staco,shock)
  !-----------------------------------------------------------------------
  !****f* outrut/reasta
  ! NAME 
  !    reasta
  ! DESCRIPTION
  !    This routine reads thge natural length calculation
  ! INPUT
  ! OUTPUT
  ! USES
  !    ecoute
  ! USED BY
  !    ***_reanut
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip),  intent(out) :: kfl_ellen
  integer(ip),  intent(out) :: kfl_stabi
  integer(ip),  intent(out) :: kfl_taust
  integer(ip),  intent(out) :: kfl_sgsno
  integer(ip),  intent(out) :: kfl_sgsti
  integer(ip),  intent(out) :: kfl_shock
  real(rp),     intent(out) :: staco(*)
  real(rp),     intent(out) :: shock  
  integer(ip)               :: istab

  staco(1) = 1.0_rp ! Diffusive term
  staco(2) = 1.0_rp ! Convective term
  staco(3) = 1.0_rp ! Reactive term
  
  do while( words(1) /= 'ENDST') 

     if( words(1)=='PARAM' ) then
        do istab = 1,3
           staco(istab) = param(istab)
        end do

     else if( words(1) == 'TAUST' ) then
        call reatau(kfl_taust)

     else if( words(1) == 'STRAT' ) then
        if( words(2) == 'SU   ' .or. words(2) == 'FIRST' ) then
           kfl_stabi = -1
        else if( words(2) == 'ASGS ' ) then
           kfl_stabi =  0
        else if( words(2) == 'FULLO' ) then 
           kfl_stabi =  1
        else if( words(2) == 'OSS  ' ) then
           kfl_stabi =  2
           if( exists('LIMIT') ) kfl_stabi = 3
        end if

     else if( words(1) == 'ELEME' ) then
        call realen(kfl_ellen)

     else if( words(1) == 'TRACK' ) then 
        if(exists('NONLI')) then
           kfl_sgsno = 1 
        end if
        if(exists('TIME ')) kfl_sgsti = 1 

     else if( words(1) == 'SHOCK' ) then
        if( exists('ISOTR') .or. exists('ON   ') ) then
           kfl_shock = 1
           shock     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
        else if( exists('ANISO') ) then
           kfl_shock = 2
           shock     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
        end if
        
     end if

     call ecoute('REASTA')

  end do

end subroutine reasta
