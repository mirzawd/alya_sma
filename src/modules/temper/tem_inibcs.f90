!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_inibcs()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_inibcs
  ! NAME
  !    tem_inibcs
  ! DESCRIPTION
  !    This routine applied boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_temper
  use mod_opebcs
  use mod_memory, only : memory_size
  use mod_communications, only : PAR_MAX
  implicit none
  integer(ip)  :: ipoin,pnodb,iboun,inodb,icomp
  integer(ip)  :: pblty

  !-------------------------------------------------------------
  !
  ! Allocate memory
  !
  !-------------------------------------------------------------
  
  call tem_membcs(1_ip)
  
  if( INOTMASTER .and. INOTEMPTY ) then
    
     do ipoin = 1,npoin
        kfl_fixno_tem(1,ipoin) = -1
     end do
     if( kfl_conbc_tem == 0 ) then
        call tem_membcs(2_ip)
     end if

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        if(kfl_conbc_tem==0) then
           iffun     =  1
           kfl_funno => kfl_funno_tem
           kfl_funtn => kfl_funtn_tem 
        else
           iffun      =  0
        end if
        ifloc     =  0
        ifbop     =  0
        kfl_fixno => kfl_fixno_tem
        bvess     => bvess_tem(:,:,1)
        tncod     => tncod_tem
        call reacod(10_ip)
        
     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        !! JMZA. 2018MAR25 
        if(kfl_conbc_tem==0.and..False.) then
           iffun     =  1
           kfl_funbo => kfl_funbo_tem
        else
           iffun      =  0
        end if 

        kfl_fixbo => kfl_fixbo_tem
        bvnat     => bvnat_tem(:,:,1)
        tbcod     => tbcod_tem
        call reacod(20_ip)
        
     end if

     !-------------------------------------------------------------
     !
     ! Put wall condition if delta_tem is negative
     !
     !-------------------------------------------------------------

     if(delta_dom>zetem.or.delta_dom<-zetem) then
        delta_tem=delta_dom
     end if
     if((delta_tem<-zetem).or.(kfl_delta==-1_ip)) then
        do iboun=1,nboun
           if(kfl_fixbo_tem(iboun)==3) then
              pblty=ltypb(iboun)
              pnodb=nnode(pblty)        
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 if(kfl_fixno_tem(1,ipoin)==3) then
                    kfl_fixno_tem(1,ipoin)  = 1
                 end if
              end do
           end if
        end do
     end if
     delta_tem = max(delta_tem,0.0_rp)
     if(kfl_delta==-1_ip) delta_tem = 0.0_rp
     
     call tem_bcntoe()

     !-------------------------------------------------------------
     !
     ! Put -1 value to 0
     !
     !-------------------------------------------------------------

     do ipoin=1,npoin
        if(kfl_fixno_tem(1,ipoin)==-1) kfl_fixno_tem(1,ipoin)=0
     end do

     !-------------------------------------------------------------
     !
     ! Non-constant BC are saved in second line of arrays
     !
     !-------------------------------------------------------------
     
     if( kfl_conbc_tem == 0 ) then
        do ipoin = 1,npoin
           bvess_tem(1,ipoin,2) = bvess_tem(1,ipoin,1)
        end do
        do iboun = 1,nboun
           bvnat_tem(:,iboun,2) = bvnat_tem(:,iboun,1)
        end do
     end if

     !------------------------------------------------------------------
     !
     ! Non-constant boundary conditions
     !
     !------------------------------------------------------------------

     if( kfl_regim_tem == 4 ) then 
        do icomp = 1,min(memory_size(bvess_tem,3_ip),memory_size(bvtem_tem,3_ip))
           do ipoin = 1,npoin
              bvtem_tem(1,ipoin,icomp) = bvess_tem(1,ipoin,icomp)
           end do
        end do
     end if
     
  end if


  if( kfl_conbc_tem == 0 ) then
     !
     ! exists fixity 7 ? - outside if master because i need master to enter for parari max
     !
     kfl_exist_fixi7_tem = 0
     do ipoin = 1,npoin
        if ( kfl_fixno_tem(1,ipoin) == 7 ) kfl_exist_fixi7_tem = 1
     end do

     call PAR_MAX(kfl_exist_fixi7_tem)

  end if

end subroutine tem_inibcs
