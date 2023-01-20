!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_detect_ourflows.f90
!> @date    26/01/2018
!> @author  Guillaume Houzeaux
!> @brief   Detect outflows 
!> @details Detect outflows with only one single free node
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_detect_outflows()

  use def_master
  use def_domain
  use def_nastin
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_SUM
  use mod_messages,       only : messages_live
  use mod_std
  implicit none
  
  integer(ip)          :: ipoin,izdom,jpoin,jbopo
  integer(ip), pointer :: node_count(:,:)

  nullify(node_count)

  call memory_alloca(mem_modul(1:2,modul),'NODE_COUNT','nsi_detect_outflows',node_count,2_ip,npoin)

  do ipoin = 1,npoin

     if( lpoty(ipoin) > 0 .and. maxval(kfl_fixno_nsi(:,ipoin)) <= 0 ) then
        
        do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
           jpoin = c_dom(izdom)
           if( lpoty(jpoin) > 0 .and. ipoin /= jpoin ) then
              node_count(1,ipoin) = node_count(1,ipoin) + 1
              jbopo = lpoty(jpoin)
              if(    count(kfl_fixno_nsi(:,jpoin) >= 1 ,KIND=ip) == ndime .or. &
                   & ( kfl_fixno_nsi(1,jpoin) >= 1 .and. kfl_fixrs_nsi(jpoin) /= 0 ) ) then
                 node_count(2,ipoin) = node_count(2,ipoin) + 1
              end if
           end if
        end do
       
     end if
  end do
  
  call PAR_INTERFACE_NODE_EXCHANGE(node_count,'SUM')

  jpoin = 0
  do ipoin = 1,npoin_own
     if( node_count(1,ipoin) == node_count(2,ipoin) .and. node_count(1,ipoin) /= 0 ) then
        call messages_live('!!!WARNING: Found two outflows connected by the node'//intost(ipoin))
        jpoin = jpoin + 1
     end if
  end do
  call PAR_SUM(jpoin)
  
  if( jpoin > 0 ) then
     call messages_live(&
          trim(intost(jpoin))//' OUTFLOWS WITH ONLY ONE NODES HAVE BEEN DETECT IN '&
          //trim(namod(modul)),'WARNING') 
  end if
  
  call memory_deallo(mem_modul(1:2,modul),'NODE_COUNT','nsi_detect_outflows',node_count)     

end subroutine nsi_detect_outflows
