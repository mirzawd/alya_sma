!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain  
!> @{
!> @file    open_close.f90
!> @author  Herbert Owen
!> @date    30/11/2015
!> @brief   Indicates if nodes must be opened or closed deppending on the angle between a vector that can be a field or the velocity
!> @details It uses the some ideas from geonor
!!          \verbatim
!!
!!          - Initially it looked at the boundaries with KFL_GEOBO(1:NBOUN): 20 freestream
!!          - Now, I will avoid any use of things from geometrical. In order to do so in the loop over iboun
!!            I will take(and do calculations) a boundary if all its nodes have abs(kfl_fixno_nsi(1,ipoin)) == 7
!!            In this way I will find the normal to the boundary without needing any codes on boundaries.
!!
!!          - In order to decide whether to open or close for nsi we use bvess_nsi (itask=1)
!!            For the other problems veloc(1:ndime,ipoin,2). 
!!      
!!
!!          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine open_close(itask,kfl_fixno,bvess,ndprb)
  use def_parame
  use def_master, only         : INOTMASTER,gisca,veloc,npoin_type
!  use def_master, only         : kfl_paral
  use def_domain, only         : nboun,MEMOR_DOM,nboun,lnodb,ltypb,lelbo,nnode,tolan,lpoty,ndime,npoin
  use mod_memory, only         : memory_alloca,memory_deallo
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE

  implicit none
  real(rp),intent(in)     :: bvess(ndime,npoin,2) ! beware bvess_nsi is used only if itask==1  else veloc(,2)
  integer(ip),intent(in)  :: itask
  integer(ip),intent(in)  :: ndprb
  integer(ip),intent(in)  :: kfl_fixno(ndprb,npoin)

  real(rp),    pointer    :: bouno(:,:)

  integer(ip)             :: inodb,ipoin,nnodb,kount,ninve, idime
  integer(ip)             :: iboun
  real(rp)                :: vnorm,xangl(3)
  nullify(bouno)
 
  if( INOTMASTER ) then    
!     call bougra()   Perhaps I may need it for bounor - Initial check says I don't
     !
     ! Allocate temporal memory pointer bouno
     !
     call memory_alloca(memor_dom,'BOUNO','open_close',bouno,ndime,max(1_ip,nboun))

     !----------------------------------------------------------------------
     !
     ! Compute the normals to the boundaries
     ! I have eliminated the control warning because it has already been done in geonor
     ! Doing it here involved calling to parari(sum) and that is a problem becuase the master does not enter open_close
     !
     !----------------------------------------------------------------------
     ninve=0
     if( nboun > 0 ) then         ! guillaume says we could save it instead of recalculating
        call bounor(nboun,lnodb,ltypb,lelbo,ninve,bouno) ! outputs the normal vector bouno for each boundary
     end if
     !
     gisca = 0_ip !initialization
     do iboun = 1,nboun
        nnodb = nnode(ltypb(iboun))
        kount = 0
        do inodb = 1,nnodb
           ipoin = lnodb(inodb,iboun)
           loop_idime: do idime =1, ndprb
              if( abs(kfl_fixno(idime,ipoin)) == 7) then 
                 kount = kount + 1 
                 exit loop_idime 
              end if
           end do loop_idime
!           if( abs(kfl_fixno(1,ipoin)) == 7) kount = kount + 1
        end do
        if( kount == nnodb ) then  ! all nodes in the boundary element have a fixno == 7 
           do inodb = 1,nnodb
              ipoin = lnodb(inodb,iboun) 
              if( ipoin <= npoin .and.lpoty(ipoin)  > 0 ) then
                 if (itask==1)then
                    !                 do idime=1,ndime
                    !                    xangl(idime) = xfiel(kfl_tran7) % a(idime+(k_transient_bcs-1)*ndime,ipoin) * x_transient_bcs + &
                    !                         xfiel(kfl_tran7) % a(idime+(k_transient_bcs)*ndime,ipoin) * (1.0_rp-x_transient_bcs)
                    !                   xangl(1:ndime) =bvess_nsi(idime,ipoin,1)
                    !                 end do
                    xangl(1:ndime) = bvess(1:ndime,ipoin,1)
                 else ! itask /=1
                    xangl(1:ndime) = veloc(1:ndime,ipoin,2)
                 end if
                 
                 call vecuni(ndime,xangl,vnorm) ! creates unitary vector (its norm stored in vnorm)
                 vnorm        = dot_product(bouno(1:ndime,iboun),xangl(1:ndime)) 
                 
                 if( vnorm <= epsilon(1.0_rp) + cos(0.5_rp*pi-tolan) ) & ! inflow
                      gisca(ipoin) = 1   ! inflow
              end if
           end do
        end if
     end do
     !     call parari('SLX',NPOIN_TYPE,npoin,gisca)   ! now i ma using PAR_INTERFACE_NODE_EXCHANGE  .. see below.
     ! si uso slx suma los valores y entonces en las fronteras le queda gisca=2   por lo cual la condicion a usar en xxx_updbcs es gisca(ipoin) >0
     ! SLX has the advantage that only the slaves enter - for parari(MAX  the master also enters.
     
     ! guillaume sugested  call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')  .. supongo anda tambien con gisca
     ! looking inside the subroutine it says  ! Exchange between slaves. the master does not enter. I do not know the diference between in my code & in my zone??
     ! I leave PAR_INTERFACE_NODE_EXCHANGE because it is newer and it is what guillaume suggested but it should also work with parari SLX
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
!     do ipoin=1,npoin
!           write(kfl_paral+2200,*)'ipoin,gisca(ipoin)',ipoin,gisca(ipoin)
!     end do

     call memory_deallo(memor_dom,'BOUNO','open_close',bouno)

  end if

end subroutine open_close
