!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_matrix()
  !------------------------------------------------------------------------
  !
  ! Compute elemental matrix and RHS of the following problem:
  !   _                             _
  !  |                             |
  !  |  rho*cp/(theta*dt) T*v dw + |  rho*cp*[u.grad(T)] dw
  ! _|W                           _|W
  !
  !    _                         _               _
  !   |                         |               |
  ! + |  k*grad(T).grad(v) dw + |  ar*Tr*v ds = |  Q*v dw
  !  _|W                       _|S             _|
  !
  !   _                               _
  !  |                               |
  !  |  rho*cp/(theta*dt) T^n*v dw - |  (qr-ar*Tr) v ds
  ! _|W                             _|S
  ! 
  ! All terms are calculated at current time step, i.e. n+theta
  !
  !------------------------------------------------------------------------

#include "def_vector_size.inc"
  use def_parame
  use def_master
  use def_temper
  use def_domain
  use def_coupli
  use def_kermod,                only : kfl_waexl_ker, temel_ker, kfl_temper_vect
  use mod_interpolation
  use mod_couplings
  use mod_communications, only : PAR_MIN
  use mod_ker_timeline,   only : ker_timeline
  use mod_messages,       only : messages_live
  use mod_timings,        only : timings_assembly
  use mod_wall_exchange,  only : ker_waexlo_getval

  implicit none
  integer(ip) :: ipoin
  real(rp)    :: time1,time2,time3,time_element,time_boundary
  
  call ker_timeline('INI_ASSEMBLY')
  !
  ! Initializations
  !
  call inisol()
  resgs_tem(1) = 0.0_rp
  resgs_tem(2) = 0.0_rp


  call cputim(time1) 
  if( kfl_waexl_ker == 1_ip ) &
       call ker_waexlo_getval(tempe,temel_ker)

  if( INOTMASTER ) then 

     if( kfl_discr_tem == NODAL_SCHEME ) then
        !
        ! FE: Element assembly
        !

       if(kfl_temper_vect == 1_ip) then ! temper vectorization ON in kermod

          call tem_elmope_all_fast(1_ip)
       else
           call tem_elmope_new(1_ip)
       end if
       !
        ! FE: Boundary assembly
        !
        call tem_bouope(0_ip)

        !
        ! FE: Assmebly coupling terms 
        !
        call tem_partis()

     else if( kfl_discr_tem == CELL_CENTERED_SCHEME ) then
        !
        ! FV: Assembly
        !
        call tem_finite_volume()
     end if

  end if
  call cputim(time2)

  !! AB what? if( INOTMASTER ) then 

  !! AB what?    if( kfl_discr_tem == NODAL_SCHEME ) then
  !! AB what?       !
  !! AB what?       ! FE: Boundary assembly
  !! AB what?       !
  !! AB what?       call tem_bouope()
  !! AB what?    end if
  !! AB what? end if
  call cputim(time3)
  !
  ! Impose non-uniform heat flux (from FIELDS or SPARE mesh)
  !    
  if(      kfl_sourc_tem == SOURCE_TERM_NODAL_FIELD ) then
     do ipoin = 1,npoin_own
        rhsid(ipoin) = rhsid(ipoin) + heat_source(1,ipoin)
     end do
  else if( kfl_sourc_tem == SOURCE_TERM_SPARE_MESH ) then
     do ipoin = 1,npoin
        rhsid(ipoin) = rhsid(ipoin) + heat_source(1,ipoin)
     end do
  end if
  !
  ! Timings
  !
  time_element  = time2 - time1
  time_boundary = time3 - time2
  call timings_assembly(time_element,time_boundary)  

  call ker_timeline('END_ASSEMBLY')

end subroutine tem_matrix

subroutine tem_test()
  use def_master
  use mod_parall
  use def_domain
  use def_coupli
  use mod_couplings
  implicit none

  integer(ip)         :: kcoup,icoup,ndofn,kpoin,ipoin
  real(rp),   pointer :: xx_send(:)
  real(rp),   pointer :: xx_recv(:)

  nullify(xx_send,xx_recv)
  
  if( INOTMASTER ) then
     allocate( xx_send(npoin) )
     allocate( xx_recv(npoin) )
     xx_send = 1.0_rp
     xx_recv = 0.0_rp
     ndofn    = 1
     do kcoup = 1,ncoup_implicit_n
        icoup = lcoup_implicit_n(kcoup)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,xx_recv,xx_send)
        color_target = coupling_type(icoup) % color_target
        if(I_AM_IN_COLOR(color_target)) then
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              write(*,'(a,i1,1x,i3,100(1x,e12.6))')'N: ',kfl_paral,lninv_loc(ipoin),xx_recv(ipoin)
           end do
        end if
     end do    
     deallocate( xx_send )
     deallocate( xx_recv )
  end if
    
  if( INOTMASTER ) then
     allocate( xx_send(npoin) )
     allocate( xx_recv(npoin) )
     xx_send = 1.0_rp
     xx_recv = 0.0_rp
     ndofn    = 1
     do kcoup = 1,ncoup_implicit_d
        icoup = lcoup_implicit_d(kcoup)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,xx_recv,xx_send)
        color_target = coupling_type(icoup) % color_target
        if(I_AM_IN_COLOR(color_target)) then
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              write(*,'(a,i1,1x,i3,100(1x,e12.6))')'D: ',kfl_paral,lninv_loc(ipoin),xx_recv(ipoin)
           end do
        end if
     end do    
     deallocate( xx_send )
     deallocate( xx_recv )
  end if
  
  
end subroutine tem_test
