!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module  mod_ker_polynomial
  use def_kintyp, only :  ip,rp,r1p
  use def_kermod, only :  typ_valpr_ker
  use def_master, only :  INOTMASTER
  use def_kermod, only : number_space_time_function,space_time_function
  use mod_ker_space_time_function  
  implicit none

  interface
    function func_template( T, params ) result(res) 
      use def_parame, only: rp
      implicit none
      real(rp), intent(in) :: T, params(:) 
      real(rp)             :: res  
    end function  
  end interface

  logical(ip)  ::               debug    = .false.
  character(5) :: ker_polynomial_name    = 'POLYN'

  private 
  public :: ker_polynomial_name
  public :: ker_polynomial_readat 
  public :: ker_polynomial_allaws
  public :: ker_polynomial_proper
 

contains 

  subroutine ker_polynomial_allaws(prope_ker, ilaws, id_module, iwhere, id_gradi)
  !
  !!mod_ker_polynomial, only : ker_polynomial_allaws, ker_polynomial_name  
  ! ker_allaws.f90 +286 
  !   call ker_polynomial_allaws(condu_ker, 9, ID_TEMPER, 'IELEM', 0)
  !
  implicit none
  type(typ_valpr_ker)  :: prope_ker
  integer(ip)          :: ilaws, id_module, id_gradi 
  character(5)         :: iwhere

  prope_ker % llaws(ilaws) % wname     =  ker_polynomial_name
  prope_ker % llaws(ilaws) % lresp(1)  =  id_module
  prope_ker % llaws(ilaws) % where     =  iwhere       ! -> 'IPOIN'   -> prope_ker % value_ipoin 
                                                       ! -> 'IELEM'   -> prope_ker % value_ielem   
                                                       ! -> 'CONST'   -> prope_ker % value_const 
  prope_ker % llaws(ilaws) % kfl_gradi =  id_gradi     ! -> kfl_grele -> prope_ker % grval_ielem 
  prope_ker % llaws(ilaws) % kfl_deriv =  0            ! -> kfl_drele -> prope_ker % drval_ielem 
  prope_ker % llaws(ilaws) % kfl_grder =  0            ! -> kfl_gdele -> prope_ker % gdval_ielem 

  if(debug) print *, "[ker_polynomial_allwas]", ker_n_functions
  end subroutine


  subroutine ker_polynomial_readat(prope_ker, imate)
  !
  !!mod_ker_polynomial, only : ker_polynomial_readat, ker_polynomial_name  
  ! ker_readat.f90 +347 
  ! else if (exists(ker_polynomial_name)) then
  !   call ker_polynomial_readat( condu_ker, imate )
  !
  implicit none
  type(typ_valpr_ker)  :: prope_ker
  integer(ip)          :: imate

  prope_ker % wlaws(imate) = ker_polynomial_name  

  if(debug) print *, "[ker_polynomial_readat]", ker_n_functions  
  end subroutine 


  subroutine ker_polynomial_proper(prope_ker, imate, which_time) 
  !
  !!mod_ker_polynomial, only : ker_polynomial_proper, ker_polynomial_name  
  ! mod_ker_proper.f90 +3257 
  ! else if( prope_ker % wlaws(imate) == ker_polynomial_name ) then
  !   call ker_polynomial_proper(prope_ker, imate, which_time)
  ! 
  implicit none
  type(typ_valpr_ker)  :: prope_ker
  integer(ip)          :: imate, which_time

  call ker_polynomial_assemble(prope_ker, imate, which_time, func02)  

  if(debug) print *, "[ker_polynomial_proper]", ker_n_functions, int( prope_ker%rlaws(1,imate) ), prope_ker%rlaws(2,imate)  
  end subroutine


  subroutine ker_polynomial_assemble(prope_ker, imate, which_time, FunctionT)
    use def_domain
    use def_master
    use def_kermod
    implicit none
    type(typ_valpr_ker)         :: prope_ker
    integer(ip), intent(in)     :: imate
    integer(ip), intent(in)     :: which_time
    procedure(func_template)    :: FunctionT    

    !
    integer(ip)          :: pnode,ielem,ipoin, igaus
    integer(ip)          :: pgaus,inode,pelty,inodb,pblty
    integer(ip)          :: pnodb,iboun,pgaub,igaub, kpoin 
    real(rp)                     :: gbtem
    real(rp)                     :: gpcar(ndime,mnode,mgaus)
    real(rp)                     :: gpvol(mgaus)
    real(rp)                     :: dummr(ndime*mgaus)
    real(rp)                     :: elcod(ndime,mnode)
    real(rp)                     :: eltem(mnode)
    real(rp)                     :: gp_temperature(mgaus)
    real(rp)                     :: params(2), T0, T 
    integer(ip)          :: ilaws  

    ilaws = prope_ker % ilaws(imate) 
    if( prope_ker % llaws(ilaws) % where == 'IELEM') then 
       !call runend("ERROR:  where /= 'IELEM' !!") 
       ! 
       if(IMASTER) return
       ! 
       params(1) = prope_ker % rlaws(1,imate)
       params(2) = prope_ker % rlaws(3,imate)
       T0        = prope_ker % rlaws(2,imate)

       do ielem = 1,nelem
          !
          if( lmate(ielem) == imate ) then
             pelty = ltype(ielem)
             if( pelty > 0 ) then
                pgaus = ngaus(pelty)
                pnode = lnnod(ielem)
                ! 
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   eltem(        inode) = tempe(ipoin,which_time)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do
                ! 
                call elmcar(&
                     pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                     elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                     dummr,ielem)
                ! 
                gp_temperature=0.0_rp
                do igaus = 1,pgaus
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      gp_temperature(igaus) = gp_temperature(igaus) + eltem(inode) * elmar(pelty) % shape(inode,igaus)
                   end do
                   if (gp_temperature(igaus) <= epsilon(1.0_rp) ) then
                      ! For initialization we default to reference temperature
                      gp_temperature(igaus) = T0
                   endif
                   prope_ker % value_ielem(ielem) % a(igaus) = FunctionT( gp_temperature(igaus), params) 
                end do
                ! 
             end if
          end if
          !
       end do
       prope_ker % kfl_nedsm = 1
       ! 
       do iboun = 1,nboun
          pblty = abs(ltypb(iboun))
          pnodb = nnode(pblty)
          ielem = lelbo(iboun)
          ! 
          if( lmate(ielem) == imate ) then
             pelty = ltype(ielem)
             if( pelty > 0 ) then
                pgaub = ngaus(pblty)
                do inodb = 1,pnodb
                   eltem(inodb) = tempe(lnodb(inodb,iboun),1)
                end do
                ! 
                do igaub = 1,pgaub
                   gbtem = 0.0_rp
                   do inodb = 1,pnodb
                      gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                   end do
                   if (gbtem <= epsilon(1.0_rp) ) then
                      ! For initialization we default to reference temperature
                      gbtem = T0 
                   endif
                   prope_ker % value_iboun(iboun) % a(igaub) = FunctionT( gbtem, params ) 
                end do
                ! 
             end if
          end if
          ! 
       end do
       ! 
    else & 
         if( prope_ker % llaws(ilaws) % where == 'IPOIN') then
       !
       call runend("ERROR: 'IPOIN'-> 'IELEM', ESTO FALLA Y NO TENGO IDEA PORQUE...")
       ! 
       do kpoin = 1,nmatn(imate)
          ipoin = lmatn(imate) % l(kpoin)
          T     = tempe(ipoin,which_time)
          prope_ker % value_ipoin(ipoin) = FunctionT( T, params )  
       end do
       ! 
    endif
    ! 
  end subroutine ker_polynomial_assemble


  function func01(T, params) result(kappa) 
  implicit none
  real(rp), intent(in) :: T, params(:) 
  real(rp)             :: B1, S1, kappa 
  B1 = params(1) ! B1 = 2.6438e-3 J/s/m/K**0.5 
  S1 = params(2) ! S1 = 254.4 K
  !
  kappa = B1 * T**1.5_rp / ( T + S1*(1e-12_rp)**(1.0_rp/T) )  
 !print *, B1, S1, T, kappa 
  !
  end function 


  function func02(t_in, params) result(y_out)
  implicit none
  real(rp), intent(in) :: t_in, params(:) 
  real(rp)             :: y_out
  real(rp)             :: x_in=0.0_rp, y_in=0.0_rp, z_in=0.0_rp  
  integer(ip) :: ifun 
    !  
    ifun = int( params(1) ) 
    call ker_space_time_function(ifun, x_in, y_in, z_in, t_in, y_out)
    !
  end function


end module 
