!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_boumat.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Element gather
!> @details Element gather for the boundary elements
!> @}
!-----------------------------------------------------------------------

! subroutine neu_boumat(pnodb,pgaub,pnode,iboun,lboel,gbsha,gbsur,baloc,borad,bvnat,elmat,elrhs,dotpr_retorna)
subroutine neu_boumat(pnodb,pgaub,pnode,iboun,lboel,gbsha,gbsur,baloc,borad,elmat,elrhs,dotpr_retorna)

  use def_kintyp, only :  ip,rp
  use def_parame, only :  pi,Stefan_Boltzmann
  use def_domain, only :  ndime
  use def_neutro
!  , only :  ADR_neu
!  use def_neutro, only :  num_energies_neu
!  use def_neutro, only :  num_directions_neu
!  use def_neutro, only :  current_energy_neu
!  use def_neutro, only :  current_direction_neu
!  use def_neutro, only :  direc_neu
!  use def_neutro, only :  weigd_neu
!  use def_neutro, only :  scattering_neu
!  use def_neutro, only :  nitsche_neu
!  use def_neutro, only :  kfl_smobo_neu,source_bound

  implicit none
  integer(ip), intent(in)  :: pnodb !> Number of nodes in the  boundary element
  integer(ip), intent(in)  :: pgaub !> Number of gauss points in the boundary element
  integer(ip), intent(in)  :: pnode !> Number of nodes in the element (Domain)
  integer(ip), intent(in)  :: iboun !> numer of BC
  integer(ip), intent(in)  :: lboel(pnodb)  ! connectivity matrix in the contour (relates number of nodes with elements)
  real(rp),    intent(in)  :: gbsha(pnodb,pgaub)  !> test/shape function, it is saved in elmar % shape
  real(rp),    intent(in)  :: gbsur(pgaub) !> surface , jacobian times weight of the Gauss point
  real(rp),    intent(in)  :: baloc(ndime,ndime,pgaub)   !> Local directions, normal to boundary element (baloc(1:ndime, ndime)) and
                                                         ! two tangencial (baloc(1:ndime,1) and baloc(1:ndime,2) for 3D 
  real(rp),    intent(in)  :: borad(num_energies_neu,num_directions_neu,pnodb) !> Boundary radiation for each energy, direction and 
                                                                               ! boundary node
!   real(rp),    intent(in)  :: bvnat !> Natural boundary conditions
  real(rp),    intent(out) :: elmat(pnode,pnode) !> Element matrix
  real(rp),    intent(out) :: elrhs(pnode) !> Element right-hand side
  real(rp),    intent(out) :: dotpr_retorna
  integer(ip)              :: idofn,jdofn,jdir,nout !> Indexes
  integer(ip)              :: outgo(num_directions_neu) !> Outward direction
  integer(ip)              :: igaub,inodb,inode !> Indexes for boundary gauss point, boundary node and element node
!   real(rp)                 :: btemp !,bolt4!> Used to be the radiation at the boundary (boundary temperature)
  real(rp)                 :: ndots(num_directions_neu) !> producto escalar de direccion con normal hacia fuera (n dot s)
  real(rp)                 :: tract,right(pnode)!> variables intermedias para computar condiciones de contorno
  real(rp)                 :: work1,dotpr,mmwat !> variables intermedias
  real(rp)                 :: gprad(num_energies_neu,num_directions_neu)!> Radiacion interpolada a los puntos de gauss.
  
  integer(ip)              :: refle,influx, isour

  real(rp)                 :: emiss,bsour,albedo2

!   integer(ip)              :: idime,ii,num_real(num_directions_neu),direc_reflex(ndime,num_directions_neu)
!   real(rp)                 :: distancia, disaux


  if(kfl_fixbo_neu(current_energy_neu,current_direction_neu) % l(iboun) ==3 .or. &
     kfl_fixbo_neu(current_energy_neu,current_direction_neu) % l(iboun)==2) then
  
     refle = kfl_fixbo_neu(current_energy_neu,current_direction_neu) % l(iboun)
     albedo2 = albedo_neu
  else

     albedo2 = 0.0_rp
     refle=0
  endif


  influx=kfl_fixbo_neu(current_energy_neu,current_direction_neu) % l(iboun)


  emiss = 1.0_rp
  bsour = 0.0_rp
  if(influx==4) then
     isour = kfl_funbo_neu(current_energy_neu,current_direction_neu) % l(iboun)
     bsour = source_bound_neu(isour,current_energy_neu) !5.67E+4 !100000.0_rp
  endif
  
  !
  ! Exiting direction, no boundary condition
  !
  igaub = 1
  dotpr = dot_product(direc_neu(1:ndime,current_direction_neu),baloc(1:ndime,ndime,igaub))
  dotpr_retorna=dotpr
  if( dotpr > 0.0_rp ) return  ! if outgoing direction then return
  
  ndots = 0.0_rp
  ndots(current_direction_neu) = dotpr_retorna
  tract = -ndots(current_direction_neu) * nitsche_neu  * bsour/ pi ! * emiss 

  gauss: do igaub = 1,pgaub  ! impose b.c.

   !   ndots(current_direction_neu) = dotpr_retorna
     !btemp        = btemp_neu(iboun, igaub)
     !bolt4        = btemp*btemp*btemp*btemp*boltz   ! For us this is bvnat
     !bolt4        = bvnat   ! For us this is bvnat
     !bolt4 = (tempe**4.0_rp)*Stefan_Boltzmann !>We no longer need to calculate the radiation
     !
     ! if nodal matrices, add contibution of nitsche's smooth term to matrix
     !
     do inodb = 1,pnodb
        inode              = lboel(inodb) !> We find to which node the boundary point corresponds to
        mmwat              = -ndots(current_direction_neu) * nitsche_neu * gbsha(inodb,igaub) !> Gauss matrix contribution
        elmat(inode,inode) = elmat(inode,inode) + mmwat * gbsur(igaub)!> We add the contribution to the elemental matrix
     end do
     !
     ! rhs  emmission contribution: u_in
     !
     !tract = -ndots(current_direction_neu) * nitsche_neu * bvnat   ! emiss*bolt4/pi
   !   tract = -ndots(current_direction_neu) * nitsche_neu  * bsour/ pi ! * emiss 
     right =  0.0_rp
     !
     !   If reflection: rho/pi*(...) : DE MOMENTO = ZERO
     !
     if( refle == 3 .or. refle==2 ) then
        !
        ! store exiting directions : this ro/2*pi * int_beta wb*ub* (s.n) (beta is outflow)
        ! outgo = list of directions that go out
        !
        outgo = 0
        nout = 0
        do idofn = 1,num_directions_neu
           dotpr = dot_product(direc_neu(1:ndime,idofn),baloc(1:ndime,ndime,igaub))
           if( dotpr > 0.0_rp ) then    ! exiting directions
              ndots(idofn) = dotpr      ! n·s (s:direction of propagation)
              nout          = nout + 1    ! This stores the number of directions pointing outwards
              outgo(nout)   = idofn       ! store exiting directions
           end if
           
         !   direc_reflex(1:ndime,idofn) = 0.0_rp
         !   do idime=1,ndime
         !      direc_reflex(idime,idofn) = direc_neu(idime,idofn)-2.0_rp*baloc(idime,ndime,igaub)*dotpr
         !   end do
         !   distancia=1000.0_rp
         !   do ii=1,num_directions_neu
         !      disaux=0.0
         !      do idime=1,ndime
         !         disaux = disaux + (direc_reflex(idime,idofn)-direc_neu(idime,ii))**2
         !      end do
         !      disaux = sqrt(disaux)
         !      if(disaux<distancia) then
         !         distancia=disaux
         !         num_real(idofn)=ii
         !      end if
         !   end do


        end do

        if( kfl_smobo_neu == 0 ) then   ! =1 always smooth
           !
           ! Nitsche's boundary conditions (reflective terms)
           !
           gprad = 0.0_rp
           work1 = 0.0_rp
           do jdir = 1,nout
              jdofn = outgo(jdir)
              do inodb = 1,pnodb
                 gprad(current_energy_neu,jdofn) = gprad(current_energy_neu,jdofn) + &
                                                      borad(current_energy_neu,jdofn,inodb) * gbsha(inodb,igaub)
              end do
              work1  = work1 + albedo2 * weigd_neu(jdofn) * ndots(jdofn) * gprad(current_energy_neu,jdofn)
           end do
           tract = tract - work1  / pi * ndots(current_direction_neu) * nitsche_neu

        else
           !
           ! Smooth nitsche's boundary conditions
           !
           do inodb = 1,pnodb
             
              do jdir = 1,nout
                 jdofn        = outgo(jdir)
                 !> integral n.s u ds we do not multiply by gpsur or gpsha because we integrate over the directions, and 
                 !> this integral is already done with the weights of the directions
                 right(inodb) = right(inodb) + albedo2 * weigd_neu(jdofn) * ndots(jdofn) * borad(current_energy_neu,jdofn,inodb)
              end do
              right(inodb) = - right(inodb) * ndots(current_direction_neu) * nitsche_neu/ pi   
              ! integral por reflectividad/ pi nitsche *ndots z
           end do
        end if
     end if
     !
     ! Assembly
     !
     do inodb = 1,pnodb
        inode = lboel(inodb)
       ! elemental rhs  vector  due to emissivity (tract) and reflectivity (right)
        elrhs(inode) = elrhs(inode) + gbsha(inodb,igaub) * (tract+right(inodb)) * gbsur(igaub) 
     end do

  end do gauss

!stop
end subroutine neu_boumat
