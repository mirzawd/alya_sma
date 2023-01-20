!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_cross_wall.f90
!> @author  Guillaume Houzeaux
!> @date    05/11/2015
!> @brief   Check if a particle has crossed the wall
!> @details Check if a particle has crossed the wall
!>                                                     
!>                                                      |// wall
!>           lagrtyp(ilagr) % coord           coord_kp1 |//
!>                      o--------------------------->o  |//
!>                                 dista                |//
!>                                                      |//
!>           lagrtyp(ilagr) % coord                     |//  p2
!>                      o-------------------------------*---->o
!>                      <....................................>
!>                            dista + d/2 + hleng*toler |//
!>                                                      |//
!>                 p1                                   |//  p2
!>                  o---|-------------------------------*---->o
!>                  <...>
!>                d/2 + hleng*toler                     |//
!>                                                      |//
!>           * is the intersection point
!
!> @} 
!-----------------------------------------------------------------------

subroutine pts_cross_wall(&
     ielem,pnode,ilagr,diame,hleng,toler,shapf,deriv,iwall,&
     iboun,coord_kp1,veloc_kp1,accel_kp1,xinte,t,dt_k)

  use def_kintyp, only : ip,rp
  use def_master, only : zeror 
  use def_domain, only : coord,lnodb,ndime
  use def_domain, only : mnodb,lnnob,lnods 
  use def_domain, only : ltypb,nbset,lbset
  use def_domain, only : walld,mnode
  use def_domain, only : nboun
  use mod_elmgeo, only : elmgeo_segfac
  use mod_elmgeo, only : elmgeo_projection_on_a_face
  use mod_elmgeo, only : elmgeo_intersection_segment_TRI06
  use mod_elmgeo, only : elmgeo_cartesian_derivatives
  use mod_memory, only : memory_size
  use def_partis, only : leleboun_pts
  use def_partis, only : dimin_pts
  use def_partis, only : kfl_walld_pts
  use def_partis, only : kfl_fixbo_pts
  use def_partis, only : lboue_pts
  use def_partis, only : bouno_pts
  use def_partis, only : bvnat_pts  
  use def_partis, only : PTS_OUTFLOW_CONDITION
  use def_partis, only : PTS_WALL_CONDITION
  use def_partis, only : PTS_SLIP_CONDITION
  use def_partis, only : PTS_BOUNCING_CONDITION
  use def_partis, only : PTS_PARTICLE_HITS_WALL 
  use def_partis, only : PTS_PARTICLE_OUTFLOW
  use def_partis, only : lagrtyp
  implicit none

  integer(ip), intent(in)       :: ielem
  integer(ip), intent(in)       :: pnode
  integer(ip), intent(in)       :: ilagr
  real(rp),    intent(in)       :: diame
  real(rp),    intent(in)       :: hleng
  real(rp),    intent(in)       :: toler
  real(rp),    intent(in)       :: shapf(pnode)
  real(rp),    intent(in)       :: deriv(ndime,mnode)
  integer(ip), intent(out)      :: iwall
  integer(ip), intent(out)      :: iboun
  real(rp),    intent(inout)    :: coord_kp1(ndime)
  real(rp),    intent(inout)    :: veloc_kp1(ndime)
  real(rp),    intent(inout)    :: accel_kp1(ndime)
  real(rp),    intent(out)      :: xinte(ndime) 
  real(rp),    intent(out)      :: t 
  real(rp),    intent(out)      :: dt_k 
  integer(ip)                   :: pnodb,kboun,iinte,pblty,idime
  real(rp)                      :: facod(ndime,mnodb),p2(ndime),dista
  real(rp)                      :: s1,s2,udotn,vect(ndime),p1(ndime)
  real(rp)                      :: norma,vect_normal(ndime)
  real(rp)                      :: grwal(ndime),elwal(mnode),alpha,grnor
  real(rp)                      :: gpcar(ndime,mnode),elcod(ndime,mnode)
  real(rp)                      :: dista_min,xboun(3)
  
  iwall = 0
  
  !-------------------------------------------------------------------
  !
  ! Loop over boundaries connected to IELEM to check intersection
  !
  !-------------------------------------------------------------------
  
  if( lboue_pts(ielem) > 0 ) then

     iinte         = 0
     vect          = 0.0_rp
     vect(1:ndime) = coord_kp1(1:ndime)-lagrtyp(ilagr) % coord(1:ndime)
     dista         = sqrt(dot_product(vect,vect))
     s1            = (         0.5_rp * diame + toler * hleng ) / ( dista + zeror )
     s2            = ( dista + 0.5_rp * diame + toler * hleng ) / ( dista + zeror )
     !s1            = (         0.5_rp * diame ) / ( dista + zeror )
     !s2            = ( dista + 0.5_rp * diame ) / ( dista + zeror )
     p1(1:ndime)   = lagrtyp(ilagr) % coord(1:ndime) - vect(1:ndime) * s1
     p2(1:ndime)   = lagrtyp(ilagr) % coord(1:ndime) + vect(1:ndime) * s2

     loop_kboun: do kboun = 1,memory_size(leleboun_pts(ielem) % l)
        iboun                  = leleboun_pts(ielem) % l(kboun)
        pnodb                  = lnnob(iboun)
        pblty                  = ltypb(iboun)
        facod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
        call elmgeo_segfac(ndime,pnodb,facod,p1,p2,iinte,xinte,toler)
        udotn = dot_product(vect(1:ndime),bouno_pts(1:ndime,iboun))
        if( udotn < 0.0_rp ) iinte = 0
        if( iinte /= 0 ) exit loop_kboun
     end do loop_kboun

     !-------------------------------------------------------------------
     !
     ! Trajectory crosses boundary IBOUN at coordinate XINTE
     !
     !-------------------------------------------------------------------

     if( iinte /= 0 ) then

        if(      kfl_fixbo_pts(iboun) == PTS_OUTFLOW_CONDITION ) then
           !
           ! Outflow
           !
           iwall = 1
           lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_OUTFLOW 
           coord_kp1(1:ndime) = xinte(1:ndime)
           if( nbset > 0 ) lagrtyp(ilagr) % boundary_set = lbset(iboun)
           return
           
        else if( kfl_fixbo_pts(iboun) == PTS_WALL_CONDITION ) then
           !
           ! Wall
           !
           if( bvnat_pts(1,iboun) > zeror ) then
              call runend('PTS_SOLITE: PARTICLE BOUNCING NOT CODED')
           else
              iwall = 1
              lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_HITS_WALL
              coord_kp1(1:ndime) = xinte(1:ndime)
              if( nbset > 0 ) lagrtyp(ilagr) % boundary_set = lbset(iboun)
           end if
           return
           
        end if

     end if

  end if

  !-------------------------------------------------------------------
  !
  ! Wall test: particle will not be deposited on boundary deposition
  ! map DEPOB_PTS
  !
  !-------------------------------------------------------------------

  if( kfl_walld_pts > 0 ) then

     elwal(1:pnode) = walld(lnods(1:pnode,ielem))
     dista          = dot_product(shapf(1:pnode),elwal(1:pnode))
     
     if( dista < (0.5_rp-toler) * diame + dimin_pts ) then
        iwall = 1
        !
        ! Wall distance gradient to guess direction of wall
        !
        do idime = 1,ndime
           elcod(idime,1:pnode) = coord(idime,lnods(1:pnode,ielem))
        end do
        call elmgeo_cartesian_derivatives(ndime,pnode,elcod,deriv,gpcar)
        do idime = 1,ndime
           grwal(idime) = -dot_product(gpcar(idime,1:pnode),elwal(1:pnode))
        end do
        grnor = sqrt(dot_product(grwal,grwal))
        if( grnor > 0.0_rp ) grwal = grwal / grnor
        !
        ! Guess deposition time and intersection
        !
        lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_HITS_WALL
        vect(1:ndime)  = coord_kp1(1:ndime)-lagrtyp(ilagr) % coord(1:ndime)
        norma          = sqrt(dot_product(vect,vect))
        if( norma == 0.0_rp ) then
           t              = t
           xinte(1:ndime) = lagrtyp(ilagr) % coord(1:ndime)
        else
           vect_normal    = vect / norma
           alpha          = dot_product(grwal,vect_normal)
           if( alpha /= 0.0_rp ) dista = dista / alpha           
           xinte(1:ndime) = lagrtyp(ilagr) % coord(1:ndime) + dista * vect_normal(1:ndime)
           t              = t - dt_k
           dt_k           = dt_k * dista / norma
           t              = t + dt_k
        end if
        coord_kp1(1:ndime) = xinte(1:ndime)
        !
        ! Try to guess the boundary... only if we are in the first layer of elements
        !
        lagrtyp(ilagr) % boundary_set = -1
        
        if( nbset > 0 ) then
           
           if( memory_size(leleboun_pts(ielem) % l) > 0 ) then
              iboun = leleboun_pts(ielem) % l(1)
              lagrtyp(ilagr) % boundary_set = lbset(iboun)
           else 
              dista_min = huge(1.0_rp)
              iboun = 0
              do kboun = 1,nboun
                 pnodb = lnnob(kboun)
                 do idime = 1,ndime                 
                    xboun(idime) = sum(coord(idime,lnodb(1:pnodb,kboun)))
                    xboun(idime) = xboun(idime) / real(pnodb,rp)
                 end do
                 dista = dot_product(xboun(1:ndime)-xinte(1:ndime),xboun(1:ndime)-xinte(1:ndime))
                 if( dista < dista_min ) then
                    dista_min = dista
                    iboun = kboun
                 end if
              end do
              if( iboun > 0 ) lagrtyp(ilagr) % boundary_set = lbset(iboun)
           end if
           
        end if
        
     end if
  end if

end subroutine pts_cross_wall
