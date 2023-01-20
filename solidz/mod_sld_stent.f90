!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_stent.f90
!> @author  Eva Casoni
!> @date    June, 2019
!>          - Subroutine written
!> @brief   ToolBox for stent crimping and expansion simulation with balloon 
!>          expandable 
!>
!>          \verbatim
!>             To be integrated in a general manner, but fro testing pursposes
!>             now it is written in a module appart. 
!>
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------------

module mod_sld_stent

  use def_kintyp, only : ip, rp, lg

  implicit none

  public :: sld_calculation_contact_stent
  public :: sld_set_boundaries
  public :: sld_local_basis_stent
  public :: sld_calc_reaction_stent
  public :: sld_release_nodes_stent
  public :: sld_nodes_release_stent
  public :: sld_save_rhsid_stent
  public :: sld_csys_rotuni_stent
  private:: sld_project_sphere_on_plane
  private:: sld_spline_projection
contains



  !----------------------------------------------------------------------------
  !> @author  Eva Casoni
  !> @date    June, 2019
  !> @brief   Steps for crimping and expansion
  !> @details 
  !>         
  !----------------------------------------------------------------------------

  subroutine sld_calculation_contact_stent(ipoin,ibopo,ifunc,Rfunc)

    use def_master, only : displ
    use def_master, only : TIME_N, ITER_K, cutim
    use def_domain, only : ndime, coord
    use def_solidz, only : kfl_fixno_sld
    use def_solidz, only : jacrot_du_dq_sld, jacrot_dq_du_sld,fcont_sld
    use def_solidz, only : r_fin_stent, kfl_conta_stent
    use def_kermod, only : space_time_function
    !use mod_sld_csys, only: sld_csys_rotuni
 
    implicit none

    integer(ip), intent(in)  :: ipoin,ibopo,ifunc
    real(rp),    intent(in)  :: Rfunc(ndime)

    integer(ip)         :: idime, itott
    real(rp)            :: Rcurr, g, Rinit, Rnew
    real(rp)            :: dummy_matrix(ndime,ndime),dinew(ndime)
    real(rp), parameter :: tolc=0.0_rp
    real(rp), parameter :: tol=1.0E-14_rp
    !
    !
    ! Distances (projections)
    !
    Rinit = sqrt(coord(1,ipoin)**2 + coord(ndime,ipoin)**2)
    Rcurr = sqrt((coord(1,ipoin)+displ(1,ipoin,TIME_N))**2 + (coord(ndime,ipoin)+displ(ndime,ipoin,TIME_N))**2)
!!$    !
!!$    ! Initialize non-contact nodes
!!$    !
    do idime = 1, ndime
       if ( kfl_fixno_sld(idime,ipoin) == 3_ip ) kfl_fixno_sld(idime,ipoin) = 0_ip
    end do

    !
    ! Existance of contact: Gap is calculated taken the undeformed mesh.
    !
    itott = (ipoin-1) * ndime
    !
    !

  if ((space_time_function(ifunc) % name == 'CRIMP' .and. kfl_conta_stent == 1_ip)) then 


    if (Rfunc(1) >= r_fin_stent) then
       !
       ! Crimping
       ! 

       g = (1.0_rp + tolc)*Rfunc(1) - Rcurr
       if ( g < 0 ) then

 
         kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
          !
          ! Projection to the surface (in local axes)
          !
          dinew(1) = Rfunc(1) - Rcurr
          dinew(2) = 0.0_rp
          if ( ndime == 3 ) dinew(3) = 0.0_rp
          ! Local axes: Local --> Global (contact reads in global so rotate)
          ! dinew is in local normal-tangent basis, displ is in global framework.
          ! 1. rotate dinew to global
          call sld_csys_rotuni_stent(2_ip,ndime,ipoin,dinew(1:ndime))
          ! 2. impose it to displ, which is in global
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)
    
       else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
          ! i used to be a contact node... i am out of the contact radius
          kfl_fixno_sld(1      ,ipoin) = - 3_ip 
          kfl_fixno_sld(2:ndime,ipoin) =   0_ip  
       end if

     else if (Rfunc(1) < r_fin_stent) then
       kfl_conta_stent = -1_ip

    end if
    

  else if ((space_time_function(ifunc) % name == 'CRIMP' .and. kfl_conta_stent == -1_ip)) then
       !
       ! Free expansion
       !      
       g = (1.0_rp + tolc)*Rfunc(2) - Rinit

      ! if  ((g < 0) .and. fcont_sld(itott+1) < 0.0_rp) then!SI ESTO SIEMPRE SE CUMPLE PARA TODOS LOS NODOS DEL BOUNDARY,
       if  (g < 0) then!SI ESTO SIEMPRE SE CUMPLE PARA TODOS LOS NODOS DEL BOUNDARY,
                                                              !ENTONCES SIEMPRE ESTARAN EN CONTACTO

          kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
          !
          ! Projection to the surface (in local axes)
          !
          dinew(1) = Rfunc(2) - Rcurr
          dinew(2) = 0.0_rp
          if ( ndime == 3 ) dinew(3) = 0.0_rp
            ! Local axes: Local --> Global (contact reads in global so rotate)
          call sld_csys_rotuni_stent(2_ip,ndime,ipoin,dinew(1:ndime))
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)


         ! if (fcont_sld(itott+1) > 0.0_rp ) then   !ESTA FUERZA CREO QUE TIENE EL SIGNO CAMBIADO
         !    ! contact force is opposite to normal (negative), then, release the node
         !    kfl_fixno_sld(1      ,ipoin) =   0_ip
         !    kfl_fixno_sld(2:ndime,ipoin) =   0_ip  
         ! end if
      end if

     else if ((space_time_function(ifunc) % name == 'EXPAN'.and. kfl_conta_stent == 2_ip) .and. Rfunc(1) <= r_fin_stent) then   
   
       ! Forced expansion
       !
       g = (1.0_rp + tolc)*Rfunc(1) - Rcurr
       if ( g > 0 ) then


          ! 
          ! Compute spline: when the projection is not over a cilinder, we define the curve as a spline
          ! as a first aproximation
          ! 
          call sld_spline_projection(coord(2,ipoin), Rfunc(1), Rnew)        

          Rnew = Rnew/(5.0_rp - 4.0_rp*cutim)

          kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
          !
          ! Projection to the surface (in local axes)
          !
          !dinew(1) = (Rfunc(1) - Rcurr)*(1.0_rp + 0.01_rp*coord(2,ipoin))
          dinew(1) = Rfunc(1) - Rcurr
          dinew(2) = 0.0_rp
          if ( ndime == 3 ) dinew(3) = 0.0_rp
          ! Local axes: Local --> Global (contact reads in global so rotate)
          call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
               jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)

       else if (g<=0) then

          dinew(1) = 0.0_rp
          dinew(2) = 0.0_rp
          dinew(3) = 0.0_rp
          call sld_csys_rotuni_stent(2_ip,ndime,ipoin,dinew(1:ndime))
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)
          

       else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
          ! i used to be a contact node... i am out of the contact radius
          kfl_fixno_sld(1      ,ipoin) = - 3_ip 
          kfl_fixno_sld(2:ndime,ipoin) =   0_ip  

       end if

     else if ((space_time_function(ifunc) % name == 'EXPAN'.and. kfl_conta_stent == -2_ip) ) then   
       !
       ! Expansion unloading
       !
       g = (1.0_rp + tolc)*Rfunc(2) - Rinit
       if ( g > 0 ) then

          kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
          !
          ! Projection to the surface (in local axes)
          !
          dinew(1) = Rfunc(2) - Rcurr
          dinew(2) = 0.0_rp
          if ( ndime == 3 ) dinew(3) = 0.0_rp
          ! Local axes: Local --> Global (contact reads in global so rotate)
          call sld_csys_rotuni_stent(2_ip,ndime,ipoin,dinew(1:ndime))
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)

       else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
          ! i used to be a contact node... 

          if (fcont_sld(1,ipoin) <= 0.0_rp ) then
             ! contact force is opposite to normal (negative), then, release the node
             kfl_fixno_sld(1      ,ipoin) = - 3_ip
             kfl_fixno_sld(2:ndime,ipoin) =   0_ip  

          end if

       end if

     else if (space_time_function(ifunc) % name == 'CHARG') then
       !
       ! Charge expanded stent into cylinder
       !
       !if  ( (coord(2,ipoin)+displ(2,ipoin,TIME_N)) <= 0.0_rp .and. (coord(2,ipoin)+displ(2,ipoin,TIME_N)) > -4.0_rp ) then  
       if  ( (coord(2,ipoin)+displ(2,ipoin,TIME_N)) < 0.00_rp .and. (coord(2,ipoin)+displ(2,ipoin,TIME_N)) > -4.0_rp ) then  

          g = (1.0_rp + tolc)*Rfunc(1) - Rcurr

          if ( g < 0 ) then  

            kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
            !
            ! Projection to the surface (in local axes)
            !
            dinew(1) = (Rfunc(1) - Rcurr)
            dinew(2) = 0.0_rp
            if ( ndime == 3 ) dinew(3) = 0.0_rp
            
            ! Local axes: Local --> Global (contact reads in global so rotate)
            ! dinew is in local normal-tangent basis, displ is in global framework.
            ! 1. rotate dinew to global
            !call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
            !     jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
            call sld_csys_rotuni_stent(2_ip,ndime,ipoin,dinew(1:ndime))

            ! 2. impose it to displ, which is in global
            displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)
            !displ(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,ITER_K) + dinew(1:ndime)

         else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
            ! i used to be a contact node... i am out of the contact radius
            kfl_fixno_sld(1      ,ipoin) = - 3_ip 
            kfl_fixno_sld(2:ndime,ipoin) =   0_ip  
         end if

      !else if ( (coord(2,ipoin)+displ(2,ipoin,TIME_N)) <= -4.0_rp) then
      else if ( (coord(2,ipoin)+displ(2,ipoin,TIME_N)) <= -4.0_rp) then
            !g = (1.0_rp + tolc)*3.0_rp - Rcurr
            g = (1.0_rp + tolc)*3.00_rp - Rcurr

            if ( g < 0) then
              kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
              !
              ! Projection to the surface (in local axes)
              !
              !dinew(1) = (3.0_rp - Rcurr)
              dinew(1) = 3.00_rp - Rcurr
              dinew(2) = 0.0_rp
              if ( ndime == 3 ) dinew(3) = 0.0_rp
              ! Local axes: Local --> Global (contact reads in global so rotate)
              ! dinew is in local normal-tangent basis, displ is in global framework.
              ! 1. rotate dinew to global
              !call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
              !   jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
              call sld_csys_rotuni_stent(2_ip,ndime,ipoin,dinew(1:ndime))
              ! 2. impose it to displ, which is in global
              displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)
              !displ(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,ITER_K) + dinew(1:ndime)

           else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
              ! i used to be a contact node... i am out of the contact radius
              kfl_fixno_sld(1      ,ipoin) = - 3_ip 
              kfl_fixno_sld(2:ndime,ipoin) =   0_ip  
           end if
      end if

    else if (space_time_function(ifunc) % name == 'MOVEV') then
    !
    ! Move the cylinder instead of the stent
    !  
 
       if  ( (coord(2,ipoin)-cutim*12.5_rp) <= 0.0_rp .and. (coord(2,ipoin)-cutim*12.5_rp) > -4.0_rp ) then  
       
          g = (1.0_rp + tolc)*Rfunc(1) - Rcurr
    
          if ( g < 0 ) then  
      
           
            kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
            !
            ! Projection to the surface (in local axes)
            !
            dinew(1) = (Rfunc(1) - Rcurr)
            dinew(2) = 0.0_rp
            if ( ndime == 3 ) dinew(3) = 0.0_rp
            ! Local axes: Local --> Global (contact reads in global so rotate)
            ! dinew is in local normal-tangent basis, displ is in global framework.
            ! 1. rotate dinew to global
            call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
                 jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
            ! 2. impose it to displ, which is in global
            displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)
            displ(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,ITER_K) + dinew(1:ndime)

         else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
            ! i used to be a contact node... i am out of the contact radius
            kfl_fixno_sld(1      ,ipoin) = - 3_ip 
            kfl_fixno_sld(2:ndime,ipoin) =   0_ip  
         end if

      else if ( (coord(2,ipoin)-cutim*12.5_rp) <= -4.0_rp) then
            g = (1.0_rp + tolc)*3.0_rp - Rcurr

            if ( g < 0) then

              kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
              !
              ! Projection to the surface (in local axes)
              !
              dinew(1) = (3.0_rp - Rcurr)
              !dinew(1) = 3.49_rp - Rcurr
              dinew(2) = 0.0_rp
              if ( ndime == 3 ) dinew(3) = 0.0_rp
              ! Local axes: Local --> Global (contact reads in global so rotate)
              ! dinew is in local normal-tangent basis, displ is in global framework.
              ! 1. rotate dinew to global
              call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
                  jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
              ! 2. impose it to displ, which is in global
              displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)
              displ(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,ITER_K) + dinew(1:ndime)

           else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
              ! i used to be a contact node... i am out of the contact radius
              kfl_fixno_sld(1      ,ipoin) = - 3_ip 
              kfl_fixno_sld(2:ndime,ipoin) =   0_ip  
         end if
      end if
   else

       call runend('WE TRY CONTACT USING FUNCTIONS AND IT CAN BE POSSIBLE: '//trim(space_time_function(ifunc) % name))

    end if

  end subroutine sld_calculation_contact_stent


  subroutine sld_set_boundaries()
     use def_solidz, only : conbo_sld, contactbou_sld
     use def_solidz, only : kfl_funno_sld, kfl_conta_stent
     use def_solidz, only : kfl_contn_stent
     use def_solidz, only : kfl_fixrs_sld
     use def_domain, only : kfl_codno
     use def_domain, only : npoin

     implicit none

     integer(ip)          :: ipoin
          

          do ipoin = 1,npoin

            if (abs(kfl_conta_stent) == 1_ip) then
              !
              ! Crimping
              !
              contactbou_sld = conbo_sld(1)
              if ( kfl_codno(1,ipoin) == conbo_sld(1) ) then
                  kfl_funno_sld(ipoin) = -1_ip
                  kfl_contn_stent(ipoin) = 1_ip
              end if
              if ( kfl_codno(1,ipoin) == conbo_sld(2) ) kfl_funno_sld(ipoin) = 0_ip

            else if (abs(kfl_conta_stent) == 2_ip) then
              !
              ! Expansion
              !
              contactbou_sld = conbo_sld(2)
              if ( kfl_codno(1,ipoin) == conbo_sld(1) ) kfl_funno_sld(ipoin) = 0_ip
              if ( kfl_codno(1,ipoin) == conbo_sld(2) ) then
                  kfl_funno_sld(ipoin) = -1_ip
                  kfl_contn_stent(ipoin) = 1_ip   !nodes on the contact surface of the stent
              end if

            else if (abs(kfl_conta_stent) == 3_ip) then
              !
              ! Charge stent on cylinder
              !
              contactbou_sld = conbo_sld(1)
              if ( kfl_codno(1,ipoin) == conbo_sld(1) ) then
                  kfl_funno_sld(ipoin) = -1_ip
                  kfl_contn_stent(ipoin) = 1_ip   !nodes on the contact surface of the stent
                !  kfl_fixrs_sld(ipoin) = 10_ip
              end if
              !if ( kfl_codno(1,ipoin) == conbo_sld(2) ) kfl_funno_sld(ipoin) = 0_ip

            else if (abs(kfl_conta_stent) == 4_ip) then
              !
              ! Move cylinder
              !
              contactbou_sld = conbo_sld(1)
              if ( kfl_codno(1,ipoin) == conbo_sld(1) ) then
                  kfl_funno_sld(ipoin) = -1_ip
                  kfl_fixrs_sld(ipoin) = -1_ip
              end if
              if ( kfl_codno(1,ipoin) == conbo_sld(2) ) then 
                  kfl_funno_sld(ipoin) = 0_ip
                  kfl_fixrs_sld(ipoin) = 0_ip
              end if
              if ( kfl_codno(1,ipoin) == conbo_sld(2) )  kfl_funno_sld(ipoin) = 0_ip
            end if

          end do

        !


  end subroutine sld_set_boundaries

  subroutine sld_local_basis_stent(ipoin,ibopo,iroty,rotma)

     use def_kintyp,      only : ip, rp, lg
     use mod_maths,       only : maths_normalize_vector
     use mod_maths,       only : maths_vectorial_product
     use def_domain,      only : ndime,coord
     use def_solidz,      only : kfl_local_sld
     use def_solidz,      only : kfl_csysl_sld
     use def_solidz,      only : csysl_sld
     use def_solidz,      only : SLD_CSYS_CYLINDRICAL
 
     implicit none

     integer(ip), intent(in)  :: ipoin
     integer(ip), intent(in)  :: ibopo
     integer(ip), intent(in)  :: iroty
     real(rp),    intent(out) :: rotma(ndime,ndime)     

     real(rp)     :: rho
     real(rp)     :: R(ndime), A(ndime), T(ndime)
     real(rp)     :: localbasis(ndime,ndime)     

 
     if (kfl_local_sld == 1_ip) then

        localbasis = 0.0_rp

        if (kfl_csysl_sld == SLD_CSYS_CYLINDRICAL) then
           !
           ! Only 3D
           !
           if (ndime == 2_ip) then

           else if (ndime == 3_ip) then 
                   !
                   ! Radial axis
                   !
                   if( abs(csysl_sld(3*ndime)) > 0.0_rp ) then
                      !
                      ! X-Y Plane (Axial=Z)
                      !
                      R(1) = coord(1,ipoin) - csysl_sld(1)
                      R(2) = coord(2,ipoin) - csysl_sld(2)
                      R(3) = 0.0_rp
                      call maths_normalize_vector(ndime,R(:),rho)

                   else if( abs(csysl_sld(3*ndime-1)) > 0.0_rp  ) then
                      !
                      ! X-Z Plane (Axial=Y)
                      !
                      R(1) = coord(1,ipoin) - csysl_sld(1)
                      R(2) = 0.0_rp
                      R(3) = coord(3,ipoin) - csysl_sld(3)
                      call maths_normalize_vector(ndime,R(:),rho)

                   else if( abs(csysl_sld(3*ndime-2)) > 0.0_rp ) then
                      !
                      ! Y-Z Plane (Axial=X)
                      !
                      R(1) = 0.0_rp
                      R(2) = coord(2,ipoin) - csysl_sld(2)
                      R(3) = coord(3,ipoin) - csysl_sld(3)
                      call maths_normalize_vector(ndime,R(:),rho)

                   end if
                   !
                   ! Axial axis
                   A(1) = csysl_sld(3*ndime-2)
                   A(2) = csysl_sld(3*ndime-1)
                   A(3) = csysl_sld(3*ndime)
                   call maths_normalize_vector(ndime,A(:),rho)

                   ! Local basis
                   localbasis(:,1) = R(:)
                   localbasis(:,3) = A(:)
                   ! T = A x R
                   call maths_vectorial_product(A(:), R(:), T(:), 3_ip)
                   localbasis(:,2) = T(:)

           end if

           rotma = localbasis

        end if

     end if    

  end subroutine sld_local_basis_stent
   
  subroutine sld_nodes_release_stent()
    use mod_communications, only : PAR_SUM
    use def_solidz,         only : kfl_goite_sld
    use def_solidz,         only : new_to_release

    implicit none

    !
    ! Check existance of sticky nodes
    !
    if ( kfl_goite_sld /= 1_ip ) then
 
      call sld_calc_reaction_stent()

      if (new_to_release == 1_ip) then
        new_to_release = 0_ip
        call sld_release_nodes_stent()    !EVA: creo que para Asto no hace falta funcion
        
        kfl_goite_sld = 1_ip
      end if

    end if
    !
    ! Check if all subdomains are converged to the solution
    ! If the sum is greater than 1 (not converged) repeat iteration
    call PAR_SUM(kfl_goite_sld,'IN MY CODE')
    if (kfl_goite_sld >=1_ip) kfl_goite_sld = 1_ip  

  end subroutine sld_nodes_release_stent



  subroutine sld_calc_reaction_stent()

    use def_domain, only : lpoty
    use def_master, only : INOTMASTER, cutim
    use def_master, only : momod, modul, soltyp
    use def_domain, only : npoin, ndime
!    use def_domain, only : kfl_codno, npoin, ndime
!    use def_solidz, only : contactbou_sld, release_nodes, new_to_release
    use def_solidz, only : fcont_sld, frxid_sld
    use def_solidz, only :  kfl_fixno_sld
    use def_solidz, only : jacrot_du_dq_sld, jacrot_dq_du_sld

    implicit none

    integer(ip)           :: ipoin, itott, ibopo
    real(rp)              :: dummy_matrix(ndime,ndime)
!    real(rp)              :: fcont(ndime*npoin)
    type(soltyp), pointer :: solve(:)

    solve => momod(modul) % solve(1:)

    if ( INOTMASTER ) then
       do ipoin = 1, npoin  !Creo que tiene que ser bucle en boundary nodes
          itott = (ipoin-1) * ndime
          ibopo = lpoty(ipoin)
 
       if (cutim > 1.0_rp .and. kfl_fixno_sld(1,ipoin) == 3_ip) then
 
        !  if (any(kfl_codno(:,ipoin) == contactbou_sld) ) then

            !fcont(itott+1:itott+ndime) = frxid_sld(itott+1:itott+ndime) !este es el q tiene Gerard  
            fcont_sld(1:ndime,ipoin) = frxid_sld(itott+1:itott+ndime)   

          !  if (fcont_sld(itott+1) > 0.0_rp ) then
          !    kfl_fixno_sld(1,ipoin) = 0_ip
          !  end if

        !    if (fcont(itott+1) >= 0.0_rp .and. release_nodes(ipoin) == 0_ip) then
        !      release_nodes(ipoin) = 1_ip
        !      new_to_release = 1_ip
        !    end if
            !
            ! Rotate contact force (Local --> Global)
            !
            call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,fcont_sld(:,ipoin), &
                 jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
        !  end if
       end if

       end do
    end if

  end subroutine sld_calc_reaction_stent



  subroutine sld_release_nodes_stent()

    use def_master, only : INOTMASTER
    use def_domain, only : kfl_codno, npoin, ndime
    use def_solidz, only : release_nodes, kfl_fixno_sld, contactbou_sld

    implicit none

    integer(ip)          :: ipoin

    if ( INOTMASTER ) then
       do ipoin = 1, npoin
          if ( kfl_codno(1,ipoin) == contactbou_sld .and. release_nodes(ipoin) == 1 ) then
             kfl_fixno_sld(1:ndime,ipoin) = 0_ip
          end if
       end do
    end if

  end subroutine sld_release_nodes_stent


  subroutine sld_save_rhsid_stent()

    use def_master, only : INOTMASTER, rhsid
    use def_master, only : momod, modul, soltyp
    use def_domain, only : npoin_own, ndime, lpoty
!    use def_domain, only : kfl_codno, npoin
    use def_solidz, only : contact_friction_sld
!    use def_solidz, only : saved_rhsid, contactbou_sld
    use def_solidz, only : kfl_fixno_sld
!    use def_solidz, only : fcont_sld

    implicit none

    integer(ip)           :: ipoin, idime, ibopo, itott
    type(soltyp), pointer :: solve(:)
    real(rp)              :: rhsid_ipoin(ndime),rhsid_ipoin_tangent_norm,rhsid_ipoin_tangent_normalized(ndime),friction_force

    solve => momod(modul) % solve(1:)

    if ( INOTMASTER ) then
       !if (.not. associated(saved_rhsid)) then
       !   allocate(saved_rhsid(ndime,npoin))
       !else
       !   print *,"ERROR ALLOCATING SAVE_RHSID. ALLOCATING WHEN ASSOCIATED"
       !end if
  !     saved_rhsid = 0.0_rp
  !     fcont_sld = 0.0_rp
  !     do ipoin = 1, npoin
  !        ibopo = lpoty(ipoin)
  !        itott = (ipoin-1)*ndime
  !        if ( kfl_codno(1,ipoin) == contactbou_sld ) then
  !           do idime = 1, ndime
  !              fcont_sld(itott+idime) =  rhsid(itott+idime)
  !           end do
 
 !           if ( ndime == 2_ip ) then
  !              saved_rhsid(1,ipoin) = -rhsid(2*ipoin-1)
  !              saved_rhsid(2,ipoin) = -rhsid(2*ipoin  )
  !           else
  !              saved_rhsid(1,ipoin) = -rhsid(3*ipoin-2)
  !              saved_rhsid(2,ipoin) = -rhsid(3*ipoin-1)
  !              saved_rhsid(3,ipoin) = -rhsid(3*ipoin )
  !           end if
  !        end if
  !     end do

       if (contact_friction_sld > 0.0_rp .and. ndime > 2_ip) then

          !npoin_own: nodos de frontera propios, que no pertenecen a otro subdominio

          do ipoin = 1, npoin_own ! do only for own nodes, because it modify the residual 
             ibopo = lpoty(ipoin)
             itott = (ipoin-1)*ndime
             !
             ! Friction
             !
             ! 1. check the fixno contact condition
             if (kfl_fixno_sld(1,ipoin) == 3_ip) then             
                                
                ! store the contact force in a temporary local vector to be used later
                do idime=1,ndime
                   rhsid_ipoin(idime) = rhsid(itott+idime)
                end do
                
                if (contact_friction_sld > 1.0_rp) then

                   ! sticky always
                   kfl_fixno_sld(2,ipoin)= 3_ip
                   kfl_fixno_sld(3,ipoin)= 3_ip
                    
                else
                   
                   ! 2. compute a normalized tangential force (as fixno==3, it is already in the local framework because
                   ! it was previously rotated in elmope)
                   rhsid_ipoin_tangent_norm= sqrt(rhsid_ipoin(2)*rhsid_ipoin(2) + rhsid_ipoin(3)*rhsid_ipoin(3))
                   rhsid_ipoin_tangent_normalized(1) = 0.0_rp               !normal direction, not used
                   if (rhsid_ipoin_tangent_norm > 0.0_rp) then
                      rhsid_ipoin_tangent_normalized(2) = rhsid_ipoin(2)/rhsid_ipoin_tangent_norm     !first tangent
                      rhsid_ipoin_tangent_normalized(3) = rhsid_ipoin(3)/rhsid_ipoin_tangent_norm     !second tangent
                   end if
                   ! 3. compute friction force (mu * N)
                   ! check if the node force is negative (positive rhs), i.e. acting against the surface
                   ! if positive (negative rhs), the normal force is opposing the contact, so there should not be friction
                   friction_force = 0.0_rp
                   if (rhsid_ipoin(1) > 0.0_rp ) friction_force = contact_friction_sld * abs(rhsid_ipoin(1)) 
                   
                   ! 4. check if stick (mu * N > tangential force) or slip (the opposite)
                   if (rhsid_ipoin_tangent_norm > 0.0_rp) then
                      if (rhsid_ipoin_tangent_norm > friction_force) then
                         ! slip: local tangential force is larger than friction
                         ! 5. release tangent degrees of freedom
                         kfl_fixno_sld(2,ipoin)= 0_ip
                         kfl_fixno_sld(3,ipoin)= 0_ip
                         ! 6. correct tangent forces ( f_t = f_t - mu * N)
                         rhsid_ipoin(2) = rhsid_ipoin(2) - rhsid_ipoin_tangent_normalized(2) * abs(friction_force)     !first tangent
                         rhsid_ipoin(3) = rhsid_ipoin(3) - rhsid_ipoin_tangent_normalized(3) * abs(friction_force)    !second tangent
                         ! 7. reassign tangent rhsid for this node
                         rhsid((itott+2)) = rhsid_ipoin(2)                     
                         rhsid((itott+3)) = rhsid_ipoin(3)                     
                      else
                         ! stick: otherwise
                         ! 5. fix tangent boundary conditions
                         kfl_fixno_sld(2,ipoin)= 3_ip  !ESTO ES REDUNDANTE, NO HARIA FALTA PONERLO
                         kfl_fixno_sld(3,ipoin)= 3_ip
                      end if                      
                   end if
                end if

             end if
             
          end do

       end if
! we don't need this because rhsid is recomputed with local information
!       call rhsmod(ndime,fcont_sld) ! forces
!       call rhsmod(ndime,rhsid    ) ! rhsid

    end if

  end subroutine sld_save_rhsid_stent


  subroutine sld_csys_rotuni_stent(itask,ndofn,ipoin,unrot)

    use def_master, only : INOTMASTER
    use def_domain, only : lpoty
    use def_solidz, only : jacrot_du_dq_sld

    implicit none

    integer(ip), intent(in)    :: itask          !< Transformation
    integer(ip), intent(in)    :: ndofn          !< Dimensions
    integer(ip), intent(in)    :: ipoin          !< Point id
    real(rp),    intent(inout) :: unrot(ndofn)   !< Vector transformed

    integer(ip)                :: ibopo
    real(rp)                   :: worma(ndofn), worve(ndofn),romatt(ndofn,ndofn)

    if ( INOTMASTER ) then

       select case (itask)

       case( 1_ip )
          !
          ! Global to local A = B^T * C
          !
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             !
             ! Boundary conditions in the user local system
             !
             worve(1:ndofn) = unrot
             romatt = transpose(jacrot_du_dq_sld(:,:,ipoin))
             worma = matmul(romatt,worve)
             unrot(1:ndofn) = worma(1:ndofn)
          else
             call runend('MOD_SLD_CSYS: THIS POINT DOES NOT HAVE A ROTATION MATRIX FOR LOCAL AXES')
          end if

       case( 2_ip )
          !
          ! Local to global: A = B * C
          !
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             !
             ! Boundary conditions in the local system
             !
             worve(1:ndofn) = unrot(1:ndofn)
             worma = matmul(jacrot_du_dq_sld(:,:,ipoin),worve)
             unrot(1:ndofn) = worma(1:ndofn)
          else
             call runend('MOD_SLD_CSYS: THIS POINT DOES NOT HAVE A ROTATION MATRIX FOR LOCAL AXES')
          end if

       end select

    end if

end subroutine sld_csys_rotuni_stent

  !----------------------------------------------------------------------------
  !> @author  Eva Casoni
  !> @date    June, 2019
  !> @brief   Steps for crimping and expansion
  !> @details 
  !>         
  !----------------------------------------------------------------------------

  subroutine sld_calculation_contact_identer(ipoin,ibopo,ifunc,Rfunc)


! In this case Rfunc is only the displacement in vertical direction: x and z position do not vary. Always in global coordinates

    use def_master, only : displ
    use def_master, only : TIME_N, cutim
    use def_domain, only : ndime, coord
    use def_solidz, only : kfl_fixno_sld
    use def_solidz, only : jacrot_du_dq_sld, jacrot_dq_du_sld
!    use def_solidz, only : fcont_sld

    implicit none

    integer(ip), intent(in)  :: ipoin,ibopo,ifunc
    real(rp),    intent(in)  :: Rfunc(ndime)

    integer(ip)         :: itott
!    integer(ip)         :: idime
    real(rp)            :: Rcurr, g, Rinit
    real(rp)            :: dummy_matrix(ndime,ndime),dinew(ndime)
    real(rp)            :: coord_init(ndime), coord_curr(ndime), y_sphere
    real(rp), parameter :: tolc=0.0_rp
    real(rp), parameter :: tol=1.0E-14_rp
    !
    !
    ! Distances (projections)
    !
    coord_init = coord(:,ipoin)
    coord_curr = coord(:,ipoin) + displ(:,ipoin,TIME_N)

    call sld_project_sphere_on_plane(coord_init, coord_curr, dinew(2), y_sphere)

!!$    !
!!$    ! Initialize non-contact nodes
!!$    !
  !  do idime = 1, ndime
  !     if ( kfl_fixno_sld(idime,ipoin) == 3_ip ) kfl_fixno_sld(idime,ipoin) = 0_ip
  !  end do

    !
    ! Existance of contact: Gap is calculated taken the undeformed mesh.
    !
    itott = (ipoin-1) * ndime
    !
    if ( cutim <= 1.0_rp + tol ) then
       !
       ! Compress
       ! 
       g = (1.0_rp + tolc)*y_sphere - coord_curr(2)
       if ( g < 0 ) then

          kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
          !
          ! Projection to the surface (in local axes)
          !
          dinew(1) = 0.0_rp
          if ( ndime == 3 ) dinew(3) = 0.0_rp
         ! Local axes: Local --> Global (contact reads in global so rotate)
          ! dinew is in local normal-tangent basis, displ is in global framework.
          ! 1. rotate dinew to global
          call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
               jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
          ! 2. impose it to displ, which is in global
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)


       else if (abs(kfl_fixno_sld(1,ipoin)) == 3) then
          ! i used to be a contact node... i am out of the contact radius
          kfl_fixno_sld(1      ,ipoin) = - 3_ip 
          kfl_fixno_sld(2:ndime,ipoin) =   0_ip  

       end if


    else if ( cutim <= 2.0_rp + tol .and. cutim > 1.0_rp + tol ) then

       !
       ! Free expansion
       !      

       g = (1.0_rp + tolc)*Rfunc(2) - Rinit
!       g = (1.0_rp + tolc)*Rfunc(2) - Rcurr

       !if  ((g < 0) .and. fcont_sld(itott+1) < 0.0_rp) then!SI ESTO SIEMPRE SE CUMPLE PARA TODOS LOS NODOS DEL BOUNDARY,
       if  (g < 0) then!SI ESTO SIEMPRE SE CUMPLE PARA TODOS LOS NODOS DEL BOUNDARY,
                                                              !ENTONCES SIEMPRE ESTARAN EN CONTACTO

          kfl_fixno_sld(1,ipoin) = 3_ip ! flag contact
          !
          ! Projection to the surface (in local axes)
          !
          dinew(1) = Rfunc(2) - Rcurr
          dinew(2) = 0.0_rp
          if ( ndime == 3 ) dinew(3) = 0.0_rp
            ! Local axes: Local --> Global (contact reads in global so rotate)
          call sld_rotsys(2_ip,1_ip,1_ip,ndime,ndime,dummy_matrix,dinew(1:ndime), &
               jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
          displ(1:ndime,ipoin,TIME_N) = displ(1:ndime,ipoin,TIME_N) + dinew(1:ndime)

         ! if (fcont_sld(itott+1) < 0.0_rp) then
         !     kfl_fixno_sld(1,ipoin) = 0_ip
         ! end if

       end if

    else



       call runend('WE TRY CONTACT USING FUNCTIONS AND IT CAN BE POSSIBLE')

    end if

  end subroutine sld_calculation_contact_identer

  !-------------------------------------------------------------------------
  !
  !> @author  Eva Casoni
  !> @date    June, 2019
  !> @brief   Subroutine to projects points of a sphere in a plane
  !> @details 
  !
  !-------------------------------------------------------------------------


  subroutine sld_project_sphere_on_plane(coord_init, coord_curr, disp_y, y_sphere)


    use def_parame, only : pi
    use def_domain, only : ndime

    implicit none

    real(rp),    intent(in)  :: coord_init(ndime), coord_curr(ndime)
    real(rp),    intent(in)  :: disp_y
    real(rp),    intent(out) :: y_sphere

    real(rp)            :: Phi, thita, alpha

    ! Having coordinate x and coordinate z of the point in the plane, they are the same
    ! for the sphere, hence compute thita and phi for the parametrization
    !
    ! x = R cos(thita)sin(Phi)
    ! y = alpha + R sin(thita)sin(Phi)
    ! z = R cos(Phi)   

    Phi = acos(coord_curr(3))*180_rp/pi
    thita = acos(coord_curr(1)/sin(Phi))*180_rp/pi
   
    alpha = 1.0_rp - disp_y                   !new center of sphere, y coordinate
    y_sphere = alpha + sin(thita)*sin(Phi)

  end subroutine sld_project_sphere_on_plane

  !-------------------------------------------------------------------------
  !
  !> @author  Eva Casoni
  !> @date    June, 2019
  !> @brief   Subroutine to projects points of a sphere in a plane
  !> @details 
  !
  !-------------------------------------------------------------------------
  subroutine  sld_spline_projection(y, Rmax, Rnew)

     
    ! The spline is defined by 4 points, is the revolution surface of the mandril sent 
    ! by iVascular for the Nitinol stent deformation
    ! I am approximating it by straight splines as a first trial
    ! The points are (y,R): (0,10), (2,14), (3,12), (4, 14)
  
    real(rp),    intent(in)  :: y, Rmax
    real(rp),    intent(out) :: Rnew

    if (y <= 25.0_rp) then
       Rnew = (4.0_rp/25.0_rp)*y + 10.0_rp

    else if (y > 25.0_rp .and. y <= 32.0_rp) then
       Rnew = (123.0_rp/7.0_rp) - (1.0_rp/7.0_rp)*y

    else if (y > 32.0_rp) then
       Rnew = 9.8_rp + 0.1_rp*y

    end if

  end subroutine sld_spline_projection
end module mod_sld_stent
