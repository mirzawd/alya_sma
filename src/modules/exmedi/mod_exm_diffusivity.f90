!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_diffusivity.f90
!> @author  Adria Quintanas-Corominas
!> @author  Alfonso Santiago
!> @author  Mariano Vázquez
!> @date    02/2021
!> @brief   Diffusivity tensor calculations
!> @details Diffusivity tensor calculations
!-----------------------------------------------------------------------

module mod_exm_diffusivity

   use      def_parame
   use      def_master
   use      def_domain
   use      def_exmedi
   use      def_kermod,     only : kfl_fiber_long_fun, kfl_fiber_tang_fun, kfl_fiber_norm_fun
   use      mod_messages,   only : messages_live 
   use      mod_eccoupling
   use      mod_biofibers,  only : biofibers

   implicit none

   public :: exm_diffusivity_init
   public :: exm_diffusivity_at_gp_comcnd 
   public :: exm_diffusivity_at_gp_coucnd

   contains

      subroutine exm_diffusivity_init( )
         use mod_elmgeo,                    only : elmgeo_jacobian_matrix
         USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE, ieee_is_nan 

         !----------------------------------
         integer(ip)                      :: ielem, imate, inode, igaus, ipoin, idime, pelty, pnode, pgaus
         real(rp)                         :: gp_fbr(3,mgaus), gp_sht(3,mgaus), gp_nrm(3,mgaus) 
         real(rp)                         :: xbalo(3,3,mgaus),xnorm(3),xlong,xtra1,xtra2
         real(rp)                         :: xmean, xauxi
         real(rp)                         :: el_coord(ndime,mnode)
         real(rp)                         :: gp_factor, gp_W(mgaus), gp_N(mnode,mgaus), gp_dNde(ndime,mnode,mgaus)
         real(rp)                         :: gp_J_det
         !----------------------------------
         ! Check consistency
         if( ndime == 2_ip ) call runend('exm_diften_initialise: EXMEDI NOT PREPARED FOR 2D PROBLEMS')

         ! Initalise
         cedif_exm(:,:,:) = 0.0_rp         

         ! Compute nodal-diffusivity tensor
         elements: do ielem = 1, nelem

            ! Recover element inforamtion 
            pelty = ltype(ielem)

            active_elem: if( pelty > 0_ip )then
                imate = lmate(ielem)
                pnode = nnode(pelty)
                pgaus = ngaus(pelty)

                ! Gather coordinates
                do inode = 1, pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1, ndime
                       el_coord(idime,inode) = coord(idime,ipoin)
                    end do
                 end do

                ! Get interpolation
                gp_W(1:pgaus) = elmar(pelty) % weigp(1:pgaus)
                gp_N(1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
                gp_dNde(1:ndime,1:pnode,1:pgaus) = elmar(pelty) % deriv(1:ndime,1:pnode,1:pgaus)

                ! Compute reference fibers at GP
                call biofibers % get_reference_fibers_at_gp(&
                    ielem,pnode,pgaus,gp_N(1:pnode,1:pgaus),gp_fbr(1:3,pgaus),gp_sht(1:3,pgaus),gp_nrm(1:3,pgaus))


                ! Compute local basis matrix according to fibers
                if( kfl_ortho_diffusion_exm(imate) )then
                    ! - Orthotropic basis
                   
                    do igaus = 1, pgaus
                        ! Normalise fiber vector xbalo(:,1)
                        xauxi=sqrt(dot_product(gp_fbr(1:ndime,igaus),gp_fbr(1:ndime,igaus)))
                        if( xauxi > 1.0e-16_rp )then
                           xbalo(1:ndime,1,igaus) = gp_fbr(1:ndime,igaus)/xauxi
                        else
                           xbalo(1:ndime,1,igaus) = 0.0_rp 
                        endif

                        ! Sheet vector is read-in
                        xauxi=sqrt(dot_product(gp_sht(1:ndime,igaus),gp_sht(1:ndime,igaus)))
                        if( xauxi > 1.0e-16_rp )then
                           xbalo(1:ndime,2,igaus) = gp_sht(1:ndime,igaus)/xauxi
                        else
                           xbalo(1:ndime,2,igaus) = 0.0_rp 
                        endif

                        ! Sheet-normal vector is read-in
                        xauxi=sqrt(dot_product(gp_nrm(1:ndime,igaus),gp_nrm(1:ndime,igaus)))
                        if( xauxi > 1.0e-16_rp )then
                           xbalo(1:ndime,3,igaus) = gp_nrm(1:ndime,igaus)/xauxi
                        else
                           xbalo(1:ndime,3,igaus) = 0.0_rp 
                        endif
                    enddo

                    ! Set material parameters
                    xlong= gdiff_exm(1,1,imate)
                    xtra1= gdiff_exm(1,2,imate)
                    xtra2= gdiff_exm(1,3,imate)

                else
                    ! - Transversely isotropic basis

                    do igaus = 1, pgaus
                        ! Normalise fiber vector xbalo(:,1)
                        xauxi=sqrt(dot_product(gp_fbr(1:ndime,igaus),gp_fbr(1:ndime,igaus)))
                        if( xauxi > 1.0e-16_rp )then
                           xbalo(1:ndime,1,igaus) = gp_fbr(1:ndime,igaus)/xauxi
                        else
                           xbalo(1:ndime,1,igaus) = 0.0_rp 
                        endif

                        ! Evaluate arbitrary cross-fibre vector xbalo(:,2)
                        xnorm(1)=0.0_rp; xnorm(2)=1.0_rp; xnorm(3)=0.0_rp
                        xbalo(1,2,igaus) = xnorm(2)*xbalo(3,1,igaus)-xnorm(3)*xbalo(2,1,igaus)
                        xbalo(2,2,igaus) = xnorm(3)*xbalo(1,1,igaus)-xnorm(1)*xbalo(3,1,igaus)
                        xbalo(3,2,igaus) = xnorm(1)*xbalo(2,1,igaus)-xnorm(2)*xbalo(1,1,igaus)
                        xauxi=sqrt(dot_product(xbalo(1:3,2,igaus),xbalo(1:3,2,igaus)))

                        ! - in the case of parallel vectors
                        if( xauxi < 0.000001_rp )then
                            xnorm(1) = 1.0_rp; xnorm(2) = 0.0_rp; xnorm(3) = 0.0_rp
                            xbalo(1,2,igaus) = xnorm(2)*xbalo(3,1,igaus)-xnorm(3)*xbalo(2,1,igaus)
                            xbalo(2,2,igaus) = xnorm(3)*xbalo(1,1,igaus)-xnorm(1)*xbalo(3,1,igaus)
                            xbalo(3,2,igaus) = xnorm(1)*xbalo(2,1,igaus)-xnorm(2)*xbalo(1,1,igaus)
                            xauxi=sqrt(dot_product(xbalo(1:3,2,igaus),xbalo(1:3,2,igaus)))
                        end if
                        if( xauxi > 1.0e-16_rp )then
                           xbalo(1:ndime,2,igaus) = xbalo(1:ndime,2,igaus)/xauxi
                        else
                           xbalo(1:ndime,2,igaus) = 0.0_rp 
                        endif

                        ! axial vector
                        xbalo(1,3,igaus) = xbalo(2,1,igaus)*xbalo(3,2,igaus) - xbalo(3,1,igaus)*xbalo(2,2,igaus)
                        xbalo(2,3,igaus) = xbalo(3,1,igaus)*xbalo(1,2,igaus) - xbalo(1,1,igaus)*xbalo(3,2,igaus)
                        xbalo(3,3,igaus) = xbalo(1,1,igaus)*xbalo(2,2,igaus) - xbalo(2,1,igaus)*xbalo(1,2,igaus)
                        xauxi=sqrt(dot_product(xbalo(1:3,3,igaus),xbalo(1:3,3,igaus)))
                        if( xauxi > 1.0e-16_rp )then
                           xbalo(1:ndime,3,igaus) = xbalo(1:ndime,3,igaus)/xauxi
                        else
                           xbalo(1:ndime,3,igaus) = 0.0_rp 
                        endif

                    enddo

                    ! Set material parameters
                    xlong= gdiff_exm(1,1,imate)
                    xtra1= gdiff_exm(1,2,imate)
                    xtra2= xtra1

               endif

               ! From GP to NODES 
               do igaus = 1, pgaus
                   call elmgeo_jacobian_matrix(ndime,pnode,el_coord,gp_dNde(1:ndime,1:pnode,igaus),gp_J_det)
                   xmean = gp_W(igaus) * gp_J_det
                   do inode = 1, pnode
                       ipoin = lnods(inode,ielem)
                       gp_factor = gp_N(inode,igaus) * xmean

                       cedif_exm(1,1,ipoin) = cedif_exm(1,1,ipoin) + gp_factor * (&
                           xlong * xbalo(1,1,igaus) * xbalo(1,1,igaus) + &
                           xtra1 * xbalo(1,2,igaus) * xbalo(1,2,igaus) + &
                           xtra2 * xbalo(1,3,igaus) * xbalo(1,3,igaus) )
                       cedif_exm(2,2,ipoin) = cedif_exm(2,2,ipoin) + gp_factor * (&
                           xlong * xbalo(2,1,igaus) * xbalo(2,1,igaus) + &
                           xtra1 * xbalo(2,2,igaus) * xbalo(2,2,igaus) + &
                           xtra2 * xbalo(2,3,igaus) * xbalo(2,3,igaus) )
                       cedif_exm(3,3,ipoin) = cedif_exm(3,3,ipoin) + gp_factor * (&
                           xlong * xbalo(3,1,igaus) * xbalo(3,1,igaus) + &
                           xtra1 * xbalo(3,2,igaus) * xbalo(3,2,igaus) + &
                           xtra2 * xbalo(3,3,igaus) * xbalo(3,3,igaus) ) 
                       cedif_exm(1,2,ipoin) = cedif_exm(1,2,ipoin) + gp_factor * (&
                           xlong * xbalo(1,1,igaus) * xbalo(2,1,igaus) + &
                           xtra1 * xbalo(1,2,igaus) * xbalo(2,2,igaus) + &
                           xtra2 * xbalo(1,3,igaus) * xbalo(2,3,igaus) )
                       cedif_exm(1,3,ipoin) = cedif_exm(1,3,ipoin) + gp_factor * (&
                           xlong * xbalo(1,1,igaus) * xbalo(3,1,igaus) + &
                           xtra1 * xbalo(1,2,igaus) * xbalo(3,2,igaus) + &
                           xtra2 * xbalo(1,3,igaus) * xbalo(3,3,igaus) ) 
                       cedif_exm(2,3,ipoin) = cedif_exm(2,3,ipoin) + gp_factor * (&
                           xlong * xbalo(2,1,igaus) * xbalo(3,1,igaus) + & 
                           xtra1 * xbalo(2,2,igaus) * xbalo(3,2,igaus) + &
                           xtra2 * xbalo(2,3,igaus) * xbalo(3,3,igaus) )

                       cedif_exm(3,2,ipoin) = cedif_exm(2,3,ipoin)
                       cedif_exm(3,1,ipoin) = cedif_exm(1,3,ipoin)
                       cedif_exm(2,1,ipoin) = cedif_exm(1,2,ipoin)

                   enddo
               enddo 

            endif active_elem
         enddo elements

         ! Interchange with other workers
         call rhsmod(ndime*ndime,cedif_exm)

         ! Average
         do ipoin=1, npoin
             cedif_exm(1:ndime,1:ndime,ipoin) = cedif_exm(1:ndime,1:ndime,ipoin) / vmass(ipoin)
             if (any(ieee_is_nan(cedif_exm(1:ndime,1:ndime,ipoin)))) &
                  call runend("Some of cedif_exm (diffusivity tensor) values are NaN at the global point "//trim(intost(lninv_loc(ipoin))))
             if (.not. all(IEEE_IS_FINITE(cedif_exm(1:ndime,1:ndime,ipoin)))) &
                  call runend("Some of cedif_exm (diffusivity tensor) values are infinite at the global point "//trim(intost(lninv_loc(ipoin))))
         end do
         !----------------------------------
      end subroutine exm_diffusivity_init

      subroutine exm_diffusivity_get_current_tensor_at_gp( ielem, imate, difin, noion )
         implicit none
         integer(ip), intent(in)  :: ielem, imate
         real(rp),    intent(out) :: difin(ndime,ndime,mgaus)
         integer(ip), intent(out) :: noion
         integer(ip)              :: pelty, pnode, pgaus, jelem, igaus
         real(rp)                 :: gp_fbr(3,mgaus), gp_sht(3,mgaus), gp_nrm(3,mgaus)
         real(rp)                 :: xauxi, xbalo(3,3,mgaus),xnorm(3),xlong,xtra1,xtra2
         real(rp)                 :: gp_N(mnode,mgaus)
         ! Recover element information
         jelem = abs(ielem)
         pelty = ltype(jelem)
         pnode = nnode(pelty)
         pgaus = ngaus(pelty)
         noion = 0_ip
         ! 
         ! Get interpolation
         gp_N(1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)

         ! Compute reference fibers at GP
         call biofibers % get_current_fibers_at_gp( &
             ielem,pnode,pgaus,gp_N(1:pnode,1:pgaus),gp_fbr(1:3,pgaus),gp_sht(1:3,pgaus),gp_nrm(1:3,pgaus))


         ! Compute local basis matrix according to fibers
         if( kfl_ortho_diffusion_exm(imate) )then
             ! - Orthotropic basis
                   
             do igaus = 1, pgaus
                 ! Normalise fiber vector xbalo(:,1)
                 xauxi=sqrt(dot_product(gp_fbr(1:ndime,igaus),gp_fbr(1:ndime,igaus)))
                 if( xauxi > 1.0e-16_rp )then
                     xbalo(1:ndime,1,igaus) = gp_fbr(1:ndime,igaus)/xauxi
                 else
                     xbalo(1:ndime,1,igaus) = 0.0_rp 
                 endif

                 ! Sheet vector is read-in
                 xauxi=sqrt(dot_product(gp_sht(1:ndime,igaus),gp_sht(1:ndime,igaus)))
                 if( xauxi > 1.0e-16_rp )then
                    xbalo(1:ndime,2,igaus) = gp_sht(1:ndime,igaus)/xauxi
                 else
                    xbalo(1:ndime,2,igaus) = 0.0_rp 
                 endif

                 ! Sheet-normal vector is read-in
                 xauxi=sqrt(dot_product(gp_nrm(1:ndime,igaus),gp_nrm(1:ndime,igaus)))
                 if( xauxi > 1.0e-16_rp )then
                    xbalo(1:ndime,3,igaus) = gp_nrm(1:ndime,igaus)/xauxi
                 else
                    xbalo(1:ndime,3,igaus) = 0.0_rp 
                 endif
             enddo

             ! Set material parameters
             xlong= gdiff_exm(1,1,imate)
             xtra1= gdiff_exm(1,2,imate)
             xtra2= gdiff_exm(1,3,imate)

         else
             ! - Transversely isotropic basis

             do igaus = 1, pgaus
                 ! Normalise fiber vector xbalo(:,1)
                 xauxi=sqrt(dot_product(gp_fbr(1:ndime,igaus),gp_fbr(1:ndime,igaus)))
                 if( xauxi > 1.0e-16_rp )then
                    xbalo(1:ndime,1,igaus) = gp_fbr(1:ndime,igaus)/xauxi
                 else
                    xbalo(1:ndime,1,igaus) = 0.0_rp 
                 endif

                 ! Evaluate arbitrary cross-fibre vector xbalo(:,2)
                 xnorm(1)=0.0_rp; xnorm(2)=1.0_rp; xnorm(3)=0.0_rp
                 xbalo(1,2,igaus) = xnorm(2)*xbalo(3,1,igaus)-xnorm(3)*xbalo(2,1,igaus)
                 xbalo(2,2,igaus) = xnorm(3)*xbalo(1,1,igaus)-xnorm(1)*xbalo(3,1,igaus)
                 xbalo(3,2,igaus) = xnorm(1)*xbalo(2,1,igaus)-xnorm(2)*xbalo(1,1,igaus)
                 xauxi=sqrt(dot_product(xbalo(1:3,2,igaus),xbalo(1:3,2,igaus)))

                 ! - in the case of parallel vectors
                 if( xauxi < 0.000001_rp )then
                     xnorm(1) = 1.0_rp; xnorm(2) = 0.0_rp; xnorm(3) = 0.0_rp
                     xbalo(1,2,igaus) = xnorm(2)*xbalo(3,1,igaus)-xnorm(3)*xbalo(2,1,igaus)
                     xbalo(2,2,igaus) = xnorm(3)*xbalo(1,1,igaus)-xnorm(1)*xbalo(3,1,igaus)
                     xbalo(3,2,igaus) = xnorm(1)*xbalo(2,1,igaus)-xnorm(2)*xbalo(1,1,igaus)
                     xauxi=sqrt(dot_product(xbalo(1:3,2,igaus),xbalo(1:3,2,igaus)))
                 end if
                 if( xauxi > 1.0e-16_rp )then
                    xbalo(1:ndime,2,igaus) = xbalo(1:ndime,2,igaus)/xauxi
                 else
                    xbalo(1:ndime,2,igaus) = 0.0_rp 
                 endif

                 ! axial vector
                 xbalo(1,3,igaus) = xbalo(2,1,igaus)*xbalo(3,2,igaus) - xbalo(3,1,igaus)*xbalo(2,2,igaus)
                 xbalo(2,3,igaus) = xbalo(3,1,igaus)*xbalo(1,2,igaus) - xbalo(1,1,igaus)*xbalo(3,2,igaus)
                 xbalo(3,3,igaus) = xbalo(1,1,igaus)*xbalo(2,2,igaus) - xbalo(2,1,igaus)*xbalo(1,2,igaus)
                 xauxi=sqrt(dot_product(xbalo(1:3,3,igaus),xbalo(1:3,3,igaus)))
                 if( xauxi > 1.0e-16_rp )then
                    xbalo(1:ndime,3,igaus) = xbalo(1:ndime,3,igaus)/xauxi
                 else
                    xbalo(1:ndime,3,igaus) = 0.0_rp 
                 endif

             enddo

             ! Set material parameters
             xlong= gdiff_exm(1,1,imate)
             xtra1= gdiff_exm(1,2,imate)
             xtra2= xtra1

         endif

         do igaus = 1, pgaus
            difin(1,1,igaus) = &
                 xlong * xbalo(1,1,igaus) * xbalo(1,1,igaus) + &
                 xtra1 * xbalo(1,2,igaus) * xbalo(1,2,igaus) + &
                 xtra2 * xbalo(1,3,igaus) * xbalo(1,3,igaus)
             difin(2,2,igaus) = & 
                 xlong * xbalo(2,1,igaus) * xbalo(2,1,igaus) + &
                 xtra1 * xbalo(2,2,igaus) * xbalo(2,2,igaus) + &
                 xtra2 * xbalo(2,3,igaus) * xbalo(2,3,igaus)
             difin(3,3,igaus) = &
                 xlong * xbalo(3,1,igaus) * xbalo(3,1,igaus) + &
                 xtra1 * xbalo(3,2,igaus) * xbalo(3,2,igaus) + &
                 xtra2 * xbalo(3,3,igaus) * xbalo(3,3,igaus)
             difin(1,2,igaus) = & 
                 xlong * xbalo(1,1,igaus) * xbalo(2,1,igaus) + &
                 xtra1 * xbalo(1,2,igaus) * xbalo(2,2,igaus) + &
                 xtra2 * xbalo(1,3,igaus) * xbalo(2,3,igaus)
             difin(1,3,igaus) = & 
                 xlong * xbalo(1,1,igaus) * xbalo(3,1,igaus) + &
                 xtra1 * xbalo(1,2,igaus) * xbalo(3,2,igaus) + &
                 xtra2 * xbalo(1,3,igaus) * xbalo(3,3,igaus)
             difin(2,3,igaus) = &
                 xlong * xbalo(2,1,igaus) * xbalo(3,1,igaus) + &
                 xtra1 * xbalo(2,2,igaus) * xbalo(3,2,igaus) + &
                 xtra2 * xbalo(2,3,igaus) * xbalo(3,3,igaus)

             difin(3,2,igaus) = difin(2,3,igaus)
             difin(3,1,igaus) = difin(1,3,igaus)
             difin(2,1,igaus) = difin(1,2,igaus)
         enddo
       
      end subroutine exm_diffusivity_get_current_tensor_at_gp 

      subroutine exm_diffusivity_at_gp_comcnd(ielem, imate, difin, noion) 
         implicit none
         integer(ip), intent(in)  :: ielem, imate
         real(rp),    intent(out) :: difin(ndime,ndime,mgaus)
         integer(ip), intent(out) :: noion
         real(rp)                 :: xshap, difin_norm, aux
         integer(ip)              :: jelem,inode,igaus,ipoin,pelty,pnode,pgaus,idime, jdime

         if(ndime== 2)then
            call runend ('EXMEDI COMCOND not coded for 2D')
         end if
            
         if (kfl_cemod_exm == 2) then
            call runend('EXM_COMCND: NO BIDOMAIN MODEL WITH FIBERS PROGRAMMED')
         end if
            
         ! Recover element information
         jelem= ielem
         if (ielem < 0) jelem= -ielem
         pelty = ltype(jelem)
         pnode = nnode(pelty)
         pgaus = ngaus(pelty)

         ! Gather operations     
         difin_norm = 0.0_rp
         do igaus=1,pgaus
            difin(:,:,igaus) = 0.0_rp
            do inode = 1,pnode
               ipoin = lnods(inode,jelem)
               xshap = elmar(pelty)%shape(inode,igaus)
               difin(:,:,igaus) = difin(:,:,igaus) + xshap * cedif_exm(:,:,ipoin)
            end do

            do idime=1,ndime
               do jdime=1,ndime
                  difin_norm = difin_norm &
                     + difin(jdime,idime,igaus) * difin(jdime,idime,igaus)
               end do
            end do

         end do

         noion= 0
         if (difin_norm < 1.0e-10_rp) then
            aux= (gdiff_exm(1,1,imate)+gdiff_exm(1,2,imate))/2.0_rp
            do idime=1,ndime
               difin(idime,idime,imate)= aux
            end do
         end if

      end subroutine exm_diffusivity_at_gp_comcnd


      subroutine exm_diffusivity_at_gp_coucnd(ielem, imate, difin, nfibe, noion)

         implicit none

         integer(ip), intent(in)  :: ielem, imate
         real(rp),    intent(out) :: difin(ndime,ndime)
         real(rp),    intent(in)  :: nfibe(ndime)
         integer(ip), intent(out) :: noion
         integer(ip)              :: idime
         real(rp)                 :: xauxi,xrefi(12),xbalo(3,3),xnorm(3),xlong,xtra1,xtra2

         noion= 0_ip
         if (kfl_cellmod(imate) == EXMSLD_CELL_NOMODEL) then
            noion= 1_ip
            difin(1,1)= 0.01_rp*gdiff_exm(1,1,1)
            difin(2,2)= 0.01_rp*gdiff_exm(1,1,1)
            difin(ndime,ndime)= 0.01_rp*gdiff_exm(1,1,1)     
            return
         end if
         
         xrefi= 0.0_rp
         do idime= 1,ndime
            xrefi(idime)= nfibe(idime)
         end do
       
         xrefi(ndime+1)       = gdiff_exm(1,1,imate)  ! longitudinal intra
         xrefi(ndime+2)       = gdiff_exm(1,2,imate)  ! cross-wise (1) intra
         xrefi(ndime+ndime)   = gdiff_exm(1,3,imate)  ! cross-wise (2) intra
         
         xauxi = xrefi(1) * xrefi(1) + xrefi(2) * xrefi(2)                 
         if (ndime==3) xauxi= xauxi+xrefi(3)*xrefi(3)
         xauxi= sqrt(xauxi)
         
         xbalo= 0.0_rp
         if (xauxi > 1.0e-8_rp) then              
       
            xbalo(1,1) = xrefi(1)/xauxi
            xbalo(2,1) = xrefi(2)/xauxi
            if (ndime==3) xbalo(ndime,1) = xrefi(ndime)/xauxi
            
            if( ndime == 2 ) then
               
               xbalo(1,2)=   - xbalo(2,1)
               xbalo(2,2)=     xbalo(1,1)
       
               difin(1,1) = xrefi(ndime+1)*xbalo(1,1)*xbalo(1,1) &
                    + xrefi(ndime+2)*xbalo(2,1)*xbalo(2,1)
       
               difin(2,2) = xrefi(ndime+1)*xbalo(1,2)*xbalo(1,2) &
                    + xrefi(ndime+2)*xbalo(2,2)*xbalo(2,2)
       
               difin(1,2) = xrefi(ndime+1)*xbalo(1,2)*xbalo(1,1) &
                    + xrefi(ndime+2)*xbalo(2,1)*xbalo(2,2)
       
               difin(2,1) = difin(1,2)
                                 
            else if (ndime == 3) then
       
               xnorm(1)=0.0_rp
               xnorm(2)=1.0_rp
               xnorm(3)=0.0_rp
               
               xbalo(1,2) = xnorm(2)*xbalo(3,1)-xnorm(3)*xbalo(2,1)
               xbalo(2,2) = xnorm(3)*xbalo(1,1)-xnorm(1)*xbalo(3,1)
               xbalo(3,2) = xnorm(1)*xbalo(2,1)-xnorm(2)*xbalo(1,1)
               
               xauxi=xbalo(1,2)*xbalo(1,2)+xbalo(2,2)*xbalo(2,2)+xbalo(3,2)*xbalo(3,2)
               xauxi=sqrt(xauxi)
               
               if (xauxi < 0.000001_rp) then
                  xnorm(1) = 1.0_rp
                  xnorm(2) = 0.0_rp
                  xnorm(3) = 0.0_rp
                  
                  xbalo(1,2) = xnorm(2)*xbalo(3,1)-xnorm(3)*xbalo(2,1)
                  xbalo(2,2) = xnorm(3)*xbalo(1,1)-xnorm(1)*xbalo(3,1)
                  xbalo(3,2) = xnorm(1)*xbalo(2,1)-xnorm(2)*xbalo(1,1)
                  
                  xauxi = xbalo(1,2)*xbalo(1,2)+xbalo(2,2)*xbalo(2,2)+xbalo(3,2)*xbalo(3,2)
                  xauxi = sqrt(xauxi)                                  
               end if
               
               xbalo(1,2) = xbalo(1,2)/xauxi
               xbalo(2,2) = xbalo(2,2)/xauxi
               xbalo(3,2) = xbalo(3,2)/xauxi
               
               ! vector axial
               
               xbalo(1,3) = xbalo(2,1)*xbalo(3,2)-xbalo(3,1)*xbalo(2,2)
               xbalo(2,3) = xbalo(3,1)*xbalo(1,2)-xbalo(1,1)*xbalo(3,2)
               xbalo(3,3) = xbalo(1,1)*xbalo(2,2)-xbalo(2,1)*xbalo(1,2)
               
               xlong= xrefi(ndime+1)
               xtra1= xrefi(ndime+2)
               xtra2= xrefi(ndime+2)
               
               difin(1,1) = xlong*xbalo(1,1)*xbalo(1,1) &
                    + xtra1*xbalo(1,2)*xbalo(1,2) &
                    + xtra2*xbalo(1,3)*xbalo(1,3)
               difin(2,2) = xlong*xbalo(2,1)*xbalo(2,1) &
                    + xtra1*xbalo(2,2)*xbalo(2,2) &
                    + xtra2*xbalo(2,3)*xbalo(2,3)
               difin(3,3) = xlong*xbalo(3,1)*xbalo(3,1) &
                    + xtra1*xbalo(3,2)*xbalo(3,2) &
                    + xtra2*xbalo(3,3)*xbalo(3,3)
               difin(1,2) = xlong*xbalo(1,1)*xbalo(2,1) &
                    + xtra1*xbalo(1,2)*xbalo(2,2) &
                    + xtra2*xbalo(1,3)*xbalo(2,3)
               difin(1,3) = xlong*xbalo(1,1)*xbalo(3,1) &
                    + xtra1*xbalo(1,2)*xbalo(3,2) &
                    + xtra2*xbalo(1,3)*xbalo(3,3)
               difin(2,3) = xlong*xbalo(2,1)*xbalo(3,1) &
                    + xtra1*xbalo(2,2)*xbalo(3,2) &
                    + xtra2*xbalo(2,3)*xbalo(3,3)
               difin(3,2) = difin(2,3)
               difin(3,1) = difin(1,3)
               difin(2,1) = difin(1,2)
               
            end if
            
         end if

      end subroutine exm_diffusivity_at_gp_coucnd

end module mod_exm_diffusivity

