!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup AMR
!> @{
!> @file    mod_error_estimate.f90
!> @author  houzeaux
!> @date    2020-05-05
!> @brief   Error estimate
!> @details Some error estimation techniques
!-----------------------------------------------------------------------

module mod_error_estimate

  use def_kintyp_basic,        only : ip,rp,lg
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use def_master,              only : optional_argument
  use def_master,              only : zeror
  use def_master,              only : npart
  use def_master,              only : optional_argument
  use def_master,              only : INOTEMPTY
  use def_AMR,                 only : AMR_LAPLACIAN
  use def_AMR,                 only : AMR_SHARP_LAPLACIAN
  use def_AMR,                 only : AMR_DISTANCE_POINT
  use def_AMR,                 only : AMR_HESSIAN
  use def_AMR,                 only : amrp0
  use mod_matrix,              only : matrix_CSR_parallel_SpMV
  use mod_matrices,            only : matrices_laplacian
  use def_elmtyp
  use def_domain
  use mod_elmgeo
  use mod_communications
  use mod_memory
  use mod_maths
  use mod_gradie
  use mod_postpr
  implicit none

  public :: error_estimate_mesh_size
  public :: error_estimate_mesh_size_multipleSol
  public :: error_estimate_initialization
  public :: set_mesh_size
  private

  real(rp), pointer :: lapla(:)
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estimate_initialization()

    nullify(lapla)
    
  end subroutine error_estimate_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Mesh size
  !> @details Compute the error and then the optimum mesh size
  !> 
  !-----------------------------------------------------------------------

   subroutine error_estimate_mesh_size(solut,hh_opt,nelem_opt,METHOD,maxMeshSize,npoin_opt,isNodalField)
    
     use mod_metricComputation, only: compute_sizeField_from_sol_viaHessian
!     use def_master, only: kfl_paral
     

    real(rp),    pointer,  intent(in)    :: solut(:)
    real(rp),    pointer,  intent(inout) :: hh_opt(:)
    integer(ip), optional, intent(in)    :: NELEM_OPT
    integer(ip), optional, intent(in)    :: NPOIN_OPT
    integer(ip), optional, intent(in)    :: METHOD
    real(rp),    optional, intent(in)    :: maxMeshSize
    logical(lg), optional, intent(in)    :: isNodalField
    real(rp),    pointer                 :: hh(:)
    real(rp),    pointer                 :: ee(:)
    integer(ip)                          :: my_method
    
    integer(ip)                          :: ielem, pelty, pnode, ipoin, inode
    
    integer(ip)                          :: numTargetNodes
    
    my_method = optional_argument(AMR_SHARP_LAPLACIAN,METHOD)
   
    if(my_method == AMR_HESSIAN) then
       !print*,'compute_sizeField_from_sol_viaHessian, kfl_paral: ',kfl_paral
       
       if(.not.present(maxMeshSize)) then
         call runend('A maximum mesh size must be defined.. I will compute a default as maximum current mesh size * 10 (or some const)')
       end if
       
       if(present(npoin_opt)) then
         numTargetNodes = npoin_opt
       else if(present(nelem_opt)) then
         numTargetNodes = get_numNodes_fromNumElems_simplices(nelem_opt)
       else
         call runend('requires target elems/nodes')
         !call compute_sizeField_from_sol_viaHessian(hh,solut,ndime)
       end if
       
       call compute_sizeField_from_sol_viaHessian(hh,solut,ndime,maxMeshSize,numTargetNodes)
       
       if(npoin>0_ip) then
         
         if(isNodalField) then
           if( .not. associated(hh_opt) ) call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,npoin)
           hh_opt = hh
         else
           if( .not. associated(hh_opt) ) & 
                call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,nelem)
           hh_opt = 0.0_rp
           do ielem = 1,nelem
              pelty = ltype(ielem) 
              pnode = nnode(pelty)
              do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                hh_opt(ielem) = hh_opt(ielem) + hh(ipoin)
              end do
              hh_opt(ielem) = hh_opt(ielem) / pnode
           end do
         end if
         call memory_deallo(memor_dom,'HH','mod_error_estimate',hh)
      end if
       
    else if( my_method /= AMR_DISTANCE_POINT ) then 
       nullify(hh,ee)
       call memory_alloca(memor_dom,'EE','mod_error_estimate',ee,nelem)
       call memory_alloca(memor_dom,'HH','mod_error_estimate',hh,nelem)

       call error_estime_circumradius(hh)

       if(      my_method == AMR_SHARP_LAPLACIAN ) then
          call error_estime_hessian(solut,ee)
       else if( my_method == AMR_LAPLACIAN ) then
          call error_estime_laplacian_matrix(solut,ee)
       end if
       call error_estime_e_to_h(ee,hh,hh_opt,FACTOR=1.0_rp,NELEM_OPT=nelem_opt)

       !block
       !  use mod_postpr
       !  call postpr(ee,(/'ERROR','SCALA','NELEM'/),1_ip,1.0_rp)
       !end block

       call memory_deallo(memor_dom,'EE','mod_error_estimate',ee)
       call memory_deallo(memor_dom,'HH','mod_error_estimate',hh)
    else
       call error_estimate_distance_point(hh_opt,amrp0)
    endif
    
  end subroutine error_estimate_mesh_size
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gargallo
  !> @date    2022-05-16
  !> @brief   Mesh size
  !> @details Compute the error and then the optimum mesh size
  !> 
  !-----------------------------------------------------------------------

   subroutine error_estimate_mesh_size_multipleSol(solut,hh_opt,nelem_opt,METHOD,maxMeshSize,npoin_opt,isNodalField)
    
     use mod_metricComputation, only: compute_sizeField_from_sol_viaHessian
     use mod_metricComputation, only: compute_sizeField_from_multiSol_viaHessian
     

    real(rp),    pointer,  intent(in)    :: solut(:,:)
    real(rp),    pointer,  intent(inout) :: hh_opt(:)
    integer(ip), optional, intent(in)    :: NELEM_OPT
    integer(ip), optional, intent(in)    :: NPOIN_OPT
    integer(ip), optional, intent(in)    :: METHOD
    real(rp),    optional, intent(in)    :: maxMeshSize
    logical(lg), optional, intent(in)    :: isNodalField
    real(rp),    pointer                 :: hh(:)
    integer(ip)                          :: my_method
    
    integer(ip)                          :: numTargetNodes
    
    my_method = optional_argument(AMR_SHARP_LAPLACIAN,METHOD)
   
    if(my_method == AMR_HESSIAN) then
       
       if(.not.present(maxMeshSize)) then
         call runend('A maximum mesh size must be defined.. I will compute a default as maximum current mesh size * 10 (or some const)')
       end if
       
       if(present(npoin_opt)) then
         numTargetNodes = npoin_opt
       else if(present(nelem_opt)) then
         numTargetNodes = get_numNodes_fromNumElems_simplices(nelem_opt)
       else
         call runend('requires target elems/nodes')
       end if
       
       call compute_sizeField_from_multiSol_viaHessian(hh,solut,ndime,maxMeshSize,numTargetNodes)
       
       if(npoin>0_ip) then
         if( .not. associated(hh_opt) ) call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,npoin)
         hh_opt = hh
         call memory_deallo(memor_dom,'HH','mod_error_estimate',hh)
       end if
       
    else
      call runend('Mutiple solutions only available with hessian metric computation')
    endif
    
  end subroutine error_estimate_mesh_size_multipleSol
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gargallo
  !> @date    2022-02-16
  !> @brief   Elems to nodes estimation
  !> @details Elems to nodes estimation (based on structured simplicial grid, sufficiently large mesh)
  !> 
  !-----------------------------------------------------------------------
  
  function get_numNodes_fromNumElems_simplices(numElems) result(numNodes)
    integer(ip), intent(in) :: numElems
    integer(ip) :: numNodes
    
    if(ndime.eq.2_ip) then
      numNodes = ceiling(numElems/2.0_rp)
    else if(ndime.eq.3_ip) then
      numNodes = ceiling(numElems/6.0_rp)
    else
      call runend('error_estimate_mesh_size: not implemented for other dimensions')
    end if
    
  end function get_numNodes_fromNumElems_simplices
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gargallo
  !> @date    2022-02-16
  !> @brief   Mesh size
  !> @details Set optimum mesh size
  !> 
  !-----------------------------------------------------------------------
  subroutine set_mesh_size(solut,hh_opt,nelem_opt,METHOD,maxMeshSize,npoin_opt,isNodalField)
   
    use mod_metricComputation, only: compute_smoothSizeField
!    use def_master, only: kfl_paral
  
    use mod_out_paraview, only: out_paraview_inp
    use def_kintyp_mesh_basic
    use mod_strings,               only : integer_to_string
    use def_kermod

   real(rp),    pointer,  intent(in)    :: solut(:)
   real(rp),    pointer,  intent(inout) :: hh_opt(:)
   integer(ip), optional, intent(in)    :: NELEM_OPT
   integer(ip), optional, intent(in)    :: NPOIN_OPT
   integer(ip), optional, intent(in)    :: METHOD
   real(rp),    optional, intent(in)    :: maxMeshSize
   logical(lg), optional, intent(in)    :: isNodalField
   
   real(rp),    pointer                 :: hh(:)

   integer(ip)                          :: numTargetNodes
  
!   type(mesh_type_basic)        :: mesh_out
!   character(len=100) :: mesh_name_out
   
   nullify(hh)
   hh => solut
   
   if(present(npoin_opt)) then
     numTargetNodes = npoin_opt
   else if(present(nelem_opt)) then
     if(ndime.eq.2_ip) then
       numTargetNodes = ceiling(nelem_opt/2.0_rp)
     else if(ndime.eq.3_ip) then
       numTargetNodes = ceiling(nelem_opt/6.0_rp)
     else
       call runend('error_estimate_mesh_size: not implemented for other dimensions')
     end if
   else
     call runend('requires target elems/nodes')
   end if
   !
   call compute_smoothSizeField(hh_opt,hh,ndime,maxMeshSize,numTargetNodes)
   !
   
!    call mesh_out      % init      ('MESH_NEW')
!    call mesh_out      % copy      (meshe(ndivi))
!    mesh_name_out = './unitt/mesh_smoothSizeField'//integer_to_string(kfl_paral) !/unitt_adapti_hessian_alya
!    if(npoin>0_ip) call out_paraview_inp(mesh_out,filename=TRIM(mesh_name_out),nodeField=hh_opt)
!    call mesh_out      % deallo()
   !
!    nullify(hh_opt)
!    if(npoin>0_ip) then
!      call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,npoin)
!      hh_opt(:) = solut(:)
!    end if

 end subroutine set_mesh_size
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Compute the circumradius of the elements
  !>          Circumradius is only coded for TRI03, QUA04, TET04. For
  !>          other elements, the average edge length is taken
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estime_circumradius(hh)

    real(rp), pointer, intent(inout) :: hh(:)
    integer(ip)                      :: ielem
    integer(ip)                      :: pelty,pnode
    
    if( .not. associated(hh) ) & 
         call memory_alloca(memor_dom,'HH'    ,'mod_error_estimate',hh,nelem)
    
    do ielem = 1,nelem
       pelty     = ltype(ielem)
       pnode     = lnnod(ielem)
       hh(ielem) = elmgeo_circumradius(pelty,coord(1:ndime,lnods(1:pnode,ielem)))
   end do

  end subroutine error_estime_circumradius
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Compute the circumradius of the elements
  !>          Circumradius is only coded for TRI03, QUA04, TET04. For
  !>          other elements, the average edge length is taken
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estimate_distance_point(hh_opt,amrp0)

    real(rp),              pointer, intent(inout) :: hh_opt(:)
    real(rp),              intent(in)             :: amrp0(ndime)
    real(rp)                                      :: auxr(ndime)
    integer(ip)                                   :: ielem, pnode
    integer(ip)                                   :: inode, ipoin
    
    if( .not. associated(hh_opt) ) & 
         call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,nelem)
    
    do ielem = 1,nelem
       auxr = 0.0_rp
       pnode     = lnnod(ielem)
       do inode = 1,pnode
          ipoin                = lnods(inode,ielem)
          auxr(1:ndime)        = auxr(1:ndime) + ((0.1_rp/real(pnode)) * coord(1:ndime,ipoin))
       end do
       auxr(1:ndime) = amrp0(1:ndime)-auxr(1:ndime)
       hh_opt(ielem) = sqrt(dot_product(auxr,auxr))-0.45_rp
    end do


  end subroutine error_estimate_distance_point

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-08
  !> @brief   Laplacian matrix
  !> @details Error estimate based onn Laplacian matrix
  !>          Error = L u
  !> 
  !-----------------------------------------------------------------------
  
  subroutine error_estime_laplacian_matrix(uu,ee)

    real(rp),   pointer, intent(in)    :: uu(:)
    real(rp),   pointer, intent(inout) :: ee(:)
    integer(ip)                        :: ielem,pelty,inode,ipoin
    real(rp),   pointer                :: ee_node(:)

    nullify(ee_node)
    call memory_alloca(memor_dom,'EE_NODE','mod_error_estimate',ee_node,npoin)
    
    call matrices_laplacian(lapla,DEALLOCATE_MATRIX=.true.,NAME='LAPLA')
    if( INOTEMPTY ) then
       call matrix_CSR_parallel_SpMV(npoin,1_ip,1_ip,r_dom,c_dom,lapla,uu,ee_node)
       do ielem = 1,nelem
          pelty = abs(ltype(ielem))
          ee(ielem) = 0.0_rp
          do inode = 1,element_type(pelty) % number_nodes
             ipoin = lnods(inode,ielem)
             ee(ielem) = ee(ielem) + abs(ee_node(ipoin)) * elmar(pelty) % shacg(inode)
          end do
       end do       
    end if
    call memory_deallo(memor_dom,'LAPLA'  ,'mod_error_estimate',lapla)
    call memory_deallo(memor_dom,'EE_NODE','mod_error_estimate',ee_node)
    
  end subroutine error_estime_laplacian_matrix
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Convert error
  !> @details Convert the error into a mesh size
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estime_hessian(uu,ee)

    real(rp), pointer, intent(in)    :: uu(:)
    real(rp), pointer, intent(inout) :: ee(:)
    real(rp), pointer                :: gradu(:,:)
    integer(ip)                      :: ipoin,idime,inode,ielem
    integer(ip)                      :: pnode,pelty
    real(rp)                         :: gpcar(ndime,mnode,1)              ! dN/dxi
    real(rp)                         :: elgra(ndime,mnode)
    real(rp)                         :: elcod(ndime,mnode)
    real(rp)                         :: gplap

    nullify(gradu)
    call memory_alloca(memor_dom,'GRADU','mod_error_estimate',gradu,ndime,npoin)
    if( npoin > 0 ) call gradie(uu,gradu)

    if( .not. associated(ee) ) & 
         call memory_alloca(memor_dom,'EE','mod_error_estimate',ee,nelem)

    do ielem = 1,nelem
       pelty = ltype(ielem)
       pnode = lnnod(ielem)
       do inode = 1,pnode
          ipoin                = lnods(inode,ielem)
          elgra(1:ndime,inode) = gradu(1:ndime,ipoin)
          elcod(1:ndime,inode) = coord(1:ndime,ipoin)
       end do

       call elmgeo_cartesian_derivatives(ndime,pnode,elcod,elmar(pelty) % dercg,gpcar)
      
       gplap = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             gplap = gplap + gpcar(idime,inode,1) * elgra(idime,inode)
          end do
       end do
       ee(ielem) = abs(gplap)
       
    end do
    
  end subroutine error_estime_hessian
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Convert error
  !> @details Convert the error into a mesh size
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estime_e_to_h(ee,hh,hh_opt,FACTOR,NELEM_OPT)
    
    real(rp),              pointer, intent(in)    :: ee(:)
    real(rp),              pointer, intent(in)    :: hh(:)
    real(rp),              pointer, intent(inout) :: hh_opt(:)
    real(rp),    optional,          intent(in)    :: FACTOR
    integer(ip), optional,          intent(in)    :: NELEM_OPT
    real(rp)                                      :: rr,kk,nn,dd,xfact
    real(rp)                                      :: beta,gamma,alpha
    integer(ip)                                   :: ielem
    real(rp)                                      :: xx
    
    if( .not. associated(hh_opt) ) & 
         call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,nelem)
    !
    ! Target number of elements
    !
    xx = optional_argument(1.0_rp,FACTOR)
    if( present(NELEM_OPT) ) then
       nn = real(nelem_opt,rp)*xx
    else
       nn = real(nelem,rp)*xx
       call PAR_SUM(nn)
    end if
    
    dd    = real(ndime,rp)
    kk    = 2.0_rp
    alpha = 2.0_rp * kk / dd
    xfact = 1.0_rp/(dd*(1.0_rp+alpha))

    beta  = 0.0_rp
    gamma = 2.0_rp/(1.0_rp+alpha)
    !
    ! beta = Sum_{i=1}^N ei^{2/(1+alpha)}
    !
    do ielem = 1,nelem
       beta = beta + abs(ee(ielem)) ** gamma
    end do
    call PAR_SUM(beta)
     
    beta = beta * ( alpha ** ((2.0_rp+alpha)/(1.0_rp+alpha)) + alpha ** (1.0_rp/(1.0_rp+alpha)) )
    beta = ( (1.0_rp + alpha ) * nn / beta ) ** (1.0_rp/dd)
    beta = beta * ( alpha ** xfact )
   
    do ielem = 1,nelem      
       !print*, 'Checking hh...'
       !if (hh(ielem) < 0.0_rp) then
       !    print*, 'By Sigmar...', ielem, hh(ielem)
       !end if
       !print*, 'Done!'
       rr            = ( abs(ee(ielem)) ** (2.0_rp*xfact) ) * beta
       hh_opt(ielem) = hh(ielem) / (rr+zeror)
       !hh_opt(ielem) = abs(hh_opt(ielem))
    end do

  end subroutine error_estime_e_to_h

end module mod_error_estimate
!> @}
 
