!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_outvar.f90
!> @author  Guillaume Houzeaux
!> @date    20/02/2013
!> @brief   Output a postprocess variable
!> @details Output a postprocess variable for gusano.
!> @} 
!-----------------------------------------------------------------------
subroutine gus_outvar(jvari,imesh)
  
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_gusano
  use mod_arrays
  use mod_outvar
  use mod_maths_basic, only : maths_normalize_vector
  use mod_elmgeo,      only : element_type
  use mod_memory_basic
  use mod_witness
  use mod_maths
  use mod_communications_point_to_point
  implicit none
  
  integer(ip), intent(in) :: jvari   !< 1 for veloc , 2 press , etc...
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ivari
  real(rp)                :: nn(3),xx(3)
  real(rp)                :: radiu,R,U
  real(rp),    pointer    :: vecto(:,:)
  integer(ip)             :: inode,ielem,ipoin,pelty
  integer(ip)             :: kpoin,iline,icont
  !
  ! Define postprocess variable
  !
  ivari = abs(jvari)

  select case ( arrays_name(ivari) )  

  case ( 'VEL1D' )
     !
     ! VEL1D: Velocity
     !
     gesca => vel1d(:,1)
     call arrays(ivari,'POSTPROCESS',gesca,MESH_ID=imesh)
     
  case ( 'VELOC' )
     !
     ! VELOC: Velocity
     !
     do ielem = 1,nelem
        pelty = ltype(ielem)
        if( abs(pelty) < element_end ) then
           nn(1:ndime) = (coord(1:ndime,lnods(2,ielem))-coord(1:ndime,lnods(1,ielem)))
           call maths_normalize_vector(ndime,nn)
           do inode = 1,element_type(pelty) % number_nodes
              ipoin = lnods(inode,ielem)
              veloc(1:ndime,ipoin,1) = vel1d(ipoin,1) * nn(1:ndime) 
           end do
        end if
     end do
     gevec => veloc(:,:,1)
     call arrays(ivari,'POSTPROCESS',gevec,MESH_ID=imesh)
     
  case ( 'PRESS' )
     !
     ! PRESS: Pressure
     !
     gesca => press(:,1) 
     call arrays(ivari,'POSTPROCESS',gesca,MESH_ID=imesh)

  case ( 'PROJM' )
     !
     ! PROJM_GUS: Projection
     !
     gesca => projm_gus(:,1) 
     call arrays(ivari,'POSTPROCESS',gesca,MESH_ID=imesh)

  case ( 'ANGLE' )
     !
     ! ANGLE_GUS
     !
     call arrays(ivari,'POSTPROCESS',angle_gus,MESH_ID=imesh)

  case ( 'AREAS' )
     !
     ! AREAS: this variable can belong to kermod if it comes from a function
     !
     gesca => areas(:,1)
     call arrays(ivari,'POSTPROCESS',gesca,MESH_ID=imesh,FORCE=.true.)

  case ( 'VEL3D' ) 
     !
     ! VEL3D_GUS
     !
     if( imesh /= 0 ) then
        nullify(vecto)
        call memory_alloca(mem_modul(1:2,modul),'VECTO','gus_outvar',vecto,3_ip,npoin)
        call arrays(arrays_number('VEL3D'),'ALLOCATE',vel3d_gus,3_ip,witness_mesh(imesh) % mesh % npoin)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           if( pelty > 0 ) then
              nn(1:ndime) = (coord(1:ndime,lnods(2,ielem))-coord(1:ndime,lnods(1,ielem)))
              do inode = 1,element_type(pelty) % number_nodes
                 ipoin = lnods(inode,ielem)
                 vecto(1:ndime,ipoin) = vecto(1:ndime,ipoin) + nn(1:ndime)
              end do
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(vecto,'SUM')
        do ipoin = 1,npoin
           call maths_normalize_vector(ndime,vecto(1:ndime,ipoin))           
        end do
        do kpoin = 1,witness_mesh(imesh) % mesh % npoin
           ipoin                = witness_mesh(imesh) % perm(kpoin)
           xx(1:3)              = coord(1:3,ipoin)-witness_mesh(imesh) % mesh % coord(1:3,kpoin)
           radiu                = sqrt(dot_product(xx,xx))
           R                    = sqrt(areas(ipoin,1)/pi)
           U                    = vel1d(ipoin,1)*(1.0_rp-(radiu/R)**2)
           vel3d_gus(1:3,kpoin) = vecto(1:3,ipoin)*U
        end do
        call arrays(ivari,'POSTPROCESS',vel3d_gus,MESH_ID=imesh)
        
        call memory_deallo(mem_modul(1:2,modul),'VECTO','gus_outvar',vecto)
        call arrays(arrays_number('VEL3D'),'DEALLOCATE',vel3d_gus)
     end if
     
  case ( 'PRE3D' )
     !
     ! PRE3D_GUS
     !
     if( imesh /= 0 ) then
        call arrays(arrays_number('PRE3D'),'ALLOCATE',pre3d_gus,witness_mesh(imesh) % mesh % npoin)
        do kpoin = 1,witness_mesh(imesh) % mesh % npoin
           ipoin                = witness_mesh(imesh) % perm(kpoin)
           pre3d_gus(kpoin)     = press(ipoin,1)
        end do
        call arrays(ivari,'POSTPROCESS',pre3d_gus,MESH_ID=imesh)
        call arrays(arrays_number('PRE3D'),'DEALLOCATE',pre3d_gus)
     end if

  case ( 'FLOWR' )
     !
     ! FLOWR: Flow rate
     !
     gesca => flowr(:,1)
     call arrays(ivari,'POSTPROCESS',gesca,MESH_ID=imesh)
     
  case ( 'PROJC' )
     !
     ! PROJC_GUS: Projection
     !
     gesca => projc_gus(:,1) 
     call arrays(ivari,'POSTPROCESS',gesca,MESH_ID=imesh)

  case ( 'LINEL' )
     !
     ! LINEL: Linelets of preconditioner 
     !
     icont = 0
     do ipoin = 1,npoin
        rhsid(ipoin) = 0.0_rp
     end do
     if( INOTEMPTY ) then
        do iline = 1,solve(1) % nline
           icont = icont+1
           do ipoin = solve(1) % lline(iline),solve(1) % lline(iline+1)-1
              kpoin = solve(1) % lrenup(ipoin)
              rhsid(kpoin) = real(icont,rp)
           end do
        end do
     end if
     call arrays(ivari,'POSTPROCESS',rhsid,MESH_ID=imesh,FORCE=.true.)

  case default

     return

  end select

  nullify(gesca,gevec,ger3p)
  
end subroutine gus_outvar
