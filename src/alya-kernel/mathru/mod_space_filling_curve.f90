!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup SFC_Toolbox
!> @{
!> @name    ToolBox for Space filling curve
!> @file    mod_space_filling_curve.f90
!> @author  houzeaux
!> @date    21/02/2018 
!> @brief   ToolBox for SFC
!> @details ToolBox for SFC
!-----------------------------------------------------------------------

module mod_space_filling_curve

  use def_kintyp, only : ip,rp,lg
  implicit none
  real(rp),  parameter :: epsil = epsilon(1.0_rp)

  private 

  public :: space_filling_curve_output 

contains 

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    08/11/2016
  !> @brief   SFC output
  !> @details SFC output. If variable WEIGHT is preent, it is 
  !>          postprocessed. In the case of partitioning woth SFC, WEIGHT
  !>          can be the subdomain the bin belongs to. 
  !>          If message argument 'LEXICAL', it means that WEIGHT is
  !>          given in the lexical order.
  !
  !----------------------------------------------------------------------

  subroutine space_filling_curve_output(ndime_sfc,n,deltx,delty,deltz,weight,wmess)

    use def_kintyp, only : ip,rp
    use mod_maths,  only : maths_sfc_1d_to2d3d_tab
    use mod_maths,  only : maths_mapping_3d_to_1d

    implicit none
    integer(ip),  intent(in)                    :: ndime_sfc         !< Dimension of the SFC
    integer(ip),  intent(in)                    :: n                 !< Size of SFC in each direction
    real(rp),     intent(in)                    :: deltx             !< Bin size in x-direction
    real(rp),     intent(in)                    :: delty             !< Bin size in y-direction
    real(rp),     intent(in)                    :: deltz             !< Bin size in z-direction
    integer(ip),  intent(in), pointer, optional :: weight(:)         !< Weight of elements
    character(*), intent(in),          optional :: wmess
    integer(ip)                                 :: d,x,y,z,ii,ii0
    integer(ip)                                 :: jj,ntot,ntot1,n1
    integer(ip)                                 :: imate,dlex
    real(rp)                                    :: xcoor(3)
    real(rp),     pointer                       :: coord_sfc(:,:)
    integer(ip),  pointer                       :: weight_loc(:)

    dlex  = 0
    n1    = n+1
    ntot  = n**ndime_sfc
    ntot1 = n1**ndime_sfc

    if( present(weight) ) then
       if( size(weight,kind=ip) /= ntot ) call runend('SPACE_FILLING_CURVE_OUTPUT: WRONG SIZE')
    end if

    if( present(weight) .and. present(wmess) ) then
       if( trim(wmess) == 'LEXICAL' ) then
          allocate( weight_loc(ntot) )
          if( ndime_sfc == 2 ) then
             do d = 1,ntot         
                call maths_sfc_1d_to2d3d_tab(n,d,x,y)
                dlex = maths_mapping_3d_to_1d(n,n,1_ip,x,y,1_ip)
                weight_loc(d) = weight(dlex)
             end do
          else
             do d = 1,ntot         
                call maths_sfc_1d_to2d3d_tab(n,d,x,y,z)
                dlex = maths_mapping_3d_to_1d(n,n,n,x,y,z)
                weight_loc(d) = weight(dlex)
             end do
          end if
       end if
    end if
    if( dlex == 0 .and. present(weight) ) weight_loc => weight

    nullify(coord_sfc)
    allocate(coord_sfc(ndime_sfc,n1**ndime_sfc))

    open(unit=90,file='space.post.msh',status='unknown')
    open(unit=91,file='space.post.res',status='unknown')
    !
    ! GiD header
    !
    if( ndime_sfc == 2 ) then
       write(90,'(a)') 'MESH SFC dimension 2 Elemtype Quadrilateral Nnode 4'       
    else
       write(90,'(a)') 'MESH SFC dimension 3 Elemtype Hexahedra Nnode 8'       
    end if
    !
    ! Bin coordinates
    !
    write(90,'(a)') 'coordinates'  

    ii = 0

    if( ndime_sfc == 2 ) then
       do y = 1,n1
          do x = 1,n1
             ii = ii + 1
             coord_sfc(1,ii) = real(x-1,rp) * deltx
             coord_sfc(2,ii) = real(y-1,rp) * delty
             write(90,1) ii,coord_sfc(1:2,ii)
          end do
       end do
    else
       do z = 1,n1
          do y = 1,n1
             do x = 1,n1
                ii = ii + 1
                coord_sfc(1,ii) = real(x-1,rp) * deltx
                coord_sfc(2,ii) = real(y-1,rp) * delty
                coord_sfc(3,ii) = real(z-1,rp) * deltz
                write(90,1) ii,coord_sfc(1:3,ii)
             end do
          end do
       end do
    end if
    !
    ! Bin center coordinates
    !
    jj = ii 

    if( ndime_sfc == 2 ) then

       do d = 1,ntot
          call maths_sfc_1d_to2d3d_tab(n,d,x,y)
          ii = x + n1*( y-1 )
          jj = jj + 1
          xcoor(1:2) = 0.25_rp * ( coord_sfc(1:2,ii)     + coord_sfc(1:2,ii+1) &
               &                 + coord_sfc(1:2,ii+n+2) + coord_sfc(1:2,ii+n+1) )
          write(90,1) jj,xcoor(1:2)
       end do

    else

       do d = 1,ntot
          call maths_sfc_1d_to2d3d_tab(n,d,x,y,z)
          ii  = x  + n1*( y-1 + (z-1)*n1 )
          ii0 = ii + n1*n1
          jj  = jj + 1
          xcoor(1:3) = 0.125_rp * ( coord_sfc(1:3,ii)      + coord_sfc(1:3,ii+1)    &
               &                  + coord_sfc(1:3,ii+n+2)  + coord_sfc(1:3,ii+n+1)  &
               &                  + coord_sfc(1:3,ii0)     + coord_sfc(1:3,ii0+1)   &
               &                  + coord_sfc(1:3,ii0+n+2) + coord_sfc(1:3,ii0+n+1) )
          write(90,1) jj,xcoor(1:3)
       end do

    end if

    write(90,'(a)') 'end coordinates'  
    !
    ! Bin elements
    !
    write(90,'(a)') 'elements'  

    imate = 1

    if( ndime_sfc == 2 ) then

       do d = 1,ntot
          if( present(weight) ) imate = min(1_ip,weight_loc(d))
          call maths_sfc_1d_to2d3d_tab(n,d,x,y)
          ii = x + n1*( y-1 )
          write(90,2) d,ii,ii+1,ii+n+2,ii+n+1,imate
       end do

    else

       do d = 1,ntot
          if( present(weight) ) imate = min(1_ip,weight_loc(d))
          call maths_sfc_1d_to2d3d_tab(n,d,x,y,z)
          ii  = x  + n1*( y-1 + (z-1)*n1 )
          ii0 = ii + n1*n1
          write(90,2) d,ii,ii+1,ii+n+2,ii+n+1,ii0,ii0+1,ii0+n+2,ii0+n+1,imate
       end do

    end if

    write(90,'(a)') 'end elements'  
    !
    ! Bar elements between consecutive bins of the SFC
    !
    if( ndime_sfc == 2 ) then
       write(90,'(a)') 'MESH SFC dimension 2 Elemtype Linear Nnode 2'       
    else
       write(90,'(a)') 'MESH SFC dimension 3 Elemtype Linear Nnode 2'       
    end if
    write(90,'(a)') 'elements'

    do d = 1,ntot-1
       write(90,2) d+ntot,ntot1+d,ntot1+d+1
    end do

    write(90,'(a)') 'end elements'  
    !
    ! Results
    !
    write(91,*) 'GiD Post Results File 1.0'
    write(91,*) ' '
    !
    ! 1D coordinate
    !
    if( ndime_sfc == 2 ) then
       write(91,*) 'GaussPoints GP_QUAD4 Elemtype Quadrilateral'
    else
       write(91,*) 'GaussPoints GP_QUAD4 Elemtype Hexahedra'
    end if
    write(91,*) 'Number of Gauss Points: 1'
    write(91,*) 'Natural Coordinates: Internal'
    write(91,*) 'End GaussPoints'

    write(91,*) 'Result 1D_COORD 1D_COORD 0 Scalar OnGaussPoints GP_QUAD4'
    write(91,*) 'ComponentNames 1D_COORD'
    write(91,*) 'values'
    do d = 1,ntot
       write(91,2) d,d
    end do
    write(91,*) 'end values' 

    write(91,*) 'GaussPoints GP_LINEAR Elemtype Linear'
    write(91,*) 'Number of Gauss Points: 1'
    write(91,*) 'Natural Coordinates: Internal'
    write(91,*) 'End GaussPoints'

    write(91,*) 'Result 1D_COORD 1D_COORD_NODE 0 Scalar OnGaussPoints GP_LINEAR'
    write(91,*) 'ComponentNames 1D_COORD_NODE'
    write(91,*) 'values'
    do d = 1,ntot-1
       write(91,2) d+ntot,d
    end do
    write(91,*) 'end values' 
    !
    ! Weight
    !
    if( present(weight) ) then
       write(91,*) 'Result WEIGHT WEIGHT 0 Scalar OnGaussPoints GP_QUAD4'
       write(91,*) 'ComponentNames WEIGHT'
       write(91,*) 'values'
       do d = 1,ntot
          write(91,2) d,weight_loc(d)
       end do
       write(91,*) 'end values' 
    end if
    if( present(weight) ) then
       write(91,*) 'Result WEIGHT_NODE WEIGHT_NODE 0 Scalar OnNodes'
       write(91,*) 'ComponentNames WEIGHT_NODE'
       write(91,*) 'values'
       do d = 1,ntot
          write(91,2) d+ntot,weight_loc(d)
       end do
       write(91,*) 'end values' 
    end if
    !
    ! Deallocate
    !
    deallocate(coord_sfc)
    if( dlex /= 0 ) deallocate( weight_loc ) 

    close(90)
    close(91)

1   format(i9,3(1x,e12.6))
2   format(10(1x,i9))

  end subroutine space_filling_curve_output

end module mod_space_filling_curve
!> @}
