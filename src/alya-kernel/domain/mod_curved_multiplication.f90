!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



MODULE mod_curved_multiplication
  use def_kintyp, only: ip, rp, lg
  use def_elmtyp, only: TRI03, QUA04, HEX08

  IMPLICIT NONE
  SAVE

  !-Module variables
  integer(ip), parameter, private :: maxOrder = 2_ip
  integer(ip), parameter, private :: maxElementType = HEX08

  type shapeFunctionType
    real(rp), allocatable :: values(:,:)
  end type shapeFunctionType

  type(shapeFunctionType), private, target :: shapeFunctions(maxElementType,maxOrder)

  !-Module functions
  public  :: curve_subdivided_element
  private :: curve_subdivided_tri03, curve_subdivided_qua04, curve_subdivided_hex08, &
             shapeFunctions_tri, shapeFunctions_qua, shapeFunctions_hex, &
             getMonomials_trip2, FINDInv

CONTAINS

subroutine curve_subdivided_element(elementType,numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
  curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
  ! -----------------------------------------------------------------------
  !****f* domain/curve_subdivided_element
  ! NAME
  !    curve_subdivided_element
  ! DESCRIPTION
  !    This subroutine curves a subdivided element according to a curved geometry field
  ! OUTPUT
  ! USED BY
  !    submsh
  !***
  ! -----------------------------------------------------------------------
  implicit none
  integer(ip),  intent(in)  :: elementType, numNewCoordinates, numNewElements, orderGeometry, numNodes, numDimPhys
  real(rp),     intent(in)  :: curved_geometry(numNodes,numDimPhys)
  real(rp),     intent(out) :: subdivided_coordinates(numNewCoordinates,numDimPhys)
  real(rp),     intent(out) :: subdivided_curved_geometry(numNodes,numDimPhys,numNewElements)

  select case(elementType)
    case(TRI03)
      call curve_subdivided_tri03(numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
        curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
    case(QUA04)
      call curve_subdivided_qua04(numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
        curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
    case(HEX08)
      call curve_subdivided_hex08(numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
        curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
    case default
      call runend('CURVED_MULTIPLICATION: Not implemented curved subdivision for this kind of element')
  end select

return
end subroutine curve_subdivided_element
!
!
!
subroutine curve_subdivided_tri03(numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
  curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
  implicit none

  integer(ip),  intent(in)  :: numNewCoordinates, numNewElements, orderGeometry, numNodes, numDimPhys
  real(rp),     intent(in)  :: curved_geometry(numNodes,numDimPhys)
  real(rp),     intent(out) :: subdivided_coordinates(numNewCoordinates,numDimPhys)
  real(rp),     intent(out) :: subdivided_curved_geometry(numNodes,numDimPhys,numNewElements)

  integer(ip), parameter :: referenceDimension = 2_ip
  integer(ip), parameter :: numEvalPoints = 9_ip !=numNewCoordinates
  !
  real(rp)    :: refPoints(referenceDimension,numNodes)
  real(rp)    :: evalPoints(referenceDimension,numEvalPoints)
  real(rp), pointer   :: shapeF(:,:)!(numEvalPoints,numNodes)

  real(rp)    :: newCurvedGeometryPoints(numEvalPoints,numDimPhys)!size(curved_geometry,2,KIND=ip))
  real(rp)    :: newCurvedGeometry(numNodes+numEvalPoints,numDimPhys)!size(curved_geometry,2,KIND=ip))

  integer(ip) :: newConnectivity(numNodes,numNewElements)

  integer(ip) :: idime,inode,ielem

  ! 0. Initialize variables and check consistency of the geometry data
  if(numNewCoordinates.ne.3) then
    call runend('CURVED_MULTIPLICATION: Check consistency in the number of new created coordinates.')
  end if
  if(numNewElements.ne.4) then
    call runend('CURVED_MULTIPLICATION: Not correct number of subdivided elements')
  end if
  if(orderGeometry.gt.2) then
    call runend('CURVED_MULTIPLICATION: Curving only allowed for order 2.')
  end if

  ! A. New curved coordinates: direct from the high-order field
  subdivided_coordinates(:,:) = curved_geometry(4:6,:)

  ! B. Multiply curved field for the subelements
  if( .not.allocated(shapeFunctions(TRI03,orderGeometry)%values) ) then
    ! B.0. Define reference element
    refPoints(1,:) = (/ -1.0_rp,  1.0_rp, -1.0_rp,  0.0_rp,  0.0_rp, -1.0_rp /)
    refPoints(2,:) = (/ -1.0_rp, -1.0_rp,  1.0_rp, -1.0_rp,  0.0_rp,  0.0_rp /)
    ! B.1. Define points in the reference space where we want new curvature
    evalPoints(1,:) = (/ -0.5_rp,  0.5_rp,  0.5_rp, -0.5_rp, -1.0_rp, -1.0_rp, -0.5_rp,  0.0_rp, -0.5_rp/) !xi
    evalPoints(2,:) = (/ -1.0_rp, -1.0_rp, -0.5_rp,  0.5_rp,  0.5_rp, -0.5_rp, -0.5_rp, -0.5_rp,  0.0_rp/) !eta
    ! B.2. Evaluate shape functions in the points where we want curvature
    allocate( shapeFunctions(TRI03,orderGeometry)%values(numEvalPoints,numNodes) )
    shapeFunctions(TRI03,orderGeometry)%values = &
      shapeFunctions_tri(orderGeometry,numNodes,numEvalPoints,refPoints,evalPoints)
  end if

  shapeF => shapeFunctions(TRI03,orderGeometry)%values

  ! B.3. Compute new curved points using isoparametric mapping
  newCurvedGeometryPoints = 0_ip
  do inode=1,numNodes
    do idime=1,numDimPhys
      newCurvedGeometryPoints(:,idime) = newCurvedGeometryPoints(:,idime)+shapeF(:,inode)*curved_geometry(inode,idime)
    end do
  end do
  newCurvedGeometry(1:6,:)  = curved_geometry(1:6,:)
  newCurvedGeometry(7:15,:) = newCurvedGeometryPoints

  ! B.4. Construct new curved geometry
  newConnectivity(:,1) = (/ 1, 4, 6,  7, 13, 12 /)
  newConnectivity(:,2) = (/ 4, 5, 6, 14, 15, 13 /)
  newConnectivity(:,3) = (/ 4, 2, 5,  8,  9, 14 /)
  newConnectivity(:,4) = (/ 6, 5, 3, 15, 10, 11 /)
  do ielem=1,numNewElements
    subdivided_curved_geometry(:,:,ielem) = newCurvedGeometry( newConnectivity(:,ielem) ,:)
  end do
return
end subroutine curve_subdivided_tri03
!
!
!
subroutine curve_subdivided_qua04(numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
  curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
  implicit none

  integer(ip),  intent(in)  :: numNewCoordinates, numNewElements, orderGeometry, numNodes, numDimPhys
  real(rp),     intent(in)  :: curved_geometry(numNodes,numDimPhys)
  real(rp),     intent(out) :: subdivided_coordinates(numNewCoordinates,numDimPhys)
  real(rp),     intent(out) :: subdivided_curved_geometry(numNodes,numDimPhys,numNewElements)

  integer(ip), parameter :: referenceDimension = 2_ip
  integer(ip), parameter :: numEvalPoints = 16_ip !=numNewCoordinates
  !
  real(rp)    :: refPoints(referenceDimension,numNodes)
  real(rp)    :: evalPoints(referenceDimension,numEvalPoints)
  real(rp), pointer   :: shapeF(:,:)!(numEvalPoints,numNodes)

  real(rp)    :: newCurvedGeometryPoints(numEvalPoints,numDimPhys)!size(curved_geometry,2,KIND=ip))
  real(rp)    :: newCurvedGeometry(numNodes+numEvalPoints,numDimPhys)!size(curved_geometry,2,KIND=ip))

  integer(ip) :: newConnectivity(numNodes,numNewElements)

  integer(ip) :: idime,inode,ielem

  real(rp), parameter :: h =  0.5_rp

  ! 0. Initialize variables and check consistency of the geometry data
  if(numNewCoordinates.ne.5) then
    call runend('CURVED_MULTIPLICATION: Check consistency in the number of new created coordinates.')
  end if
  if(numNewElements.ne.4) then
    call runend('CURVED_MULTIPLICATION: Not correct number of subdivided elements')
  end if
  if(orderGeometry.gt.2) then
    call runend('CURVED_MULTIPLICATION: Curving only allowed for order 2.')
  end if

  ! A. New curved coordinates: direct from the high-order field
  subdivided_coordinates(:,:) = curved_geometry(5:9,:)

  ! B. Multiply curved field for the subelements
  if( .not.allocated(shapeFunctions(QUA04,orderGeometry)%values) ) then
    ! B.0. Define reference element
    refPoints(1,:) = (/ -1.0_rp,  1.0_rp,  1.0_rp, -1.0_rp,  0.0_rp,  1.0_rp,  0.0_rp, -1.0_rp,  0.0_rp /)
    refPoints(2,:) = (/ -1.0_rp, -1.0_rp,  1.0_rp,  1.0_rp, -1.0_rp,  0.0_rp,  1.0_rp,  0.0_rp,  0.0_rp /)
    ! B.1. Define points in the reference space where we want new curvature
    evalPoints(1,:) =  (/-h,h,1.0_rp,1.0_rp,h,-h,-1.0_rp,-1.0_rp,-h,0.0_rp,h,-h,h,-h,0.0_rp,h /) !xi
    evalPoints(2,:) =  (/-1.0_rp,-1.0_rp,-h,h,1.0_rp,1.0_rp,h,-h,-h,-h,-h,0.0_rp,0.0_rp,h,h,h /) !eta
    ! B.2. Evaluate shape functions in the points where we want curvature
    allocate( shapeFunctions(QUA04,orderGeometry)%values(numEvalPoints,numNodes) )
    shapeFunctions(QUA04,orderGeometry)%values = &
      shapeFunctions_qua(orderGeometry,numNodes,numEvalPoints,refPoints,evalPoints)
  end if
  shapeF => shapeFunctions(QUA04,orderGeometry)%values

  ! B.3. Compute new curved points using isoparametric mapping
  newCurvedGeometryPoints = 0_ip
  do inode=1,numNodes
    do idime=1,numDimPhys
      newCurvedGeometryPoints(:,idime) = newCurvedGeometryPoints(:,idime)+shapeF(:,inode)*curved_geometry(inode,idime)
    end do
  end do
  newCurvedGeometry(1:9,:)   = curved_geometry(1:9,:)
  newCurvedGeometry(10:25,:) = newCurvedGeometryPoints

  ! B.4. Construct new curved geometry
  newConnectivity(:,1) = (/ 1,5,9,8,10,19,21,17,18 /)
  newConnectivity(:,2) = (/ 5,2,6,9,11,12,22,19,20 /)
  newConnectivity(:,3) = (/ 9,6,3,7,22,13,14,24,25 /)
  newConnectivity(:,4) = (/ 8,9,7,4,21,24,15,16,23 /)
  do ielem=1,numNewElements
    subdivided_curved_geometry(:,:,ielem) = newCurvedGeometry( newConnectivity(:,ielem) ,:)
  end do
return
end subroutine curve_subdivided_qua04
!
!
!
subroutine curve_subdivided_hex08(numNewCoordinates,numNewElements,orderGeometry,numNodes,numDimPhys,&
  curved_geometry,subdivided_coordinates,subdivided_curved_geometry)
  implicit none

  integer(ip),  intent(in)  :: numNewCoordinates, numNewElements, orderGeometry, numNodes, numDimPhys
  real(rp),     intent(in)  :: curved_geometry(numNodes,numDimPhys)
  real(rp),     intent(out) :: subdivided_coordinates(numNewCoordinates,numDimPhys)
  real(rp),     intent(out) :: subdivided_curved_geometry(numNodes,numDimPhys,numNewElements)

  integer(ip), parameter :: referenceDimension = 3
  integer(ip), parameter :: numEvalPoints = 27*8!numNodes*numNewElements!19 !=numNewCoordinates
  !
  real(rp)    :: refPoints(referenceDimension,numNodes)
  real(rp)    :: evalPoints(referenceDimension,numEvalPoints)
  real(rp), pointer   :: shapeF(:,:)!(numEvalPoints,numNodes)
  real(rp)    :: xinew(numNodes),etanew(numNodes),zetanew(numNodes)

  real(rp)    :: newCurvedGeometryPoints(numEvalPoints,numDimPhys)!size(curved_geometry,2,KIND=ip))
  !real(rp)    :: newCurvedGeometry(numNodes+numEvalPoints,numDimPhys)!size(curved_geometry,2,KIND=ip))

  !integer(ip) :: newConnectivity(numNodes,numNewElements)

  integer(ip) :: idime,inode,ielem

  real(rp), parameter :: h =  0.5_rp

  ! 0. Initialize variables and check consistency of the geometry data
  if(numNewElements.ne.8_ip) then
    call runend('CURVED_MULTIPLICATION: Not correct number of subdivided elements')
  end if
  if(orderGeometry.gt.2_ip) then
    call runend('CURVED_MULTIPLICATION: Curving only allowed for order 2.')
  end if

  ! A. New curved coordinates: direct from the high-order field
  subdivided_coordinates(:,:) = curved_geometry(9:27,:)

  ! B. Multiply curved field for the subelements
  if( .not.allocated(shapeFunctions(HEX08,orderGeometry)%values) ) then
    ! B.0. Define reference element
    refPoints(1,:) = (/ -1.0_rp,  1.0_rp,  1.0_rp, -1.0_rp,  -1.0_rp,  1.0_rp, 1.0_rp, -1.0_rp, 0.0_rp,  1.0_rp, 0.0_rp, -1.0_rp,&
                         0.0_rp,  1.0_rp,  0.0_rp, -1.0_rp,  -1.0_rp,  1.0_rp, 1.0_rp, -1.0_rp, 0.0_rp,  1.0_rp, 0.0_rp, -1.0_rp,&
                         0.0_rp,  0.0_rp,  0.0_rp  /)
    refPoints(2,:) = (/ -1.0_rp, -1.0_rp,  1.0_rp,  1.0_rp,  -1.0_rp, -1.0_rp, 1.0_rp,  1.0_rp,-1.0_rp,  0.0_rp, 1.0_rp,  0.0_rp,&
                        -1.0_rp,  0.0_rp,  1.0_rp,  0.0_rp,  -1.0_rp, -1.0_rp, 1.0_rp,  1.0_rp, 0.0_rp,  0.0_rp, 0.0_rp,  0.0_rp,&
                        -1.0_rp,  1.0_rp,  0.0_rp  /)
    refPoints(3,:) = (/ -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp,   1.0_rp,  1.0_rp, 1.0_rp,  1.0_rp,-1.0_rp, -1.0_rp,-1.0_rp, -1.0_rp,&
                         1.0_rp,  1.0_rp,  1.0_rp,  1.0_rp,   0.0_rp,  0.0_rp, 0.0_rp , 0.0_rp,-1.0_rp,  0.0_rp, 1.0_rp,  0.0_rp,&
                         0.0_rp,  0.0_rp,  0.0_rp  /)

    ! B.1. Define evaluation points on reference element to get new curvature field
    xinew   = refPoints(1,:)/2.0_rp
    etanew  = refPoints(2,:)/2.0_rp
    zetanew = refPoints(3,:)/2.0_rp

    evalPoints(1,(0*numNodes+1):(1*numNodes)) = xinew   - h
    evalPoints(2,(0*numNodes+1):(1*numNodes)) = etanew  - h
    evalPoints(3,(0*numNodes+1):(1*numNodes)) = zetanew - h

    evalPoints(1,(1*numNodes+1):(2*numNodes)) = xinew   + h
    evalPoints(2,(1*numNodes+1):(2*numNodes)) = etanew  - h
    evalPoints(3,(1*numNodes+1):(2*numNodes)) = zetanew - h

    evalPoints(1,(2*numNodes+1):(3*numNodes)) = xinew   + h
    evalPoints(2,(2*numNodes+1):(3*numNodes)) = etanew  + h
    evalPoints(3,(2*numNodes+1):(3*numNodes)) = zetanew - h

    evalPoints(1,(3*numNodes+1):(4*numNodes)) = xinew   - h
    evalPoints(2,(3*numNodes+1):(4*numNodes)) = etanew  + h
    evalPoints(3,(3*numNodes+1):(4*numNodes)) = zetanew - h

    evalPoints(1,(4*numNodes+1):(5*numNodes)) = evalPoints(1,(0*numNodes+1):(1*numNodes))
    evalPoints(2,(4*numNodes+1):(5*numNodes)) = evalPoints(2,(0*numNodes+1):(1*numNodes))
    evalPoints(3,(4*numNodes+1):(5*numNodes)) = zetanew + h

    evalPoints(1,(5*numNodes+1):(6*numNodes)) = evalPoints(1,(1*numNodes+1):(2*numNodes))
    evalPoints(2,(5*numNodes+1):(6*numNodes)) = evalPoints(2,(1*numNodes+1):(2*numNodes))
    evalPoints(3,(5*numNodes+1):(6*numNodes)) = zetanew + h

    evalPoints(1,(6*numNodes+1):(7*numNodes)) = evalPoints(1,(2*numNodes+1):(3*numNodes))
    evalPoints(2,(6*numNodes+1):(7*numNodes)) = evalPoints(2,(2*numNodes+1):(3*numNodes))
    evalPoints(3,(6*numNodes+1):(7*numNodes)) = zetanew + h

    evalPoints(1,(7*numNodes+1):(8*numNodes)) = evalPoints(1,(3*numNodes+1):(4*numNodes))
    evalPoints(2,(7*numNodes+1):(8*numNodes)) = evalPoints(2,(3*numNodes+1):(4*numNodes))
    evalPoints(3,(7*numNodes+1):(8*numNodes)) = zetanew + h

    allocate( shapeFunctions(HEX08,orderGeometry)%values(numEvalPoints,numNodes) )
    shapeFunctions(HEX08,orderGeometry)%values = &
      shapeFunctions_hex(orderGeometry,numNodes,numEvalPoints,refPoints,evalPoints)
  end if
  shapeF => shapeFunctions(HEX08,orderGeometry)%values

  ! B.3. Compute new curved points using isoparametric mapping
  newCurvedGeometryPoints = 0.0_rp
  do inode=1,numNodes
    do idime=1,numDimPhys
      newCurvedGeometryPoints(:,idime) = newCurvedGeometryPoints(:,idime)+shapeF(:,inode)*curved_geometry(inode,idime)
    end do
  end do
!
  ! B.4. Construct new curved geometry
!   newConnectivity(1:8,1) = (/  1,9,21,12,17,25,27,24 /)
!   newConnectivity(1:8,2) = (/  9,2,10,21,25,18,22,27 /)
!   newConnectivity(1:8,3) = (/ 21,10,3,11,27,22,19,26 /)
!   newConnectivity(1:8,4) = (/ 12,21,11,4,24,27,26,20 /)
!   newConnectivity(1:8,5) = (/ 17,25,27,24,5,13,23,16 /)
!   newConnectivity(1:8,6) = (/ 25,18,22,27,13,6,14,23 /)
!   newConnectivity(1:8,7) = (/ 27,22,19,26,23,14,7,15 /)
!   newConnectivity(1:8,8) = (/ 24,27,26,20,16,23,15,8 /)
  do ielem=1,numNewElements
    subdivided_curved_geometry(:,:,ielem) = newCurvedGeometryPoints( ((ielem-1)*numNodes+1):(ielem*numNodes) ,:)
    !The following would save to calculate 8*8 points.. but it is not relevant compared with the 6*(81-8) points that would have to be computed,
    !taking into accont that now the code is really readable, whereas the other way it could not be understood
!     subdivided_curved_geometry(1:8,:,ielem) = curved_geometry(newConnectivity(1:8,ielem),:)
  end do
return
end subroutine curve_subdivided_hex08
!
!
!
function shapeFunctions_tri(order,numNodes,numEvalPoints,refPoints,evalPoints) result(shapeF)
  implicit none

  integer(ip), intent(in) :: order,numNodes,numEvalPoints
  real(rp) ,intent(in)  :: refPoints(2,numNodes), evalPoints(2,numEvalPoints)
  real(rp) :: shapeF(numEvalPoints,numNodes)

  real(rp) :: SHref(numNodes,numNodes)
  real(rp) :: SHrefInv(numNodes,numNodes)
  real(rp) :: SHeval(numEvalPoints,numNodes)

  integer(ip) :: errorFlag

  if(order.gt.2_ip) then
    call runend('CURVED_MULTIPLICATION: Curving only allowed for order 2.')
  end if

  SHref  = getMonomials_trip2(numNodes,transpose(refPoints))
  SHeval = getMonomials_trip2(numEvalPoints,transpose(evalPoints))

  call FINDInv(transpose(SHref), SHrefInv, numNodes, errorflag)
  if(errorFlag.eq.1) then
    call runend('CURVED_MULTIPLICATION: Cannot generate high-order shape functions for mesh multiplication')
  end if

  shapeF = transpose( matmul(SHrefInv,transpose(SHeval)) )

  return
end function shapeFunctions_tri
!
!
!
function shapeFunctions_qua(order,numNodes,numEvalPoints,refPoints,evalPoints) result(shapeF)
  implicit none

  integer(ip), intent(in) :: order,numNodes,numEvalPoints
  real(rp) ,intent(in)  :: refPoints(2,numNodes), evalPoints(2,numEvalPoints)
  real(rp) :: shapeF(numEvalPoints,numNodes)

  real(rp) :: xi(numEvalPoints),eta(numEvalPoints)

  if(order.gt.2) then
    call runend('CURVED_MULTIPLICATION: Curving only allowed for order 2.')
  end if

  if((refPoints(1,1).ne.(-1.0_rp)).or.(refPoints(1,3).ne.(1.0_rp)).or.(refPoints(2,1).ne.(-1.0_rp)).or.(refPoints(2,3).ne.(1.0_rp))) then
    call runend('CURVED_MULTIPLICATION: Quad shape functions only implemented for reference element [-1,1]x[-1,1]')
  end if

  xi = evalPoints(1,:)
  eta = evalPoints(2,:)

  shapeF(:,1) = xi*(xi-1.0_rp)*eta*(eta-1.0_rp)*0.25_rp
  shapeF(:,2) = xi*(xi+1.0_rp)*eta*(eta-1.0_rp)*0.25_rp
  shapeF(:,3) = xi*(xi+1.0_rp)*eta*(eta+1.0_rp)*0.25_rp
  shapeF(:,4) = xi*(xi-1.0_rp)*eta*(eta+1.0_rp)*0.25_rp
  shapeF(:,5) = (1.0_rp-xi*xi)*eta*(eta-1.0_rp)*0.5_rp
  shapeF(:,6) = xi*(xi+1.0_rp)*(1.0_rp-eta*eta)*0.5_rp
  shapeF(:,7) = (1.0_rp-xi*xi)*eta*(eta+1.0_rp)*0.5_rp
  shapeF(:,8) = xi*(xi-1.0_rp)*(1.0_rp-eta*eta)*0.5_rp
  shapeF(:,9) = (1.0_rp-xi*xi)*(1.0_rp-eta*eta)

  return
end function shapeFunctions_qua
!
!
!
function shapeFunctions_hex(order,numNodes,numEvalPoints,refPoints,evalPoints) result(shapeF)
  implicit none

  integer(ip), intent(in) :: order,numNodes,numEvalPoints
  real(rp) ,intent(in)  :: refPoints(3,numNodes), evalPoints(3,numEvalPoints)
  real(rp) :: shapeF(numEvalPoints,numNodes)

  real(rp) :: xi(numEvalPoints),eta(numEvalPoints),zeta(numEvalPoints)

  if(order.gt.2) then
    call runend('CURVED_MULTIPLICATION: Curving only allowed for order 2.')
  end if

  if((refPoints(1,1).ne.(-1.0_rp)).or.(refPoints(1,7).ne.(1.0_rp)).or.(refPoints(2,1).ne.(-1.0_rp)).or.(refPoints(2,7).ne.(1.0_rp))) then
    call runend('CURVED_MULTIPLICATION: Quad shape functions only implemented for reference element [-1,1]^3')
  end if

  xi   = evalPoints(1,:)
  eta  = evalPoints(2,:)
  zeta = evalPoints(3,:)

  shapeF(:,1)  = xi*(xi-1.0_rp)*eta*(eta-1.0_rp)*0.25_rp *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,2)  = xi*(xi+1.0_rp)*eta*(eta-1.0_rp)*0.25_rp *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,3)  = xi*(xi+1.0_rp)*eta*(eta+1.0_rp)*0.25_rp *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,4)  = xi*(xi-1.0_rp)*eta*(eta+1.0_rp)*0.25_rp *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,9)  = (1.0_rp-xi*xi)*eta*(eta-1.0_rp)*0.5_rp  *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,10) = xi*(xi+1.0_rp)*(1.0_rp-eta*eta)*0.5_rp  *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,11) = (1.0_rp-xi*xi)*eta*(eta+1.0_rp)*0.5_rp  *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,12) = xi*(xi-1.0_rp)*(1.0_rp-eta*eta)*0.5_rp  *zeta*(zeta-1.0_rp)/2.0_rp
  shapeF(:,21) = (1.0_rp-xi*xi)*(1.0_rp-eta*eta)      *zeta*(zeta-1.0_rp)/2.0_rp

  shapeF(:,17) = xi*(xi-1.0_rp)*eta*(eta-1.0_rp)*0.25_rp *(1.0_rp-zeta*zeta)
  shapeF(:,18) = xi*(xi+1.0_rp)*eta*(eta-1.0_rp)*0.25_rp *(1.0_rp-zeta*zeta)
  shapeF(:,19) = xi*(xi+1.0_rp)*eta*(eta+1.0_rp)*0.25_rp *(1.0_rp-zeta*zeta)
  shapeF(:,20) = xi*(xi-1.0_rp)*eta*(eta+1.0_rp)*0.25_rp *(1.0_rp-zeta*zeta)
  shapeF(:,25) = (1.0_rp-xi*xi)*eta*(eta-1.0_rp)*0.5_rp  *(1.0_rp-zeta*zeta)
  shapeF(:,22) = xi*(xi+1.0_rp)*(1.0_rp-eta*eta)*0.5_rp  *(1.0_rp-zeta*zeta)
  shapeF(:,26) = (1.0_rp-xi*xi)*eta*(eta+1.0_rp)*0.5_rp  *(1.0_rp-zeta*zeta)
  shapeF(:,24) = xi*(xi-1.0_rp)*(1.0_rp-eta*eta)*0.5_rp  *(1.0_rp-zeta*zeta)
  shapeF(:,27) = (1.0_rp-xi*xi)*(1.0_rp-eta*eta)      *(1.0_rp-zeta*zeta)

  shapeF(:,5)  = xi*(xi-1.0_rp)*eta*(eta-1.0_rp)*0.25_rp *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,6)  = xi*(xi+1.0_rp)*eta*(eta-1.0_rp)*0.25_rp *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,7)  = xi*(xi+1.0_rp)*eta*(eta+1.0_rp)*0.25_rp *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,8)  = xi*(xi-1.0_rp)*eta*(eta+1.0_rp)*0.25_rp *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,13) = (1.0_rp-xi*xi)*eta*(eta-1.0_rp)*0.5_rp  *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,14) = xi*(xi+1.0_rp)*(1.0_rp-eta*eta)*0.5_rp  *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,15) = (1.0_rp-xi*xi)*eta*(eta+1.0_rp)*0.5_rp  *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,16) = xi*(xi-1.0_rp)*(1.0_rp-eta*eta)*0.5_rp  *zeta*(zeta+1.0_rp)/2.0_rp
  shapeF(:,23) = (1.0_rp-xi*xi)*(1.0_rp-eta*eta)      *zeta*(zeta+1.0_rp)/2.0_rp

  return
end function shapeFunctions_hex
!
!
!
function getMonomials_trip2(numPoints,x) result(M)
  implicit none
  integer(ip) ,intent(in)  :: numPoints
  real(rp)    ,intent(in)  :: x(numPoints,6)

  real(rp) :: M(size(x,1,KIND=ip),6)

  M(:,1) = 1.0_rp
  M(:,2) = x(:,1)
  M(:,3) = x(:,2)
  M(:,4) = x(:,1)*x(:,2)
  M(:,5) = x(:,1)*x(:,1)
  M(:,6) = x(:,2)*x(:,2)

  return
end function
!
!
!
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
  !Subroutine to find the inverse of a square matrix
  !Author : Louisda16th a.k.a Ashwith J. Rego
  !Reference : Algorithm has been well explained in:
  !http://math.uww.edu/~mcfarlat/inverse.htm
  !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
  IMPLICIT NONE
  !Declarations
  INTEGER(ip), INTENT(IN) :: n
  INTEGER(ip), INTENT(OUT) :: errorflag !Return error status. -1 for error, 0 for normal
  REAL(rp), INTENT(IN), DIMENSION(n,n) :: matrix !Input matrix
  REAL(rp), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

  LOGICAL :: FLAG = .TRUE.
  INTEGER(ip) :: i, j, k
  REAL(rp) :: m
  REAL(rp), DIMENSION(n,2*n) :: augmatrix !augmented matrix

  !Augment input matrix with an identity matrix
  DO i = 1, n
  DO j = 1, 2*n
  IF (j <= n ) THEN
  augmatrix(i,j) = matrix(i,j)
  ELSE IF ((i+n) == j) THEN
  augmatrix(i,j) = 1
  Else
  augmatrix(i,j) = 0
  ENDIF
  END DO
  END DO

  !Reduce augmented matrix to upper traingular form
  DO k =1, n-1
  IF (augmatrix(k,k) == 0) THEN
  FLAG = .FALSE.
  DO i = k+1, n
  IF (augmatrix(i,k) /= 0) THEN
  DO j = 1,2*n
  augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
  END DO
  FLAG = .TRUE.
  EXIT
  ENDIF
  IF (FLAG .EQV. .FALSE.) THEN
  PRINT*, "Matrix is non - invertible"
  inverse = 0
  errorflag = -1
  return
  ENDIF
  END DO
  ENDIF
  DO j = k+1, n
  m = augmatrix(j,k)/augmatrix(k,k)
  DO i = k, 2*n
  augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
  END DO
  END DO
  END DO

  !Test for invertibility
  DO i = 1, n
  IF (augmatrix(i,i) == 0) THEN
  PRINT*, "Matrix is non - invertible"
  inverse = 0
  errorflag = -1
  return
  ENDIF
  END DO

  !Make diagonal elements as 1
  DO i = 1 , n
  m = augmatrix(i,i)
  DO j = i , (2 * n)
  augmatrix(i,j) = (augmatrix(i,j) / m)
  END DO
  END DO

  !Reduced right side half of augmented matrix to identity matrix
  DO k = n-1, 1, -1
  DO i =1, k
  m = augmatrix(i,k+1)
  DO j = k, (2*n)
  augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
  END DO
  END DO
  END DO

  !store answer
  DO i =1, n
  DO j = 1, n
  inverse(i,j) = augmatrix(i,j+n)
  END DO
  END DO
  errorflag = 0
  RETURN
END SUBROUTINE FINDInv
!
!
!
END MODULE mod_curved_multiplication




!SUBROUTINE M66INV (A, AINV, OK_FLAG)
!
!
!


! function getMonomials_trip2(x,element) result(M)
!
!     switch element.type
!     case 'tri'
!         switch element.order
!             case 2
!             M = [ ones(size(x,1,KIND=ip),1) x(:,1) x(:,2)...
!                   x(:,1).*x(:,2) x(:,1).*x(:,1) x(:,2).*x(:,2) ]';
!             case 3
!             M = [ ones(size(x,1,KIND=ip),1) x(:,1) x(:,2)...
!                   x(:,1).*x(:,2) x(:,1).*x(:,1) x(:,2).*x(:,2) ...
!                   x(:,1).^3 x(:,1).^2.*x(:,2) x(:,2).^2.*x(:,1) x(:,2).^3 ]';
!         end
!     case 'tet'
!         if(element.order ==1)
!             M = [ ones(size(x,1,KIND=ip),1) x(:,1) x(:,2) x(:,3)  ]';
!         elseif(element.order ==2)
!             M = [ ones(size(x,1,KIND=ip),1) x(:,1) x(:,2) x(:,3) ...
!                   x(:,1).*x(:,2) x(:,1).*x(:,3)  x(:,2).*x(:,3) ...
!                   x(:,1).*x(:,1) x(:,2).*x(:,2)  x(:,3).*x(:,3)  ]';
!         end
!     end
! end
! !




