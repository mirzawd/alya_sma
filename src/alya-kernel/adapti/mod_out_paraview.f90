!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



MODULE mod_out_paraview
  use def_kintyp_basic,       only: ip,rp,lg
  use def_kintyp_mesh_basic,  only: mesh_type_basic
  
  implicit none
  
  private
  
  public:: out_paraview_inp
  
CONTAINS
  
subroutine out_paraview_inp(mesh,filename,elemField,nodeField,nodeFields)
!*************************************************
!*
!*     Output for Paraview of the surf mesh
!*
!*************************************************
  implicit none
  !
  !*** input variables
  !
  type(mesh_type_basic) , intent(in)   :: mesh
  character(len=*),optional,intent(in) :: filename
  real(rp),  optional   , intent(in)   :: elemField(:)
  real(rp),  optional   , intent(in)   :: nodeField(:)
  real(rp),  optional   , intent(in)   :: nodeFields(:,:)

  character(len=200) :: fname
  integer(ip)        , pointer    :: lnods(:,:)
  real(rp)           , pointer    :: x(:,:)
  integer(ip) :: npoin,nnode,nelem
  !
  !*** local variables
  !
  integer(ip) :: ielem,inode
  character(len=4) :: elemType
  
  integer(4)                                      :: iunit4
  integer(ip)                                     :: ioerr
  logical(lg)                                     :: opened
  !
  if(present(filename)) then
    fname = trim(filename)
  else
    fname = 'MESH'
  end if

  nullify(lnods)
  lnods=>mesh%lnods
  nullify(x)
  x    =>mesh%coord
  nnode=size(mesh%lnods,1)
  nelem=size(mesh%lnods,2)
  npoin=size(x,2)
  !
  !*** file to write
  !
  do iunit4 = 90_4,1000_4
     inquire(unit=iunit4,opened=opened,iostat=ioerr)
     if( ioerr /= 0 )  cycle
     if( .not. opened ) exit
  end do
  open(iunit4,file=trim(fname)//'.inp',status='unknown')
  !
  !*** header
  !
  if((nnode.eq.4).or.(nnode.eq.3)) then
    if(present(elemField).and.(.not.present(nodeField))) then
      write(iunit4,'(2(i8,1x),3(1x,i8))')  npoin,nelem,0,1,0
    else if(present(nodeField).and.(.not.present(elemField))) then
      write(iunit4,'(2(i8,1x),3(1x,i8))')  npoin,nelem,1,0,0
    else if(present(nodeField).and.present(elemField)) then
      write(iunit4,'(2(i8,1x),3(1x,i8))')  npoin,nelem,1,1,0
    else if(present(nodeFields)) then
      write(iunit4,'(2(i8,1x),3(1x,i8))')  npoin,nelem,size(nodeFields,2),0,0
    else
      write(iunit4,'(2(i8,1x),3(1x,i8))')  npoin,nelem,0,0,0
    end if
  else
     call runend(' out_inp: incorrect number of nodes')
  end if
  !
  !*** coordinates
  !
  do inode = 1,npoin
    if(size(x,1)==2) then
      write(iunit4,'(i8,3(1x,E15.7))') inode,x(1,inode),x(2,inode),0.0
    else
      write(iunit4,'(i8,3(1x,E15.7))') inode,x(1,inode),x(2,inode),x(3,inode)
    end if
  end do
  !
  !*** elements
  !
  do ielem = 1,nelem
    if(nnode==3) then
      elemType = 'tri '
      write(iunit4,'(i8,1x,i1,1x,A,3(1x,i8))')  ielem,0,elemType,(lnods(inode,ielem),inode=1,nnode)
    else if(nnode==4) then
      elemType = 'tet'
      write(iunit4,'(i8,1x,i1,1x,A,4(1x,i8))')  ielem,0,elemType,(lnods(inode,ielem),inode=1,nnode)
    end if
  end do
  !
  !*** write nodal/elemental fields
  !
  if(present(nodeFields)) then
    if(size(nodeFields,2)==2) then
      write(iunit4,'(A)') '2 1 1'
      write(iunit4,'(A)') 'nodeField1, None'
      write(iunit4,'(A)') 'nodeField2, None'

      do inode = 1,npoin
        write(iunit4,'(i8,1x,F9.3,1x,E15.7)') inode,nodeFields(inode,:)
      end do
    else
      call runend('progtram paraview output of more nodal fields')
    end if
  else if(present(nodeField)) then
    write(iunit4,'(2(1x,i1))') 1,1
    write(iunit4,'(A)') 'nodeField, None'

    do inode = 1,npoin
      write(iunit4,'(i8,1x,E15.7)') inode,nodeField(inode)
    end do
  end if
  if(present(elemField)) then
    write(iunit4,'(2(1x,i1))') 1,1
    write(iunit4,'(A)') 'elemField, None'

    do ielem = 1,nelem
      write(iunit4,'(i8,1x,E15.7)') ielem,elemField(ielem)
    end do
  end if

  close(iunit4)
end subroutine out_paraview_inp



END MODULE mod_out_paraview
