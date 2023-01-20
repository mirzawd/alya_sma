!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_elmsmo(                          &
     pnode,pgaus,pelty,pevat,lelch,kfl_defor,   &
     kfl_elfix,gpcar,gpvol,elcod,elbve,hleng,   &
     voele,asmin,asmax,asele,vomin,vomax,elmat, &
     elrhs)
  !----------------------------------------------------------------------
  !****f* domain/ale_elmsmo
  ! NAME 
  !    ale_elmsmo
  ! DESCRIPTION
  !    Compute the elemental weighted Laplacian matrix
  !    ( alpha * grad X , grad v )
  ! USES
  ! USED BY
  !    smooth
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp 
  use def_elmtyp, only       :  ELEXT
  use def_domain, only       :  mnode,ndime
  use def_elmtyp, only       :  QUA04,HEX08,PEN06
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: pelty
  integer(ip), intent(in)    :: lelch
  integer(ip), intent(in)    :: kfl_defor
  integer(ip), intent(in)    :: kfl_elfix(ndime,*)  
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elbve(ndime,pnode)
  real(rp),    intent(in)    :: hleng(ndime)
  real(rp),    intent(in)    :: voele
  real(rp),    intent(in)    :: asmin
  real(rp),    intent(in)    :: asmax
  real(rp),    intent(inout) :: asele
  real(rp),    intent(in)    :: vomin
  real(rp),    intent(in)    :: vomax
  real(rp),    intent(out)   :: elmat(ndime*pnode,ndime*pnode)
  real(rp),    intent(out)   :: elrhs(ndime*pnode)
  integer(ip)                :: inode,jnode,kdime,igaus,idime
  integer(ip)                :: jdime,pevat,idofn,jdofn
  real(rp)                   :: fact1,gpder(ndime,ndime)
  real(rp)                   :: f1,f2,gpdif(3)

  !----------------------------------------------------------------------
  !
  ! Determine stiffness
  !
  !----------------------------------------------------------------------

  if( kfl_defor == 1 ) then
     !
     ! Only smoothing
     !
     do idime = 1,ndime
        gpdif(idime) = hleng(idime)
     end do 

  else if( kfl_defor == 2 ) then
     !
     ! Conserve small elements only
     !
     do idime = 1,ndime
        gpdif(idime) = (1.0_rp+(1.0_rp-vomin/vomax)/(voele/vomax))
     end do

  else if( kfl_defor == 4 ) then
     !
     ! Sort hleng: hleng(1)=max; hleng(ndime)=min
     ! HLENG(1)     = Max length
     ! HLENG(NDIME) = Min length
     !     
     do idime = 1,ndime
        gpdif(idime) = asele
     end do

  else if( kfl_defor == 5 ) then
     !
     ! Isotropic and uniform Laplacian
     !
     do idime = 1,ndime
        gpdif(idime) = 1.0_rp
     end do

  else if( kfl_defor == 6 ) then
     !
     ! Aspect ratio and element volume
     !     
     if( pelty /= QUA04 .and. pelty /= HEX08 .and. pelty /= PEN06 ) then
        asele = 1.0_rp
     end if
     
     do idime = 1,ndime
        !f1           = (1.0_rp+(1.0_rp-asmin/asmax)/(asele/asmax))**1.0_rp
        !f2           = (1.0_rp+(1.0_rp-vomin/vomax)/(voele/vomax))**1.0_rp
        f1           = asele
        f2           = vomax/voele
        gpdif(idime) = 1.0_rp * f1 * f2
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  pevat = pnode * ndime
  do jdofn = 1,pevat
     do idofn = 1,pevat
        elmat(idofn,jdofn) = 0.0_rp
     end do
     elrhs(jdofn) = 0.0_rp
  end do

  !----------------------------------------------------------------------
  !
  ! Laplacian matrix: ( alpha * grad X , grad v )
  !
  !----------------------------------------------------------------------

  do igaus = 1,pgaus

     gpder = 0.0_rp
     do inode = 1,pnode
        do jdime = 1,ndime
           do idime = 1,ndime
              gpder(idime,jdime) = gpder(idime,jdime) &
                   + elcod(jdime,inode) * gpcar(idime,inode,igaus)
           end do
        end do
     end do

     do inode = 1,pnode
        do jnode = inode+1,pnode
           fact1 = 0.0_rp
           do kdime = 1,ndime
              fact1 = fact1 + gpcar(kdime,inode,igaus) &
                   &        * gpcar(kdime,jnode,igaus)
           end do
           fact1 = fact1 * gpvol(igaus)
           idofn = ( inode - 1 ) * ndime 
           jdofn = ( jnode - 1 ) * ndime 
           do idime = 1,ndime
              idofn = idofn + 1
              jdofn = jdofn + 1
              elmat(idofn,jdofn) = elmat(idofn,jdofn) + fact1 * gpdif(idime) 
              elmat(jdofn,idofn) = elmat(jdofn,idofn) + fact1 * gpdif(idime) 
           end do
        end do
        fact1 = 0.0_rp
        do kdime = 1,ndime
           fact1 = fact1 + gpcar(kdime,inode,igaus) &
                &        * gpcar(kdime,inode,igaus)
        end do
        fact1 = fact1 * gpvol(igaus)
        idofn = ( inode - 1 ) * ndime 
        do idime = 1,ndime
           idofn = idofn + 1
           elmat(idofn,idofn) = elmat(idofn,idofn) + fact1 * gpdif(idime) 
        end do

     end do

  end do


  !----------------------------------------------------------------------
  !
  ! Extension elements
  !
  !----------------------------------------------------------------------

  if( lelch == ELEXT ) then
     do idofn = ndime+1,pevat
        do jdofn = 1,pevat
           elmat(idofn,jdofn) = 0.0_rp
        end do
        elrhs(idofn) = 0.0_rp
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Prescribe Laplacian on walls
  !
  !----------------------------------------------------------------------

  do inode = 1,pnode
     do idime = 1,ndime
        if( kfl_elfix(idime,inode) > 0 ) then
           idofn = ( inode - 1 ) * ndime + idime
           fact1 = elmat(idofn,idofn)
           do jdofn = 1,pevat
              elrhs(jdofn)       = elrhs(jdofn) - elmat(jdofn,idofn) * elbve(idime,inode)
              elmat(jdofn,idofn) = 0.0_rp
              elmat(idofn,jdofn) = 0.0_rp
           end do
           elmat(idofn,idofn) = fact1
           elrhs(idofn)       = fact1 * elbve(idime,inode)
        end if
     end do
  end do

end subroutine ale_elmsmo
