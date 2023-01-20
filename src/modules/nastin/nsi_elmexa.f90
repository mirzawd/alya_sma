!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmexa(&
     pgaus,pnode,gpsha,elcod,gpden,gpvis,gppor,gpgvi,&
     cutim,baloc,gprhs,gprhc,gprh2)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmexa
  ! NAME 
  !    nsi_elmexa
  ! DESCRIPTION
  !    Compute RHS for exact solution
  ! USES
  !    nsi_exacso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_nastin, only     :  kfl_exacs_nsi
  implicit none
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: gpsha(pnode,*)
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: gpden(*)
  real(rp),    intent(in)  :: gpvis(*)
  real(rp),    intent(in)  :: gppor(*)
  real(rp),    intent(in)  :: gpgvi(ndime,*)
  real(rp),    intent(in)  :: cutim
  real(rp),    intent(in)  :: baloc(ndime)
  real(rp),    intent(out) :: gprhs(ndime,*)
  real(rp),    intent(out) :: gprhc(*)
  real(rp),    intent(out) :: gprh2(*)
  integer(ip)              :: igaus,inode,idime
  real(rp)                 :: dummr,gpcod(3)

  if( kfl_exacs_nsi /= 0 ) then

     if( pgaus == -1 ) then
        !
        ! Contribution on a a particular boundary Gauss point
        !
        igaus = 1
        do idime = 1,ndime
           gpcod(idime) = 0.0_rp
        end do
        do inode = 1,pnode           
           do idime = 1,ndime
              gpcod(idime) = gpcod(idime) + elcod(idime,inode) * gpsha(inode,1)
           end do
        end do
        call nsi_exacso(&
             3_ip,gpcod,gpden(igaus),gpvis(igaus),&
             gppor(igaus),gpgvi(1,igaus),&
             dummr,dummr,dummr,dummr,baloc,&
             gprhs(1,igaus),gprhc(igaus),gprh2(igaus))

     else
        !
        ! Contribution on an element
        !
        do igaus = 1,pgaus
           do idime = 1,ndime
              gpcod(idime) = 0.0_rp
           end do
           do inode = 1,pnode           
              do idime = 1,ndime
                 gpcod(idime) = gpcod(idime) + elcod(idime,inode) * gpsha(inode,igaus)
              end do
           end do
           call nsi_exacso(&
                2_ip,gpcod,gpden(igaus),gpvis(igaus),&
                gppor(igaus),gpgvi(1,igaus),&
                dummr,dummr,dummr,dummr,baloc,&
                gprhs(1,igaus),gprhc(igaus),gprh2(igaus))
        end do

     end if
     
  end if

end subroutine nsi_elmexa

