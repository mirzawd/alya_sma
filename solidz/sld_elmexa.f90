!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_elmexa(&
     itask,pgaus,pnode,elcod,gpsha,gpcar,gpdds,gpvol,&
     eldis,gpdis,gpgdi,elrhs)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmexa
  ! NAME 
  !    sld_elmexa
  ! DESCRIPTION
  !    Compute RHS for exact solution
  ! USES
  !    sld_exacso
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  use def_domain, only : ndime,mnode
  use def_solidz, only : kfl_exacs_sld
  use def_solidz, only : err01_sld,err02_sld,err0i_sld
  use def_solidz, only : err11_sld,err12_sld,err1i_sld
  
  implicit none
  
  integer(ip), intent(in)    :: itask,pgaus,pnode
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: eldis(ndime,pnode)
  real(rp),    intent(out)   :: gpdis(ndime,pgaus)
  real(rp),    intent(out)   :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: elrhs(ndime,pnode)
  
  integer(ip)                :: igaus,inode,idime,jdime
  real(rp)                   :: gpcod(3)
  real(rp)                   :: u(3)            
  real(rp)                   :: dudx(ndime,ndime)   
  real(rp)                   :: d2udx2(27) 
  real(rp)                   :: gprhs(3)
  real(rp)                   :: diffe,absou

  if( kfl_exacs_sld /= 0 ) then

     if( itask == 1 ) then
        !
        ! Assemble RHS
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
           call sld_exacso(2_ip,gpcod,u,dudx,d2udx2,gpdds(1,1,1,1,igaus),gprhs)
           do inode = 1,pnode
              do idime = 1,ndime
                 elrhs(idime,inode) = elrhs(idime,inode) &
                      - gprhs(idime) * gpvol(igaus) * gpsha(inode,igaus)
              end do
           end do
        end do

     else
        !
        ! Compute solution at Gauss points
        !
        ! GPDIS(I)   = ui
        ! GPGDI(I,J) = dui/dxj
        !
        do igaus = 1,pgaus
           do idime = 1,ndime
              gpdis(idime,igaus) = 0.0_rp
              do jdime = 1,ndime
                 gpgdi(idime,jdime,igaus) = 0.0_rp
              end do
           end do
        end do
        do igaus = 1,pgaus
           do idime = 1,ndime
              do inode = 1,pnode
                 gpdis(idime,igaus) = gpdis(idime,igaus) &
                      + eldis(idime,inode) * gpsha(inode,igaus)
              end do
              do jdime = 1,ndime
                 do inode = 1,pnode
                    gpgdi(idime,jdime,igaus) = gpgdi(idime,jdime,igaus) &
                         + eldis(idime,inode) * gpcar(jdime,inode,igaus)
                 end do
              end do
           end do
        end do

        do igaus = 1,pgaus
           do idime = 1,ndime
              gpcod(idime) = 0.0_rp
           end do
           do inode = 1,pnode           
              do idime = 1,ndime
                 gpcod(idime) = gpcod(idime) + elcod(idime,inode) * gpsha(inode,igaus)
              end do
           end do
           call sld_exacso(1_ip,gpcod,u,dudx,d2udx2,gpdds(1,1,1,1,igaus),gprhs)
           !
           ! Errors
           ! 
           do idime = 1,ndime
              diffe        = abs(gpdis(idime,igaus)-u(idime))
              absou        = abs(u(idime))
              err01_sld(1) = err01_sld(1) + diffe * gpvol(igaus)
              err02_sld(1) = err02_sld(1) + diffe * diffe * gpvol(igaus)
              err0i_sld(1) = max(err0i_sld(1),diffe)
              err01_sld(2) = err01_sld(2) + absou * gpvol(igaus)
              err02_sld(2) = err02_sld(2) + absou * absou * gpvol(igaus)
              err0i_sld(2) = max(err0i_sld(2),absou)
              do jdime = 1,ndime
                 diffe        = abs(gpgdi(jdime,idime,igaus)-dudx(idime,jdime))
                 absou        = abs(dudx(idime,jdime))
                 err11_sld(1) = err11_sld(1) + diffe * gpvol(igaus)
                 err12_sld(1) = err12_sld(1) + diffe * diffe * gpvol(igaus)
                 err1i_sld(1) = max(err1i_sld(1),diffe)
                 err11_sld(2) = err11_sld(2) + absou * gpvol(igaus)
                 err12_sld(2) = err12_sld(2) + absou * absou * gpvol(igaus)
                 err1i_sld(2) = max(err1i_sld(2),absou)
              end do
           end do
        end do

     end if

  end if

end subroutine sld_elmexa

