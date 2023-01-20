!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_elmpno(&
     pnode,pelty,lelch,kfl_elfix,kfl_iffix,&
     elcod,eless,asele,elmat,elrhs)
  
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_elmpno
  ! NAME
  !    elmpno
  ! DESCRIPTION
  !    Assemble the Laplacian-like matrix
  ! USED BY
  !    Turnon 
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp, only     :  ip,rp
  use def_elmtyp, only     :  QUA04,HEX08,PEN06,ELEXT
  use def_domain, only     :  ndime
!  use def_master, only     :  kfl_fixno_ale
  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pelty
  integer(ip), intent(in)  :: lelch
  integer(ip), intent(in)  :: kfl_elfix(ndime,*)  
  integer(ip), intent(in)  :: kfl_iffix
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: eless(ndime,pnode)
  real(rp),    intent(in)  :: asele
  real(rp),    intent(out) :: elmat(ndime*pnode,ndime*pnode)
  real(rp),    intent(out) :: elrhs(ndime,pnode)
  integer(ip)              :: inode,jnode,idime
!  integer(ip)              :: ipoin
  integer(ip)              :: idofn,jdofn,jdime
  real(rp)                 :: fact1

  elmat = 0.0_rp
  elrhs = 0.0_rp
  
  do inode = 1,pnode
     do idime = 1,ndime
        idofn = (inode-1)*ndime+idime
        do jnode = 1,pnode
           jdofn = (jnode-1)*ndime+idime
           elmat(idofn,jdofn) = -1.0_rp / real(pnode-1,rp)
        end do
        elmat(idofn,idofn) = 1.0_rp
     end do
  end do

  if( asele > 5.0_rp ) then

     if(       ( ndime == 2 .and. pelty == QUA04 ) &
          .or. ( ndime == 3 .and. ( pelty == HEX08 .or. pelty == PEN06 ) ) ) then
        !
        ! Fix node
        !
        !do inode = 1,pnode
        !   do jnode = 1,pnode
        !      elmat(inode,jnode) = 0.0_rp
        !   end do
        !   elmat(inode,inode) = 1.0e5_rp
        !   elrhs(inode)       = 1.0e5_rp
        !end do
        !do inode = 1,pnode
        !   ipoin = lnods(inode)
        !   kfl_fixno_ale(1,ipoin) = 1
        !end do
        !
        ! Relax stiffness
        !
        do inode = 1,pnode
           do idime = 1,ndime
              idofn = (inode-1)*ndime+idime             
              elmat(idofn,idofn) = elmat(idofn,idofn) + asele
              elrhs(idime,inode) = elrhs(idime,inode) + asele
           end do
        end do
     end if

  end if
  !
  ! Assemble RHS
  !
  do inode = 1,pnode
     do idime = 1,ndime
        elrhs(idime,inode) = elrhs(idime,inode) * elcod(idime,inode)
     end do
  end do  

  !----------------------------------------------------------------------
  !
  ! Extension elements
  !
  !----------------------------------------------------------------------
  
  if( lelch == ELEXT ) then
     do idofn = ndime+1,ndime*pnode
        do jdofn = 1,ndime*pnode
           elmat(idofn,jdofn) = 0.0_rp
        end do
     end do
     do inode = 1,pnode
        elrhs(1:ndime,2:pnode) = 0.0_rp
     end do
  end if
  
  !----------------------------------------------------------------------
  !
  ! Prescribe Laplacian on walls
  !
  !----------------------------------------------------------------------
  
  if( kfl_iffix == 0 ) then
     
     do inode = 1,pnode
        do idime = 1,ndime
           if( kfl_elfix(idime,inode) > 0 ) then
              idofn = ( inode - 1 ) * ndime + idime
              fact1 = elmat(idofn,idofn)
              do jnode = 1,pnode
                 do jdime = 1,ndime
                    jdofn = (jnode-1)*ndime+jdime
                    elrhs(jdime,jnode) = elrhs(jdime,jnode) - elmat(jdofn,idofn) * eless(idime,inode)
                 end do
              end do
              do jdofn = 1,ndime*pnode
                 elmat(jdofn,idofn) = 0.0_rp
                 elmat(idofn,jdofn) = 0.0_rp
              end do
              if( fact1 == 0.0_rp ) then
                 elmat(idofn,idofn) = 1.0_rp
              else
                 elmat(idofn,idofn) = fact1 
              end if
              elrhs(idime,inode)    = elmat(idofn,idofn) * eless(idime,inode)
           end if
        end do
     end do    
  end if

end subroutine ale_elmpno
