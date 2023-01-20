!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmdi3(&
     pnode,pevat,lnods,elauu,elaup,elapp,elapu,elrbu,elrbp,elcmm)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmdi3
  ! NAME 
  !    nsi_elmdi3
  ! DESCRIPTION
  ! USES
  !    nsi_rotma3
  ! USED BY
  !    nsi_elmop3
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,exnor,skcos,lpoty
  use def_nastin, only       :  kfl_confi_nsi,nodpr_nsi,kfl_local_nsi,&
       &                        kfl_fixno_nsi,kfl_fixrs_nsi,bvess_nsi,&
       &                        skcos_nsi,valpr_nsi,kfl_imppr_nsi,&
       &                        kfl_fixpr_nsi,kfl_stabi_nsi

  use def_kermod, only       : kfl_adj_prob
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elauu(pevat,pevat)
  real(rp),    intent(inout) :: elaup(pevat,pnode)
  real(rp),    intent(inout) :: elapp(pnode,pnode)
  real(rp),    intent(inout) :: elapu(pnode,pevat)
  real(rp),    intent(inout) :: elrbu(pevat)
  real(rp),    intent(inout) :: elrbp(pnode)
  real(rp),    intent(inout) :: elcmm(pevat,pevat)
  real(rp)                   :: adiag,worma(9),vimpr_nsi,adia2
  integer(ip)                :: inode,ipoin,ievav,idime
  integer(ip)                :: jevav,iroty,ibopo,jnode

  if( kfl_stabi_nsi /= -1 ) then
     !
     ! Prescribe one pressure degree of freedom if the flow is confined
     !
     if( kfl_confi_nsi >= 0 .and. nodpr_nsi > 0 ) then

        do inode = 1,pnode
           ipoin = lnods(inode)
           if( ipoin == nodpr_nsi ) then
              adiag = elapp(inode,inode)
              do jevav = 1,pnode * ndime
                 elrbu(jevav)       = elrbu(jevav)-elaup(jevav,inode) * valpr_nsi
                 elaup(jevav,inode) = 0.0_rp 
                 elapu(inode,jevav) = 0.0_rp
              end do
              do jnode = 1,pnode
                 elrbp(jnode)       = elrbp(jnode)-elapp(jnode,inode) * valpr_nsi
                 elapp(jnode,inode) = 0.0_rp 
                 elapp(inode,jnode) = 0.0_rp                 
              end do
              elapp(inode,inode) = adiag
              elrbp(inode)       = valpr_nsi*adiag
           end if
        end do

     end if
     !
     ! Impose pressure at nodes with kfl_fixpr_nsi - We are modifying the monolithic problem
     !
     if ( kfl_imppr_nsi > 0 ) then
        vimpr_nsi = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode)
           if ( kfl_fixpr_nsi(1,ipoin) > 0 ) then 
              adiag = elapp(inode,inode)
              do jevav = 1,pnode * ndime
                 elrbu(jevav)       = elrbu(jevav)-elaup(jevav,inode) * vimpr_nsi
                 elaup(jevav,inode) = 0.0_rp 
                 elapu(inode,jevav) = 0.0_rp
              end do
              do jnode = 1,pnode
                 elrbp(jnode)       = elrbp(jnode)-elapp(jnode,inode) * vimpr_nsi
                 elapp(jnode,inode) = 0.0_rp 
                 elapp(inode,jnode) = 0.0_rp                 
              end do
              elapp(inode,inode) = adiag
              elrbp(inode)       = vimpr_nsi*adiag
           end if
        end do

     end if

  else

     if( kfl_confi_nsi >= 0 .and. nodpr_nsi > 0 ) then

        do inode = 1,pnode
           ipoin = lnods(inode)
           if( ipoin == nodpr_nsi ) then
              adiag = 1.0_rp
              do jevav = 1,pnode * ndime
                 elrbu(jevav)       = elrbu(jevav)-elaup(jevav,inode) * valpr_nsi
                 elaup(jevav,inode) = 0.0_rp 
                 elapu(inode,jevav) = 0.0_rp
              end do
              do jnode = 1,pnode
                 elrbp(jnode)       = elrbp(jnode)-elapp(jnode,inode) * valpr_nsi
                 elapp(jnode,inode) = 0.0_rp 
                 elapp(inode,jnode) = 0.0_rp                 
              end do
              elapp(inode,inode) = adiag
              elrbp(inode)       = valpr_nsi*adiag
           end if
        end do

     end if

  end if
  !
  ! Rotate the nodal matrices, the element nodal velocities and RHS to 
  ! prescribe boundary conditions in a skew-system, either the tangent
  ! one or another prescribed by the user. Also, periodical boundary 
  ! conditions are accounted for.
  !
  if( kfl_local_nsi == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        ibopo = lpoty(ipoin)
        if( ibopo > 0 ) then
           iroty = kfl_fixrs_nsi(ipoin)
           if( iroty == -1 ) then                                    ! Tangent system
              call nsi_rotma3(&
                   &      inode,pnode,pevat,elauu,elaup,&
                   &      elapu,elrbu,elcmm,exnor(1,1,ibopo),worma)
           else if( iroty >= 1 ) then                                ! Given system
              call nsi_rotma3(&
                   &      inode,pnode,pevat,elauu,elaup,&
                   &      elapu,elrbu,elcmm,skcos(1,1,iroty),worma)
           else if( iroty == -2 ) then                               ! Given system
              call nsi_rotma3(&
                   &      inode,pnode,pevat,elauu,elaup,&
                   &      elapu,elrbu,elcmm,skcos_nsi(1,1,ibopo),worma)
           else if( iroty == -3 ) then                               ! Geometrical normal
              call nsi_rotma3(&
                   &      inode,pnode,pevat,elauu,elaup,&
                   &      elapu,elrbu,elcmm,skcos(1,1,ibopo),worma)
           end if
        end if
     end do
  end if
  !  
  ! Velocity Dirichlet boundary conditions
  !
  do inode = 1,pnode
     ievav = (inode-1) * ndime
     ipoin = lnods(inode)
     do idime = 1,ndime
        ievav = ievav+1
        if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
           adiag = elauu(ievav,ievav)
           adia2 = elcmm(ievav,ievav)
           do jevav = 1,pevat
              elauu(ievav,jevav) = 0.0_rp
              elcmm(ievav,jevav) = 0.0_rp
           end do
           do jevav = 1,pevat
              elrbu(jevav)       = elrbu(jevav) - elauu(jevav,ievav) * bvess_nsi(idime,ipoin,1) 
              elauu(jevav,ievav) = 0.0_rp
              ! Since nothing is sent to the rhs - This implies that the Consistent mass matrix is only ready for cases were
              ! the value to be prescribed is 0. This is what happens for the end of step velocity since one solves for an increment
              elcmm(ievav,jevav) = 0.0_rp
           end do
           do jnode = 1,pnode
              elaup(ievav,jnode) = 0.0_rp
              if (kfl_adj_prob == 1) elapu(jnode,ievav) = 0.0_rp     ! added for the adjoint problem
           end do
           elauu(ievav,ievav) = adiag
           elrbu(ievav)       = adiag * bvess_nsi(idime,ipoin,1)
           elcmm(ievav,ievav) = adiag
        end if
     end do
  end do

end subroutine nsi_elmdi3
