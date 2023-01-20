!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmma2(&
     pnode,gpsha,gpcar,gpvel,gplev,gple0,gpvol,&
     taust,fact1,elmat,elrhs)

  !------------------------------------------------------------------------
  !****f* Levels/lev_elmma2
  ! NAME 
  !    lev_elmmat
  ! DESCRIPTION
  !    ORDER=1:
  ! Level Set convection equation for redistanciation
  !      1. Compute elemental matrix 
  ! 
  ! USES
  ! USED BY
  !    lev_elmope(3_ip)
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp  
  use def_domain, only     :  ndime
!  use def_master, only     :  thicl
  use def_levels

  implicit none
  integer(ip), intent(in)  :: pnode  
  real(rp),    intent(in)  :: taust,fact1,gpvol
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode),gplev(ncomp_lev),gple0
  real(rp),    intent(in)  :: gpvel(ndime)
  real(rp),    intent(out) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)              :: inode,jnode
  real(rp)                 :: fact2,fact3,tleve,gpsig

!     kfl_modif = 0

     if(ndime==2) then

        do inode=1,pnode
           tleve=gpsha(inode)&                                         ! v
                & +taust*(   gpcar(1,inode)*gpvel(1)&                  ! tau*u.grad(v)
                &          + gpcar(2,inode)*gpvel(2) ) 
           do jnode=1,pnode
              ! Residual
              fact2=fact1*gpsha(jnode)&                                ! phi/(theta*dt)
                   & +gpvol*(   gpcar(1,jnode)*gpvel(1)&
                   &          + gpcar(2,jnode)*gpvel(2) )
              ! Assembly
              fact3=tleve*fact2
              elmat(inode,jnode)=elmat(inode,jnode)+fact3

           end do
           elrhs(inode)=elrhs(inode)+tleve*fact1*gplev(2)
           if(gple0==0_rp) then
              gpsig=1_rp
           else
              gpsig=gple0/abs(gple0)
           endif   
           elrhs(inode)=elrhs(inode)+gpvol*tleve*gpsig
        end do

!        if(kfl_modif==1) then
!           !
!           ! Add Diffusion term to stabilice:  ( k*grad(Nj) , grad(Ni) )
!           ! k is determined such that peclet=1.  Therefore k=|u|*h. In order to avoid bringing h to this subroutine 
!           ! we can also obtain k=2*tau*|u|**2  since tau=h/(2*|u|)  for the pure advection problem.
!           ! Moreover we shall only add diffusion far from the interface that is if |LS|>2*thicl
!           if(abs(gple0) > 3.0_rp*thicl) then
!              gpdif = 2.0_rp*taust*((gpvel(1)*gpvel(1))+(gpvel(2)*gpvel(2)))
!           else
!              gpdif = 0.0_rp
!           end if
!           fact0=gpdif*gpvol
!           do inode=1,pnode
!              do jnode=1,pnode
!                 elmat(inode,jnode)=elmat(inode,jnode)+fact1&
!                      *gpcar(1,inode)&
!                      *gpcar(1,jnode)    
!                 elmat(inode,jnode)=elmat(inode,jnode)+fact1&
!                      *gpcar(2,inode)&
!                      *gpcar(2,jnode)
!              end do
!           end do
!        end if
     else

        do inode=1,pnode

           tleve=gpsha(inode)&                                         ! v
                & +taust*(   gpcar(1,inode)*gpvel(1)&                  ! tau*u.grad(v)
                &          + gpcar(2,inode)*gpvel(2)& 
                &          + gpcar(3,inode)*gpvel(3) ) 
           do jnode=1,pnode
              ! Residual
              fact2=fact1*gpsha(jnode)&                                ! phi/(theta*dt)
                   & +gpvol*(   gpcar(1,jnode)*gpvel(1)&
                   &          + gpcar(2,jnode)*gpvel(2)&
                   &          + gpcar(3,jnode)*gpvel(3) )
              ! Assembly
              fact3=tleve*fact2
              elmat(inode,jnode)=elmat(inode,jnode)+fact3

           end do
           elrhs(inode)=elrhs(inode)+tleve*fact1*gplev(2)
           if(gple0==0_rp) then
              gpsig=1_rp
           else
              gpsig=gple0/abs(gple0)
           endif   
           elrhs(inode)=elrhs(inode)+gpvol*tleve*gpsig
        end do

     end if

end subroutine lev_elmma2
