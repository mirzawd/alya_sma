!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmmat(&
     pnode,gpsha,gpcar,gpvel,gplev,gpvol,&
     taust,fact1,elmat,elrhs)

  !------------------------------------------------------------------------
  !****f* Levels/lev_elmmat
  ! NAME 
  !    lev_elmmat
  ! DESCRIPTION
  !    ORDER=1:
  ! Level Set convection equation
  !      1. Compute elemental matrix 
  ! 
  ! USES
  ! USED BY
  !    lev_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp  
  use def_domain, only     :  ndime
  use def_levels

  implicit none
  integer(ip), intent(in)  :: pnode  
  real(rp),    intent(in)  :: taust,fact1,gpvol
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode),gpvel(ndime),gplev(ncomp_lev)
  real(rp),    intent(out) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)              :: inode,jnode,itime
  real(rp)                 :: fact2,fact3,tleve

  if(kfl_tisch_lev==1) then

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
        end do

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
        end do

     end if

  else if(kfl_tisch_lev==2) then

     if(ndime==2) then

        do inode=1,pnode

           tleve=gpsha(inode)&                                         ! v
                & +taust*(   gpcar(1,inode)*gpvel(1)&                  ! tau*u.grad(v)
                &          + gpcar(2,inode)*gpvel(2) ) 
           do jnode=1,pnode
              ! Residual
              fact2=fact1*pabdf_lev(1)*gpsha(jnode)&                   ! phi/(theta*dt)
                   & +gpvol*(   gpcar(1,jnode)*gpvel(1)&
                   &          + gpcar(2,jnode)*gpvel(2) )
              ! Assembly
              fact3=tleve*fact2
              elmat(inode,jnode)=elmat(inode,jnode)+fact3

           end do

           do itime=2,kfl_tiacc_lev+1
              elrhs(inode)=elrhs(inode)-tleve*fact1*gplev(itime)*pabdf_lev(itime)
           end do

        end do

     else

        do inode=1,pnode

           tleve=gpsha(inode)&                                         ! v
                & +taust*(   gpcar(1,inode)*gpvel(1)&                  ! tau*u.grad(v)
                &          + gpcar(2,inode)*gpvel(2)& 
                &          + gpcar(3,inode)*gpvel(3) ) 
           do jnode=1,pnode
              ! Residual
              fact2=fact1*pabdf_lev(1)*gpsha(jnode)&                                ! phi/(theta*dt)
                   & +gpvol*(   gpcar(1,jnode)*gpvel(1)&
                   &          + gpcar(2,jnode)*gpvel(2)&
                   &          + gpcar(3,jnode)*gpvel(3) )
              ! Assembly
              fact3=tleve*fact2
              elmat(inode,jnode)=elmat(inode,jnode)+fact3

           end do
           do itime=2,kfl_tiacc_lev+1
              elrhs(inode)=elrhs(inode)-tleve*fact1*gplev(itime)*pabdf_lev(itime)
           end do
        end do

     end if

  end if

end subroutine lev_elmmat
