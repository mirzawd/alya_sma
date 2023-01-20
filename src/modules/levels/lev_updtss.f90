!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_updtss
  !-----------------------------------------------------------------------
  !****f* Levels/lev_updtss
  ! NAME 
  !    lev_updtss
  ! DESCRIPTION
  !    This routine computes the time step size for the level set 
  !    convection equation.
  ! USED BY
  !    lev_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  use mod_communications, only : PAR_MIN
  implicit none 
  integer(ip)      :: ielem,idime,inode,ipoin
  integer(ip)      :: pnode,pelty      
  real(rp)         :: dtcri,dtmin,rdinv,gpdet,h,c,v
  real(rp)         :: elcod(ndime,mnode),elvel(ndime,mnode)
  real(rp)         :: gpcar(ndime,mnode),xjacm(9),xjaci(9)
  !
  ! Compute minimum element time step
  !
  dtmin = dtime

  if(ittim==0) then
     dtinv =1.0_rp/dtime  
  else
     if(kfl_paral/=0) then

        c=0.0_rp
        v=0.0_rp
        rdinv=1.0_rp/real(ndime)
        do ielem = 1,nelem
           pelty=ltype(ielem)
           if( pelty > 0 ) then
              pnode=nnode(pelty)
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do
              rdinv = 1.0_rp/real(ndime)
              ! Look for maximum velocity 
              call lev_velfun(kfl_advec_lev,ndime,pnode,lnods(1,ielem),elcod,elvel)
              do inode=1,pnode
                 do idime=1,ndime
                    v=abs(elvel(idime,inode))
                    c=max(v,c)
                 end do
              end do

              call elmder(pnode,ndime,elmar(pelty)%dercg,elcod,gpcar,gpdet,xjacm,xjaci)
              h     = (elmar(pelty)%weicg*gpdet)**rdinv    ! h
              if(c>1.e-6) then
                 dtcri = h/c                               ! dt=h/c
              else 
                 dtcri = dtime
              endif
              dtmin=min(dtmin,dtcri)
           end if
        end do
     end if
     !
     ! Look for minimum over whole mesh
     !
     if( IPARALL ) call PAR_MIN(dtmin)

     dtcri_lev = dtmin
     dtinv_lev = 1.0_rp/(dtcri_lev*safet_lev)
     if((kfl_timco==1).or.(kfl_timco_lev==1)) dtinv=max(dtinv,dtinv_lev) 

  endif

end subroutine lev_updtss
