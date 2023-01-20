!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_assresdiff
    
  implicit none
  
  public :: assresdiff

contains


  subroutine assresdiff(ndofn_glo,ndofn_ele,idesvar,pnode,lnods,elresdiff,resdiff,kdofn)
    !------------------------------------------------------------------------
    !****f* mathru/assresdiff
    ! NAME
    !    assresdiff
    ! DESCRIPTION
    !    Assembly of the RHS
    ! USES
    ! USED BY
    !    *_elmope
    !***
    !------------------------------------------------------------------------
    use def_kintyp, only        :  ip,rp
    use def_domain, only        :  npoin
    use def_kermod, only        :  kfl_ndvars_opt
    implicit none


    integer(ip),  intent(in), optional :: kdofn
    integer(ip),  intent(in)    :: ndofn_glo,ndofn_ele,pnode,idesvar
    integer(ip),  intent(in)    :: lnods(pnode)
    real(rp),     intent(in)    :: elresdiff(*)
    real(rp),     intent(inout) :: resdiff(kfl_ndvars_opt,*)
    integer(ip)                 :: inode,ipoin,idofl,idofg,idofn,ndof2

    if( present(kdofn) ) then

       if( ndofn_ele == 1 ) then
          do inode = 1,pnode
             ipoin = lnods(inode)
             idofg = (ipoin-1) * ndofn_glo + kdofn
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             resdiff(idesvar,idofg) = resdiff(idesvar,idofg) + elresdiff(inode)
          end do
       else
          call runend('MATRIX_ASSRHS: VERY STRANGE RHS ASSEMBLY INDEED...')
       end if

    else
       if(ndofn_glo==1) then
          !
          ! 1 DOF
          !
          do inode=1,pnode
             ipoin=lnods(inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             resdiff(idesvar,ipoin)=resdiff(idesvar,ipoin)+elresdiff(inode)
          end do
       else if(ndofn_glo==2) then
          !
          ! 2 DOF's
          !
          do inode=1,pnode
             ipoin=lnods(inode)
             idofg=2*ipoin-1
             idofl=2*inode-1
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             resdiff(idesvar,idofg)=resdiff(idesvar,idofg)+elresdiff(idofl)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             resdiff(idesvar,idofg+1)=resdiff(idesvar,idofg+1)+elresdiff(idofl+1)
          end do
       else if(ndofn_glo>2) then
          !
          ! >2 DOF's
          !
          do inode=1,pnode
             ipoin=lnods(inode)
             idofg=(ipoin-1)*ndofn_glo
             idofl=(inode-1)*ndofn_glo
             do idofn=1,ndofn_glo
                idofg=idofg+1
                idofl=idofl+1
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                resdiff(idesvar,idofg)=resdiff(idesvar,idofg)+elresdiff(idofl)
             end do
          end do
          
       else if(ndofn_glo<0) then
          !
          ! >2 DOF's
          !
          ndof2=abs(ndofn_glo)
          do inode=1,pnode
             ipoin=lnods(inode)
             idofg=(ipoin-1)*ndof2
             idofl=(inode-1)*ndof2
             do idofn=1,ndof2
                idofg=(idofn-1)*npoin+ipoin
                idofl=idofl+1
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                resdiff(idesvar,idofg)=resdiff(idesvar,idofg)+elresdiff(idofl)
             end do
          end do
       end if

    end if

  end subroutine assresdiff

end module mod_assresdiff
