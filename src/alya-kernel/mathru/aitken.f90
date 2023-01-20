!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    aitken.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Compute the Aitken parameter for a certain mask of the domain.
!> @details Compute the Aitken parameter for a certain mask of ghe domain.
!>          The ITASK input parameter sets the mask:
!>          itask= 0: go through all the domain's nodes (default)
!>          itask= 1: go through all the contact nodes on the fuid side
!>          itask= 2: ...
!> @} 
!-----------------------------------------------------------------------
subroutine aitken(itask,nn1,nn2,nn3,v1,v2,v3,kcom1,kcom2,kcom3,kdime,relpa)
  use def_kintyp
  use def_domain
  use def_master
  use def_elmtyp, only : NODE_CONTACT_SOLID, NODE_CONTACT_FLUID
  implicit none
  integer(ip), intent(in)    :: nn1 !> nodal dimension of the first vector
  integer(ip), intent(in)    :: nn2 !> nodal dimension of the second vector
  integer(ip), intent(in)    :: nn3 !> nodal dimension of the third vector
  integer(ip), intent(in)    :: kcom1 !> mesh dimension of the first vector
  integer(ip), intent(in)    :: kcom2 !> mesh dimension of the second
  integer(ip), intent(in)    :: kcom3 !> mesh dimension of the third vector
  integer(ip), intent(in)    :: kdime !> local dimension for the Aitken parameter
  integer(ip), intent(in)    :: itask !> mask 
  real(rp),    intent(in)    :: v1(*),v2(*)
  real(rp),    intent(inout) :: v3(*)
  real(rp),    intent(inout) :: relpa
  real(rp),    target        :: xfact(2)
  integer(ip)                :: ii,i1,i2,i3,idime
  real(rp)                   :: deltu,rdiff

  xfact(1) = 0.0_rp
  xfact(2) = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential case
     !
     do ii=1,npoin
        if (itask == 1) then
           if (lnoch(ii) .ne. NODE_CONTACT_FLUID) cycle
        end if
        i1=(ii-1)*nn1+kcom1
        i2=(ii-1)*nn2+kcom2
        i3=(ii-1)*nn3+kcom3
        do idime=0,kdime-1
           deltu        = v1(i1+idime)-v2(i2+idime)
           rdiff        = v3(i3+idime)-deltu
           xfact(1)     = xfact(1) + rdiff*deltu
           xfact(2)     = xfact(2) + rdiff*rdiff
           v3(i3+idime) = deltu
        end do
     end do

     if (itask == 1) then
        ! wall suggest to do this correction, que te da directamente el w
        if(xfact(2)/=0.0_rp) then
           relpa = - relpa * xfact(1)/xfact(2)
        else
           relpa = 1.0_rp
        end if
     else
        if(xfact(2)/=0.0_rp) then
           relpa = relpa + (relpa-1.0_rp)*xfact(1)/xfact(2)
        else
           relpa = 0.0_rp
        end if
     end if

  else if(kfl_paral>=0) then
     ! 
     ! Parallel case
     !
     if(kfl_paral>=1) then
        do ii=1,npoi1
           if (itask == 1) then
              if (lnoch(ii) .ne. NODE_CONTACT_FLUID) cycle
           end if
           i1=(ii-1)*nn1+kcom1
           i2=(ii-1)*nn2+kcom2
           i3=(ii-1)*nn3+kcom3
           do idime=0,kdime-1  
              deltu        = v1(i1+idime)-v2(i2+idime)
              rdiff        = v3(i3+idime)-deltu
              xfact(1)     = xfact(1) + rdiff*deltu
              xfact(2)     = xfact(2) + rdiff*rdiff
              v3(i3+idime) = deltu
           end do
        end do
        do ii=npoi2,npoi3
           if (itask == 1) then
              if (lnoch(ii) .ne. NODE_CONTACT_FLUID) cycle
           end if
           i1=(ii-1)*nn1+kcom1
           i2=(ii-1)*nn2+kcom2
           i3=(ii-1)*nn3+kcom3
           do idime=0,kdime-1  
              deltu        = v1(i1+idime)-v2(i2+idime)
              rdiff        = v3(i3+idime)-deltu
              xfact(1)     = xfact(1) + rdiff*deltu
              xfact(2)     = xfact(2) + rdiff*rdiff
              v3(i3+idime) = deltu
           end do
        end do
     end if
     nparr =  2
     parre => xfact
     call par_operat(3_ip)
     if(xfact(2)/=0.0_rp) then
        relpa = relpa + (relpa-1.0_rp)*xfact(1)/xfact(2)
     else
        relpa = 0.0_rp
     end if

  end if

end subroutine aitken
