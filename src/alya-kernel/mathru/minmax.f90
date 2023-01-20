!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine minmax(ndim1,ndim2,ndim3,vecto,vemin,vemax)
  !------------------------------------------------------------------------
  !****f* mathru/minmax
  ! NAME 
  !    minmax
  ! DESCRIPTION
  !    Compute the minimum and maximum of a vector 
  ! USES
  ! USED BY
  !    Modules: *_cvgunk
  !***
  !------------------------------------------------------------------------
  use      def_kintyp
  use      def_master
  use      mod_communications, only : PAR_MIN
  use      mod_communications, only : PAR_MAX
  implicit none
  integer(ip), intent(in)  :: ndim1,ndim2,ndim3
  real(rp),    intent(in)  :: vecto(ndim1,ndim2)
  real(rp),    intent(out) :: vemin,vemax
  integer(ip)              :: idim1,idim2,ndim4
  real(rp)                 :: uvalu

  vemin= huge(1.0_rp)
  vemax=-huge(1.0_rp)

  if(kfl_paral==-1.or.(kfl_paral==0.and.kfl_ptask==0)) then
     ! 
     ! Sequential case
     !
     if(ndim1==1) then
        do idim2=1,ndim2
           if(vecto(1,idim2)>vemax) vemax=vecto(1,idim2)
           if(vecto(1,idim2)<vemin) vemin=vecto(1,idim2)
        end do

     else if(ndim3<0) then
        ndim4=-ndim3
        do idim2=1,ndim2
           if(vecto(ndim4,idim2)>vemax) vemax=vecto(ndim4,idim2)
           if(vecto(ndim4,idim2)<vemin) vemin=vecto(ndim4,idim2)
        end do

     else
        ndim4=min(ndim1,ndim3)
        do idim2=1,ndim2
           uvalu=0.0_rp
           do idim1=1,ndim4
              uvalu=uvalu+vecto(idim1,idim2)*vecto(idim1,idim2)
           end do
           uvalu=sqrt(uvalu)
           if(uvalu>vemax) vemax=uvalu
           if(uvalu<vemin) vemin=uvalu
        end do
     end if

  else
     !
     ! Parallel case
     !
     if(kfl_paral>=1) then
        if(ndim1==1) then
           do idim2=1,ndim2
              if(vecto(1,idim2)>vemax) vemax=vecto(1,idim2)
              if(vecto(1,idim2)<vemin) vemin=vecto(1,idim2)
           end do
        else
           ndim4=min(ndim1,ndim3)
           do idim2=1,ndim2
              uvalu=0.0_rp
              do idim1=1,ndim4
                 uvalu=uvalu+vecto(idim1,idim2)*vecto(idim1,idim2)
              end do
              uvalu=sqrt(uvalu)
              if(uvalu>vemax) vemax=uvalu
              if(uvalu<vemin) vemin=uvalu
           end do
        end if
     end if
     !
     ! Minimum
     !
     call PAR_MIN(vemin)
     !
     ! Maximum
     ! 
     call PAR_MAX(vemax)
    
  end if

end subroutine minmax
