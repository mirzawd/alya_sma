!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ibmdef(itask)
  !-----------------------------------------------------------------------
  !****f* master/ibmdef
  ! NAME
  !    ibmdef
  ! DESCRIPTION
  !    This routines creates the IB structure
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,ii,jj

  real(rp),    pointer    :: Fv(:), Fp(:)
  real(rp),    pointer    :: Tv(:), Tp(:)
  real(rp),    pointer    :: pp_pf(:,:),pp_vf(:,:),pp_pt(:,:),pp_vt(:,:)

  if( itask == 6 ) then

     !-------------------------------------------------------------------
     !
     ! Reduce sum on particles pressure+viscous forces and torques
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        call memgen(zero,nimbo*12,zero)
        ii = 0
        do iimbo = 1,nimbo
           Fv        => imbou(iimbo)%vforce
           Fp        => imbou(iimbo)%pforce
           Tv        => imbou(iimbo)%vtorqu
           Tp        => imbou(iimbo)%ptorqu
           ii        =  ii + 1
           gesca(ii) =  Fv(1)
           ii        =  ii + 1
           gesca(ii) =  Fv(2)
           ii        =  ii + 1
           gesca(ii) =  Fv(3)
           ii        =  ii + 1
           gesca(ii) =  Fp(1)
           ii        =  ii + 1
           gesca(ii) =  Fp(2)
           ii        =  ii + 1
           gesca(ii) =  Fp(3)
           ii        =  ii + 1
           gesca(ii) =  Tv(1)
           ii        =  ii + 1
           gesca(ii) =  Tv(2)
           ii        =  ii + 1
           gesca(ii) =  Tv(3)
           ii        =  ii + 1
           gesca(ii) =  Tp(1)
           ii        =  ii + 1
           gesca(ii) =  Tp(2)
           ii        =  ii + 1
           gesca(ii) =  Tp(3)
        end do

        nparr =  nimbo*12
        parre => gesca
        call par_operat(3_ip)

        ii = 0
        do iimbo = 1,nimbo
           Fv    => imbou(iimbo)%vforce
           Fp    => imbou(iimbo)%pforce
           Tv    => imbou(iimbo)%vtorqu
           Tp    => imbou(iimbo)%ptorqu
           ii    =  ii + 1
           Fv(1) =  gesca(ii)
           ii    =  ii + 1
           Fv(2) =  gesca(ii) 
           ii    =  ii + 1
           Fv(3) =  gesca(ii) 
           ii    =  ii + 1
           Fp(1) =  gesca(ii)
           ii    =  ii + 1
           Fp(2) =  gesca(ii) 
           ii    =  ii + 1
           Fp(3) =  gesca(ii) 
           ii    =  ii + 1
           Tv(1) =  gesca(ii) 
           ii    =  ii + 1
           Tv(2) =  gesca(ii) 
           ii    =  ii + 1
           Tv(3) =  gesca(ii) 
           ii    =  ii + 1
           Tp(1) =  gesca(ii) 
           ii    =  ii + 1
           Tp(2) =  gesca(ii) 
           ii    =  ii + 1
           Tp(3) =  gesca(ii)           
        end do
        call memgen(two,nimbo*12,zero)
     end if

  end if

  if( itask == 7 ) then  ! the same as itask 6 but for rigid bodies

     !-------------------------------------------------------------------
     !
     ! Reduce sum on rigid body pressure+viscous forces and torques
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        call memgen(zero,nrbod*(12+120),zero)
        ii = 0
        do iimbo = 1,nrbod
           Fv        => rbbou(iimbo)%vforce
           Fp        => rbbou(iimbo)%pforce
           Tv        => rbbou(iimbo)%vtorqu
           Tp        => rbbou(iimbo)%ptorqu
           pp_pf => rbbou(iimbo) % pp_pf   ! to postprocess forces separatelly by sets
           pp_vf => rbbou(iimbo) % pp_vf
           pp_pt => rbbou(iimbo) % pp_pt
           pp_vt => rbbou(iimbo) % pp_vt
   
           ii        =  ii + 1
           gesca(ii) =  Fv(1)
           ii        =  ii + 1
           gesca(ii) =  Fv(2)
           ii        =  ii + 1
           gesca(ii) =  Fv(3)
           ii        =  ii + 1
           gesca(ii) =  Fp(1)
           ii        =  ii + 1
           gesca(ii) =  Fp(2)
           ii        =  ii + 1
           gesca(ii) =  Fp(3)
           ii        =  ii + 1
           gesca(ii) =  Tv(1)
           ii        =  ii + 1
           gesca(ii) =  Tv(2)
           ii        =  ii + 1
           gesca(ii) =  Tv(3)
           ii        =  ii + 1
           gesca(ii) =  Tp(1)
           ii        =  ii + 1
           gesca(ii) =  Tp(2)
           ii        =  ii + 1
           gesca(ii) =  Tp(3)

           do jj=1,10
              ii        =  ii + 1
              gesca(ii) =  pp_vf(1,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_vf(2,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_vf(3,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_pf(1,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_pf(2,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_pf(3,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_vt(1,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_vt(2,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_vt(3,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_pt(1,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_pt(2,jj)
              ii        =  ii + 1
              gesca(ii) =  pp_pt(3,jj)
           end do

        end do

        nparr =  nrbod*(12+120)
        parre => gesca
        call par_operat(3_ip)

        ii = 0
        do iimbo = 1,nrbod
           Fv        => rbbou(iimbo)%vforce
           Fp        => rbbou(iimbo)%pforce
           Tv        => rbbou(iimbo)%vtorqu
           Tp        => rbbou(iimbo)%ptorqu
           pp_pf => rbbou(iimbo) % pp_pf   ! to postprocess forces separatelly by sets
           pp_vf => rbbou(iimbo) % pp_vf
           pp_pt => rbbou(iimbo) % pp_pt
           pp_vt => rbbou(iimbo) % pp_vt

           ii    =  ii + 1
           Fv(1) =  gesca(ii)
           ii    =  ii + 1
           Fv(2) =  gesca(ii) 
           ii    =  ii + 1
           Fv(3) =  gesca(ii) 
           ii    =  ii + 1
           Fp(1) =  gesca(ii)
           ii    =  ii + 1
           Fp(2) =  gesca(ii) 
           ii    =  ii + 1
           Fp(3) =  gesca(ii) 
           ii    =  ii + 1
           Tv(1) =  gesca(ii) 
           ii    =  ii + 1
           Tv(2) =  gesca(ii) 
           ii    =  ii + 1
           Tv(3) =  gesca(ii) 
           ii    =  ii + 1
           Tp(1) =  gesca(ii) 
           ii    =  ii + 1
           Tp(2) =  gesca(ii) 
           ii    =  ii + 1
           Tp(3) =  gesca(ii)           

           do jj=1,10
              ii    =  ii + 1
              pp_vf(1,jj) =  gesca(ii)
              ii    =  ii + 1
              pp_vf(2,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_vf(3,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_pf(1,jj) =  gesca(ii)
              ii    =  ii + 1
              pp_pf(2,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_pf(3,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_vt(1,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_vt(2,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_vt(3,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_pt(1,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_pt(2,jj) =  gesca(ii) 
              ii    =  ii + 1
              pp_pt(3,jj) =  gesca(ii)
           end do
        end do
        call memgen(two,nrbod*(12+120),zero)
     end if

  end if

end subroutine ibmdef
