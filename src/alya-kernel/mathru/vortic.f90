!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine vortic(itask)
  !------------------------------------------------------------------------
  !***** mathru/vortic
  ! NAME 
  !    vortic
  ! DESCRIPTION
  !    This routine computes the vorticity and the invariants of the velocity 
  !    gradient.
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master, only    : veloc,vorti,INOTMASTER
  use def_domain, only    : ndime,npoin,memor_dom
  use mod_memory, only    : memory_alloca
  use mod_memory, only    : memory_deallo
  use mod_gradie
  use mod_maths

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,jdime,jtask,kdime
  integer(ip)             :: ndofn
  integer(4)              :: i
  real(rp)                :: omega(ndime,ndime), strain(ndime,ndime), grad_u(ndime,ndime,npoin)
  real(rp)                :: eig(ndime),auxL2(ndime,ndime)
  real(rp)                :: omega2(ndime,ndime), strain2(ndime,ndime)
  real(rp)                :: mod_omega, mod_strain, disto

  jtask = abs(itask)

  if( ndime == 1 ) return

  if( INOTMASTER ) then
     !
     ! Initialization
     !
     if( jtask == 1 ) then
        ndofn = ndime + 1
     else if( jtask == 4 ) then
        ndofn = ndime + 1
     else if( jtask == 2 ) then
        if( ndime == 3 ) then
           ndofn = ndime
        else
           ndofn = 1
        end if
     else if( jtask == 3 ) then
        call memory_deallo(memor_dom,'VORTI','vortic',vorti)
        return
     end if
     if( itask  < 0 .and. .not. associated(vorti) ) then
        call memory_alloca(memor_dom,'VORTI','vortic',vorti,ndofn,npoin)
     end if
     vorti = 0.0_rp

     !
     ! grad(u) using gradie
     !

     call gradie(veloc(:,:,1),grad_u)

     !
     ! Loop over elements
     !
     do ipoin = 1,npoin

        ! Vorticity
        if (ndime == 3) then
           vorti(1,ipoin) = (grad_u(2,3,ipoin) - grad_u(3,2,ipoin))
           vorti(2,ipoin) = (grad_u(3,1,ipoin) - grad_u(1,3,ipoin))
           vorti(3,ipoin) = (grad_u(1,2,ipoin) - grad_u(2,1,ipoin))
        else if (ndime == 2) then
           vorti(1,ipoin) = (grad_u(1,2,ipoin) - grad_u(2,1,ipoin))
        end if

        if( jtask == 1 ) then ! Q-criterion 
           mod_strain = 0.0_rp
           mod_omega  = 0.0_rp
           do idime = 1,ndime
              do jdime = 1,ndime  
                 omega(idime,jdime)  = 0.5_rp * (grad_u(idime,jdime,ipoin)-grad_u(jdime,idime,ipoin))
                 strain(idime,jdime) = 0.5_rp * (grad_u(idime,jdime,ipoin)+grad_u(jdime,idime,ipoin))
                 mod_omega  = mod_omega  + omega(idime,jdime)**2_ip
                 mod_strain = mod_strain + strain(idime,jdime)**2_ip
              end do
           end do
           vorti(ndime+1,ipoin) = 0.5_rp *(mod_omega - mod_strain)
        else if (jtask == 4 ) then ! Lambda-2 method
           do idime = 1,ndime
              do jdime = 1,ndime  
                 omega(idime,jdime)  = 0.5_rp * (grad_u(idime,jdime,ipoin)-grad_u(jdime,idime,ipoin))
                 strain(idime,jdime) = 0.5_rp * (grad_u(idime,jdime,ipoin)+grad_u(jdime,idime,ipoin))
              end do
           end do
           omega2  = 0.0_rp
           strain2 = 0.0_rp
           auxL2 = 0.0_rp
           do idime = 1,3
              do jdime = 1,3
                 do kdime = 1,ndime
                    omega2(idime,jdime)  = omega2(idime,jdime)  + omega(idime,kdime)*omega(kdime,jdime)  
                    strain2(idime,jdime) = strain2(idime,jdime) + strain(idime,kdime)*strain(kdime,jdime)  
                 end do
                 auxL2(idime,jdime) = omega2(idime,jdime) + strain2(idime,jdime)
              end do
           end do
           call maths_eigen_3x3_symmetric_matrix(auxL2,eig)
           vorti(ndime+1,ipoin) = eig(2)
        else if (jtask == 5 ) then ! Distortion time
           !  http://dx.doi.org/10.1016/j.jbiomech.2015.11.049
           strain2 = 0.0_rp
           disto=0.0_rp
           do idime = 1,ndime
              do jdime = 1,ndime  
                 strain(idime,jdime) = 0.5_rp * (grad_u(idime,jdime,ipoin)+grad_u(jdime,idime,ipoin))
              end do
           end do

           strain2=matmul(strain,strain)
           
           do i=1,ndime
              disto=disto+strain2(i,i)
           enddo

           if(disto>0.0_rp) disto=1/sqrt(disto/2.0_rp)

           vorti(1,ipoin) = disto


            
        endif
     end do
  end if

end subroutine vortic
