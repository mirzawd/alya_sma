!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outtra(itask)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_outtra
  ! NAME 
  !    nsi_outtra
  ! DESCRIPTION
  !    This routine computes the traction on boundaries 
  ! USES
  ! USED BY
  !    nsi_output
  !  itask = 0  only on boundary nodes
  !  itask = 1  on all nodes, =0 in interior nodes
  !  outputs gevec (does not work in complex geometries)
  !   
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_memchk
  use def_kermod
  use mod_gradie
  use mod_memory
  use mod_ker_proper
  use mod_frivel,     only : frivel
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ibopo,idime,dummi
  integer(8)              :: metmp(2)=0_8
  integer(4)              :: istat
  real(rp)                :: n1,n2,n3,s11,s12,s22,s13,s23,s33,us2
  real(rp)                :: sn1,sn2,sn3,rho,nu,u,ustar
  real(rp),    pointer    :: visco_tmp(:) 

  if( delta_nsi > 0.0_rp .or. kfl_delta == 1 ) then 
     !
     ! Wall law: no velocity gradient needed
     !
  else
     ! 
     ! Velocity gradients GRADV_NSI
     ! 
     call memory_alloca(mem_modul(1:2,modul),'GRADV_NSI','nsi_outtra',gradv_nsi,ntens,npoin)
     call gradie(veloc(1:ndime,1:npoin,1),gradv_nsi)
     nullify(  visco_tmp )
     allocate( visco_tmp(npoin) )
     call ker_proper('VISCO','NPOIN',dummi,dummi,visco_tmp)
     do ipoin = 1,npoin  
        gradv_nsi(1:ntens,ipoin) = gradv_nsi(1:ntens,ipoin) * visco_tmp(ipoin) 
     end do
     deallocate( visco_tmp )
  end if 

  if( itask == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Tangential traction on boundary nodes only
     !
     !-------------------------------------------------------------------

     if( delta_nsi > 0.0_rp .or. kfl_delta == 1 ) then 
        !
        ! Wall law
        !
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if(ibopo>=1) then
              rho = prope_nsi(1,ipoin)
              nu  = prope_nsi(2,ipoin)/rho
              u   = 0.0_rp
              do idime=1,ndime
                 u = u + veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
              end do
              u = sqrt(u)
              if( kfl_rough >  0 ) rough_dom = rough(ipoin) 
              if( kfl_delta == 1 ) then
                 call frivel(kfl_ustar,ywalp(ibopo),rough_dom,u,nu,ustar)
              else
                 call frivel(kfl_ustar,delta_nsi,rough_dom,u,nu,ustar)
              end if
              if(u/=0.0_rp) then
                 us2 = rho*ustar*ustar/u
                 do idime=1,ndime
                    gevec(idime,ibopo) = -us2*veloc(idime,ipoin,1)
                 end do
              end if
           end if
        end do
     else
        !
        ! Up to the wall
        !
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if(ibopo>=1) then
              n1  = exnor(1,1,ibopo)
              n2  = exnor(2,1,ibopo)
              s11 = gradv_nsi(1,ipoin)      ! ( du/dx + du/dx )
              s22 = gradv_nsi(2,ipoin)      ! ( dv/dy + dv/dy )
              s12 = gradv_nsi(3,ipoin)      ! ( du/dy + dv/dx )
              sn1 = s11*n1+s12*n2
              sn2 = s12*n1+s22*n2
              if(ndime==3) then
                 n3   = exnor(3,1,ibopo) 
                 s33  = gradv_nsi(4,ipoin)  ! ( dw/dz + dw/dz )
                 s13  = gradv_nsi(5,ipoin)  ! ( du/dz + dw/dx )
                 s23  = gradv_nsi(6,ipoin)  ! ( dv/dz + dw/dy )
                 sn1  = sn1+s13*n3
                 sn2  = sn2+s23*n3
                 sn3  = s13*n1+s23*n2+s33*n3
                 gevec(1,ibopo)=sn1        
                 gevec(2,ibopo)=sn2       
                 gevec(3,ibopo)=sn3      
              else                 
                 gevec(1,ibopo)=sn1     
                 gevec(2,ibopo)=sn2       
              end if
              do idime = 1,ndime
                 gevec(idime,ibopo) = -press(ipoin,1) * exnor(idime,1,ibopo) + gevec(idime,ibopo)
              end do
           end if
        end do
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 do idime = 1,ndime
                    gevec(idime,ibopo) = -press(ipoin,1) * exnor(idime,1,ibopo) + gevec(idime,ibopo)
                 end do
              end if
           end do
        end if

     end if

  else  ! itask .ne. 0

     !-------------------------------------------------------------------
     !
     ! Tangential traction on all nodes. Set to 0 on interior nodes
     !
     !-------------------------------------------------------------------

     if( delta_nsi > 0.0_rp .or. kfl_delta == 1 ) then
        !
        ! Wall law (uses entire velocity, not tangent)
        !
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo >= 1 ) then
              rho = prope_nsi(1,ipoin)
              nu  = prope_nsi(2,ipoin)/rho
              u   = 0.0_rp
              do idime=1,ndime
                 u = u + veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
              end do
              u = sqrt(u)
              if( kfl_rough >  0 ) rough_dom = rough(ipoin) 
              if( kfl_delta == 1 ) then
                 call frivel(kfl_ustar,ywalp(ibopo),rough_dom,u,nu,ustar)
              else
                 call frivel(kfl_ustar,delta_nsi,rough_dom,u,nu,ustar)
              end if
              if( u /= 0.0_rp ) then
                 us2 = rho*ustar*ustar/u
                 do idime=1,ndime  ! it should be tangent velocity
                    gevec(idime,ipoin) = -us2*veloc(idime,ipoin,1)
                 end do
              end if
           end if
        end do

     else
        !
        ! Up to the wall
        !
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo >= 1 ) then
              n1  = exnor(1,1,ibopo)
              n2  = exnor(2,1,ibopo)
              s11 = gradv_nsi(1,ipoin)
              s22 = gradv_nsi(2,ipoin)
              s12 = gradv_nsi(3,ipoin)
              sn1 = s11*n1+s12*n2
              sn2 = s12*n1+s22*n2
              if( ndime == 3 ) then
                 n3   = exnor(3,1,ibopo) 
                 s33  = gradv_nsi(4,ipoin)
                 s13  = gradv_nsi(5,ipoin)
                 s23  = gradv_nsi(6,ipoin)
                 sn1  = sn1+s13*n3
                 sn2  = sn2+s23*n3
                 sn3  = s13*n1+s23*n2+s33*n3
                 gevec(1,ipoin)=sn1        
                 gevec(2,ipoin)=sn2       
                 gevec(3,ipoin)=sn3      
              else
                 gevec(1,ipoin)=sn1         
                 gevec(2,ipoin)=sn2       
              end if
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then
              do idime = 1,ndime
                 gevec(idime,ipoin) = -press(ipoin,1) * exnor(idime,1,ibopo) + gevec(idime,ipoin)
                 gevec(idime,ipoin) = -gevec(idime,ipoin)
              end do
           end if
        end do
     end if
  end if
  !
  ! Deallocate memory
  !
  if( delta_nsi > 0.0_rp .or. kfl_delta == 1 ) then        
  else
     call memchk(two,istat,metmp,'GRADV_NSI','nsi_outtra',gradv_nsi)
     deallocate(gradv_nsi,stat=istat)
     if(istat/=0) call memerr(two,'GRADV_NSI','nsi_outtra',0_ip)
  end if

end subroutine nsi_outtra
