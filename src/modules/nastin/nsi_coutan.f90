!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_coutan()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_coutan
  ! NAME 
  !    nsi_coutan
  ! DESCRIPTION
  !    This routine computes the tangential component of the traction
  !    for the RANS/LES two layer wall model
  ! USES
  ! USED BY
  !    nsi_plugin
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_memory
  use def_kermod
  use mod_gradie
  use mod_ker_proper,     only : ker_proper

  implicit none
  integer(ip)             :: ipoin,ibopo,idime,dummi
  real(rp)                :: n1,n2,n3,s11,s12,s22,s13,s23,s33,stnn
  real(rp)                :: sn1,sn2,sn3,stx,sty,stz
  real(rp),    pointer    :: visco_tmp(:)

  nullify(  visco_tmp )
  
  if( INOTEMPTY ) then
     ! 
     ! Velocity gradients GRADV_NSI
     ! 
     call memory_alloca(mem_modul(1:2,modul),'GRADV_NSI','nsi_coutan',gradv_nsi,ntens,npoin)
     call gradie(veloc(1:ndime,1:npoin,1),gradv_nsi)
     call memory_alloca(mem_modul(1:2,modul),'VISCO_TMP','nsi_coutan',visco_tmp,npoin)
     call ker_proper('VISCO','NPOIN',dummi,dummi,visco_tmp)
     do ipoin = 1,npoin  
        gradv_nsi(1:ntens,ipoin) = gradv_nsi(1:ntens,ipoin) * visco_tmp(ipoin) 
     end do
     call memory_deallo(mem_modul(1:2,modul),'VISCO_TMP','nsi_coutan',visco_tmp)

     !-------------------------------------------------------------------
     !
     ! Tangential traction on all nodes. Set to 0 on interior nodes
     !
     !-------------------------------------------------------------------

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
              stnn = sn1*n1+sn2*n2+sn3*n3
              stx  = sn1-stnn*n1
              sty  = sn2-stnn*n2
              stz  = sn3-stnn*n3
              tracr_nsi(1,ipoin)=stx         
              tracr_nsi(2,ipoin)=sty       
              tracr_nsi(3,ipoin)=stz      
           else
              stnn = sn1*n1+sn2*n2
              stx  = sn1-stnn*n1
              sty  = sn2-stnn*n2
              tracr_nsi(1,ipoin)=stx         
              tracr_nsi(2,ipoin)=sty       
           end if
        else
           do idime = 1,ndime
              tracr_nsi(idime,ipoin) = 0.0_rp
           end do
        end if
     end do

     !
     ! Deallocate memory
     !
     call memory_deallo(mem_modul(1:2,modul),'GRADV_NSI','nsi_coutan',gradv_nsi)

  end if

end subroutine nsi_coutan

!------------------------------------------------------------------------ 
!
! - Velocity strain rates gravb(ntens,npoin)
!   gradv_nsi(1,ipoin)=mu*(du/dx+du/dx)     
!   gradv_nsi(2,ipoin)=mu*(dv/dy+dv/dy)     
!   gradv_nsi(3,ipoin)=mu*(du/dy+dv/dx)     
!   gradv_nsi(4,ipoin)=mu*(dw/dz+dw/dz)     
!   gradv_nsi(5,ipoin)=mu*(du/dz+dw/dz)     
!   gradv_nsi(6,ipoin)=mu*(dv/dz+dw/dy)
! 
! In the case of the law of the wall, on walls:
! 
! (sig.n).t = rho*(U*)^2
! 
! Otherwise:      
!  
! Let sij=1/2(ui,j+uj,i), the traction is given by      
!         +-                                               -+
!         | -p*n1 + 2*mu*s11*n1 + 2*mu*s12*n2 + 2*mu*s13*n3 |
! sig.n = | -p*n2 + 2*mu*s12*n1 + 2*mu*s22*n2 + 2*mu*s23*n3 |
!         | -p*n3 + 2*mu*s13*n1 + 2*mu*s23*n2 + 2*mu*s33*n3 |
!         +-                                               -+
!       = (sn1,sn2,sn3)
! 
! The tangential component of the traction is computed as:
! (sig.n).t = sig.n-(n.sig.n)n
!           = (sn1,sn2,sn3)-(sn1*n1+sn2*n2+sn3*n3)(n1,n2,n3)
!             +-                                       -+
!             | sn1 - sn1*n1*n1 - sn2*n2*n1 - sn3*n3*n1 |
!           = | sn2 - sn1*n1*n2 - sn2*n2*n2 - sn3*n3*n2 |
!             | sn3 - sn1*n1*n3 - sn2*n2*n3 - sn3*n3*n3 |
!             +-                                       -+
! 
! NB: the pressure does not intervene in the tangential stress
!     and grave(j,i)=mu*(ui,j+uj,i) so that sni=grave(j,i)*nj      
!***
!-----------------------------------------------------------------------
