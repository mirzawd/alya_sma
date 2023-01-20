!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_boumat(&
     pnode,pnodb, igaub,iboun,lboel,lelbo,xmmat,xmrhs,gbsha,&
     gbsur,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* temper/tem_boumat
  ! NAME 
  !    tem_boumat
  ! DESCRIPTION
  !    Assemble boundary contribution
  ! USES
  ! USED BY
  !    tem_bouope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ltype, lnods
  use def_kermod, only       :  kfl_waexl_ker, kfl_waexl_imp_ker, lexlo_ker, shape_waexl
  use def_master, only       :  tempe
!  use def_temper, only       :  kfl_fixbo_tem
  
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb, igaub, iboun
  integer(ip), intent(in)    :: lboel(pnodb)
  integer(ip), intent(in)    :: lelbo
  real(rp),    intent(in)    :: xmmat,xmrhs,gbsha(pnodb),gbsur
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inodb,jnodb,inode,jnode, ielem, pelty, index
  real(rp)                   :: xmuit, wetem

  if ( kfl_waexl_ker == 0.or.kfl_waexl_imp_ker==0) then 
     do inodb=1,pnodb  
        inode=lboel(inodb)
        elrhs(inode)=elrhs(inode)+gbsha(inodb)*xmrhs*gbsur
        xmuit=xmmat*gbsha(inodb)*gbsur
        do jnodb=1,pnodb
           jnode=lboel(jnodb)
           elmat(inode,jnode)=elmat(inode,jnode)&
                +xmuit*gbsha(jnodb)
        end do
     end do
  else      ! WALL EXCHANGHE LOCATION implicit 
     ielem = lelbo
     pelty = ltype(ielem)
     index = lexlo_ker(igaub,iboun)
     do jnode =1, pnode                
        xmuit=xmmat* shape_waexl(jnode, index)*gbsur
        wetem = wetem +  shape_waexl(jnode, index)*tempe(lnods(jnode, ielem), 1)
        do inodb = 1,pnodb
           inode=lboel(inodb) 
           elmat(inode,jnode)=elmat(inode,jnode)&
                +xmuit*gbsha(inodb)           
        end do
     end do
  end if



!! AB check   if (  kfl_waexl_ker == 1.and.kfl_waexl_imp_ker==1.and. kfl_fixbo_tem(iboun)==3) then 
!! AB check      ! WALL EXCHANGHE LOCATION implicit, ensamble h in matrix (nonsymmetric)
!! AB check      ielem = lboel(pnodb+1)
!! AB check      index = lexlo_ker(igaub,iboun)
!! AB check      do jnode =1, pnode                
!! AB check         xmuit=xmmat* shape_waexl(jnode, index)*gbsur
!! AB check         do inodb = 1,pnodb
!! AB check            inode=lboel(inodb) 
!! AB check            elmat(inode,jnode)=elmat(inode,jnode)&
!! AB check                 +xmuit*gbsha(inodb)           
!! AB check         end do
!! AB check      end do
!! AB check      ! RHS
!! AB check      do inodb =1, pnodb
!! AB check         inode=lboel(inodb)
!! AB check         elrhs(inode)=elrhs(inode)+gbsha(inodb)*xmrhs*gbsur
!! AB check      end do
!! AB check 
!! AB check   else ! classical
!! AB check      do inodb=1,pnodb  
!! AB check         inode=lboel(inodb)
!! AB check         elrhs(inode)=elrhs(inode)+gbsha(inodb)*xmrhs*gbsur
!! AB check         xmuit=xmmat*gbsha(inodb)*gbsur
!! AB check         do jnodb=1,pnodb
!! AB check            jnode=lboel(jnodb)
!! AB check            elmat(inode,jnode)=elmat(inode,jnode)&
!! AB check                 +xmuit*gbsha(jnodb)
!! AB check         end do
!! AB check      end do
!! AB check <<<<<<< .working
!! AB check =======
!! AB check   else      ! WALL EXCHANGHE LOCATION implicit 
!! AB check      ielem = lelbo
!! AB check      pelty = ltype(ielem)
!! AB check      index = lexlo_ker(igaub,iboun)
!! AB check      do jnode =1, pnode                
!! AB check         xmuit=xmmat* shape_waexl(jnode, index)*gbsur
!! AB check         wetem = wetem +  shape_waexl(jnode, index)*tempe(lnods(jnode, ielem), 1)
!! AB check         do inodb = 1,pnodb
!! AB check            inode=lboel(inodb) 
!! AB check            elmat(inode,jnode)=elmat(inode,jnode)&
!! AB check                 +xmuit*gbsha(inodb)           
!! AB check         end do
!! AB check      end do
!! AB check >>>>>>> .merge-right.r12343
!! AB check   end if

end subroutine tem_boumat
