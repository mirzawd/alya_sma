!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine bouder_der(pnode,nnodb,ndime,ndimb,deriv,bocod,elcod,baloc_der,eucta_der)
  !-----------------------------------------------------------------------
  !****f* Domain/bouder_der
  ! NAME
  !    bouder_der
  ! DESCRIPTION
  !    This routine calculates the baloc and eucta derivatives w.r.t. coordinates
  ! USES
  !    vecnor
  !    vecpro
  ! USED BY
  !    nsm_bouope
  !    tem_bouope
  ! SOURCE
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_kermod, only     :  kfl_adj_prob
  implicit none
  integer(ip), intent(in)  :: nnodb,ndime,ndimb,pnode
  real(rp),    intent(in)  :: deriv(max(1_ip,ndimb),nnodb)
  real(rp),    intent(in)  :: bocod(ndime,nnodb),elcod(ndime,pnode)
  real(rp),    intent(out) :: eucta_der(ndime,nnodb),baloc_der(ndime,ndime,ndime,nnodb)
  
  real(rp)                    eucta,baloc(ndime,ndime),produ,cocog(ndime),dummr,produ_der(ndime,nnodb)
  real(rp)                    baloc_copy(ndime,ndime),baloc_copy_der(ndime,ndime,ndime,nnodb),baloc_copy0(ndime,ndime)
  real(rp)                    baloc_copy0_der(ndime,ndime,ndime,nnodb)
  integer(ip)                 idime,inode
  
  if (kfl_adj_prob == 0) return
  
  baloc_der = 0.0_rp
  eucta_der = 0.0_rp
  baloc_copy = 0.0_rp
  baloc_copy0 = 0.0_rp
  baloc_copy_der = 0.0_rp
  baloc_copy0_der = 0.0_rp

  dummr=1.0_rp/real(pnode,rp)
  
  if( ndime == 1 ) then
     !
     ! 1D
     !
!      baloc(1,1) =   1.0_rp
!      eucta      =   baloc(1,1)*baloc(1,1)
!      eucta      =   sqrt(eucta)
     
  else if( ndime == 2 ) then
     !
     ! 2D
     !
     call mbmabt(baloc_copy,bocod,deriv,ndime,ndimb,nnodb)  ! Evaluates the tangent vectors
     
     baloc_copy(1,2) =   baloc_copy(2,1)
     baloc_copy(2,2) = - baloc_copy(1,1)
     
     eucta      =   baloc_copy(1,2) * baloc_copy(1,2) &
          &       + baloc_copy(2,2) * baloc_copy(2,2)
     eucta      =   sqrt(eucta)
     
     baloc_copy_der(1,1,1,1) = deriv(1,1)
     baloc_copy_der(1,1,1,2) = deriv(1,2)
     
     baloc_copy_der(2,1,2,1) = deriv(1,1)
     baloc_copy_der(2,1,2,2) = deriv(1,2)

     baloc_copy_der(1,2,2,1) = baloc_copy_der(2,1,2,1)
     baloc_copy_der(1,2,2,2) = baloc_copy_der(2,1,2,2)
     
     baloc_copy_der(2,2,1,1) = -baloc_copy_der(1,1,1,1)
     baloc_copy_der(2,2,1,2) = -baloc_copy_der(1,1,1,2)
     
     do idime = 1, ndime
       do inode = 1, nnodb
         eucta_der(idime,inode) = 0.5_rp*(2.0_rp*baloc_copy_der(1,2,idime,inode)*baloc_copy(1,2) + &
                                          2.0_rp*baloc_copy_der(2,2,idime,inode)*baloc_copy(2,2))/eucta
       enddo
     enddo
     
     do idime = 1, ndime
       do inode = 1, nnodb
         produ_der(idime,inode) = -0.5_rp*(2.0_rp*baloc_copy_der(1,2,idime,inode)*baloc_copy(1,2) + &
                                          2.0_rp*baloc_copy_der(2,2,idime,inode)*baloc_copy(2,2))/(eucta*eucta*eucta)
       enddo
     enddo
         
     
     cocog(1)=0.0_rp
     cocog(2)=0.0_rp
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
        cocog(2)=cocog(2)+elcod(2,inode)
     end do
     cocog(1)=cocog(1)*dummr
     cocog(2)=cocog(2)*dummr

     produ=      baloc_copy(1,1)*baloc_copy(1,1)
     produ=produ+baloc_copy(2,1)*baloc_copy(2,1)     
     
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,1)=produ*baloc_copy(1,1)
        baloc(2,1)=produ*baloc_copy(2,1)
        baloc(1,2)=produ*baloc_copy(1,2)
        baloc(2,2)=produ*baloc_copy(2,2)
     end if

     if(produ/=0.0_rp) then
        do idime = 1, ndime
          do inode = 1, nnodb
            baloc_der(1,1,idime,inode) = produ_der(idime,inode)*baloc_copy(1,1) + produ*baloc_copy_der(1,1,idime,inode)
            baloc_der(2,1,idime,inode) = produ_der(idime,inode)*baloc_copy(2,1) + produ*baloc_copy_der(2,1,idime,inode)
            baloc_der(1,2,idime,inode) = produ_der(idime,inode)*baloc_copy(1,2) + produ*baloc_copy_der(1,2,idime,inode)
            baloc_der(2,2,idime,inode) = produ_der(idime,inode)*baloc_copy(2,2) + produ*baloc_copy_der(2,2,idime,inode)
          enddo
        enddo
     end if

     produ=(cocog(1)-bocod(1,1))*baloc(1,2)&
          +(cocog(2)-bocod(2,1))*baloc(2,2)

     if(produ>0.0_rp) then
       do idime = 1, ndime
         do inode = 1, nnodb
           baloc_der(1,1,idime,inode) = -baloc_der(1,1,idime,inode)     ! t1=-t1                      
           baloc_der(1,2,idime,inode) = -baloc_der(1,2,idime,inode)     ! n =-n
           baloc_der(2,1,idime,inode) = -baloc_der(2,1,idime,inode)     ! t1=-t1                      
           baloc_der(2,2,idime,inode) = -baloc_der(2,2,idime,inode)     ! n =-n
         enddo
       enddo
     end if
     
  else if( ndime == 3 ) then
     !
     ! 3D
     !
     call mbmabt(baloc_copy0,bocod,deriv,ndime,ndimb,nnodb)  ! Evaluates the tangent vectors
     
     do inode = 1, nnodb  
       baloc_copy0_der(1,1,1,inode) = deriv(1,inode) !dbaloc_copy0(1,1)/dx_inode
       baloc_copy0_der(1,1,2,inode) = 0.0_rp         !dbaloc_copy0(1,1)/dy_inode
       baloc_copy0_der(1,1,3,inode) = 0.0_rp         !dbaloc_copy0(1,1)/dz_inode
       
       baloc_copy0_der(1,2,1,inode) = deriv(2,inode) !dbaloc_copy0(1,2)/dx_inode
       baloc_copy0_der(1,2,2,inode) = 0.0_rp         !dbaloc_copy0(1,2)/dy_inode
       baloc_copy0_der(1,2,3,inode) = 0.0_rp         !dbaloc_copy0(1,2)/dz_inode
       
       baloc_copy0_der(2,1,1,inode) = 0.0_rp         !dbaloc_copy0(2,1)/dx_inode
       baloc_copy0_der(2,1,2,inode) = deriv(1,inode) !dbaloc_copy0(2,1)/dy_inode
       baloc_copy0_der(2,1,3,inode) = 0.0_rp         !dbaloc_copy0(2,1)/dz_inode
       
       baloc_copy0_der(2,2,1,inode) = 0.0_rp         !dbaloc_copy0(2,2)/dx_inode
       baloc_copy0_der(2,2,2,inode) = deriv(2,inode) !dbaloc_copy0(2,2)/dy_inode
       baloc_copy0_der(2,2,3,inode) = 0.0_rp         !dbaloc_copy0(2,2)/dz_inode
       
       baloc_copy0_der(3,1,1,inode) = 0.0_rp         !dbaloc_copy0(3,1)/dx_inode
       baloc_copy0_der(3,1,2,inode) = 0.0_rp         !dbaloc_copy0(3,1)/dy_inode
       baloc_copy0_der(3,1,3,inode) = deriv(1,inode) !dbaloc_copy0(3,1)/dz_inode
       
       baloc_copy0_der(3,2,1,inode) = 0.0_rp         !dbaloc_copy0(3,2)/dx_inode
       baloc_copy0_der(3,2,2,inode) = 0.0_rp         !dbaloc_copy0(3,2)/dy_inode
       baloc_copy0_der(3,2,3,inode) = deriv(2,inode) !dbaloc_copy0(3,2)/dz_inode  
       
     enddo
     
     baloc_copy0(1,3) =   baloc_copy0(2,1) * baloc_copy0(3,2) - baloc_copy0(3,1) * baloc_copy0(2,2)
     baloc_copy0(2,3) =   baloc_copy0(3,1) * baloc_copy0(1,2) - baloc_copy0(1,1) * baloc_copy0(3,2)
     baloc_copy0(3,3) =   baloc_copy0(1,1) * baloc_copy0(2,2) - baloc_copy0(2,1) * baloc_copy0(1,2)
          
     do idime = 1, ndime
       do inode = 1, nnodb
         baloc_copy0_der(1,3,idime,inode) = baloc_copy0_der(2,1,idime,inode) * baloc_copy0(3,2) - & !dbaloc_copy0(1,3)/dx_inode
                                            baloc_copy0_der(3,1,idime,inode) * baloc_copy0(2,2) + &
                                            baloc_copy0(2,1) * baloc_copy0_der(3,2,idime,inode) - &
                                            baloc_copy0(3,1) * baloc_copy0_der(2,2,idime,inode)
                                            
         baloc_copy0_der(2,3,idime,inode) = baloc_copy0_der(3,1,idime,inode) * baloc_copy0(1,2) - & !dbaloc_copy0(2,3)/dx_inode
                                            baloc_copy0_der(1,1,idime,inode) * baloc_copy0(3,2) + &
                                            baloc_copy0(3,1) * baloc_copy0_der(1,2,idime,inode) - &
                                            baloc_copy0(1,1) * baloc_copy0_der(3,2,idime,inode)
                                            
         baloc_copy0_der(3,3,idime,inode) = baloc_copy0_der(1,1,idime,inode) * baloc_copy0(2,2) - & !dbaloc_copy0(3,3)/dx_inode
                                            baloc_copy0_der(2,1,idime,inode) * baloc_copy0(1,2) + &
                                            baloc_copy0(1,1) * baloc_copy0_der(2,2,idime,inode) - &
                                            baloc_copy0(2,1) * baloc_copy0_der(1,2,idime,inode)
        
       enddo
     enddo
     
     eucta      =   baloc_copy0(1,3) * baloc_copy0(1,3) &
          &       + baloc_copy0(2,3) * baloc_copy0(2,3) &
          &       + baloc_copy0(3,3) * baloc_copy0(3,3)
     eucta      =   sqrt(eucta)

     do idime = 1, ndime
       do inode = 1, nnodb
         eucta_der(idime,inode)      = 0.5_rp * ( 2.0_rp*baloc_copy0_der(1,3,idime,inode) * baloc_copy0(1,3) + &
                                                  2.0_rp*baloc_copy0_der(2,3,idime,inode) * baloc_copy0(2,3) + &
                                                  2.0_rp*baloc_copy0_der(3,3,idime,inode) * baloc_copy0(3,3)   )/eucta                
       enddo
     enddo
       
     ! recalculate t1 so that it is orthogonal to t2
     baloc_copy(1,1) =   baloc_copy0(2,3) * baloc_copy0(3,2) - baloc_copy0(3,3) * baloc_copy0(2,2)
     baloc_copy(2,1) =   baloc_copy0(3,3) * baloc_copy0(1,2) - baloc_copy0(1,3) * baloc_copy0(3,2)
     baloc_copy(3,1) =   baloc_copy0(1,3) * baloc_copy0(2,2) - baloc_copy0(2,3) * baloc_copy0(1,2)
     
     baloc_copy(1,2) = baloc_copy0(1,2)
     baloc_copy(2,2) = baloc_copy0(2,2)
     baloc_copy(3,2) = baloc_copy0(3,2)
     
     baloc_copy(1,3) = baloc_copy0(1,3)
     baloc_copy(2,3) = baloc_copy0(2,3)
     baloc_copy(3,3) = baloc_copy0(3,3)
     
     do idime = 1, ndime
       do inode = 1, nnodb
         baloc_copy_der(1,1,idime,inode) = baloc_copy0_der(2,3,idime,inode) * baloc_copy0(3,2) - &
                                           baloc_copy0_der(3,3,idime,inode) * baloc_copy0(2,2) + &
                                           baloc_copy0(2,3) * baloc_copy0_der(3,2,idime,inode) - & 
                                           baloc_copy0(3,3) * baloc_copy0_der(2,2,idime,inode)
                                           
         baloc_copy_der(2,1,idime,inode) = baloc_copy0_der(3,3,idime,inode) * baloc_copy0(1,2) - &
                                           baloc_copy0_der(1,3,idime,inode) * baloc_copy0(3,2) + &
                                           baloc_copy0(3,3) * baloc_copy0_der(1,2,idime,inode) - &
                                           baloc_copy0(1,3) * baloc_copy0_der(3,2,idime,inode)
                                           
         baloc_copy_der(3,1,idime,inode) = baloc_copy0_der(1,3,idime,inode) * baloc_copy0(2,2) - &
                                           baloc_copy0_der(2,3,idime,inode) * baloc_copy0(1,2) + &
                                           baloc_copy0(1,3) * baloc_copy0_der(2,2,idime,inode) - &
                                           baloc_copy0(2,3) * baloc_copy0_der(1,2,idime,inode)
                                           
         baloc_copy_der(1,2,idime,inode) = baloc_copy0_der(1,2,idime,inode)
         baloc_copy_der(2,2,idime,inode) = baloc_copy0_der(2,2,idime,inode)
         baloc_copy_der(3,2,idime,inode) = baloc_copy0_der(3,2,idime,inode)
     
         baloc_copy_der(1,3,idime,inode) = baloc_copy0_der(1,3,idime,inode)
         baloc_copy_der(2,3,idime,inode) = baloc_copy0_der(2,3,idime,inode)
         baloc_copy_der(3,3,idime,inode) = baloc_copy0_der(3,3,idime,inode)
                                           
       enddo
     enddo
     
     
     cocog(1)=0.0_rp
     cocog(2)=0.0_rp
     cocog(3)=0.0_rp
     do inode=1,pnode
        cocog(1)=cocog(1)+elcod(1,inode)
        cocog(2)=cocog(2)+elcod(2,inode)
        cocog(3)=cocog(3)+elcod(3,inode)
     end do
     cocog(1)=cocog(1)*dummr
     cocog(2)=cocog(2)*dummr
     cocog(3)=cocog(3)*dummr

     !!!!!!!!!!!!!!!!!!!!!!
     
     produ =       baloc_copy(1,1)*baloc_copy(1,1)
     produ = produ+baloc_copy(2,1)*baloc_copy(2,1)
     produ = produ+baloc_copy(3,1)*baloc_copy(3,1)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,1)=produ*baloc_copy(1,1)
        baloc(2,1)=produ*baloc_copy(2,1)
        baloc(3,1)=produ*baloc_copy(3,1)
     end if
     
     do idime = 1, ndime
       do inode = 1, nnodb
         produ_der(idime,inode) = -0.5_rp*(2.0_rp*baloc_copy_der(1,1,idime,inode)*baloc_copy(1,1) + &
                                           2.0_rp*baloc_copy_der(2,1,idime,inode)*baloc_copy(2,1) + &
                                           2.0_rp*baloc_copy_der(3,1,idime,inode)*baloc_copy(3,1))*(produ*produ*produ)
       enddo
     enddo
     if(produ/=0.0_rp) then
        do idime = 1, ndime
          do inode = 1, nnodb
            baloc_der(1,1,idime,inode) = produ_der(idime,inode)*baloc_copy(1,1) + produ*baloc_copy_der(1,1,idime,inode)
            baloc_der(2,1,idime,inode) = produ_der(idime,inode)*baloc_copy(2,1) + produ*baloc_copy_der(2,1,idime,inode)
            baloc_der(3,1,idime,inode) = produ_der(idime,inode)*baloc_copy(3,1) + produ*baloc_copy_der(3,1,idime,inode)
          enddo
        enddo
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!
     
     produ =       baloc_copy(1,2)*baloc_copy(1,2)
     produ = produ+baloc_copy(2,2)*baloc_copy(2,2)
     produ = produ+baloc_copy(3,2)*baloc_copy(3,2)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,2)=produ*baloc_copy(1,2)
        baloc(2,2)=produ*baloc_copy(2,2)
        baloc(3,2)=produ*baloc_copy(3,2)
     end if
     
     do idime = 1, ndime
       do inode = 1, nnodb
         produ_der(idime,inode) = -0.5_rp*(2.0_rp*baloc_copy_der(1,2,idime,inode)*baloc_copy(1,2) + &
                                           2.0_rp*baloc_copy_der(2,2,idime,inode)*baloc_copy(2,2) + &
                                           2.0_rp*baloc_copy_der(3,2,idime,inode)*baloc_copy(3,2))*(produ*produ*produ)
       enddo
     enddo
     if(produ/=0.0_rp) then
        do idime = 1, ndime
          do inode = 1, nnodb
            baloc_der(1,2,idime,inode) = produ_der(idime,inode)*baloc_copy(1,2) + produ*baloc_copy_der(1,2,idime,inode)
            baloc_der(2,2,idime,inode) = produ_der(idime,inode)*baloc_copy(2,2) + produ*baloc_copy_der(2,2,idime,inode)
            baloc_der(3,2,idime,inode) = produ_der(idime,inode)*baloc_copy(3,2) + produ*baloc_copy_der(3,2,idime,inode)
          enddo
        enddo
     end if
      
     !!!!!!!!!!!!!!!!!!!!!!!
     
     produ =       baloc_copy(1,3)*baloc_copy(1,3)
     produ = produ+baloc_copy(2,3)*baloc_copy(2,3)
     produ = produ+baloc_copy(3,3)*baloc_copy(3,3)
     if(produ/=0.0_rp) then
        produ=1.0_rp/sqrt(produ)
        baloc(1,3)=produ*baloc_copy(1,3)
        baloc(2,3)=produ*baloc_copy(2,3)
        baloc(3,3)=produ*baloc_copy(3,3)
     end if
     
     do idime = 1, ndime
       do inode = 1, nnodb
         produ_der(idime,inode) = -0.5_rp*(2.0_rp*baloc_copy_der(1,3,idime,inode)*baloc_copy(1,3) + &
                                           2.0_rp*baloc_copy_der(2,3,idime,inode)*baloc_copy(2,3) + &
                                           2.0_rp*baloc_copy_der(3,3,idime,inode)*baloc_copy(3,3))*(produ*produ*produ)
       enddo
     enddo
     if(produ/=0.0_rp) then
        do idime = 1, ndime
          do inode = 1, nnodb
            baloc_der(1,3,idime,inode) = produ_der(idime,inode)*baloc_copy(1,3) + produ*baloc_copy_der(1,3,idime,inode)
            baloc_der(2,3,idime,inode) = produ_der(idime,inode)*baloc_copy(2,3) + produ*baloc_copy_der(2,3,idime,inode)
            baloc_der(3,3,idime,inode) = produ_der(idime,inode)*baloc_copy(3,3) + produ*baloc_copy_der(3,3,idime,inode)
          enddo
        enddo
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!!!!
     
     produ=(cocog(1)-bocod(1,1))*baloc(1,3)&
          +(cocog(2)-bocod(2,1))*baloc(2,3)&
          +(cocog(3)-bocod(3,1))*baloc(3,3)

     if(produ>0.0_rp) then
       do idime = 1, ndime
         do inode = 1, nnodb
        baloc_der(1,1,idime,inode) = -baloc_der(1,1,idime,inode)     ! t1=-t1                      
        baloc_der(1,3,idime,inode) = -baloc_der(1,3,idime,inode)     ! n =-n
        baloc_der(2,1,idime,inode) = -baloc_der(2,1,idime,inode)     ! t1=-t1                      
        baloc_der(2,3,idime,inode) = -baloc_der(2,3,idime,inode)     ! n =-n
        baloc_der(3,1,idime,inode) = -baloc_der(3,1,idime,inode)     ! t1=-t1                      
        baloc_der(3,3,idime,inode) = -baloc_der(3,3,idime,inode)     ! n =-n
         enddo
       enddo
     end if

  end if

end subroutine bouder_der

