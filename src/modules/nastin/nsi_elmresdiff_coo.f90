!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_elmresdiff_coo.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute the derivation of GPRES r. t. design variables
!
!-------------------------------------------------------------------
!        +-           -+    +-                      -+   +- -+   +-  -+
!        | elmresudiff  |   | elauudiff    elaupdiff |   | u |   | elrbudiff |
!        |              | = |                        | * |   | - |           |
!        | elmresbdiff  |   | elapudiff    elappdiff |   | p |   | elrbpdiff |
!        +-           -+    +-                      -+   +- -+   +-  -+
!-------------------------------------------------------------------
!
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_elmresdiff_coo(pnode,pgaus,pevat,gpvol,gpvol_der,gpsha,lnods, &
                              p1vec,p2vec,p2sca,elvel,wgrgr,h,h_der, &
                              gpsp1,gpsp2,rmomu,gpden,gpadv,gpvis,gpcar,gpcar_der,rcont,elresucoodiff,elrespcoodiff,ielem)

  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
!  use def_nastin, only     :  dtinv_nsi, pabdf_nsi
  use def_master

  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: p1vec(pnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime, pnode)
  real(rp),    intent(in)    :: p2vec(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: p2sca(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol_der(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: wgrgr(pnode,pnode,pgaus)
  real(rp),    intent(in)    :: h
  real(rp),    intent(in)    :: h_der(ndime,pnode)
  real(rp),    intent(in)    :: gpsp1(pgaus)                      ! tau1'
  real(rp),    intent(in)    :: gpsp2(pgaus)                      ! tau2'
  real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpcar_der(ndime,mnode,ndime,mnode,pgaus)
  real(rp),    intent(in)    :: rcont(ndime,pnode,pgaus)
  real(rp),    intent(out)   :: elresucoodiff(ndime,pnode,ndime,pnode)
  real(rp),    intent(out)   :: elrespcoodiff(pnode,ndime,pnode)

  integer(ip)                :: idime,inode,jnode,igaus,jdime,idofv,ipoin
  integer(ip)                :: knode
  integer(ip)                :: idof1,idof2,idof3,jdof1,jdof2,jdof3,ind
  real(rp)                   :: elauudiff(pnode*ndime,pnode*ndime,ndime,pnode)
  real(rp)                   :: elaupdiff(pnode*ndime,pnode,ndime,pnode)
  real(rp)                   :: elappdiff(pnode,pnode,ndime,pnode)
  real(rp)                   :: elapudiff(pnode,pnode*ndime,ndime,pnode)
  real(rp)                   :: elrbudiff(ndime,pnode,ndime,pnode)
  real(rp)	             :: elrbpdiff(pnode,ndime,pnode)
  real(rp)                   :: elraudiff(ndime,pnode,ndime,pnode)
  real(rp)	             :: elrapdiff(pnode,ndime,pnode)
  real(rp)	             :: fact0, fact1,fact2,fact3,fact4,fact6,elvel1(ndime*pnode),elpre1(pnode)
  real(rp)                   :: auuprodu(ndime*pnode,ndime,pnode),aupprodp(ndime*pnode,ndime,pnode)
  real(rp)                   :: apuprodu(pnode,ndime,pnode),appprodp(pnode,ndime,pnode)
  real(rp)                   :: gpsp1diff(ndime,pnode,pgaus)                      ! d(tau1')/d(x)
  real(rp)    		     :: gpsp2diff(ndime,pnode,pgaus)                      ! d(tau2')/d(x)
  real(rp)                   :: p1vecdiff(ndime,pnode,pnode,pgaus)
  real(rp)                   :: p2vecdiff(ndime,pnode,ndime,pnode,pgaus)
  real(rp)                   :: rmomu1(pnode,pgaus),gpvno,rmomu1_der(ndime,pnode,pnode,pgaus)
  real(rp)                   :: rcont_der(ndime,mnode,ndime,mnode,pgaus),p2sca_der(ndime,pnode,pnode,pgaus)
  real(rp)                   :: wgrgr_der(ndime,pnode,pnode,pnode,pgaus)
  real(rp)                   :: fact1_der,fact2_der,fact3_der,fact3_de1,fact4_der,fact6_der,fact0_der
  
  
  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  elrbpdiff = 0.0_rp
  elrespcoodiff = 0.0_rp
  elrbudiff = 0.0_rp
  elresucoodiff = 0.0_rp
  elappdiff = 0.0_rp
  elauudiff = 0.0_rp
  elaupdiff = 0.0_rp
  elapudiff = 0.0_rp  
  apuprodu = 0.0_rp
  appprodp = 0.0_rp
  auuprodu = 0.0_rp
  aupprodp = 0.0_rp
  rmomu1_der = 0.0_rp
  
  !
  ! elvel ---> elvel1 
  !
  do inode = 1,pnode
    do idime = 1,ndime
      ind = (inode-1)*ndime + idime
      elvel1(ind) = elvel(idime,inode)
    enddo
  enddo
  !
  ! press_nsi ---> elpre1 
  !
  do inode = 1,pnode
    ipoin = lnods(inode)
    elpre1(inode) = press_forw(ipoin,1)
  end do
  !
  !  gpsp1diff and gpsp2diff
  !
  do igaus = 1,pgaus
    call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
    do idime = 1,ndime
      do inode = 1,pnode
        gpsp1diff(idime,inode,igaus) = h_der(idime,inode)/h*( 4.0_rp*gpvis(igaus)/(h*h) + 1/gpsp1(igaus) )*gpsp1(igaus)*gpsp1(igaus)
        gpsp2diff(idime,inode,igaus) = 2.0_rp*gpden(igaus)*gpvno*h_der(idime,inode)
      enddo
    enddo
  enddo        
 
  !
  !  rmomu --->  rmomu1
  !      
  do inode = 1,pnode
    do igaus = 1,pgaus
      rmomu1(inode,igaus) = rmomu(inode,igaus)
    enddo
  enddo
  
  !----------------------------------------------------------------------
  !
  ! Test functions derivatives
  !
  !----------------------------------------------------------------------
  !
  ! P1VEC =  v * [ (tau1'/tau1) - tau1' * sig ] + tau1' * rho * (uc.grad)v
  ! P2VEC =  tau1' * grad(q)  ( and * rho if in Low-Mach, but not yet)
  !
  if( ndime == 2 ) then

    do igaus = 1,pgaus
      
      do jdime = 1, ndime
        do jnode = 1,pnode
          
          fact1 = gpsp1diff(jdime,jnode,igaus)*gpden(igaus)*gpvol(igaus)
          fact2 = gpsp1(igaus)*gpden(igaus)*gpvol(igaus)
          fact3     = gpsp1(igaus) * gpvol(igaus)
          fact3_der = gpsp1diff(jdime,jnode,igaus) * gpvol(igaus)
          fact3_de1 = gpsp1(igaus) * gpvol_der(jdime,jnode,igaus)
          do inode = 1,pnode
            p1vecdiff(jdime,jnode,inode,igaus)   = fact1 *    & 
                                     ( gpadv(1,igaus) * gpcar(1,inode,igaus) &
                                     + gpadv(2,igaus) * gpcar(2,inode,igaus) ) &
                                     + fact2 * &
                                     ( gpadv(1,igaus) * gpcar_der(1,inode,jdime,jnode,igaus) &
                                     + gpadv(2,igaus) * gpcar_der(2,inode,jdime,jnode,igaus) ) &                                     
                                     + p1vec(inode,igaus)/gpvol(igaus)*gpvol_der(jdime,jnode,igaus)
            
            p2vecdiff(jdime,jnode,1,inode,igaus) = (fact3_der+fact3_de1) * gpcar(1,inode,igaus)
            p2vecdiff(jdime,jnode,2,inode,igaus) = (fact3_der+fact3_de1) * gpcar(2,inode,igaus)
            p2vecdiff(jdime,jnode,1,inode,igaus) = p2vecdiff(jdime,jnode,1,inode,igaus) + fact3 * gpcar_der(1,inode,jdime,jnode,igaus)
            p2vecdiff(jdime,jnode,2,inode,igaus) = p2vecdiff(jdime,jnode,2,inode,igaus) + fact3 * gpcar_der(2,inode,jdime,jnode,igaus)
            
            p2sca_der(jdime,jnode,inode,igaus)   = gpvol_der(jdime,jnode,igaus) * gpsha(inode,igaus)
            
            do knode = 1,pnode
              wgrgr_der(jdime,jnode,inode,knode,igaus) = &
                   &   gpcar_der(1,inode,jdime,jnode,igaus)*gpcar(1,knode,igaus) +&
                   &   gpcar(1,inode,igaus)*gpcar_der(1,knode,jdime,jnode,igaus) +&
                   &   gpcar_der(2,inode,jdime,jnode,igaus)*gpcar(2,knode,igaus) +&        
                   &   gpcar(2,inode,igaus)*gpcar_der(2,knode,jdime,jnode,igaus)             
            end do
            
          enddo !inode
        enddo !jnode
      end do !jdime
    end do !igaus

  else

    do igaus = 1,pgaus
      
      do jdime = 1, ndime
        do jnode = 1,pnode
          fact1 = gpsp1diff(jdime,jnode,igaus)*gpden(igaus)*gpvol(igaus)
          fact2 = gpsp1(igaus)*gpden(igaus)*gpvol(igaus)
          fact3     = gpsp1(igaus) * gpvol(igaus)
          fact3_der = gpsp1diff(jdime,jnode,igaus) * gpvol(igaus)
          fact3_de1 = gpsp1(igaus) * gpvol_der(jdime,jnode,igaus)
          do inode = 1,pnode
            p1vecdiff(jdime,jnode,inode,igaus)   = fact1 *    & 
                                     ( gpadv(1,igaus) * gpcar(1,inode,igaus) &
                                     + gpadv(2,igaus) * gpcar(2,inode,igaus) &
                                     + gpadv(3,igaus) * gpcar(3,inode,igaus) ) &
                                     + fact2 * &
                                     ( gpadv(1,igaus) * gpcar_der(1,inode,jdime,jnode,igaus) &
                                     + gpadv(2,igaus) * gpcar_der(2,inode,jdime,jnode,igaus) &
                                     + gpadv(3,igaus) * gpcar_der(3,inode,jdime,jnode,igaus) ) &                                     
                                     + p1vec(inode,igaus)/gpvol(igaus)*gpvol_der(jdime,jnode,igaus)                                
            
            p2vecdiff(jdime,jnode,1,inode,igaus) = (fact3_der+fact3_de1) * gpcar(1,inode,igaus)
            p2vecdiff(jdime,jnode,2,inode,igaus) = (fact3_der+fact3_de1) * gpcar(2,inode,igaus)
            p2vecdiff(jdime,jnode,3,inode,igaus) = (fact3_der+fact3_de1) * gpcar(3,inode,igaus)
            p2vecdiff(jdime,jnode,1,inode,igaus) = p2vecdiff(jdime,jnode,1,inode,igaus) + fact3 * gpcar_der(1,inode,jdime,jnode,igaus)
            p2vecdiff(jdime,jnode,2,inode,igaus) = p2vecdiff(jdime,jnode,2,inode,igaus) + fact3 * gpcar_der(2,inode,jdime,jnode,igaus)
            p2vecdiff(jdime,jnode,3,inode,igaus) = p2vecdiff(jdime,jnode,3,inode,igaus) + fact3 * gpcar_der(3,inode,jdime,jnode,igaus)
            
            p2sca_der(jdime,jnode,inode,igaus)   = gpvol_der(jdime,jnode,igaus) * gpsha(inode,igaus)
            do knode = 1,pnode
              wgrgr_der(jdime,jnode,inode,knode,igaus) = &
                   &   gpcar_der(1,inode,jdime,jnode,igaus)*gpcar(1,knode,igaus) +&
                   &   gpcar(1,inode,igaus)*gpcar_der(1,knode,jdime,jnode,igaus) +&
                   &   gpcar_der(2,inode,jdime,jnode,igaus)*gpcar(2,knode,igaus) +&        
                   &   gpcar(2,inode,igaus)*gpcar_der(2,knode,jdime,jnode,igaus) +&            
                   &   gpcar_der(3,inode,jdime,jnode,igaus)*gpcar(3,knode,igaus) +&        
                   &   gpcar(3,inode,igaus)*gpcar_der(3,knode,jdime,jnode,igaus)             
            end do
        
          enddo !inode
        enddo !jnode
      end do !jdime
    end do !igaus

  endif
  !
  ! Derivatives of RMOMU1 = RMOMU1 - rho/(dt*theta)*u (since gprhs is zero we have to cancle the time temrs from RMOMU1)
  !
  
  do igaus = 1,pgaus
    do jdime = 1, ndime
      do jnode = 1,pnode
        do inode = 1,pnode
          do idime = 1,ndime
            rmomu1_der(jdime,jnode,inode,igaus) = rmomu1_der(jdime,jnode,inode,igaus) + &
                                                  gpden(igaus) * gpadv(idime,igaus) * gpcar_der(idime,inode,jdime,jnode,igaus)
            rcont_der(jdime,jnode,idime,inode,igaus) = gpcar_der(idime,inode,jdime,jnode,igaus)
          enddo
        end do
      end do
    enddo
  end do
    
!   do igaus = 1,pgaus
!     fact1 = dtinv_nsi * pabdf_nsi(1) * gpden(igaus)
!     do inode = 1,pnode
! 	rmomu1(inode,igaus) = rmomu1(inode,igaus) - fact1 * gpsha(inode,igaus)
!     end do
!   end do

  !----------------------------------------------------------------------
  !
  ! elauudiff
  !
  !----------------------------------------------------------------------
  if( ndime == 2 ) then
    do igaus = 1,pgaus
    
      do jdime = 1,ndime
        do knode = 1,pnode    
          fact0     = gpsp2(igaus) * gpvol(igaus)  ! (1) + (5)
          fact0_der = gpsp2diff(jdime,knode,igaus) * gpvol(igaus) + gpsp2(igaus) * gpvol_der(jdime,knode,igaus)  ! (1) + (5)
          fact6     = gpvis(igaus) * gpvol(igaus)
          fact6_der = gpvis(igaus) * gpvol_der(jdime,knode,igaus)
          do inode = 1,pnode
            idof1 = 2*inode-1
            idof2 = idof1+1
            fact1     = fact0 * gpcar(1,inode,igaus)
            fact1_der = fact0_der * gpcar(1,inode,igaus) + fact0 * gpcar_der(1,inode,jdime,knode,igaus)
            fact2     = fact0 * gpcar(2,inode,igaus)
            fact2_der = fact0_der * gpcar(2,inode,igaus) + fact0 * gpcar_der(2,inode,jdime,knode,igaus)
            do jnode = 1,pnode
              jdof1              =   2*jnode-1
              jdof2              =   jdof1+1

              fact4              =   p1vec(inode,igaus) * rmomu1(jnode,igaus) &                  ! (2): ( rmomu1(u) , p1vec(v) ) 
                                   + fact6 * wgrgr(inode,jnode,igaus)                            ! (7): ( mu dui/dxj , dvi/dxj )
                                   
              fact4_der          =   p1vecdiff(jdime,knode,inode,igaus) * rmomu1(jnode,igaus) &  ! (2): ( rmomu1(u) , p1vec(v) )
                                   + p1vec(inode,igaus) * rmomu1_der(jdime,knode,jnode,igaus) &  ! (2): ( rmomu1(u) , p1vec(v) ) 
                                   + fact6_der * wgrgr(inode,jnode,igaus)                     &  ! (7): ( mu dui/dxj , dvi/dxj )
                                   + fact6 * wgrgr_der(jdime,knode,inode,jnode,igaus)            ! (7): ( mu dui/dxj , dvi/dxj )
                                   
              elauudiff(idof1,jdof1,jdime,knode) = elauudiff(idof1,jdof1,jdime,knode) + fact1_der * rcont(1,jnode,igaus) &
                                                                          + fact1 * rcont_der(jdime,knode,1,jnode,igaus) &
                                                                          + fact4_der                                       ! Auu_xx
                                                                          
              elauudiff(idof2,jdof1,jdime,knode) = elauudiff(idof2,jdof1,jdime,knode) + fact2_der * rcont(1,jnode,igaus) &
                                                                          + fact2 * rcont_der(jdime,knode,1,jnode,igaus)    ! Auu_yx
                                                                          
              elauudiff(idof1,jdof2,jdime,knode) = elauudiff(idof1,jdof2,jdime,knode) + fact1_der * rcont(2,jnode,igaus) &
                                                                          + fact1 * rcont_der(jdime,knode,2,jnode,igaus)    ! Auu_xy
                                                                          
              elauudiff(idof2,jdof2,jdime,knode) = elauudiff(idof2,jdof2,jdime,knode) + fact2_der * rcont(2,jnode,igaus) & 
                                                                          + fact2 * rcont_der(jdime,knode,2,jnode,igaus) & 
                                                                          + fact4_der                                       ! Auu_yy
            end do
          end do
          
        enddo ! knode
      enddo ! jdime

    end do
    
  else
  
    do igaus = 1,pgaus
    
      do jdime = 1,ndime
        do knode = 1,pnode    
          fact0     = gpsp2(igaus) * gpvol(igaus)  ! (1) + (5)
          fact0_der = gpsp2diff(jdime,knode,igaus) * gpvol(igaus) + gpsp2(igaus) * gpvol_der(jdime,knode,igaus)  ! (1) + (5)
          fact6     = gpvis(igaus) * gpvol(igaus)
          fact6_der = gpvis(igaus) * gpvol_der(jdime,knode,igaus)
          do inode = 1,pnode
            idof1 = 3*inode-2
            idof2 = idof1+1
            idof3 = idof2+1
            fact1     = fact0 * gpcar(1,inode,igaus)
            fact1_der = fact0_der * gpcar(1,inode,igaus) + fact0 * gpcar_der(1,inode,jdime,knode,igaus)
            fact2     = fact0 * gpcar(2,inode,igaus)
            fact2_der = fact0_der * gpcar(2,inode,igaus) + fact0 * gpcar_der(2,inode,jdime,knode,igaus)
            fact3     = fact0 * gpcar(3,inode,igaus)
            fact3_der = fact0_der * gpcar(3,inode,igaus) + fact0 * gpcar_der(3,inode,jdime,knode,igaus)
            do jnode = 1,pnode
              jdof1              =   3*jnode-2
              jdof2              =   jdof1+1
              jdof3              =   jdof2+1
              fact4              =   p1vec(inode,igaus) * rmomu1(jnode,igaus) &                  ! (2): ( rmomu1(u) , p1vec(v) ) 
                                   + fact6 * wgrgr(inode,jnode,igaus)                            ! (7): ( mu dui/dxj , dvi/dxj )
              fact4_der          =   p1vecdiff(jdime,knode,inode,igaus) * rmomu1(jnode,igaus) &  ! (2): ( rmomu1(u) , p1vec(v) )
                                   + p1vec(inode,igaus) * rmomu1_der(jdime,knode,jnode,igaus) &  ! (2): ( rmomu1(u) , p1vec(v) ) 
                                   + fact6_der * wgrgr(inode,jnode,igaus)                     &  ! (7): ( mu dui/dxj , dvi/dxj )
                                   + fact6 * wgrgr_der(jdime,knode,inode,jnode,igaus)            ! (7): ( mu dui/dxj , dvi/dxj )
                                   
              elauudiff(idof1,jdof1,jdime,knode) = elauudiff(idof1,jdof1,jdime,knode) + fact1_der * rcont(1,jnode,igaus) &
                                                                          + fact1 * rcont_der(jdime,knode,1,jnode,igaus) &
                                                                          + fact4_der                                    ! Auu_xx
                                                                          
              elauudiff(idof2,jdof1,jdime,knode) = elauudiff(idof2,jdof1,jdime,knode) + fact2_der * rcont(1,jnode,igaus) &
                                                                          + fact2 * rcont_der(jdime,knode,1,jnode,igaus) ! Auu_yx
                                                                          
              elauudiff(idof3,jdof1,jdime,knode) = elauudiff(idof3,jdof1,jdime,knode) + fact3_der * rcont(1,jnode,igaus) &
                                                                          + fact3 * rcont_der(jdime,knode,1,jnode,igaus) ! Auu_zx
                                                                          
              elauudiff(idof1,jdof2,jdime,knode) = elauudiff(idof1,jdof2,jdime,knode) + fact1_der * rcont(2,jnode,igaus) &
                                                                          + fact1 * rcont_der(jdime,knode,2,jnode,igaus) ! Auu_xy
                                                                          
              elauudiff(idof2,jdof2,jdime,knode) = elauudiff(idof2,jdof2,jdime,knode) + fact2_der * rcont(2,jnode,igaus) & 
                                                                          + fact2 * rcont_der(jdime,knode,2,jnode,igaus) & 
                                                                          + fact4_der                                    ! Auu_yy
                                                                          
              elauudiff(idof3,jdof2,jdime,knode) = elauudiff(idof3,jdof2,jdime,knode) + fact3_der * rcont(2,jnode,igaus) & 
                                                                          + fact3 * rcont_der(jdime,knode,2,jnode,igaus) ! Auu_zy
                                                                          
              elauudiff(idof1,jdof3,jdime,knode) = elauudiff(idof1,jdof3,jdime,knode) + fact1_der * rcont(3,jnode,igaus) &
                                                                          + fact1 * rcont_der(jdime,knode,3,jnode,igaus) ! Auu_xz
                                                                          
              elauudiff(idof2,jdof3,jdime,knode) = elauudiff(idof2,jdof3,jdime,knode) + fact2_der * rcont(3,jnode,igaus) & 
                                                                          + fact2 * rcont_der(jdime,knode,3,jnode,igaus) ! Auu_yz
                                                                          
              elauudiff(idof3,jdof3,jdime,knode) = elauudiff(idof3,jdof3,jdime,knode) + fact3_der * rcont(3,jnode,igaus) & 
                                                                          + fact3 * rcont_der(jdime,knode,3,jnode,igaus) & 
                                                                          + fact4_der                                    ! Auu_zz
                                                                          
            end do
          end do
          
        enddo ! knode
      enddo ! jdime

    end do

  end if
  
  
  !----------------------------------------------------------------------
  !
  ! elaupdiff
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: - ( p , div(v) ) + ( grad(p) , p1vec(v)-v )
  !  
  if( ndime == 2 ) then
    do igaus = 1,pgaus
    
      do jdime = 1,ndime
        do knode = 1,pnode    
           do jnode = 1,pnode
              fact2 = gpvol(igaus)*gpsha(jnode,igaus)
              fact2_der = gpvol_der(jdime,knode,igaus)*gpsha(jnode,igaus)
              do inode = 1,pnode
                 fact1              = -gpvol(igaus)*gpsha(inode,igaus)+p1vec(inode,igaus)
                 fact1_der          = -gpvol_der(jdime,knode,igaus)*gpsha(inode,igaus)+p1vecdiff(jdime,knode,inode,igaus)
                 
                 idofv              = (inode-1)*ndime+1
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2*gpcar_der(1,inode,jdime,knode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2_der*gpcar(1,inode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1_der*gpcar(1,jnode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1*gpcar_der(1,jnode,jdime,knode,igaus)
                 idofv              = idofv+1
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2*gpcar_der(2,inode,jdime,knode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2_der*gpcar(2,inode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1_der*gpcar(2,jnode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1*gpcar_der(2,jnode,jdime,knode,igaus)
                                  
                 
              end do !inode
           end do !jnode
        enddo ! knode
      enddo ! jdime
	  
    end do ! igaus
    
  else

    do igaus = 1,pgaus
    
      do jdime = 1,ndime
        do knode = 1,pnode    
           do jnode = 1,pnode
              fact2 = gpvol(igaus)*gpsha(jnode,igaus)
              fact2_der = gpvol_der(jdime,knode,igaus)*gpsha(jnode,igaus)
              do inode = 1,pnode
                 fact1              = -gpvol(igaus)*gpsha(inode,igaus)+p1vec(inode,igaus)
                 fact1_der          = -gpvol_der(jdime,knode,igaus)*gpsha(inode,igaus)+p1vecdiff(jdime,knode,inode,igaus)
                 
                 idofv              = (inode-1)*ndime+1
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2*gpcar_der(1,inode,jdime,knode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2_der*gpcar(1,inode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1_der*gpcar(1,jnode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1*gpcar_der(1,jnode,jdime,knode,igaus)
                 idofv              = idofv+1
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2*gpcar_der(2,inode,jdime,knode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2_der*gpcar(2,inode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1_der*gpcar(2,jnode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1*gpcar_der(2,jnode,jdime,knode,igaus)
                 idofv              = idofv+1
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2*gpcar_der(3,inode,jdime,knode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) - fact2_der*gpcar(3,inode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1_der*gpcar(3,jnode,igaus)
                 elaupdiff(idofv,jnode,jdime,knode) = elaupdiff(idofv,jnode,jdime,knode) + fact1*gpcar_der(3,jnode,jdime,knode,igaus)
                                  
                 
              end do !inode
           end do !jnode
        enddo ! knode
      enddo ! jdime
	  
    end do ! igaus
  
  
  end if
  !----------------------------------------------------------------------
  !
  ! elapudiff
  !
  !----------------------------------------------------------------------
  !
  ! ( div(u) , (tau2^{-1}*tau2')*q )   --- or ( div (rho u), q) if Low-Mach 
  ! + ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )  (only laplacian form of div[2 mu eps(u)] )
  !
  
  if( ndime == 2 ) then
    
    do igaus = 1,pgaus
      do jdime = 1,ndime
        do knode = 1,pnode    
    
          do inode = 1,pnode
            idof1 = 2*inode-1
            idof2 = idof1+1
            do jnode = 1,pnode
            
               elapudiff(jnode,idof1,jdime,knode) = elapudiff(jnode,idof1,jdime,knode) &     
                                                      + rcont(1,inode,igaus) * p2sca_der(jdime,knode,jnode,igaus) & 
                                                      + rcont_der(jdime,knode,1,inode,igaus) * p2sca(jnode,igaus) & 
                                                      + rmomu1(inode,igaus)   * p2vecdiff(jdime,knode,1,jnode,igaus) &
                                                      + rmomu1_der(jdime,knode,inode,igaus)   * p2vec(1,jnode,igaus)
                   
               elapudiff(jnode,idof2,jdime,knode) = elapudiff(jnode,idof2,jdime,knode) &
                                                      + rcont(2,inode,igaus) * p2sca_der(jdime,knode,jnode,igaus) & 
                                                      + rcont_der(jdime,knode,2,inode,igaus) * p2sca(jnode,igaus) & 
                                                      + rmomu1(inode,igaus)   * p2vecdiff(jdime,knode,2,jnode,igaus) &
                                                      + rmomu1_der(jdime,knode,inode,igaus)   * p2vec(2,jnode,igaus)
            end do
          end do
          
        enddo
      enddo
    end do    
    
  else
  
    do igaus = 1,pgaus
      do jdime = 1,ndime
        do knode = 1,pnode    
    
          do inode = 1,pnode
            idof1 = 3*inode-2
            idof2 = idof1+1
            idof3 = idof2+1
            do jnode = 1,pnode
            
               elapudiff(jnode,idof1,jdime,knode) = elapudiff(jnode,idof1,jdime,knode) &     
                                                      + rcont(1,inode,igaus) * p2sca_der(jdime,knode,jnode,igaus) & 
                                                      + rcont_der(jdime,knode,1,inode,igaus) * p2sca(jnode,igaus) & 
                                                      + rmomu1(inode,igaus)   * p2vecdiff(jdime,knode,1,jnode,igaus) &
                                                      + rmomu1_der(jdime,knode,inode,igaus)   * p2vec(1,jnode,igaus)
                   
               elapudiff(jnode,idof2,jdime,knode) = elapudiff(jnode,idof2,jdime,knode) &
                                                      + rcont(2,inode,igaus) * p2sca_der(jdime,knode,jnode,igaus) & 
                                                      + rcont_der(jdime,knode,2,inode,igaus) * p2sca(jnode,igaus) & 
                                                      + rmomu1(inode,igaus)   * p2vecdiff(jdime,knode,2,jnode,igaus) &
                                                      + rmomu1_der(jdime,knode,inode,igaus)   * p2vec(2,jnode,igaus)
               
               elapudiff(jnode,idof3,jdime,knode) = elapudiff(jnode,idof3,jdime,knode) &
                                                      + rcont(3,inode,igaus) * p2sca_der(jdime,knode,jnode,igaus) & 
                                                      + rcont_der(jdime,knode,3,inode,igaus) * p2sca(jnode,igaus) & 
                                                      + rmomu1(inode,igaus)   * p2vecdiff(jdime,knode,3,jnode,igaus) &
                                                      + rmomu1_der(jdime,knode,inode,igaus)   * p2vec(3,jnode,igaus)
                                                      
            end do
          end do
          
        enddo
      enddo
    end do       
  
  end if
  
  !----------------------------------------------------------------------
  !
  ! elappdiff
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: ( grad(p) , tau1' grad(q) ) jdime,jnode
  ! 
  if( ndime == 2 ) then
  
    do igaus=1,pgaus
      do jdime = 1,ndime
        do knode = 1,pnode    
          do inode=1,pnode
            do jnode=inode+1,pnode
              fact1 =  p2vecdiff(jdime,knode,1,jnode,igaus) * gpcar(1,inode,igaus) &
                   & + p2vecdiff(jdime,knode,2,jnode,igaus) * gpcar(2,inode,igaus) &
                   & + p2vec(1,jnode,igaus) * gpcar_der(1,inode,jdime,knode,igaus) &
                   & + p2vec(2,jnode,igaus) * gpcar_der(2,inode,jdime,knode,igaus)

              elappdiff(jnode,inode,jdime,knode) = elappdiff(jnode,inode,jdime,knode) + fact1
              elappdiff(inode,jnode,jdime,knode) = elappdiff(inode,jnode,jdime,knode) + fact1
            end do
            fact1 =  p2vecdiff(jdime,knode,1,inode,igaus) * gpcar(1,inode,igaus) &
                 & + p2vecdiff(jdime,knode,2,inode,igaus) * gpcar(2,inode,igaus) &
                 & + p2vec(1,inode,igaus) * gpcar_der(1,inode,jdime,knode,igaus) &
                 & + p2vec(2,inode,igaus) * gpcar_der(2,inode,jdime,knode,igaus)
            elappdiff(inode,inode,jdime,knode) = elappdiff(inode,inode,jdime,knode) + fact1
          end do 
        end do
      end do
    end do
        
  else
  
    do igaus=1,pgaus
      do jdime = 1,ndime
        do knode = 1,pnode    
          do inode=1,pnode
            do jnode=inode+1,pnode
              fact1 =  p2vecdiff(jdime,knode,1,jnode,igaus) * gpcar(1,inode,igaus) &
                   & + p2vecdiff(jdime,knode,2,jnode,igaus) * gpcar(2,inode,igaus) &
                   & + p2vecdiff(jdime,knode,3,jnode,igaus) * gpcar(3,inode,igaus) &
                   & + p2vec(1,jnode,igaus) * gpcar_der(1,inode,jdime,knode,igaus) &
                   & + p2vec(2,jnode,igaus) * gpcar_der(2,inode,jdime,knode,igaus) &
                   & + p2vec(3,jnode,igaus) * gpcar_der(3,inode,jdime,knode,igaus)

              elappdiff(jnode,inode,jdime,knode) = elappdiff(jnode,inode,jdime,knode) + fact1
              elappdiff(inode,jnode,jdime,knode) = elappdiff(inode,jnode,jdime,knode) + fact1
            end do
            fact1 =  p2vecdiff(jdime,knode,1,inode,igaus) * gpcar(1,inode,igaus) &
                 & + p2vecdiff(jdime,knode,2,inode,igaus) * gpcar(2,inode,igaus) &
                 & + p2vecdiff(jdime,knode,3,inode,igaus) * gpcar(3,inode,igaus) &
                 & + p2vec(1,inode,igaus) * gpcar_der(1,inode,jdime,knode,igaus) &
                 & + p2vec(2,inode,igaus) * gpcar_der(2,inode,jdime,knode,igaus) &
                 & + p2vec(3,inode,igaus) * gpcar_der(3,inode,jdime,knode,igaus)
            elappdiff(inode,inode,jdime,knode) = elappdiff(inode,inode,jdime,knode) + fact1
          end do 
        end do
      end do
    end do
    
  end if
  
  !----------------------------------------------------------------------
  !
  ! elraudiff = elauudiff * u  + elaupdiff * p 
  ! elrapdiff = elapudiff * u  + elappdiff * p 
  !
  !----------------------------------------------------------------------
  
  !elauudiff * u
  do idime = 1, ndime
    do inode = 1, pnode
      call mbvab0(auuprodu(:,idime,inode),elauudiff(:,:,idime,inode),elvel1,pevat,pevat)
    enddo
  enddo
  
  !elaupdiff * p
  do idime = 1, ndime
    do inode = 1, pnode
      call mbvab0(aupprodp(:,idime,inode),elaupdiff(:,:,idime,inode),elpre1,pevat,pnode)
    enddo
  enddo
  
  !elapudiff * u
  do idime = 1, ndime
    do inode = 1, pnode
      call mbvab0(apuprodu(:,idime,inode),elapudiff(:,:,idime,inode),elvel1,pnode,pevat)
    enddo
  enddo
  
  !elappdiff * p
  do idime = 1, ndime
    do inode = 1, pnode
      call mbvab0(appprodp(:,idime,inode),elappdiff(:,:,idime,inode),elpre1,pnode,pnode)
    enddo
  enddo
  
  
  ! auuprodu + aupprodp --> elraudiff
  ! apuprodu + appprodp --> elrapdiff
  do jdime = 1, ndime
    do knode = 1, pnode
      do inode = 1,pnode
        elrapdiff(inode,jdime,knode) = apuprodu(inode,jdime,knode)+appprodp(inode,jdime,knode)
        do idime = 1,ndime
          ind = (inode-1)*ndime + idime
          elraudiff(idime,inode,jdime,knode) = auuprodu(ind,jdime,knode)+aupprodp(ind,jdime,knode)
        enddo
      enddo
    enddo !knode
  enddo !jdime

  
  ! elauudiff --> 0.0 (for assembly at the end)
  ! elaupdiff --> 0.0 (for assembly at the end)
  ! elapudiff --> 0.0 (for assembly at the end)
  ! elappdiff --> 0.0 (for assembly at the end)
  elauudiff = 0.0_rp
  elaupdiff = 0.0_rp
  elapudiff = 0.0_rp
  elappdiff = 0.0_rp

  !----------------------------------------------------------------------
  !
  ! elresucoodiff = elraudiff + elrbudiff (0)
  ! elrespcoodiff = elrapdiff + elrbpdiff (0)
  !
  !----------------------------------------------------------------------
  
  do jdime = 1, ndime
    do knode = 1, pnode
      do inode = 1,pnode
        elrespcoodiff(inode,jdime,knode) = elrapdiff(inode,jdime,knode) - elrbpdiff(inode,jdime,knode)
        do idime = 1,ndime
          elresucoodiff(idime,inode,jdime,knode) = elraudiff(idime,inode,jdime,knode) - elrbudiff(idime,inode,jdime,knode)
        enddo
      enddo
    enddo !knode
  enddo !jdime



end subroutine nsi_elmresdiff_coo

