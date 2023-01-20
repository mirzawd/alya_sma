!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmsld(&
     pnode,pgaus,pgaub,pnodb,lboel,lnodf,gbcar,gpcar,gbsha,&
     shaga,baloc,gbvis,gbsur,elvel,bopre,gpgve,gbpre,elrhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmsld
  ! NAME 
  !    nsi_elmsld
  ! DESCRIPTION
  !   Assemble force
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  forcf
  use def_domain, only       :  ndime,mnode
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pgaub,pnodb
  integer(ip), intent(in)    :: lboel(pnodb)
  integer(ip), intent(in)    :: lnodf(pnodb)
  real(rp),    intent(out)   :: gbcar(ndime,pnode)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gbsha(pnodb,pgaub)
  real(rp),    intent(in)    :: shaga(pgaus,pnode)
  real(rp),    intent(inout) :: baloc(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gbvis(pgaub)
  real(rp),    intent(in)    :: gbsur(pgaub)
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: bopre(pnodb)
  real(rp),    intent(out)   :: gpgve(ndime,ndime) 
  real(rp),    intent(out)   :: gbpre
  real(rp),    intent(out)   :: elrhs(ndime,pnodb)
  integer(ip)                :: idime,jdime,igaus,inode,inodb,igaub,ipoin
  real(rp)                   :: fact1,fact2

  do igaub = 1,pgaub
     !
     ! Normalize baloc
     !
     fact1 = 0.0_rp
     do idime = 1,ndime
        fact1 = fact1 + baloc(idime,ndime,igaub) * baloc(idime,ndime,igaub)
     end do
     fact1 = 1.0_rp / sqrt(fact1)
     do idime = 1,ndime
        baloc(idime,ndime,igaub) = fact1 * baloc(idime,ndime,igaub)
     end do 
     !
     ! dNi/dx
     !
     do inode = 1,pnode                                       
        do idime = 1,ndime                                    
           gbcar(idime,inode) = 0.0_rp                        
           do inodb = 1,pnodb                                  
              do igaus = 1,pgaus
                 gbcar(idime,inode) = gbcar(idime,inode)&
                      + shaga(igaus,lboel(inodb)) * gbsha(inodb,igaub)&
                      * gpcar(idime,inode,igaus)
              end do
           end do
        end do
     end do
     !
     ! Pressure and Grad(u)
     !  
     gbpre = 0.0_rp
     do inodb = 1,pnodb
        gbpre = gbpre + bopre(inodb) * gbsha(inodb,igaub)
     end do
     do idime = 1,ndime
        do jdime = 1,ndime
           gpgve(jdime,idime) = 0.0_rp
        end do
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           do jdime = 1,ndime
              gpgve(jdime,idime) = gpgve(jdime,idime) &
                   + elvel(idime,inode) * gbcar(jdime,inode)
           end do
        end do
     end do
     !
     ! Assemble - sig.n = - ( -p n + 2*mu*1/2*[grad(u)+grad(u)^t].n )
     !
     fact1 = gbsur(igaub) * gbpre
     fact2 = gbsur(igaub) * gbvis(igaub)
     do inodb = 1,pnodb
        do idime = 1,ndime
           elrhs(idime,inodb) = elrhs(idime,inodb) &
                + fact1 * baloc(idime,ndime,igaub) * gbsha(inodb,igaub)
           do jdime =1 ,ndime
              elrhs(idime,inodb) = elrhs(idime,inodb) &
                   - fact2 * ( gpgve(idime,jdime) +  gpgve(jdime,idime) ) &
                   * baloc(jdime,ndime,igaub)  * gbsha(inodb,igaub)   
           end do
        end do
     end do
  end do
  !
  ! Assemble
  !
  if( ndime == 2 ) then
     do inode = 1,pnodb
        ipoin = lnodf(inode)
        forcf(1,ipoin) = forcf(1,ipoin) + elrhs(1,inode)
        forcf(2,ipoin) = forcf(2,ipoin) + elrhs(2,inode)
     end do
  else
     do inode = 1,pnodb
        ipoin = lnodf(inode)
        forcf(1,ipoin) = forcf(1,ipoin) + elrhs(1,inode)
        forcf(2,ipoin) = forcf(2,ipoin) + elrhs(2,inode)
        forcf(3,ipoin) = forcf(3,ipoin) + elrhs(3,inode)
     end do
  end if

end subroutine nsi_elmsld

