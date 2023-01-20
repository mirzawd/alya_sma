!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine cartbo(&
     itask,lboel,gbsha,shaga,gpcar,gpsha,shapp,&
     gbcar,pnodb,pnode,pgaus)
  !-----------------------------------------------------------------------
  !****f* Domain/cartbo
  ! NAME 
  !    cartbo
  ! DESCRIPTION
  !    This routine computes the following extrapolating from Gauss
  !    points to boundaries:
  !    ITASK = 1,2 ... Cartesian derivatives on the boundary Gauss points 
  !          = 2 ..... Shape function ate boundary Gauss points
  ! USES
  ! USED BY
  !    nsi_bouset
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  implicit none
  integer(ip), intent(in)  :: itask,pnodb,pnode,pgaus
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbsha(pnodb),shaga(pgaus,pnode)
  real(rp),    intent(out) :: gbcar(ndime,pnode),shapp(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  integer(ip)              :: inode,idime,inodb,igaus

  if( itask == 3 ) then
     !
     ! Shape function gpsha at boundary Gauss points
     !     
     do inode = 1,pnode                      
        shapp(inode) = 0.0_rp                 
        do inodb = 1,pnodb                    
           do igaus = 1,pgaus
              shapp(inode) = shapp(inode) &
                   + shaga(igaus,lboel(inodb)) * gbsha(inodb) &
                   * gpsha(inode,igaus)
           end do
        end do
     end do
  else
     !
     ! Cartesian derivates at boundary Gauss points 
     !            
     do inode = 1,pnode                                       
        do idime = 1,ndime                                    
           gbcar(idime,inode) = 0.0_rp                        
           do inodb = 1,pnodb                                  
              do igaus = 1,pgaus
                 gbcar(idime,inode) = gbcar(idime,inode)&
                      + shaga(igaus,lboel(inodb)) * gbsha(inodb)&
                      * gpcar(idime,inode,igaus)
              end do
           end do
        end do
     end do

     if( itask == 2 ) then
        !
        ! Shape function gpsha at boundary Gauss points
        !     
        do inode=1,pnode                      
           shapp(inode)=0.0_rp                 
           do inodb=1,pnodb                    
              do igaus=1,pgaus
                 shapp(inode)=shapp(inode)&
                      + shaga(igaus,lboel(inodb)) * gbsha(inodb)&
                      * gpsha(inode,igaus)
              end do
           end do
        end do
     end if
  end if

end subroutine cartbo
