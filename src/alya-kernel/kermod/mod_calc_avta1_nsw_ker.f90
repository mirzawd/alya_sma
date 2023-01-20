!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_calc_avta1_nsw_ker.f90
!> @author  Herbert Owen
!> @brief   Obtains some values needed to apply wall law increasing viscosity in the first element and using no slip 
!> @details - obtains: avta1_nsw_ker 
!>          - avta1_nsw_ker    : xxxxd.
!> @} 
!-----------------------------------------------------------------------
module mod_calc_avta1_nsw_ker
  use def_kintyp,                   only : ip,rp
  implicit none
  private
  public :: calc_avta1_nsw_ker

contains
  subroutine calc_avta1_nsw_ker(avwei)
    use def_master,                   only : INOTMASTER
    use def_kermod,                   only : kfl_noslw_ker, kfl_nswel_ker, avta1_nsw_ker, normal_nsw_ker, avupo_ker!  , avwei_ker
    use def_domain,                   only : lnods,nelem,nnode,ntens,ltype,mgaus,coord,elmar,ndime,mnode,ngaus
    use mod_ker_proper,               only : ker_proper
    implicit none
    real(rp),intent(in)               :: avwei

    integer(ip)   :: ielem,inode,ipoin,idime,ikoun,igaus,jdime
    integer(ip)   :: pnode,plapl,pelty,pgaus
    integer(ip)   :: dummi
    real(rp)      :: auxi,auxde,auxvi,auxmut
    real(rp)      :: elcod(ndime,mnode)

    real(rp)      :: gpgnavv(ndime,mgaus)  ! Gauss Point Gradient in Normal dir of AVerage Velocity  - employ  mgaus since pgaus is not defined
    real(rp)      :: elavv(ndime,mnode)    ! employ  mnode since pnode is not defined
    real(rp)      :: elnnsw(ndime)
    real(rp)      :: elgnavv(ndime)        ! ELement Gradient in Normal dir of AVerage Velocity 

    real(rp)      :: gpmut(mgaus)                          ! mut
    real(rp)      :: gpvis(mgaus)                          ! mu
    real(rp)      :: gpden(mgaus)                          ! Density    
    real(rp)      :: gpsha(mnode,mgaus)                    ! N
    real(rp)      :: gpder(ndime,mnode,mgaus)              ! dN/dsi                
    real(rp)      :: gpcar(ndime,mnode,mgaus)              ! dN/dxi
    real(rp)      :: gphes(ntens,mnode,mgaus)              ! d2N/dxidxj
    real(rp)      :: gpvol(mgaus)


    if ( kfl_noslw_ker /= 0_ip ) then 

       if( INOTMASTER ) then     ! pensar que pasa con el master
          !
          ! Loop over elements - this could be vectorized - but perhaps it does not introduce much benefit since it is only needed in some elements
          !
          elem0: do ielem = 1,nelem

             ikoun = kfl_nswel_ker(ielem)
             if( ikoun > 0_ip ) then
                !
                ! Element properties and dimensions
                !
                pelty = ltype(ielem) 
                pnode = nnode(pelty)
                pgaus = ngaus(pelty)
                plapl = 0_ip 
                !        porde = lorde(pelty)
                !        ptopo = ltopo(pelty)
                !        pevat = ndime * pnode
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   do idime = 1,ndime
                      elcod(idime,inode) = coord(idime,ipoin)
                      elavv(idime,inode) = avupo_ker(idime,ipoin)
                   end do
                end do
                do idime = 1,ndime
                   elnnsw(idime) = normal_nsw_ker(idime,kfl_nswel_ker(ielem)) 
                end do
                !
                ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL  ! borrowed from nsi_elmope - just for gpsha & gpcar
                !
                call elmca2(&
                     pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                     elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                     gpder,gpcar,gphes,ielem)
                !--------------------------------------------------------------------
                !
                ! Properties and local time step
                !
                !--------------------------------------------------------------------

                call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)     ! rho
                call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,gpsha,gpcar)     ! mu
                call ker_proper('TURBU','PGAUS',dummi,ielem,gpmut,pnode,pgaus,gpsha,gpcar)     ! mut

                !----------------------------------------------------------------------
                !
                ! Gauss point values
                !
                !----------------------------------------------------------------------
                !
                ! GPGNAVV = dj uav_i nj =    dj N_I_i nj  Uav_I
                !
                gpgnavv(:,:)   = 0.0_rp

                do igaus = 1,pgaus
                   do inode = 1,pnode
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            gpgnavv(idime,igaus) = gpgnavv(idime,igaus) + elavv(idime,inode) * &
                                 gpcar(jdime,inode,igaus) * elnnsw(jdime)
                         end do
                      end do
                   end do
                end do
                !
                ! Obtain the average over gauss points of the normal gradient of the average velocity for all gauss points   - I could merge this loop with the upper one and obtain elgnavv directly
                ! Also average density and viscosity
                !
                elgnavv  = 0.0_rp
                auxvi = 0.0_rp
                auxde = 0.0_rp
                auxmut = 0.0_rp
                do igaus = 1,pgaus
                   do idime = 1,ndime
                      elgnavv(idime)  = elgnavv(idime)  + gpgnavv(idime,igaus)
                   end do
                   auxvi  = auxvi  + gpvis(igaus)
                   auxde  = auxde  + gpden(igaus)
                   auxmut = auxmut + gpmut(igaus)
                end do
                do idime = 1,ndime
                   elgnavv(idime) = elgnavv(idime) / real(pgaus,rp)
                end do
                auxvi  = auxvi  / real(pgaus,rp)
                auxde  = auxde  / real(pgaus,rp)
                auxmut = auxmut / real(pgaus,rp)
                !
                ! Substract normal component to keep only tangential one
                !
                !
                auxi = 0.0_rp
                do idime=1,ndime
                   auxi = auxi + elgnavv(idime) * elnnsw(idime)
                end do
                do idime=1,ndime
                   elgnavv(idime) = elgnavv(idime) - auxi * elnnsw(idime)
                end do
                !
                ! Time average of  (mu+mut) d u_t / dni --- beware at the begining of the run avta1_nsw  might have erroneous values -- we should find a better way to initialise it 
                !
                do idime=1,ndime
                   avta1_nsw_ker(idime,ikoun) =  (1.0_rp-avwei) * avta1_nsw_ker(idime,ikoun) +  &
                        (  avwei * (auxvi + auxde * auxmut ) ) * elgnavv(idime)

!                   cambiando a esta linea logro que me3 se comporte como me5
!                   avta1_nsw_ker(idime,ikoun) =   (auxvi + auxde * auxmut ) * elgnavv(idime)
                end do

             end if

          end do elem0

       end if

    end if

  end subroutine calc_avta1_nsw_ker

end module mod_calc_avta1_nsw_ker
