!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_bouset(ibsec,ibset)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_bouset
  ! NAME
  !    chm_bouset
  ! DESCRIPTION
  !    This routine computes variables on a boundary set W.
  !    The variable are:
  !    1. setfl: set flux of class i  =  int_S k grad(Ci).n ds
  ! USES
  !    bouder
  !    chenor
  ! USED BY
  !    chm_outset
  !***
  !-----------------------------------------------------------------------
  use def_master,      only : postp, INOTMASTER, vbset, conce, advec
  use def_domain,      only : ndimb, mgaus, elmar, nboun, nmate, lbset, ltypb, nnode, ngaus,&
                              lnodb, coord, lelbo, ltype, lmate, lnods, ndime, mnodb, mnode,&
                              mgaub
  use def_chemic,      only : nclas_chm
  use def_kintyp,      only : ip, rp
  use mod_ker_proper,  only : ker_proper
  use mod_bouder,      only : bouder

  implicit none
  integer(ip), intent(in)  :: ibsec,ibset
  real(rp),    pointer     :: setsu(:)
  real(rp),    pointer     :: set_conce(:)
  real(rp),    pointer     :: set_mass_flux(:)
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb)
  real(rp)                 :: bocon(nclas_chm, mnodb)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: gbcon(nclas_chm, mgaub)
  real(rp)                 :: bovel(ndime,mnodb)
  real(rp)                 :: gbvel(ndime,mgaub)
  integer(ip)              :: ielem,inode,ipoin,nn
  integer(ip)              :: idime,igaub,iboun,inodb,pblty
  integer(ip)              :: ispec,dummi
  integer(ip)              :: pnodb,pmate,pnode,pelty,pgaus,pgaub
  real(rp)                 :: eucta,dsurf
  real(rp)                 :: gbden(mgaus),mcon

  external                 :: chenor

  if( INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! Initialization
     !
     !----------------------------------------------------------------------

     nn            =  postp(1) % nvabs + 1
     setsu         => vbset( nn:nn , ibset )             ! Surface
     set_mass_flux => vbset(  1: 1 , ibset )             ! Mass flux
     set_conce     => vbset(  2: 9 , ibset )             ! Mean species mass fraction

     setsu         = 0.0_rp
     set_conce     = 0.0_rp
     set_mass_flux = 0.0_rp
     !
     ! Loop over elements
     !
     boundaries: do iboun = 1,nboun

        if( lbset(iboun) == ibsec ) then

           !----------------------------------------------------------------
           !
           ! Element properties and dimensions and gather
           !
           !----------------------------------------------------------------

           pblty = ltypb(iboun)
           pnodb = nnode(pblty)
           pgaub = ngaus(pblty)
           pmate = 1

           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              do ispec = 1,nclas_chm
                 bocon(ispec,inodb) = conce(ipoin,ispec,1)
              enddo
              do idime = 1,ndime
                 bocod(idime,inodb) = coord(idime,ipoin)
                 bovel(idime,inodb) = advec(idime,ipoin,1) !!!DEFINE
              end do
           end do

           ielem = lelbo(iboun)
           pelty = ltype(ielem)
           if( nmate > 1 ) pmate = lmate(ielem)

           if (pelty > 0) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do

              gbcon = 0.0_rp
              gbvel = 0.0_rp
              do igaub = 1,pgaub
                 do inodb = 1,pnodb
                    do ispec = 1,nclas_chm
                       gbcon(ispec, igaub) = gbcon(ispec,igaub) + elmar(pblty)%shape(inodb,igaub) * bocon(ispec, inodb)
                    enddo
                    do idime = 1,ndime
                       gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                    end do
                 end do
              end do

              call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)

              !----------------------------------------------------------------
              !
              ! Loop over Gauss points
              !
              !----------------------------------------------------------------

              gauss_points: do igaub = 1,pgaub

                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
                      bocod,baloc,eucta)
                 call chenor(pnode,baloc,bocod,elcod)
                 dsurf = elmar(pblty)%weigp(igaub)*eucta
                 setsu = setsu + dsurf

                 !-------------------------------------------------------------
                 !
                 ! Mass flux
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(1) /= 0 .or.  postp(1) % npp_setsb(2) /= 0 ) then
                   do idime=1,ndime
                      set_mass_flux = set_mass_flux + gbden(igaub) * gbvel(idime,igaub) * baloc(idime,ndime) * dsurf
                   end do
                 end if

                 !-------------------------------------------------------------
                 !
                 ! Mass flux weghted average of concentrations
                 ! Here the numerator is summed, then in chm_oustet
                 ! it is divided by set_mass_flux
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(2) /= 0 ) then
                    do ispec = 1, min(8_ip,nclas_chm)
                       do idime=1,ndime
                          mcon = gbden(igaub)*gbcon(ispec,igaub)
                          set_conce(ispec) = set_conce(ispec) + mcon * gbvel(idime,igaub)*baloc(idime,ndime) * dsurf
                       end do
                    enddo
                 end if

              end do gauss_points
           end if
        end if

     end do boundaries

  end if

end subroutine chm_bouset
