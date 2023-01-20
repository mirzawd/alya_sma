!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_bouset.f90
!> @author  Mariano Vazquez
!> @date    April, 2010
!>          - Subroutine creation
!> @author  Gerard Guillamet
!>          August, 2018
!>          - Subroutine refactoring
!> @brief   This routine computes variables on a boundary set W.
!> @details
!>
!>          This routine computes variables on a boundary set W.
!>          The variable are:
!>          0        SETSU: surface                         =  int_W  = meas(W)
!>          1  -> 3  SETFO: Force                           =
!>          4        SETNF: Normal force                    =
!>          5  -> 7  SETDI: Averaged disp.                  =  Sum_{ui in W}/nodb (u)|_i
!>          8        SETDM: Averaged disp. (Magnitude)      =  sqrt(ux^2 + uy^2 + uz^2)
!>          9  -> 11 SETFR: Sum. Reaction force             =  Sum_{i in W} (Fext - Fint - M·a)|_i
!>          12       SETFM: Sum. reaction force (Magnitude) =  sqrt(fr_x^2 + fr_y^2 + fr_z^2)
!>          13 -> 15 SETCF: Sum. contact force              =  Sum(fc)
!>          16       SETFC: Sum. contact force (Magnitude)  =  sqrt(fc_x^2 + fc_y^2 + fc_z^2)
!>          17       SETMP: Pressure                        =
!>          18       SETNT: Normal traction                 =
!>          19 -> 21 SETCL: Sum. contact force in local     =
!>          22 -> 24 SETTR: Sum. force                      =  Sum(f)
!>          25       SETRT: Sum. force (Magnitude)          =  sqrt(f_x^2 + f_y^2 + f_z^2)
!>
!>          Values of SETDI and SETDM are averaged further on in sld_outset
!>
!>@warning  <GGU> SETFO and SETNF requires a revision
!>
!> @}
!------------------------------------------------------------------------------

subroutine sld_bouset(ibsec,ibset)

  use def_kintyp,         only : ip, rp
  use def_master,         only : ITER_K
  use def_master,         only : postp, displ, gisca, vbset, npoi1
  use def_domain,         only : ndime, npoin, nboun, ndimb
  use def_domain,         only : elmar, coord
  use def_domain,         only : mnodb, mnode, mgaus, mgaub
  use def_domain,         only : lnodb, lboel, ltype, lnods, lnnob, ltypb, lbset, lelbo
  use def_domain,         only : nnode, ngaus
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_parall,         only : commd
  use def_solidz,         only : gppio_sld, bvnat_sld
  use def_solidz,         only : kfl_fixno_sld, frxid_sld
  use def_solidz,         only : kfl_conta_sld, fcont_sld
  use mod_sld_csys,       only : sld_csys_rotuni
  use mod_bouder
  use mod_communications, only : PAR_SUM

  implicit none

  external                 :: chenor
  external                 :: memgen
  
  integer(ip), intent(in)  :: ibsec         !< Boundary number
  integer(ip), intent(in)  :: ibset         !< Set code
  real(rp),    pointer     :: setsu(:),setfo(:),setnf(:),setfr(:),setfm(:),setdi(:),setdm(:)
  real(rp),    pointer     :: setfc(:),setcf(:),setcl(:)
  real(rp),    pointer     :: setmp(:),setnt(:),settr(:),setrt(:)
  integer(ip)              :: ielem,inode,ipoin,igaus,idime,jdime,nn,jj,idofn
  integer(ip)              :: pnode,pgaus,iboun,igaub,inodb
  integer(ip)              :: pelty,pblty,pnodb,pgaub
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                 :: bodis(ndime,mnodb),eldis(ndime,mnode)
  real(rp)                 :: gppio(ndime,ndime,mgaus)              ! Piola tensor P
  real(rp)                 :: gbsur,eucta
  real(rp)                 :: gbdis(ndime,mgaub)
  real(rp)                 :: gbpre(mgaub)
  real(rp)                 :: gbtra(ndime,mgaub)
  real(rp)                 :: dummr(3),tract(3),tmatr,fcont(ndime)
  real(rp)                 :: multiplicity_loc

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  nn    =  postp(1) % nvabs + 1
  setsu => vbset( nn:nn, ibset ) ! Surface
  setfo => vbset(  1:3 , ibset ) ! Force
  setnf => vbset(  4:4 , ibset ) ! Normal Force
  setdi => vbset(  5:7 , ibset ) ! Averaged displacement
  setdm => vbset(  8:8 , ibset ) ! Averaged displacement (Magnitude)
  setfr => vbset(  9:11, ibset ) ! Sum of reaction forces
  setfm => vbset( 12:12, ibset ) ! Sum of reaction forces (Magnitude)
  setcf => vbset( 13:15, ibset ) ! Sum of contact forces
  setfc => vbset( 16:16, ibset ) ! Sum of contact forces (Magnitude)
  setmp => vbset( 17:17, ibset ) ! Pressure
  setnt => vbset( 18:18, ibset ) ! Normal traction
  setcl => vbset( 19:21, ibset ) ! Sum of contact forces (local system)
  settr => vbset( 22:24, ibset ) ! Sum of traction forces (local system)
  setrt => vbset( 25:25, ibset ) ! Sum of traction forces (Magnitude)
  setsu =  0.0_rp
  setfo =  0.0_rp
  setnf =  0.0_rp
  setdi =  0.0_rp
  setdm =  0.0_rp
  setfr =  0.0_rp
  setfm =  0.0_rp
  setcf =  0.0_rp
  setfc  = 0.0_rp
  setmp =  0.0_rp
  setnt =  0.0_rp
  setcl =  0.0_rp
  settr =  0.0_rp
  setrt =  0.0_rp
  
  boundaries: do iboun = 1,nboun

     if ( lbset(iboun) == ibsec ) then

        !----------------------------------------------------------------
        !
        ! Element properties, dimensions and gather
        !
        !----------------------------------------------------------------

        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        pgaub = ngaus(pblty)

        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
           bodis(1:ndime,inodb) = displ(1:ndime,ipoin,ITER_K)
        end do

        ielem = lelbo(iboun)
        pelty = ltype(ielem)
        if ( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              eldis(1:ndime,inode) = displ(1:ndime,ipoin,ITER_K)
           end do

           !----------------------------------------------------------------
           !
           ! Initialize values at Gauss points
           !
           !----------------------------------------------------------------

           gbdis(:,:)   = 0.0_rp
           gbpre(:)     = 0.0_rp
           gbtra(:,:)   = 0.0_rp
           gppio(:,:,:) = 0.0_rp

           !----------------------------------------------------------------
           !
           ! Loop over Gauss points
           !
           !----------------------------------------------------------------

           gauss_points: do igaub = 1,pgaub

              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&  ! Cartesian derivative
                   bocod,baloc,eucta)                                 ! and Jacobian
              gbsur = elmar(pblty)%weigp(igaub)*eucta
              setsu = setsu + gbsur
              call chenor(pnode,baloc,bocod,elcod)                    ! Check normal

              !-------------------------------------------------------------
              !
              ! Compute traction
              ! t = P . N  (1st piola . normal in the reference system)
              !
              !-------------------------------------------------------------

              if ( postp(1) % npp_setsb(1) /= 0 .or. postp(1) % npp_setsb(4) /= 0 ) then
                 !
                 ! Extrapolate from Gauss points igauss to node inode and
                 ! then interpolation at boundary Gauss points igaub
                 !
                 do igaus = 1,pgaus
                    do inodb = 1,pnodb
                       tmatr =  elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                            & * elmar(pblty)%shape(inodb,igaub)
                       do jdime = 1,ndime
                          do idime = 1,ndime
                             gppio(idime,jdime,igaub) = gppio(idime,jdime,igaub) &
                                  + gppio_sld(ielem) % a(idime,jdime,igaus)*tmatr
                          end do
                       end do
                    end do
                 end do

                 do idime = 1,ndime
                    tract(idime)= 0.0_rp
                    do jdime = 1,ndime
                       tract(idime) = tract(idime) &
                            + gppio(jdime,idime,igaub) * baloc(jdime,ndime)
                    end do
                 end do
                 !
                 ! FORCE, F_y, F_z
                 !
                 if ( postp(1) % npp_setsb(1) /= 0  ) then
                    ! OPTION 1 : (x y z) component of the force
                    setfo(1) = setfo(1) + gbsur * tract(1)
                    setfo(2) = setfo(2) + gbsur * tract(2)
                    setfo(3) = setfo(3) + gbsur * tract(3)
                 end if
                 !
                 ! NORMA (Normal force)
                 !
                 if ( postp(1) % npp_setsb(4) /= 0  ) then
                    dummr(1) = 0.0_rp
                    do idime = 1,ndime
                       dummr(1) = dummr(1) + tract(idime) * baloc(idime,ndime)
                    end do
                    setnf(1) = setnf(1) + dummr(1) * gbsur
                 end if

              end if

              !-------------------------------------------------------------
              !
              ! Averaged displacement
              !
              !-------------------------------------------------------------

              if( any(postp(1) % npp_setsb(5:8) /= 0_ip) ) then

                 do inodb = 1,pnodb
                    gbdis(1:ndime,igaub) = gbdis(1:ndime,igaub) + &
                         elmar(pblty)%shape(inodb,igaub) * bodis(1:ndime,inodb)
                 end do
                 !
                 ! DIBOX, DIBOY, DIBOZ
                 !
                 setdi(1:ndime) = setdi(1:ndime) + gbsur * gbdis(1:ndime,igaub)
                 !
                 ! DIBOU (Magnitude)
                 !
                 if ( postp(1) % npp_setsb(8) /= 0 ) setdm = sqrt(sum(setdi(1:ndime)**2))

              end if

              !-------------------------------------------------------------
              !
              ! Pressure and Normal traction
              !
              !-------------------------------------------------------------

              if( any(postp(1) % npp_setsb(17:18) /= 0_ip) .or. &
                  any(postp(1) % npp_setsb(22:25) /= 0_ip) ) then

                 do igaus = 1,pgaus
                    do inodb = 1,pnodb
                       tmatr =  elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                            & * elmar(pblty)%shape(inodb,igaub)
                       do jdime = 1,ndime
                          do idime = 1,ndime
                             gppio(idime,jdime,igaub) = gppio(idime,jdime,igaub) &
                                  + gppio_sld(ielem) % a(idime,jdime,igaus)*tmatr
                          end do
                       end do
                    end do
                 end do

                 ! Compute pressure
                 do idime = 1,ndime
                    tract(idime)= 0.0_rp
                    do jdime = 1,ndime
                       tract(idime) = tract(idime) &
                            + gppio(jdime,idime,igaub) * baloc(jdime,ndime)
                    end do
                 end do

                 do inodb = 1,pnodb
                    gbpre(igaub) = gbpre(igaub) &
                         + elmar(pblty)%shape(inodb,igaub) * bvnat_sld(1,iboun,ITER_K)
                 end do
                 !
                 ! PRESS
                 !
                 if( postp(1) % npp_setsb(17) /= 0_ip ) then
                    setmp = setmp + gbsur * gbpre(igaub)
                 end if
                 !
                 ! TRNOR: Normal traction
                 !
                 if( postp(1) % npp_setsb(18) /= 0_ip ) then
                    dummr(1) = 0.0_rp
                    do idime = 1,ndime
                       dummr(1) = dummr(1) + tract(idime) * baloc(idime,ndime)
                    end do
                    setnt = setnt + dummr(1) * gbsur
                 end if
                 !
                 ! FOBOU: Forces
                 !
                 if( any(postp(1) % npp_setsb(22:25) /= 0_ip) ) then
                    do inodb = 1,pnodb
                       gbtra(1:ndime,igaub) = gbtra(1:ndime,igaub) + &
                            elmar(pblty)%shape(inodb,igaub) * bvnat_sld(1:ndime,iboun,ITER_K)
                    end do
                    !
                    ! FOBOX, FOBOY, FOBOZ
                    !
                    settr(1:ndime) = settr(1:ndime) + gbsur * gbtra(1:ndime,igaub)
                    !
                    ! FOBOU (Magnitude)
                    !
                    if( postp(1) % npp_setsb(25) /= 0 ) setrt = sqrt(sum(settr(1:ndime)**2))

                 end if
                 
              end if

           end do gauss_points

        end if

     end if

  end do boundaries

  !----------------------------------------------------------
  !
  ! Sum of the Reaction and Contact forces
  !
  !----------------------------------------------------------

  if( any(postp(1) % npp_setsb(9:16) /= 0_ip) .or. &
      any(postp(1) % npp_setsb(19:21) /= 0_ip) ) then

     call memgen(1_ip,npoin,0_ip)
     do iboun = 1,nboun
        if ( lbset(iboun) == ibset ) then
           do inodb = 1,lnnob(iboun)
              ipoin = lnodb(inodb,iboun)
              gisca(ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY CODE')

     !
     ! FRBOX, FRBOY, FRBOZ (and FRBOU)
     !
     if( any(postp(1) % npp_setsb(9:12) /= 0_ip) ) then
        do ipoin = 1,npoin
           idofn = (ipoin-1)*ndime
           if ( gisca(ipoin) /= 0_ip ) then
              do idime = 1,ndime
                 if ( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                    setfr(idime) = setfr(idime) - frxid_sld(idofn+idime) / real(gisca(ipoin),rp)
                 end if
              end do
           end if
        end do
        ! FRBOU
        if ( postp(1) % npp_setsb(12) /= 0_ip ) setfm = sqrt(sum(setfr(1:ndime)**2))
        
     end if
     !
     ! FCONX, FCONY, FCONZ (and FCONT)
     !
     if( any(postp(1) % npp_setsb(13:16) /= 0_ip) ) then

        if( kfl_conta_sld /= 0_ip ) then
           jj = 0
           do ipoin = 1,npoin
              idofn = (ipoin-1)*ndime
              if( ipoin <= npoi1 ) then
                 multiplicity_loc = 1.0_rp
              else
                 jj = jj + 1
                 multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
              end if
              if( kfl_fixno_sld(1,ipoin) == 3_ip ) then
                 fcont = fcont_sld(1:ndime,ipoin)
                 setcf(1:ndime) = setcf(1:ndime) + fcont(1:ndime) / multiplicity_loc
              end if
           end do
           ! FCONT
           if ( postp(1) % npp_setsb(16) /= 0_ip ) setfc = sqrt(sum(setcf(1:ndime)**2))
        end if

     end if
     !
     ! FCONO, FCOT1, FCOT2
     !
     if( any(postp(1) % npp_setsb(19:21) /= 0_ip) ) then

        if( kfl_conta_sld /= 0_ip ) then
           jj = 0
           do ipoin = 1,npoin
              idofn = (ipoin-1)*ndime
              if( ipoin <= npoi1 ) then
                 multiplicity_loc = 1.0_rp
              else
                 jj = jj + 1
                 multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
              end if
              if( kfl_fixno_sld(1,ipoin) == 3_ip ) then
                 fcont = fcont_sld(1:ndime,ipoin)
                 ! Global --> Local
                 call sld_csys_rotuni(1_ip,ndime,ipoin,fcont)
                 setcl(1:ndime) = setcl(1:ndime) + fcont(1:ndime) / multiplicity_loc
              end if
           end do
        end if

     end if

     call memgen(3_ip,npoin,0_ip)

  end if

end subroutine sld_bouset
