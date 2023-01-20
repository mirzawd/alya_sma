!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_inibcs()
  !-----------------------------------------------------------------------
  !****f* partis/chm_inibcs
  ! NAME
  !    chm_inibcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! OUTPUT
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master,       only : INOTMASTER, fleve, conce
  use def_domain,       only : bvess, ifbop, iffun, ifloc, kfl_fixbo, kfl_fixno, kfl_funno, kfl_funtn, kfl_icodb, kfl_icodn, nboun,&
                               nelem, neset, npoin, tbcod, tncod, coord, leset, lnods, lnnod
  use def_kintyp,       only : ip, rp
  use def_chemic,       only : nvar_CMC_chm, kfl_allcl_chm, kfl_bc_alpha_CMC_chm, kfl_conbc_chm, kfl_model_chm, kfl_robin_chm,&
                               kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm, nclas_chm, nvar_therm_CMC_chm, kfl_fixno_chm,&
                               bvess_chm, kfl_field_chm, kfl_initi_chm, xinit_chm, kfl_fixbo_chm, tncod_chm, kfl_funno_chm,&
                               kfl_funtn_chm, tbcod_chm

  implicit none
  integer(ip)  :: ipoin,iclas,iboun,inode,ielem,ieset
  real(rp)     :: value_xxx
  real(rp)     :: Heaviside,x,y,rc,xc,yc
  integer(ip)  :: pos_fields(nvar_CMC_chm)

  external     :: chm_membcs
  external     :: reacod
  external     :: runend

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call chm_membcs(1_ip)

     if( kfl_conbc_chm == 0 ) then
        call chm_membcs(4_ip)
     end if

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        ifloc     =  0
        ifbop     =  0

        if( kfl_allcl_chm == 1 ) then

           iclas     =  1
           if(kfl_conbc_chm==0) then
              iffun     =  1
           else
              iffun      =  0
           end if
           kfl_fixno => kfl_fixno_chm(iclas:iclas,:)
           bvess     => bvess_chm(iclas:iclas,:)
           tncod     => tncod_chm(iclas:)
           call reacod(10_ip)
           do ipoin = 1,npoin
              do iclas = 2,nclas_chm
                 kfl_fixno_chm(iclas,ipoin) = kfl_fixno_chm(1,ipoin)
                 bvess_chm(iclas,ipoin)     = bvess_chm(1,ipoin)
              end do
           end do

        else

           if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then

              if (kfl_bc_alpha_CMC_chm == 0_ip) then

                 if (kfl_solve_enth_CMC_chm == 0_ip) then
                    if (kfl_field_chm(2) < kfl_field_chm(1) .and. kfl_field_chm(2) > 0) then
                       ! Temperature listed before species in the boundary conditions file but not used
                       pos_fields = (/( iclas, iclas = 2_ip, nclas_chm+1_ip )/)
                    else
                       pos_fields = (/( iclas, iclas = 1_ip, nclas_chm )/)
                    end if

                 else
                    if (kfl_field_chm(2) < kfl_field_chm(1)) then
                       ! Temperature listed before species in the boundary conditions file
                       pos_fields = (/( iclas, iclas = 2_ip, nclas_chm+1_ip ), 1_ip/)
                    else
                       ! Temperature listed after species in the boundary conditions file
                       pos_fields = (/( iclas, iclas = 1_ip, nclas_chm+1_ip )/)
                    end if

                 end if
              else

                 if (kfl_solve_enth_CMC_chm == 0_ip) then
                    if (kfl_field_chm(2) < kfl_field_chm(1) .and. kfl_field_chm(2) > 0) then
                       ! Temperature listed before species in the boundary conditions file but not used
                       pos_fields = (/( 2_ip, iclas = 1_ip, nclas_chm )/)
                    else
                       pos_fields = (/( 1_ip, iclas = 1_ip, nclas_chm )/)
                    end if

                 else
                    if (kfl_field_chm(2) < kfl_field_chm(1)) then
                       ! Temperature listed before species in the boundary conditions file
                       pos_fields = (/( 2_ip, iclas = 1_ip, nclas_chm ), 1_ip/)
                    else
                       ! Temperature listed after species in the boundary conditions file
                       pos_fields = (/( 1_ip, iclas = 1_ip, nclas_chm ), 2_ip/)
                    end if

                 end if
              end if

              do iclas = 1, nvar_therm_CMC_chm   ! Since we do not consider explicitly Dirichlet BCs for soot
                 kfl_fixno => kfl_fixno_chm(iclas:iclas,:)
                 bvess     => bvess_chm(iclas:iclas,:)
                 tncod     => tncod_chm(pos_fields(iclas):)
                 kfl_funno => kfl_funno_chm(:,iclas)
                 kfl_funtn => kfl_funtn_chm(:,iclas)
                 if(kfl_conbc_chm==0) then
                    iffun     =  1
                 else
                    iffun     =  0
                 end if
                 call reacod(10_ip)

              end do

           else
              do iclas = 1,nclas_chm
                 kfl_fixno => kfl_fixno_chm(iclas:iclas,:)
                 bvess     => bvess_chm(iclas:iclas,:)
                 tncod     => tncod_chm(iclas:)
                 kfl_funno => kfl_funno_chm(:,iclas)
                 kfl_funtn => kfl_funtn_chm(:,iclas)
                 if(kfl_conbc_chm==0) then
                    iffun     =  1
                 else
                    iffun     =  0
                 end if
                 call reacod(10_ip)

              end do
           end if
        end if

     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        do iclas = 1,nclas_chm
           kfl_fixbo => kfl_fixbo_chm(iclas,:)
           tbcod     => tbcod_chm(iclas:)
           call reacod(20_ip)
        end do

     end if

     !-------------------------------------------------------------
     !
     ! Initial solution
     !
     !-------------------------------------------------------------

     do iclas = 1,nclas_chm
        !
        ! Initialize to zero, for any case:
        !
        do ipoin=1,npoin
           if(kfl_fixno_chm(iclas,ipoin)==0) then
              bvess_chm(iclas,ipoin)=0.0_rp
           end if
        end do


        if( kfl_initi_chm(iclas) == 1 ) then

           call runend('CHM_INIBCS: kfl_initi_chm(iclas) == 1 DEPRECATED, CHECKOUT REVISION BEFORE DEC 2018')

        else if( kfl_initi_chm(iclas) == 2 ) then
           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 bvess_chm(iclas,ipoin)=xinit_chm(iclas,1)
              end if
           end do
        else if( kfl_initi_chm(iclas) == 3 ) then ! Levels
           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 if ( fleve(ipoin,1)>=0) then
                    Heaviside=1.0_rp
                 else
                    Heaviside=0.0_rp
                 endif
                 bvess_chm(iclas,ipoin)=xinit_chm(iclas,1) + (xinit_chm(iclas,2)-xinit_chm(iclas,1))*Heaviside
              endif
           end do

        else if( kfl_initi_chm(iclas) == 4 ) then
           !
           !  circle of radius rc arround xc,yc
           !
           xc=0.005_rp
           yc=0.005_rp
           rc=0.002_rp

           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 x=coord(1,ipoin)
                 y=coord(2,ipoin)
                 if ( ((x-xc)**2+(y-yc)**2) - rc**2 >= 0.0_rp) then
                    Heaviside=1.0_rp
                 else
                    Heaviside=0.0_rp
                 endif
                 bvess_chm(iclas,ipoin)=xinit_chm(iclas,1) + (xinit_chm(iclas,2)-xinit_chm(iclas,1))*Heaviside
              endif
           end do


        else if( kfl_initi_chm(iclas) < 0 ) then
           ieset = -kfl_initi_chm(iclas)
           if( neset == 0 ) call runend('CHM_INIBCS: ELEMENT SETS NOT DEFINED')
           do ielem = 1,nelem
              if( leset(ielem) == ieset ) then
                 do inode = 1,lnnod(ielem)
                    ipoin = lnods(inode,ielem)
                    if(kfl_fixno_chm(iclas,ipoin)==0) then
                       bvess_chm(iclas,ipoin)=xinit_chm(iclas,1)
                    end if
                 end do
              end if
           end do
        end if
     end do

     do iclas = 1,nclas_chm
        if( kfl_initi_chm(iclas) == 2 ) then
           do ipoin=1,npoin
              conce(ipoin,iclas,1)=xinit_chm(iclas,1)
              conce(ipoin,iclas,2)=xinit_chm(iclas,1)
              conce(ipoin,iclas,3)=xinit_chm(iclas,1)
           end do
        else if (kfl_initi_chm(iclas) == 3) then

           do ipoin=1,npoin
              if ( fleve(ipoin,1)>=0) then
                 Heaviside=1.0_rp
              else
                 Heaviside=0.0_rp
              endif
              value_xxx= xinit_chm(iclas,1)  + (xinit_chm(iclas,2)-xinit_chm(iclas,1)) * Heaviside
              conce(ipoin,iclas,1)=value_xxx
              conce(ipoin,iclas,2)=value_xxx
              conce(ipoin,iclas,3)=value_xxx
           end do

        else if (kfl_initi_chm(iclas) == 4) then

           xc= 0.005_rp
           yc= 0.005_rp
           rc= 0.002_rp
           do ipoin=1,npoin
              x=coord(1,ipoin)
              y=coord(2,ipoin)
              if ( ((x-xc)**2+(y-yc)**2) - rc**2 >= 0.0_rp) then
                 Heaviside=1.0_rp
              else
                 Heaviside=0.0_rp
              endif
              value_xxx= xinit_chm(iclas,1)  + (xinit_chm(iclas,2)-xinit_chm(iclas,1)) * Heaviside
              conce(ipoin,iclas,1)=value_xxx
              conce(ipoin,iclas,2)=value_xxx
              conce(ipoin,iclas,3)=value_xxx
           end do

        end if
     end do

     !-------------------------------------------------------------
     !
     ! KFL_ROBIN_CHM: Count number of Robin boundaries
     !
     !-------------------------------------------------------------

     kfl_robin_chm = 0
     iboun         = 0
     do while( iboun < nboun )
        iboun = iboun+1
        iclas = 0
        do while( iclas < nclas_chm )
           iclas = iclas+1
           if( kfl_fixbo_chm(iclas,iboun) == 2 ) then
              kfl_robin_chm = kfl_robin_chm+ 1
              iclas         = nclas_chm
           end if
        end do
     end do

  end if
!

end subroutine chm_inibcs
