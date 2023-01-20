!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_outvar.f90
!> @author  Mariano Vazquez
!> @brief   This routine bridges output with postpr
!> @date    16/11/1966
!> @details This routine bridges output with postpr
!> @} 
!-----------------------------------------------------------------------
subroutine exm_outvar(ivari,imesh)
  !-----------------------------------------------------------------------
  !****f* exmedi/exm_outvar
  ! NAME 
  !    exm_outvar
  ! DESCRIPTION
  !    This routine bridges output with postpr
  ! USES
  ! USED BY
  !    exm_output
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_postpr

  use      def_exmedi

  use      mod_postpr
  use      mod_iofile
  use      mod_memory

  use      mod_arrays,              only : arrays
  use      mod_arrays,              only : arrays_name
  use      mod_arrays,              only : arrays_number
  use      mod_outvar,              only : outvar
  use      mod_eccoupling   
  
  implicit none
  integer(ip), intent(in)  :: ivari
  integer(ip), intent(in)  :: imesh   !< Mesh to postprocess on
  integer(ip)              :: icurr,iconc
  integer(ip)              :: ipoin,n,kmodel_imate
  real(rp)                 :: xvalp,xvalq,discr,xauxa,xauxb,xauxc,temp1,xvap1,xvap2,xvap3,&
       vap2r,vap2i,xvalv,xvalu,phi,temp2,vauxi
  complex(rp)              :: xvaps(3)
  real(rp)                 :: xvapsr(3)

  real(rp)  , pointer      :: auxvar(:)

  nullify(auxvar)

  select case ( arrays_name(ivari) )  

  case ( 'FIBER' )
     !
     ! Fiber orientation
     !
     if(INOTMASTER) then
        call memgen(zero,ndime,npoin)
        if(kfl_fiber_long_fun.lt.0_ip) then
           gevec(1:ndime,1:npoin) = fiber_exm(1:ndime,1:npoin,1)
        elseif(kfl_fiber_long_fun.eq.1) then
           gevec=0.0_rp
           gevec(1,1:npoin) = 1.0_rp
        elseif(kfl_fiber_long_fun.eq.2) then
           gevec=0.0_rp
           gevec(2,1:npoin) = 1.0_rp
        elseif(kfl_fiber_long_fun.eq.3) then
           gevec=0.0_rp
           gevec(3,1:npoin) = 1.0_rp
        endif
     endif

  case ( 'GRAFI' )
     !
     ! TODO define what this is
     !
     if ( INOTMASTER ) then

        call exm_gravec(fiber,grafi_exm(1:ndime,1:ndime,1:npoin))

        ! calcular los valores propios

        select case ( arrays_name(ivari) )
           
        case ( 'GRAFI' )

           do ipoin=1,npoin
              xauxa = -(grafi_exm(1,1,ipoin) + grafi_exm(2,2,ipoin) + grafi_exm(3,3,ipoin))
              xauxb = grafi_exm(1,1,ipoin)*grafi_exm(2,2,ipoin) + &
                   grafi_exm(1,1,ipoin)*grafi_exm(3,3,ipoin) + &
                   grafi_exm(2,2,ipoin)*grafi_exm(3,3,ipoin) - &
                   grafi_exm(1,2,ipoin)*grafi_exm(2,1,ipoin) - &
                   grafi_exm(1,3,ipoin)*grafi_exm(3,1,ipoin) - &
                   grafi_exm(2,3,ipoin)*grafi_exm(3,2,ipoin) 
              xauxc = grafi_exm(1,1,ipoin)*grafi_exm(2,3,ipoin)*grafi_exm(3,2,ipoin) + &
                   grafi_exm(1,2,ipoin)*grafi_exm(2,1,ipoin)*grafi_exm(3,3,ipoin) + &
                   grafi_exm(1,3,ipoin)*grafi_exm(2,2,ipoin)*grafi_exm(3,1,ipoin) - &
                   grafi_exm(1,1,ipoin)*grafi_exm(2,2,ipoin)*grafi_exm(3,3,ipoin) - &
                   grafi_exm(1,2,ipoin)*grafi_exm(2,3,ipoin)*grafi_exm(3,1,ipoin) - &
                   grafi_exm(1,3,ipoin)*grafi_exm(2,1,ipoin)*grafi_exm(3,2,ipoin)

              xvalp = xauxb - xauxa**2/3.0_rp
              xvalq = (2.0_rp*xauxa**3 - 9.0_rp*xauxa*xauxb + 27.0_rp*xauxc)/27.0_rp
              discr = (xvalp/3.0_rp)**3 + (xvalq/2.0_rp)**2

              if(discr .lt. 0.0_rp)then        ! 3 real roots -- use the trigonometric formulation
                 phi = acos(-xvalq/2.0_rp/sqrt(abs(xvalp*xvalp*xvalp)/27.0_rp))
                 temp1= 2.0_rp*sqrt(abs(xvalp)/3.0_rp)
                 xvap1 =  temp1*cos(phi/3.0_rp)
                 xvap2 = -temp1*cos((phi+pi)/3.0_rp)
                 xvap3 = -temp1*cos((phi-pi)/3.0_rp)
              else                             ! 1 real root & 2 conjugate complex roots or 3 real roots
                 temp1 = -xvalq/2.0_rp + sqrt(discr)
                 temp2 = -xvalq/2.0_rp - sqrt(discr)
                 xvalu = abs(temp1)**(1.0_rp/3.0_rp)
                 xvalv = abs(temp2)**(1.0_rp/3.0_rp)
                 if(temp1 .lt. 0.0_rp) xvalu=-xvalu
                 if(temp2 .lt. 0.0_rp) xvalv=-xvalv
                 xvap1  = xvalu + xvalv
                 vap2r = -(xvalu + xvalv)/2.0_rp
                 vap2i =  (xvalu-xvalv)*sqrt(3.0_rp)/2.0_rp
              end if

              temp1 = xauxa/3.0_rp
              xvap1 = xvap1 - temp1
              xvap2 = xvap2 - temp1
              xvap3 = xvap3 - temp1
              vap2r = vap2r - temp1

              if(discr .lt. 0.0_rp)then
                 xvaps(1) = cmplx( xvap1, 0.0_rp, kind=rp)
                 xvaps(2) = cmplx( xvap2, 0.0_rp, kind=rp)
                 xvaps(3) = cmplx( xvap3, 0.0_rp, kind=rp)
              else if(discr .eq. 0.0_rp)then
                 xvaps(1) = cmplx( xvap1, 0.0_rp, kind=rp)
                 xvaps(2) = cmplx( vap2r, 0.0_rp, kind=rp)
                 xvaps(3) = cmplx( vap2r, 0.0_rp, kind=rp)
              else
                 xvaps(1) = cmplx( xvap1, 0.0_rp, kind=rp)
                 xvaps(2) = cmplx( vap2r, vap2i, kind=rp)
                 xvaps(3) = cmplx( vap2r, -vap2i, kind=rp)
              end if

              ! los guardo en grafi

              if (ndime .eq. 2)then
                 vauxi = sqrt(vmass(ipoin))
              else
                 vauxi = vmass(ipoin)**(1.0_rp/3.0_rp)
              end if

              grafi_exm(1,1,ipoin) = Real(xvaps(1),rp)*vauxi
              grafi_exm(1,2,ipoin) = Real(xvaps(2),rp)*vauxi
              grafi_exm(1,3,ipoin) = Real(xvaps(3),rp)*vauxi
           end do

        case ( 'XXXXX' )

           do ipoin= 1,npoin
              xvapsr(1) = grafi_exm(1,1,ipoin) * grafi_exm(1,1,ipoin) &
                   + grafi_exm(2,1,ipoin) * grafi_exm(2,1,ipoin)
              xvapsr(2) = grafi_exm(1,2,ipoin) * grafi_exm(1,2,ipoin) &
                   + grafi_exm(2,2,ipoin) * grafi_exm(2,2,ipoin)
              if (ndime==3) then
                 xvapsr(1) = xvapsr(1) &
                      + grafi_exm(3,1,ipoin) * grafi_exm(3,1,ipoin) 
                 xvapsr(2) = xvapsr(2) &
                      + grafi_exm(3,2,ipoin) * grafi_exm(3,2,ipoin) 
                 xvapsr(3) = grafi_exm(1,3,ipoin) * grafi_exm(1,3,ipoin) &
                      + grafi_exm(2,3,ipoin) * grafi_exm(2,3,ipoin) &
                      + grafi_exm(3,3,ipoin) * grafi_exm(3,3,ipoin) 
              end if
              grafi_exm(1,1,ipoin) = sqrt(xvapsr(1))
              grafi_exm(1,2,ipoin) = sqrt(xvapsr(2))
              grafi_exm(1,3,ipoin) = sqrt(xvapsr(3))
           end do
        end select

        do ipoin=1, npoin           
           grafi_exm(1:ndime,1:ndime,ipoin) = grafi_exm(1:ndime,1:ndime,ipoin)*vmass(ipoin)**(0.3333_rp)
        end do

        ! los postprocesas
        gevec => grafi_exm(1,1:ndime,1:npoin)

     end if

  case ( 'IOCON' )
     !
     ! All concentrations
     !

     if ( INOTMASTER ) then

        call memgen(zero, nconc_exm, npoin)
        do ipoin = 1,npoin     
           do iconc= 1,nconc_exm
              gevec(iconc,ipoin) = vconc(iconc,ipoin,1)
           end do
        end do
     else
        !
        ! This is necessary for outvar to get the right dimension
        !
        call memgen(0_ip,nconc_exm,1_ip)

     end if

  case ( 'FISOC' )

     call arrays(arrays_number('FISOC'),'POSTPROCESS',fisoc)
     kfl_reset_fisoc = .TRUE.
     return


  case ( 'VAUXI' )
     !
     ! VAUXI_EXM
     !
     call arrays(arrays_number('VAUXI'),'POSTPROCESS',vauxi_exm)
     return

  case ( 'VICEL' )
     !
     ! VICEL_EXM
     !
     call arrays(arrays_number('VICEL'),'POSTPROCESS',vicel_exm)
     return

  case ( 'CALCI' ) 
     !
     ! Calcium concentration
     !
     if ( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR','exm_outvar',auxvar,npoin)
        do ipoin= 1,npoin           
           n = nodemat(ipoin)
           kmodel_imate = kfl_cellmod(n)
           if (kmodel_imate == EXMSLD_CELL_TORORD) then
              auxvar(ipoin) = vconc(5,ipoin,1)
           else
              auxvar(ipoin) = vconc(1,ipoin,1)
           end if
        end do
        gesca => auxvar
     end if

  case ( 'POTAS' )
     !
     ! Calcium concentration
     !
     if ( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR','exm_outvar',auxvar,npoin)
        do ipoin= 1,npoin           
           auxvar(ipoin) = vconc(3,ipoin,1)
        end do
        gesca => auxvar
     end if

  case ( 'EFLUX' )
     !
     ! Electric flux
     !
     if(INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        gevec(1:ndime,1:npoin) = eflux_exm(1:ndime,1:npoin)
     end if

  case ( 'BIPOL' )
     !
     ! Bipolar leads fluxes
     !
     if(INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        gevec(1:ndime,1:npoin) = bipol_exm(1:ndime,1:npoin)
     end if

  case ( 'INTRA' )
     !
     ! ELMAG
     !
     if ( INOTMASTER ) then
        do ipoin= 1,npoin           
           unkno(ipoin) = elmag(ipoin,1)
        end do
        gesca => unkno
     end if

  case ( 'EXTRA' )
     !
     ! Nothing
     !
     call runend('EXM_OUTVAR: NO EXTRACELL POTENTIAL PROGRAMMED')

     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = elmag(ipoin,1)
        end do
        gesca => unkno
     end if

  case ( 'SHEET' )
     !
     ! Sheet fiber direction
     !
     if(INOTMASTER) then
        call memgen(zero,ndime,npoin)
        if(kfl_fiber_tang_fun.lt.0_ip) then
           gevec(1:ndime,1:npoin) = sheet_exm(1:ndime,1:npoin,1)
        elseif(kfl_fiber_tang_fun.eq.1) then
           gevec=0.0_rp
           gevec(1,1:npoin) = 1.0_rp
        elseif(kfl_fiber_tang_fun.eq.2) then
           gevec=0.0_rp
           gevec(2,1:npoin) = 1.0_rp
        elseif(kfl_fiber_tang_fun.eq.3) then
           gevec=0.0_rp
           gevec(3,1:npoin) = 1.0_rp
        endif
     endif

  case ( 'NORMA' )
     !
     ! Normal fiber direction
     !
     if(INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        if(kfl_fiber_norm_fun.lt.0_ip) then
           gevec(1:ndime,1:npoin) = normal_exm(1:ndime,1:npoin,1)
        elseif(kfl_fiber_norm_fun.eq.1) then
           gevec=0.0_rp
           gevec(1,1:npoin) = 1.0_rp
        elseif(kfl_fiber_norm_fun.eq.2) then
           gevec=0.0_rp
           gevec(2,1:npoin) = 1.0_rp
        elseif(kfl_fiber_norm_fun.eq.3) then
           gevec=0.0_rp
           gevec(3,1:npoin) = 1.0_rp
        endif
     endif

  case ( 'CECO1' )
     !
     ! CEDIF component 1
     !
     if ( INOTMASTER ) gevec => cedif_exm(1,:,:)

  case ( 'CECO2' )
     !
     ! CEDIF component 2
     !
     if ( INOTMASTER ) gevec => cedif_exm(2,:,:)

  case ( 'CECO3' )
     !
     ! CEDIF component 3
     !
     if ( INOTMASTER ) gevec => cedif_exm(3,:,:)

  case ( 'DISPL' )
     !
     ! DISPL
     !
     if(INOTEMPTY) then
        call memgen(0_ip,ndime,npoin)
        if( kfl_exmsld_ecc )then
           gevec=displ_lagr_ecc(1:ndime,1:npoin,1)
        else
           gevec= 0.0_rp
        endif
     endif

  case ( 'CURRE' )
     !
     ! Currents
     !
     if ( INOTMASTER ) then

        call memgen(zero, nicel_exm, npoin)
        do ipoin = 1,npoin     
           do icurr= 1,nicel_exm
              gevec(icurr,ipoin) = vicel_exm(icurr,ipoin)
           end do
        end do

     else
        !
        ! This is necessary for outvar to get the right dimension
        !
        call memgen(0_ip,nicel_exm,1_ip)

     end if

  case ( 'QNETO' )
     !
     ! QNETO
     !
     call arrays(arrays_number('QNETO'),'POSTPROCESS',qneto_exm)
     return
 
  case ( 'ISOCH' )
     !
     ! ISOCH
     !
     if ( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR','exm_outvar',auxvar,npoin)
        do ipoin= 1,npoin           
           auxvar(ipoin) = fisoc(1,ipoin)
        end do
        gesca => auxvar
     end if
     kfl_reset_fisoc = .TRUE.

  case ( 'REPOL' )
     !
     ! REPOL
     !
     if ( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR','exm_outvar',auxvar,npoin)
        do ipoin= 1,npoin           
           auxvar(ipoin) = fisoc(2,ipoin)
        end do
        gesca => auxvar
     end if
     kfl_reset_fisoc = .TRUE.

  case default

     return

  end select
  !
  ! Postprocess
  !
  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

  if( INOTMASTER ) then
     nullify(gesca)
     call memory_deallo(mem_modul(1:2,modul),'AUXVAR','exm_outvar',auxvar)
  end if

end subroutine exm_outvar
