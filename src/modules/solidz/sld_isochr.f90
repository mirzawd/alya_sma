!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



function sld_zonenode(ipoin)
  use def_kintyp
  use def_domain
  use def_master
  implicit none

  logical :: sld_zonenode
  integer(ip) :: izone,ipoin

  sld_zonenode = .true.
  izone= lzone(ID_SOLIDZ)

  if (nzone > 1) then
 !    if (lpoiz(izone,ipoin) == 0) sld_zonenode = .false.
  end if

end function sld_zonenode

subroutine sld_isochr()
  !-----------------------------------------------------------------------
  !****f* solid/sld_isochr
  ! NAME 
  !    sld_isochr
  ! DESCRIPTION
  !     Computes the isochrones of a variable, this means, the time in
  !     wich the variable reaches a predefined value.
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------

  use      def_parame
  use      def_master
  use      def_domain
  use      def_postpr
  use      def_solidz

  implicit none
  integer(ip)              :: ipoin
  real(rp)                 :: xauxi
  logical(lg)              :: sld_zonenode

  return
  
  if( INOTMASTER) then

     kfl_foten_sld = kfl_foten_sld+1
     if(kfl_foten_sld==1) then
        call sld_calculate_pp_tensors
     endif

     do ipoin= 1,npoin
        if (.not.(sld_zonenode(ipoin))) cycle
 

         !!!!!
         ! displacement isochrone computation
         !!!!!

          if (ndime==1) then
            xauxi = sqrt(displ(1,ipoin,1)**2)
          elseif(ndime==2) then
            xauxi = sqrt(displ(1,ipoin,1)**2+displ(2,ipoin,1)**2)
          elseif(ndime==3) then
            xauxi = sqrt(displ(1,ipoin,1)**2+displ(2,ipoin,1)**2+displ(3,ipoin,1)**2)
          else
            call runend('SLD_ISOCHR: something went wrong. ndime not equal to 1.or.2.or.3')
          endif

          if (xauxi > thiso_sld(1)) then !! Displacement isocrhone
                if (iswav_sld(ipoin,1) == 0) then
                   isoch_sld(ipoin,1) = cutim ! 1 is displacement. 2 is von misses
                   iswav_sld(ipoin,1) = 1_ip ! 1 is for the displacement wave. 2 is for von mises wave.
                end if
          end if


         !!!!!
         ! Von mises stress isochrone computation
         !!!!!
         !
         ! For 2D: sqrt(3*(0.25*(sigxx-sigyy)**2+sigxy**2))
         !
         ! For 3D: sqrt(0.5*((sigxx-sigyy)**2+((sigyy-sigzz)**2+(sigzz-sigxx)**2)+6*(sigxy**2+sigxz**2+sigzy**2))
         !
         !!

          if(ndime==1 .or. ndime==2) then
            xauxi = sqrt(3_rp*(0.25_rp*((caust_sld(1,ipoin)-caust_sld(2,ipoin))**2.0_rp)+&
                               caust_sld(3,ipoin)**2.0_rp))


          elseif(ndime==3) then
            xauxi = (sqrt(0.5_rp*((caust_sld(1,ipoin)-caust_sld(2,ipoin))**2.0_rp+&
                                  (caust_sld(2,ipoin)-caust_sld(3,ipoin))**2.0_rp+&
                                  (caust_sld(3,ipoin)-caust_sld(1,ipoin))**2.0_rp+&               
                          6.0_rp*(caust_sld(4,ipoin)**2.0_rp+caust_sld(5,ipoin)**2.0_rp+&
                                  caust_sld(6,ipoin)**2.0_rp))))
          else
            call runend('SLD_ISOCHR: something went wrong. ndime not equal to 1.or.2.or.3')
          endif


          if (xauxi > thiso_sld(2)) then !! Von mises isochrone
                if (iswav_sld(ipoin,2) == 0) then
                   isoch_sld(ipoin,2) = cutim ! 1 is displacement. 2 is von misses
                   iswav_sld(ipoin,2) = 1_ip ! 1 is for the displacement wave. 2 is for von mises wave.
                end if
          end if
     end do

  end if

end subroutine sld_isochr
