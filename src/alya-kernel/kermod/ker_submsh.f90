!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_submsh.f90
!> @author  Guillaume Houzeaux
!> @date    05/12/2012
!> @brief   Project on support geometry
!> @details Project on support geometry
!> @} 
!-----------------------------------------------------------------------
subroutine ker_submsh()
  use def_kintyp
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_kdtree
  use mod_postpr
  use mod_ker_deform
  use mod_outfor, only : outfor
  use mod_messages, only : livinf
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip)           :: iboun,ipoin,idime,ndefo,kfl_dmeth
  integer(ip)           :: itotn,kpoin,izdom,jpoin
  real(rp)              :: dummr(3),dista,chkdi,dummy,proje(3)
  real(rp)              :: displ_min,displ_max,displ_ave
  type(soltyp), pointer :: solv2(:)

  if( kfl_suppo == 1 ) then

     modul = ID_KERMOD
     call livinf(0_ip,'PROJECT MM BOUNDARY MESH ONTO SUPPORT SURFACE',0_ip)

     if( INOTMASTER ) then
        !
        ! Construct KD-Tree
        !
        chkdi = 1.0e9_rp
        !call kdtree(&
        !     1_ip,mnodb_mm,npoin_mm,nboun_mm,&
        !     coord_mm,lnodb_mm,ltypb_mm)
        !
        ! Look for minimum distance to the surface for all boundary nodes
        ! Displacement is BVESS_SUPPO_KER on these nodes
        !
        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 ) then
              call dpopar(&
                   ndime,coord(1:ndime,ipoin),npoin_mm,mnodb_mm,&
                   nboun_mm,chkdi,ltypb_mm,lnodb_mm,&
                   coord_mm,dista,dummr,proje,iboun) 
              do idime = 1,ndime
                 bvess_suppo_ker(idime,ipoin) = proje(idime)-coord(idime,ipoin)
              end do
           end if
        end do
     end if
     !
     ! Solve mesh displacement
     ! 
     solv2     => momod(modul) % solve(3:)
     ndefo     =  1
     kfl_dmeth =  6  ! Deformaiton method
     call deform_deform(& 
          ndefo,kfl_dmeth,1.0_rp,ID_KERMOD,kfl_fixno_suppo_ker,bvess_suppo_ker,&
          coord,amatr,unkno,rhsid,solv2)
     !
     ! Compute displacement
     !
     if( INOTMASTER ) then
        do ipoin = 1,solve_sol(1) % nequa
           itotn = (ipoin-1) * ndime
           do idime = 1,ndime
              itotn = itotn + 1
              displ_ker(idime,ipoin) = unkno(itotn)
              !
              ! Coordinate update formerly performed in deform_deform
              !
              coord(idime,ipoin) = coord(idime,ipoin) + displ_ker(idime,ipoin)
           end do
        end do
     end if
     !
     ! Compute statistics
     !     
     displ_min =  1.0e12_rp
     displ_max = -1.0e12_rp
     displ_ave =  0.0_rp
     kpoin     =  0
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           if( lpoty(ipoin) /= 0 ) then
              if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
                 dista = 1.0e12_rp
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if( ipoin /= jpoin ) then
                       dummy = 0.0_rp
                       do idime = 1,ndime
                          dummy = dummy + ( coord(idime,ipoin) - coord(idime,jpoin) ) ** 2
                       end do
                       if( dummy < dista ) dista = dummy
                    end if
                 end do
                 dummy = 0.0_rp
                 do idime = 1,ndime
                    dummy = dummy + displ_ker(idime,ipoin) * displ_ker(idime,ipoin)
                 end do
                 kpoin     = kpoin + 1
                 dista     = sqrt( dummy / dista )
                 displ_min = min( displ_min , dista )
                 displ_max = max( displ_max , dista )
                 displ_ave = displ_ave + dista
              end if
           end if
        end do
     end if
     call PAR_SUM(displ_min)
     call PAR_SUM(displ_max)
     call PAR_SUM(displ_ave)
     call PAR_SUM(kpoin)

     displ_ave = displ_ave / real(kpoin,rp)
     routp(1) = displ_min
     routp(2) = displ_max
     routp(3) = displ_ave
     call outfor(58_ip,lun_outpu,' ')
     !
     ! Postprocess
     !
     !wopos(1) = 'XXXXX' 
     !wopos(2) = 'VECTO'
     !wopos(3) = 'NPOIN'
     !call postpr(bvess_suppo_ker,wopos,ittim,cutim)   
     !do ipoin = 1,npoin_mm
     !   write(100,*) ipoin+1040,coord_mm(1,ipoin),coord_mm(2,ipoin)
     !end do
     !do iboun=1,nboun_mm
     !   write(101,*) iboun+1120,lnodb_mm(1,iboun)+1040,lnodb_mm(2,iboun)+1040
     !end do
     !call runend('POPO')
     !
     ! Deallocate memory as we will no longer project onto the support surface
     !
     if( INOTMASTER ) then
        !call kdtree(&
        !     2_ip,mnodb_mm,npoin_mm,nboun_mm,&
        !     coord_mm,lnodb_mm,ltypb_mm)
        call ker_memory(-5_ip) 
        call ker_memory(-6_ip) 
     end if
     kfl_suppo = -1
     !
     ! Recompute arrays depending on mesh coordinates
     !
     modul = ID_KERNEL
     call domarr(2_ip)

  end if

end subroutine ker_submsh
