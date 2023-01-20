!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Unity_Tests_Toolbox
!> Toolbox for unity tests, to be executed at runtime
!> @{
!> @name    ToolBox for unity tests
!> @file    mod_unity_tests.f90
!> @date    16/01/2017
!> @author  Guillaume Houzeaux
!> @brief   Unity tests
!> @details Unity tests
!
!-----------------------------------------------------------------------

module mod_unity_tests

  use def_kintyp, only : ip,rp,lg
  use def_elmtyp
  use def_domain, only : lrule,ngaus
  use def_domain, only : lexis,nelty
  use def_domain, only : ldime,elmar
  use def_domain, only : lquad,ltopo
  use def_domain, only : memor_dom
  use def_master, only : lun_outpu
  use def_master, only : ioutp,routp,coutp
  use def_master, only : lninv_loc
  use def_master, only : INOTMASTER,IMASTER,INOTSLAVE
  use mod_elmgeo, only : elmgeo_cartesian_derivatives_jacobian
  use mod_elmgeo, only : element_type
  use mod_elmgeo, only : elmgeo_shapf_deriv_heslo
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_outfor, only : outfor
  use mod_rulepw, only : rulepw
  implicit none
  private

  integer(ip), parameter :: max_order=12

  public :: unity_tests_polynomial_integration
  public :: unity_tests_integration_rules
  public :: unity_tests_check_halos
  public :: unity_tests_SpMV
  public :: unity_tests_mod_exchange
  
contains
 
  !-----------------------------------------------------------------------
  !
  !> @brief   Integration rule
  !> @details Test the exactness of selected integration rules
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine unity_tests_integration_rules()

    integer(ip)          :: ielem,inode,ipoin,igaus
    integer(ip)          :: pgaus,pelty,pquad
    integer(ip)          :: iorde,ntens,prule
    real(rp)             :: gpcod(3),rorde
    real(rp)             :: err02(max_order),xinte(max_order),xvolu

    real(rp),    pointer :: xjaci(:,:,:)
    real(rp),    pointer :: xjacm(:,:,:)

    real(rp),    pointer :: gpcrt(:,:,:)  
    real(rp),    pointer :: gpdet(:)  
    real(rp),    pointer :: gpvol(:)  
    real(rp),    pointer :: elcod(:,:)  

    integer(ip)          :: ndime
    integer(ip)          :: pnode
    integer(ip)          :: npoin
    integer(ip)          :: nelem
    integer(ip), pointer :: lnods(:,:)
    real(rp),    pointer :: coord(:,:)  

    call outfor(80_ip,lun_outpu,' ')

    do pelty = 1,nelty

       if( lexis(pelty) /= 0 ) then

          if(      pelty == BAR02 ) then
             call unity_tests_BAR02(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == BAR03 ) then
             call unity_tests_BAR03(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == BAR04 ) then
             call unity_tests_BAR04(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == TRI03 ) then
             call unity_tests_TRI03(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == TRI10 ) then
             call unity_tests_TRI10(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == QUA04 ) then
             call unity_tests_QUA04(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == TRI06 ) then
             call unity_tests_TRI06(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == QUA09 ) then
             call unity_tests_QUA09(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == QUA16 ) then
             cycle
          else if( pelty == TET04 ) then
             call unity_tests_TET04(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == HEX08 ) then
             call unity_tests_HEX08(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == TET10 ) then
             call unity_tests_TET10(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == HEX27 ) then
             call unity_tests_HEX27(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == HEX64 ) then
             cycle
          else if( pelty == PEN06 ) then
             call unity_tests_PEN06(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == PEN18 ) then
             call unity_tests_PEN18(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == PYR05 ) then
             call unity_tests_PYR05(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == TET20 ) then
             call unity_tests_TET20(ndime,pnode,npoin,nelem,lnods,coord)
          else if( pelty == POINT ) then
             goto 10 
          else if( pelty == BAR3D ) then
             goto 10 
          else if( pelty == SHELL ) then
             goto 10 
          else
             call runend('UNITY TEST NOT CODED!')
          end if

          pgaus = ngaus(pelty)
          prule = lrule(pelty)
          pquad = lquad(pelty)
          ndime = ldime(pelty)
          ntens = 3*ndime-3
          
          allocate(xjaci(ndime,ndime,pgaus))
          allocate(xjacm(ndime,ndime,pgaus))
          allocate(gpcrt(ndime,pnode,pgaus))
          allocate(gpdet(pgaus))
          allocate(gpvol(pgaus))
          allocate(elcod(ndime,pnode))

          xvolu = 0.0_rp
          xinte = 0.0_rp

          elements: do ielem = 1,nelem

             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do
             call elmgeo_cartesian_derivatives_jacobian(ndime,pnode,pnode,pgaus,elcod,elmar(pelty) % deriv,xjaci,gpcrt,gpdet)

             gpvol = elmar(pelty) % weigp * gpdet  
             xvolu = xvolu + sum(gpvol)               

             gauss_points: do igaus = 1,pgaus
                gpcod = 0.0_rp
                do inode = 1,pnode
                   gpcod(1:ndime) = gpcod(1:ndime) + elcod(1:ndime,inode) * elmar(pelty) % shape(inode,igaus)
                end do
                do iorde = 1,max_order
                   rorde = real(iorde,rp)
                   xinte(iorde) = xinte(iorde) + gpvol(igaus) * ( gpcod(1)**rorde + gpcod(2)**rorde + gpcod(3)**rorde ) 
                end do

             end do gauss_points

          end do elements
          ! 
          ! Compute error
          !
          iorde = 0
          do while( iorde < max_order )
             iorde        = iorde + 1
             rorde        = real(iorde,rp)
             err02(iorde) = xinte(iorde) - real(ndime,rp) / (rorde+1.0_rp)
             if( abs(err02(iorde)) > 1.0e-12_rp .or. abs(xvolu-1.0_rp) > 1.0e-12_rp ) then
                coutp(1) = trim(element_type(pelty) % name)
                if( pquad == 1 ) then
                   coutp(2) = 'CLOSE'
                else
                   coutp(2) = ' OPEN'
                end if
                ioutp(1) = pgaus
                ioutp(2) = iorde - 1
                call outfor(81_ip,lun_outpu)
                iorde = max_order
             end if
          end do

          deallocate(xjaci)
          deallocate(xjacm)
          deallocate(gpcrt)
          deallocate(gpdet)
          deallocate(gpvol)
          deallocate(elcod) 
          deallocate(lnods)
          deallocate(coord)

10        continue
          
       end if

    end do

  end subroutine unity_tests_integration_rules


  subroutine unity_tests_polynomial_integration()

    integer(ip)          :: ielem,inode,ipoin,igaus
    integer(ip)          :: pgaus,pelty,ielty,ptopo,pquad
    integer(ip)          :: iorde,ierro,ntens,prule
    real(rp)             :: gpcod(3),rorde
    real(rp)             :: err02(max_order),xinte(max_order),xvolu

    real(rp),    pointer :: shapf(:,:)
    real(rp),    pointer :: deriv(:,:,:)
    real(rp),    pointer :: heslo(:,:,:)
    real(rp),    pointer :: posgp(:,:)
    real(rp),    pointer :: weigp(:)
    real(rp),    pointer :: xjaci(:,:,:)
    real(rp),    pointer :: xjacm(:,:,:)

    real(rp),    pointer :: gpcrt(:,:,:)  
    real(rp),    pointer :: gpdet(:)  
    real(rp),    pointer :: gpvol(:)  
    real(rp),    pointer :: elcod(:,:)  

    integer(ip)          :: ndime
    integer(ip)          :: pnode
    integer(ip)          :: npoin
    integer(ip)          :: nelem
    integer(ip), pointer :: lnods(:,:)
    real(rp),    pointer :: coord(:,:)  

    call outfor(78_ip,lun_outpu,' ')

    do ielty = 1,6

       if(      ielty == 1 ) then
          call unity_tests_TET10(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 3
          pelty = TET10
       else if( ielty == 2 ) then
          call unity_tests_HEX27(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 1
          pelty = HEX27
       else if( ielty == 3 ) then
          call unity_tests_PEN06(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 5
          pelty = PEN06
       else if( ielty == 4 ) then
          call unity_tests_PYR05(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 7
          pelty = PYR05
       else if( ielty == 5 ) then
          call unity_tests_TRI06(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 3
          pelty = TRI06
        else if( ielty == 6 ) then
          call unity_tests_QUA09(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 1
          pelty = QUA09
       else if( ielty == 7 ) then
          call unity_tests_BAR02(ndime,pnode,npoin,nelem,lnods,coord)
          prule = -1
          pelty = BAR02
       else if( ielty == 8 ) then
          call unity_tests_BAR03(ndime,pnode,npoin,nelem,lnods,coord)
          prule = -1
          pelty = BAR03 
       else if( ielty == 8 ) then
          call unity_tests_TRI03(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 3
          pelty = TRI03 
       else if( ielty == 9 ) then
          call unity_tests_QUA04(ndime,pnode,npoin,nelem,lnods,coord)
          prule = 1
          pelty = QUA04
       end if 

       ntens = 3*ndime-3
       ptopo = ltopo(pelty)
       pquad = lquad(pelty)
          
       do pgaus = 1,100

          allocate(posgp(ndime,pgaus))
          allocate(weigp(pgaus))
          call rulepw(ndime,pgaus,ptopo,pquad,posgp,weigp,ierro)

          if( ierro == 0 ) then

             allocate(shapf(pnode,pgaus))
             allocate(deriv(ndime,pnode,pgaus))
             allocate(heslo(ntens,pnode,pgaus))
             allocate(xjaci(ndime,ndime,pgaus))
             allocate(xjacm(ndime,ndime,pgaus))
             allocate(gpcrt(ndime,pnode,pgaus))
             allocate(gpdet(pgaus))
             allocate(gpvol(pgaus))
             allocate(elcod(ndime,pnode))
             
             !do igaus = 1,pgaus
             !   call elmgeo_shapf_deriv_heslo(&
             !        ndime,pnode,posgp(:,igaus),shapf(:,igaus),deriv(:,:,igaus),heslo(:,:,igaus),ierro) 
             !end do
             call elmgeo_shapf_deriv_heslo(&
                  ndime,pnode,pgaus,posgp,shapf,deriv,heslo,ierro) 

             !call shafal(posgp,ndime,pnode,pgaus,ntens,shapf,deriv,heslo,ierro)    

             err02 = 0.0_rp
             xvolu = 0.0_rp
             xinte = 0.0_rp

             elements: do ielem = 1,nelem

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do
                call elmgeo_cartesian_derivatives_jacobian(ndime,pnode,pnode,pgaus,elcod,deriv,xjaci,gpcrt,gpdet)

                gpvol = weigp * gpdet  
                xvolu = xvolu + sum(gpvol)               

                gauss_points: do igaus = 1,pgaus
                   gpcod = 0.0_rp
                   do inode = 1,pnode
                      gpcod(1:ndime) = gpcod(1:ndime) + elcod(1:ndime,inode) * shapf(inode,igaus)
                   end do
                   do iorde = 1,max_order
                      rorde = real(iorde,rp)
                      xinte(iorde) = xinte(iorde) + gpvol(igaus) * ( gpcod(1)**rorde + gpcod(2)**rorde + gpcod(3)**rorde ) 
                   end do

                end do gauss_points

             end do elements
             !
             ! Compute error
             !
             do iorde = 1,max_order
                rorde        = real(iorde,rp)
                err02(iorde) = xinte(iorde) - real(ndime,rp) / (rorde+1.0_rp)
                if( abs(err02(iorde)) < 1.0e-12_rp .and. abs(xvolu-1.0_rp) < 1.0e-12_rp ) then
                   coutp(2) = 'OK'
                else
                   coutp(2) = 'FAILED'
                end if
                coutp(1)     = trim(element_type(pelty) % name)
                ioutp(1)     = iorde
                ioutp(2)     = pgaus
                routp(1)     = xvolu
                routp(2)     = err02(iorde)
                call outfor(79_ip,lun_outpu,' ')
             end do

             write(lun_outpu,*)

             deallocate(shapf)
             deallocate(deriv)
             deallocate(heslo)
             deallocate(xjaci)
             deallocate(xjacm)
             deallocate(gpcrt)
             deallocate(gpdet)
             deallocate(gpvol)
             deallocate(elcod) 

          end if

          deallocate(posgp)
          deallocate(weigp)

       end do

       deallocate(lnods)
       deallocate(coord)

    end do

  end subroutine unity_tests_polynomial_integration

  subroutine unity_tests_BAR02(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  1
    nnode =  2
    nelem =  2
    npoin =  3

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 1,2 /)
   lnods(:,2) = (/ 2,3 /)

   coord(:,1) = (/ 0.0_rp /)     
   coord(:,2) = (/ 0.5_rp /)
   coord(:,3) = (/ 1.0_rp /)

 end subroutine unity_tests_BAR02

  subroutine unity_tests_BAR03(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  1
    nnode =  3
    nelem =  1
    npoin =  3

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 1,3,2 /)

   coord(:,1) = (/ 0.0_rp /)     
   coord(:,2) = (/ 0.5_rp /)
   coord(:,3) = (/ 1.0_rp /)

 end subroutine unity_tests_BAR03

  subroutine unity_tests_BAR04(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  1
    nnode =  4
    nelem =  1
    npoin =  4

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 1,3,4,2 /)

   coord(:,1) = (/ 0.0_rp /)     
   coord(:,2) = (/ 1.0_rp/3.0_rp /)
   coord(:,3) = (/ 2.0_rp/3.0_rp /)
   coord(:,4) = (/ 1.0_rp /)

 end subroutine unity_tests_BAR04

  subroutine unity_tests_TRI03(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  2
    nnode =  3
    nelem =  2
    npoin =  4

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 1,2,4 /)
   lnods(:,2) = (/ 2,3,4 /)

   coord(:, 1) = (/ 0.0_rp , 0.0_rp /)     
   coord(:, 2) = (/ 1.0_rp , 0.0_rp /)    
   coord(:, 3) = (/ 1.0_rp , 1.0_rp /)    
   coord(:, 4) = (/ 0.0_rp , 1.0_rp /) 

  end subroutine unity_tests_TRI03

  subroutine unity_tests_QUA04(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  2
    nnode =  4
    nelem =  1
    npoin =  4

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 1,2,3,4 /)

   coord(:, 1) = (/ 0.0_rp , 0.0_rp /)     
   coord(:, 2) = (/ 1.0_rp , 0.0_rp /)    
   coord(:, 3) = (/ 1.0_rp , 1.0_rp /)    
   coord(:, 4) = (/ 0.0_rp , 1.0_rp /) 

 end subroutine unity_tests_QUA04

  subroutine unity_tests_TRI06(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  2
    nnode =  6
    nelem =  4
    npoin = 13

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 13,5 ,9 ,10,  7, 12 /)
   lnods(:,2) = (/ 8 ,5 ,13, 6, 10, 11 /)
   lnods(:,3) = (/ 1 ,5 ,8 , 2,  6,  3 /)
   lnods(:,4) = (/ 9 ,5 ,1 , 7,  2,  4 /)

   coord(:, 1) = (/  0.000000e+00_rp , 1.000000e+00_rp /)     
   coord(:, 2) = (/  2.500000e-01_rp , 7.500000e-01_rp /)    
   coord(:, 3) = (/  5.000000e-01_rp , 1.000000e+00_rp /)    
   coord(:, 4) = (/  0.000000e+00_rp , 5.000000e-01_rp /)    
   coord(:, 5) = (/  5.000000e-01_rp , 5.000000e-01_rp /)    
   coord(:, 6) = (/  7.500000e-01_rp , 7.500000e-01_rp /)    
   coord(:, 7) = (/  2.500000e-01_rp , 2.500000e-01_rp /)    
   coord(:, 8) = (/  1.000000e+00_rp , 1.000000e+00_rp /)    
   coord(:, 9) = (/  0.000000e+00_rp , 0.000000e+00_rp /)    
   coord(:,10) = (/  7.500000e-01_rp , 2.500000e-01_rp /)    
   coord(:,11) = (/  1.000000e+00_rp , 5.000000e-01_rp /)    
   coord(:,12) = (/  5.000000e-01_rp , 0.000000e+00_rp /)    
   coord(:,13) = (/  1.000000e+00_rp , 0.000000e+00_rp /) 

  end subroutine unity_tests_TRI06

  subroutine unity_tests_TRI10(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  
    real(rp)                            :: x13, x23

    ndime =  2
    nnode =  10
    nelem =  2
    npoin = 16

    x13 = 1.0_rp/3.0_rp
    x23 = 2.0_rp/3.0_rp

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

   lnods(:,1) = (/ 1, 2, 4, 5, 6, 14, 16, 11, 12, 13 /)
   lnods(:,2) = (/ 3, 4, 2, 9, 10, 16, 14, 7, 8, 15 /)

   coord(:, 1) = (/  0.0_rp , 0.0_rp /)
   coord(:, 2) = (/  0.0_rp , 1.0_rp /)
   coord(:, 3) = (/  1.0_rp , 1.0_rp /)
   coord(:, 4) = (/  0.0_rp , 1.0_rp /)

   coord(:, 5) = (/  x13 , 0.0_rp /)
   coord(:, 6) = (/  x23 , 0.0_rp /)
   coord(:, 7) = (/  1.0_rp , x13 /)
   coord(:, 8) = (/  1.0_rp , x23 /)
   coord(:, 9) = (/  x23 , 1.0_rp /)
   coord(:,10) = (/  x13 , 0.0_rp /)
   coord(:,11) = (/  0.0_rp , x23 /)
   coord(:,12) = (/  0.0_rp , x13 /)
   coord(:,13) = (/  x13 , x13 /)
   coord(:,14) = (/  x23 , x13 /)
   coord(:,15) = (/  x23 , x23 /)
   coord(:,16) = (/  x13 , x23 /)

  end subroutine unity_tests_TRI10

  subroutine unity_tests_QUA09(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  2
    nnode =  9
    nelem =  1
    npoin =  9

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 9, 6, 1, 5, 7, 3, 2, 8, 4 /)

    coord(:, 1) = (/  0.000000e+00_rp , 1.000000e+00_rp /)     
    coord(:, 2) = (/  0.000000e+00_rp , 5.000000e-01_rp /)    
    coord(:, 3) = (/  5.000000e-01_rp , 1.000000e+00_rp /)    
    coord(:, 4) = (/  5.000000e-01_rp , 5.000000e-01_rp /)    
    coord(:, 5) = (/  0.000000e+00_rp , 0.000000e+00_rp /)    
    coord(:, 6) = (/  1.000000e+00_rp , 1.000000e+00_rp /)    
    coord(:, 7) = (/  1.000000e+00_rp , 5.000000e-01_rp /)    
    coord(:, 8) = (/  5.000000e-01_rp , 0.000000e+00_rp /)    
    coord(:, 9) = (/  1.000000e+00_rp , 0.000000e+00_rp /) 

  end subroutine unity_tests_QUA09

  subroutine unity_tests_TET04(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  3
    nnode =  4
    nelem =  6
    npoin =  8

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 5 , 3 , 1 , 7 /)
    lnods(:,2) = (/ 7 , 4 , 1 , 5 /)
    lnods(:,3) = (/ 2 , 1 , 4 , 5 /)
    lnods(:,4) = (/ 5 , 2 , 6 , 4 /)
    lnods(:,5) = (/ 5 , 8 , 7 , 4 /)
    lnods(:,6) = (/ 8 , 6 , 4 , 5 /)

    coord(:, 1) = (/ 0.000000e+00_rp , 1.000000e+00_rp , 1.000000e+00_rp /)   
    coord(:, 2) = (/ 0.000000e+00_rp , 0.000000e+00_rp , 1.000000e+00_rp /)    
    coord(:, 3) = (/ 0.000000e+00_rp , 1.000000e+00_rp , 0.000000e+00_rp /)    
    coord(:, 4) = (/ 1.000000e+00_rp , 1.000000e+00_rp , 1.000000e+00_rp /)    
    coord(:, 5) = (/ 0.000000e+00_rp , 0.000000e+00_rp , 0.000000e+00_rp /)    
    coord(:, 6) = (/ 1.000000e+00_rp , 0.000000e+00_rp , 1.000000e+00_rp /)    
    coord(:, 7) = (/ 1.000000e+00_rp , 1.000000e+00_rp , 0.000000e+00_rp /)    
    coord(:, 8) = (/ 1.000000e+00_rp , 0.000000e+00_rp , 0.000000e+00_rp /)    
    
  end subroutine unity_tests_TET04

  subroutine unity_tests_TET20(ndime,nnode,npoin,nelem,lnods,coord)

     ! WIP
     !TODO: Add test for TET20

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  3
    nnode =  20
    nelem =  6
    npoin =  64

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) =  (/ 5, 4, 2, 1, 41, 42, 34, 33, 54, 53, 25, 26, 9 , 10, 16, 15, 57, 43, 55, 35 /) 
    lnods(:,2) =  (/ 6, 3, 8, 7, 45, 46, 50, 49, 38, 37, 20, 19, 21, 22, 30, 29, 58, 48, 40, 52 /) 
    lnods(:,3) =  (/ 2, 6, 8, 5, 27, 28, 37, 38, 59, 60, 53, 54, 24, 23, 17, 18, 61, 56, 62, 39 /) 
    lnods(:,4) =  (/ 2, 4, 8, 3, 33, 34, 31, 32, 59, 60, 12, 11, 50, 49, 13, 14, 63, 36, 64, 51 /) 
    lnods(:,5) =  (/ 3, 6, 8, 2, 46, 45, 37, 38, 49, 50, 11, 12, 60, 59, 27, 28, 58, 47, 64, 61 /) 
    lnods(:,6) =  (/ 4, 5, 2, 8, 42, 41, 53, 54, 33, 34, 32, 31, 59, 60, 23, 24, 57, 44, 63, 62 /) 

    coord(:,1 ) = (/ 0.0000000000000000_rp, 0.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,2 ) = (/ 1.0000000000000000_rp, 0.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,3 ) = (/ 1.0000000000000000_rp, 1.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,4 ) = (/ 0.0000000000000000_rp, 1.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,5 ) = (/ 0.0000000000000000_rp, 0.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,6 ) = (/ 1.0000000000000000_rp, 0.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,7 ) = (/ 1.0000000000000000_rp, 1.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,8 ) = (/ 0.0000000000000000_rp, 1.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,40) = (/ 0.6666666666666667_rp, 0.6666666666666666_rp, 0.9999999999999998_rp /) 
    coord(:,9 ) = (/ 0.3333333333333333_rp, 0.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,10) = (/ 0.6666666666666666_rp, 0.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,11) = (/ 1.0000000000000000_rp, 0.3333333333333333_rp, 0.0000000000000000_rp /) 
    coord(:,12) = (/ 1.0000000000000000_rp, 0.6666666666666666_rp, 0.0000000000000000_rp /) 
    coord(:,13) = (/ 0.6666666666666667_rp, 1.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,14) = (/ 0.3333333333333334_rp, 1.0000000000000000_rp, 0.0000000000000000_rp /) 
    coord(:,15) = (/ 0.0000000000000000_rp, 0.6666666666666667_rp, 0.0000000000000000_rp /) 
    coord(:,16) = (/ 0.0000000000000000_rp, 0.3333333333333334_rp, 0.0000000000000000_rp /) 
    coord(:,17) = (/ 0.3333333333333333_rp, 0.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,18) = (/ 0.6666666666666666_rp, 0.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,19) = (/ 1.0000000000000000_rp, 0.3333333333333333_rp, 1.0000000000000000_rp /) 
    coord(:,20) = (/ 1.0000000000000000_rp, 0.6666666666666666_rp, 1.0000000000000000_rp /) 
    coord(:,21) = (/ 0.6666666666666667_rp, 1.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,22) = (/ 0.3333333333333334_rp, 1.0000000000000000_rp, 1.0000000000000000_rp /) 
    coord(:,23) = (/ 0.0000000000000000_rp, 0.6666666666666667_rp, 1.0000000000000000_rp /) 
    coord(:,24) = (/ 0.0000000000000000_rp, 0.3333333333333334_rp, 1.0000000000000000_rp /) 
    coord(:,25) = (/ 0.0000000000000000_rp, 0.0000000000000000_rp, 0.3333333333333333_rp /) 
    coord(:,26) = (/ 0.0000000000000000_rp, 0.0000000000000000_rp, 0.6666666666666666_rp /) 
    coord(:,27) = (/ 1.0000000000000000_rp, 0.0000000000000000_rp, 0.3333333333333333_rp /) 
    coord(:,28) = (/ 1.0000000000000000_rp, 0.0000000000000000_rp, 0.6666666666666666_rp /) 
    coord(:,29) = (/ 1.0000000000000000_rp, 1.0000000000000000_rp, 0.3333333333333333_rp /) 
    coord(:,30) = (/ 1.0000000000000000_rp, 1.0000000000000000_rp, 0.6666666666666666_rp /) 
    coord(:,31) = (/ 0.0000000000000000_rp, 1.0000000000000000_rp, 0.3333333333333333_rp /) 
    coord(:,32) = (/ 0.0000000000000000_rp, 1.0000000000000000_rp, 0.6666666666666666_rp /) 
    coord(:,33) = (/ 0.6666666666666667_rp, 0.3333333333333333_rp, 0.0000000000000000_rp /) 
    coord(:,34) = (/ 0.3333333333333334_rp, 0.6666666666666666_rp, 0.0000000000000000_rp /) 
    coord(:,35) = (/ 0.3333333333333334_rp, 0.3333333333333334_rp, 0.0000000000000000_rp /) 
    coord(:,36) = (/ 0.6666666666666667_rp, 0.6666666666666666_rp, 0.0000000000000000_rp /) 
    coord(:,37) = (/ 0.6666666666666667_rp, 0.3333333333333333_rp, 1.0000000000000000_rp /) 
    coord(:,38) = (/ 0.3333333333333334_rp, 0.6666666666666666_rp, 1.0000000000000000_rp /) 
    coord(:,39) = (/ 0.3333333333333334_rp, 0.3333333333333334_rp, 0.9999999999999998_rp /) 
    coord(:,40) = (/ 0.6666666666666667_rp, 0.6666666666666666_rp, 0.9999999999999998_rp /) 
    coord(:,41) = (/ 0.0000000000000000_rp, 0.3333333333333333_rp, 0.6666666666666667_rp /) 
    coord(:,42) = (/ 0.0000000000000000_rp, 0.6666666666666666_rp, 0.3333333333333334_rp /) 
    coord(:,43) = (/ 0.0000000000000000_rp, 0.3333333333333334_rp, 0.3333333333333334_rp /) 
    coord(:,44) = (/ 0.0000000000000000_rp, 0.6666666666666667_rp, 0.6666666666666667_rp /) 
    coord(:,45) = (/ 1.0000000000000000_rp, 0.3333333333333333_rp, 0.6666666666666667_rp /) 
    coord(:,46) = (/ 1.0000000000000000_rp, 0.6666666666666666_rp, 0.3333333333333334_rp /) 
    coord(:,47) = (/ 0.9999999999999998_rp, 0.3333333333333333_rp, 0.3333333333333334_rp /) 
    coord(:,48) = (/ 0.9999999999999998_rp, 0.6666666666666666_rp, 0.6666666666666667_rp /) 
    coord(:,49) = (/ 0.3333333333333333_rp, 1.0000000000000000_rp, 0.6666666666666667_rp /) 
    coord(:,50) = (/ 0.6666666666666666_rp, 1.0000000000000000_rp, 0.3333333333333334_rp /) 
    coord(:,51) = (/ 0.3333333333333334_rp, 0.9999999999999998_rp, 0.3333333333333334_rp /) 
    coord(:,52) = (/ 0.6666666666666667_rp, 0.9999999999999998_rp, 0.6666666666666667_rp /) 
    coord(:,53) = (/ 0.3333333333333333_rp, 0.0000000000000000_rp, 0.6666666666666667_rp /) 
    coord(:,54) = (/ 0.6666666666666666_rp, 0.0000000000000000_rp, 0.3333333333333334_rp /) 
    coord(:,55) = (/ 0.3333333333333333_rp, 0.0000000000000000_rp, 0.3333333333333334_rp /) 
    coord(:,56) = (/ 0.6666666666666666_rp, 0.0000000000000000_rp, 0.6666666666666667_rp /) 
    coord(:,57) = (/ 0.3333333333333334_rp, 0.3333333333333333_rp, 0.3333333333333334_rp /) 
    coord(:,58) = (/ 0.6666666666666666_rp, 0.6666666666666665_rp, 0.6666666666666667_rp /) 
    coord(:,59) = (/ 0.3333333333333333_rp, 0.6666666666666667_rp, 0.6666666666666667_rp /) 
    coord(:,60) = (/ 0.6666666666666666_rp, 0.3333333333333334_rp, 0.3333333333333334_rp /) 
    coord(:,61) = (/ 0.6666666666666667_rp, 0.3333333333333334_rp, 0.6666666666666667_rp /) 
    coord(:,62) = (/ 0.3333333333333333_rp, 0.3333333333333334_rp, 0.6666666666666667_rp /) 
    coord(:,63) = (/ 0.3333333333333334_rp, 0.6666666666666667_rp, 0.3333333333333334_rp /) 
    coord(:,64) = (/ 0.6666666666666666_rp, 0.6666666666666666_rp, 0.3333333333333334_rp /) 
   
  end subroutine unity_tests_TET20

  subroutine unity_tests_HEX08(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  3
    nnode =  8
    nelem =  1
    npoin =  8

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 1,2,3,4,5,6,7,8 /)

    coord(:, 1) = (/ 0.0_rp , 0.0_rp , 0.0_rp /)
    coord(:, 2) = (/ 1.0_rp , 0.0_rp , 0.0_rp /)   
    coord(:, 3) = (/ 1.0_rp , 1.0_rp , 0.0_rp /)   
    coord(:, 4) = (/ 0.0_rp , 1.0_rp , 0.0_rp /)   
    coord(:, 5) = (/ 0.0_rp , 0.0_rp , 1.0_rp /)   
    coord(:, 6) = (/ 1.0_rp , 0.0_rp , 1.0_rp /)   
    coord(:, 7) = (/ 1.0_rp , 1.0_rp , 1.0_rp /)
    coord(:, 8) = (/ 0.0_rp , 1.0_rp , 1.0_rp /)  

  end subroutine unity_tests_HEX08

  subroutine unity_tests_PYR05(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  3
    nnode =  5
    nelem =  6
    npoin =  9

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 1,5,6,2,9 /)
    lnods(:,2) = (/ 2,6,7,3,9 /)
    lnods(:,3) = (/ 3,7,8,4,9 /)
    lnods(:,4) = (/ 1,4,8,5,9 /)
    lnods(:,5) = (/ 1,2,3,4,9 /)
    lnods(:,6) = (/ 5,8,7,6,9 /)

    coord(:, 1) = (/ 0.0_rp , 0.0_rp , 0.0_rp /)
    coord(:, 2) = (/ 1.0_rp , 0.0_rp , 0.0_rp /)   
    coord(:, 3) = (/ 1.0_rp , 1.0_rp , 0.0_rp /)   
    coord(:, 4) = (/ 0.0_rp , 1.0_rp , 0.0_rp /)   
    coord(:, 5) = (/ 0.0_rp , 0.0_rp , 1.0_rp /)   
    coord(:, 6) = (/ 1.0_rp , 0.0_rp , 1.0_rp /)   
    coord(:, 7) = (/ 1.0_rp , 1.0_rp , 1.0_rp /)
    coord(:, 8) = (/ 0.0_rp , 1.0_rp , 1.0_rp /)  
    coord(:, 9) = (/ 0.5_rp , 0.5_rp , 0.5_rp /)  

  end subroutine unity_tests_PYR05

  subroutine unity_tests_PEN06(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  3
    nnode =  6
    nelem =  2
    npoin =  8

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 1,2,4,5,6,8 /)
    lnods(:,2) = (/ 2,3,4,6,7,8 /)

    coord(:, 1) = (/ 0.0_rp , 0.0_rp , 0.0_rp /)
    coord(:, 2) = (/ 1.0_rp , 0.0_rp , 0.0_rp /)   
    coord(:, 3) = (/ 1.0_rp , 1.0_rp , 0.0_rp /)   
    coord(:, 4) = (/ 0.0_rp , 1.0_rp , 0.0_rp /)   
    coord(:, 5) = (/ 0.0_rp , 0.0_rp , 1.0_rp /)   
    coord(:, 6) = (/ 1.0_rp , 0.0_rp , 1.0_rp /)   
    coord(:, 7) = (/ 1.0_rp , 1.0_rp , 1.0_rp /)
    coord(:, 8) = (/ 0.0_rp , 1.0_rp , 1.0_rp /)  

  end subroutine unity_tests_PEN06

  ! New PEN18 subroutine

  subroutine unity_tests_PEN18(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)

    ndime =  3
    nnode =  18
    nelem =  2
    npoin =  27

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 1 , 2 , 4 , 5 , 6 , 8 ,  9 , 13 , 12 , 26 , 20 , 24 , 14 , 18 , 17 , 19 , 27 , 25 /)
    lnods(:,2) = (/ 2 , 3 , 4 , 6 , 7 , 8 , 10 , 11 , 13 , 20 , 22 , 24 , 15 , 16 , 18 , 21 , 23 , 27 /)

    coord(:, 1) = (/ 0.0_rp , 0.0_rp , 0.0_rp /)
    coord(:, 2) = (/ 1.0_rp , 0.0_rp , 0.0_rp /)
    coord(:, 3) = (/ 1.0_rp , 1.0_rp , 0.0_rp /)
    coord(:, 4) = (/ 0.0_rp , 1.0_rp , 0.0_rp /)
    coord(:, 5) = (/ 0.0_rp , 0.0_rp , 1.0_rp /)
    coord(:, 6) = (/ 1.0_rp , 0.0_rp , 1.0_rp /)
    coord(:, 7) = (/ 1.0_rp , 1.0_rp , 1.0_rp /)
    coord(:, 8) = (/ 0.0_rp , 1.0_rp , 1.0_rp /)
    coord(:, 9) = (/ 0.5_rp , 0.0_rp , 0.0_rp /)
    coord(:,10) = (/ 1.0_rp , 0.5_rp , 0.0_rp /)
    coord(:,11) = (/ 0.5_rp , 1.0_rp , 0.0_rp /)
    coord(:,12) = (/ 0.0_rp , 0.5_rp , 0.0_rp /)
    coord(:,13) = (/ 0.5_rp , 0.5_rp , 0.0_rp /)
    coord(:,14) = (/ 0.5_rp , 0.0_rp , 1.0_rp /)
    coord(:,15) = (/ 1.0_rp , 0.5_rp , 1.0_rp /)
    coord(:,16) = (/ 0.5_rp , 1.0_rp , 1.0_rp /)
    coord(:,17) = (/ 0.0_rp , 0.5_rp , 1.0_rp /)
    coord(:,18) = (/ 0.5_rp , 0.5_rp , 1.0_rp /)
    coord(:,19) = (/ 0.5_rp , 0.0_rp , 0.5_rp /)
    coord(:,20) = (/ 1.0_rp , 0.0_rp , 0.5_rp /)
    coord(:,21) = (/ 1.0_rp , 0.5_rp , 0.5_rp /)
    coord(:,22) = (/ 1.0_rp , 1.0_rp , 0.5_rp /)
    coord(:,23) = (/ 0.5_rp , 1.0_rp , 0.5_rp /)
    coord(:,24) = (/ 0.0_rp , 1.0_rp , 0.5_rp /)
    coord(:,25) = (/ 0.0_rp , 0.5_rp , 0.5_rp /)
    coord(:,26) = (/ 0.0_rp , 0.0_rp , 0.5_rp /)
    coord(:,27) = (/ 0.5_rp , 0.5_rp , 0.5_rp /)

  end subroutine unity_tests_PEN18

  subroutine unity_tests_TET10(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =  3
    nnode = 10
    nelem =  6
    npoin = 27

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/ 1 , 4 , 8 , 3 , 12 , 20 , 25 , 21 , 11 , 24 /)
    lnods(:,2) = (/ 3 , 7 , 8 , 1 , 19 , 15 , 24 , 21 , 27 , 25 /)  
    lnods(:,3) = (/ 5 , 8 , 7 , 1 , 16 , 15 , 26 , 17 , 25 , 27 /)  
    lnods(:,4) = (/ 1 , 5 , 6 , 7 , 17 , 13 , 22 , 27 , 26 , 14 /)  
    lnods(:,5) = (/ 1 , 2 , 3 , 7 ,  9 , 10 , 21 , 27 , 23 , 19 /) 
    lnods(:,6) = (/ 2 , 6 , 7 , 1 , 18 , 14 , 23 ,  9 , 22 , 27 /) 

    coord(:, 1) = (/ 0.000000e+00_rp , 0.000000e+00_rp , 0.000000e+00_rp /)
    coord(:, 2) = (/ 1.000000e+00_rp , 0.000000e+00_rp , 0.000000e+00_rp /)   
    coord(:, 3) = (/ 1.000000e+00_rp , 1.000000e+00_rp , 0.000000e+00_rp /)   
    coord(:, 4) = (/ 0.000000e+00_rp , 1.000000e+00_rp , 0.000000e+00_rp /)   
    coord(:, 5) = (/ 0.000000e+00_rp , 0.000000e+00_rp , 1.000000e+00_rp /)   
    coord(:, 6) = (/ 1.000000e+00_rp , 0.000000e+00_rp , 1.000000e+00_rp /)   
    coord(:, 7) = (/ 1.000000e+00_rp , 1.000000e+00_rp , 1.000000e+00_rp /)
    coord(:, 8) = (/ 0.000000e+00_rp , 1.000000e+00_rp , 1.000000e+00_rp /)   
    coord(:, 9) = (/ 5.000000e-01_rp , 0.000000e+00_rp , 0.000000e+00_rp /)   
    coord(:,10) = (/ 1.000000e+00_rp , 5.000000e-01_rp , 0.000000e+00_rp /)   
    coord(:,11) = (/ 5.000000e-01_rp , 1.000000e+00_rp , 0.000000e+00_rp /)   
    coord(:,12) = (/ 0.000000e+00_rp , 5.000000e-01_rp , 0.000000e+00_rp /)   
    coord(:,13) = (/ 5.000000e-01_rp , 0.000000e+00_rp , 1.000000e+00_rp /) 
    coord(:,14) = (/ 1.000000e+00_rp , 5.000000e-01_rp , 1.000000e+00_rp /)   
    coord(:,15) = (/ 5.000000e-01_rp , 1.000000e+00_rp , 1.000000e+00_rp /)   
    coord(:,16) = (/ 0.000000e+00_rp , 5.000000e-01_rp , 1.000000e+00_rp /)   
    coord(:,17) = (/ 0.000000e+00_rp , 0.000000e+00_rp , 5.000000e-01_rp /)   
    coord(:,18) = (/ 1.000000e+00_rp , 0.000000e+00_rp , 5.000000e-01_rp /)   
    coord(:,19) = (/ 1.000000e+00_rp , 1.000000e+00_rp , 5.000000e-01_rp /)
    coord(:,20) = (/ 0.000000e+00_rp , 1.000000e+00_rp , 5.000000e-01_rp /)   
    coord(:,21) = (/ 5.000000e-01_rp , 5.000000e-01_rp , 0.000000e+00_rp /)   
    coord(:,22) = (/ 5.000000e-01_rp , 0.000000e+00_rp , 5.000000e-01_rp /)   
    coord(:,23) = (/ 1.000000e+00_rp , 5.000000e-01_rp , 5.000000e-01_rp /)   
    coord(:,24) = (/ 5.000000e-01_rp , 1.000000e+00_rp , 5.000000e-01_rp /)   
    coord(:,25) = (/ 0.000000e+00_rp , 5.000000e-01_rp , 5.000000e-01_rp /)   
    coord(:,26) = (/ 5.000000e-01_rp , 5.000000e-01_rp , 1.000000e+00_rp /)   
    coord(:,27) = (/ 5.000000e-01_rp , 5.000000e-01_rp , 5.000000e-01_rp /)  

  end subroutine unity_tests_TET10

  subroutine unity_tests_HEX27(ndime,nnode,npoin,nelem,lnods,coord)

    integer(ip), intent(out)          :: ndime
    integer(ip), intent(out)          :: nnode
    integer(ip), intent(out)          :: npoin
    integer(ip), intent(out)          :: nelem
    integer(ip), intent(inout), pointer :: lnods(:,:)
    real(rp),    intent(inout), pointer :: coord(:,:)  

    ndime =   3
    nnode =  27
    nelem =   8
    npoin = 125

    allocate(lnods(nnode,nelem))
    allocate(coord(ndime,npoin))

    lnods(:,1) = (/106 , 110 , 122 , 118 , 114 , 120 , 125 , 124 , 105 , 109 , 117 , 104 , 103 , 108 , 121 , 116 , 113 , 119 , 123 , 112 , 100 , 101 , 107 , 115 , 102 , 111 , 99  /)
    lnods(:,2) = (/51 , 106 , 118 , 90 , 83 , 114 , 124 , 98 , 19 , 104 , 61 , 49 , 50 , 103 , 116 , 89 , 56 , 112 , 91 , 82 , 17 , 18 , 102 , 60 , 48 , 55 , 16 /)
    lnods(:,3) = (/28 , 64 , 110 , 106 , 80 , 95 , 120 , 114 , 27 , 31 , 105 , 15 , 26 , 62 , 108 , 103 , 79 , 77 , 113 , 54 , 13 , 25 , 29 , 101 , 14 , 53 , 12 /)
    lnods(:,4) = (/73 , 28 , 106 , 51 , 94 , 80 , 114 , 83 , 24 , 15 , 19 , 45 , 72 , 26 , 103 , 50 , 78 , 54 , 56 , 81 , 11 , 23 , 14 , 18 , 44 , 52 , 10 /)
    lnods(:,5) = (/42 , 70 , 96 , 87 , 106 , 110 , 122 , 118 , 38 , 69 , 85 , 40 , 9 , 33 , 84 , 59 , 105 , 109 , 117 , 104 , 37 , 7 , 32 , 58 , 8 , 100 , 6 /)
    lnods(:,6) = (/76 , 42 , 87 , 97 , 51 , 106 , 118 , 90 , 41 , 40 , 86 , 75 , 47 , 9 , 59 , 88 , 19 , 104 , 61 , 49 , 39 , 5 , 8 , 57 , 46 , 17 , 4 /)
    lnods(:,7) = (/67 , 92 , 70 , 42 , 28 , 64 , 110 , 106 , 65 , 68 , 38 , 36 , 22 , 63 , 33 , 9 , 27 , 31 , 105 , 15 , 34 , 21 , 30 , 7 , 3 , 13 , 2 /)
    lnods(:,8) = (/93 , 67 , 42 , 76 , 73 , 28 , 106 , 51 , 66 , 36 , 41 , 74 , 71 , 22 , 9 , 47 , 24 , 15 , 19 , 45 , 35 , 20 , 3 , 5 , 43 , 11 , 1 /)

    coord(:,1  ) = (/   7.500000e-01_rp   , 7.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,2  ) = (/   2.500000e-01_rp   , 7.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,3  ) = (/   5.000000e-01_rp   , 7.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,4  ) = (/   7.500000e-01_rp   , 7.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,5  ) = (/   7.500000e-01_rp   , 7.500000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,6  ) = (/   2.500000e-01_rp   , 7.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,7  ) = (/   2.500000e-01_rp   , 7.500000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,8  ) = (/   5.000000e-01_rp   , 7.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,9  ) = (/   5.000000e-01_rp   , 7.500000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,10 ) = (/   7.500000e-01_rp  ,  2.500000e-01_rp   ,   7.500000e-01_rp /)   
    coord(:,11 ) = (/   7.500000e-01_rp  ,  5.000000e-01_rp   ,   7.500000e-01_rp /)   
    coord(:,12 ) = (/   2.500000e-01_rp  ,  2.500000e-01_rp   ,   7.500000e-01_rp /)   
    coord(:,13 ) = (/   2.500000e-01_rp  ,  5.000000e-01_rp   ,   7.500000e-01_rp /)   
    coord(:,14 ) = (/   5.000000e-01_rp  ,  2.500000e-01_rp   ,   7.500000e-01_rp /)   
    coord(:,15 ) = (/   5.000000e-01_rp  ,  5.000000e-01_rp   ,   7.500000e-01_rp /)   
    coord(:,16 ) = (/   7.500000e-01_rp  ,  2.500000e-01_rp   ,   2.500000e-01_rp /)   
    coord(:,17 ) = (/   7.500000e-01_rp  ,  5.000000e-01_rp   ,   2.500000e-01_rp /)   
    coord(:,18 ) = (/   7.500000e-01_rp  ,  2.500000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,19 ) = (/   7.500000e-01_rp  ,  5.000000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,20 ) = (/   7.500000e-01_rp  ,  7.500000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,21 ) = (/   2.500000e-01_rp  ,  7.500000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,22 ) = (/   5.000000e-01_rp  ,  7.500000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,23 ) = (/   7.500000e-01_rp  ,  2.500000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,24 ) = (/   7.500000e-01_rp  ,  5.000000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,25 ) = (/   2.500000e-01_rp  ,  2.500000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,26 ) = (/   5.000000e-01_rp  ,  2.500000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,27 ) = (/   2.500000e-01_rp  ,  5.000000e-01_rp   ,   1.000000e+00_rp /)   
    coord(:,28 ) = (/   5.000000e-01_rp   , 5.000000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,29 ) = (/   0.000000e+00_rp   , 2.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,30 ) = (/   0.000000e+00_rp   , 7.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,31 ) = (/   0.000000e+00_rp   , 5.000000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,32 ) = (/   0.000000e+00_rp   , 7.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,33 ) = (/   0.000000e+00_rp   , 7.500000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,34 ) = (/   2.500000e-01_rp   , 1.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,35 ) = (/   7.500000e-01_rp   , 1.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,36 ) = (/   5.000000e-01_rp   , 1.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,37 ) = (/   2.500000e-01_rp   , 1.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,38 ) = (/   2.500000e-01_rp   , 1.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,39 ) = (/   7.500000e-01_rp   , 1.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,40 ) = (/   5.000000e-01_rp   , 1.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,41 ) = (/   7.500000e-01_rp   , 1.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,42 ) = (/   5.000000e-01_rp   , 1.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,43 ) = (/   1.000000e+00_rp   , 7.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,44 ) = (/   1.000000e+00_rp   , 2.500000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,45 ) = (/   1.000000e+00_rp   , 5.000000e-01_rp    ,  7.500000e-01_rp /)   
    coord(:,46 ) = (/   1.000000e+00_rp   , 7.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,47 ) = (/   1.000000e+00_rp   , 7.500000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,48 ) = (/   1.000000e+00_rp   , 2.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,49 ) = (/   1.000000e+00_rp   , 5.000000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,50 ) = (/   1.000000e+00_rp   , 2.500000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,51 ) = (/   1.000000e+00_rp   , 5.000000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,52 ) = (/   7.500000e-01_rp   , 0.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,53 ) = (/   2.500000e-01_rp   , 0.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,54 ) = (/   5.000000e-01_rp   , 0.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,55 ) = (/   7.500000e-01_rp   , 0.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,56 ) = (/   7.500000e-01_rp   , 0.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,57 ) = (/   7.500000e-01_rp   , 7.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,58 ) = (/   2.500000e-01_rp   , 7.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,59 ) = (/   5.000000e-01_rp   , 7.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,60 ) = (/   7.500000e-01_rp   , 2.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,61 ) = (/   7.500000e-01_rp   , 5.000000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,62 ) = (/   0.000000e+00_rp   , 2.500000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,63 ) = (/   0.000000e+00_rp   , 7.500000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,64 ) = (/   0.000000e+00_rp   , 5.000000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,65 ) = (/   2.500000e-01_rp   , 1.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,66 ) = (/   7.500000e-01_rp   , 1.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,67 ) = (/   5.000000e-01_rp   , 1.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,68 ) = (/   0.000000e+00_rp   , 1.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,69 ) = (/   0.000000e+00_rp   , 1.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,70 ) = (/   0.000000e+00_rp   , 1.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,71 ) = (/   1.000000e+00_rp   , 7.500000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,72 ) = (/   1.000000e+00_rp   , 2.500000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,73 ) = (/   1.000000e+00_rp   , 5.000000e-01_rp    ,  1.000000e+00_rp /)   
    coord(:,74 ) = (/   1.000000e+00_rp   , 1.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,75 ) = (/   1.000000e+00_rp   , 1.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,76 ) = (/   1.000000e+00_rp   , 1.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,77 ) = (/   0.000000e+00_rp   , 0.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,78 ) = (/   7.500000e-01_rp   , 0.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,79 ) = (/   2.500000e-01_rp   , 0.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,80 ) = (/   5.000000e-01_rp   , 0.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,81 ) = (/   1.000000e+00_rp   , 0.000000e+00_rp    ,  7.500000e-01_rp /)   
    coord(:,82 ) = (/   1.000000e+00_rp   , 0.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,83 ) = (/   1.000000e+00_rp   , 0.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,84 ) = (/   0.000000e+00_rp   , 7.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,85 ) = (/   2.500000e-01_rp   , 1.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,86 ) = (/   7.500000e-01_rp   , 1.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,87 ) = (/   5.000000e-01_rp   , 1.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,88 ) = (/   1.000000e+00_rp   , 7.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,89 ) = (/   1.000000e+00_rp   , 2.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,90 ) = (/   1.000000e+00_rp   , 5.000000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,91 ) = (/   7.500000e-01_rp   , 0.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,92 ) = (/   0.000000e+00_rp   , 1.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,93 ) = (/   1.000000e+00_rp   , 1.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,94 ) = (/   1.000000e+00_rp   , 0.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,95 ) = (/   0.000000e+00_rp   , 0.000000e+00_rp    ,  1.000000e+00_rp /)   
    coord(:,96 ) = (/   0.000000e+00_rp   , 1.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,97 ) = (/   1.000000e+00_rp   , 1.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,98 ) = (/   1.000000e+00_rp   , 0.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,99 ) = (/   2.500000e-01_rp   , 2.500000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,100) = (/   2.500000e-01_rp  ,  5.000000e-01_rp   ,   2.500000e-01_rp /)   
    coord(:,101) = (/   2.500000e-01_rp  ,  2.500000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,102) = (/   5.000000e-01_rp  ,  2.500000e-01_rp   ,   2.500000e-01_rp /)   
    coord(:,103) = (/   5.000000e-01_rp  ,  2.500000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,104) = (/   5.000000e-01_rp  ,  5.000000e-01_rp   ,   2.500000e-01_rp /)   
    coord(:,105) = (/   2.500000e-01_rp  ,  5.000000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,106) = (/   5.000000e-01_rp  ,  5.000000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,107) = (/   0.000000e+00_rp  ,  2.500000e-01_rp   ,   2.500000e-01_rp /)   
    coord(:,108) = (/   0.000000e+00_rp  ,  2.500000e-01_rp   ,   5.000000e-01_rp /)   
    coord(:,109) = (/   0.000000e+00_rp   , 5.000000e-01_rp    ,  2.500000e-01_rp /)   
    coord(:,110) = (/   0.000000e+00_rp   , 5.000000e-01_rp    ,  5.000000e-01_rp /)   
    coord(:,111) = (/   2.500000e-01_rp   , 0.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,112) = (/   5.000000e-01_rp   , 0.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,113) = (/   2.500000e-01_rp   , 0.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,114) = (/   5.000000e-01_rp   , 0.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,115) = (/   2.500000e-01_rp   , 2.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,116) = (/   5.000000e-01_rp   , 2.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,117) = (/   2.500000e-01_rp   , 5.000000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,118) = (/   5.000000e-01_rp   , 5.000000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,119) = (/   0.000000e+00_rp   , 0.000000e+00_rp    ,  2.500000e-01_rp /)   
    coord(:,120) = (/   0.000000e+00_rp   , 0.000000e+00_rp    ,  5.000000e-01_rp /)   
    coord(:,121) = (/   0.000000e+00_rp   , 2.500000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,122) = (/   0.000000e+00_rp   , 5.000000e-01_rp    ,  0.000000e+00_rp /)   
    coord(:,123) = (/   2.500000e-01_rp   , 0.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,124) = (/   5.000000e-01_rp   , 0.000000e+00_rp    ,  0.000000e+00_rp /)   
    coord(:,125) = (/   0.000000e+00_rp   , 0.000000e+00_rp    ,  0.000000e+00_rp /) 

  end subroutine unity_tests_HEX27
 
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    4/05/2017
  !> @brief   Halos unity test
  !> @details Unity tests for halos: check global numbering and
  !>          communications
  !
  !-----------------------------------------------------------------------

  subroutine unity_tests_check_halos()
  
    use def_domain,         only : nboun_2,lnodb
    use def_domain,         only : nboun,nelem_2
    use def_domain,         only : lnnob,lnnod,lmast
    use def_domain,         only : lnods,nelem
    use def_domain,         only : mnode,mnodb
    use mod_communications, only : PAR_FROM_GHOST_BOUNDARY_EXCHANGE
    use mod_communications, only : PAR_FROM_GHOST_ELEMENT_EXCHANGE
    implicit none
    
    integer(ip)          :: iboun,inodb
    integer(ip)          :: ielem,inode,ipoin
    integer(ip), pointer :: lnodb_halo(:,:)
    integer(ip), pointer :: lnods_halo(:,:)
    
    nullify(lnodb_halo)
    nullify(lnods_halo)
    
    if( INOTMASTER ) then
       !
       ! Boundaries
       !
       call memory_alloca(memor_dom,'LNODB_HALO','unity_tests_check_halos',lnodb_halo,mnodb,nboun_2)
       do iboun = 1,nboun_2
          do inodb = 1,lnnob(iboun)
             lnodb_halo(inodb,iboun) = lninv_loc(lnodb(inodb,iboun))
          end do
       end do
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE(lnodb_halo,'SUBSTITUTE','IN MY CODE')
       do iboun = 1,nboun
          do inodb = 1,lnnob(iboun)
             ipoin = lnodb(inodb,iboun)
             if( lmast(ipoin) == 0 ) then
                if( lnodb_halo(inodb,iboun) /= lninv_loc(lnodb(inodb,iboun)) ) then
                   call runend('UNITY_TESTS_CHECK_HALOS: WRONG BOUNDARY HALO')
                end if
             end if
          end do
       end do
       call memory_deallo(memor_dom,'LNODB_HALO','unity_tests_check_halos',lnodb_halo)
       !
       ! Elements
       !
       call memory_alloca(memor_dom,'LNODS_HALO','unity_tests_check_halos',lnods_halo,mnode,nelem_2)
       do ielem = 1,nelem_2
          do inode = 1,lnnod(ielem)
             lnods_halo(inode,ielem) = lninv_loc(lnods(inode,ielem))
          end do
       end do
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE(lnods_halo,'SUBSTITUTE','IN MY CODE')
       do ielem = 1,nelem
          do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             if( lmast(ipoin) == 0 ) then
                if( lnods_halo(inode,ielem) /= lninv_loc(lnods(inode,ielem)) ) then
                   call runend('UNITY_TESTS_CHECK_HALOS: WRONG ELEMENT HALO')
                end if
             end if
          end do
       end do
       call memory_deallo(memor_dom,'LNODS_HALO','unity_tests_check_halos',lnods_halo)       
    end if
    
  end subroutine unity_tests_check_halos

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    4/05/2017
  !> @brief   SpMV
  !> @details SpMV: y = Ax
  !
  !-----------------------------------------------------------------------

  subroutine unity_tests_SpMV()

    use def_domain,         only : npoin,ndime
    use def_domain,         only : coord,vmass
    use def_domain,         only : r_dom,c_dom
    use def_domain,         only : nzdom
    use mod_matrices,       only : matrices_gradient_divergence
    use mod_matrix,         only : matrix_CSR_parallel_SpMV
    use mod_communications, only : PAR_MIN
    use mod_communications, only : PAR_MAX
    
    integer(ip)       :: ipoin
    real(rp), pointer :: xx(:,:)
    real(rp), pointer :: yy(:,:)
    real(rp), pointer :: dy(:,:)
    real(rp), pointer :: dx(:,:)
    real(rp), pointer :: Grad(:,:)
    real(rp), pointer :: Div(:,:)
    real(rp)          :: ymin,ymax
    real(rp)          :: dmin,dmax
    integer(ip)       :: idime
    
    nullify(xx,yy,Grad,Div)
    
    if( INOTMASTER ) then

       allocate(xx(1_ip,npoin))
       allocate(yy(ndime,npoin))
       
       allocate(dx(ndime,npoin))
       allocate(dy(1_ip,npoin))
       
       allocate(Grad(ndime,nzdom))
       allocate(Div(ndime,nzdom))

       xx   = 0.0_rp
       yy   = 0.0_rp
       dx   = 0.0_rp
       dy   = 0.0_rp
       Grad = 0.0_rp 
       Div  = 0.0_rp
       
       xx(1,1:npoin) = 2.3_rp * coord(1,1:npoin) + 3.0_rp * coord(2,1:npoin)
       dx(1,1:npoin) = 2.3_rp * coord(1,1:npoin)
       dx(2,1:npoin) = 3.0_rp * coord(2,1:npoin)
       dx(1,1:npoin) = 1.0_rp
       dx(2,1:npoin) = 0.0_rp

       call matrices_gradient_divergence(Grad,Div)
       call matrix_CSR_parallel_SpMV(npoin,ndime,1_ip,r_dom,c_dom,Grad,xx,yy)
       call matrix_CSR_parallel_SpMV(npoin,1_ip,ndime,r_dom,c_dom,Div,dx,dy)

       do idime = 1,ndime
          yy(idime,1:npoin) = yy(idime,1:npoin) / vmass(1:npoin)
       end do
       dy(1,1:npoin) = dy(1,1:npoin) / vmass(1:npoin)
       
       ymin = minval(yy)
       ymax = maxval(yy)
       dmin = minval(dy)
       dmax = maxval(dy)
       
    end if

    do ipoin = 1,npoin
       print*,ipoin,yy(1,ipoin),yy(2,ipoin),dy(1,ipoin)
    end do
    
    call PAR_MIN(ymin)
    call PAR_MAX(ymax)
    call PAR_MIN(dmin)
    call PAR_MAX(dmax)
    if( INOTSLAVE ) print*,'MIN/MAX GRAD=',ymin,ymax
    if( INOTSLAVE ) print*,'MIN/MAX DIV =',dmin,dmax
    !call runend('O.K.!')
    
    if( INOTMASTER ) then
       deallocate(xx)
       deallocate(yy)
       deallocate(dx)
       deallocate(dy)
       deallocate(Grad)
       deallocate(Div)
    end if
    
  end subroutine unity_tests_SpMV

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-27
  !> @brief   Test for parallel exchange
  !> @details Unity test for broadcasting data
  !> 
  !-----------------------------------------------------------------------
  
  subroutine unity_tests_mod_exchange()

    use mod_exchange
    
    integer(ip)          :: ivar0,ivar1(3),ivar2(2,2),ivarr(0:2)
    integer(ip), pointer :: ipoi1(:)
    real(rp)             :: rvar0,rvar1(3),rvar2(2,2)
    character(5)         :: cvar0,cvar1(2)
    logical(lg)          :: lvar0

    integer(ip)          :: ivas0,ivas1(3),ivas2(2,2),ivass(0:2)
    real(rp)             :: rvas0,rvas1(3),rvas2(2,2)
    character(5)         :: cvas0,cvas1(2)
    logical(lg)          :: lvas0

    integer(ip)          :: ierro

    nullify(ipoi1)
    
    ivar0      = 1
    ivar1      = (/ 1_ip,2_ip,3_ip /)
    ivar2(1,1) = 1
    ivar2(1,2) = 2
    ivar2(2,1) = 3
    ivar2(2,2) = 4
    ivarr      = (/  1_ip,2_ip,3_ip /)
    rvar0      = 1
    rvar1      = (/ 1.0_rp,2.0_rp,3.0_rp /)
    rvar2(1,1) = 1.0_rp
    rvar2(1,2) = 2.0_rp
    rvar2(2,1) = 3.0_rp
    rvar2(2,2) = 4.0_rp
    lvar0      = .true.
    cvar0      = 'ABCDE'
    cvar1      = (/ 'ABCDE','FGHIK' /)

    if( IMASTER ) then
       ivas0 = ivar0
       ivas1 = ivar1
       ivas2 = ivar2
       ivass = ivarr
       rvas0 = rvar0
       rvas1 = rvar1
       rvas2 = rvar2
       cvas0 = cvar0
       cvas1 = cvar1
       lvas0 = lvar0
    end if

    call exchange_add(ivas0)
    call exchange_add(ivas1)
    call exchange_add(ivas2)
    call exchange_add(ivass)
    call exchange_add(rvas0)
    call exchange_add(rvas1)
    call exchange_add(rvas2)
    call exchange_add(cvas0)
    call exchange_add(cvas1)
    call exchange_add(lvas0)
    call exchange_end()
    
    ierro = 0
    
    if( ivas0 /= ivar0 )        ierro = 1 
    if( any( ivas1 /= ivar1 ) ) ierro = 2
    if( any( ivas2 /= ivar2 ) ) ierro = 3
    if( any( ivass /= ivarr ) ) ierro = 3
    if( rvas0 /= rvar0 )        ierro = 4
    if( any( rvas1 /= rvar1 ) ) ierro = 5
    if( any( rvas2 /= rvar2 ) ) ierro = 6
    if( .not. lvar0           ) ierro = 7
    if( cvas0 /= cvar0 )        ierro = 8
    if( any( cvar1 /= cvas1 ) ) ierro = 9

    if( ierro /= 0 ) call runend('MOD_UNITY_TESTS: EXCHANGE FAILS')

  end subroutine unity_tests_mod_exchange

end module mod_unity_tests
!> @}                                                       
