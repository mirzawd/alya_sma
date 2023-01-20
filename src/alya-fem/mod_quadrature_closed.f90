!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_quadrature_gauss.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Close quadrature
!> @details Close quadrature rule
!>   
!-----------------------------------------------------------------------

module mod_quadrature_closed

  use def_kintyp_basic, only : ip,rp
  implicit none
  private

  public :: close_poi
  public :: close_bar
  public :: close_tri
  public :: close_qua
  public :: close_pen
  public :: close_pyr

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Point integration
  !> @details Point integration, here only for consistency
  !> 
  !-----------------------------------------------------------------------

  subroutine close_poi(ngaus,weigp,ierro)

    integer(ip), intent(in)  :: ngaus
    integer(ip), intent(out) :: ierro
    real(rp),    intent(out) :: weigp(ngaus)

    ierro    = 0
    weigp(1) = 1.0_rp

  end subroutine close_poi

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Closed rule for bars
  !> @details This routine sets up the integration constants of closed rules in
  !>          one dimension.
  !>
  !-----------------------------------------------------------------------

  subroutine close_bar(ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ngaus
    integer(ip), intent(out) :: ierro
    real(rp),    intent(out) :: posgp(ngaus),weigp(ngaus)

    ierro=0

    if(ngaus==2) then
       posgp(1)=-1.0_rp
       posgp(2)= 1.0_rp
       weigp(1)= 1.0_rp
       weigp(2)= 1.0_rp

    else if(ngaus==3) then
       posgp(1)=-1.0_rp
       posgp(2)= 0.0_rp
       posgp(3)= 1.0_rp
       weigp(1)= 1.0_rp/3.0_rp
       weigp(2)= 4.0_rp/3.0_rp
       weigp(3)= 1.0_rp/3.0_rp

    else if(ngaus==4) then
       posgp(1)=-1.0_rp
       posgp(2)=-1.0_rp/3.0_rp
       posgp(3)= 1.0_rp/3.0_rp
       posgp(4)= 1.0_rp
       weigp(1)= 1.0_rp/4.0_rp
       weigp(2)= 3.0_rp/4.0_rp
       weigp(3)= 3.0_rp/4.0_rp
       weigp(4)= 1.0_rp/4.0_rp

    else

       ierro=1

    end if

  end subroutine close_bar

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Closed rule for hex
  !> @details This routine sets up the integration constants of closed 
  !>          integration rules for brick elements:
  !>
  !>               NDIME = 1             NDIME = 2             NDIME = 3
  !>
  !>           NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
  !>           -----  ----------     -----  ----------     -----  ----------  
  !>             2      q1           2 x 2     q1          2x2x2     q1   
  !>             3      q3           3 x 3     q3          3x3x3     q3
  !>             4      q3           4 x 4     q3          4x4x4     q3
  !>                                   8       p2            20      p2 
  !>
  !-----------------------------------------------------------------------

  subroutine close_qua(ndime,ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ndime,ngaus
    integer(ip), intent(out) :: ierro
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: posgl(4), weigl(4)
    integer(ip)              :: inoga(64),igaus,nlocs,ilocs,jlocs,klocs

    ierro=0
    !
    ! 2D 8-point integration
    !
    if(ndime==2.and.ngaus==8) then
       do igaus=5,8
          weigp(igaus-4)=-1.0_rp/3.0_rp
          weigp(igaus  )= 4.0_rp/3.0_rp
       end do
       posgp(1,5)= 0.0_rp 
       posgp(1,6)= 1.0_rp
       posgp(1,7)= 0.0_rp
       posgp(1,8)=-1.0_rp
       posgp(2,5)=-1.0_rp
       posgp(2,6)= 0.0_rp
       posgp(2,7)= 1.0_rp
       posgp(2,8)= 0.0_rp
       return
    end if
    !
    ! 3D 20-point integration
    !
    if(ndime==3.and.ngaus==20) then
       nlocs=nint(real(ngaus, rp)**(1.0_rp/3.0_rp))
       do igaus=1,8
          weigp(igaus)=-1.0_rp
       end do
       do igaus=9,20
          weigp(igaus)= 4.0_rp/3.0_rp
       end do
       posgl(1)=-1.0_rp
       posgl(2)= 0.0_rp
       posgl(3)= 1.0_rp
       call close_chaord(inoga,27_ip)
       igaus=0
       do ilocs=1,nlocs
          do jlocs=1,nlocs
             do klocs=1,nlocs
                igaus=igaus+1
                if(inoga(igaus) <= 20) then
                   posgp(1,inoga(igaus))=posgl(ilocs)
                   posgp(2,inoga(igaus))=posgl(jlocs)
                   posgp(3,inoga(igaus))=posgl(klocs)
                end if
             end do
          end do
       end do
       return
    end if
    !
    ! Rules obtained from one-dimensional integration
    !
    if(ndime==2) then
       nlocs=nint(sqrt(real(ngaus, rp)))
    else if(ndime==3) then
       nlocs=nint(real(ngaus, rp)**(1.0_rp/3.0_rp))
    end if

    ilocs = nlocs**ndime
    call close_chaord(inoga,ilocs)

    if(nlocs==2) then
       posgl(1)=-1.0_rp
       posgl(2)= 1.0_rp
       weigl(1)= 1.0_rp
       weigl(2)= 1.0_rp
    else if(nlocs==3) then
       posgl(1)=-1.0_rp
       posgl(2)= 0.0_rp
       posgl(3)= 1.0_rp
       weigl(1)= 1.0_rp/3.0_rp
       weigl(2)= 4.0_rp/3.0_rp
       weigl(3)= 1.0_rp/3.0_rp
    else if(nlocs==4) then
       posgl(1)=-1.0_rp
       posgl(2)=-1.0_rp/3.0_rp
       posgl(3)= 1.0_rp/3.0_rp
       posgl(4)= 1.0_rp
       weigl(1)= 1.0_rp/4.0_rp
       weigl(2)= 3.0_rp/4.0_rp
       weigl(3)= 3.0_rp/4.0_rp
       weigl(4)= 1.0_rp/4.0_rp
    else
       ierro=1
    end if

    if(ndime==2) then
       igaus=0
       do ilocs=1,nlocs
          do jlocs=1,nlocs
             igaus=igaus+1
             weigp(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)
             posgp(1,inoga(igaus))=posgl(ilocs)
             posgp(2,inoga(igaus))=posgl(jlocs)
          end do
       end do
       if(ngaus==5) then                                    !  4       3
          do igaus=1,4                                      !
             weigp(igaus)=1.0_rp/3.0_rp                     !      5
          end do                                            !
          weigp(  5)=8.0_rp/3.0_rp                          !  1       2
          posgp(1,5)=0.0_rp
          posgp(2,5)=0.0_rp
       end if
    else if(ndime==3) then
       igaus=0
       do ilocs=1,nlocs
          do jlocs=1,nlocs
             do klocs=1,nlocs
                igaus=igaus+1
                weigp(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
                posgp(1,inoga(igaus))=posgl(ilocs)
                posgp(2,inoga(igaus))=posgl(jlocs)
                posgp(3,inoga(igaus))=posgl(klocs)
             end do
          end do
       end do
       if(ngaus==9) then                                 !    8------7
          do igaus=1,8                                   !   /|     /|
             weigp(igaus)=weigp(igaus)/3.0_rp            !  5------6 |
          end do                                         !  | | 9  | |
          weigp(  9)=16.0_rp/3.0_rp                      !  | 4----|-3
          posgp(1,9)=0.0_rp                              !  |/     |/
          posgp(2,9)=0.0_rp                              !  1------2
          posgp(3,9)=0.0_rp
       end if
    end if

  end subroutine close_qua

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Close rule for tet
  !> @details This routine sets up the integration constants of closed rules for
  !>          triangles and tetrahedra 
  !> 
  !>                 NDIME = 2             NDIME = 3
  !> 
  !>              NGAUS  EXACT POL.     NGAUS  EXACT POL. 
  !>              -----  ----------     -----  ----------
  !>                3       p1            4       p1
  !>                4       p2            5       p2
  !>                6       p1           10       p2
  !>                7       p3           11       p3
  !>               10       p1           15       p3
  !>                                     20       p2 & x^4,...      
  !> 
  !-----------------------------------------------------------------------

  subroutine close_tri(ndime,ngaus,posgp,weigp,ierro)

    integer(ip), intent(out) :: ierro
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: w1,w2,w3
    !
    ! Line integral
    !
    ierro=0

    if(ndime==2) then
       !
       ! Area integral (triangles)
       !
       if(ngaus==3) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          weigp(  1)= 1.0_rp/6.0_rp
          weigp(  2)= 1.0_rp/6.0_rp
          weigp(  3)= 1.0_rp/6.0_rp
       else if(ngaus==4) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(1,4)= 1.0_rp/3.0_rp
          posgp(2,4)= 1.0_rp/3.0_rp
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 3.0_rp/8.0_rp
       else if(ngaus==6) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(1,4)= 0.5_rp
          posgp(2,4)= 0.0_rp
          posgp(1,5)= 0.5_rp
          posgp(2,5)= 0.5_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 0.5_rp
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 1.0_rp/8.0_rp
          weigp(  5)= 1.0_rp/8.0_rp
          weigp(  6)= 1.0_rp/8.0_rp
       else if(ngaus==7) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(1,4)= 0.5_rp
          posgp(2,4)= 0.0_rp
          posgp(1,5)= 0.5_rp
          posgp(2,5)= 0.5_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 0.5_rp
          posgp(1,7)= 1.0_rp/3.0_rp
          posgp(2,7)= 1.0_rp/3.0_rp
          weigp(  1)= 1.0_rp/40.0_rp
          weigp(  2)= 1.0_rp/40.0_rp
          weigp(  3)= 1.0_rp/40.0_rp
          weigp(  4)= 1.0_rp/15.0_rp
          weigp(  5)= 1.0_rp/15.0_rp
          weigp(  6)= 1.0_rp/15.0_rp
          weigp(  7)= 9.0_rp/40.0_rp
       else if(ngaus==10) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(1, 4)= 1.0_rp/3.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(1, 5)= 2.0_rp/3.0_rp
          posgp(2, 5)= 0.0_rp
          posgp(1, 6)= 2.0_rp/3.0_rp 
          posgp(2, 6)= 1.0_rp/3.0_rp
          posgp(1, 7)= 1.0_rp/3.0_rp
          posgp(2, 7)= 2.0_rp/3.0_rp
          posgp(1, 8)= 0.0_rp
          posgp(2, 8)= 2.0_rp/3.0_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 1.0_rp/3.0_rp
          posgp(1,10)= 1.0_rp/3.0_rp
          posgp(2,10)= 1.0_rp/3.0_rp
          weigp(   1)= 1.0_rp/54.0_rp
          weigp(   2)= 1.0_rp/54.0_rp
          weigp(   3)= 1.0_rp/54.0_rp
          weigp(   4)= 1.0_rp/18.0_rp
          weigp(   5)= 1.0_rp/18.0_rp
          weigp(   6)= 1.0_rp/18.0_rp
          weigp(   7)= 1.0_rp/18.0_rp
          weigp(   8)= 1.0_rp/18.0_rp 
          weigp(   9)= 1.0_rp/18.0_rp 
          weigp(  10)= 1.0_rp/9.0_rp
       else
          ierro=1
       end if

    else if(ndime==3) then
       !
       ! Volume integral ( tetrahedra )
       !
       if(ngaus==4) then
          posgp(1,1)= 0.0_rp 
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp
          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 1.0_rp/24.0_rp
       else if(ngaus==5) then
          posgp(1,1)= 0.0_rp 
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp
          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          posgp(1,5)= 1.0_rp/4.0_rp
          posgp(2,5)= 1.0_rp/4.0_rp
          posgp(3,5)= 1.0_rp/4.0_rp
          weigp(  1)= 1.0_rp/120.0_rp
          weigp(  2)= 1.0_rp/120.0_rp
          weigp(  3)= 1.0_rp/120.0_rp
          weigp(  4)= 1.0_rp/120.0_rp
          weigp(  5)= 2.0_rp/15.0_rp
       else if(ngaus==10) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(3, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 0.5_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 0.5_rp
          posgp(2, 6)= 0.5_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 0.0_rp
          posgp(2, 7)= 0.5_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 0.5_rp
          posgp(2, 8)= 0.0_rp
          posgp(3, 8)= 0.5_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 0.5_rp
          posgp(3, 9)= 0.5_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.0_rp
          posgp(3,10)= 0.5_rp
          !weigp(   1)=-1.0_rp/120.0_rp
          !weigp(   2)=-1.0_rp/120.0_rp 
          !weigp(   3)=-1.0_rp/120.0_rp 
          !weigp(   4)=-1.0_rp/120.0_rp 
          !weigp(   5)= 1.0_rp/30.0_rp 
          !weigp(   6)= 1.0_rp/30.0_rp 
          !weigp(   7)= 1.0_rp/30.0_rp 
          !weigp(   8)= 1.0_rp/30.0_rp 
          !weigp(   9)= 1.0_rp/30.0_rp 
          !weigp(  10)= 1.0_rp/30.0_rp 
          weigp(   1)=1.0_rp/240.0_rp
          weigp(   2)=1.0_rp/240.0_rp 
          weigp(   3)=1.0_rp/240.0_rp 
          weigp(   4)=1.0_rp/240.0_rp 
          weigp(   5)= 1.0_rp/40.0_rp 
          weigp(   6)= 1.0_rp/40.0_rp 
          weigp(   7)= 1.0_rp/40.0_rp 
          weigp(   8)= 1.0_rp/40.0_rp 
          weigp(   9)= 1.0_rp/40.0_rp 
          weigp(  10)= 1.0_rp/40.0_rp 
       else if(ngaus==11) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(3, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 0.5_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 0.5_rp
          posgp(2, 6)= 0.5_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 0.0_rp
          posgp(2, 7)= 0.5_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 0.5_rp
          posgp(2, 8)= 0.0_rp
          posgp(3, 8)= 0.5_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 0.5_rp
          posgp(3, 9)= 0.5_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.0_rp
          posgp(3,10)= 0.5_rp
          posgp(1,11)= 1.0_rp/4.0_rp
          posgp(2,11)= 1.0_rp/4.0_rp
          posgp(3,11)= 1.0_rp/4.0_rp
          weigp(   1)= 1.0_rp/360.0_rp
          weigp(   2)= 1.0_rp/360.0_rp
          weigp(   3)= 1.0_rp/360.0_rp
          weigp(   4)= 1.0_rp/360.0_rp
          weigp(   5)= 1.0_rp/90.0_rp
          weigp(   6)= 1.0_rp/90.0_rp
          weigp(   7)= 1.0_rp/90.0_rp
          weigp(   8)= 1.0_rp/90.0_rp
          weigp(   9)= 1.0_rp/90.0_rp
          weigp(  10)= 1.0_rp/90.0_rp
          weigp(  11)= 4.0_rp/45.0_rp
       else if(ngaus==15) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(3, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 0.5_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 0.5_rp
          posgp(2, 6)= 0.5_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 0.0_rp
          posgp(2, 7)= 0.5_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 0.0_rp
          posgp(2, 8)= 0.0_rp
          posgp(3, 8)= 0.5_rp
          posgp(1, 9)= 0.5_rp
          posgp(2, 9)= 0.0_rp
          posgp(3, 9)= 0.5_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.5_rp
          posgp(3,10)= 0.5_rp
          posgp(1,11)= 1.0_rp/3.0_rp
          posgp(2,11)= 1.0_rp/3.0_rp
          posgp(3,11)= 0.0_rp
          posgp(1,12)= 1.0_rp/3.0_rp
          posgp(2,12)= 0.0_rp
          posgp(3,12)= 1.0_rp/3.0_rp
          posgp(1,13)= 1.0_rp/3.0_rp
          posgp(2,13)= 1.0_rp/3.0_rp
          posgp(3,13)= 1.0_rp/3.0_rp
          posgp(1,14)= 0.0_rp
          posgp(2,14)= 1.0_rp/3.0_rp
          posgp(3,14)= 1.0_rp/3.0_rp
          posgp(1,15)= 1.0_rp/4.0_rp
          posgp(2,15)= 1.0_rp/4.0_rp
          posgp(3,15)= 1.0_rp/4.0_rp
          weigp(   1)= 17.0_rp/5040.0_rp
          weigp(   2)= 17.0_rp/5040.0_rp
          weigp(   3)= 17.0_rp/5040.0_rp
          weigp(   4)= 17.0_rp/5040.0_rp
          weigp(   5)=  2.0_rp/315.0_rp
          weigp(   6)=  2.0_rp/315.0_rp
          weigp(   7)=  2.0_rp/315.0_rp
          weigp(   8)=  2.0_rp/315.0_rp
          weigp(   9)=  2.0_rp/315.0_rp
          weigp(  10)=  2.0_rp/315.0_rp
          weigp(  11)=  9.0_rp/840.0_rp
          weigp(  12)=  9.0_rp/840.0_rp
          weigp(  13)=  9.0_rp/840.0_rp
          weigp(  14)=  9.0_rp/840.0_rp
          weigp(  15)= 16.0_rp/315.0_rp
       else if(ngaus==20) then                           ! Integrates P2
          posgp(1, 1)= 0.0_rp                                ! and quartic mono-
          posgp(2, 1)= 0.0_rp                                ! mials to avoid zero 
          posgp(3, 1)= 0.0_rp                                ! weights at the 
          posgp(1, 2)= 1.0_rp                                ! edges
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 1.0_rp/3.0_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 2.0_rp/3.0_rp
          posgp(2, 6)= 0.0_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 2.0_rp/3.0_rp
          posgp(2, 7)= 1.0_rp/3.0_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 1.0_rp/3.0_rp
          posgp(2, 8)= 2.0_rp/3.0_rp
          posgp(3, 8)= 0.0_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 2.0_rp/3.0_rp
          posgp(3, 9)= 0.0_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 1.0_rp/3.0_rp
          posgp(3,10)= 0.0_rp
          posgp(1,11)= 0.0_rp
          posgp(2,11)= 0.0_rp
          posgp(3,11)= 1.0_rp/3.0_rp
          posgp(1,12)= 2.0_rp/3.0_rp
          posgp(2,12)= 0.0_rp
          posgp(3,12)= 1.0_rp/3.0_rp
          posgp(1,13)= 0.0_rp
          posgp(2,13)= 2.0_rp/3.0_rp
          posgp(3,13)= 1.0_rp/3.0_rp
          posgp(1,14)= 0.0_rp
          posgp(2,14)= 0.0_rp
          posgp(3,14)= 2.0_rp/3.0_rp
          posgp(1,15)= 1.0_rp/3.0_rp
          posgp(2,15)= 0.0_rp
          posgp(3,15)= 2.0_rp/3.0_rp
          posgp(1,16)= 0.0_rp
          posgp(2,16)= 1.0_rp/3.0_rp
          posgp(3,16)= 2.0_rp/3.0_rp
          posgp(1,17)= 1.0_rp/3.0_rp
          posgp(2,17)= 1.0_rp/3.0_rp
          posgp(3,17)= 0.0_rp
          posgp(1,18)= 1.0_rp/3.0_rp
          posgp(2,18)= 0.0_rp
          posgp(3,18)= 1.0_rp/3.0_rp
          posgp(1,19)= 1.0_rp/3.0_rp
          posgp(2,19)= 1.0_rp/3.0_rp
          posgp(3,19)= 1.0_rp/3.0_rp
          posgp(1,20)= 0.0_rp
          posgp(2,20)= 1.0_rp/3.0_rp
          posgp(3,20)= 1.0_rp/3.0_rp
          w1 = 383.0_rp/2400.0_rp
          w2 = 1.0_rp/240.0_rp - w1
          w3 = 3.0_rp/80.0_rp  - w2
          weigp(   1)= w1
          weigp(   2)= w1
          weigp(   3)= w1
          weigp(   4)= w1
          weigp(   5)= w2
          weigp(   6)= w2
          weigp(   7)= w2
          weigp(   8)= w2
          weigp(   9)= w2
          weigp(  10)= w2
          weigp(  11)= w2
          weigp(  12)= w2
          weigp(  13)= w2
          weigp(  14)= w2
          weigp(  15)= w2
          weigp(  16)= w2
          weigp(  17)= w3
          weigp(  18)= w3
          weigp(  19)= w3
          weigp(  20)= w3
       else
          ierro=1
       end if
    end if
    !
    ! Errors
    !
    !if(ierro==1) call runend('CLOSE_TRI: NOT AVAILABLE INTEGRATION RULE')

  end subroutine close_tri

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Closed rule for prisms
  !> @details This routine sets up the integration constants of closed rules for
  !           PRISMS
  ! 
  !                NDIME = 3    
  ! 
  !             NGAUS  EXACT POL.
  !             -----  ----------  
  !               3       p1       
  ! 
  !-----------------------------------------------------------------------

  subroutine close_pen(ndime,ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ndime,ngaus
    integer(ip), intent(out) :: ierro
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)

    ierro=0
    !
    ! Area integral
    !
    if(ndime==3) then
       if(ngaus==6) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp
          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          posgp(1,5)= 1.0_rp
          posgp(2,5)= 0.0_rp
          posgp(3,5)= 1.0_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 1.0_rp
          posgp(3,6)= 1.0_rp
          weigp(  1)= 1.0_rp/12.0_rp
          weigp(  2)= 1.0_rp/12.0_rp
          weigp(  3)= 1.0_rp/12.0_rp
          weigp(  4)= 1.0_rp/12.0_rp
          weigp(  5)= 1.0_rp/12.0_rp
          weigp(  6)= 1.0_rp/12.0_rp

       else if(ngaus==18) then

          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp

          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          posgp(1,5)= 1.0_rp
          posgp(2,5)= 0.0_rp
          posgp(3,5)= 1.0_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 1.0_rp
          posgp(3,6)= 1.0_rp

          ! 2nd order nodes follow GiD,
          ! just in case...

          posgp(1,7)= 0.5_rp
          posgp(2,7)= 0.0_rp
          posgp(3,7)= 0.0_rp
          posgp(1,8)= 0.5_rp
          posgp(2,8)= 0.5_rp
          posgp(3,8)= 0.0_rp
          posgp(1,9)= 0.0_rp
          posgp(2,9)= 0.5_rp
          posgp(3,9)= 0.0_rp

          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.0_rp
          posgp(3,10)= 0.5_rp
          posgp(1,11)= 1.0_rp
          posgp(2,11)= 0.0_rp
          posgp(3,11)= 0.5_rp
          posgp(1,12)= 0.0_rp
          posgp(2,12)= 1.0_rp
          posgp(3,12)= 0.5_rp

          posgp(1,13)= 0.5_rp
          posgp(2,13)= 0.0_rp
          posgp(3,13)= 1.0_rp
          posgp(1,14)= 0.5_rp
          posgp(2,14)= 0.5_rp
          posgp(3,14)= 1.0_rp
          posgp(1,15)= 0.0_rp
          posgp(2,15)= 0.5_rp
          posgp(3,15)= 1.0_rp

          posgp(1,16)= 0.5_rp
          posgp(2,16)= 0.0_rp
          posgp(3,16)= 0.5_rp
          posgp(1,17)= 0.5_rp
          posgp(2,17)= 0.5_rp
          posgp(3,17)= 0.5_rp
          posgp(1,18)= 0.0_rp
          posgp(2,18)= 0.5_rp
          posgp(3,18)= 0.5_rp

          weigp(  1)= 1.0_rp/48.0_rp
          weigp(  2)= 1.0_rp/48.0_rp
          weigp(  3)= 1.0_rp/48.0_rp
          weigp(  4)= 1.0_rp/48.0_rp
          weigp(  5)= 1.0_rp/48.0_rp
          weigp(  6)= 1.0_rp/48.0_rp

          weigp(  7)= 1.0_rp/16.0_rp
          weigp(  8)= 1.0_rp/16.0_rp
          weigp(  9)= 1.0_rp/16.0_rp
          weigp( 13)= 1.0_rp/16.0_rp
          weigp( 14)= 1.0_rp/16.0_rp
          weigp( 15)= 1.0_rp/16.0_rp

          weigp( 10)= 0.0_rp
          weigp( 11)= 0.0_rp
          weigp( 12)= 0.0_rp
          weigp( 16)= 0.0_rp
          weigp( 17)= 0.0_rp
          weigp( 18)= 0.0_rp

       else
          ierro=1
       end if
    else if(ndime==2.and.ngaus==0) then
    else
       ierro=1
    end if
    !
    ! Errors
    !
    !  if(ierro==1) call runend('CLOSE_PEN: NOT AVAILABLE INTEGRATION RULE')


  end subroutine close_pen

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Closed rule for pyramid
  !> @details There is no closed quadrature rule for the pyramid as the Jacobian
  !>          is zero at the apex!
  !> 
  !-----------------------------------------------------------------------

  subroutine close_pyr(ndime,ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ndime,ngaus
    integer(ip), intent(out) :: ierro
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    integer(ip)              :: ii
    real(rp)                 :: g1
    real(rp)                 :: j,k
    real(rp)                 :: jk(2,20)

    jk(1,1)  = -1.0_rp 
    jk(2,1)  = -1.0_rp
    jk(1,2)  =  1.0_rp
    jk(2,2)  = -1.0_rp
    jk(1,3)  =  1.0_rp
    jk(2,3)  =  1.0_rp
    jk(1,4)  = -1.0_rp
    jk(2,4)  =  1.0_rp

    ierro=1
    posgp=0.0_rp
    weigp=0.0_rp

    if( ngaus == 5 ) then
       g1 = 8.0_rp*sqrt(2.0_rp/15.0_rp)/5.0_rp
       do ii = 1,4
          j = jk(1,ii)
          k = jk(2,ii)
          posgp(1,ii) = j
          posgp(2,ii) = k
          posgp(3,ii) =-2.0_rp/3.0_rp
          weigp(  ii) = 81.0_rp/100.0_rp           
       end do
       posgp(1,5) = 0.0_rp
       posgp(2,5) = 0.0_rp
       posgp(3,5) = 2.0_rp/5.0_rp
       weigp(  5) = 125.0_rp/27.0_rp  
    end if


  end subroutine close_pyr

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-07-27
  !> @brief   Change the ordering of integration points
  !> @details Change the ordering of integration points for closed integration
  !>          rules, according to : INOGA(IGAUS)=INODE
  !> 
  !-----------------------------------------------------------------------

  subroutine close_chaord(inoga,ninte)

    !-----------------------------------------------------------------------
    !
    !     Change the ordering of integration points for closed integration
    !     rules, according to : INOGA(IGAUS)=INODE
    !
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: ninte
    integer(ip), intent(out) :: inoga(*)

    if(ninte==4) then                                   ! 2 x 2 
       inoga(1)= 1
       inoga(2)= 4
       inoga(3)= 2
       inoga(4)= 3
    else if(ninte==9) then                              ! 3 x 3
       inoga(1)= 1
       inoga(2)= 8
       inoga(3)= 4
       inoga(4)= 5
       inoga(5)= 9
       inoga(6)= 7
       inoga(7)= 2
       inoga(8)= 6
       inoga(9)= 3
    else if(ninte==16) then                             ! 4 x 4
       inoga( 1)= 1
       inoga( 2)=12
       inoga( 3)=11
       inoga( 4)= 4
       inoga( 5)= 5                                     ! 4  10   9   3
       inoga( 6)=13                                     !
       inoga( 7)=16                                     ! 11 16   15  8
       inoga( 8)=10                                     ! 
       inoga( 9)= 6                                     ! 12 13   14  7
       inoga(10)=14                                     ! 
       inoga(11)=15                                     ! 1   5   6   2
       inoga(12)= 9                                        
       inoga(13)= 2
       inoga(14)= 7
       inoga(15)= 8
       inoga(16)= 3
    else if(ninte==8) then                              ! 2 x 2 x 2 
       inoga(1)= 1
       inoga(2)= 5
       inoga(3)= 4
       inoga(4)= 8
       inoga(5)= 2
       inoga(6)= 6
       inoga(7)= 3
       inoga(8)= 7
    else if(ninte==27) then                             ! 3 x 3 x 3 
       inoga( 1)= 1
       inoga( 2)=13
       inoga( 3)= 5
       inoga( 4)=12
       inoga( 5)=25
       inoga( 6)=20
       inoga( 7)= 4
       inoga( 8)=16
       inoga( 9)= 8
       inoga(10)= 9
       inoga(11)=22
       inoga(12)=17
       inoga(13)=21
       inoga(14)=27
       inoga(15)=26
       inoga(16)=11
       inoga(17)=24
       inoga(18)=19
       inoga(19)= 2
       inoga(20)=14
       inoga(21)= 6
       inoga(22)=10
       inoga(23)=23
       inoga(24)=18
       inoga(25)= 3
       inoga(26)=15
       inoga(27)= 7
    else if(ninte==64) then
       inoga( 1)= 1
       inoga( 2)=12
       inoga( 3)=21
       inoga( 4)= 5
       inoga( 5)=16
       inoga( 6)=44    
       inoga( 7)=52
       inoga( 8)=32
       inoga( 9)=15
       inoga(10)=43
       inoga(11)=51
       inoga(12)=31     
       inoga(13)= 4
       inoga(14)=20
       inoga(15)=24
       inoga(16)= 8
       inoga(17)= 9
       inoga(18)=37      
       inoga(19)=45
       inoga(20)=25
       inoga(21)=33
       inoga(22)=57
       inoga(23)=61
       inoga(24)=53     
       inoga(25)=36
       inoga(26)=60
       inoga(27)=64
       inoga(28)=56
       inoga(29)=14
       inoga(30)=42     
       inoga(31)=50
       inoga(32)=30
       inoga(33)=10
       inoga(34)=38
       inoga(35)=46
       inoga(36)=26     
       inoga(37)=34
       inoga(38)=58
       inoga(39)=62
       inoga(40)=54
       inoga(41)=35
       inoga(42)=59    
       inoga(43)=63
       inoga(44)=55
       inoga(45)=13
       inoga(46)=41
       inoga(47)=49
       inoga(48)=29     
       inoga(49)= 2
       inoga(50)=18
       inoga(51)=22
       inoga(52)= 6
       inoga(53)=11
       inoga(54)=39     
       inoga(55)=47
       inoga(56)=27
       inoga(57)=12
       inoga(58)=40
       inoga(59)=48
       inoga(60)=28     
       inoga(61)= 3
       inoga(62)=19
       inoga(63)=23
       inoga(64)= 7
    end if

  end subroutine close_chaord

end module mod_quadrature_closed
