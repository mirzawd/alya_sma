!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain/bouele
!> @{
!> @file    bouele.f90
!> @author  houzeaux
!> @date    2022-10-06
!> @brief   Defines the boundary Gauss point from the volume
!>          Gauss point as well as the type of quadrature
!> @details Define the following variables
!>          LEXIB ... If type exists
!>          NGAUB ... Boundary Gauss points (if IELTY > 0)
!>          LQUAB ... Boundary element quadrature (if IELTY > 0)
!>    
!> @} 
!-----------------------------------------------------------------------

subroutine bouele(nelty,ngaus,lquad,ielty,ngaub,lquab,lexib)

  use def_kintyp, only : rp
  use def_elmtyp  
  implicit none
  
  integer(ip), intent(in)    :: nelty,ngaus,lquad,ielty
  integer(ip), intent(inout) :: ngaub(nelty),lquab(nelty),lexib(nelty)

  select case ( ielty )

  case ( :9   ) 
     !
     ! Point
     !
     lexib(POINT) = 1
     lquab(POINT) = lquad
     ngaub(POINT) = 1

  case( TRI03 )  
     !
     ! TRI03
     !
     lexib(BAR02) = 1
     lquab(BAR02) = lquad

     if( ngaub(BAR02) == 0 ) then

        if( lquad == 0 ) then               ! open rule

           select case ( ngaus )
           case (  :1  ) ; ngaub(BAR02) = 1                 ! exact for P1 
           case ( 2:4  ) ; ngaub(BAR02) = 2                 ! exact for P3
           case ( 5:7  ) ; ngaub(BAR02) = 3                 ! exact for P5
           case ( 8:13 ) ; ngaub(BAR02) = 4                 ! exact for P7
           end select

        else if( lquad == 1 ) then          ! closed rule

           select case ( ngaus )
           case (  :4  ) ; ngaub(BAR02) = 2
           case ( 5:7  ) ; ngaub(BAR02) = 3
           case ( 8:10 ) ; ngaub(BAR02) = 4
           end select

        end if
     end if

  case ( TRI06 )    
     !
     ! TRI06
     !
     lexib(BAR03) = 1
     lquab(BAR03) = lquad

     if( ngaub(BAR03) == 0 ) then 

        if( lquad == 0 ) then               ! open rule

           select case ( ngaus ) 
           case (  :1 ) ; ngaub(BAR03) = 1                 ! exact for P1 
           case ( 2:4 ) ; ngaub(BAR03) = 2                 ! exact for P3
           case ( 5:7 ) ; ngaub(BAR03) = 3                 ! exact for P5
           case (8:13 ) ; ngaub(BAR03) = 4                 ! exact for P7
           end select

        else if( lquad == 1 ) then          ! closed rule

           select case ( ngaus ) 
           case (  :4  ) ; ngaub(BAR03) = 2
           case ( 5:7  ) ; ngaub(BAR03) = 3
           case ( 8:10 ) ; ngaub(BAR03) = 4
           end select

        end if
     end if

  case ( TRI10 )    
     !
     ! TRI10
     !
     lexib(BAR04) = 1
     lquab(BAR04) = lquad

     if( ngaub(BAR04) == 0 ) then

        if( lquad == 0 ) then               ! open rule

           select case ( ngaus )
           case (  :1  ) ; ngaub(BAR04) = 1                 ! exact for P1 
           case ( 2:4  ) ; ngaub(BAR04) = 2                 ! exact for P3
           case ( 5:7  ) ; ngaub(BAR04) = 3                 ! exact for P5
           case ( 8:13 ) ; ngaub(BAR04) = 4                 ! exact for P7
           end select

        else if( lquad == 1 ) then          ! closed rule

           select case ( ngaus )
           case ( :4   ) ; ngaub(BAR04) = 2
           case ( 5:7  ) ; ngaub(BAR04) = 3
           case ( 8:10 ) ; ngaub(BAR04) = 4
           end select

        end if
     end if

  case ( QUA04 )    
     !
     ! QUA04
     !
     lexib(BAR02) = 1
     lquab(BAR02) = lquad
     if( ngaub(BAR02) == 0 ) ngaub(BAR02) = int(sqrt(real(ngaus,rp)))

  case ( QUA08 )    
     !
     ! QUA08
     !
     lexib(BAR03) = 1
     lquab(BAR03) = lquad
     if( ngaub(BAR03) == 0 ) ngaub(BAR03) = int(sqrt(real(ngaus,rp)))

  case ( QUA09 )    
     !
     ! QUA09
     !
     lexib(BAR03) = 1
     lquab(BAR03) = lquad
     if( ngaub(BAR03) == 0 ) ngaub(BAR03) = int(sqrt(real(ngaus,rp)))

  case ( QUA16 )   
     !
     ! QUA16
     !
     lexib(BAR04) = 1
     lquab(BAR04) = lquad
     if( ngaub(BAR04) == 0 ) ngaub(BAR04) = int(sqrt(real(ngaus,rp)))

  case ( TET04 )  
     !
     ! TET04
     !
     lexib(TRI03) = 1
     lquab(TRI03) = lquad

     if( ngaub(TRI03) == 0 ) then
        if( lquad == 0 ) then                     ! open rule

           select case ( ngaus )
           case ( 1   ) ; ngaub(TRI03) = 1        ! exact for P1
           case ( 2:5 ) ; ngaub(TRI03) = ngaus-1  ! exact for P2/P3
           case ( 8   ) ; ngaub(TRI03) = 6        ! exact for P4
           case ( 11  ) ; ngaub(TRI03) = 6        ! exact for P4
           case ( 14  ) ; ngaub(TRI03) = 7        ! exact for P5
           case ( 15  ) ; ngaub(TRI03) = 7        ! exact for P5
           case ( 29  ) ; ngaub(TRI03) = 13       ! 
           end select

        else if( lquad == 1 ) then                ! closed rule

           select case ( ngaus )
           case ( 1:5  ) ; ngaub(TRI03) = 3
           case ( 6:11 ) ; ngaub(TRI03) = 6
           case ( 15   ) ; ngaub(TRI03) = 7
           case ( 20   ) ; ngaub(TRI03) = 10
           end select

        end if
     end if

  case ( TET10 ) 
     !
     ! TET10
     !
     lexib(TRI06) = 1
     lquab(TRI06) = lquad

     if( ngaub(TRI06) == 0 ) then

        if( lquad == 0 ) then                  ! open rule

           select case ( ngaus )
           case ( 1   ) ; ngaub(TRI06) = 1                     ! exact for P1
           case ( 2:5 ) ; ngaub(TRI06) = ngaus-1               ! exact for P2/P3
           case ( 8   ) ; ngaub(TRI06) = 6                     ! exact for P4
           case ( 11  ) ; ngaub(TRI06) = 6                     ! exact for P4
           case ( 14  ) ; ngaub(TRI06) = 7                     ! exact for P5
           case ( 15  ) ; ngaub(TRI06) = 7                     ! exact for P5
           case ( 45  ) ; ngaub(TRI06) = 7                     ! exact for P5
           case ( 65  ) ; ngaub(TRI06) = 7                     ! exact for P5
           end select

        else if( lquad == 1 ) then             ! closed rule

           select case ( ngaus )               
           case ( :5   ) ; ngaub(TRI06) = 3
           case ( 6:11 ) ; ngaub(TRI06) = 6
           case ( 15   ) ; ngaub(TRI06) = 7
           case ( 20   ) ; ngaub(TRI06) = 10
           end select

        end if
     end if

  case ( TET20 ) 
     !
     ! TET20
     !
     lexib(TRI10) = 1
     lquab(TRI10) = lquad

     if( ngaub(TRI10) == 0 ) then
        
        if( lquad == 0 ) then                  ! open rule

           select case ( ngaus ) 
           case ( 1   ) ; ngaub(TRI10) = 1                     ! exact for P1
           case ( 2:5 ) ; ngaub(TRI10) = ngaus-1               ! exact for P2/P3
           case ( 8   ) ; ngaub(TRI10) = 6                     ! exact for P4
           case ( 11  ) ; ngaub(TRI10) = 6                     ! exact for P4
           case ( 14  ) ; ngaub(TRI10) = 7                     ! exact for P5
           case ( 15  ) ; ngaub(TRI10) = 7                     ! exact for P5
           case ( 29  ) ; ngaub(TRI10) = 13                    ! exact for P6 (Must use this one!!!!)
           end select
           
        else if( lquad == 1 ) then              ! closed rule (don't use!!!!)

           select case ( ngaus ) 
           
           case ( 1:5  ) ; ngaub(TRI10) = 3
           case ( 6:11 ) ; ngaub(TRI10) = 6
           case ( 15   ) ; ngaub(TRI10) = 7
           case ( 20   ) ; ngaub(TRI10) = 10
           end select
           
        end if
     end if

  case ( PYR05 ) 
     !
     ! PYR05
     !
     lexib(TRI03) = 1
     lexib(QUA04) = 1

     select case ( ngaus )
     case ( 1 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 1
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 1
     case ( 5,6,8,9 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 3
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 4        
     case ( 13,18 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 6
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 9                  
     case default
        call runend('BOUELE: WRONG NUMBER OF GAUSS POINTS FOR PYRA_5 ')
     end select

  case ( PYR14 )   
     !
     ! PYR14
     !
     lexib(TRI06) = 1
     lexib(QUA08) = 1
     call runend('BOUELE: PYRA_14 ELEMENT IS NOT READY')

  case ( PEN06 )   
     !
     ! PEN06
     !
     lexib(TRI03) = 1
     lexib(QUA04) = 1
     lquab(TRI03) = lquad
     lquab(QUA04) = lquad

     select case ( ngaus ) 
     case ( 1  )
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 1
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 1
     case ( 6  ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 3
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 4
     case ( 8  ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 4
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 4
     case ( 11 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 6
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 9
     case ( 16 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 7
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 9
     case ( 24 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 13
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 16
     case ( 29 ) 
        if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 13
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 16
     case default
        call runend('PEN06: NOT CODED')
     end select

  case ( PEN15 ) 
     !
     ! PEN15
     !
     lexib(TRI06 ) = 1
     lexib(QUA08) = 1

     select case ( ngaus )
     case ( 1 ) 
        if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 1
        if( ngaub(QUA08) == 0 ) ngaub(QUA08) = 1
     case ( 6 ) 
        if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 3
        if( ngaub(QUA08) == 0 ) ngaub(QUA08) = 4
     end select

     call runend('BOUELE: PENTA_15 ELEMENT IS NOT READY')

  case ( PEN18 ) 
     !
     ! PEN18
     !
     lexib(TRI06) = 1
     lexib(QUA09) = 1

     select case ( ngaus )
     case ( 1  ) 
        if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 1
        if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 1
     case ( 6  ) 
        if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 3
        if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 4
     case ( 18 ) 
        if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 6
        if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 9
     case ( 21 ) 
        if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7
        if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 9
     end select

     ! call runend('BOUELE: PENTA_18 ELEMENT IS NOT READY')

  case ( HEX08 )    
     !
     ! HEX08
     !
     lexib(QUA04) = 1
     lquab(QUA04) = lquad
     if( ngaub(QUA04) == 0 ) ngaub(QUA04) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))

     !case ( HEX20 ) then  

     !   lexib(QUA08) = 1
     !   lquab(QUA08) = lquad
     !   if(ngaus==20) then
     !      if( ngaub(QUA08) == 0 ) ngaub(QUA08) = 9
     !   else
     !      if( ngaub(QUA08) == 0 ) ngaub(QUA08) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))
     !   end if

  case ( HEX27 )  
     !
     ! HEX27
     !
     lexib(QUA09) = 1
     lquab(QUA09) = lquad
     if( ngaub(QUA09) == 0 ) ngaub(QUA09) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))

  case ( HEX64 )  
     !
     ! HEX64
     !
     lexib(QUA16) = 1
     lquab(QUA16) = lquad
     if( ngaub(QUA16) == 0 ) ngaub(QUA16) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))

  case ( SHELL )    
     !
     ! SHELL
     !
     lexib(BAR02) = 1
     lquab(BAR02) = lquad
     if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 1

  case ( BAR3D )    
     !
     ! BAR3D
     !
     lexib(POINT) = 1
     lquab(POINT) = lquad
     ngaub(POINT) = 1

  case default

     call runend('BOUELE: UNDEFINED ELEMENT')

  end select

end subroutine bouele

