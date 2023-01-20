!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine hsort2(itask,nrows,ivi1,ivi2,ivou)
  !------------------------------------------------------------------------
  !****f* mathru/hsort2
  ! NAME
  !    hsort2
  ! DESCRIPTION
  !    Quick sorting of IVOU using IVI1. The element in ivi1 are sorting in:
  !    ITASK = 1 ... Decreasing value, i.e., ivi1(1) > ivi1(2) > ...
  !    ITASK = 2 ... Increasing value, i.e., ivi1(1) < ivi1(2) < ...
  ! INPUT
  !    ITASK ... 1,2 for decreasing, increasing order
  !    NROWS ... Size of IVI1
  !    IVI1 .... Array to be ordered
  !    IVOU .... Array to be ordered
  ! OUTPUT
  !    IVI1 .... Ordered array
  !    IVOU .... Ordered array
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg  
  implicit none
  integer(ip), intent(in)    :: itask,nrows
  integer(ip), intent(inout) :: ivi1(*) 
  integer(ip), intent(inout) :: ivi2(*) 
  integer(ip), intent(inout) :: ivou(*) 
  integer(ip)                :: leni, ir, ii, jj, iau1, iau2, jaux

  select case(itask)

  case(1_ip)

     !-------------------------------------------------------------------
     !
     ! Decreasing order
     !
     !-------------------------------------------------------------------
     
     if(nrows<2) then
        return
     end if

     leni = (nrows/2) + 1
     ir  = nrows

100  continue

     if (leni>1) then
        leni = leni - 1
        iau1 = ivi1(leni)
        iau2 = ivi2(leni)
        jaux = ivou(leni)
     else
        iau1     = ivi1(ir)
        ivi1(ir) = ivi1(1)
        iau2     = ivi2(ir)
        ivi2(ir) = ivi2(1)

        jaux     = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           ivi1(1) = iau1
           ivi2(1) = iau2
           ivou(1) = jaux
           return
        endif
     end if

     ii = leni
     jj = leni + leni

200  if (jj<=ir) then
        if (jj<ir) then
           if ( ivi1(jj)>ivi1(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iau1>ivi1(jj) ) then
           ivi1(ii) = ivi1(jj)
           ivi2(ii) = ivi2(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 200
     end if

     ivi1(ii) = iau1
     ivi2(ii) = iau2
     ivou(ii) = jaux

     goto 100

  case (2_ip)

     !-------------------------------------------------------------------
     !
     ! Increasing order
     !
     !-------------------------------------------------------------------

     if(nrows<2) then
        return
     end if

     leni = (nrows/2) + 1
     ir  = nrows

300  continue

     if ( leni > 1 ) then
        leni = leni - 1
        iau1 = ivi1(leni)
        iau2 = ivi2(leni)
        jaux = ivou(leni)
     else
        iau1     = ivi1(ir)
        ivi1(ir) = ivi1(1)
        iau2     = ivi2(ir)
        ivi2(ir) = ivi2(1)

        jaux     = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           ivi1(1) = iau1
           ivi2(1) = iau2
           ivou(1) = jaux
           goto 301
        endif
     end if

     ii = leni
     jj = leni + leni

400  if (jj<=ir) then
        if (jj<ir) then
           if ( ivi1(jj)<ivi1(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iau1<ivi1(jj) ) then
           ivi1(ii) = ivi1(jj)
           ivi2(ii) = ivi2(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 400
     end if

     ivi1(ii) = iau1
     ivi2(ii) = iau2
     ivou(ii) = jaux

     goto 300

301  continue

     do ii = 1,nrows-1
        do jj = ii+1,nrows
           if( ivi1(ii) == ivi1(jj) ) then
              if( ivi2(ii) > ivi2(jj) ) then
                 iau1     = ivi1(ii)
                 ivi1(ii) = ivi1(jj) 
                 ivi1(jj) = iau1
                 iau2     = ivi2(ii)
                 ivi2(ii) = ivi2(jj) 
                 ivi2(jj) = iau2
                 jaux     = ivou(ii)
                 ivou(ii) = ivou(jj) 
                 ivou(jj) = jaux
              end if
           end if
        end do
     end do
 

  case(3)

     if(nrows<2) then
        return
     end if

     do jj=2,nrows
        iau1=ivi1(jj)
        iau2=ivi2(jj)
        jaux=ivou(jj)
        do ii=jj-1,1,-1
           if(ivi1(ii)<=iau1) exit
           ivi1(ii+1)=ivi1(ii)
           ivi2(ii+1)=ivi2(ii)
           ivou(ii+1)=ivou(ii)
        end do
        ivi1(ii+1)=iau1
        ivi2(ii+1)=iau2
        ivou(ii+1)=jaux
     end do

  end select

end subroutine hsort2
