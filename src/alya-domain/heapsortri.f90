!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine heapsortri(itask,nrows,rvin,ivou)
  !------------------------------------------------------------------------
  !****f* mathru/heapsortri
  ! NAME
  !    heapsorti2
  ! DESCRIPTION
  !    Quick sorting of an integer IVOU using a real rVIN. 
  !    The element in rvin are sorting in:
  !    ITASK = 1 ... Decreasing value, i.e., rvin(1) > rvin(2) > ...
  !    ITASK = 2 ... Increasing value, i.e., rvin(1) < rvin(2) < ...
  ! INPUT
  !    ITASK ... 1,2 for decreasing, increasing order
  !    NROWS ... Size of RVIN
  !    RVIN .... Array to be ordered (input)
  !    IVOU .... Array to be ordered (output)
  ! OUTPUT
  !    RVIN .... Ordered array
  !    IVOU .... Ordered array
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg  
  implicit none
  integer(ip), intent(in)    :: itask,nrows
  real(rp),    intent(inout) :: rvin(*) 
  integer(ip), intent(inout) :: ivou(*) 
  integer(ip)                :: leni, ir, ii, jj, jaux
  real(rp)                   :: raux

  select case(itask)

  case(1)
     !
     ! Decreasing order
     !
     if(nrows<2) then
        return
     end if

     leni = (nrows/2) + 1
     ir  = nrows

100  continue

     if (leni>1) then
        leni = leni - 1
        raux = rvin(leni)
        jaux = ivou(leni)
     else
        raux = rvin(ir)
        rvin(ir) = rvin(1)

        jaux = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           rvin(1) = raux
           ivou(1) = jaux
           return
        endif
     end if

     ii = leni
     jj = leni + leni

200  if (jj<=ir) then
        if (jj<ir) then
           if ( rvin(jj)>rvin(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (raux>rvin(jj) ) then
           rvin(ii) = rvin(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 200
     end if

     rvin(ii) = raux
     ivou(ii) = jaux

     goto 100

  case(2)
     !
     ! Increasing order
     !
     if(nrows<2) then
        return
     end if

     leni = (nrows/2) + 1
     ir  = nrows

300  continue

     if (leni>1) then
        leni = leni - 1
        raux = rvin(leni)
        jaux = ivou(leni)
     else
        raux = rvin(ir)
        rvin(ir) = rvin(1)
        jaux = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           rvin(1) = raux
           ivou(1) = jaux
           return
        endif
     end if

     ii = leni
     jj = leni + leni

400  if (jj<=ir) then
        if (jj<ir) then
           if ( rvin(jj)<rvin(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (raux<rvin(jj) ) then
           rvin(ii) = rvin(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 400
     end if

     rvin(ii) = raux
     ivou(ii) = jaux

     goto 300

  case(3)

     if(nrows<2) then
        return
     end if

     do jj=2,nrows
        raux=rvin(jj)
        jaux=ivou(jj)
        do ii=jj-1,1,-1
           if(rvin(ii)<=raux) exit
           rvin(ii+1)=rvin(ii)
           ivou(ii+1)=ivou(ii)
        end do
        rvin(ii+1)=raux
        ivou(ii+1)=jaux
     end do

  end select

end subroutine heapsortri
