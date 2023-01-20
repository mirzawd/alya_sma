!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine heapsorti2(itask,nrows,ivin,ivou)
  !------------------------------------------------------------------------
  !****f* mathru/heapsorti2
  ! NAME
  !    heapsorti2
  ! DESCRIPTION
  !    Quick sorting of IVOU using IVIN. The element in ivin are sorting in:
  !    ITASK = 1 ... Decreasing value, i.e., ivin(1) > ivin(2) > ...
  !    ITASK = 2 ... Increasing value, i.e., ivin(1) < ivin(2) < ...
  ! INPUT
  !    ITASK ... 1,2 for decreasing, increasing order
  !    NROWS ... Size of IVIN
  !    IVIN .... Array to be ordered
  !    IVOU .... Array to be ordered
  ! OUTPUT
  !    IVIN .... Ordered array
  !    IVOU .... Ordered array
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg  
  implicit none
  integer(ip), intent(in)    :: itask,nrows
  integer(ip), intent(inout) :: ivin(*) 
  integer(ip), intent(inout) :: ivou(*) 
  integer(ip)                :: leni, ir, ii, jj, iaux, jaux

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
        iaux = ivin(leni)
        jaux = ivou(leni)
     else
        iaux = ivin(ir)
        ivin(ir) = ivin(1)

        jaux = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           ivin(1) = iaux
           ivou(1) = jaux
           return
        endif
     end if

     ii = leni
     jj = leni + leni

200  if (jj<=ir) then
        if (jj<ir) then
           if ( ivin(jj)>ivin(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iaux>ivin(jj) ) then
           ivin(ii) = ivin(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 200
     end if

     ivin(ii) = iaux
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
        iaux = ivin(leni)
        jaux = ivou(leni)
     else
        iaux = ivin(ir)
        ivin(ir) = ivin(1)
        jaux = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           ivin(1) = iaux
           ivou(1) = jaux
           return
        endif
     end if

     ii = leni
     jj = leni + leni

400  if (jj<=ir) then
        if (jj<ir) then
           if ( ivin(jj)<ivin(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iaux<ivin(jj) ) then
           ivin(ii) = ivin(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 400
     end if

     ivin(ii) = iaux
     ivou(ii) = jaux

     goto 300

  case(3)

     if(nrows<2) then
        return
     end if

     do jj=2,nrows
        iaux=ivin(jj)
        jaux=ivou(jj)
        do ii=jj-1,1,-1
           if(ivin(ii)<=iaux) exit
           ivin(ii+1)=ivin(ii)
           ivou(ii+1)=ivou(ii)
        end do
        ivin(ii+1)=iaux
        ivou(ii+1)=jaux
     end do

  end select

end subroutine heapsorti2
