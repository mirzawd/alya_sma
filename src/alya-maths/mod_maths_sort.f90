!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical solvers
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for algebraic solvers
!> @details ToolBox for algebraic solvers
!
!-----------------------------------------------------------------------

module mod_maths_sort

  use def_maths
  use mod_maths_basic
  implicit none

  private
 
  interface maths_heap_sort
     module procedure maths_heap_sort_I1,&
          &           maths_heap_sort_I2,&
          &           maths_heap_sort_IP_2,&
          &           maths_heap_sort_RP1,&
          &           maths_heap_sort_RP2
  end interface maths_heap_sort
  interface maths_quick_sort
     module procedure maths_quick_sort_rp,&
          &           maths_quick_sort_rp_1,&
          &           maths_quick_sort_ip,&
          &           maths_quick_sort_ip_1
  end interface maths_quick_sort

  public :: maths_heap_sort
  public :: maths_quick_sort
  public :: maths_geometrical_sort_using_coordinates
  public :: maths_geometrical_sort_using_sfc

contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   heap sort
  !> @details Equalize some arrays
  !>          Quick sorting. The element in ivin are sorting in:
  !>          ITASK = 1 ... Decreasing value, i.e., ivin(1) > ivin(2) > ...
  !>          ITASK = 2 ... Increasing value, i.e., ivin(1) < ivin(2) < ...
  !
  !----------------------------------------------------------------------

  subroutine maths_heap_sort_I1(itask,nrows,ivin,message,ivo1,ivo2,PERMUTATION,ELIMINATE_REPLICATES)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(inout)           :: nrows
    integer(ip),  intent(inout)           :: ivin(*)
    character(*), intent(in),    optional :: message
    integer(ip),  intent(inout), optional :: ivo1(*)
    integer(ip),  intent(inout), optional :: ivo2(*)
    integer(ip),  intent(inout), optional :: PERMUTATION(*)
    logical(lg),  intent(in),    optional :: ELIMINATE_REPLICATES
    integer(ip)                           :: leni,ir,ii,jj,iaux
    integer(ip)                           :: iau1,iau2,kk
    logical(lg)                           :: if_eliminate_replicates

    if_eliminate_replicates = .false.
    if( present(ELIMINATE_REPLICATES) ) if_eliminate_replicates = ELIMINATE_REPLICATES
    if( present(message) ) then
       if( trim(message) == 'ELIMINATE DUPLICATES') if_eliminate_replicates = .true.
    end if

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

100    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
       else
          iaux = ivin(ir)
          if( present(ivo1) ) iau1 = ivo1(ir)
          if( present(ivo2) ) iau2 = ivo2(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) ivo1(ir) = ivo1(1)
          if( present(ivo2) ) ivo2(ir) = ivo2(1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) then
             iau1     = ivo1(ir)
             ivo1(ir) = ivo1(1)
          end if
          if( present(ivo2) ) then
             iau2     = ivo2(ir)
             ivo2(ir) = ivo2(1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

       if( present(PERMUTATION) ) then
401       if( jj <= ir ) then
             if( jj < ir ) then
                if ( PERMUTATION(ivin(jj)) < PERMUTATION(ivin(jj+1)) ) then
                   jj = jj + 1
                end if
             end if

             if( PERMUTATION(iaux) < PERMUTATION(ivin(jj)) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(ii) = ivo1(jj)
                if( present(ivo2) ) ivo2(ii) = ivo2(jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 401
          end if
       else
400       if( jj <= ir ) then
             if( jj < ir ) then
                if ( ivin(jj) < ivin(jj+1) ) then
                   jj = jj + 1
                end if
             end if

             if( iaux < ivin(jj) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(ii) = ivo1(jj)
                if( present(ivo2) ) ivo2(ii) = ivo2(jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 400
          end if
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 300

    end select
    !
    ! Eliminate duplicates
    !
500 if( if_eliminate_replicates ) then
       if( itask == 1 ) stop
       jj = 0
       kk = 0
       do ii = 1,nrows
          if( ivin(ii) <= jj ) then
             ivin(ii) = 0
             if( present(ivo1) ) ivo1(ii) = 0
             if( present(ivo2) ) ivo2(ii) = 0
          else
             kk       = kk + 1
             ivin(kk) = ivin(ii)
             jj       = ivin(ii)
             if( present(ivo1) ) ivo1(kk) = ivo1(ii)
             if( present(ivo2) ) ivo2(kk) = ivo2(ii)
          end if
       end do
       if( kk < nrows ) then
          ivin(kk+1:nrows) = 0
          if( present(ivo1) ) ivo1(kk+1:nrows) = 0
          if( present(ivo2) ) ivo2(kk+1:nrows) = 0
       end if
       nrows = kk
    end if

  end subroutine maths_heap_sort_I1

  pure subroutine maths_heap_sort_I2(itask,nrows,ivin,ivo1,ivo2,message)

    integer(ip),  intent(in)                        :: itask
    integer(ip),  intent(inout)                     :: nrows
    integer(ip),  intent(inout)                     :: ivin(*)
    integer(ip),  intent(inout),  pointer           :: ivo1(:,:)
    integer(ip),  intent(inout),  pointer, optional :: ivo2(:,:)
    character(*), intent(in),              optional :: message
    integer(ip)                                     :: leni,ir,ii,jj,iaux,ncol1,ncol2
    integer(ip),                  pointer           :: iau1(:)
    integer(ip),                  pointer           :: iau2(:)

    nullify(iau1)
    nullify(iau2)

    ncol1 = size(ivo1,2,KIND=ip)
    if( ncol1 < 1 ) return
    !if( size(ivo1,1,KIND=ip) < nrows ) stop
    allocate(iau1(ncol1))

    if( present(ivo2) ) then
       ncol2 = size(ivo2,2,KIND=ip)
       if( ncol2 < 1 ) return
       !if( size(ivo2,1,KIND=ip) < nrows ) stop
       allocate(iau2(ncol2))
    end if

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

100    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          iau1(1:ncol1) = ivo1(leni,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2) = ivo2(leni,1:ncol2)
          end if
       else
          iaux = ivin(ir)
          iau1(1:ncol1) = ivo1(ir,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2) = ivo2(ir,1:ncol2)
          end if
          ivin(ir) = ivin(1)
          ivo1(ir,1:ncol1) = ivo1(1,1:ncol1)
          if( present(ivo2) ) then
             ivo2(ir,1:ncol2) = ivo2(1,1:ncol2)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivo1(1,1:ncol1) = iau1(1:ncol1)
             if( present(ivo2) ) then
                ivo2(1,1:ncol2) = iau2(1:ncol2)
             end if
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivo1(ii,1:ncol1) = ivo1(jj,1:ncol1)
             if( present(ivo2) ) then
                ivo2(ii,1:ncol2) = ivo2(jj,1:ncol2)
             end if
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       ivo1(ii,1:ncol1) = iau1(1:ncol1)
       if( present(ivo2) ) then
          ivo2(ii,1:ncol2) = iau2(1:ncol2)
       end if

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          iau1(1:ncol1) = ivo1(leni,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2) = ivo2(leni,1:ncol2)
          end if
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          iau1(1:ncol1)    = ivo1(ir,1:ncol1)
          ivo1(ir,1:ncol1) = ivo1(1,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2)    = ivo2(ir,1:ncol2)
             ivo2(ir,1:ncol2) = ivo2(1,1:ncol2)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivo1(1,1:ncol1) = iau1(1:ncol1)
             if( present(ivo2) ) then
                ivo2(1,1:ncol2) = iau2(1:ncol2)
             end if
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivo1(ii,1:ncol1) = ivo1(jj,1:ncol1)
             if( present(ivo2) ) then
                ivo2(ii,1:ncol2) = ivo2(jj,1:ncol2)
             end if
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       ivin(ii) = iaux
       ivo1(ii,1:ncol1) = iau1(1:ncol1)
       if( present(ivo2) ) then
          ivo2(ii,1:ncol2) = iau2(1:ncol2)
       end if

       goto 300

    end select

    if( associated(iau1) ) deallocate(iau1)
    if( associated(iau2) ) deallocate(iau2)

    return
    !
    ! Eliminate duplicates
    !
500 continue
    if( present(message) ) then
       if( trim(message) == 'ELIMINATE DUPLICATES') then
          !stop
       end if
    end if

  end subroutine maths_heap_sort_I2

  pure subroutine maths_heap_sort_IP_2(itask,ivou,NROWS)

    integer(ip),                     intent(in)    :: itask
    integer(ip),            pointer, intent(inout) :: ivou(:,:)
    integer(ip),  optional,          intent(in)    :: NROWS
    integer(ip)                                    :: leni,ir,ii,jj,iaux
    integer(ip)                                    :: kk,kdime,nn
    integer(ip),  pointer                          :: ivin(:)
    integer(ip),  pointer                          :: iau1(:)
    
    nullify(ivin,iau1)

    if( associated(ivou) ) then

       kdime = size(ivou,1_ip)
       if( present(NROWS) ) then
          nn = NROWS
       else
          nn = size(ivou,2_ip)
       end if
       allocate(ivin(nn))
       allocate(iau1(kdime))

       do kk = 1,kdime

          do ii = 1,nn
             ivin(ii) = ivou(kk,ii)
          end do

          select case ( itask )

          case ( 1_ip )
             !
             ! Decreasing order
             !
             if( nn < 2 ) then
                goto 500
             end if

             leni = (nn/2) + 1
             ir  = nn

100          continue

             if( leni > 1 ) then
                leni = leni - 1
                iaux = ivin(leni)
                iau1 = ivou(:,leni)
             else
                iaux = ivin(ir)
                iau1 = ivou(:,ir)
                ivin(ir) = ivin(1)
                ivou(:,ir) = ivou(:,1)

                ir = ir - 1

                if( ir == 1 ) then
                   ivin(1) = iaux
                   ivou(:,1) = iau1
                   goto 500
                end if
             end if

             ii = leni
             jj = leni + leni

200          if( jj <= ir ) then
                if( jj < ir ) then
                   if( ivin(jj) > ivin(jj+1) ) then
                      jj = jj + 1
                   end if
                end if

                if( iaux > ivin(jj) ) then
                   ivin(ii) = ivin(jj)
                   ivou(:,ii) = ivou(:,jj)
                   ii = jj
                   jj = jj + jj
                else
                   jj = ir + 1
                endif

                goto 200
             end if

             ivin(ii) = iaux
             ivou(:,ii) = iau1

             goto 100

          case ( 2_ip )
             !
             ! Increasing order
             !
             if( nn < 2 ) then
                goto 500
             end if

             leni = (nn/2) + 1
             ir  = nn

300          continue

             if( leni > 1 ) then
                leni       = leni - 1
                iaux       = ivin(leni)
                iau1       = ivou(:,leni)
             else
                iaux       = ivin(ir)
                ivin(ir)   = ivin(1)
                iau1       = ivou(:,ir)
                ivou(:,ir) = ivou(:,1)
                ir         = ir - 1

                if( ir == 1 ) then
                   ivin(1)    = iaux
                   ivou(:,1)  = iau1
                   goto 500
                end if
             end if

             ii = leni
             jj = leni + leni

400          if( jj <= ir ) then
                if( jj < ir ) then
                   if ( ivin(jj) < ivin(jj+1) ) then
                      jj = jj + 1
                   end if
                end if

                if( iaux < ivin(jj) ) then
                   ivin(ii) = ivin(jj)
                   ivou(:,ii) = ivou(:,jj)
                   ii = jj
                   jj = jj + jj
                else
                   jj = ir + 1
                end if

                goto 400
             end if

             ivin(ii)   = iaux
             ivou(:,ii) = iau1

             goto 300


          end select

500       continue

       end do

       deallocate(iau1)
       deallocate(ivin)

    end if

  end subroutine maths_heap_sort_IP_2

  pure subroutine maths_heap_sort_RP1(itask,nrows,rvin,ivo1,ivo2,rvo1,rvo2)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(in)              :: nrows
    real(rp),     intent(inout)           :: rvin(*)
    integer(ip),  intent(inout), optional :: ivo1(*)
    integer(ip),  intent(inout), optional :: ivo2(*)
    real(rp),     intent(inout), optional :: rvo1(*)
    real(rp),     intent(inout), optional :: rvo2(*)
    real(rp)                              :: iaux
    integer(ip)                           :: leni,ir,ii,jj
    integer(ip)                           :: iau1,iau2
    real(rp)                              :: rau1,rau2
    
    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

100    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = rvin(leni)
          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
          if( present(rvo1) ) rau1 = rvo1(leni)
          if( present(rvo2) ) rau2 = rvo2(leni)
       else
          iaux = rvin(ir)
          if( present(ivo1) ) iau1 = ivo1(ir)
          if( present(ivo2) ) iau2 = ivo2(ir)
          if( present(rvo1) ) rau1 = rvo1(ir)
          if( present(rvo2) ) rau2 = rvo2(ir)
          rvin(ir) = rvin(1)
          if( present(ivo1) ) ivo1(ir) = ivo1(1)
          if( present(ivo2) ) ivo2(ir) = ivo2(1)
          if( present(rvo1) ) rvo1(ir) = rvo1(1)
          if( present(rvo2) ) rvo2(ir) = rvo2(1)

          ir = ir - 1

          if( ir == 1 ) then
             rvin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             if( present(rvo1) ) rvo1(1) = rau1
             if( present(rvo2) ) rvo2(1) = rau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( rvin(jj) > rvin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > rvin(jj) ) then
             rvin(ii) = rvin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             if( present(rvo1) ) rvo1(ii) = rvo1(jj)
             if( present(rvo2) ) rvo2(ii) = rvo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       rvin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2
       if( present(rvo1) ) rvo1(ii) = rau1
       if( present(rvo2) ) rvo2(ii) = rau2

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = rvin(leni)

          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
          if( present(rvo1) ) rau1 = rvo1(leni)
          if( present(rvo2) ) rau2 = rvo2(leni)
       else
          iaux     = rvin(ir)
          rvin(ir) = rvin(1)
          if( present(ivo1) ) then
             iau1     = ivo1(ir)
             ivo1(ir) = ivo1(1)
          end if
          if( present(ivo2) ) then
             iau2     = ivo2(ir)
             ivo2(ir) = ivo2(1)
          end if
          if( present(rvo1) ) then
             rau1     = rvo1(ir)
             rvo1(ir) = rvo1(1)
          end if
          if( present(rvo2) ) then
             rau2     = rvo2(ir)
             rvo2(ir) = rvo2(1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             rvin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             if( present(rvo1) ) rvo1(1) = rau1
             if( present(rvo2) ) rvo2(1) = rau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( rvin(jj) < rvin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < rvin(jj) ) then
             rvin(ii) = rvin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             if( present(rvo1) ) rvo1(ii) = rvo1(jj)
             if( present(rvo2) ) rvo2(ii) = rvo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       rvin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2
       if( present(rvo1) ) rvo1(ii) = rau1
       if( present(rvo2) ) rvo2(ii) = rau2

       goto 300

    end select

500 continue

  end subroutine maths_heap_sort_RP1

  pure subroutine maths_heap_sort_RP2(itask,ndofn,nrows,ivin,ivo1,ivo2,PERMUTATION)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(in)              :: ndofn
    integer(ip),  intent(inout)           :: nrows
    integer(ip),  intent(inout)           :: ivin(*)
    real(rp),     intent(inout), optional :: ivo1(ndofn,*)
    real(rp),     intent(inout), optional :: ivo2(ndofn,*)
    integer(ip),  intent(inout), optional :: PERMUTATION(*)
    integer(ip)                           :: iaux
    integer(ip)                           :: leni,ir,ii,jj
    real(rp)                              :: iau1(ndofn),iau2(ndofn)

    select case ( itask )

    case ( 1_ip )
       !
       ! Oups
       !
       return

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)

          if( present(ivo1) ) iau1 = ivo1(:,leni)
          if( present(ivo2) ) iau2 = ivo2(:,leni)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) then
             iau1       = ivo1(:,ir)
             ivo1(:,ir) = ivo1(:,1)
          end if
          if( present(ivo2) ) then
             iau2       = ivo2(:,ir)
             ivo2(:,ir) = ivo2(:,1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(:,1) = iau1
             if( present(ivo2) ) ivo2(:,1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

       if( present(PERMUTATION) ) then
401       if( jj <= ir ) then
             if( jj < ir ) then
                if ( PERMUTATION(ivin(jj)) < PERMUTATION(ivin(jj+1)) ) then
                   jj = jj + 1
                end if
             end if

             if( PERMUTATION(iaux) < PERMUTATION(ivin(jj)) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(:,ii) = ivo1(:,jj)
                if( present(ivo2) ) ivo2(:,ii) = ivo2(:,jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 401
          end if
       else
400       if( jj <= ir ) then
             if( jj < ir ) then
                if ( ivin(jj) < ivin(jj+1) ) then
                   jj = jj + 1
                end if
             end if

             if( iaux < ivin(jj) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(:,ii) = ivo1(:,jj)
                if( present(ivo2) ) ivo2(:,ii) = ivo2(:,jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 400
          end if
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(:,ii) = iau1
       if( present(ivo2) ) ivo2(:,ii) = iau2

       goto 300

    end select

500 continue

  end subroutine maths_heap_sort_RP2

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   heap sort
  !> @details Quick sort
  !>          Quick sorting. The element in ivin are sorting in:
  !>          ITASK = 1 ... Decreasing value, i.e., ivin(1) > ivin(2) > ...
  !>          ITASK = 2 ... Increasing value, i.e., ivin(1) < ivin(2) < ...
  !
  !----------------------------------------------------------------------

  recursive subroutine maths_quick_sort_rp_1(itask,ivin)

    integer(ip),          intent(in)    :: itask
    real(rp),    pointer, intent(inout) :: ivin(:)
    integer(ip)                         :: iq,nrows

    nrows = size(ivin)
    if( nrows > 1) then
       call maths_partition_rp(itask,nrows,ivin,iq)
       call maths_quick_sort_rp(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_rp(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_rp_1

 recursive subroutine maths_quick_sort_rp(itask,nrows,ivin)

    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    real(rp),    intent(inout) :: ivin(*)
    integer(ip)                :: iq

    if( nrows > 1) then
       call maths_partition_rp(itask,nrows,ivin,iq)
       call maths_quick_sort_rp(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_rp(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_rp

  recursive subroutine maths_quick_sort_ip_1(itask,ivin)

    integer(ip),          intent(in)    :: itask
    integer(ip), pointer, intent(inout) :: ivin(:)
    integer(ip)                         :: iq,nrows

    nrows = size(ivin)
    if( nrows > 1) then
       call maths_partition_ip(itask,nrows,ivin,iq)
       call maths_quick_sort_ip(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_ip(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_ip_1

 recursive subroutine maths_quick_sort_ip(itask,nrows,ivin)

    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    integer(ip), intent(inout) :: ivin(*)
    integer(ip)                :: iq

    if( nrows > 1) then
       call maths_partition_ip(itask,nrows,ivin,iq)
       call maths_quick_sort_ip(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_ip(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_ip

  subroutine maths_partition_rp(itask,nrows,ivin,marker)
    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    real(rp),    intent(inout) :: ivin(*)
    integer(ip), intent(out)   :: marker
    integer(ip)                :: i,j
    real(rp)                   :: temp
    real(rp)                   :: x      ! pivot point

    x = ivin(1)
    i= 0
    j= nrows + 1

    if( itask == 1 ) then
       do
          j = j-1
          do
             if( ivin(j) >= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) <= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    else
       do
          j = j-1
          do
             if( ivin(j) <= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) >= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    end if

  end subroutine maths_partition_rp

  subroutine maths_partition_ip(itask,nrows,ivin,marker)
    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    integer(ip), intent(inout) :: ivin(*)
    integer(ip), intent(out)   :: marker
    integer(ip)                :: i,j
    integer(ip)                :: temp
    integer(ip)                :: x      ! pivot point

    x = ivin(1)
    i= 0
    j= nrows + 1

    if( itask == 1 ) then
       do
          j = j-1
          do
             if( ivin(j) >= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) <= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    else
       do
          j = j-1
          do
             if( ivin(j) <= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) >= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    end if

  end subroutine maths_partition_ip
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux and gguillamet
  !> @date    2019-09-02
  !> @brief   Geometrical ordering
  !> @details Uses a geometrical ordering based on coordinates
  !>          It is not required to use stable sort method.
  !>
  !>          \verbatim
  !>          The following example shows how the algorithm sort a list
  !>          of nodes using same coordinates:
  !>
  !>          Unsorted list:
  !>          node      x     y     z
  !>             1   -3.0  -2.0   3.0
  !>             2   -3.0  -2.0   4.0
  !>             3    0.0  -2.0  -7.0
  !>             4    0.0  -2.0   2.0
  !>             5   -1.0  -2.0  -1.0
  !>             6   -1.0   1.0   2.0
  !>             7   -3.0   2.0  12.0
  !>             8    0.0  -2.0   8.0
  !>             9    0.0  -2.0   9.0
  !>
  !>          Sorted list:
  !>          node      x     y     z
  !>             1   -3.0  -2.0   3.0
  !>             2   -3.0  -2.0   4.0
  !>             7   -3.0   2.0  12.0
  !>             5   -1.0  -2.0  -1.0
  !>             6   -1.0   1.0   2.0
  !>             3    0.0  -2.0  -7.0
  !>             4    0.0  -2.0   2.0
  !>             8    0.0  -2.0   8.0
  !>             9    0.0  -2.0   9.0
  !>          \endverbatim
  !>
  !>          Unit test:
  !>
  !>          integer(ip)          :: lsize,ii,nn,nz,iboun
  !>          integer(ip), pointer :: kpoin(:)
  !>          real(rp),    pointer :: xcoor(:,:)
  !>          integer(ip), pointer :: ia(:)
  !>          integer(ip), pointer :: ja(:)
  !>          real(rp),    pointer :: aa(:)
  !>          real(rp),    pointer :: diag(:)
  !>          lsize = 9
  !>          allocate(kpoin(lsize))
  !>          allocate(xcoor(ndime,lsize))
  !>          do ii = 1,lsize
  !>             kpoin(ii) = ii
  !>          end do
  !>          xcoor(1:3,1) = (/ -3.0_rp, -2.0_rp,  3.0_rp /)
  !>          xcoor(1:3,2) = (/ -3.0_rp, -2.0_rp,  4.0_rp /)
  !>          xcoor(1:3,3) = (/  0.0_rp, -2.0_rp, -7.0_rp /)
  !>          xcoor(1:3,4) = (/  0.0_rp, -2.0_rp,  2.0_rp /)
  !>          xcoor(1:3,5) = (/ -1.0_rp, -2.0_rp, -1.0_rp /)
  !>          xcoor(1:3,6) = (/ -1.0_rp,  1.0_rp,  2.0_rp /)
  !>          xcoor(1:3,7) = (/ -3.0_rp,  2.0_rp, 12.0_rp /)
  !>          xcoor(1:3,8) = (/  0.0_rp, -2.0_rp,  8.0_rp /)
  !>          xcoor(1:3,9) = (/  0.0_rp, -2.0_rp,  9.0_rp /)
  !>          print*,'ORIGINAL'
  !>          do ii = 1,lsize
  !>             print*,kpoin(ii),xcoor(1:3,ii)
  !>          end do
  !>          call maths_geometrical_sort_using_coordinates(2_ip,ndime,lsize,xcoor,kpoin)
  !>          print*,''
  !>          print*,'FULL SORTED'
  !>          do ii = 1,lsize
  !>             print*,kpoin(ii),xcoor(1:3,ii)
  !>          end do
  !>          deallocate(kpoin)
  !>          deallocate(xcoor)
  !>
  !-----------------------------------------------------------------------

  subroutine maths_geometrical_sort_using_coordinates(itask,pdime,nrows,xcoor,ivin,TOLERANCE)

    integer(ip), intent(in)             :: itask               !< Ordering
    integer(ip), intent(in)             :: pdime               !< Dimensions
    integer(ip), intent(inout)          :: nrows               !< Size list nodes
    real(rp),    intent(inout)          :: xcoor(pdime,nrows)  !< Coordinates
    integer(ip), intent(inout)          :: ivin(*)             !<
    real(rp),    intent(in),   optional :: TOLERANCE           !< Tolerance
    integer(ip)                         :: kk,mm,irows,jrows
    integer(ip)                         :: idime,jdime
    real(rp),    allocatable            :: xtmp(:,:)
    real(rp)                            :: my_toler
    integer(ip)                         :: perm(2,3)

    perm(1,1) = 2
    perm(2,1) = 3

    perm(1,2) = 1
    perm(2,2) = 3

    perm(1,3) = 1
    perm(2,3) = 2

    if( nrows < 2 ) return

    if( present(TOLERANCE) ) then
       my_toler = TOLERANCE
    else
       my_toler = epsil
    end if

    allocate(xtmp(nrows,pdime))
    do irows = 1,nrows
       xtmp(irows,:) = xcoor(:,irows)
    end do
    if( pdime == 2 ) then
       call maths_heap_sort_RP1(itask,nrows,xtmp(:,1),ivin,rvo1=xtmp(:,2))
    else
       call maths_heap_sort_RP1(itask,nrows,xtmp(:,1),ivin,rvo1=xtmp(:,2),rvo2=xtmp(:,3))
    end if

    do idime = 1,pdime-1
       irows = 0
       do while( irows < nrows )
          
          irows = irows + 1
          jrows = irows 
          kk    = 0
          do_kk: do while( jrows < nrows )
             jrows = jrows + 1
             if( idime == 1 ) then
                if( abs(xtmp(irows,1)-xtmp(jrows,1)) <= my_toler ) then
                   kk = kk + 1
                   if( jrows >= nrows ) exit do_kk
                else
                   exit do_kk
                end if
             else if( idime == 2 ) then
                if(  abs(xtmp(irows,1)-xtmp(jrows,1)) <= my_toler .and. &
                     abs(xtmp(irows,2)-xtmp(jrows,2)) <= my_toler ) then
                   kk = kk + 1
                  if( jrows >= nrows ) exit do_kk
                else
                   exit do_kk
                end if                
             end if
          end do do_kk
          
          if( kk > 0 ) then
             jrows = irows + kk
             mm    = kk + 1
             if( mm > 1 ) then
                jdime = idime + 1
                if( pdime == 2 ) then
                   call maths_heap_sort_RP1(itask,mm,&
                        xtmp(irows:jrows,jdime),&
                        ivin(irows:jrows),&
                        RVO1=xtmp(irows:jrows,perm(1,jdime)))
                else
                   call maths_heap_sort_RP1(itask,mm,&
                        xtmp(irows:jrows,jdime),&
                        ivin(irows:jrows),&
                        RVO1=xtmp(irows:jrows,perm(1,jdime)),&
                        RVO2=xtmp(irows:jrows,perm(2,jdime)))              
                end if
             end if
          end if
          irows = irows + kk 
       end do
    end do

1   continue 
    do irows = 1,nrows
       xcoor(:,irows) = xtmp(irows,:) 
    end do

    deallocate(xtmp)

  end subroutine maths_geometrical_sort_using_coordinates

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Geometrical ordering
  !> @details Uses a geometrical ordering based on SFC
  !>
  !-----------------------------------------------------------------------

  subroutine maths_geometrical_sort_using_sfc(itask,pdime,nrows,xcoor,ivin,NUMBER_BOXES)

    integer(ip), intent(in)           :: itask
    integer(ip), intent(in)           :: pdime
    integer(ip), intent(inout)        :: nrows
    real(rp),    intent(inout)        :: xcoor(pdime,nrows)
    integer(ip), intent(inout)        :: ivin(*)
    integer(ip), intent(in), optional :: NUMBER_BOXES
    integer(ip)                       :: boxes(3),ii,jj,kk,idime
    integer(ip)                       :: nsfc,nd,d,isame,irows
    real(rp)                          :: xleng(3)
    real(rp)                          :: xnorm,rsfc
    real(rp)                          :: xmini(3)
    real(rp)                          :: xmaxi(3)
    real(rp)                          :: toler
    integer(ip), allocatable          :: i2dto1d(:,:,:)
    integer(ip), allocatable          :: coord1d(:)

    if( nrows <= 1 ) then

       return

    else

       if( present(NUMBER_BOXES) ) then
          nsfc = NUMBER_BOXES
       else
          rsfc = real(nrows,rp)**(1.0_rp/real(pdime,rp))  ! Number of nodes in each direction
          nsfc = min(64_ip,2_ip**int(rsfc+1.0_rp,ip))
       end if
       !
       ! Sizes
       !
       toler = 1.0e-6_rp
       xleng = 0.0_rp
       xmini = 0.0_rp
       xmaxi = 0.0_rp
       do idime = 1,pdime
          xmini(idime) = minval(xcoor(idime,:))
          xmaxi(idime) = maxval(xcoor(idime,:))
       end do
       xleng(1:pdime) = xmaxi(1:pdime)-xmini(1:pdime)
       xnorm          = maxval(xleng)
       !
       ! Assign minimum size in each direction
       !
       xmini = xmini - 0.01_rp*xnorm
       xmaxi = xmaxi + 0.01_rp*xnorm
       !
       ! Look for 1D coordinates
       !
       allocate(coord1d(nrows))
       isame = 1
       do while( isame == 1 )

          boxes = nsfc
          if( pdime == 2 ) boxes(3) = 1
          nd = nsfc**pdime
          allocate(i2dto1d(boxes(1),boxes(2),boxes(3)))

          if( pdime == 2 ) then
             do d = 1,nd
                call maths_sfc_1d_to2d3d_tab(nsfc,d,ii,jj)
                i2dto1d(ii,jj,1) = d
             end do
             do irows = 1,nrows
                call maths_mapping_coord_to_3d(pdime,boxes,xmaxi,xmini,xcoor(:,irows),ii,jj,kk)
                coord1d(irows) =  i2dto1d(ii,jj,1)
             end do
          else
             do d = 1,nd
                call maths_sfc_1d_to2d3d_tab(nsfc,d,ii,jj,kk)
                i2dto1d(ii,jj,kk) = d
             end do
             do irows = 1,nrows
                call maths_mapping_coord_to_3d(pdime,boxes,xmaxi,xmini,xcoor(:,irows),ii,jj,kk)
                coord1d(irows) =  i2dto1d(ii,jj,kk)
             end do
          end if
          !
          ! Check if some nodes have the same 1D coordinate
          !
          isame = 0
          loop_irows: do irows = 1,nrows-1
             d     = coord1d(irows)
             isame = count(coord1d(irows+1:nrows)==d,KIND=ip)
             if( isame > 0 ) exit loop_irows
          end do loop_irows
          deallocate(i2dto1d)
          if( isame /= 0 ) nsfc = 2_ip * nsfc

       end do
       !
       ! Order according to 1D coordinate
       !
       call maths_heap_sort(itask,nrows,coord1d,' ',ivin)
       deallocate(coord1d)

    end if

  end subroutine maths_geometrical_sort_using_sfc
  
end module mod_maths_sort
!> @}
