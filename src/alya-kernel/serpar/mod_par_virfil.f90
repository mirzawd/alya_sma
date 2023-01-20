!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_par_virfil

  use def_parame
  use def_master
  use def_parall
  use mod_iofile
  use mod_memchk
  use mod_parall, only : par_memor
  integer(4) :: istat
  integer(8) :: mlenc,mleni,mlenr,mlenl

  public :: par_inibuf, par_wribuf, par_dumbuf, par_deabuf
  private

  type Buffer
     character(150)        :: file_path
     character(1), pointer :: buf_c(:)
     integer(ip),  pointer :: buf_i(:)
     real(rp),     pointer :: buf_r(:)
     logical(lg),  pointer :: buf_l(:)
     integer(ip)           :: npari(1000) ! '1000' max number of writes for any domain
     integer(ip)           :: nparr(1000)
     integer(ip)           :: nparc(1000)
     integer(ip)           :: nparl(1000)
     integer(ip)           :: count   = 1
     integer(ip)           :: c_start = 1
     integer(ip)           :: i_start = 1
     integer(ip)           :: r_start = 1
     integer(ip)           :: l_start = 1
     integer(ip)           :: mlenc
     integer(ip)           :: mleni
     integer(ip)           :: mlenr
     integer(ip)           :: mlenl
  end type Buffer

  type(Buffer), allocatable :: myBuf (:)


  contains

  subroutine par_inibuf(itask)
    implicit none
    integer(ip), intent(in) :: itask
    integer(ip)             :: i,i1,i2,maxin,maxre,maxlo
    real(rp)                :: xfact

    maxin = 5000
    maxre = 5000
    maxlo = 5000

    if( itask == -1 ) then
       i1    = 0
       i2    = npart_par
       allocate( myBuf(0:npart_par) )
   else
       i1    = itask
       i2    = itask
       maxin = max(maxin,npari)
       maxre = max(maxre,nparr)
       maxlo = max(maxlo,nparl)
    end if
    
    do i = i1,i2

       call par_filnam(1_ip,i,fil_rstar_par,myBuf(i)%file_path)
       ! 
       ! Max memory
       !               
       xfact            = real(rmbyt_par,rp)*1024.0_rp*1024.0_rp*1024.0_rp/real(npart_par*(ip+rp),rp)
       myBuf(i) % mlenc = 2000
       myBuf(i) % mleni = max(maxin,int(xfact,ip))
       myBuf(i) % mlenr = max(maxre,int(xfact,ip))
       myBuf(i) % mlenl = max(maxlo,int(xfact,ip))

       nullify(myBuf(i) % buf_c)
       nullify(myBuf(i) % buf_i)
       nullify(myBuf(i) % buf_r)
       nullify(myBuf(i) % buf_l)

       allocate( myBuf(i) % buf_c( myBuf(i) % mlenc ), stat = istat )
       call memchk(zero,istat,par_memor,'myBuf(i) % buf_c','par_inibuf',myBuf(i) % buf_c)
       allocate( myBuf(i) % buf_i( myBuf(i) % mleni ), stat = istat )
       call memchk(zero,istat,par_memor,'myBuf(i) % buf_i','par_inibuf',myBuf(i) % buf_i)
       allocate( myBuf(i) % buf_r( myBuf(i) % mlenr ), stat = istat )
       call memchk(zero,istat,par_memor,'myBuf(i) % buf_r','par_inibuf',myBuf(i) % buf_r)
       allocate( myBuf(i) % buf_l( myBuf(i) % mlenl ), stat = istat )
       call memchk(zero,istat,par_memor,'myBuf(i) % buf_l','par_inibuf',myBuf(i) % buf_l)

       myBuf(i) % count   = 1
       myBuf(i) % c_start = 1
       myBuf(i) % i_start = 1
       myBuf(i) % r_start = 1
       myBuf(i) % l_start = 1

    end do

  end subroutine par_inibuf

  subroutine par_wribuf(subDom)
    implicit none

    integer(ip)            :: subDom, w_count
    integer(ip)            :: c_start, i_start, r_start, l_start
    integer(ip)            :: c_final, i_final, r_final, l_final

    w_count = myBuf(subDom) % count
    c_start = myBuf(subDom) % c_start
    i_start = myBuf(subDom) % i_start
    r_start = myBuf(subDom) % r_start
    l_start = myBuf(subDom) % l_start
    !
    ! Length
    !
    mlenc = int(myBuf(subDom) % mlenc,8)
    mleni = int(myBuf(subDom) % mleni,8)
    mlenr = int(myBuf(subDom) % mlenr,8) 
    mlenl = int(myBuf(subDom) % mlenl,8)
    !
    ! Save sizes
    !
    c_final = c_start + nparc - 1
    i_final = i_start + npari - 1
    r_final = r_start + nparr - 1
    l_final = l_start + nparl - 1
    !
    ! Reallocate if we are too short in memory
    !
    if( 1 == 1 ) then
       if( int(c_final,8) > mlenc .or. int(i_final,8) > mleni .or. int(r_final,8) > mlenr ) then
          call par_dumbuf(subDom)
          call par_inibuf(subDom)
          w_count = myBuf(subDom) % count
          c_start = myBuf(subDom) % c_start
          i_start = myBuf(subDom) % i_start
          r_start = myBuf(subDom) % r_start
          l_start = myBuf(subDom) % l_start
          c_final = c_start + nparc - 1
          i_final = i_start + npari - 1
          r_final = r_start + nparr - 1
          l_final = l_start + nparl - 1
       end if
       myBuf(subDom) % nparc(w_count) = nparc
       myBuf(subDom) % npari(w_count) = npari
       myBuf(subDom) % nparr(w_count) = nparr
       myBuf(subDom) % nparl(w_count) = nparl
    else
       myBuf(subDom) % nparc(w_count) = nparc
       myBuf(subDom) % npari(w_count) = npari
       myBuf(subDom) % nparr(w_count) = nparr
       myBuf(subDom) % nparl(w_count) = nparl
       if( int(c_final,8) > mlenc ) then
          !myBuf( subDom ) % buf_c = memrea(i_final,par_memor,'buf_c','par_wribuf',myBuf( subDom ) % buf_s)
          call runend('PAR_WRIBUF: INCREASE CHARACTER MEMORY FOR VIRTUAL FILE')    
       end if
       if( int(i_final,8) > mleni ) then
          call memrea(i_final,par_memor,'buf_i','par_wribuf',myBuf( subDom ) % buf_i)
       end if
       if( int(r_final,8) > mlenr ) then
          call memrea(r_final,par_memor,'buf_r','par_wribuf',myBuf( subDom ) % buf_r)
       end if
       if( int(l_final,8) > mlenl ) then
          call memrea(l_final,par_memor,'buf_l','par_wribuf',myBuf( subDom ) % buf_l)
       end if
    end if
    !
    ! Save content
    !
    if( nparc > 0 ) myBuf( subDom ) % buf_c( c_start:c_final ) = parch( 1:nparc )
    if( npari > 0 ) myBuf( subDom ) % buf_i( i_start:i_final ) = parin( 1:npari )
    if( nparr > 0 ) myBuf( subDom ) % buf_r( r_start:r_final ) = parre( 1:nparr )
    if( nparl > 0 ) myBuf( subDom ) % buf_l( l_start:l_final ) = parlo( 1:nparl )
    !
    ! Count++
    !
    myBuf(subDom) % count = myBuf(subDom) % count + 1
    !
    ! Starts++
    !
    myBuf(subDom) % c_start = c_final + 1
    myBuf(subDom) % i_start = i_final + 1
    myBuf(subDom) % r_start = r_final + 1
    myBuf(subDom) % l_start = l_final + 1

  end subroutine par_wribuf

  subroutine par_dumbuf(itask)
    implicit none
    !
    ! Dump buffer to restart files
    !
    integer(ip), intent(in) :: itask
    integer(ip)             :: ii, jj, kk, iunit, ipari, i1, i2
    integer(ip)             :: c_start, i_start, r_start, l_start
    integer(ip)             :: c_stop,  i_stop,  r_stop,  l_stop
    integer(4)              :: iunit4
    integer(4), pointer     :: parin4(:)
    character(20)           :: cdum1,cdum2

    if( itask == -1 ) then
       i1 = 0
       i2 = npart_par
    else
       i1 = itask
       i2 = itask
    end if

    do ii = i1,i2

       c_start = 1
       i_start = 1
       r_start = 1 
       l_start = 1 

       iunit   = lun_aonlp_par + ii 
       iunit4  = int(iunit,4)
       cdum1   = intost(ii)
       cdum2   = 'PARALL RESTART '//trim(cdum1)

       if( kfl_filio_par == 1 ) &           
            call iofile(zero,iunit,adjustl(trim(myBuf(ii)%file_path)),trim(cdum2),'old','unformatted','append')

       do jj=1, myBuf(ii) % count-1

          c_stop = myBuf(ii) % nparc(jj)
          i_stop = myBuf(ii) % npari(jj)
          r_stop = myBuf(ii) % nparr(jj)
          l_stop = myBuf(ii) % nparl(jj)

          if( kfl_bytes_par == 4 .and. ip /= 4 ) then
             !
             ! Output in INTEGER(4)
             !
             allocate(parin4(npari),stat=istat) 
             call memchk(zero,istat,par_memor,'parin4','par_dumbuf',parin4)
             kk = 0
             do ipari = i_start,i_start+i_stop-1
                kk = kk + 1
                parin4(ipari) = int( myBuf(kk) % buf_i(ipari),4)
             end do
             write(iunit4) int(myBuf(ii) % npari(jj),4), int(myBuf(ii) % nparr(jj),4), int(myBuf(ii) % nparc(jj),4)
             if( i_stop > 0 ) write(iunit4) parin4( 1 : kk ) 
          else
             ! 
             ! Output in INTEGER(IP)
             !
             write(iunit4) myBuf(ii) % npari(jj), myBuf(ii) % nparr(jj), myBuf(ii) % nparc(jj)
             if( i_stop > 0 ) write(iunit4) myBuf(ii) % buf_i( i_start : i_start+i_stop-1 ) 
          end if

          if( r_stop > 0 ) write(iunit4) myBuf(ii) % buf_r( r_start : r_start+r_stop-1 )
          if( c_stop > 0 ) write(iunit4) myBuf(ii) % buf_c( c_start : c_start+c_stop-1 )
          if( l_stop > 0 ) write(iunit4) myBuf(ii) % buf_l( l_start : l_start+l_stop-1 )
           
          c_start = c_start + c_stop
          i_start = i_start + i_stop
          r_start = r_start + r_stop
          l_start = l_start + l_stop
           
          if( kfl_bytes_par == 4 .and. ip /= 4 ) then
             call memchk(two,istat,par_memor,'parin4','par_dumbuf',parin4)
             deallocate(parin4,stat=istat) 
             if(istat/=0) call memerr(two,'PARIN4','par_dumbuf',0_ip)
          end if

       end do

       if( kfl_filio_par == 1 ) close(iunit4)
       
    end do
    !
    ! Deallocate buffer
    !
    call par_deabuf(itask)

  end subroutine par_dumbuf

  subroutine par_deabuf(itask)
    implicit none
    integer(ip), intent(in) :: itask
    integer(ip)             :: ii,i1,i2

    if( itask == -1 ) then
       i1 = 0
       i2 = npart_par
    else
       i1 = itask
       i2 = itask
    end if

    do ii = i1,i2
       !deallocate( myBuf(ii) % npari, stat = istat )
       !deallocate( myBuf(ii) % nparr, stat = istat )
       !deallocate( myBuf(ii) % nparc, stat = istat )
       deallocate( myBuf(ii) % buf_i, stat = istat )
       deallocate( myBuf(ii) % buf_r, stat = istat )
       deallocate( myBuf(ii) % buf_c, stat = istat )
       deallocate( myBuf(ii) % buf_l, stat = istat )
    end do
    if( itask == -1 ) deallocate( myBuf, stat=istat )

  end subroutine par_deabuf

end module mod_par_virfil
