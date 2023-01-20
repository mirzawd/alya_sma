!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine suplot(ndofn,unkno,luspl)
  !-----------------------------------------------------------------------
  !****f* kernel/outrut
  ! NAME 
  !    nsi_matrix
  ! DESCRIPTION
  !    This routine writes a file with a scalar function to be plotted
  !    by gnuplot using SPLOT. The domain MUST BE [a,b] x [a,b]. 
  !
  !    - In gnuplot, draw the 3d elevation with:
  !      splot 'problem.spl' w l
  !
  ! USES
  ! USED BY
  !    ***_output
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_domain
  use      mod_memchk
  implicit none
  integer(ip), intent(in)  :: ndofn,luspl
  real(rp),    intent(in)  :: unkno(ndofn,npoin) 
  integer(ip)              :: ndiv,idivx,idivy,ipoin,idofn
  real(rp)                 :: y,dy,yact,ymin,ymax,ynext,x,dx
  real(rp)                 :: xact,xmin,xnext,xmax,d2,zeroc
  integer(ip)              :: found
  integer(4)               :: istat
  integer(ip), allocatable :: list(:,:)

  ndiv = int(sqrt(real(npoin,rp)+0.5_rp))
  allocate(list(ndiv,ndiv),stat=istat) 
  call memchk(zero,istat,memor_dom,'LIST','suplot',list)

  xmin=1.0e10_rp
  xmax=0.0_rp
  ymin=1.0e10_rp
  ymax=0.0_rp
  zeroc=1.0e-8_rp
  found=0

  do ipoin=1,npoin
     x=coord(1,ipoin)
     y=coord(2,ipoin)
     xmin=min(xmin,x)
     xmax=max(xmax,x)
     ymin=min(ymin,y)
     ymax=max(ymax,y)
  end do
  xnext=xmax
  ynext=ymax

  yact =ymin
  do idivy=1,ndiv
     xact =xmin
     do idivx=1,ndiv
        ipoin=0
        do while((found==0).and.(ipoin.lt.npoin))
           ipoin=ipoin+1
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           dx = (x-xact)*(x-xact)
           dy = (y-yact)*(y-yact)
           d2 = dx + dy
           if(d2.lt.zeroc) then
              found=1
              list(idivx,idivy) = ipoin
           end if
        end do
        if(found==0) then
           write(luspl,*) 'Something goes wrong, point not found'
           stop
        else if (found==1) then
           do ipoin=1,npoin
              dy=yact-coord(2,ipoin)
              dy=dy*dy
              if(dy.lt.zeroc) then
                 x=coord(1,ipoin)
                 if(x>xact+zeroc) xnext=min(xnext,x)
              end if
           end do
        end if
        found=0
        xact=xnext
        xnext=xmax
     end do
     do ipoin=1,npoin
        y=coord(2,ipoin)
        if(y>yact+zeroc) ynext=min(ynext,y)
     end do
     yact=ynext
     ynext=ymax
  end do

  do idivy=1,ndiv
     do idivx=1,ndiv
        ipoin=list(idivx,idivy)
        write(luspl,'(20(3x,f15.5))')&
             coord(1,ipoin),coord(2,ipoin),&
             (unkno(idofn,ipoin),idofn=1,ndofn)
        !             unkno(ipoin)
     end do
     write(luspl,*)
  end do
  call memchk(two,istat,memor_dom,'LIST','suplot',list)
  deallocate(list,stat=istat) 
  if(istat.ne.0)  call memerr(two,'LIST','suplot',0_ip)


end subroutine suplot


