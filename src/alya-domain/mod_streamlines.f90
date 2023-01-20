!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup CPU_Time 
!> @{
!> @file    mod_timings.f90
!> @author  houzeaux
!> @date    2019-03-01
!> @brief   Compute assembly timings
!> @details Assembly timings
!>          
!-----------------------------------------------------------------------

module mod_streamlines

  use def_kintyp_basic, only : ip,rp,lg
  use def_kintyp_dims,  only : nelty
  use def_domain,       only : npoin
  use def_domain,       only : nelem
  use def_domain,       only : lnods  
  use def_domain,       only : ltype  
  use def_domain,       only : coord  
  use def_domain,       only : ltopo  
  use def_domain,       only : nnode  
  use def_domain,       only : ndime  
  use def_domain,       only : mnode  
  use def_domain,       only : lexis  
  use def_domain,       only : memor_dom  
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
    
  implicit none

  private

  public :: strfun

contains

  subroutine strfun(vecto,strea)
    !-----------------------------------------------------------------------
    !****f* Mathru/strfun
    ! NAME 
    !    strfun
    ! DESCRIPTION
    !    This routine computes the stream function
    ! USES
    !    strloc
    ! USED BY
    !    nsi_output
    !***
    !-----------------------------------------------------------------------

    real(rp),    pointer, intent(in)     :: vecto(:,:,:)
    real(rp),    pointer, intent(inout)  :: strea(:)
    real(rp)                             :: elvel(ndime,mnode),elcod(ndime,mnode)
    integer(ip)                          :: ielem,kpoin,itouc,jnode,iesta,iesto
    integer(ip)                          :: jpoin,ipoin,inode,pelty,pnode
    integer(ip), pointer                 :: mpoin(:)
    integer(ip), pointer                 :: nocha(:,:),noinv(:,:)
    real(rp),    pointer                 :: westr(:,:)

    nullify(mpoin,nocha,noinv,westr)
    
    if(npoin/=0) then
       !
       ! Allocate memory
       !
       call memory_alloca(memor_dom,'NOCHA','strfun',nocha,16_ip,nelty)
       call memory_alloca(memor_dom,'NOINV','strfun',noinv,16_ip,nelty)
       call memory_alloca(memor_dom,'WESTR','strfun',westr,16_ip,nelty)
       call memory_alloca(memor_dom,'MPOIN','strfun',mpoin,npoin)
       if(ndime==2) then
          iesta=10
          iesto=29
       else if(ndime==3) then
          iesta=30
          iesto=50
       end if
       do pelty=iesta,iesto
          if(lexis(pelty)/=0)&
               call chanum(&
               ltopo(pelty),nnode(pelty),nocha(:,pelty),&
               noinv(:,pelty),westr(:,pelty))
       end do
       !
       ! Compute stream function
       !    
       ielem=0
       kpoin=1
       mpoin(1)=1
       strea(1)=0.0_rp
       do while(kpoin<npoin)
          ielem=mod(ielem+1,nelem)
          if(ielem==0) ielem=nelem
          pelty=ltype(ielem)
          pnode=nnode(pelty)     
          itouc=0
          do jnode=1,pnode
             jpoin=lnods(jnode,ielem)
             if(mpoin(jpoin)>=1) itouc=itouc+1 
          end do
          if(itouc>0.and.itouc<pnode) then 
             do inode=1,pnode
                ipoin=lnods(inode,ielem)
                elvel(1:ndime,inode) = vecto(1:ndime,ipoin,1)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do
             call strloc(pnode,ndime,npoin,&
                  lnods(:,ielem),mpoin,kpoin,elcod,elvel,strea,&
                  nocha(:,pelty),noinv(:,pelty),westr(:,pelty))
          end if
       end do
       do ipoin=1,npoin
          strea(ipoin)=strea(ipoin)/real(mpoin(ipoin),rp)
       end do
       !
       ! Deallocate memory
       !
       call memory_deallo(memor_dom,'NOCHA','strfun',nocha)
       call memory_deallo(memor_dom,'NOINV','strfun',noinv)
       call memory_deallo(memor_dom,'WESTR','strfun',westr)
       call memory_deallo(memor_dom,'MPOIN','strfun',mpoin)

    end if

  end subroutine strfun

  subroutine strloc(pnode,ndime,npoin,&
       lnods,lpoin,kpoin,elcod,elvel,strea,&
       nocha,noinv,westr)
    !-----------------------------------------------------------------------
    !****f* Mathru/strloc
    ! NAME 
    !    strloc
    ! DESCRIPTION
    !    This routine computes the stream-function for each element
    ! USES
    ! USED BY
    !    stream 
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)    :: pnode,ndime,npoin
    real(rp),    intent(in)    :: elcod(ndime,pnode),elvel(ndime*pnode)
    integer(ip), intent(in)    :: lnods(pnode)
    integer(ip), intent(in)    :: nocha(pnode),noinv(pnode)
    real(rp),    intent(in)    :: westr(pnode-1)
    integer(ip), intent(inout) :: lpoin(npoin)
    real(rp),    intent(inout) :: strea(npoin)
    integer(ip)                :: kpoin,inoco,itouc,ipoco,mnodb,mnoda
    integer(ip)                :: icoun,ivxin
    integer(ip)                :: ivxfi,jpoin,ipini,ipfin,inode,lapoi
    real(rp)                   :: xcomp,ycomp,unori,unorf,strfa,preno
    real(rp)                   :: strep
    logical(lg)                :: lofin
    !
    ! Identify first node where the streamfunction is known
    !
    inoco=0
    itouc=0
    do while(itouc==0)
       inoco=inoco+1
       ipoco=lnods(inoco)
       if(lpoin(ipoco)>=1) itouc=1
    end do
    !
    ! Compute the streamfunction at the nodes on the edges of the element
    !
    inoco=nocha(inoco)
    do icoun=1,pnode
       mnodb=inoco+icoun
       lofin=.false.
       do while((.not.lofin).and.(mnodb<=pnode))
          if (nocha(mnodb)==0) then
             mnodb=mnodb+1
          else
             lofin=.true.
          end if
       end do
       if(mnodb>pnode) mnodb=mod(mnodb,pnode)
       if(mnodb==0) mnodb=pnode
       mnoda=mnodb-1
       if(mnoda==0) then
          mnoda=pnode
          do while(nocha(mnoda)==0)
             mnoda=mnoda-1
          end do
       end if
       mnodb=noinv(mnodb)
       mnoda=noinv(mnoda)
       xcomp=elcod(1,mnodb)-elcod(1,mnoda)
       ycomp=elcod(2,mnodb)-elcod(2,mnoda)
       ivxin=(mnoda-1)*2+1
       ivxfi=(mnodb-1)*2+1
       unori=-elvel(ivxin)*ycomp+elvel(ivxin+1)*xcomp
       unorf=-elvel(ivxfi)*ycomp+elvel(ivxfi+1)*xcomp
       ipini=lnods(mnoda)
       ipfin=lnods(mnodb)
       strfa=strea(ipini)/real(lpoin(ipini),rp)
       strea(ipfin)=strea(ipfin)+strfa-0.5_rp*(unori+unorf)
       if(lpoin(ipfin)==0) kpoin=kpoin+1
       lpoin(ipfin)=lpoin(ipfin)+1
    end do
    ! 
    ! Compute the streamfunction at the interior nodes. For the cubic elements
    ! P3 and Q3 (10 and 16 nodes), a special routine is called. If the same
    ! algorithm as for the rest of elements is to be used, only the corner
    ! nodes are employed to compute the streamfunction (see routine CHANUM)
    !
    if((pnode==10).or.(pnode==16)) then
       call strcub(pnode,npoin,kpoin,lpoin,lnods,strea)
    else 
       if(nocha(pnode)==0) then
          preno=0.0_rp
          do inode=1,min(pnode-1_ip,8_ip)
             jpoin=lnods(inode)
             strep=strea(jpoin)/real(lpoin(jpoin),rp)
             preno=preno+westr(inode)*strep
          end do
          lapoi=lnods(pnode)
          strea(lapoi)=preno  
          lpoin(lapoi)=lpoin(lapoi)+1
          kpoin=kpoin+1
       end if
    end if

  end subroutine strloc

  subroutine strcub(pnode,npoin,kpoin,lpoin,lnods,strea)
    !-----------------------------------------------------------------------
    !****f* Mathru/strcub
    ! NAME 
    !    strcub
    ! DESCRIPTION
    !    This routine computes the stream-function at the interior points of
    !    the 16-noded cubic quadrilateral and the 10-noded cubic triangle
    ! USES
    ! USED BY
    !    strloc
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)    :: pnode,npoin
    integer(ip), intent(in)    :: lnods(pnode)
    integer(ip), intent(inout) :: lpoin(npoin),kpoin
    real(rp),    intent(inout) :: strea(npoin)
    integer(ip)                :: inode,ipoin,ipo10
    real(rp)                   :: contr
    real(rp)                   :: strlo(12)
    !
    ! 16-noded element
    !
    if(pnode==16) then
       do inode=1,12
          ipoin=lnods(inode)
          strlo(inode)=strea(ipoin)/real(lpoin(ipoin),rp)
       end do
       strea(lnods(13))=strlo(12)+strlo( 5)-strlo(1)
       strea(lnods(14))=strlo( 6)+strlo( 7)-strlo(2)      
       strea(lnods(15))=strlo( 8)+strlo( 9)-strlo(3)
       strea(lnods(16))=strlo(10)+strlo(11)-strlo(4)      
       lpoin(lnods(13))=lpoin(lnods(13))+1
       lpoin(lnods(14))=lpoin(lnods(14))+1
       lpoin(lnods(15))=lpoin(lnods(15))+1 
       lpoin(lnods(16))=lpoin(lnods(16))+1
       kpoin = kpoin + 4
       !
       ! 10-noded element
       !
    else if(pnode==10) then
       ipo10=lnods(10)
       strea(ipo10)=0.0_rp
       do inode=1,3
          ipoin=lnods(inode)
          contr=strea(ipoin)/real(lpoin(ipoin),rp)
          strea(ipo10)=strea(ipo10)-contr
       end do
       do inode=4,9
          ipoin=lnods(inode)
          contr=strea(ipoin)/real(lpoin(ipoin),rp)
          strea(ipo10)=strea(ipo10)+contr
       end do
       strea(ipo10)=strea(ipo10)/3.0_rp
       lpoin(ipo10)=lpoin(ipo10)+1
       kpoin = kpoin + 1
    end if

  end subroutine strcub

  subroutine chanum(ptopo,pnode,nocha,noinv,westr)

    !-----------------------------------------------------------------------
    !****f* Mathru/chanum
    ! NAME 
    !    chanum
    ! DESCRIPTION
    !    This routine changes nodal numbering for stream-function algorithm
    !    and sets weigths for central node
    ! USES
    ! USED BY
    !    stream 
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: ptopo,pnode
    integer(ip), intent(out) :: nocha(pnode), noinv(pnode)
    real(rp),    intent(out) :: westr(pnode-1)
    integer(ip)              :: istop,inode
    !
    ! Initialization
    !
    do inode=1,pnode
       nocha(inode)=inode
       noinv(inode)=inode
    end do
    istop=0
    !
    ! Quads & bricks
    !
    if(ptopo==0) then
       if(pnode==4) then
          continue
       else if(pnode==5) then
          nocha(pnode)=0
          noinv(pnode)=0
          do inode=1,pnode-1
             westr(inode)=0.25_rp
          end do
       else if((pnode==8).or.(pnode==9)) then
          nocha(2)=3
          nocha(3)=5
          nocha(4)=7
          nocha(5)=2
          nocha(6)=4
          nocha(7)=6
          noinv(2)=5
          noinv(3)=2
          noinv(4)=6
          noinv(5)=3
          noinv(6)=7
          noinv(7)=4
          if(pnode==9) then
             nocha(pnode)=0
             noinv(pnode)=0
             do inode=1,4
                westr(inode  )=-0.25_rp
                westr(inode+4)= 0.50_rp
             end do
          end if
       else if(pnode==16) then
          nocha(2 )=4
          nocha(3 )=7
          nocha(4 )=10
          nocha(5 )=2
          nocha(6 )=3
          nocha(7 )=5
          nocha(8 )=6
          nocha(9 )=8
          nocha(10)=9
          nocha(11)=11
          nocha(12)=12
          nocha(13)=0
          nocha(14)=0
          nocha(15)=0
          nocha(16)=0
          noinv(2 )=5
          noinv(3 )=6
          noinv(4 )=2
          noinv(5 )=7
          noinv(6 )=8
          noinv(7 )=3
          noinv(8 )=9
          noinv(9 )=10
          noinv(10)=4
          noinv(11)=11
          noinv(12)=12
          noinv(13)=0
          noinv(14)=0
          noinv(15)=0
          noinv(16)=0
          do inode=1,4
             westr(inode)=0.25_rp
          end do
          do inode=5,8
             westr(inode)=0.0_rp
          end do
       else
          istop=1
       end if
       !
       ! Triangles & tetrahedra
       !
    else if(ptopo==1) then
       if(pnode==3) then
          continue
       else if(pnode==4) then
          nocha(pnode)=0
          noinv(pnode)=0
          do inode=1,pnode-1
             westr(inode)=1.0_rp/3.0_rp
          end do
       else if((pnode==6).or.(pnode==7)) then
          nocha(2)=3
          nocha(3)=5
          nocha(4)=2
          nocha(5)=4
          noinv(2)=4
          noinv(3)=2
          noinv(4)=5
          noinv(5)=3
          if(pnode==7) then
             nocha(pnode)=0
             noinv(pnode)=0
             do inode=1,3
                westr(inode  )=-1.0_rp/9.0_rp
                westr(inode+3)= 4.0_rp/9.0_rp
             end do
          end if
       else if(pnode==10) then
          nocha(2 )=4
          nocha(3 )=7
          nocha(4 )=2
          nocha(5 )=3
          nocha(6 )=5
          nocha(7 )=6
          nocha(10)=0
          noinv(2 )=4
          noinv(3 )=5
          noinv(4 )=2
          noinv(5 )=6
          noinv(6 )=7
          noinv(7 )=3
          noinv(10)=0
          do inode=1,3
             westr(inode)=1.0_rp/3.0_rp
          end do
          do inode=4,8
             westr(inode)=0.0_rp
          end do
       else
          istop=1
       end if
    end if
    !
    ! Error
    !
    if(istop==1)&
         call runend('CHANUM: STREAMFUNCTION FACILITY NOT AVAILABLE')

  end subroutine chanum

end module mod_streamlines
!> @}
