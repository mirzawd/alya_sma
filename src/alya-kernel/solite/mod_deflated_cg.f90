!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solver
!> @{
!> @file    mod_deflated_cg.f90
!> @author  houzeaux
!> @date    2019-12-13
!> @brief   Deflated CG
!> @details Operations required by the Delfated CG
!-----------------------------------------------------------------------

module mod_deflated_cg

  use def_kintyp,         only : ip,rp
  use def_master,         only : kfl_paral,parre,nparr,INOTMASTER
  use def_master,         only : IPARALL
  use def_domain,         only : r_dom,c_dom
  use def_solver,         only : solve_sol
  use mod_communications, only : PAR_SUM
  
  implicit none
  private

  public :: matgro
  public :: matgr2
  public :: wtvect
  public :: wvect
  public :: wtvect_without_all_reduce
  
contains

  !------------------------------------------------------------------------
  !
  ! Operations for the deflated CG
  !
  !------------------------------------------------------------------------

  subroutine matgro(ngrou,npopo,nskyl,ndofn,ia,ja,an,askyl)
    !
    ! ASKYL: Factorize group matrix
    !
    integer(ip), intent(in)            :: ngrou,npopo,nskyl,ndofn
    integer(ip), intent(in)            :: ia(*),ja(*)
    real(rp),    intent(in)            :: an(ndofn,ndofn,*)
    real(rp),    intent(inout), target :: askyl(nskyl)
    integer(ip), pointer               :: iskyl(:),lgrou(:)
    integer(ip), pointer               :: iagro(:),jagro(:)
    integer(ip)                        :: igrou,jgrou,ipoin,izsym,jpoin,izgro
    integer(ip)                        :: kskyl,idofn,jdofn,igrou1,jgrou1

    if( ngrou == 0 ) return

    if( INOTMASTER ) then

       if( solve_sol(1) % kfl_defas == 1 ) then

          !----------------------------------------------------------------
          !
          ! Fill in sparse matrix ASKYL 
          !       
          !----------------------------------------------------------------

          lgrou => solve_sol(1) % lgrou
          iagro => solve_sol(1) % iagro
          jagro => solve_sol(1) % jagro

          if( ndofn == 1 ) then

             do ipoin = 1,npopo
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izsym = ia(ipoin),ia(ipoin+1)-1
                      jpoin = ja(izsym)
                      if( lgrou(jpoin) > 0 ) then
                         jgrou = lgrou(jpoin)
                         izgro = iagro(igrou)
                         iifzgro1: do while( izgro <= iagro(igrou+1)-1 )
                            if( jagro(izgro) == jgrou ) exit iifzgro1
                            izgro = izgro + 1
                         end do iifzgro1
                         askyl(izgro) = askyl(izgro) + an(1,1,izsym)
                      end if
                   end do
                end if
             end do

          else

             call runend('MATGRO: NOT CODED')

          end if

       else

          !----------------------------------------------------------------
          !
          ! Fill in skyline matrix ASKYL 
          !       
          !----------------------------------------------------------------

          lgrou => solve_sol(1) % lgrou
          iskyl => solve_sol(1) % iskyl

          if(ndofn==1)then

             do ipoin=1,npopo
                if(lgrou(ipoin)>0) then
                   igrou=lgrou(ipoin)
                   do izsym=ia(ipoin),ia(ipoin+1)-2
                      jpoin=ja(izsym)
                      if(lgrou(jpoin)>0) then
                         jgrou=lgrou(jpoin)
                         if(igrou>jgrou) then
                            kskyl=iskyl(igrou+1)-1-(igrou-jgrou)
                            askyl(kskyl)=askyl(kskyl)+an(1,1,izsym)
                         else if(igrou<jgrou) then
                            kskyl=iskyl(jgrou+1)-1-(jgrou-igrou)
                            askyl(kskyl)=askyl(kskyl)+an(1,1,izsym)
                         else
                            kskyl=iskyl(igrou+1)-1
                            askyl(kskyl)=askyl(kskyl)+2.0_rp*an(1,1,izsym)
                         end if
                      end if
                   end do
                   izsym=ia(ipoin+1)-1
                   kskyl=iskyl(igrou+1)-1
                   askyl(kskyl)=askyl(kskyl)+an(1,1,izsym)  
                end if
             end do

          else
             
             do ipoin=1,npopo
                if(lgrou(ipoin)>0) then
                   igrou=lgrou(ipoin)

                   do izsym=ia(ipoin),ia(ipoin+1)-2
                      jpoin=ja(izsym)
                      if(lgrou(jpoin)>0) then
                         jgrou=lgrou(jpoin)

                         do idofn=1,ndofn
                            do jdofn=1,idofn

                               igrou1=(igrou-1)*ndofn+idofn
                               jgrou1=(jgrou-1)*ndofn+jdofn

                               if(igrou1>jgrou1) then

                                  kskyl=iskyl(igrou1+1)-1-(igrou1-jgrou1)
                                  askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izsym)

                               else if(igrou1<jgrou1) then

                                  kskyl=iskyl(jgrou1+1)-1-(jgrou1-igrou1)
                                  askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izsym)

                               else

                                  kskyl=iskyl(igrou1+1)-1
                                  askyl(kskyl)=askyl(kskyl)+2.0_rp*an(idofn,jdofn,izsym)

                               endif
                            enddo
                         enddo
                      end if
                   enddo

                   izsym=ia(ipoin+1)-1

                   do idofn=1,ndofn
                      do jdofn=1,idofn

                         igrou1=(igrou-1)*ndofn+idofn
                         jgrou1=(igrou-1)*ndofn+jdofn
                         kskyl=iskyl(igrou1+1)-1-(igrou1-jgrou1)
                         askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izsym) 

                      enddo
                   enddo

                end if
             end do



          end if

       end if

    end if

    if( IPARALL ) then
       !
       ! Parallel: reduce sum
       !
       call PAR_SUM(nskyl,askyl,'IN MY CODE')
       !nparr =  nskyl
       !parre => askyl
       !call par_operat(3_ip)
    end if

  end subroutine matgro

  subroutine matgr2(ngrou,npopo,nskyl,ndofn,ia,ja,an,askyl)
    !
    ! ASKYL: Factorize group matrix
    !
    integer(ip), intent(in)            :: ngrou,npopo,nskyl,ndofn
    integer(ip), intent(in)            :: ia(*),ja(*)
    real(rp),    intent(in)            :: an(ndofn,ndofn,*)
    real(rp),    intent(inout)         :: askyl(*)
    integer(ip), pointer               :: iskyl(:),lgrou(:)
    integer(ip)                        :: igrou,jgrou,ipoin,izdom,jpoin,izgro
    integer(ip)                        :: kskyl,idofn,jdofn,igrou1,jgrou1
    integer(ip), pointer               :: iagro(:),jagro(:)

    if( ngrou == 0 ) return

    if( INOTMASTER ) then

       if( solve_sol(1) % kfl_defas == 1 ) then

          !----------------------------------------------------------------
          !
          ! Fill in sparse matrix ASKYL 
          !       
          !----------------------------------------------------------------

          lgrou => solve_sol(1) % lgrou
          iagro => solve_sol(1) % iagro
          jagro => solve_sol(1) % jagro

          if( ndofn == 1 ) then

             do ipoin = 1,npopo
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izdom = ia(ipoin),ia(ipoin+1)-1
                      jpoin = ja(izdom)
                      if( lgrou(jpoin) > 0 ) then
                         jgrou = lgrou(jpoin)
                         izgro = iagro(igrou)
                         iifzgro1: do while( jagro(izgro) /= jgrou )
                            izgro = izgro + 1
                         end do iifzgro1
                         askyl(izgro) = askyl(izgro) + an(1,1,izdom)
                      end if
                   end do
                end if
             end do

          else

             call assgr2(ngrou,npopo,nskyl,ndofn,ia,ja,an,askyl)

          end if

       else

          !----------------------------------------------------------------
          !
          ! Fill in skyline matrix ASKYL 
          !       
          !----------------------------------------------------------------

          lgrou => solve_sol(1) % lgrou
          iskyl => solve_sol(1) % iskyl

          if( ndofn==1 )then

             do ipoin = 1,npopo
                if( lgrou(ipoin) > 0 ) then
                   igrou = lgrou(ipoin)
                   do izdom = ia(ipoin),ia(ipoin+1)-1
                      jpoin = ja(izdom)
                      if( jpoin < ipoin ) then
                         if( lgrou(jpoin) > 0 ) then
                            jgrou = lgrou(jpoin)
                            if( igrou > jgrou ) then
                               kskyl        = iskyl(igrou+1) - 1 - (igrou-jgrou)
                               askyl(kskyl) = askyl(kskyl) + an(1,1,izdom)
                            else if( igrou < jgrou ) then
                               kskyl        = iskyl(jgrou+1) - 1 - (jgrou-igrou)
                               askyl(kskyl) = askyl(kskyl) + an(1,1,izdom)
                            else
                               kskyl        = iskyl(igrou+1) - 1
                               askyl(kskyl) = askyl(kskyl) + 2.0_rp*an(1,1,izdom)
                            end if
                         end if
                      else if( ipoin == jpoin ) then
                         kskyl        = iskyl(igrou+1) - 1
                         askyl(kskyl) = askyl(kskyl) + an(1,1,izdom)  
                      end if
                   end do
                end if
             end do

          else
             
             call runend('DEFALTED: SKYLINE NOT WORKING FOR NDOFN > 1 ')

             do ipoin=1,npopo
                if(lgrou(ipoin)>0) then
                   igrou=lgrou(ipoin)

                   do izdom=ia(ipoin),ia(ipoin+1)-2
                      jpoin=ja(izdom)
                      if(lgrou(jpoin)>0) then
                         jgrou=lgrou(jpoin)

                         do idofn=1,ndofn
                            do jdofn=1,idofn

                               igrou1=(igrou-1)*ndofn+idofn
                               jgrou1=(jgrou-1)*ndofn+jdofn

                               if(igrou1>jgrou1) then

                                  kskyl=iskyl(igrou1+1)-1-(igrou1-jgrou1)
                                  askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izdom)

                               else if(igrou1<jgrou1) then

                                  kskyl=iskyl(jgrou1+1)-1-(jgrou1-igrou1)
                                  askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izdom)

                               else

                                  kskyl=iskyl(igrou1+1)-1
                                  askyl(kskyl)=askyl(kskyl)+2.0_rp*an(idofn,jdofn,izdom)

                               endif
                            enddo
                         enddo
                      end if
                   enddo

                   izdom=ia(ipoin+1)-1

                   do idofn=1,ndofn
                      do jdofn=1,idofn

                         igrou1=(igrou-1)*ndofn+idofn
                         jgrou1=(igrou-1)*ndofn+jdofn
                         kskyl=iskyl(igrou1+1)-1-(igrou1-jgrou1)
                         askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izdom) 

                      enddo
                   enddo

                end if
             end do

          end if

       end if

    end if
    !
    ! Parallel: reduce sum
    !
    if( IPARALL ) then
       call PAR_SUM(nskyl,askyl,'IN MY CODE')
    end if

  end subroutine matgr2

  subroutine wtvect(npopo,ngrou,ndofn,xsmall,xbig)
    !
    ! XSMALL= W^T.XBIG 
    !
    use def_kintyp, only             :  ip,rp
    use def_master, only             :  IPARALL,parre,nparr,npoi1,npoi2,npoi3,&
         &                              NPOIN_REAL_12DI,ISLAVE,ISEQUEN
!    use def_master, only             :  icoml
    use def_solver, only             :  solve_sol
    implicit none
    integer(ip), intent(in)          :: npopo,ngrou,ndofn
    real(rp),    intent(in)          :: xbig(*)
    real(rp),    intent(out), target :: xsmall(ngrou*ndofn)
    integer(ip)                      :: ipoin,igrou,ipoin1,igrou1,idofn

    do igrou=1,solve_sol(1) % ngrou*ndofn
       xsmall(igrou)=0.0_rp
    end do

    if( ISEQUEN ) then

       if(ndofn==1)then  
          do ipoin=1,npopo
             igrou=solve_sol(1) % lgrou(ipoin)
             if(igrou>0) xsmall(igrou)=xsmall(igrou)+xbig(ipoin)
          end do
       else
          do ipoin=1,npopo
             igrou=solve_sol(1) % lgrou(ipoin)
             do idofn=1,ndofn
                igrou1=(igrou-1)*ndofn+idofn
                ipoin1=(ipoin-1)*ndofn+idofn
                if(igrou>0) xsmall(igrou1)=xsmall(igrou1)+xbig(ipoin1)
             enddo
          end do
       end if

    else if( ISLAVE ) then

       if(ndofn==1)then  

          do ipoin=1,npoi1
             igrou=solve_sol(1) % lgrou(ipoin)
             if(igrou>0) xsmall(igrou)=xsmall(igrou)+xbig(ipoin)
          end do

          do ipoin=npoi2,npoi3
             igrou=solve_sol(1) % lgrou(ipoin)
             if(igrou>0) xsmall(igrou)=xsmall(igrou)+xbig(ipoin)
          end do

       else

          do ipoin=1,npoi1
             igrou=solve_sol(1) % lgrou(ipoin)
             igrou1=(igrou-1)*ndofn
             ipoin1=(ipoin-1)*ndofn
             do idofn=1,ndofn 
                igrou1=igrou1+1
                ipoin1=ipoin1+1
                if(igrou>0) xsmall(igrou1)=xsmall(igrou1)+xbig(ipoin1)
             enddo
          end do

          do ipoin=npoi2,npoi3
             igrou=solve_sol(1) % lgrou(ipoin)
             igrou1=(igrou-1)*ndofn
             ipoin1=(ipoin-1)*ndofn
             do idofn=1,ndofn 
                igrou1=igrou1+1
                ipoin1=ipoin1+1
                if(igrou>0) xsmall(igrou1)=xsmall(igrou1)+xbig(ipoin1)
             enddo
          end do

       end if

    end if

    if( solve_sol(1) % kfl_defso == 1 .and. solve_sol(1) % kfl_defas == 2 ) return

    if( IPARALL ) then

       nparr =  solve_sol(1) % ngrou*ndofn
       parre => xsmall 

       if( solve_sol(1) % kfl_gathe == 1 ) then
          !
          ! Parallel: all gather v
          !
          !icoml =  solve_sol(1) % icoml
          call par_gatsen()

       else if( solve_sol(1) % kfl_gathe == 0 ) then
          !
          ! Parallel: reduce sum
          !
          call par_operat(3_ip)   

       else if( solve_sol(1) % kfl_gathe == 2 ) then
          !
          ! Parallel: send/receive
          !
          call par_cregrp(3_ip)

       end if

    end if

  end subroutine wtvect

  subroutine wvect(npopo,ndofn,xsmall,xbig)
    !
    ! XBIG= W.XSMALL 
    !
    use def_kintyp, only     :  ip,rp
    use def_solver, only     :  solve_sol
    use def_elmtyp, only     :  NOHOL
    use def_domain, only     :  lnoch
    implicit none
    integer(ip), intent(in)  :: npopo,ndofn
    real(rp),    intent(in)  :: xsmall(*)
    real(rp),    intent(out) :: xbig(*)
    integer(ip)              :: ipoin,igrou,ipoin1,igrou1,idofn

    if(ndofn==1)then

       do ipoin=1,npopo
          igrou=solve_sol(1) % lgrou(ipoin)
          if(igrou>0) then
             xbig(ipoin)=xsmall(igrou)
             if( lnoch(ipoin) == NOHOL ) xbig(ipoin)=0.0_rp
          else
             xbig(ipoin)=0.0_rp
          end if
       end do

    else

       do ipoin=1,npopo
          igrou=solve_sol(1) % lgrou(ipoin)
          igrou1=(igrou-1)*ndofn
          ipoin1=(ipoin-1)*ndofn
          do idofn=1,ndofn 
             igrou1=igrou1+1
             ipoin1=ipoin1+1
             if(igrou>0) then
                xbig(ipoin1)=xsmall(igrou1)
                if( lnoch(ipoin) == NOHOL ) xbig(ipoin1)=0.0_rp
             else
                xbig(ipoin1)=0.0_rp
             end if
          enddo
       end do

    endif

  end subroutine wvect

  subroutine wtvect_without_all_reduce(npopo,ngrou,ndofn,xsmall,xbig)
    !
    ! XSMALL= W^T.XBIG 
    !
    use def_kintyp, only             :  ip,rp
    use def_master, only             :  npoi1,npoi2,npoi3,&
         &                              NPOIN_REAL_12DI,ISLAVE,ISEQUEN
    use def_solver, only             :  solve_sol
    implicit none
    integer(ip), intent(in)          :: npopo,ngrou,ndofn
    real(rp),    intent(in)          :: xbig(*)
    real(rp),    intent(out), target :: xsmall(ngrou*ndofn)
    integer(ip)                      :: ipoin,igrou,ipoin1,igrou1,idofn

    do igrou=1,solve_sol(1) % ngrou*ndofn
       xsmall(igrou)=0.0_rp
    end do

    if( ISEQUEN ) then

       if(ndofn==1)then  
          do ipoin=1,npopo
             igrou=solve_sol(1) % lgrou(ipoin)
             if(igrou>0) xsmall(igrou)=xsmall(igrou)+xbig(ipoin)
          end do
       else
          do ipoin=1,npopo
             igrou=solve_sol(1) % lgrou(ipoin)
             do idofn=1,ndofn
                igrou1=(igrou-1)*ndofn+idofn
                ipoin1=(ipoin-1)*ndofn+idofn
                if(igrou>0) xsmall(igrou1)=xsmall(igrou1)+xbig(ipoin1)
             enddo
          end do
       end if

    else if( ISLAVE ) then

       if(ndofn==1)then  

          do ipoin=1,npoi1
             igrou=solve_sol(1) % lgrou(ipoin)
             if(igrou>0) xsmall(igrou)=xsmall(igrou)+xbig(ipoin)
          end do

          do ipoin=npoi2,npoi3
             igrou=solve_sol(1) % lgrou(ipoin)
             if(igrou>0) xsmall(igrou)=xsmall(igrou)+xbig(ipoin)
          end do

       else

          do ipoin=1,npoi1
             igrou=solve_sol(1) % lgrou(ipoin)
             igrou1=(igrou-1)*ndofn
             ipoin1=(ipoin-1)*ndofn
             do idofn=1,ndofn 
                igrou1=igrou1+1
                ipoin1=ipoin1+1
                if(igrou>0) xsmall(igrou1)=xsmall(igrou1)+xbig(ipoin1)
             enddo
          end do

          do ipoin=npoi2,npoi3
             igrou=solve_sol(1) % lgrou(ipoin)
             igrou1=(igrou-1)*ndofn
             ipoin1=(ipoin-1)*ndofn
             do idofn=1,ndofn 
                igrou1=igrou1+1
                ipoin1=ipoin1+1
                if(igrou>0) xsmall(igrou1)=xsmall(igrou1)+xbig(ipoin1)
             enddo
          end do

       end if

    end if

    if( solve_sol(1) % kfl_defso == 1 .and. solve_sol(1) % kfl_defas == 2 ) return

  end subroutine wtvect_without_all_reduce

end module mod_deflated_cg
!> @}
