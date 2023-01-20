!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!
! Operations for the deflated CG
!
!------------------------------------------------------------------------

subroutine matgrx2(ngrou,npopo,nskyl,ndofn,ia,ja,an,askyl)
  !
  ! ASKYL: Factorize group matrix
  !
  use def_kintyp, only               :  ip,rp
  use def_master, only               :  INOTMASTER,parcx,nparx,IPARALL
  use def_solver, only               :  solve_sol
  use mod_memchk
  implicit none
  integer(ip), intent(in)            :: ngrou,npopo,nskyl,ndofn
  integer(ip), intent(in)            :: ia(*),ja(*)
  complex(rp), intent(in)            :: an(ndofn,ndofn,*)
  complex(rp), intent(inout), target :: askyl(nskyl)
  integer(ip), pointer               :: iskyl(:),lgrou(:)
  integer(ip)                        :: igrou,jgrou,ipoin,izdom,jpoin,izgro,izgro1
  integer(ip)                        :: info,kskyl,idofn,jdofn,igrou1,jgrou1
  integer(ip), pointer               :: iagro(:),jagro(:)

  if( ngrou==0 ) return

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
      
            do ipoin = 1,npopo
              if( lgrou(ipoin) > 0 ) then
                 igrou = lgrou(ipoin)
                 do izdom = ia(ipoin),ia(ipoin+1)-1
                    jpoin = ja(izdom)
                    if( lgrou(jpoin) > 0 ) then
                       jgrou = lgrou(jpoin)
                       izgro = iagro(igrou)
                       iifzgro2: do while( jagro(izgro) /= jgrou )
                          izgro = izgro + 1
                       end do iifzgro2
                       do jdofn = 1,ndofn
                         do idofn = 1,ndofn
                           izgro1 = (izgro-1)*ndofn*ndofn+(jdofn-1)*ndofn+idofn
                           askyl(izgro1) = askyl(izgro1) + an(idofn,jdofn,izdom)
                         end do
                       end do  
                    end if
                 end do
              end if
           end do         
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

  if( IPARALL ) then
     !
     ! Parallel: reduce sum
     !
     nparx =  nskyl
     parcx => askyl
     call par_operat(3_ip)
  end if
  !
  ! Inverse matrix ASKYL
  !
  if( INOTMASTER .and. solve_sol(1) % kfl_defas == 0 ) then
     !call chofac(ngrou*ndofn,nskyl,iskyl,askyl,info)
     call runend('MATGR2: ERROR WHILE DOING CHOLESKY FACTORIZATION')
     !if(info/=0) call runend('MATGR2: ERROR WHILE DOING CHOLESKY FACTORIZATION')
  end if

end subroutine matgrx2

subroutine wtvectx(npopo,ngrou,ndofn,xsmall,xbig)
  !
  ! XSMALL= W^T.XBIG 
  !
  use def_kintyp, only             :  ip,rp
  use def_master, only             :  IPARALL,parcx,nparx,npoi1,npoi2,npoi3,&
       &                              NPOIN_REAL_12DI,icoml,ISLAVE,ISEQUEN
!  use def_master, only             :  kfl_paral
  use def_solver, only             :  solve_sol
  implicit none
  integer(ip), intent(in)          :: npopo,ngrou,ndofn
  complex(rp), intent(in)          :: xbig(*)
  complex(rp), intent(out), target :: xsmall(ngrou*ndofn)
  integer(ip)                      :: ipoin,igrou,ipoin1,igrou1,idofn
!  real(rp)                         :: cpu_smal1,cpu_smal2,cpu_smal

  do igrou=1,solve_sol(1) % ngrou*ndofn
     xsmall(igrou)=(0.0_rp,0.0_rp)
  end do

  if( ISEQUEN ) then
  !!call cputim(cpu_smal1)
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
  !!call cputim(cpu_smal2)
  !!cpu_smal = cpu_smal2 - cpu_smal1
  !!write(*,*)'Time to small1:',cpu_smal,kfl_paral

  else if( ISLAVE ) then
  !!call cputim(cpu_smal1)
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
  !!call cputim(cpu_smal2)
  !!cpu_smal = cpu_smal2 - cpu_smal1
  !!write(*,*)'Time to small1:',cpu_smal,kfl_paral
  end if

  if( IPARALL ) then
  !!call cputim(cpu_smal1)
     nparx =  solve_sol(1) % ngrou*ndofn
     parcx => xsmall 

     if( solve_sol(1) % kfl_gathe == 1 ) then
        !
        ! Parallel: all gather v
        !
        icoml =  solve_sol(1) % icoml
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
  !!call cputim(cpu_smal2)
  !!cpu_smal = cpu_smal2 - cpu_smal1
  !!write(*,*)'Time to small exch:',cpu_smal,kfl_paral
  end if

end subroutine wtvectx

subroutine wvectx(npopo,ndofn,xsmall,xbig)
  !
  ! XBIG= W.XSMALL 
  !
  use def_kintyp, only     :  ip,rp
  use def_solver, only     :  solve_sol
  use def_domain, only     :  nimbo
  use def_master, only     :  ID_DODEME,lntib,kfl_coibm
  implicit none
  integer(ip), intent(in)  :: npopo,ndofn
  complex(rp), intent(in)  :: xsmall(*)
  complex(rp), intent(out) :: xbig(*)
  integer(ip)              :: ipoin,igrou,ipoin1,igrou1,idofn

  if(ndofn==1)then
     
     if( nimbo > 0 ) then

        call runend('DEFOPE: UNIFY DODEME AND IMMBOU')

     else if( nimbo > 0 .and. kfl_coibm == 1 ) then

        do ipoin=1,npopo
           igrou=solve_sol(1) % lgrou(ipoin)
           if(igrou>0) then
              xbig(ipoin)=xsmall(igrou)
              if( lntib(ipoin) > 0 ) xbig(ipoin)=(0.0_rp,0.0_rp)
           else
              xbig(ipoin)=(0.0_rp,0.0_rp)
           end if
        end do

     else

        do ipoin=1,npopo
           igrou=solve_sol(1) % lgrou(ipoin)
           if(igrou>0) then
              xbig(ipoin)=xsmall(igrou)
           else
              xbig(ipoin)=(0.0_rp,0.0_rp)
           end if
        end do

     end if

  else

     if( nimbo > 0 ) then
        do ipoin=1,npopo
           igrou=solve_sol(1) % lgrou(ipoin)
           igrou1=(igrou-1)*ndofn
           ipoin1=(ipoin-1)*ndofn
           do idofn=1,ndofn 
              igrou1=igrou1+1
              ipoin1=ipoin1+1
              if(igrou>0) then
                 xbig(ipoin1)=xsmall(igrou1)
                 if( lntib(ipoin) > 0 ) xbig(ipoin1)=(0.0_rp,0.0_rp)
              else
                 xbig(ipoin1)=(0.0_rp,0.0_rp)
              end if
           enddo
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
              else
                 xbig(ipoin1)=(0.0_rp,0.0_rp)
              end if
           enddo
        end do       
     end if

  endif

end subroutine wvectx
