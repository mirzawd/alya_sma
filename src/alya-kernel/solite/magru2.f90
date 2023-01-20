!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine magru2(ngrou,npopo,nskyl,nbvar,ia,ja,an,askyl,idprecon,wa2)
  !
  ! ASKYL: Factorize group matrix
  !
  use def_kintyp, only               :  ip,rp
  use def_master, only               :  kfl_paral,parre,nparr
  use def_solver, only               :  solve_sol,SOL_NO_PRECOND,SOL_LINELET,&
       &                        SOL_SQUARE,SOL_DIAGONAL,SOL_MATRIX
  use mod_memchk
  use mod_skyline, only   :  lufact
  implicit none
  integer(ip), intent(in)            :: ngrou,npopo,nskyl,nbvar,idprecon
  integer(ip), intent(in)            :: ia(*),ja(*)
  real(rp),    intent(in)            :: an(nbvar,nbvar,*),wa2(*)
  real(rp),    intent(inout), target :: askyl(nskyl)
  integer(ip), pointer               :: iskyl(:),lgrou(:),idiag(:)
  integer(ip)                        :: igrou,jgrou,ipoin,izdom,jpoin,igrou1,jgrou1,idofn,jdofn
  integer(ip)                        :: info,kskyl
  !  real(rp)                           :: lsky(ngrou,ngrou),usky(ngrou,ngrou),asky2(ngrou,ngrou)
  !  real(rp)                           :: asky3(ngrou,ngrou),mu(ngrou),mu2(ngrou) 
!  integer(ip)                        :: iplace,kgrou


  if(ngrou==0) return

  if(kfl_paral/=0) then
     !
     ! Fill in skyline matrix ASKYL 
     !       
     lgrou => solve_sol(1)%lgrou
     iskyl => solve_sol(1)%iskyl
     idiag => solve_sol(1)%idiag

     if(nbvar==1)then

        if(idprecon==SOL_NO_PRECOND)then

            do ipoin=1,npopo
             if(lgrou(ipoin)>0) then
               igrou=lgrou(ipoin)
               do izdom=ia(ipoin),ia(ipoin+1)-1
                  jpoin=ja(izdom)
 
                  if(lgrou(jpoin)>0) then
                     jgrou=lgrou(jpoin)
                     if(igrou<jgrou) then
                         kskyl=iskyl(jgrou+1)-(jgrou-igrou)
                        askyl(kskyl)=askyl(kskyl)+an(1,1,izdom)
                     else     
                        kskyl=idiag(igrou)-(igrou-jgrou)
                        askyl(kskyl)=askyl(kskyl)+an(1,1,izdom)
                     endif
                  end if
               end do
            end if
         end do
       
        else if(idprecon==SOL_DIAGONAL)then 

         do ipoin=1,npopo
             if(lgrou(ipoin)>0) then
               igrou=lgrou(ipoin)
               do izdom=ia(ipoin),ia(ipoin+1)-1
                  jpoin=ja(izdom)
 
                  if(lgrou(jpoin)>0) then
                     jgrou=lgrou(jpoin)
                     if(igrou<jgrou) then
                         kskyl=iskyl(jgrou+1)-(jgrou-igrou)
                        askyl(kskyl)=askyl(kskyl)+wa2(igrou)*an(1,1,izdom)
                     else     
                        kskyl=idiag(igrou)-(igrou-jgrou)
                        askyl(kskyl)=askyl(kskyl)+wa2(igrou)*an(1,1,izdom)
                     endif
                  end if
               end do
            end if
         end do

      else if(idprecon==SOL_SQUARE)then 

         do ipoin=1,npopo
             if(lgrou(ipoin)>0) then
               igrou=lgrou(ipoin)
               do izdom=ia(ipoin),ia(ipoin+1)-1
                  jpoin=ja(izdom)
 
                  if(lgrou(jpoin)>0) then
                     jgrou=lgrou(jpoin)
                     if(igrou<jgrou) then
                         kskyl=iskyl(jgrou+1)-(jgrou-igrou)
                        askyl(kskyl)=askyl(kskyl)+wa2(igrou)*an(1,1,izdom)*wa2(jgrou)
                     else     
                        kskyl=idiag(igrou)-(igrou-jgrou)
                        askyl(kskyl)=askyl(kskyl)+wa2(igrou)*an(1,1,izdom)*wa2(jgrou)
                     endif
                  end if
               end do
            end if
         end do

       else
        call runend('MAGRU2: PRECONDITIONER NOT READY')
       endif   
 
     else

        if(idprecon==SOL_NO_PRECOND)then

        do ipoin=1,npopo
           if(lgrou(ipoin)>0) then
              igrou=lgrou(ipoin)

              do izdom=ia(ipoin),ia(ipoin+1)-1
                 jpoin=ja(izdom)

                 if(lgrou(jpoin)>0) then
                    jgrou=lgrou(jpoin)

                    do idofn=1,nbvar
                       do jdofn=1,nbvar 

                          igrou1=(igrou-1)*nbvar+idofn
                          jgrou1=(jgrou-1)*nbvar+jdofn

                          if(igrou1<jgrou1) then
                             kskyl=iskyl(jgrou1+1)-(jgrou1-igrou1)
                             askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izdom)
                          else     
                             kskyl=idiag(igrou1)-(igrou1-jgrou1)
                             askyl(kskyl)=askyl(kskyl)+an(idofn,jdofn,izdom)
                          endif

                       end do
                    end do


                 end if
              end do
           end if
        end do

        else if(idprecon==SOL_DIAGONAL)then 


         do ipoin=1,npopo
           if(lgrou(ipoin)>0) then
              igrou=lgrou(ipoin)

              do izdom=ia(ipoin),ia(ipoin+1)-1
                 jpoin=ja(izdom)

                 if(lgrou(jpoin)>0) then
                    jgrou=lgrou(jpoin)

                    do idofn=1,nbvar
                       do jdofn=1,nbvar 

                          igrou1=(igrou-1)*nbvar+idofn
                          jgrou1=(jgrou-1)*nbvar+jdofn

                          if(igrou1<jgrou1) then
                             kskyl=iskyl(jgrou1+1)-(jgrou1-igrou1)
                             askyl(kskyl)=askyl(kskyl)+wa2(igrou1)*an(idofn,jdofn,izdom)
                          else     
                             kskyl=idiag(igrou1)-(igrou1-jgrou1)
                             askyl(kskyl)=askyl(kskyl)+wa2(igrou1)*an(idofn,jdofn,izdom)
                          endif

                       end do
                    end do


                 end if
              end do
           end if
        end do
    else if(idprecon==SOL_SQUARE)then 


         do ipoin=1,npopo
           if(lgrou(ipoin)>0) then
              igrou=lgrou(ipoin)

              do izdom=ia(ipoin),ia(ipoin+1)-1
                 jpoin=ja(izdom)

                 if(lgrou(jpoin)>0) then
                    jgrou=lgrou(jpoin)

                    do idofn=1,nbvar
                       do jdofn=1,nbvar 

                          igrou1=(igrou-1)*nbvar+idofn
                          jgrou1=(jgrou-1)*nbvar+jdofn

                          if(igrou1<jgrou1) then
                             kskyl=iskyl(jgrou1+1)-(jgrou1-igrou1)
                             askyl(kskyl)=askyl(kskyl)+wa2(igrou1)*an(idofn,jdofn,izdom)*wa2(jgrou1)
                          else     
                             kskyl=idiag(igrou1)-(igrou1-jgrou1)
                             askyl(kskyl)=askyl(kskyl)+wa2(igrou1)*an(idofn,jdofn,izdom)*wa2(jgrou1)
                          endif

                       end do
                    end do


                 end if
              end do
           end if
        end do    
    else
        call runend('MAGRU2: PRECONDITIONER NOT READY')
       endif   
 
     endif

  end if

  !
  !     Copy askyl in asky2
  !  
  !   do igrou=1,ngrou
  !       do jgrou=1,ngrou
  !           asky2(igrou,jgrou)=0.0d+00
  !       enddo
  !   enddo 

  !   do igrou=1,ngrou
  !       do jgrou=1,ngrou
  !           asky3(igrou,jgrou)=0.0d+00
  !       enddo
  !   enddo
  !    do igrou=1,ngrou
  !!       do jgrou=1,ngrou
  !           lsky(igrou,jgrou)=0.0d+00
  !       enddo
  !   enddo
  !   do igrou=1,ngrou
  !       do jgrou=1,ngrou
  !           usky(igrou,jgrou)=0.0d+00
  !       enddo
  !   enddo


  ! do igrou=1,ngrou
  !       jgrou=igrou-(idiag(igrou)-iskyl(igrou))
  !       do iplace=iskyl(igrou),idiag(igrou)-1
  !           asky2(igrou,jgrou)=askyl(iplace)          
  !           jgrou=jgrou+1
  !       enddo
  !
  !     Diagonal U
  !
  !       asky2(igrou,igrou)=askyl(idiag(igrou)) 

  !         jgrou=igrou-(iskyl(igrou+1)-1-idiag(igrou))
  !       do iplace=idiag(igrou)+1,iskyl(igrou+1)-1
  !           asky2(jgrou,igrou)=askyl(iplace)
  !           jgrou=jgrou+1
  !       enddo
  !  enddo








  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     nparr =  nskyl
     parre => askyl
     call par_operat(3_ip)
  end if
  !
  ! Inverse matrix ASKYL
  !
  if(kfl_paral/=0) then
     call lufact(ngrou*nbvar,nskyl,iskyl,askyl,idiag,info)
     if(info/=0) call runend('MATGRO: ERROR WHILE DOING LU FACTORIZATION')
  end if

  !
  !    DBG FATORIZATION
  !  


  ! do igrou=1,ngrou
  !       jgrou=igrou-(idiag(igrou)-iskyl(igrou))
  !       do iplace=iskyl(igrou),idiag(igrou)-1
  !           lsky(igrou,jgrou)=askyl(iplace)          
  !           jgrou=jgrou+1
  !       enddo
  !
  !     Diagonal  L
  !     
  !       lsky(igrou,igrou)=1.0d+00
  !
  !     Diagonal U
  !
  !       usky(igrou,igrou)=askyl(idiag(igrou)) 

  !         jgrou=igrou-(iskyl(igrou+1)-1-idiag(igrou))
  !       do iplace=idiag(igrou)+1,iskyl(igrou+1)-1
  !           usky(jgrou,igrou)=askyl(iplace)
  !           jgrou=jgrou+1
  !       enddo
  !  enddo

  !  do igrou=1,ngrou
  !     do jgrou=1,ngrou
  !       do kgrou=1,ngrou

  !        asky3(igrou,jgrou)=asky3(igrou,jgrou)+lsky(igrou,kgrou)*usky(kgrou,jgrou)

  !       enddo
  !     enddo
  !  enddo

  ! do igrou=1,ngrou
  !    do jgrou=1,ngrou

  !      if(abs(asky3(igrou,jgrou)-asky2(igrou,jgrou))>1.0d-8)then
  !          write(*,*)asky3(igrou,jgrou),asky2(igrou,jgrou)
  !          call abort 
  !      endif
  !      write(*,*)asky2(igrou,jgrou)
  !    enddo
  ! enddo


  ! do igrou=1,ngrou
  !     mu(igrou)=2.0d+00
  ! enddo

  ! call LU_solve(&                                      ! A'.mu = W^T.r_{-1} 
  !             ngrou,nskyl,&
  !             iskyl,1_ip,askyl,&
  !             mu,ngrou,idiag,info)

  ! do igrou=1,ngrou
  !    do jgrou=1,ngrou
  !      mu2(igrou)=mu2(igrou)+asky3(igrou,jgrou)*mu(jgrou)
  !    enddo
  ! enddo 



end subroutine magru2

