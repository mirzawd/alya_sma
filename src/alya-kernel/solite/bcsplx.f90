!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine bcsplx(nbnodes,nbvar,an,ja,ia,xx,yy) 

  !----------------------------------------------------------------------------------   
  ! Sources/kernel/solite/bcsplx.f90
  ! NAME 
  !     bcsplx
  ! DESCRIPTION
  !     This routine multiplies a non-symmetric complex matrix stored in BCSR format
  !     with a complex vector:
  !                         YY = A XX 
  ! INPUTPUT ARGUMENTS
  !    NBNODES .... Number of nodes
  !    NBVAR ...... Number of variables in each node
  !    AN ......... Complex matrix
  !    JA ......... Compressed Sparse format: index vector for column numbers
  !    IA ......... Compressed Sparse format: index vector for beginning of a row block
  !    XX ......... Complex vector
  ! OUTPUT ARGUMENTS
  !    YY ......... Complex result vector
  ! USES
  ! USED BY
  !-------------------------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER,NPOIN_TYPE
  !Declaration statements
  implicit none

  !Dummy arguments
  integer(ip), intent(in)          :: nbnodes,nbvar
  complex(rp), intent(in)          :: an(*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  complex(rp), intent(in)          :: xx(*)
  complex(rp), intent(out), target :: yy(*)

  !Local variables
  integer(ip)                      :: ii,ih,jj,jh,kk,kh,ll,col,colh
  complex(rp)                      :: temp

  !------------------------------------------------------------
  ! y = Ax
  !------------------------------------------------------------
  if (INOTMASTER) then
     if (nbvar == 4_ip) then
        ! !$omp parallel do                                &
        ! !$omp            default(shared)                 &
        ! !$omp            private(ii,ih,jj,jh,col,colh)   &
        ! !$omp            schedule(static)                
        do ii = 1,nbnodes
           ih = 4_ip * (ii - 1_ip)
           yy(ih+1) = (0.0_rp,0.0_rp)
           yy(ih+2) = (0.0_rp,0.0_rp)
           yy(ih+3) = (0.0_rp,0.0_rp)
           yy(ih+4) = (0.0_rp,0.0_rp)
           do jj = ia(ii),ia(ii+1)-1_ip
              jh = 16_ip * (jj - 1_ip)
              col = ja(jj)
              colh = 4_ip * (col - 1_ip)
              yy(ih+1) = yy(ih+1) + (an(jh+1) * xx(colh+1) &
                   +  an(jh+2) * xx(colh+2) &
                   +  an(jh+3) * xx(colh+3) &
                   +  an(jh+4) * xx(colh+4))
              yy(ih+2) = yy(ih+2) + (an(jh+5) * xx(colh+1) &
                   +  an(jh+6) * xx(colh+2) &
                   +  an(jh+7) * xx(colh+3) &
                   +  an(jh+8) * xx(colh+4))
              yy(ih+3) = yy(ih+3) + (an(jh+9)  * xx(colh+1) &
                   +  an(jh+10) * xx(colh+2) &
                   +  an(jh+11) * xx(colh+3) &
                   +  an(jh+12) * xx(colh+4))
              yy(ih+4) = yy(ih+4) + (an(jh+13) * xx(colh+1) &
                   +  an(jh+14) * xx(colh+2) &
                   +  an(jh+15) * xx(colh+3) &
                   +  an(jh+16) * xx(colh+4))                                   
           end do
        end do
        ! !$omp end parallel do
     else
        ! !$omp parallel                                                      &
        ! !$omp            default(shared)                                    &
        ! !$omp            private(temp,ii,ih,kk,kh,jj,jh,col,colh,ll)        
        temp = (0.0_rp,0.0_rp)
        ! !$omp do schedule(static)                         
        do ii = 1,nbnodes
           ih = nbvar * (ii - 1_ip)
           do kk = 1,nbvar
              kh = nbvar * (kk - 1_ip)
              yy(ih+kk) = (0.0_rp,0.0_rp)
              do jj = ia(ii),ia(ii+1)-1_ip
                 jh = nbvar * nbvar * (jj - 1_ip)
                 col = ja(jj)
                 colh = nbvar * (col - 1_ip)
                 do ll = 1,nbvar
                    temp = temp + an(jh+kh+ll) * xx(colh+ll)
                 end do
                 yy(ih+kk) = yy(ih+kk) + temp
              end do
           end do
        end do
        ! !$omp end do
        ! !$omp end parallel
     end if
     call pararx('SLX',NPOIN_TYPE,nbnodes*nbvar,yy)
     !call rhsmox(nbvar,yy)
  end if

end subroutine bcsplx

subroutine bcsplx2(nbnodes,nbvar,an,ja,ia,xx,yy) 

  !----------------------------------------------------------------------------------   
  ! Sources/kernel/solite/bcsplx.f90
  ! NAME 
  !     bcsplx2
  ! DESCRIPTION
  !     This routine multiplies a non-symmetric complex matrix stored in BCSR format
  !     with a complex vector:
  !                         YY = A XX - YY
  ! INPUTPUT ARGUMENTS
  !    NBNODES .... Number of nodes
  !    NBVAR ...... Number of variables in each node
  !    AN ......... Complex matrix
  !    JA ......... Compressed Sparse format: index vector for column numbers
  !    IA ......... Compressed Sparse format: index vector for beginning of a row block
  !    XX ......... Complex vector
  ! OUTPUT ARGUMENTS
  !    YY ......... Complex result vector
  ! USES
  ! USED BY
  !-------------------------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER,NPOIN_TYPE
  !Declaration statements
  implicit none

  !Dummy arguments
  integer(ip), intent(in)            :: nbnodes,nbvar
  complex(rp), intent(in)            :: an(*)
  integer(ip), intent(in)            :: ja(*),ia(*)
  complex(rp), intent(in)            :: xx(*)
  complex(rp), intent(inout), target :: yy(*)

  !Local variables
  integer(ip)                      :: ii,ih,jj,jh,col,colh
!  integer(ip)                      :: kk,kh,ll
  complex(rp)                      :: yy1,yy2,yy3,yy4
!  complex(rp)                      :: temp

  !------------------------------------------------------------
  ! y = Ax - y
  !------------------------------------------------------------
  if (INOTMASTER) then
     if (nbvar == 4_ip) then
        !$OMP parallel do                                                &
        !$OMP            default(shared)                                 &
        !$OMP            private(ii,ih,jj,jh,col,colh,yy1,yy2,yy3,yy4)   &
        !$OMP            schedule(guided)          
        do ii = 1,nbnodes
           ih = 4_ip * (ii - 1_ip)
           !yy(ih+1) = (0.0_rp,0.0_rp)
           !yy(ih+2) = (0.0_rp,0.0_rp)
           !yy(ih+3) = (0.0_rp,0.0_rp)
           !yy(ih+4) = (0.0_rp,0.0_rp)
           yy1 = (0.0_rp,0.0_rp)
           yy2 = (0.0_rp,0.0_rp)
           yy3 = (0.0_rp,0.0_rp)
           yy4 = (0.0_rp,0.0_rp)
           do jj = ia(ii),ia(ii+1)-1_ip
              jh = 16_ip * (jj - 1_ip)
              col = ja(jj)
              colh = 4_ip * (col - 1_ip)
              yy1 = yy1 + (an(jh+1) * xx(colh+1) &
                   +  an(jh+2) * xx(colh+2) &
                   +  an(jh+3) * xx(colh+3) &
                   +  an(jh+4) * xx(colh+4))
              yy2 = yy2 + (an(jh+5) * xx(colh+1) &
                   +  an(jh+6) * xx(colh+2) &
                   +  an(jh+7) * xx(colh+3) &
                   +  an(jh+8) * xx(colh+4))
              yy3 = yy3 + (an(jh+9)  * xx(colh+1) &
                   +  an(jh+10) * xx(colh+2) &
                   +  an(jh+11) * xx(colh+3) &
                   +  an(jh+12) * xx(colh+4))
              yy4 = yy4 + (an(jh+13) * xx(colh+1) &
                   +  an(jh+14) * xx(colh+2) &
                   +  an(jh+15) * xx(colh+3) &
                   +  an(jh+16) * xx(colh+4))                                 
           end do
           yy(ih+1)=yy1 - yy(ih+1)
           yy(ih+2)=yy2 - yy(ih+2)
           yy(ih+3)=yy3 - yy(ih+3)
           yy(ih+4)=yy4 - yy(ih+4)
        end do
        !$OMP end parallel do
     else
        print *,'ERROR: Not implemented'
        !		  temp = (0.0_rp,0.0_rp)
        !		  do ii = 1,nbnodes
        !  		  ih = nbvar * (ii - 1_ip)
        !			  do kk = 1,nbvar
        !				  kh = nbvar * (kk - 1_ip)
        !				  yy(ih+kk) = (0.0_rp,0.0_rp)
        !				  do jj = ia(ii),ia(ii+1)-1_ip
        !					  jh = nbvar * nbvar * (jj - 1_ip)
        !					  col = ja(jj)
        !				          colh = nbvar * (col - 1_ip)
        !					  do ll = 1,nbvar
        !						  temp = temp + an(jh+kh+ll) * xx(colh+ll)
        !					  end do
        !					  yy(ih+kk) = yy(ih+kk) + temp
        !				  end do
        !			  end do
        !		  end do
     end if
     call pararx('SLX',NPOIN_TYPE,nbnodes*nbvar,yy)
     !call rhsmox(nbvar,yy)
     !call rhsmod(nbvar,yy)
  end if

end subroutine bcsplx2

subroutine bcsplx2opt(nbnodes,nbvar,an,ja,ia,xx,yy) 

  !----------------------------------------------------------------------------------   
  ! Sources/kernel/solite/bcsplx.f90
  ! NAME 
  !     bcsplx2
  ! DESCRIPTION
  !     This routine multiplies a non-symmetric complex matrix stored in BCSR format
  !     with a complex vector:
  !                         YY = A XX - YY
  ! INPUTPUT ARGUMENTS
  !    NBNODES .... Number of nodes
  !    NBVAR ...... Number of variables in each node
  !    AN ......... Complex matrix
  !    JA ......... Compressed Sparse format: index vector for column numbers
  !    IA ......... Compressed Sparse format: index vector for beginning of a row block
  !    XX ......... Complex vector
  ! OUTPUT ARGUMENTS
  !    YY ......... Complex result vector
  ! USES
  ! USED BY
  !-------------------------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER,NPOIN_TYPE
 !Declaration statements
  implicit none

 !Dummy arguments
  integer(ip), intent(in)          :: nbnodes,nbvar
  complex(rp), intent(in)          :: an(*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  complex(rp), intent(in)          :: xx(*)
  complex(rp), intent(out), target :: yy(*)

 !Local variables
  integer(ip)                      :: ii,ih,jj,jh,col,colh
!  integer(ip)                      :: kk,kh,ll
  integer(ip)                      :: jjstart, jjend
!  integer(ip)                      :: sizean
  complex(rp)                      :: yy1,yy2,yy3,yy4
!  complex(rp)                      :: temp
  complex(rp)                      :: xx1,xx2,xx3,xx4

  real(rp) :: normL2

 !------------------------------------------------------------
 ! y = Ax - y
 !------------------------------------------------------------
  if (INOTMASTER) then
     if (nbvar == 4_ip) then
!!$OMP parallel do                                                                &
!!$OMP            default(shared)                                                 &
!!$OMP            private(ii,ih,jj,jh,col,colh,yy1,yy2,yy3,yy4,xx1,xx2,xx3,xx4)   &
!!$OMP            schedule(guided)          

        !sizean=size(an)
        !normL2 = dot_product(an(1:sizean),an(1:sizean))
        !if(normL2>0) then


        do ii = 1,nbnodes
           ih = 4_ip * (ii - 1_ip)
           !yy(ih+1) = (0.0_rp,0.0_rp)
           !yy(ih+2) = (0.0_rp,0.0_rp)
           !yy(ih+3) = (0.0_rp,0.0_rp)
           !yy(ih+4) = (0.0_rp,0.0_rp)
           yy1 = (0.0_rp,0.0_rp)
           yy2 = (0.0_rp,0.0_rp)
           yy3 = (0.0_rp,0.0_rp)
           yy4 = (0.0_rp,0.0_rp)
           
           jjstart = ia(ii)
           jjend = ia(ii+1)-1_ip


           ! This is not the norm
           normL2 = real(dot_product(an((16*jjstart-15):(16*jjend)),an((16*jjstart-15):(16*jjend))),rp)

           if(normL2>0.0_rp) then

           !do jj = ia(ii),ia(ii+1)-1_ip
           do jj = jjstart,jjend
              jh = 16_ip * (jj - 1_ip)
              col = ja(jj)
              colh = 4_ip * (col - 1_ip)

              xx1 = xx(colh+1)  
              xx2 = xx(colh+2)  
              xx3 = xx(colh+3)  
              xx4 = xx(colh+4)  

              yy1 = yy1     + (an(jh+1) * xx1 &
                         +  an(jh+2) * xx2 &
                         +  an(jh+3) * xx3 &
                         +  an(jh+4) * xx4)
              yy2 = yy2     + (an(jh+5) * xx1 &
                         +  an(jh+6) * xx2 &
                         +  an(jh+7) * xx3 &
                         +  an(jh+8) * xx4)
              yy3 = yy3     + (an(jh+9)  * xx1 &
                         +  an(jh+10) * xx2 &
                         +  an(jh+11) * xx3 &
                         +  an(jh+12) * xx4)
              yy4 = yy4     + (an(jh+13) * xx1 &
                         +  an(jh+14) * xx2 &
                         +  an(jh+15) * xx3 &
                         +  an(jh+16) * xx4)                                 
           end do

           end if

           yy(ih+1)=yy1 - yy(ih+1)
           yy(ih+2)=yy2 - yy(ih+2)
           yy(ih+3)=yy3 - yy(ih+3)
           yy(ih+4)=yy4 - yy(ih+4)
        end do
!!$OMP end parallel do
        !else
        !   yy(:) = -yy(:)
        !end if



     else
                 print *,'ERROR: Not implemented'
!        	  temp = (0.0_rp,0.0_rp)
!        	  do ii = 1,nbnodes
!        	  ih = nbvar * (ii - 1_ip)
!        		  do kk = 1,nbvar
!        			  kh = nbvar * (kk - 1_ip)
!        			  yy(ih+kk) = (0.0_rp,0.0_rp)
!        			  do jj = ia(ii),ia(ii+1)-1_ip
!        				  jh = nbvar * nbvar * (jj - 1_ip)
!        				  col = ja(jj)
!        			          colh = nbvar * (col - 1_ip)
!        				  do ll = 1,nbvar
!        					  temp = temp + an(jh+kh+ll) * xx(colh+ll)
!        				  end do
!        				  yy(ih+kk) = yy(ih+kk) + temp
!        			  end do
!        		  end do
!        	  end do
     end if
     call pararx('SLX',NPOIN_TYPE,nbnodes*nbvar,yy)
     !call rhsmox(nbvar,yy)
     !call rhsmod(nbvar,yy)
   end if

end subroutine bcsplx2opt



