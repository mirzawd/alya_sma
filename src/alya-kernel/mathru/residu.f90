!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine residu(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)
  !------------------------------------------------------------------------
  !****f* mathru/residu
  ! NAME 
  !    residu
  ! DESCRIPTION
  !    Compute the L2, L1, Linf difference (with relaxation) between
  !    two vectors: 
  !    redif = || [r*v1+(1-r)*v2] - v2|| / ||r*v1+(1-r)*v2||  
  ! INPUT
  !    V1(NN1,*)
  !    V2(NN2,*)
  !    residual over the KDIME dimensions, starting from KCOM1 and KCOM2.
  ! USES
  ! USED BY
  !    Modules: *_cvgunk
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_domain
  use def_master
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip), intent(in)  :: nn1,nn2,norm,kcom1,kcom2,kdime
  real(rp),    intent(in)  :: v1(*),v2(*),relax
  real(rp),    intent(out) :: redif
  integer(ip)              :: ii,i1,i2,idime
  real(rp)                 :: va,vo,vd,ra,ro,resin(2)
  real(rp),    target      :: resid(2)

  ra       = relax
  ro       = 1.0_rp-ra
  resid(1) = 0.0_rp
  resid(2) = 0.0_rp

  if(kfl_paral==-1) then

     select case(norm)

     case(0)

        if( ra == 1.0_rp ) then
           do ii = 1,npoin
              i1    = (ii-1)*nn1+kcom1
              i2    = (ii-1)*nn2+kcom2
              resin = 0.0_rp
              do idime = 0,kdime-1
                 va       = v1(i1+idime)
                 vd       = va-v2(i2+idime)
                 resin(1) = resin(1) + vd * vd
                 resin(2) = resin(2) + va * va
              end do
              resin    = sqrt(resin)
              resid(1) = max(resid(1),resin(1))
              resid(2) = max(resid(2),resin(2))
           end do
        else
           do ii = 1,npoin
              i1    = (ii-1)*nn1+kcom1
              i2    = (ii-1)*nn2+kcom2
              resin = 0.0_rp
              do idime = 0,kdime-1
                 vo       = v2(i2+idime)
                 va       = ra*v1(i1+idime)+ro*vo
                 vd       = va-vo
                 resin(1) = resin(1) + vd * vd
                 resin(2) = resin(2) + va * va                 
              end do
              resin    = sqrt(resin)
              resid(1) = max(resid(1),resin(1))
              resid(2) = max(resid(2),resin(2))              
           end do
        end if
        if(resid(2)>zeror) then
           redif = resid(1)/resid(2)
        else
           redif = resid(1)
        end if

     case(1)

        if(ra==1.0_rp) then
           do ii=1,npoin
              i1=(ii-1)*nn1+kcom1
              i2=(ii-1)*nn2+kcom2
              do idime=0,kdime-1
                 va = v1(i1+idime)
                 vd = va-v2(i2+idime)
                 resid(1) = resid(1) + abs(vd)
                 resid(2) = resid(2) + abs(va)
              end do
           end do
        else
           do ii=1,npoin
              i1=(ii-1)*nn1+kcom1
              i2=(ii-1)*nn2+kcom2
              do idime=0,kdime-1
                 vo = v2(i2+idime)
                 va = ra*v1(i1+idime)+ro*vo
                 resid(1) = resid(1) + abs(va-vo)
                 resid(2) = resid(2) + abs(va)
              end do
           end do
        end if
        if(resid(2)>zeror) then
           redif = resid(1)/resid(2)
        else
           redif = resid(1)
        end if
        
     case(2)

        if(ra==1.0_rp) then
           do ii=1,npoin
              i1=(ii-1)*nn1+kcom1
              i2=(ii-1)*nn2+kcom2
              do idime=0,kdime-1
                 va = v1(i1+idime)
                 vd = va-v2(i2+idime)
                 resid(1) = resid(1) + vd*vd
                 resid(2) = resid(2) + va*va
              end do
           end do
        else
           do ii=1,npoin
              i1=(ii-1)*nn1+kcom1
              i2=(ii-1)*nn2+kcom2
              do idime=0,kdime-1
                 vo = v2(i2+idime)
                 va = ra*v1(i1+idime)+ro*vo
                 vd = va-vo
                 resid(1) = resid(1) + vd*vd
                 resid(2) = resid(2) + va*va
              end do
           end do
        end if

        if(resid(2)>zeror) then
           redif = sqrt(resid(1)/resid(2))
        else
           redif = sqrt(resid(1))
        end if

     end select

  else if(kfl_paral>=0) then

     select case(norm)

     case(0)

        do ii = 1,npoin
           i1    = (ii-1)*nn1+kcom1
           i2    = (ii-1)*nn2+kcom2
           resin = 0.0_rp
           do idime = 0,kdime-1
              vo       = v2(i2+idime)
              va       = ra*v1(i1+idime)+ro*vo
              vd       = va-vo
              resin(1) = resin(1) + vd * vd
              resin(2) = resin(2) + va * va                 
           end do
           resin    = sqrt(resin)
           resid(1) = max(resid(1),resin(1))
           resid(2) = max(resid(2),resin(2))              
        end do
        
        call PAR_MAX(2_ip,resid)

        if(resid(2)>zeror) then
           redif = resid(1)/resid(2)
        else
           redif = resid(1)
        end if

     case(1)

        if(kfl_paral>=1) then
           do ii=1,npoi1
              i1=(ii-1)*nn1+kcom1
              i2=(ii-1)*nn2+kcom2
              do idime=0,kdime-1                 
                 vo = v2(i2+idime)
                 va = ra*v1(i1+idime)+ro*vo                              
                 resid(1) = resid(1) + abs(va-vo)
                 resid(2) = resid(2) + abs(va)
              end do
           end do
           do ii=npoi2,npoi3
              i1=(ii-1)*nn1+kcom1
              i2=(ii-1)*nn2+kcom2
              do idime=0,kdime-1   
                 vo = v2(i2+idime)
                 va = ra*v1(i1+idime)+ro*vo                                           
                 resid(1) = resid(1) + abs(va-vo)
                 resid(2) = resid(2) + abs(va)
              end do
           end do
        end if
        call PAR_SUM(2_ip,resid) 

        if(resid(2)>zeror) then
           redif = resid(1)/resid(2)
        else
           redif = resid(1)
        end if

     case(2)

        if(kfl_paral>=1) then
           if(ra==1.0_rp) then
              do ii=1,npoi1
                 i1=(ii-1)*nn1+kcom1
                 i2=(ii-1)*nn2+kcom2
                 do idime=0,kdime-1 
                    va = v1(i1+idime)
                    vd = va-v2(i2+idime)
                    resid(1) = resid(1) + vd*vd
                    resid(2) = resid(2) + va*va
                 end do
              end do
              do ii=npoi2,npoi3
                 i1=(ii-1)*nn1+kcom1
                 i2=(ii-1)*nn2+kcom2
                 do idime=0,kdime-1 
                    va = v1(i1+idime)     
                    vd = va-v2(i2+idime)
                    resid(1) = resid(1) + vd*vd
                    resid(2) = resid(2) + va*va
                 end do
              end do
           else
              do ii=1,npoi1
                 i1=(ii-1)*nn1+kcom1
                 i2=(ii-1)*nn2+kcom2
                 do idime=0,kdime-1                 
                    vo = v2(i2+idime)
                    va = ra*v1(i1+idime)+ro*vo    
                    vd = va-vo
                    resid(1) = resid(1) + vd*vd
                    resid(2) = resid(2) + va*va
                 end do
              end do
              do ii=npoi2,npoi3
                 i1=(ii-1)*nn1+kcom1
                 i2=(ii-1)*nn2+kcom2
                 do idime=0,kdime-1                 
                    vo = v2(i2+idime)
                    va = ra*v1(i1+idime)+ro*vo        
                    vd = va-vo
                    resid(1) = resid(1) + vd*vd
                    resid(2) = resid(2) + va*va
                 end do
              end do
           end if
        end if        
        call PAR_SUM(2_ip,resid) 

        if(resid(2)>zeror) then
           redif = sqrt(resid(1)/resid(2))
        else
           redif = sqrt(resid(1))
        end if

     end select

  end if

end subroutine residu
