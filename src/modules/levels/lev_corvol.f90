!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_corvol
  !-----------------------------------------------------------------------
  !****f* Levels/lev_corvol
  ! NAME 
  !    lev_corvol
  ! DESCRIPTION
  !    Volume correction
  !                                    
  ! USES
  !    
  ! USED BY
  !    lev_updunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use def_solver
  use mod_messages, only : livinf
  implicit none

  integer(ip) :: itcor,ipoin
  real(rp)    :: corvol, errvol, ervolr 
  real(rp)    :: eps, ervold

  call livinf(59_ip,'VOLUME CORRECTION',0_ip)

  eps=1e-06_rp
  corvol=0.0_rp
  itcor=0_ip
  call lev_calvol()
  errvol=volit_lev-volrf_lev
  ervolr=errvol/(volrf_lev+zeror)

  do while((abs(ervolr)>1.0e-08_rp).and.(itcor<10_ip))
     itcor=itcor+1_ip
     if(INOTMASTER) then
        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)+eps
        enddo
     endif
     call lev_calvol()
     ervold=(volit_lev-volrf_lev-errvol)/eps

     if(ervold/=0.0_rp) then 

        if(INOTMASTER) then
           do ipoin=1,npoin
              if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)-eps-corvol
           enddo
        endif

        corvol=corvol-errvol/ervold
        if(INOTMASTER) then
           do ipoin=1,npoin
              if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)+corvol
           enddo
        endif
        call lev_calvol()
        errvol=volit_lev-volrf_lev
        ervolr=errvol/volrf_lev
     else
        call runend('LEV_CORVOL: VOLUME CORRECTION FAILED')

     endif
  enddo


end subroutine lev_corvol


subroutine lev_corvolhhh
  !-----------------------------------------------------------------------
  !****f* Levels/lev_corvol
  ! NAME 
  !    lev_corvol
  ! DESCRIPTION
  !    Volume correction  - my implementation similar to what I had in Zephyr - being tested
  !                                    
  ! USES
  !    
  ! USED BY
  !    lev_updunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use def_solver
  use mod_messages, only : livinf
  implicit none

  integer(ip)   :: ipoin
  real(rp)      :: dells(3),resik(2),toler,fivo0,auxde,denom
  integer(ip)   :: miter,kiter

  call livinf(59_ip,'lev_corvol_hhh:VOLUME CORRECTION',0_ip)

  toler=1.0e-5_rp
  miter=20

  dells(3)=0.0_rp
  call lev_calvol()
  resik(2)=volrf_lev-volit_lev
  fivo0=volit_lev   ! for use at the stopping criteria
  write(735,*)'resik(2),volrf_lev,fivo0 ',resik(2),volrf_lev,fivo0 

  dells(2)=1.0e-6_rp   ! buscar un mejor iniguess , ej: paper de orlando
  dells(2)=sign(dells(2),resik(2))

  if(INOTMASTER) then
     do ipoin=1,npoin
        if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)+dells(2)
     enddo
  endif

  call lev_calvol()
  resik(1)=volrf_lev-volit_lev

  iterat: do kiter=1,miter
    denom = resik(1)-resik(2)  ! to avoid divide by zero - perhaps a better solution can be found
    if(abs(denom)>1.0e-15_rp) then
        dells(1) = dells(2) - ( resik(1)* (dells(2)-dells(3)) / denom )
    else
        dells(1) = 0.0_rp 
    end if
    resik(2)=resik(1)
    dells(3)=dells(2)
    dells(2)=dells(1)

    auxde=dells(1)-dells(3)

    if(INOTMASTER) then
       do ipoin=1,npoin
          if(kfl_fixno_lev(1,ipoin)<1) fleve(ipoin,1)=fleve(ipoin,1)+auxde
       enddo
    endif

    call lev_calvol()
    resik(1)=volrf_lev-volit_lev
    write(*,*)'kiter,resik(1),dells(1),volit_lev',kiter,resik(1),dells(1),volit_lev
    write(738,*)'kiter,resik(1),dells(1),volit_lev',kiter,resik(1),dells(1),volit_lev
    if (abs(resik(1))<toler*fivo0) exit iterat
  end do iterat

end subroutine lev_corvolhhh


