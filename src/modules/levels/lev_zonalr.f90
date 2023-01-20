!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_zonalr
  !-----------------------------------------------------------------------
  !****f* Levels/lev_zonalr
  ! NAME 
  !    lev_zonalr
  ! DESCRIPTION
  !    Apply redistanciation only in one zone
  ! USES
  !    
  ! USED BY
  !    lev_redist
  !    lev_redieq
  !***
  !-----------------------------------------------------------------------
  use def_parame, only       :  ip,rp
  use def_master, only       :  fleve,ittim,gisca
  use def_levels, only       :  dista_lev,kfl_zonal_lev,kfl_fixno_lev,tyred_lev
  use def_domain, only       :  coord,npoin
  implicit none
  integer(ip)             :: ipoin,nfrou_lev
  real(rp)                :: abs_y,facto

  nfrou_lev=10   ! for the moment I leave it fixed but it might be included in lev.dat file

  select case(kfl_zonal_lev)

  case(1_ip)

     do ipoin=1,npoin
        if(kfl_fixno_lev(1,ipoin)<1) then
           if(abs(coord(2,ipoin))<3.0_rp) then
              if(fleve(ipoin,1)>=0) then
                 fleve(ipoin,1)=dista_lev(ipoin)
              else
                 fleve(ipoin,1)=-dista_lev(ipoin)
              endif
           endif
        end if
     end do

  case(2_ip)

     do ipoin=1,npoin
        if(kfl_fixno_lev(1,ipoin)<1) then
           abs_y=abs(coord(2,ipoin))
           if(abs_y<2.7_rp) then
              if(fleve(ipoin,1)>=0) then
                 fleve(ipoin,1)=dista_lev(ipoin)
              else
                 fleve(ipoin,1)=-dista_lev(ipoin)
              endif
           else if(abs_y<3.3_rp) then
              facto=(abs_y-2.7_rp)/0.6_rp
              if(fleve(ipoin,1)>=0) then
                 fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
              else
                 fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
              endif
           endif
        end if
     end do

  case(3_ip)

     if(tyred_lev==1) then

        if(mod(ittim,nfrou_lev)==0) then  !do it in all cut nodes at a certain frequency given in nfrou_lev

           do ipoin=1,npoin
              if (kfl_fixno_lev(1,ipoin)<1) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              end if
           end do

        else   !do it only on the cut nodes close to the boat at the frequency determined by nfred_lev that should preferably be=1
           do ipoin=1,npoin
              if (kfl_fixno_lev(1,ipoin)<1) then
                 if(abs(coord(2,ipoin))<3.0_rp) then
                    if(fleve(ipoin,1)>=0) then
                       fleve(ipoin,1)=dista_lev(ipoin)
                    else
                       fleve(ipoin,1)=-dista_lev(ipoin)
                    endif
                 endif
              end if
           end do

        end if


     else  ! tyred

        if(mod(ittim,nfrou_lev)==0) then  !do it in all cut nodes at a certain frequency given in nfrou_lev

           do ipoin=1,npoin
              if ( (kfl_fixno_lev(1,ipoin)<1).and.(gisca(ipoin)>0) ) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              end if
           end do

        else   !do it only on the cut nodes close to the boat at the frequency determined by nfred_lev that should preferably be=1
           do ipoin=1,npoin
              if( (kfl_fixno_lev(1,ipoin)<1).and.(gisca(ipoin)>0) ) then
                 if(abs(coord(2,ipoin))<3.0_rp) then
                    if(fleve(ipoin,1)>=0) then
                       fleve(ipoin,1)=dista_lev(ipoin)
                    else
                       fleve(ipoin,1)=-dista_lev(ipoin)
                    endif
                 endif
              end if
           end do

        end if

     end if! tyred


  case(4_ip)

     if(tyred_lev==1) then

        if(mod(ittim,nfrou_lev)==0) then  !do it in all cut nodes at a certain frequency given in nfrou_lev

           do ipoin=1,npoin
              if (kfl_fixno_lev(1,ipoin)<1) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              end if
           end do

        else   !do it only on the cut nodes close to the boat at the frequency determined by nfred_lev that should preferably be=1
           do ipoin=1,npoin
              if (kfl_fixno_lev(1,ipoin)<1) then
                 abs_y=abs(coord(2,ipoin))
                 if(abs_y<2.7_rp) then
                    if(fleve(ipoin,1)>=0) then
                       fleve(ipoin,1)=dista_lev(ipoin)
                    else
                       fleve(ipoin,1)=-dista_lev(ipoin)
                    endif
                 else if(abs_y<3.3_rp) then
                    facto=(abs_y-2.7_rp)/0.6_rp
                    if(fleve(ipoin,1)>=0) then
                       fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
                    else
                       fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
                    endif
                 endif
              end if
           end do

        end if


     else  ! tyred

        if(mod(ittim,nfrou_lev)==0) then  !do it in all cut nodes at a certain frequency given in nfrou_lev

           do ipoin=1,npoin
              if( (kfl_fixno_lev(1,ipoin)<1).and.(gisca(ipoin)>0) ) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              end if
           end do

        else   !do it only on the cut nodes close to the boat at the frequency determined by nfred_lev that should preferably be=1

           do ipoin=1,npoin
              if( (kfl_fixno_lev(1,ipoin)<1).and.(gisca(ipoin)>0) ) then
                 abs_y=abs(coord(2,ipoin))
                 if(abs_y<2.7_rp) then
                    if(fleve(ipoin,1)>=0) then
                       fleve(ipoin,1)=dista_lev(ipoin)
                    else
                       fleve(ipoin,1)=-dista_lev(ipoin)
                    endif
                 else if(abs_y<3.3_rp) then
                    facto=(abs_y-2.7_rp)/0.6_rp
                    if(fleve(ipoin,1)>=0) then
                       fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
                    else
                       fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
                    endif
                 endif
              end if
           end do

        end if

     end if! tyred

  case(5_ip)  ! update only patially outside but at each time step

     facto=0.1_rp

     if(tyred_lev==1) then


        do ipoin=1,npoin

           if(kfl_fixno_lev(1,ipoin)<1) then
              if(abs(coord(2,ipoin))<3.0_rp) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 end if
              else
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1) = facto*dista_lev(ipoin) + (1.0_rp-facto)*fleve(ipoin,1)
                 else
                    fleve(ipoin,1)= -facto*dista_lev(ipoin) + (1.0_rp-facto)*fleve(ipoin,1)
                 end if
              end if
           end if
        end do

     else  ! tyred

        do ipoin=1,npoin

           if( (kfl_fixno_lev(1,ipoin)<1).and.(gisca(ipoin)>0) ) then
              if(abs(coord(2,ipoin))<3.0_rp) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 end if
              else
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1) = facto*dista_lev(ipoin) + (1.0_rp-facto)*fleve(ipoin,1)
                 else
                    fleve(ipoin,1)= -facto*dista_lev(ipoin) + (1.0_rp-facto)*fleve(ipoin,1)
                 end if
              end if
           end if
        end do

     end if! tyred


  case(6_ip)   ! the same idea as 5 but with a smooth transition between 2.7 and 3.3 

     facto=0.1_rp

     if(tyred_lev==1) then

        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)<1) then
              abs_y=abs(coord(2,ipoin))
              if(abs_y<2.7_rp) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              else if(abs_y<3.3_rp) then
                 facto =  (101.0_rp/20.0_rp) - (1.5_rp*abs_y)    ! 1 at 2,7 and 0,1 at 3.3
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
                 endif
              else                       !    abs_y>3.3_rp
                 facto =  0.1_rp
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
                 endif
              endif
           end if
        end do


     else  ! tyred

        do ipoin=1,npoin
           if( (kfl_fixno_lev(1,ipoin)<1).and.(gisca(ipoin)>0) ) then
              abs_y=abs(coord(2,ipoin))
              if(abs_y<2.7_rp) then
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=-dista_lev(ipoin)
                 endif
              else if(abs_y<3.3_rp) then
                 facto =  (101.0_rp/20.0_rp) - (1.5_rp*abs_y)    ! 1 at 2,7 and 0,1 at 3.3
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
                 endif
              else                       !    abs_y>3.3_rp
                 facto =  0.1_rp
                 if(fleve(ipoin,1)>=0) then
                    fleve(ipoin,1)=facto*fleve(ipoin,1) + (1.0_rp-facto)*dista_lev(ipoin)
                 else
                    fleve(ipoin,1)=facto*fleve(ipoin,1) - (1.0_rp-facto)*dista_lev(ipoin)
                 endif
              endif
           end if
        end do

     end if! tyred

  end select


end subroutine lev_zonalr
