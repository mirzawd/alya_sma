!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_funcre(&
     itask,idofn,ndime,dinew,bvess,kfixn,palaw,npara,ifuge,kfixi,rtini,rtifi,timev)
  !------------------------------------------------------------------------
  !
  ! This subroutine computes transient boundary conditions
  !
  ! CAVEAT: THIS IS DIFFERENT THAN KERNEL'S FUNCRE!!! 
  !
  !------------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip), intent(in)  :: npara,ifuge,kfixi,ndime,itask,idofn
  real(rp),    intent(in)  :: palaw(20),timev
  integer(ip)              :: iftbc,kfixn(ndime),idime
!  integer(ip)              :: ipara
  real(rp)                 :: timea,zerom
!  real(rp)                 :: timeb,funca,funcb,timec
  real(rp)                 :: rtini,rtifi,tirel,tista,tiend,tirep
!  real(rp)                 :: timei,timef
  real(rp)                 :: dinew,bvess(ndime),xx_funcre,bvref(ndime)



  zerom=epsilon(1.0_rp)

  tista= palaw(11)
  tiend= palaw(12)
  tirep= palaw(13)
  bvref =palaw(15:15+ndime-1)


  if (itask == 0_ip) then                    ! initial fixities
     
     ! nada
     
  else if (itask == 1_ip) then               ! check only fixities 

     call runend('SLD_FUNCRE: TIME VARYING DISPLACEMENT MUST BE DONE NOW THROUGH A DATA FILE!!')

     if ((timev+zerom) > rtini) then   
        if (timev < rtifi) then                    ! within the time lapse
           if (kfixi == 1 .or. kfixi==2 .or. kfixi==3 .or. kfixi==4) then 
              do idime=1,ndime
                 if (kfixn(idime) > 0) then
                    bvess(idime)= bvref(idime)
                 end if
              end do

              if (kfixi == 3) then    ! impose to all dimensions, even when fixity code is 0
                 do idime=1,ndime
                    if (kfixn(idime) == 0) kfixn(idime) = 1                    !   release the fixed fixities
                 end do                 
              end if

!              kfixn(1:ndime)= 1                    !   displacement is reference value
!              bvess(1:ndime)= bvref(1:ndime)
           end if
        else                                       ! outside the time lapse
           bvess(1:ndime)= bvref(1:ndime)
           if (kfixi == 1) then                    ! 
!              kfixn(1:ndime)= 1                    !   don't move
           else if (kfixi == 2) then               
!              kfixn(1:ndime)= 1                    !   don't move
           else if (kfixi == 3) then               
!              kfixn(1:ndime)= 1                    !   don't move
           else if (kfixi == 4) then               
              do idime=1,ndime
                 if (kfixn(idime) == 1 .or. kfixn(idime)==2) kfixn(idime) = 0                    !   release the fixed fixities
              end do
           end if
        end if
     end if

  else if (itask == 2_ip) then

     call runend('SLD_FUNCRE: TIME VARYING DISPLACEMENT MUST BE DONE NOW THROUGH A DATA FILE!!')

     dinew=bvess(idofn)
     iftbc= 0
     if ((timev+zerom) > rtini) then   
        if (timev < rtifi) then                    ! within the time lapse
           tirel= timev - rtini
           iftbc= 1
        else
           tirel= rtifi - rtini
           iftbc= 1
        end if
     end if
     
     if (iftbc == 0) return
     
     xx_funcre= 1.0_rp
     
     if(ifuge==0) then
        !
        ! No time dependence 
        !
        xx_funcre=1.0_rp
        
     else if(ifuge==1) then
        !
        ! Polynomial evolution
        !
        xx_funcre= palaw(1) &
             + palaw(2) * tirel &
             + palaw(3) * tirel * tirel &
             + palaw(4) * tirel * tirel * tirel

     else if(ifuge==2) then
        !
        ! Periodic evolution
        !
        if(palaw(1)-zerom<=timev.and.timev<=palaw(2)+zerom) then 
           xx_funcre=palaw(3)*cos(palaw(4)*timev+palaw(5))+palaw(6)
        else if (timev>palaw(2)+zerom) then
           timea=palaw(2)
           xx_funcre=palaw(3)*cos(palaw(4)*timea+palaw(5))+palaw(6)
        else if (timev<palaw(1)-zerom) then
           timea=palaw(1)
           xx_funcre=palaw(3)*cos(palaw(4)*timea+palaw(5))+palaw(6)
        end if
        
     else if(ifuge==3) then
        !
        ! Discrete evolution
        !

        call runend('SLD_FUNCRE: DEPRECATED OPTION!!!!')


!        timei=palaw(1)
!        timef=palaw((npara/2-1)*2+1)
!        
!        if(timev<=timei) then
!           xx_funcre=palaw(2)
!        else
!           if(timev>=timef) then            ! Look for the time inside the period
!              timec=timev
!              do while(timec>timef)
!                 timec=timec-(timef-timei)
!              end do
!           else
!              timec=timev
!           end if
!           ipara=0
!           do while(ipara<npara/2)
!              ipara=ipara+1
!              if(timec<palaw((ipara-1)*2+1)) then
!                 timea=palaw((ipara-1)*2+1)
!                 funca=palaw((ipara-1)*2+2)
!                 timeb=palaw((ipara-2)*2+1)
!                 funcb=palaw((ipara-2)*2+2)
!                 xx_funcre=(funcb-funca)/(timeb-timea)*(timec-timea)+funca
!                 ipara=npara/2
!              end if
!           end do
!        end if
        
     else if(ifuge==4) then
        !
        ! Special function to change boundary values
        !
        xx_funcre=palaw(2)
        
     else if(ifuge==5) then
        !
        ! Marek Prymon's function
        !
        
        
     end if

     !
     !  Compute new value
     !     

     dinew= bvess(idofn)*xx_funcre

     
  else if (itask == 10_ip) then                    ! update time lapses

     call runend('SLD_FUNCRE: TIME VARYING DISPLACEMENT MUST BE DONE NOW THROUGH A DATA FILE!!')

     
     if ((timev+zerom) > rtini) then   
        if (timev < rtifi) then                    ! within the time lapse
           ! ok, do nothing
        else                                       ! outside the time lapse
           rtini= rtifi + tirep
           rtifi= rtini + (tiend - tista)        
        end if
     end if
     
  end if
  
end subroutine sld_funcre
