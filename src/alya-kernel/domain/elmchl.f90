!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmchl(&
     tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,porde,&
     hnatu,kfl_advec,kfl_ellen)
  !-----------------------------------------------------------------------
  !****f* Domain/elmchl
  ! NAME
  !   elmchl
  ! DESCRIPTION
  !   This routine computes the characteristic element lengths CHALE 
  !   according to a given strategy. CHALE is divided by two for
  !   quadratic elements:
  !   KFL_ELLEN = 0 ... CHALE(1) = Minimum element length
  !                 ... CHALE(2) = Minimum element length
  !   KFL_ELLEN = 1 ... CHALE(1) = Maximum element length
  !                 ... CHALE(2) = Maximum element length
  !   KFL_ELLEN = 2 ... CHALE(1) = Average element length
  !                 ... CHALE(2) = Average element length
  !   KFL_ELLEN = 3 ... IF KFL_ADVEC = 1:
  !                     CHALE(1) = Flow direction
  !                     CHALE(2) = Flow direction
  !                     ELSE IF KFL_ADVEC =0:
  !                     CHALE(1) = Minimum element length
  !                     CHALE(2) = Minimum element length
  !   KFL_ELLEN = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
  !                 ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
  !   KFL_ELLEN = 5 ... CHALE(1) = Length in flow direction
  !                 ... CHALE(2) = Minimum element kength
  !
  !   KFL_ELLEN = 7 ... CHALE(1) = sqrt(area) or (volum)**(1/3) 
  !                 ... CHALE(2) = sqrt(area) or (volum)**(1/3) 
  !   KFL_ELLEN = 8 ... CHALE(1) = volum / largest face area
  !                 ... CHALE(2) = volum / largest face area
  !
  ! OUTPUT
  !   CHALE
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     : ip,rp
  use def_domain, only     : ndime
  use def_domain, only     : mnode
  use def_domain, only     : elmar
  use def_domain, only     : ltopo
  
  implicit none
  
  external                 :: mbvab1
  external                 :: velchl
  external                 :: elmder
  external                 :: elmafa
  
  integer(ip), intent(in)  :: pnode,porde,kfl_ellen
  integer(ip), intent(in)  :: kfl_advec
  integer(ip), intent(in)  :: pelty
  real(rp),    intent(in)  :: hnatu
  real(rp),    intent(out) :: chale(2)
  real(rp),    intent(in)  :: tragl(ndime,ndime),hleng(ndime)
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(out) :: chave(ndime,2)
  integer(ip)              :: idime,inode
  real(rp)                 :: elno1,elno2
  real(rp)                 :: volum,hleng_aux
  real(rp)                 :: xjaci(9),xjacm(9),gpdet
  real(rp)                 :: gpcar(ndime,mnode)
  real(rp)                 :: aface(3)
  
  if(kfl_ellen==0) then 
     !
     ! Minimum element length
     !
     chale(1)=hleng(ndime) 
     chale(2)=chale(1)

  else if(kfl_ellen==1) then   
     !
     ! Maximum element length
     !     
     chale(1)=hleng(1) 
     chale(2)=chale(1)

  else if(kfl_ellen==2) then 
     !
     ! Average length
     !
     chale(1)=0.0_rp
     do idime=1,ndime
        chale(1)=chale(1)+hleng(idime)
     end do
     chale(1)=chale(1)/real(ndime,rp) 
     chale(2)=chale(1)

  else if(kfl_ellen==3) then 
     !
     ! Length in flow direction
     !
     if( kfl_advec/=0 ) then 
        !
        ! Characteristic element velocity (average)
        !
        chave=0.0_rp
        do idime=1,ndime
           do inode=1,pnode
              chave(idime,1)=chave(idime,1)+elvel(idime,inode)
           end do
           chave(idime,1)=chave(idime,1)/real(pnode,rp)
        end do
        !
        ! Characteristic element length u^l = J^(-t) u^g
        !
        call mbvab1(chave(1,2),tragl,chave(1,1),ndime,ndime,elno2,elno1)
        if(elno2>1.0e-16_rp.and.elno1>1.0e-16_rp) then
           chale(1)=hnatu*elno1/elno2
        else
           chale(1)=hleng(ndime)
        end if
        chale(2)=chale(1)
        chale(2)=hleng(ndime)
        if (ndime ==3 ) then
           chale(2)=(hleng(ndime)*hleng(2)*hleng(1))**(1.0_rp/3.0_rp)
        else if (ndime==2) then
           chale(2)=sqrt(hleng(2)*hleng(1))
        end if
     else
        chale(1)=hleng(ndime)       
        chale(2)=chale(1)
     end if

  else if(kfl_ellen==4) then 
     !
     ! sqrt(hmin*hmax)
     !
     chale(1)=sqrt(hleng(1)*hleng(ndime))     
     chale(2)=chale(1)

  else if(kfl_ellen==5) then 
     !
     ! Along velocity direction
     !
     call velchl(pnode,elcod,elvel,chale,hleng)

  else if(kfl_ellen==6) then 
     !
     ! Mixed element length - hmin for tau1, hmax for tau2 - here we only obtain the values for tau1 - tau2 directly in nsi_elmsgs
     !
     chale(1)=hleng(ndime) 
     chale(2)=chale(1)
     
  else if(kfl_ellen==7) then
     !
     ! From length/area/volume for bar, quadrilateral and hexahedral elements
     !
     if( ltopo(pelty) == -1 .or. ltopo(pelty) == 0 ) then
        volum = 0.0_rp
        do inode = 1,pnode
           call elmder(&
                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                elcod,gpcar,gpdet,xjacm,xjaci)
           volum = volum + elmar(pelty)%weigc(inode)*gpdet
        end do
        if(      ndime == 1_ip ) then
           hleng_aux = hleng(ndime) 
        else if( ndime == 2_ip ) then
           hleng_aux = sqrt(volum)
        else if( ndime == 3_ip ) then
           hleng_aux = volum**(1.0_rp/3.0_rp)
        end if
        chale(1) = hleng_aux
        chale(2) = chale(1)
     else
        chale(1)=hleng(ndime)       
        chale(2)=chale(1)
     end if

  else if(kfl_ellen==8) then
     !
     ! Volume/largest face (Hexahedrons)
     !
     if( ltopo(pelty) == 0_ip ) then
        volum = 0.0_rp
        do inode = 1,pnode
           call elmder(&
                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                elcod,gpcar,gpdet,xjacm,xjaci)
           volum = volum + elmar(pelty)%weigc(inode)*gpdet
        end do
        if( ndime == 3_ip ) then
           call elmafa(pelty,pnode,elcod,aface)
        end if
        chale(1) = volum/aface(1)
        chale(2) = chale(1)
     else
        chale(1)=hleng(ndime)       
        chale(2)=chale(1)
     end if
    
  end if
  !
  ! Divide h by 2 for quadratic elements and 3 for cubic elements
  !
  chale(1) = chale(1)/real(porde,rp)
  chale(2) = chale(2)/real(porde,rp)

end subroutine elmchl
