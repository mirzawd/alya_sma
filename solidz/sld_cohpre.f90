!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_cohpre(itask,ielem,iboun,igaub,xdisp,tract,stiff)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_cohpre
  ! NAME
  !    sld_cohpre
  ! DESCRIPTION
  !    Compute some Gauss values
  ! INPUT
  !    XDISP ... Jump in displacement at interface Gauss point ... delta
  ! OUTPUT
  !    TRACT ... Cohesive traction at interface Gauss point ...... t(delta)
  !    STIFF ... Cohesive tangent stiffness for implicit ......... dt/ddelta
  ! USES
  !    sld_cohesive_law_xxx
  ! USED BY
  !    sld_cohesi
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,nmate,lelch
  use def_solidz, only       :  cranx_sld,lawch_sld,lmate_sld
  use def_solidz, only       :  lecoh_sld,dfric_sld,dceff_sld,dslip_sld
  use def_solidz, only       :  cohnx_sld
  use def_elmtyp
  use def_master 
  implicit none
  integer(ip), intent(in)    :: igaub,iboun,itask,ielem
  real(rp),    intent(in)    :: xdisp(ndime)
  real(rp),    intent(out)   :: tract(ndime),stiff(ndime,ndime)
  real(rp)                   :: cnorm(ndime),ctant(ndime)
  real(rp)                   :: delta(2),tcbar(2),yybar(2,2)
  real(rp)                   :: stemp(2,ndime)
  integer(ip)                :: idime,jdime,pmate
  
  logical                    :: debugging
  real(rp), parameter        :: pert = 1.0e-8_rp 
  
  debugging = .false. !(ielem == 916 .and. ittim == 30)
  
  pmate = 1
  if ( nmate > 1 ) then
    pmate = lmate_sld(ielem)
  end if
 
     !
     ! DELTA(1): opening in normal direction and CNORM: crack normal vector
     !
     delta(1) = 0.0_rp
     if (lelch(ielem) == ELCUT ) then
        cnorm = cranx_sld(:,ielem)
     else if (lelch(ielem)==ELCOH) then
        cnorm = cohnx_sld(:,ielem)
     end if 
     do idime = 1,ndime
        delta(1)  = delta(1)  + xdisp(idime)*cnorm(idime)
     end do
     !
     ! DELTA(2): opening in tangential direction and CTANT: crack tangential direction
     !
     delta(2) = 0.0_rp
     ctant    = 0.0_rp
     do idime = 1,ndime
        ctant(idime) = xdisp(idime) - delta(1)*cnorm(idime)
        delta(2) = delta(2) + ctant(idime)*ctant(idime)
     end do
     delta(2) = sqrt(delta(2))
     if (delta(2) > 1.0e-12_rp) ctant = ctant/delta(2)

     !
     ! Get traction TCBAR and tangent stiffness YYBAR in local coordinate system
     !
     if (delta(1) > 0.0_rp) then
        !
        ! Case 1: Crack opening (cohesive or tractino-free)
        !
        if (debugging) then
           write(*,*)' '
           write(*,*)'<<<crack opening>>>'
        end if
        if (lecoh_sld(ielem) /= 0) then
           !
           ! Compute traction from selected cohesive law
           !
           if (lawch_sld(pmate)==900) then
              call sld_cohesive_law_900(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)
           else if (lawch_sld(pmate)==901) then
              call sld_cohesive_law_901(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)
           else if (lawch_sld(pmate)==902) then
              call sld_cohesive_law_902(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)
           else if (lawch_sld(pmate)==903) then
              call sld_cohesive_law_903(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)
           end if
           
        else
           !
           ! No traction upon crack opening (traction-free)
           !
           tcbar = 0.0_rp
           yybar = 0.0_rp
        end if
        dfric_sld(ielem,iboun,igaub,1) = 0.0_rp
        
     else 
        !
        ! Case 2: Crack closing (contact and friction)
        !
        if (debugging) then
           write(*,*)' '
           write(*,*)'<<<crack closing>>>'
        end if
        call sld_contac(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)
        dceff_sld(ielem,iboun,igaub,1) = 0.0_rp

     end if
     dslip_sld(ielem,iboun,igaub,1) = delta(2)
     
     !
     ! TRACT: surface traction in global coordinate system
     !
     tract = 0.0_rp
     do idime = 1,ndime
        tract(idime) = tract(idime) + tcbar(1)*cnorm(idime) + tcbar(2)*ctant(idime)
     end do
     
     !
     ! STIFF: tangent stiffness in global coordinate system (in case implicit)
     !
     if (itask == 2) then
        stemp = 0.0_rp
        do idime = 1,ndime
           stemp(1,idime) = stemp(1,idime) + yybar(1,1)*cnorm(idime) + yybar(1,2)*ctant(idime)
           stemp(2,idime) = stemp(2,idime) + yybar(2,1)*cnorm(idime) + yybar(2,2)*ctant(idime)
        end do
        stiff = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              stiff(idime,jdime) = stiff(idime,jdime) + cnorm(idime)*stemp(1,jdime) &
                                                     + ctant(idime)*stemp(2,jdime)
           end do
        end do
     end if
     
  
  if (debugging) then
    write(*,*)'ielem = ',ielem,' ** iboun = ',iboun,' ** igaub = ',igaub
    write(*,*)'tract =     ** stiff =                  ** xdisp = '
    do idime = 1,ndime
      write(*,'(1x,e12.5,2x,2(1x,e12.5),2x,1x,e12.5)')tract(idime),(stiff(idime,jdime),jdime =1,ndime),xdisp(idime)
    end do
  end if
  
100 format (2(F16.8,','))
end subroutine sld_cohpre

