!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updrel()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_updrel
  ! NAME 
  !    nsi_updrel
  ! DESCRIPTION
  !    This routine updates the relaxation factor
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  implicit none
  integer(ip) :: ipoin,idime,idofn,jdofn

  if( kfl_relax_nsi == 2 .or. kfl_relap_nsi == 2 ) then
     !
     ! Aitken's method
     !
     if(itinn(modul)/=1) then 

        if( NSI_MONOLITHIC ) then
           !
           ! Monolithic
           !
           if(kfl_relax_nsi==2)&
                call aitken(&
                0_ip,ndime,solve(1)%ndofn,ndime,veloc,unkno,dunkn_nsi,&
                one,one,one,ndime,relpa_nsi(1))
           if(kfl_relap_nsi==2)&
                call aitken(&
                 0_ip,1_ip,solve(1)%ndofn,1_ip,press,unkno,dunkp_nsi,&
                 one,solve(1)%ndofn,one,one,relpa_nsi(2)) 
             
        else 
           !
           ! BGS
           !
           if(ivari_nsi==ivari_nsi_mom.and.kfl_relax_nsi==2)&
                call aitken(&
                0_ip,ndime,solve(1)%ndofn,ndime,veloc,unkno,dunkn_nsi,&
                one,one,one,ndime,relpa_nsi(1))
           if(ivari_nsi==ivari_nsi_cont.and.kfl_relap_nsi==2)&
                call aitken(&
                0_ip,1_ip,1_ip,1_ip,press,unkno(ndbgs_nsi+1),dunkp_nsi,&
                one,one,one,one,relpa_nsi(2))

        end if

     else

        if(kfl_paral/=0) then

           if( NSI_MONOLITHIC ) then
              !
              ! Monolithic
              !
              if(kfl_relax_nsi==2) then
                 do ipoin=1,npoin
                    idofn=(ipoin-1)*ndime
                    jdofn=(ipoin-1)*solve(1)%ndofn
                    do idime=1,ndime
                       idofn=idofn+1
                       jdofn=jdofn+1
                       dunkn_nsi(idime,ipoin)=veloc(idime,ipoin,1)-unkno(jdofn)
                    end do
                 end do
              end if
              if(kfl_relap_nsi==2) then
                 do ipoin=1,npoin
                    jdofn=ipoin*solve(1)%ndofn
                    dunkp_nsi(ipoin)=press(ipoin,1)-unkno(jdofn)
                 end do
              end if
              
           else
              !
              ! BGS
              !
              if(ivari_nsi==ivari_nsi_mom.and.kfl_relax_nsi==2) then
                 do ipoin=1,npoin
                    idofn=(ipoin-1)*ndime
                    do idime=1,ndime
                       idofn=idofn+1
                       dunkn_nsi(idime,ipoin)=veloc(idime,ipoin,1)-unkno(idofn)
                    end do
                 end do
              end if
              if(ivari_nsi==ivari_nsi_cont.and.kfl_relap_nsi==2) then
                 do ipoin=1,npoin
                    dunkp_nsi(ipoin)=press(ipoin,1)-unkno(ndbgs_nsi+ipoin)
                 end do
              end if

           end if
        end if
        relpa_nsi(1) = 0.0_rp
        relpa_nsi(2) = 0.0_rp
     end if
     relax_nsi = 1.0_rp-relpa_nsi(1)
     relap_nsi = 1.0_rp-relpa_nsi(2)
  end if

end subroutine nsi_updrel
