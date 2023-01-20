!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_intbcs()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_intbcs
  ! NAME 
  !    tem_intbcs
  ! DESCRIPTION
  !    This routine interpolates boundary conditions
  ! USED BY
  !    tem_updbcs
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use def_temper
  use mod_memchk
  use mod_ecoute, only      : ecoute
  use mod_memory, only      : memory_alloca, memory_deallo
  use mod_ecoute, only      : ecoute_set_read_unit
  use mod_ecoute, only      : ecoute_set_write_unit
  use mod_bouder
  implicit none
  integer(ip)          :: ipoin,idime,knodb(mnodb),ibopo,pnodb,pblty
  integer(ip)          :: iboun,inodb
  integer(ip), pointer :: itbck(:)
  real(rp),    pointer :: rtbck(:),vbmas(:)
  real(rp)             :: bocod(ndime,mnodb),baloc(ndime,ndime)
  real(rp)             :: eucta,teold,tenew,gbdet

  if(kfl_paral==-1) then

     select case(kfl_intbc_tem)

     case(2)
        !
        ! Per boundary (e.g. Marek Prymon)
        !
        call memgen(zero,nbopo,zero)
        call memory_alloca(mem_modul(1:2,modul),'VBMAS','tem_intbcs',vbmas,nbopo)
        call memory_alloca(mem_modul(1:2,modul),'RTBCK','tem_intbcs',rtbck,nboun)
        call memory_alloca(mem_modul(1:2,modul),'ITBCK','tem_intbcs',itbck,nboun)
        !
        ! Reads boundary conditions from file
        !
        call ecoute_set_read_unit (lun_intbc_tem)            ! Reading file
        call ecoute_set_write_unit(momod(modul) % lun_outpu) ! Writing file
        words(1)=''
        param(1)=-999999.0_rp

        do while(words(1)/='TIME ')
           call ecoute('tem_intbcs')
        end do
        !do while(abs(param(1)-cutim)>1.0e-3)
        !   do while(words(1)/='TIME ')
        !      call ecoute('tem_intbcs')
        !   end do
        !end do
        call ecoute('tem_intbcs')
        if(iknbo_tem==1) then                   ! Unknown boundary
           do while(words(1)/='ENDTI')
              pnodb=nnpar-1
              knodb(1:pnodb)=int(param(1:pnodb))
              call finbou(pnodb,knodb,iboun)
              if(iboun<1.or.iboun>nboun) &
                   call runend('TEM_INTBCS: ERROR WHILE READING BOUNDARY')
              itbck(iboun)=1
              rtbck(iboun)=param(nnpar)
              call ecoute('tem_intbcs')
           end do
        else                                    ! Known boundary
           do while(words(1)/='ENDTI')
              iboun=int(param(1))
              if(iboun<1.or.iboun>nboun) &
                   call runend('TEM_INTBCS: ERROR WHILE READING BOUNDARY')
              itbck(iboun)=1
              rtbck(iboun)=param(2)
              call ecoute('tem_intbcs')
           end do
        end if
        !
        ! Initialization
        !
        do ibopo=1,nbopo
           gesca(ibopo)=0.0_rp
        end do
        !
        ! Project values on boundary nodes
        !
        boundaries: do iboun=1,nboun

           if(itbck(iboun)/=0) then
              pblty=ltypb(iboun)
              pnodb=nnode(pblty)
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 do idime=1,ndime
                    bocod(idime,inodb)=coord(idime,ipoin)
                 end do
              end do
              !
              ! Loop over Gauss points (which are nodes)
              !
              gauss_points: do inodb=1,pnodb
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deric(:,:,inodb),&
                      bocod,baloc,eucta)
                 gbdet        = elmar(pblty)%weigc(inodb)*eucta
                 !if(kfl_naxis==1) then
                 !   gbdet=gbdet*bocod(1,inodb)*twopi
                 !end if
                 ipoin        = lnodb(inodb,iboun)
                 ibopo        = lpoty(ipoin)
                 gesca(ibopo) = gesca(ibopo)+gbdet*rtbck(iboun)
                 vbmas(ibopo) = vbmas(ibopo)+gbdet
              end do gauss_points

           end if
        end do boundaries
        !
        ! Update fixity boundary condition (and solve diagonal system)
        !
        do ipoin=1,npoin
           if(kfl_fixno_tem(1,ipoin)==5) then
              ibopo=lpoty(ipoin)
              if(ibopo/=0) then
                 if(vbmas(ibopo)==0.0_rp) then
                    call runend('TEM_INTBCS: ZERO BOUNDARY MASS MATRIX')
                 else
                 tenew = gesca(ibopo)/vbmas(ibopo)
                 if(kfl_timei_tem/=0.and.kfl_tiacc_tem==2.and.kfl_tisch_tem==1) then
                    teold=therm(ipoin,3)
                    bvess_tem(1,ipoin,1)=0.50_rp*tenew+0.50_rp*teold
                 else
                    bvess_tem(1,ipoin,1)=tenew                 
                 end if
              end if
              end if
           end if
        end do
        !
        ! Deallocate memory
        !
        call memgen(two,one,zero)
        call memory_deallo(mem_modul(1:2,modul),'VBMAS','tem_intbcs',vbmas)
        call memory_deallo(mem_modul(1:2,modul),'RTBCK','tem_intbcs',rtbck)
        call memory_deallo(mem_modul(1:2,modul),'ITBCK','tem_intbcs',itbck)

     end select

  end if

end subroutine tem_intbcs
