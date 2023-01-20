!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine reacod(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/reacod
  ! NAME
  !    reacod
  ! DESCRIPTION
  !    Reads codes on nodes and boundaries
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use def_inpout
  use mod_iofile
  use mod_chktyp, only : check_type
  use mod_ecoute, only : ecoute
  use mod_outfor, only : outfor
  use mod_memory, only : memory_size
  use mod_windk,  only : mod_windk_create_system
  use mod_codes
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,iboun,icode,idofn,iparb,ibopo
  integer(ip)             :: ndofn,ncode,ibves,kpoin,ierro
  integer(ip)             :: jcode,icono,nmcod,kcode(mcono),ivcod,ivcob
  integer(ip)             :: ntotn,mcono_tmp
  character(20)           :: messa
  integer(ip), pointer    :: kfl_codes_tmp(:,:)

  ierro = 0

  if( itask == 1 .or. itask == -1 .or. itask == 100 ) then

     !-------------------------------------------------------------------
     !
     ! Master reads node codes
     !
     !-------------------------------------------------------------------

     call codes_read_on_nodes(itask)

  else if( itask == 2 .or. itask == 200 ) then

     !-------------------------------------------------------------------
     !
     ! Master reads boundary codes
     !
     !-------------------------------------------------------------------

     call codes_read_on_boundaries(itask)

  else if ( itask == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Master reads geometrical node codes
     !
     !-------------------------------------------------------------------

     call ecoute('reacod')
     ncode = 0
     ndofn = tgcod(1) % ndofn

     do while(words(1)/='ENDCO')

        if(words(2)/='OFF  ') then

           ncode = ncode + 1

           if( param(1) < 0 ) then
              tgcod(1) % l(ncode) % lcode(1)  = int(abs(param(1)) )
              tgcod(1) % l(ncode) % kfl_value = int( param(2) )
           else
              if( exists('FIELD') ) then
                 tgcod(1) % l(ncode) % lcode(1)  = int(abs(param(1)) )
                 tgcod(1) % l(ncode) % kfl_value = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
              else
                 tgcod(1) % l(ncode) % lcode(1)  = int( param(1) )
                 tgcod(1) % l(ncode) % kfl_value = 0
                 do idofn = 1,ndofn
                    tgcod(1) % l(ncode) % bvess(idofn) = param(1+idofn)
                 end do
              end if
           end if

        end if

        call ecoute('reacod')
     end do

     tgcod(1) % ncode = ncode

  else if ( itask == IMPOSE_NODE_CODES .or. itask == IMPOSE_EDGE_CODES ) then

     !-------------------------------------------------------------------
     !
     ! Impose node or edge codes (done only for untagged nodes)
     !
     !-------------------------------------------------------------------

     if( itask == IMPOSE_NODE_CODES ) then
        ntotn         =  npoin
        kfl_codes_tmp => kfl_codno
        mcono_tmp     =  mcono
     else if( itask == IMPOSE_EDGE_CODES ) then
        ntotn         =  meshe(ndivi) % nedge
        kfl_codes_tmp => kfl_coded
        mcono_tmp     =  2
     end if

     ndofn = memory_size(kfl_fixno,1_ip)

     if( itask == IMPOSE_NODE_CODES ) then
        if( ifbop == 1 .and. memory_size(kfl_fixno,2_ip) /= nbopo ) &
             call runend('REACOD: WRONG DIMENSIONS FOR FIXITY ARRAY')
     end if

!!!!!!!!!!!!!!!!!!!!
     ! ESTA LINEA ESTABA ASI PERO ES RARISIMA Y NO SE PARA QUE
     ! ME JODE EN UN PROBLEMA CUANDO HAY UN SOLO ELEMENTO
!!!!     if( ndofn == npoin ) ndofn = 1

     ! Y ME SIGUE JODIENDO ESTA MIERDA.. CHEQUEA LUEGO SI IBVES ES = A NPOIN
     if( ifbes == 1 ) ibves = memory_size(bvess,1_ip)
     ! ASI QUE LO CAGO Y QUE SE CAGUE
     ibves = ntotn + 1
!!!!!!!!!!!!!!!!!!!!

     do ncode = 1,tncod(1) % ncode

        !if (tncod(1) % l(ncode) % cotag  == '') then     !!! esta cosa hay que revisarla, no se puede descomentar!!!

        !
        ! Do it only when the code is untagged
        !
        if( itask == IMPOSE_NODE_CODES ) then
           if(      tncod(1) % l(ncode) % lcode(1) == mcodb+1 ) then
              nmcod = 0
           else if( tncod(1) % l(ncode) % lcode(2) == mcodb+1 ) then
              nmcod = 1
           else if( tncod(1) % l(ncode) % lcode(3) == mcodb+1 ) then
              nmcod = 2
           else
              nmcod = 3
           end if
        else
           if(      tncod(1) % l(ncode) % lcode(1) == mcodb+1 ) then
              nmcod = 0
           else if( tncod(1) % l(ncode) % lcode(2) == mcodb+1 ) then
              nmcod = 1
           else
              nmcod = 2
           end if
        end if
        !
        ! Check if value function exist
        !
        if( itask == IMPOSE_NODE_CODES ) then
           if( tncod(1) % l(ncode) % kfl_value > 0 ) then
              ivcod = tncod(1) % l(ncode) % kfl_value
              !call check_type(bvcod,ivcod,ndofn,npoin)
              call check_type(xfiel,ivcod,ndofn,npoin)
           end if
        end if
        !
        ! Order codes
        !
        do jcode = 1,mcono_tmp
           kcode(jcode) = tncod(1) % l(ncode) % lcode(jcode)
        end do
        call heapsorti1(2_ip,mcono_tmp,kcode)
        icode = tncod(1) % l(ncode) % lcode(1)

        do ipoin = 1,ntotn
           icono = 0
           do jcode = 1,mcono_tmp
              if( kfl_codes_tmp(jcode,ipoin) == abs(kcode(jcode)) ) icono = icono + 1
           end do

           if( icono == mcono_tmp ) then

              kpoin = ipoin
              if( itask == IMPOSE_NODE_CODES ) then
                 ibopo = lpoty(ipoin)
                 if( ifbop == 1 ) kpoin = ibopo
              end if

              if( kpoin == 0 ) then

                 messa = intost(ipoin)
                 ierro = ierro + 1
                 call outfor(2_ip,lun_outpu,&
                      'BOUNDARY CONDITION CANNOT BE IMPOSED ON INTERIOR NODE: '//trim(messa))
              else

                 kfl_fixno(1,kpoin) = tncod(1) % l(ncode) % kfl_fixno

                 call codfix(ndofn,kfl_fixno(1,kpoin))

                 if( ifbes == 1 ) then
                    if( ibves == ntotn ) then
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          bvess(kpoin,1) = tncod(1) % l(ncode) % bvess(1)
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          bvess(kpoin,1) = xfiel(ivcod) % a(1,ipoin,1)
                       end if
                    else
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = tncod(1) % l(ncode) % bvess(idofn)
                          end do
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = xfiel(ivcod) % a(idofn,ipoin,1)
                          end do
                       end if
                    end if
                 end if
                 if( iffun /= 0 ) then
                    kfl_funno(kpoin) = tncod(1) % l(ncode) % kfl_funno
                    kfl_funtn(kpoin) = tncod(1) % l(ncode) % kfl_funtyp
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(kpoin) = tncod(1) % l(ncode) % kfl_fixrs
                          !if( kfl_fixrs(ibopo) == -2 .and. nskew > 0 ) call geofix(kpoin,ibopo)
                       end if
                    end if
                 else
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(kpoin) = tncod(1) % l(ncode) % kfl_fixrs
                          !if(kfl_fixrs(ibopo)==-2.and.nskew>0) call geofix(kpoin,ibopo)
                       end if
                    end if
                 end if
              end if

           end if

        end do

        !end if

     end do

  else if ( itask == 20 ) then

     !-------------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------------

     do ncode = 1,tbcod(1) % ncode
        ivcob = tbcod(1) % l(ncode) % kfl_value
        if( ivcob == 0 ) then
           do iboun = 1,nboun
              if( kfl_codbo(iboun) == tbcod(1) % l(ncode) % lcode ) then
                 kfl_fixbo(iboun) = tbcod(1) % l(ncode) % kfl_fixbo
                 do iparb = 1,tbcod(1) % ndofn
                    bvnat(iparb,iboun) = tbcod(1) % l(ncode) % bvnat(iparb)
                 end do
                 if( iffun == 1 ) then
                    kfl_funbo(iboun) = tbcod(1) % l(ncode) % kfl_funbo
                    kfl_funtb(iboun) = tbcod(1) % l(ncode) % kfl_funtyp
                 end if
              end if
           end do
        else
           do iboun = 1,nboun
              if( kfl_codbo(iboun) == abs(tbcod(1) % l(ncode) % lcode) ) then
                 kfl_fixbo(iboun) = tbcod(1) % l(ncode) % kfl_fixbo
                 do iparb = 1,tbcod(1) % ndofn
                    bvnat(iparb,iboun) = xfiel(ivcob) % a(iparb,iboun,1)
                 end do
                 if( iffun == 1 ) then
                    kfl_funbo(iboun) = tbcod(1) % l(ncode) % kfl_funbo
                    kfl_funtb(iboun) = tbcod(1) % l(ncode) % kfl_funtyp
                 end if
              end if
           end do
        end if
     end do

  else if( itask == 30 ) then

     !-------------------------------------------------------------------
     !
     ! Initial conditions on nodes, only for the codes tagged as INITI
     !
     !-------------------------------------------------------------------


     ! OJO QUE NO ESTA LISTO AUN

     ndofn = memory_size(kfl_fixno,1_ip)
     if( memory_size(kfl_fixno,2_ip) /= nbopo ) &
          call runend('REACOD: WRONG DIMENSIONS FOR INITIAL CONDITIONS FIXITY ARRAY')
     if( ndofn == npoin ) ndofn = 1
     if( ifbes == 1 )     ibves = memory_size(bvess,1_ip)

     do ncode = 1,tncod(1) % ncode

        if(      tncod(1) % l(ncode) % lcode(1) == mcodb+1 ) then
           nmcod = 0
        else if( tncod(1) % l(ncode) % lcode(2) == mcodb+1 ) then
           nmcod = 1
        else if( tncod(1) % l(ncode) % lcode(3) == mcodb+1 ) then
           nmcod = 2
        else
           nmcod = 3
        end if
        !
        ! Order codes
        !
        do jcode = 1,mcono
           kcode(jcode) = tncod(1) % l(ncode) % lcode(jcode)
        end do
        call heapsorti1(2_ip,mcono,kcode)
        icode = tncod(1) % l(ncode) % lcode(1)

        do ipoin = 1,npoin
           icono = 0
           do jcode = 1,mcono
              if( kfl_codes_tmp(jcode,ipoin) == abs(kcode(jcode)) ) icono = icono + 1
           end do

           if( icono == mcono ) then

              ibopo = lpoty(ipoin)
              kpoin = ipoin
              if( ifbop == 1 ) kpoin = ibopo

              if( kpoin == 0 ) then

                 messa = intost(ipoin)
                 ierro = ierro + 1
                 call outfor(2_ip,lun_outpu,&
                      'BOUNDARY CONDITION CANNOT BE IMPOSED ON INTERIOR NODE: '//trim(messa))
              else

                 kfl_fixno(1,kpoin) = tncod(1) % l(ncode) % kfl_fixno
                 call codfix(ndofn,kfl_fixno(1,kpoin))

                 if( ifbes == 1 ) then
                    if( ibves == npoin ) then
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          bvess(kpoin,1) = tncod(1) % l(ncode) % bvess(1)
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          bvess(kpoin,1) = xfiel(ivcod) % a(1,ipoin,1)
                       end if
                    else
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = tncod(1) % l(ncode) % bvess(idofn)
                          end do
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = xfiel(ivcod) % a(idofn,ipoin,1)
                          end do
                       end if
                    end if
                 end if

                 if( iffun /= 0 ) then
                    kfl_funno(kpoin) = tncod(1) % l(ncode) % kfl_funno
                    kfl_funtn(kpoin) = tncod(1) % l(ncode) % kfl_funtyp
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(kpoin) = tncod(1) % l(ncode) % kfl_fixrs
                          !if( kfl_fixrs(ibopo) == -2 .and. nskew > 0 ) call geofix(kpoin,ibopo)
                       end if
                    end if
                 else
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(kpoin) = tncod(1) % l(ncode) % kfl_fixrs
                          !if(kfl_fixrs(ibopo)==-2.and.nskew>0) call geofix(kpoin,ibopo)
                       end if
                    end if
                 end if
              end if

           end if
        end do

     end do

  end if

  messa=intost(ierro)
  if(ierro==1) then
     call runend('REACOD: '//trim(messa)//' ERROR HAS BEEN FOUND')
  else if(ierro>=2) then
     call runend('REACOD: '//trim(messa)//' ERRORS HAVE BEEN FOUND')
  end if
  !
  ! Recover original values
  !
  iffun = 0
  ifloc = 0
  ifbop = 0
  ifbes = 1

end subroutine reacod

subroutine geofix(kpoin,ibopo)
  use def_kintyp, only    :  ip
  use def_domain, only    :  kfl_fixno,ndime,lpoin
  implicit none
  integer(ip), intent(in) :: kpoin,ibopo
  integer(ip)             :: idime,kfl_fixn2(2)

  kfl_fixn2(1)=kfl_fixno(    1,kpoin)
  kfl_fixn2(2)=kfl_fixno(ndime,kpoin)

  if(lpoin(ibopo)==0) then
     !
     ! 0 patch
     !
     kfl_fixno(1,kpoin)=kfl_fixn2(1)
     do idime=2,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(2)
     end do

  else if(lpoin(ibopo)==1) then
     !
     ! 1 patch
     !
     kfl_fixno(1,kpoin)=kfl_fixn2(1)
     do idime=2,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(2)
     end do

  else if(lpoin(ibopo)==2) then
     !
     ! 2 patches
     !
     do idime=1,min(2_ip,ndime)
        kfl_fixno(idime,kpoin)=kfl_fixn2(1)
     end do
     if(ndime==3) kfl_fixno(ndime,kpoin)=kfl_fixn2(2)

  else if(lpoin(ibopo)==3) then
     !
     ! 3 patches: corner of step type
     !
     do idime=1,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(1)
     end do

  else if(lpoin(ibopo)==-3) then
     !
     ! 3 patches: corner of bottom type
     !
     do idime=1,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(1)
     end do

  end if

end subroutine geofix
