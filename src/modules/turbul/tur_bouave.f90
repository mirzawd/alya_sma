!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_bouave(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_bouave
  ! NAME 
  !    tur_bouave
  ! DESCRIPTION
  !    This routine computes the mean velocity on boundaries which have
  !    an inlet condition for k 
  ! USES
  ! USED BY
  !    tur_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_bouder
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: baloc(ndime,ndime)
  real(rp)                :: bocod(ndime,mnodb)
  real(rp)                :: bovel(ndime,mnodb) 
  integer(ip)             :: ipoin,fixn6,idime,igaub,iboun,inodb
  integer(ip)             :: pblty,pnodb,fixn8
  real(rp)                :: eucta,gpsur,surfa,gbvel(3),gbvno,gbkin,vevel
  real(rp),    target     :: dummr(2)

  if(itask==1) avvel_tur=0.0_rp
  if(itask==2) avkin_tur=0.0_rp
  surfa=0.0_rp
  !
  ! Loop over elements
  !
  if(itask==1.and.kfl_avvel_tur==1) then

     if( INOTMASTER ) then

        if(nboun==0) then
           do ipoin=1,npoin
              if(lpoty(ipoin)/=0) then
                 if(kfl_fixno_tur(1,ipoin,1)==6.or.kfl_fixno_tur(1,ipoin,1)==8) then
                    vevel=0.0_rp
                    do idime=1,ndime
                       vevel=vevel+advec(idime,ipoin,1)*advec(idime,ipoin,1)
                    end do
                    avvel_tur = avvel_tur+sqrt(vevel)
                    surfa     = surfa+1.0_rp
                 end if
              end if
           end do
        else
           boundaries: do iboun=1,nboun

              pblty=ltypb(iboun)
              pnodb=nnode(pblty)
              fixn6=0
              fixn8=0
              boundary_nodes:  do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 if(kfl_fixno_tur(1,ipoin,1)==6) fixn6=fixn6+1
                 if(kfl_fixno_tur(1,ipoin,1)==8) fixn8=fixn8+1
              end do boundary_nodes

              if(fixn6==pnodb.or.fixn8==pnodb) then
                 !
                 ! Gather operations
                 !
                 bocod(1:ndime,1:pnodb)=coord(1:ndime,lnodb(1:pnodb,iboun))
                 bovel(1:ndime,1:pnodb)=advec(1:ndime,lnodb(1:pnodb,iboun),1)

                 gauss_points: do igaub=1,ngaus(pblty)

                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&
                         bocod,baloc,eucta)
                    gpsur=elmar(pblty)%weigp(igaub)*eucta 
                    surfa=surfa+gpsur

                    if( itask==1 ) then
                       !
                       ! Averaged velocity module: AVVEL
                       !
                       gbvel=0.0_rp
                       do inodb=1,pnodb
                          ipoin=lnodb(inodb,iboun)
                          do idime=1,ndime
                             gbvel(idime)=gbvel(idime)&
                                  +elmar(pblty)%shape(inodb,igaub)*advec(idime,ipoin,1)
                          end do
                       end do
                       call vecnor(gbvel,ndime,gbvno,2_ip)
                       avvel_tur=avvel_tur+gbvno*gpsur

                    else if( itask==2 ) then
                       !
                       ! Averaged turbulence kinetic energy: AVKIN
                       !
                       gbkin=0.0_rp
                       do inodb=1,pnodb
                          ipoin=lnodb(inodb,iboun)
                          gbkin=gbkin&
                               +elmar(pblty)%shape(inodb,igaub)*untur(1,ipoin,1)
                       end do
                       avkin_tur=avkin_tur+gbkin*gpsur
                    end if

                 end do gauss_points

              end if

           end do boundaries
        end if

     end if

     if( itask==1 .and. IPARALL ) then
        !
        ! Parallel: reduce sum
        !
        nparr     =  2
        dummr(1)  =  avvel_tur
        dummr(2)  =  surfa
        parre     => dummr
        call par_operat(3_ip)    
        avvel_tur =  dummr(1)  
        surfa     =  dummr(2)
     end if
     !
     ! Averaged values
     !
     if( surfa>0.0_rp ) then     
        if( itask==1 ) then
           avvel_tur=avvel_tur/surfa
        else if( itask==2 ) then
           avkin_tur=avkin_tur/surfa
        end if
     end if

  end if

end subroutine tur_bouave

