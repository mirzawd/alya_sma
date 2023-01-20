!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    geonor.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Computes geometrical boundary conditions
!> @details The arrays computed here are:
!!          \verbatim
!!
!!          - KFL_GEOBO(1:NBOUN): 10,20,30,40,50,60
!!      
!!            Prescribed ......... 10
!!            Freestream ......... 20 
!!            Wall law ........... 30 
!!            Symmetry / Splip ... 40
!!            No slip / Wall ..... 50
!!            Free / Outflow ..... 60
!!
!!          - KFL_GEONO(1:NBOPO) takes the following values:
!!
!!            Prescribed ......... 10
!!                                 11 Freestream becomes inflow
!!            Freestream ......... 20 Freestream is an outflow (=11 if this is an inflow)
!!            Wall law ........... 30 Local basis n, t1, t2: SKCOS along a plane
!!                                 31 Local basis n1, t, n2: SKCOS along a line
!!            Symmetry / Splip ... 40 Local basis n, t1, t2: SKCOS along a line
!!                                 41 Local basis n1, t, n2: SKCOS along a plane
!!            No slip / Wall ..... 50
!!            Free / Outflow ..... 60
!!
!!          - SKCOS(1:NDIME,1:NDIME,1:NBOPO): local basis
!!
!!          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine geonor(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain 
  use mod_memory
  use mod_chktyp,         only : check_type
  use mod_communications, only : PAR_GHOST_BOUNDARY_EXCHANGE
  use mod_communications, only : PAR_SUM,PAR_MAX
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_maths,          only : maths_local_orthonormal_basis
  use mod_maths,          only : maths_heap_sort
  use mod_outfor,         only : outfor
  use mod_messages,       only : messages_live
  use mod_domain,         only : domain_memory_allocate
  use mod_domain,         only : domain_memory_deallocate
  implicit none
  integer(ip), intent(in)  :: itask
  integer(ip), pointer     :: ledbo(:,:)
  integer(ip), pointer     :: kpoin(:)
  integer(ip), pointer     :: kfl_geon2(:)
  real(rp),    pointer     :: bouno(:,:)
  integer(ip)              :: inodb,ipoin,jnodb,jpoin,ista0,ista1
  integer(ip)              :: iboun,icolu,pnodb,ielem,icode,igeon
  integer(ip)              :: ibott,ibopo,jboun,kboun,idime,jdime
  integer(ip)              :: icoun(100),ierro,ipoi1,ipoi2,ierr2
  integer(ip)              :: ibou1,ninve,kfl_outfl,kfl30,kfl40
  integer(ip)              :: nboun_tmp,npoin_tmp,pblty,jerro
!  integer(ip)              :: iposi
  real(rp)                 :: vnorm,cosim,extno,xangl(3)
  real(rp)                 :: tange(3,2),vecau(3),exwor(ndime,2)
  integer(ip), allocatable :: lninv_tmp(:)
  !
  ! Nullify pointers
  !
  nullify(ledbo)
  nullify(kpoin)
  nullify(kfl_geon2)  
  nullify(bouno)

  ibou1 = -1_ip
  ipoi1 = -1_ip
  ipoi2 = -1_ip

  if( kfl_geome == 1 ) then

     call messages_live('COMPUTE GEOMETRICAL NORMALS')

     ierro = 0
     ierr2 = 0 
     jerro = 0
     !
     ! Periodic ndoes
     !
     if( nperi > 0 ) call runend('GEONOR: CANNOT BE USED WITH PERIODICITY')

     if( INOTMASTER ) then 
        !
        ! Save geometry dimensions: we are going to consider fringe geometry
        !
        nboun_tmp = nboun
        npoin_tmp = npoin
        !
        ! Boundary graph required C_BOU, R_BOU: must include fringe geometry
        !
        nboun = nboun_2
        npoin = npoin_2
        call domain_memory_deallocate('BOUNDARY GRAPH')
        call bougra()
        nboun = nboun_tmp
        npoin = npoin_tmp

        if( itask == 1 ) then
           !
           ! Exchange KFL_GEOBO on fringe nodes
           !
           call PAR_GHOST_BOUNDARY_EXCHANGE(kfl_geobo,'SUBSTITUTE','IN MY CODE')
           !
           ! SKCOS, KFL_GEONO: Allocate memory
           !
           call domain_memory_allocate('SKCOS')
           call domain_memory_allocate('KFL_GEONO') ! KGL_GEONO
        end if
        !
        ! Initialize SKCOS
        !
        do ibopo = 1,nbopo
           do idime = 1,ndime
              do jdime = 1,ndime
                 skcos(jdime,idime,ibopo) = exnor(jdime,idime,ibopo) 
              end do
           end do
        end do
        !
        ! Allocate temporal memory
        !
        call memory_alloca(memor_dom,'LEDBO'    ,'geonor',ledbo,2_ip,nzbou)
        call memory_alloca(memor_dom,'BOUNO'    ,'geonor',bouno,ndime,max(1_ip,nboun_2))
        call memory_alloca(memor_dom,'KPOIN'    ,'geonor',kpoin,npoin_2)
        call memory_alloca(memor_dom,'KFL_GEON2','geonor',kfl_geon2,npoin_2)
        !
        ! Order C_BOU to ensure same result in sequential and in parallel
        !
        if( INOTMASTER ) then
           allocate(lninv_tmp(npoin_2))
           do ipoin = 1,npoin_2
              icolu = r_bou(ipoin+1)-r_bou(ipoin)
              if( icolu > 0 ) then
                 ista0              = r_bou(ipoin)
                 ista1              = r_bou(ipoin+1)-1
                 lninv_tmp(1:icolu) = lninv_loc(c_bou(ista0:ista1))
                 call maths_heap_sort(2_ip,icolu,lninv_tmp,ivo1=c_bou(ista0:))
              end if
           end do
           deallocate(lninv_tmp)        
        end if
     
        !----------------------------------------------------------------------
        !
        ! Compute boundary edges IEDGE
        ! LEDBO(1,iedge) = boundary 1
        ! LEDBO(2,iedge) = boundary 2 / 0
        !
        !----------------------------------------------------------------------

        !*OMP PARALLEL                                                      &
        !*OMP DEFAULT   (NONE)                                              &
        !*OMP PRIVATE   (iboun,pblty,pnodb,ielem,inodb,ipoin,jnodb,jpoin,   &
        !*OMP           ista0,ista1,icolu,ierro,ipoi1,ipoi2,ibou1,iposi)    &
        !*OMP SHARED    (nboun_2,ltypb,nnode,lelbo,lnodb,r_bou,c_bou,ledbo) 
        !*OMP DO

        do iboun = 1,nboun_2
           pblty = abs(ltypb(iboun))
           pnodb = nnode(pblty)
           ielem = lelbo(iboun)

           if( ielem >= 1 ) then  ! Volume element
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 do jnodb = 1,pnodb
                    if( inodb /= jnodb ) then
                       jpoin = lnodb(jnodb,iboun)                
                       ista0 = r_bou(ipoin  )-1
                       ista1 = r_bou(ipoin+1)-1
                       icolu = ista0
                       do while( icolu < ista1 )
                          icolu = icolu + 1
                          if( c_bou(icolu) == jpoin ) then
                             !
                             ! Check with OpenMP
                             !
                             if( ledbo(1,icolu) == 0 ) then
                                ledbo(1,icolu) = iboun
                             else if( ledbo(2,icolu) == 0 ) then
                                ledbo(2,icolu) = iboun
                             else
                                ierro = ierro + 1
                                ipoi1 = ipoin
                                ipoi2 = jpoin
                                ibou1 = iboun
                                print*,'error in geonor'
                                print*,ipoin,jpoin,ledbo(1,icolu),ledbo(2,icolu),iboun,lninv_loc(ipoin),lninv_loc(jpoin)
                                stop
                             end if
                             icolu = ista1
                          end if
                       end do
                    end if
                 end do
              end do
           end if
        end do

        !*OMP END DO
        !*OMP END PARALLEL

     end if

     !----------------------------------------------------------------------
     !
     ! Check errors
     !
     !----------------------------------------------------------------------

     call PAR_MAX(ierro)
     if( ierro /= 0 ) then
        ioutp(1) = ierro
        ioutp(2) = ipoi1
        ioutp(3) = ipoi2
        ioutp(4) = ibou1
        call outfor(45_ip,lun_outpu,' ')
     end if

     !----------------------------------------------------------------------
     !
     ! Compute the normals to the boundaries
     !
     !----------------------------------------------------------------------

     ninve = 0
     if( INOTMASTER .and. nboun_2 > 0 ) then
        call bounor(nboun_2,lnodb,ltypb,lelbo,ninve,bouno)
     end if
     call PAR_SUM(ninve)
     if( ninve > 0 ) then
        call outfor(2_ip,lun_outpu,'GEOMETRICAL BCS: SOME BOUNDARIES NORMALS HAVE BEEN INVERTED BECAUSE THEY WERE POINTING INWARDS')
     end if

     if( INOTMASTER ) then

        !----------------------------------------------------------------------
        !
        ! Initialization
        !
        !----------------------------------------------------------------------

        !
        ! Maximum angle and initialization
        !
        cosim = cos(geoan * pi / 180.0_rp)
        do ipoin = 1,npoin_2
           kpoin(ipoin)     = -1
           kfl_geon2(ipoin) =  0
        end do
        !
        ! High priority conditions
        !
        do iboun = 1,nboun_2
           if(  kfl_geobo(iboun) == 50 ) then  ! No slip wall
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 if( kfl_geon2(ipoin) == 0 ) then
                    kpoin(ipoin)     = -2
                    kfl_geon2(ipoin) = 50
                 end if
              end do
           end if
        end do

        do iboun = 1,nboun_2
           if(  kfl_geobo(iboun) == 10 ) then  ! Inflow
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 if( kfl_geon2(ipoin) == 0 ) then
                    kpoin(ipoin)     = -2
                    kfl_geon2(ipoin) = 10
                 end if
              end do
           end if
        end do

        do iboun = 1,nboun_2
           if(  kfl_geobo(iboun) == 30 .or.&   ! Wall law
                kfl_geobo(iboun) == 40 ) then  ! Slip wall or Symmetry
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 if( kpoin(ipoin) == -1 ) kpoin(ipoin) = 0
              end do
           end if
        end do

        if( ndime == 2 ) then

           !-------------------------------------------------------------------    
           !
           ! 2D
           !
           !-------------------------------------------------------------------

           do ipoin = 1,npoin
              ibopo = lpoty(ipoin)
              if( ibopo /= 0 ) then
                 if( kpoin(ipoin) == 0 ) then
                    !
                    ! Check if point is on a slip-wall or wall-law boundary
                    !
                    kfl_outfl = 0
                    kfl30     = 0
                    kfl40     = 0
                    exwor     = 0.0_rp

                    do iboun = 1,nboun_2

                       if(&
                            kfl_geobo(iboun) == 30 .or.&   ! Wall law
                            kfl_geobo(iboun) == 40 .or.&   ! Slip law or Symmetry
                            kfl_geobo(iboun) == 20 .or.&   ! Freestream
                            kfl_geobo(iboun) == 60 ) then  ! Outflow
                          pnodb = nnode(ltypb(iboun))
                          nodes_2d: do inodb = 1,pnodb
                             if( lnodb(inodb,iboun) == ipoin ) then
                                if( kfl_geobo(iboun) == 30 ) then
                                   kfl30 = kfl30 + 1
                                   kpoin(ipoin)                = kpoin(ipoin) + 1
                                   exwor(1:ndime,kpoin(ipoin)) = bouno(1:ndime,iboun)
                                else if( kfl_geobo(iboun) == 40 ) then
                                   kfl40 = kfl40 + 1
                                   kpoin(ipoin)                = kpoin(ipoin) + 1
                                   exwor(1:ndime,kpoin(ipoin)) = bouno(1:ndime,iboun)
                                else if( kfl_geobo(iboun) == 20 ) then
                                   kfl_outfl = 1
                                else if( kfl_geobo(iboun) == 60 ) then
                                   kfl_outfl = 1
                                end if
                                exit nodes_2d
                             end if
                          end do nodes_2d
                       end if
                    end do

                    if( kpoin(ipoin) /= 0 ) then
                       !
                       ! Prescribe in local system
                       !
                       if( kfl_outfl == 1 ) then
                          skcos(1:ndime,1,ibopo) =  exwor(1:ndime,1)
                          skcos(1,2,ibopo)       = -skcos(2,1,ibopo)
                          skcos(2,2,ibopo)       =  skcos(1,1,ibopo)  
                       end if
                       if( kfl30 /= 0 ) then
                          kfl_geon2(ipoin) = 30
                       else
                          kfl_geon2(ipoin) = 40
                       end if
                       !
                       ! Check if point is on an obtuse angle or if it is a corner
                       !
                       if( kpoin(ipoin) == 2 ) then
                          vnorm = dot_product(exwor(1:ndime,1),exwor(1:ndime,2))
                          exwor(1:ndime,1)=exwor(1:ndime,1)+exwor(1:ndime,2)
                          call vecuni(ndime,exwor,extno)
                          if( vnorm < cosim ) then                          ! Point is a corner
                             call corner(exwor,ipoin,ibott)

                             if(ibott==0) then                              ! Obtuse and step type: take exnor
                                if( kfl_convx == 0 ) then
                                   kfl_geon2(ipoin) = 60
                                else if( kfl_convx == 1 ) then
                                   if( kfl30 /= 0 ) then
                                      kfl_geon2(ipoin) = 30
                                   else
                                      kfl_geon2(ipoin) = 40
                                   end if
                                else if( kfl_convx == 2 ) then
                                   kfl_geon2(ipoin) = 50
                                else
                                   if( kfl30 /= 0 ) then
                                      kfl_geon2(ipoin) = 30
                                   else
                                      kfl_geon2(ipoin) = 40
                                   end if
                                   skcos(1:ndime,1,ibopo) =  exwor(1:ndime,1)
                                   skcos(1,2,ibopo)       = -skcos(2,1,ibopo)
                                   skcos(2,2,ibopo)       =  skcos(1,1,ibopo)
                                end if
                             else                                           ! Bottom type: u=0
                                kfl_geon2(ipoin) = 50
                             end if

                          end if
                       end if
                    end if
                 end if
              end if
           end do

        else

           !-------------------------------------------------------------------
           !
           ! 3D
           !
           !-------------------------------------------------------------------
 
           !$OMP PARALLEL                                                        &
           !$OMP DEFAULT   (NONE)                                                &
           !$OMP PRIVATE   (ipoin,ibopo,kfl_outfl,kfl30,kfl40,icolu,jpoin,jboun, &
           !$OMP           kboun,vnorm,tange,vecau,jerro,ibott,iboun)            &
           !$OMP SHARED    (npoin,kpoin,r_bou,c_bou,ledbo,kfl_geobo,bouno,       &
#ifndef NDIMEPAR
           !$OMP            ndime,                                               &
#endif
           !$OMP            cosim,kfl_geon2,skcos,exnor,kfl_convx,lpoty) 
           !$OMP DO

           do ipoin = 1,npoin

              ibopo = lpoty(ipoin)

              if( ibopo /= 0 ) then

                 if( kpoin(ipoin) == 0 ) then

                    kfl_outfl = 0
                    kfl30     = 0
                    kfl40     = 0

                    loop_segment_ipoin_jpoin: do icolu = r_bou(ipoin),r_bou(ipoin+1)-1
                       jpoin = c_bou(icolu)

                       if( jpoin /= ipoin ) then
                          ! Segment ipoin-jpoin

                          jboun = ledbo(1,icolu)
                          kboun = ledbo(2,icolu)
                          if( jboun > kboun ) then
                             iboun = kboun
                             kboun = jboun
                             jboun = iboun
                          end if
                          
                          if( jboun /= 0 .and. kboun /= 0 ) then

                             !-------------------------------------------------
                             !
                             ! Angle between jboun and kboun of ipoin-jpoin
                             !
                             !-------------------------------------------------

                             if(    ( kfl_geobo(jboun) == 30 .or. kfl_geobo(jboun) == 40 )  .and. &
                                  & ( kfl_geobo(kboun) == 30 .or. kfl_geobo(kboun) == 40 ) ) then  

                                if( kfl_geobo(jboun) == 30 .or. kfl_geobo(kboun) == 30 ) then
                                   kfl30 = kfl30 + 1
                                end if
                                if( kfl_geobo(jboun) == 40 .or. kfl_geobo(kboun) == 40 ) then
                                   kfl40 = kfl40 + 1
                                end if
                                ! n_jboun . n_kboun
                                
                                vnorm = dot_product(bouno(1:ndime,jboun),bouno(1:ndime,kboun))
                                
                                if( vnorm < cosim ) then
                                   ! n_jboun x n_kboun    
                                   call vecpro(bouno(1,jboun),bouno(1,kboun),vecau,ndime)
                                   ! kpoin(ipoin)= # of time the tangent has been calculated 
                                   ! at ipoin. Max 2 times.

                                   if( kpoin(ipoin) == 0 .or. kpoin(ipoin) == -1 ) then
                                      tange(1:ndime,1) = vecau(1:ndime)
                                      kpoin(ipoin)     = 1

                                   else if( kpoin(ipoin) == 1 ) then
                                      ! This point is a corner in 2D.
                                      tange(1:ndime,2) = vecau(1:ndime)
                                      kpoin(ipoin)     = 2

                                   else if( kpoin(ipoin) == 2 ) then
                                      ! This point is a corner in 3D.
                                      kpoin(ipoin)     = 3

                                   end if

                                end if

                             else if(    kpoin(ipoin) ==  0 .and. (&
                                  &  kfl_geobo(jboun) == 30 .or.   &
                                  &  kfl_geobo(jboun) == 40 )) then 

                                kpoin(ipoin)   = -1
                                vecau(1:ndime) = bouno(1:ndime,jboun)

                                if(  kfl_geobo(kboun) == 20 .or. &
                                     kfl_geobo(kboun) == 60 ) kfl_outfl = 1

                             else if(    kpoin(ipoin) ==  0 .and. (&
                                  &  kfl_geobo(kboun) == 30 .or.   &
                                  &  kfl_geobo(kboun) == 40 )) then  

                                kpoin(ipoin)   = -1
                                vecau(1:ndime) = bouno(1:ndime,kboun)
                                
                                if(  kfl_geobo(jboun) == 20 .or. &
                                     kfl_geobo(jboun) == 60 ) kfl_outfl = 1

                             end if

                          end if

                       end if

                    end do loop_segment_ipoin_jpoin

                    if( kpoin(ipoin) == -1 ) then 
                       !
                       ! Only one boundary
                       !                       
                       if( kfl30 /= 0 ) then
                          kfl_geon2(ipoin) = 30
                       else
                          kfl_geon2(ipoin) = 40
                       end if
                       if( kfl_outfl == 1 ) then
                          skcos(1:ndime,1,ibopo) = vecau(1:ndime)
                          call maths_local_orthonormal_basis(ndime,skcos(:,:,ibopo),jerro)
                       end if

                    else if( kpoin(ipoin) == 0 ) then 

                       if( kfl30 /= 0 ) then
                          kfl_geon2(ipoin) = 30
                       else
                          kfl_geon2(ipoin) = 40
                       end if
                       if( kfl_outfl == 1 ) then
                          skcos(1:ndime,1,ibopo) = vecau(1:ndime)
                          call maths_local_orthonormal_basis(ndime,skcos(:,:,ibopo),jerro)
                       end if

                    else if( kpoin(ipoin) == 1 .or. kpoin(ipoin) == 2 ) then

                       if( kpoin(ipoin) == 1 ) then

                          skcos(1,2,ibopo) = tange(1,1)
                          skcos(2,2,ibopo) = tange(2,1)
                          skcos(3,2,ibopo) = tange(3,1)

                       else

                          call vecuni(ndime,tange(1,1),vnorm)
                          call vecuni(ndime,tange(1,2),vnorm)
                          vnorm=dot_product(tange(1:ndime,1),tange(1:ndime,2))
                          ! If t_1,ipoin . t_2,ipoin < 0 invert oneS
                          if( vnorm < 1.0e-12_rp ) then
                             tange(1,2) = -tange(1,2)
                             tange(2,2) = -tange(2,2)
                             tange(3,2) = -tange(3,2)
                          end if
                          ! Volcar el promedio en skcos (base al punto de contorno)
                          skcos(1,2,ibopo) = tange(1,1) + tange(1,2)
                          skcos(2,2,ibopo) = tange(2,1) + tange(2,2)
                          skcos(3,2,ibopo) = tange(3,1) + tange(3,2)

                       end if

                       call vecuni(ndime,skcos(1,2,ibopo),vnorm)
                       ! Calcular la tangente 2 usando la normal inicial
                       call vecpro(skcos(1,1,ibopo),skcos(1,2,ibopo),vecau,ndime)
                       ! Calcular la normal de nuevo
                       call vecpro(skcos(1,2,ibopo),vecau,skcos(1,1,ibopo),ndime)
                       call vecuni(ndime,skcos(1,1,ibopo),vnorm)

                       if( vnorm < 1.0e-12_rp ) then
                          if( kfl30 /= 0 ) then
                             kfl_geon2(ipoin) = 30
                          else
                             kfl_geon2(ipoin) = 40
                          end if

                       else

                          call vecpro(skcos(1,1,ibopo),skcos(1,2,ibopo),skcos(1,3,ibopo),ndime)
                          call vecuni(ndime,skcos(1,3,ibopo),vnorm)
                          if( kfl30 /= 0 ) then
                             kfl_geon2(ipoin) = 31
                          else
                             kfl_geon2(ipoin) = 41
                          end if
                       end if

                    else if( kpoin(ipoin) == 3 ) then
                       ! Comprobar si es una esquina tipo de marcha o tipo de fondo      
                       call corner(exnor(1,1,ibopo),ipoin,ibott)

                       if(ibott==0) then                                ! Step type

                          if( kfl_convx == 0 ) then
                             kfl_geon2(ipoin) = 60
                          else if( kfl_convx == 1 ) then
                             if( kfl30 /= 0 ) then
                                kfl_geon2(ipoin) = 30
                             else
                                kfl_geon2(ipoin) = 40
                             end if
                          else if( kfl_convx == 2 ) then
                             kfl_geon2(ipoin) = 50
                          else
                             if( kfl30 /= 0 ) then
                                kfl_geon2(ipoin) = 30
                             else
                                kfl_geon2(ipoin) = 40
                             end if
                          end if
                       else                                             ! Bottom type: u=0
                          kfl_geon2(ipoin) = 50
                       end if
                    end if

                 end if
              end if            
           end do

           !$OMP END DO
           !$OMP END PARALLEL

        end if

        !----------------------------------------------------------------------
        !
        ! Put condition with low priority
        !
        !----------------------------------------------------------------------

        do iboun = 1,nboun_2
           if(  kfl_geobo(iboun) == 60 ) then  ! Outflow
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 if( ipoin <= npoin ) then 
                    ibopo = lpoty(ipoin)          
                    if( ibopo > 0 ) then
                       if( kfl_geon2(ipoin) == 0 ) then
                          kfl_geon2(ipoin) = 60
                       end if
                    end if
                 end if
              end do
           end if
        end do

        do iboun = 1,nboun_2
           if( kfl_geobo(iboun) == 20 ) then  ! Free stream
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 if( ipoin <= npoin ) then 
                    ibopo = lpoty(ipoin)          
                    if( ibopo > 0 ) then
                       if( kfl_geon2(ipoin) == 0 ) then
                          kfl_geon2(ipoin) = 20
                       end if
                    end if
                 end if
              end do
           end if
        end do

        !----------------------------------------------------------------------
        !
        ! When freestream should be converted to inflow
        !
        !----------------------------------------------------------------------
        
        if( kfl_frees == 0 ) then
           !
           ! Wind angle
           !
           do iboun = 1,nboun_2
              if( kfl_geobo(iboun) == 20 ) then 
                 xangl(1) = sin(awind)
                 xangl(2) = cos(awind)
                 xangl(3) = 0.0_rp               
                 vnorm    = dot_product(bouno(1:ndime,iboun),xangl(1:ndime))
                 if( vnorm <= epsilon(1.0_rp) + cos(0.5_rp*pi-tolan) ) then
                    do inodb = 1,nnode(ltypb(iboun))
                       ipoin = lnodb(inodb,iboun)
                       if( ipoin <= npoin ) then 
                          ibopo = lpoty(ipoin) 
                          if( ibopo > 0 ) then
                             kfl_geon2(ipoin) = 11
                          end if
                       end if
                    end do
                 end if
              end if
           end do
        else
           !
           ! Field
           !           
           call check_type(xfiel,kfl_frees,ndime,npoin) ! Check if value function exist
           call memgen(1_ip,npoin,0_ip)
           do iboun = 1,nboun
              if( kfl_geobo(iboun) == 20 ) then 
                 do inodb = 1,nnode(ltypb(iboun))
                    ipoin = lnodb(inodb,iboun)
                    
                    xangl(    1:ndime) = xfiel(kfl_frees) % a(    1:ndime,ipoin,1)
                    call vecuni(ndime,xangl,vnorm)
                    vnorm        = dot_product(bouno(1:ndime,iboun),xangl(1:ndime))

                    if( vnorm <= epsilon(1.0_rp) + cos(0.5_rp*pi-tolan) ) then
                       if( ipoin <= npoin ) then 
                          ibopo = lpoty(ipoin) 
                          if( ibopo > 0 ) then
                             gisca(ipoin) = 1
                          end if
                       end if
                    end if
                 end do
              end if
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
           do ipoin = 1,npoin
              if( gisca(ipoin) >= 1 ) kfl_geon2(ipoin) = 11
           end do
           call memgen(3_ip,npoin,0_ip) 
        end if

        !----------------------------------------------------------------------
        !
        ! Copy KFL_GEON2(NPOIN) to KFL_GEONO(NBOPO)
        !
        !----------------------------------------------------------------------

        ierr2 = 0
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              igeon = kfl_geono(ibopo)
              if( kfl_geon2(ipoin) /= igeon ) ierr2 = ierr2 + 1
              kfl_geono(ibopo) = kfl_geon2(ipoin)
           end if
        end do

        !----------------------------------------------------------------------
        !
        ! Parall: SKCOS can be different. Take that of master's node
        !
        !----------------------------------------------------------------------

        call par_sleskc() 

        !----------------------------------------------------------------------
        !
        ! Count conditions
        !
        !----------------------------------------------------------------------

        if( itask == 1 ) then
           icoun = 0
           do ibopo = 1,nbopo
              icode = kfl_geono(ibopo)
              if( icode /= 0 ) icoun(icode) = icoun(icode) + 1
           end do
           ioutp(1) = icoun(10)
           ioutp(2) = icoun(20)
           ioutp(3) = icoun(30)
           ioutp(4) = icoun(31)
           ioutp(5) = icoun(40)
           ioutp(6) = icoun(41)
           ioutp(7) = icoun(50)
           ioutp(8) = icoun(60)
           call outfor(44_ip,lun_outpu,' ')
        end if
        !
        ! Deallocate volatile memory
        !
        call domain_memory_deallocate('BOUNDARY GRAPH') ! Deallocate boundary graph R_BOU and C_BOU
        call memory_deallo(memor_dom,'LEDBO','geonor',ledbo)
        call memory_deallo(memor_dom,'BOUNO','geonor',bouno)
        call memory_deallo(memor_dom,'KPOIN','geonor',kpoin)
        call memory_deallo(memor_dom,'KFL_GEON2','geonor',kfl_geon2)

     end if

     !----------------------------------------------------------------------
     !
     ! Check tangent vector
     !
     !----------------------------------------------------------------------

     call PAR_SUM(jerro)
     if( jerro > 0 ) then
        call runend('GEONOR:SOME TANGENT VECTORS COULD NOT BE COMPUTED')
     end if
     
     !----------------------------------------------------------------------
     !
     ! Check that, if mesh is changing, the type of conditions remains
     ! the same: if not, then the modules would have to change their 
     ! conditions    
     !
     !----------------------------------------------------------------------

     if( itask == 2 ) then
        call PAR_SUM(ierr2)
        if( ierr2 /= 0 ) then
           call runend('GEONOR: CONDITIONS HAVE CHANGED')
        end if
     end if

  end if

end subroutine geonor
