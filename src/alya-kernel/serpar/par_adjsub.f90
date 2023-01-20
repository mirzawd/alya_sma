!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_adjsub()
  !-----------------------------------------------------------------------
  !****f* Parall/par_adjsub
  ! NAME
  !    Send and receive adjacent graph from adjacent subdomains
  ! DESCRIPTION
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  implicit none
  integer(ip) :: dom_i

  if( INOTSLAVE ) return

     !do ii = 1,nneig

        !dom_i = commd % neights(ii)
        !ini   = commd % bound_size(ii)
        !bsize = commd % bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
        if( kfl_paral > dom_i ) then
           !call par_sengra(ii,dom_i)
           !call par_recgra(ii,dom_i)
        else
           !call par_recgra(ii,dom_i)
           !call par_sengra(ii,dom_i)           
        end if
        !
        ! Send my additional graph from subdomain DOM_I
        !

        !
        ! Receive additional graph from subdomain DOM_I
        !
        
        !bsize4 = int(bsize,4)
        !   call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
        !        PAR_INTEGER,  dom_i, 0_4,     &
        !        loc_rpari1(ini:), bsize4,              &
        !        PAR_INTEGER, dom_i, 0_4,      &
        !        PAR_COMM_MYCODE, status, istat )
#endif
     !end do

end subroutine par_adjsub

subroutine par_sengra(ii,dom_i)
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  implicit none
  integer(ip),  intent(in) :: ii,dom_i
  !
  ! KFL_PARALL sends to DOM_I
  !
  !ini     = commd % bound_size(ii)
  !bsize   = commd % bound_size(ii+1) - ini
  !npoin_s = 0
  !nzdom_s = 0
  !kpoin   = 0
  !
  ! Look for neighbors of adjacent nodes
  !
  !allocate( lnode_s(npoin), stat=istat )
  !do ipoin = 1,npoin
  !   lnode_s(ipoin) = 0
  !end do  
  !do jj = 1,commd % bound_dim
  !   ipoin   = commd % bound_perm(jj)
  !   nzdom_s = nzdom_s + (r_dom(ipoin+1)-r_dom(ipoin))
  !   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
  !      jpoin = c_dom(izdom)
  !      if( lnode_s(jpoin) == 0 ) then
  !         npoin_s        = npoin_s + 1
  !         lnode_s(jpoin) = npoin_s
  !      end do
  !   end do
  !end do
  !
  ! Inverse permutation
  !
  !allocate( liode_s(npoin_s), stat=istat )
  !do ipoin = 1,npoin
  !   kpoin = lnode_s(ipoin)
  !   if( kpoin /= 0 ) then
  !      liode_s(kpoin) = ipoin
  !   end if
  !end do
  !
  ! Fill in coordinates
  !
  !allocate( coord_s(ndime,npoin_s), stat=istat )   
  !do ipoin = 1,npoin
  !   kpoin = lnode_s(ipoin)
  !   if( kpoin /= 0 ) then
  !      do idime = 1,ndime
  !         coord_s(idime,kpoin) = coord(idime,ipoin)
  !      end do
  !   end do
  !end do
  !
  ! Fill in graph R_DOM_S
  !
  !allocate( r_dom_s(npoin_s+1),     stat=istat )

  !do jj = 1,commd % bound_dim
  !   ipoin          = commd % bound_perm(jj)
  !   kpoin          = lnode_s(ipoin)
  !   r_dom_s(kpoin) = r_dom(ipoin+1)-r_dom(ipoin)
  !end do

  !do kpoin = 1,npoin_s
  !   r_dom_s(kpoin+1) = r_dom_s(kpoin) + r_dom_s(kpoin+1)
  !end do
  !r_dom_s(1) = 1
  !nzdom_s    = r_dom_s(npoin_s+1)-1
  !
  ! Fill in graph C_DOM_S
  !
  !allocate( c_dom_s(nzdom_s), stat=istat )

  !do jj = 1,commd % bound_dim

  !   ipoin = commd % bound_perm(jj)
  !   kpoin = liode_s(ipoin)
  !   kzdom = r_dom(ipoin+1)-r_dom(ipoin)

  !   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
  !      nzdom_s = nzdom_s + 1
  !      jpoin   = c_dom(izdom)
  !      kpoin   = lnode_s(jpoin)
  !      !r_dom_s(kpoin) = 
  !      if( lnode_s(jpoin) /= 0 ) then
  !         npoin_s        = npoin_s + 1
  !         lnode_s(jpoin) = npoin_s
  !      end if
  !   end do

  !   r_dom_s(kpoin+1) = r_dom_s(kpoin+1) kzdom

  !end do


end subroutine par_sengra
