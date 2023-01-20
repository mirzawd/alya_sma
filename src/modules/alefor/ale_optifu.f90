!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_optifu()
  !-----------------------------------------------------------------------
  !****f* domain/ale_optifu
  ! NAME
  !    domain
  ! DESCRIPTION
  !    Laplacian ale_optifuing using Gauss-Seidel
  ! USED BY
  !    Turnon 
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_elmtyp
  use def_parame
  use def_elmtyp
  use def_alefor
  use def_domain
  use def_alefor
  use mod_memchk
  implicit none
  integer(ip)           :: ipoin,inode,izdom,telec,kk,cnode,kelem
  integer(ip)           :: ii,kpoin,iters
  integer(ip), pointer  :: conec(:,:)
  integer(4)            :: istat

  if( INOTMASTER ) then

     telec = 0
     do ipoin = 1,npoin
        telec = max(telec,pelpo_2(ipoin+1)-pelpo_2(ipoin))
     end do     
     allocate(conec(mnode,telec),stat=istat)
     call memchk(zero,istat,memor_dom,'CONEC','dod_optifu',conec)

     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then
           kfl_fixno_ale(1,ipoin) = 1
        else
           kfl_fixno_ale(1,ipoin) = 0
        end if
        if( lnoch(ipoin) == NOHOL ) kfl_fixno_ale(1,ipoin) = 1 
     end do

     do iters = 1,nsmoo_ale
        do ipoin = 1,npoin !global

           if( kfl_fixno_ale(1,ipoin) == 0 )then
              telec = 0
              do izdom = pelpo_2(ipoin), pelpo_2(ipoin+1)-1
                 kelem = lelpo_2(izdom)
                 if( kelem <= nelem ) then
                    cnode = nnode(ltype(kelem))-1 
                    if( ltype(kelem) > 0 ) telec = telec + 1
                 end if
              end do

              kk = 0 
              do izdom = pelpo_2(ipoin), pelpo_2(ipoin+1)-1
                 kelem = lelpo_2(izdom)
                 if( kelem <= nelem ) then
                    if( ltype(kelem) > 0 )then
                       kk = kk + 1
                       inode = 0
                       do ii = 1,nnode(ltype(kelem))
                          kpoin = lnods(ii,kelem)
                          if( kpoin /= ipoin ) then
                             inode = inode + 1
                             conec ( inode ,kk ) = kpoin
                          end if
                       end do
                    end if
                 end if
              end do

              call ale_optele ( ipoin, cnode , telec ,conec )

           end if

        end do
     end do

     call memchk(two,istat,memor_dom,'CONEC','dod_optifu',conec)
     deallocate(conec,stat=istat)
     if(istat/=0) call memerr(two,'CONEC','dod_optifu',0_ip)

  end if

end subroutine ale_optifu

subroutine ale_optele ( x , k , M , neigh )

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  implicit none  
  integer(ip),intent(in)         :: x , k , M
  integer(ip),intent(in)         :: neigh(mnode,*)
  integer(ip)                    :: idime,tetra,maxiter,iter,ipoin,cont
  real(rp)                       :: normgradiF,efici_min,efici_max,efici_med,efici_med_old
  real(rp)                       :: kappaS(M),kappaS_ini(M),kappaS_iter(M),efici(M)
  real(rp)                       :: gradiF(ndime),eps,coord_new(ndime),coord_ini(ndime),coord_aux(ndime)

  !print*,'nuevopuntooooooooo'

  normgradiF = 10000.0_rp
  maxiter    = 1000_ip
  iter       = 0_ip
  eps        = 0.001_rp
  ipoin      = x
  !print*,'coords inicialessssssssss',coord_ale(1:3,ipoin)
  do idime = 1,ndime
     coord_ini(idime) = coord_ale(idime,ipoin,1)
     coord_new(idime) = coord_ale(idime,ipoin,1) !!!lo inicializo al original!!!!???
  end do
  efici_max=-10000.0_rp
  efici_min=10000.0_rp
  efici_med = 0.0_rp
  efici_med_old=-10000.0_rp

  do while( normgradiF > eps .and. iter < maxiter )
     !print*,'-----------------------------------llamo a evagra con iter=',iter
     !print*,'entro a evagra'

     call ale_evagra( x , k ,M, neigh,gradiF,normgradiF,kappaS)
     !print*,'salgo de evagra'

    if(iter==0)then
        do tetra=1,M
           kappaS_ini(tetra) =kappaS (tetra)
           !print*,kappaS_ini(tetra)
        end do
     else
        efici_med = 0.0_rp
        do tetra=1,M
           kappaS_iter(tetra) = kappaS (tetra)
           efici(tetra) =(kappaS_ini(tetra) - kappaS_iter(tetra))/kappaS_ini(tetra)
           efici_max = max(efici_max,efici(tetra))
           efici_min = min(efici_min,efici(tetra))

           !print*,efici(tetra)
           efici_med = efici_med + efici(tetra)
        end do
        efici_med = efici_med/M
        if( efici_med_old>=efici_med)then
           !print*,'PAROOOOOOOOOO-efici_med_old'
           iter = maxiter
           normgradiF = 0.0_rp
        else
           efici_med_old=efici_med
        end if
     end if

     ! normgradiF = 0.0_rp
     ! do idime=1,ndime
     !    normgradiF  = normgradiF + gradiF(idime)*gradiF(idime)
     ! end do
     !print*,gradiF
     ! normgradiF = sqrt(normgradiF)
     !print*,'NORMGRADIFFFFFFFFFFF-PARAMOS?',normgradiF 
     if(normgradiF < eps)then
        !print*,'PAROOOOOOOOOO-normagradif',iter
        iter = maxiter
     else
        do idime=1,ndime
           coord_aux(idime) = coord_new(idime)
        end do
        ! print*,'llamo a evalfa'

        call ale_evalfa( x ,k,M,neigh ,gradiF,coord_new)

        do idime=1,ndime
           coord_ale(idime,ipoin,1) = coord_new(idime)
        end do
        !print*,'coordssssssss para la iter',iter
        !print*,'es igual a=',coord_ale(1:3,ipoin)
        !print*,coord_ale(1:3,ipoin)
        cont = 0
        do idime=1,ndime
           if(coord_aux(idime)==coord_new(idime))cont =  1 !!!???
        end do
        if(cont ==ndime)then
           iter = maxiter
           !print*,'entro en lo de cont'
        end if
        iter = iter + 1
        !print*,iter
     end if

  end do
end subroutine ale_optele
