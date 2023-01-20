!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_inicoa()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_solgro
  ! NAME 
  !    nsi_solcoa
  ! DESCRIPTION
  !    Initial solution using coarse grid agglomeration
  ! USED BY
  !    nsi_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_master
  use def_domain
  use def_nastin 
  use mod_nsi_schur_operations,  only : nsi_solini
  implicit none
  real(rp) :: dummr(2)

  if( IMASTER ) then
     call nsi_solini(3_ip)
     call nsi_matrix()
     call nsi_solcoa(&
          dummr,dummr,dummr,dummr,dummr,dummr)
  else
     dtinv_nsi = 0.0_rp 
     call nsi_updunk(1100_ip)
     call nsi_solini(3_ip)
     call nsi_matrix()
     call nsi_solcoa(&
          amatr(poauu_nsi:),amatr(poaup_nsi:),amatr(poapu_nsi:),amatr(poapp_nsi:),&
          rhsid,rhsid(ndbgs_nsi+1:))
  end if
  
end subroutine nsi_inicoa

subroutine nsi_solcoa(Auu,Aup,Apu,App,bu,bp)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_solcoa
  ! NAME 
  !    nsi_solcoa
  ! DESCRIPTION
  !    This routine sovles the coarse grid system
  ! USED BY
  !    nsi_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_master
  use def_domain
  use def_nastin
  use mod_memchk
  use mod_nsi_schur_operations
  use mod_messages,              only : livinf
  use mod_local_basis,           only : local_basis_global_to_local
  use mod_local_basis,           only : local_basis_local_to_global
  use mod_communications_global, only : PAR_SUM
  use mod_communications_global, only : PAR_MIN
  use mod_skyline,               only : lufact
  use mod_skyline,               only : lusolv
  implicit none
  real(rp),    intent(in) :: Auu(ndime,ndime,*)
  real(rp),    intent(in) :: Aup(ndime,*) 
  real(rp),    intent(in) :: Apu(ndime,*)
  real(rp),    intent(in) :: App(*)
  real(rp),    intent(in) :: bu(*)
  real(rp),    intent(in) :: bp(*)
  integer(ip)             :: igrou,ipoin,kskyl,kgrou,idofn,jdofn,jpoin
  integer(ip)             :: igrou1,jgrou1,idime,jdime,izdom,info,jgrou
  integer(ip)             :: mgrou,ipoin1,ii,ibopo,jj,jbopo,ifixx,jfixx
  integer(ip)             :: ngrou,ndofn,nskyl
  integer(4)              :: istat
  real(rp)                :: alpha
  integer(ip), pointer    :: iskyl(:),idiag(:),lgrou_nsi(:),gigro(:)
  real(rp),    pointer    :: Ac(:),mu(:)
  real(rp),    pointer    :: uu(:),rr(:),ww(:)
  integer(ip), pointer    :: kmatr(:)
  !
  ! Initialization
  !
  if( ngrou_dom <= 0 ) &
       call runend('NSI_INICO: CANNOT SOLVE COARSE GRID SYSTEM: DEFINE GROUPS FIRST')
  call livinf(59_ip,'SOLVE COARSE GRID SYSTEM',0_ip)
  ndofn = ndime + 1
  ngrou = ngrou_dom 
  !
  ! Parallel exchange
  !
  if( INOTMASTER ) then
     call rhsmod(ndime,bu)
     call rhsmod( 1_ip,bp)
  end if
  !
  ! Rotate global to local
  !
  if( INOTEMPTY ) call local_basis_global_to_local(kfl_fixrs_nsi,veloc,LAST_COMPONENT=1_ip)
  !
  ! Define groups
  !
  if( INOTMASTER ) then
     allocate(lgrou_nsi(npoin),stat=istat)
     call memchk(0_ip,istat,mem_modul(1:2,modul),'LGROU_NSI','nsi_solcoa',lgrou_nsi)
     do ipoin = 1,npoin
        lgrou_nsi(ipoin) = lgrou_dom(ipoin)
     end do
  end if
  !
  ! Redefine groups if some are empty
  !
  allocate(gigro(ngrou))
  do igrou = 1,ngrou
     gigro(igrou) = 0
  end do
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        igrou = lgrou_nsi(ipoin)  
        if( igrou > 0 ) gigro(igrou) = gigro(igrou) + 1
     end do
  end if
  call PAR_SUM(ngrou,gigro)
  jgrou = 0
  do igrou = 1,ngrou
     if( gigro(igrou) > 0 ) then
        jgrou = jgrou + 1
        gigro(igrou) = jgrou
     end if
  end do
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        igrou = lgrou_nsi(ipoin)
        if( igrou > 0 ) lgrou_nsi(ipoin) = gigro(igrou)
     end do
  end if
  ngrou = jgrou
  deallocate(gigro)
  !
  ! Allocate memory to depict zero lines in matrix
  !
  allocate( kmatr(ngrou*ndofn) )
  do igrou1 = 1,ngrou*ndofn
     kmatr(igrou1) = 0
  end do
  !
  ! Allocate memory for system
  !
  allocate(iskyl(ndofn*ngrou+1),stat=istat)
  call memchk(0_ip,istat,mem_modul(1:2,modul),'ISKYL','nsi_solcoa',iskyl)
  if( INOTMASTER ) then
     do igrou = 1,ngrou*ndofn+1
        iskyl(igrou) = ngrou*ndofn
     end do
     do ipoin = 1,npoin
        igrou = lgrou_nsi(ipoin)
        if( igrou > 0 ) then
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin = c_dom(izdom)           
              jgrou = lgrou_nsi(jpoin)
              if( jgrou > 0 .and. igrou >= jgrou ) then
                 do idofn = 1,ndofn 
                    kgrou = (igrou-1) * ndofn + idofn + 1  
                    do jdofn = 1,ndofn
                       mgrou = (jgrou-1) * ndofn + jdofn
                       if( mgrou < iskyl(kgrou) ) iskyl(kgrou) = mgrou
                    end do
                 end do
              end if
           end do
        end if
     end do
  end if
  call PAR_MIN(ndofn*ngrou+1,iskyl)
  nskyl    =  1
  iskyl(1) =  1
  !
  ! Skyline: velocity
  !
  allocate(idiag(ngrou*ndofn),stat=istat)
  call memchk(0_ip,istat,mem_modul(1:2,modul),'IDIAG','cregrp',idiag)
  kgrou = 0_ip
  do igrou = 1,ngrou
     do idofn = 1,ndofn 
        kgrou          = kgrou + 1
        kskyl          = kgrou - iskyl(kgrou+1)
        idiag(kgrou)   = nskyl + kskyl
        kskyl          = 2 * kskyl + 1  
        nskyl          = nskyl + kskyl
        iskyl(kgrou+1) = nskyl
     enddo
  end do
  nskyl = nskyl - 1 
  !
  ! Allocate memory for system of equations
  !
  allocate(mu(ngrou*ndofn),stat=istat)
  call memchk(0_ip,istat,mem_modul(1:2,modul),'MU','nsi_solcoa',mu)

  allocate(Ac(nskyl),stat=istat)
  call memchk(0_ip,istat,mem_modul(1:2,modul),'AC','nsi_solcoa',Ac)
  !
  ! Assemble matrix by agglomeration
  !
  if( INOTMASTER ) then

     do ipoin = 1,npoin 
        igrou = lgrou_nsi(ipoin)
        do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
           jpoin = c_dom(izdom)
           jgrou = lgrou_nsi(jpoin)

           do idime = 1,ndime
              if( kfl_fixno_nsi(idime,ipoin) == 0 ) then
                 !
                 ! Auu
                 !
                 do jdime = 1,ndime 
                    if( kfl_fixno_nsi(jdime,jpoin) == 0 ) then
                       igrou1 = (igrou-1) * ndofn + idime
                       jgrou1 = (jgrou-1) * ndofn + jdime
                       kmatr(igrou1) = 1
                       if( igrou1 < jgrou1 ) then
                          kskyl       = iskyl(jgrou1+1) - (jgrou1-igrou1)
                          Ac(kskyl) = Ac(kskyl)     + Auu(jdime,idime,izdom)
                       else     
                          kskyl       = idiag(igrou1)   - (igrou1-jgrou1)
                          Ac(kskyl) = Ac(kskyl)     + Auu(jdime,idime,izdom)
                       end if
                    end if
                 end do
                 !
                 ! Aup
                 !
                 jdime = ndofn
                 jbopo = lpoty(jpoin)
                 jfixx = 0
                 if( jbopo > 0 ) then
                    jfixx = kfl_fixpr_nsi(1,jpoin)
                 end if
                 if( jfixx == 0 ) then
                    igrou1 = (igrou-1) * ndofn + idime
                    jgrou1 = (jgrou-1) * ndofn + jdime
                    kmatr(igrou1) = 1
                    if( igrou1 < jgrou1 ) then
                       kskyl       = iskyl(jgrou1+1) - (jgrou1-igrou1)
                       Ac(kskyl) = Ac(kskyl)     + Aup(idime,izdom)
                    else     
                       kskyl       = idiag(igrou1)   - (igrou1-jgrou1)
                       Ac(kskyl) = Ac(kskyl)     + Aup(idime,izdom)
                    end if
                 end if
              end if
           end do

           idime = ndofn
           ibopo = lpoty(ipoin)
           ifixx = 0
           if( ibopo > 0 ) then
              ifixx = kfl_fixpr_nsi(1,ipoin)
           end if

           if( ifixx == 0 ) then
              !
              ! Apu
              !
              igrou1 = (igrou-1) * ndofn + idime
              do jdime = 1,ndime   
                 if( kfl_fixno_nsi(jdime,jpoin) == 0 ) then
                    jgrou1 = (jgrou-1) * ndofn + jdime
                    kmatr(igrou1) = 1
                    if( igrou1 < jgrou1 ) then
                       kskyl       = iskyl(jgrou1+1) - (jgrou1-igrou1)
                       Ac(kskyl) = Ac(kskyl)     + Apu(jdime,izdom)
                    else     
                       kskyl       = idiag(igrou1)   - (igrou1-jgrou1)
                       Ac(kskyl) = Ac(kskyl)     + Apu(jdime,izdom)
                    end if
                 end if
              end do
              !
              ! App
              !
              jdime = ndofn
              jbopo = lpoty(jpoin)
              jfixx = 0
              if( jbopo > 0 ) then
                 jfixx = kfl_fixpr_nsi(1,jpoin)
              end if
              if( jfixx == 0 ) then
                 jgrou1 = (jgrou-1) * ndofn + jdime
                 kmatr(igrou1) = 1
                 if( igrou1 < jgrou1 ) then
                    kskyl       = iskyl(jgrou1+1) - (jgrou1-igrou1)
                    Ac(kskyl) = Ac(kskyl)     + App(izdom)
                 else     
                    kskyl       = idiag(igrou1)   - (igrou1-jgrou1)
                    Ac(kskyl) = Ac(kskyl)     + App(izdom)
                 end if
              end if
           end if
        end do
     end do
  end if
  call PAR_SUM(nskyl,Ac)
  !
  ! Check matrix can be inverted
  !
  call PAR_SUM(ngrou*ndofn,kmatr)
  do igrou1 = 1,ngrou*ndofn
     if( kmatr(igrou1) == 0 ) then
        kskyl     = idiag(igrou1)
        Ac(kskyl) = 1.0_rp
     end if
  end do
  deallocate( kmatr ) 
  !
  ! Factorize matrix
  !
  if( INOTMASTER ) then
     call lufact(ngrou*ndofn,nskyl,iskyl,Ac,idiag,info)
     if( info /= 0 ) call runend('MATGRU: ERROR WHILE DOING CHOLESKY FACTORIZATION')
  end if
  !
  ! Initial solution
  !
  if( INOTMASTER ) then
     allocate( uu(ndofn*npoin) )
     allocate( ww(ndime*npoin) )
     allocate( rr(ndofn*npoin) )
  end if

  if( INOTMASTER ) then
     !
     ! Equation residual (Momentum): r = b - A x
     !
     call nsi_auuvec(1_ip,Auu,veloc(:,:,1),ww)
     call nsi_aupvec(2_ip,Aup,press(:,1),ww)     
     call rhsmod(ndime,ww)
     do ipoin = 1,npoin
        ii = ( ipoin - 1 ) * ndime
        jj = ( ipoin - 1 ) * ndofn
        do idime = 1,ndime
           ii     = ii + 1
           jj     = jj + 1
           rr(jj) = bu(ii) - ww(ii) 
        end do
     end do
     !
     ! Equation residual (Schur)
     !
     call nsi_apuvec(1_ip,Apu,veloc(:,:,1),ww)
     call nsi_appvec(2_ip,App,press(:,1),ww)
     call rhsmod(1_ip ,ww)
     do ipoin = 1,npoin
        jj     = ipoin * ndofn
        rr(jj) = bp(ipoin) - ww(ipoin)
     end do
     !
     ! Compute coarse residual r'= W^T ( b - A x )
     !
     do ipoin = 1,npoin_own
        igrou  = lgrou_nsi(ipoin)
        igrou1 = ( igrou - 1 ) * ndofn
        ipoin1 = ( ipoin - 1 ) * ndofn
        do idime = 1,ndofn
           igrou1 = igrou1 + 1
           ipoin1 = ipoin1 + 1
           mu(igrou1) = mu(igrou1) + rr(ipoin1)
        end do
     end do
  end if
  call PAR_SUM(ndofn*ngrou,mu)
 
  if( INOTMASTER ) then
     !
     ! Solve system: A' x' = r'
     !
     info = 0
     call LUsolv(&                
          ngrou*ndofn,nskyl,iskyl,1_ip,Ac,&
          mu,ngrou*ndofn,idiag,info)     
  end if
 
  call PAR_SUM(info)
  if( info /= 0 ) call runend('MATGRU: COULD NOT SOLVE INITR_DOML SYSTEM')

  if( INOTMASTER ) then
     !
     ! Update residual r = W r'
     !
     do ipoin = 1,npoin
        igrou  = lgrou_nsi(ipoin)
        igrou1 = (igrou-1) * ndofn
        ipoin1 = (ipoin-1) * ndofn
        do idime  = 1,ndime
           igrou1 = igrou1 + 1 
           ipoin1 = ipoin1 + 1
           if( kfl_fixno_nsi(idime,ipoin) == 0 ) then
              rr(ipoin1) = mu(igrou1)
           end if
        end do
        ibopo = lpoty(ipoin)
        ifixx = 0
        if( ibopo > 0 ) then
           ifixx = kfl_fixpr_nsi(1,ipoin)
        end if
        if( ifixx == 0 ) then
           igrou1 = igrou1 + 1
           ipoin1 = ipoin1 + 1
           rr(ipoin1) = mu(igrou1) 
        end if
     end do
  end if
  !
  ! Acceleration: alpha = < r , A z > // || A z ||^2
  ! z = ( W^t Ac^-1 W ) r
  ! 
  !if( INOTMASTER ) call nsi_allvec(1_ip,Auu,Aup,App,Apu,rr,uu)
  !call prodxy(ndofn,npoin,uu,vv,alpha) ! vv = rr before solving coarse system
  !call norm2x(ndofn,uu,denom)
  !alpha = alpha / ( denom  * denom )
  
  alpha = 1.0_rp

  if( INOTMASTER ) then
     !
     ! Update solution
     !
     idofn = 0
     do ipoin = 1,npoin
        idofn = ( ipoin - 1 ) * ndofn
        do idime = 1,ndime
           idofn = idofn + 1
           if( kfl_fixno_nsi(idime,ipoin) == 0 ) then
              veloc(idime,ipoin,1) = veloc(idime,ipoin,1) + alpha * rr(idofn)
           end if
        end do
        idofn = idofn + 1
        press(ipoin,1) = press(ipoin,1) + alpha * rr(idofn)
     end do
  end if
  !
  ! Rotate local to global
  !
  if( INOTEMPTY ) call local_basis_local_to_global(kfl_fixrs_nsi,veloc,LAST_COMPONENT=1_ip)
  !
  ! Deallocate memory
  !
  if( INOTMASTER ) then

     call memchk(two,istat,mem_modul(1:2,modul),'LGROU_NSI','nsi_inicoa',lgrou_nsi)
     deallocate(lgrou_nsi,stat=istat)
     if( istat /= 0 ) call memerr(two,'LGROU_NSI','nsi_inicoa',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'RR','nsi_inicoa',rr)
     deallocate(rr,stat=istat)
     if( istat /= 0 ) call memerr(two,'RR','nsi_inicoa',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'WW','nsi_inicoa',ww)
     deallocate(ww,stat=istat)
     if( istat /= 0 ) call memerr(two,'WW','nsi_inicoa',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'UU','nsi_inicoa',uu)
     deallocate(uu,stat=istat)
     if( istat /= 0 ) call memerr(two,'UU','nsi_inicoa',0_ip)

  end if

  call memchk(2_ip,istat,mem_modul(1:2,modul),'AC','nsi_solcoa',Ac)
  deallocate(Ac,stat=istat)
  if(istat/=0) call memerr(2_ip,'AC','nsi_solcoa',0_ip)

  call memchk(2_ip,istat,mem_modul(1:2,modul),'MU','nsi_solcoa',mu)
  deallocate(mu,stat=istat)
  if(istat/=0) call memerr(2_ip,'MU','nsi_solcoa',0_ip)

  call memchk(2_ip,istat,mem_modul(1:2,modul),'IDIAG','nsi_solcoa',idiag)
  deallocate(idiag,stat=istat)
  if(istat/=0) call memerr(2_ip,'IDIAG','nsi_solcoa',0_ip)

  call memchk(2_ip,istat,mem_modul(1:2,modul),'ISKYL','nsi_solcoa',iskyl)
  deallocate(iskyl,stat=istat)
  if(istat/=0) call memerr(2_ip,'ISKYL','nsi_solcoa',0_ip)

end subroutine nsi_solcoa

