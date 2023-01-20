!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmust()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmust
  ! NAME 
  !    tur_elmust
  ! DESCRIPTION
  !    Compute the generalized distance to the wall via a 
  !    Poisson equation:
  !    1. Solve Lapl(f)=-1, with f=0 on wall
  !    2. d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !    See the following references:
  !    P.G. Tucker, Differential equation-based wall distance computation for
  !         DES and RANS, J. Comp. Phys. 190 (2003) 229-248.
  !    P.G. Tucker, Int. J. Numer. Fluids 33 (2000) 869.
  !    P.G. Tucker, Appl. Math. Model. 22 (1998) 293.
  ! USES
  ! USED BY
  !    Turbul
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_gradie
  use mod_memchk
  use mod_postpr
  implicit none
  integer(ip)  :: ielem,igaus,idime,inode,ipoin,pnode,pgaus,pelty
  integer(ip)  :: izmat,izrhs,itotn
  integer(4)   :: istat
  real(rp)     :: elmal(mnode,mnode),elrhs(ndime,mnode)
  real(rp)     :: elcod(ndime,mnode),elvel(ndime,mnode)
  real(rp)     :: gpvel(ndime,mnode)
  real(rp)     :: gpcar(ndime,mnode,mgaus),gpvol(mgaus),gpdet
  real(rp)     :: xjaci(9),xjacm(9),usmax,gpdif
  character(5) :: wnada(2)

  allocate(ustar_tur(npoin),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'USTAR_TUR','tur_memarr',ustar_tur)

  if(kfl_walgo_tur==1.and.ittot_tur==0) then

     if(kfl_paral/=0) then
        do izmat=1,solve(1)%nzmat
           amatr(izmat)=0.0_rp
        end do
        do izrhs=1,2_ip*solve(1)%nzrhs
           rhsid(izrhs)=0.0_rp
        end do        
        call memgen(zero,ndime,npoin)
        do ielem=1,nelem
           !
           ! Element properties and dimensions
           !
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           pgaus=ngaus(pelty)
           !
           ! Gather operations: ELCOD
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
              end do
           end do
           !
           ! 1st order Cartesian derivatives GPCAR and GPVOL=dV=|J|*wgx
           !
           do igaus=1,pgaus     
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(:,:,igaus),& 
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
              gpvol(igaus)=elmar(pelty)%weigp(igaus)*gpdet  
           end do
           !
           ! Compute element matrix ELMAL and assemble LAPLA_TUR
           !
           call tur_elmla2(&
                pnode,pgaus,lnods(:,ielem),gpcar,elmar(pelty)%shape,&
                gpvol,elmal,elrhs)
           call assmat(&
                1_ip,pnode,pnode,nunkn_tur,solve(1)%kfl_algso,&
                ielem,lnods(:,ielem),elmal,amatr)
           call assrhs(-ndime,pnode,lnods(:,ielem),elrhs,rhsid)

        end do
     end if
     !
     ! Solve system: Lapl(f)=-1, with f=0 on wall
     !
     call tur_inisol(1_ip)
     itotn=1
     do idime=1,ndime
        do izrhs=1,solve(1)%nzrhs
           unkno(izrhs)=0.0_rp
        end do
        call solver(rhsid(itotn:),unkno,amatr,pmatr) 
        do ipoin=1,npoin
           gevec(idime,ipoin)=unkno(ipoin)
        end do
        itotn=itotn+npoin
     end do
     wnada(1) = 'GEVEC'
     wnada(2) = 'VECTO'
     call postpr(gevec,wnada,ittim,cutim)
     !
     ! Compute maximum advection
     !
     if(kfl_paral/=0) then
        usmax=0.0_rp
        do ipoin=1,npoin
           if(kfl_fixno_tur(1,ipoin,1)==3.or.kfl_fixno_tur(1,ipoin,1)==4) then
              if(ustar_tur(ipoin)>usmax) usmax=ustar_tur(ipoin)
           end if
        end do
        gpdif=1.0e-6_rp*usmax
     end if
     !
     ! Advect ustar
     !
     if(kfl_paral/=0) then
        do izmat=1,solve(1)%nzmat
           amatr(izmat)=0.0_rp
        end do
        do izrhs=1,solve(1)%nzrhs
           rhsid(izrhs)=0.0_rp
        end do        
        do ipoin=1,npoin
           unkno(ipoin)=ustar_tur(ipoin)
        end do        
        call memgen(zero,npoin,zero)
        do ielem=1,nelem
           !
           ! Element properties and dimensions
           !
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           pgaus=ngaus(pelty)
           !
           ! Gather operations: ELCOD
           !
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)
                 elvel(idime,inode)=gevec(idime,ipoin)
              end do
           end do
           !
           ! 1st order Cartesian derivatives GPCAR and GPVOL=dV=|J|*wgx
           !
           do igaus=1,pgaus     
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(:,:,igaus),& 
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
              gpvol(igaus)=elmar(pelty)%weigp(igaus)*gpdet  
           end do
           !
           ! Compute element matrix ELMAL and assemble LAPLA_TUR
           !
           call tur_elmla3(&
                pnode,pgaus,lnods(:,ielem),elvel,gpcar,elmar(pelty)%shape,&
                gpvol,gpvel,gpdif,elmal,elrhs)
           call assmat(&
                1_ip,pnode,pnode,nunkn_tur,solve(1)%kfl_algso,&
                ielem,lnods(:,ielem),elmal,amatr)
           call assrhs(1_ip,pnode,lnods(:,ielem),elrhs,rhsid)

        end do
        call solver(rhsid,unkno,amatr,pmatr) 
        do ipoin=1,npoin
           ustar_tur(ipoin)=unkno(ipoin)
        end do  
        wnada(1) = 'GESCA'
        wnada(2) = 'SCALA'
        call postpr(ustar_tur,wnada,ittim,cutim)      
        stop
     end if

     
  end if

end subroutine tur_elmust
