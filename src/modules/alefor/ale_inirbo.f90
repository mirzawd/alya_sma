!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_inirbo
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_inirbo
  ! NAME
  !    ale_inirbo
  ! DESCRIPTION
  !    This routine is similar to ibm_iniunk but for the case where rigid body is treated inside ale
  !    This routines computes for each particle:
  !    1. The volume 
  !    2. The center of gravity
  !    3. The inertia tensor
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_elmtyp
  use mod_memchk
  use def_alefor
  use mod_messages, only : livinf
  use mod_communications_global, only : PAR_SUM
  use mod_bouder

  implicit none
  integer(ip)          :: iboib,igaub,pblty,pgaub,iimbo,ipoib,inoib,isize
  integer(ip)          :: idime,iiner,nsize,pnodb,ipoin,kpoin
  integer(4)           :: istat
  real(rp)             :: baloc(ndime,ndime),bocod(ndime,max(mnoib,mnodb)),xfact
  real(rp)             :: gbsur,eucta,bouno(3),gpcib(3),rx,ry,rz,nx,ny,nz
  real(rp)             :: rx3,ry3,rz3,r3,onov3,onovn
  real(rp),    pointer :: inerc(:,:)
  integer(ip), pointer :: kfl_calgr(:)
  character(20)        :: messa

  !----------------------------------------------------------------------
  !
  ! Initial coordinates
  !
  !---------------------------------------------------------------------- 
  
  do iimbo = 1,nrbod
     do kpoin = 1,rbbou(iimbo) % npoib
        ipoin = rbbou(iimbo) % lninv(kpoin) 
        do idime = 1,ndime
           rbbou(iimbo) % cooin(idime,kpoin) = coord_ori(idime,ipoin)
           rbbou(iimbo) % cooib(idime,kpoin) = coord(idime,ipoin)
        end do
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------
  
  call livinf(98_ip,' ',0_ip)
  
  nsize = 5 * ndime - 9                                           ! =1 in 2D, 6 in 3D
  onov3 = 1.0_rp / 3.0_rp
  onovn = 1.0_rp / real(ndime,rp)
  allocate(inerc(nsize,nrbod),stat=istat)
  call memchk(zero,istat,memor_dom,'INERC','ale_inirbo',inerc)
  allocate(kfl_calgr(nrbod),stat=istat)
  call memchk(zero,istat,memor_dom,'CALGR','ale_inirbo',kfl_calgr)
  do iimbo = 1,nrbod   ! flag to decide if posgr is calculated or read, it is set based on the initial posgr
     kfl_calgr(iimbo) = 1_ip
  end do

  do iimbo = 1,nrbod
     if(rbbou(iimbo) % posgr(1) < -0.5e12_rp) then
        do idime = 1,ndime
           rbbou(iimbo) % posgr(idime) = 0.0_rp
        end do
     else
        kfl_calgr(iimbo) = 0_ip
     end if
     rbbou(iimbo) % volum = 0.0_rp
  end do

  !----------------------------------------------------------------------
  !
  ! Center of gravity and volume
  !
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod
     
     if( rbbou(iimbo) % volum == 0.0_rp ) then
        
        do iboib = 1,rbbou(iimbo) % nboib
           pblty = rbbou(iimbo) % ltyib(iboib) 
           pnodb = nnode(pblty)
           pgaub = ngaib(pblty)
           !
           ! BOCOD: Gather
           !
           do inoib = 1,pnodb
              ipoib = rbbou(iimbo) % lnoib(inoib,iboib)
              do idime = 1,ndime
                 bocod(idime,inoib) = rbbou(iimbo) % cooib(idime,ipoib)                 
              end do
           end do
           !
           ! Loop over Gauss points
           !
           do igaub = 1,pgaub
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty) % derib(:,:,igaub),&  
                   bocod,baloc,eucta)                 
              gbsur = elmar(pblty) % weiib(igaub)*eucta 
              !
              ! GPCIB: Coordinates of Gauss point
              !
              do idime = 1,ndime
                 gpcib(idime) = 0.0_rp
                 do inoib = 1,pnodb
                    gpcib(idime) = gpcib(idime) + &
                         bocod(idime,inoib) * elmar(pblty) % shaib(inoib,igaub)
                 end do
              end do
              !
              ! BOUNO: Exterior normal
              !
              !bouno(1:ndime) = baloc(1:ndime,ndime)
              call exterb(iimbo,pnodb,iboib,bouno)
              !
              ! Center of gravity: Integrate xi^2.ni
              ! Volume: Integrate xi.ni
              !
              do idime = 1,ndime
                 xfact                       = gbsur * gpcib(idime) * bouno(idime)
                 if (kfl_calgr(iimbo) == 1_ip)   rbbou(iimbo) % posgr(idime) = rbbou(iimbo) % posgr(idime) + xfact * gpcib(idime) 
                 rbbou(iimbo) % volum        = rbbou(iimbo) % volum        + xfact
              end do
           end do

        end do

     end if
     
  end do
  !
  ! For body fitted type bodies, sum in parallel
  ! Take absolute value of volume as for body fitted boundaries, the normal points outwards
  !
  do iimbo = 1,nrbod
     if (kfl_calgr(iimbo) == 1_ip) call PAR_SUM(ndime,rbbou(iimbo) % posgr)
     call PAR_SUM(rbbou(iimbo) % volum)
     rbbou(iimbo) % volum = -rbbou(iimbo) % volum
     if (kfl_calgr(iimbo) == 1_ip)  then
        do idime = 1,ndime
           rbbou(iimbo) % posgr(idime) = -rbbou(iimbo) % posgr(idime)
        end do
     end if
  end do
  !
  ! Correct cog and volume
  !
  do iimbo = 1,nrbod
     if (kfl_calgr(iimbo) == 1_ip)  then
        do idime = 1,ndime
           rbbou(iimbo) % posgr(idime) = rbbou(iimbo) % posgr(idime) * 0.5_rp
        end do
     end if
     rbbou(iimbo) % volum = rbbou(iimbo) % volum * onovn
  end do
  !
  ! Set posil = posgr if posil has not been read in ibm.dat. This is the logical choice for body fitted cases
  !
  do iimbo = 1,nrbod
     if(rbbou(iimbo) % posil(1,1) < -0.5e12_rp) then
        do idime = 1,ndime
           rbbou(iimbo) % posil(idime,1) = rbbou(iimbo) % posgr(idime)
           rbbou(iimbo) % posil(idime,2) = rbbou(iimbo) % posgr(idime)
           rbbou(iimbo) % posil(idime,3) = rbbou(iimbo) % posgr(idime)
           rbbou(iimbo) % posil(idime,4) = rbbou(iimbo) % posgr(idime)
        end do
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! Tensor of inertia
  !
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod

     if(rbbou(iimbo) % momin(1) < -0.5_rp) then  ! momin has been initialized to -1.0 to indicate that it must be calculated if 
        ! it is not read 

        do iboib = 1,rbbou(iimbo) % nboib

           pblty = rbbou(iimbo) % ltyib(iboib) 
           pnodb = nnode(pblty)
           pgaub = ngaib(pblty)
           !
           ! BOCOD: Gather
           !
           do inoib = 1,pnodb
              ipoib = rbbou(iimbo) % lnoib(inoib,iboib)
              do idime = 1,ndime
                 bocod(idime,inoib) = rbbou(iimbo) % cooib(idime,ipoib)
              end do
           end do
           !
           ! Loop over Gauss points
           !
           do igaub = 1,pgaub

              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty) % derib(:,:,igaub),&  
                   bocod,baloc,eucta)                 
              gbsur = elmar(pblty) % weiib(igaub)*eucta 
              !
              ! GPCIB: Coordinates of Gauss point
              !
              do idime = 1,ndime
                 gpcib(idime) = 0.0_rp
                 do inoib = 1,pnodb
                    gpcib(idime) = gpcib(idime) + &
                         bocod(idime,inoib) * elmar(pblty) % shaib(inoib,igaub)
                 end do
              end do
              !
              ! BOUNO: Exterior normal
              !
              call exterb(iimbo,pnodb,iboib,bouno)
              !
              ! Inertia tensor
              !
              if ( ndime == 2 ) then

                 rx = gpcib(1) - rbbou(iimbo) % posgr(1)
                 ry = gpcib(2) - rbbou(iimbo) % posgr(2)
                 nx = bouno(1) 
                 ny = bouno(2) 

                 inerc(1,iimbo) = inerc(1,iimbo) + gbsur * ( rx**3 * nx + ry**3 * ny ) * onov3

              else if ( ndime == 3 ) then

                 rx  = gpcib(1) - rbbou(iimbo) % posgr(1)
                 ry  = gpcib(2) - rbbou(iimbo) % posgr(2)
                 rz  = gpcib(3) - rbbou(iimbo) % posgr(3)
                 rx3 = rx * rx * rx * onov3
                 ry3 = ry * ry * ry * onov3
                 rz3 = rz * rz * rz * onov3
                 nx  = bouno(1) * gbsur
                 ny  = bouno(2) * gbsur
                 nz  = bouno(3) * gbsur
                 r3  = rx3 * nx + ry3 * ny + rz3 * nz 

                 inerc(1,iimbo) = inerc(1,iimbo) + ry3 * ny + rz3 * nz     ! I11
                 inerc(2,iimbo) = inerc(2,iimbo) + rx3 * nx + rz3 * nz     ! I22
                 inerc(3,iimbo) = inerc(3,iimbo) + rx3 * nx + ry3 * ny     ! I33
                 inerc(4,iimbo) = inerc(4,iimbo) - rx*rx*ry * 0.5_rp * nx  ! I12
                 inerc(5,iimbo) = inerc(5,iimbo) - rz*rz*rx * 0.5_rp * nz  ! I13
                 inerc(6,iimbo) = inerc(6,iimbo) - ry*ry*rz * 0.5_rp * ny  ! I23

              end if

           end do

        end do

     end if

  end do
  !
  ! For body fitted type bodies, sum in parallel
  !
  do iimbo = 1,nrbod
     if(rbbou(iimbo) % momin(1) < -0.5_rp) then  ! momin has been initialized to -1.0 to indicate that it must be calculated if 
        ! it is not read
        call PAR_SUM(nsize,inerc(:,iimbo))
        do isize = 1,nsize
           inerc(isize,iimbo) = -inerc(isize,iimbo)
        end do
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! Compute density form mass or mass from density
  !
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod
     if( rbbou(iimbo) % densi < 0.0_rp .and. rbbou(iimbo) % massa < 0.0_rp ) then
        rbbou(iimbo) % densi = 1.0_rp
        rbbou(iimbo) % massa = rbbou(iimbo) % densi * rbbou(iimbo) % volum
     else if( rbbou(iimbo) % densi < 0.0_rp ) then
        rbbou(iimbo) % densi = rbbou(iimbo) % massa / rbbou(iimbo) % volum
     else if( rbbou(iimbo) % massa < 0.0_rp ) then
        rbbou(iimbo) % massa = rbbou(iimbo) % densi * rbbou(iimbo) % volum
     else
        !           now we allow to prescribe both mas and density even if they are not consistent with the volume. A higher mass is a 
        !           way of relaxing the movement to obtain the same stationary solution.
     end if

     if( rbbou(iimbo) % volum <= 0.0_rp ) then
        messa = intost(iimbo)
        print*,rbbou(iimbo) % densi,rbbou(iimbo) % massa,rbbou(iimbo) % volum
        call runend('ALE_INIRBO: NEGATIVE VOLUME FOR PARTICLE '//trim(messa)//'. CHECK BOUNDARY ORIENTATION')
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! Multiply inertia tensor by density
  !
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod
     if(rbbou(iimbo) % momin(1) < -0.5_rp) then  ! momin has been initialized to -1.0 to indicate that it must be calculated if 
        ! it is not read
        do iiner = 1,nsize
           rbbou(iimbo) % momin(iiner) = inerc(iiner,iimbo) * rbbou(iimbo) % densi
        end do
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! Modify COOIB so that center of gravity is at prescribed position
  ! 
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod
     do ipoib = 1,rbbou(iimbo) % npoib
        do idime = 1,ndime
           ! rbbou(iimbo) % cooin(idime,ipoib) = rbbou(iimbo) % cooib(idime,ipoib) - rbbou(iimbo) % posgr(idime)
           rbbou(iimbo) % cooib(idime,ipoib) = rbbou(iimbo) % cooin(idime,ipoib) + rbbou(iimbo) % posil(idime,1)
           !              rbbou(iimbo) % cooi2(idime,ipoib) = rbbou(iimbo) % cooib(idime,ipoib)              
        end do
     end do
  end do
  
  call memchk(two,istat,memor_dom,'INERC','ale_inirbo',inerc)
  deallocate(inerc,stat=istat)
  if(istat/=0) call memerr(two,'INERC','ale_inirbo',0_ip)
  call memchk(two,istat,memor_dom,'CALGR','ale_inirbo',kfl_calgr)
  deallocate(kfl_calgr,stat=istat)
  if(istat/=0) call memerr(two,'CALGR','ale_inirbo',0_ip)

  npoib = 0
  nboib = 0
  do iimbo = 1,nrbod
     npoib = npoib + rbbou(iimbo) % npoib
     nboib = nboib + rbbou(iimbo) % nboib
  end do
  
  !----------------------------------------------------------------------
  !
  ! Compute fixity bvess_ale - using % cooib coord
  ! 
  !----------------------------------------------------------------------
  
  call ale_fixirb()     


end subroutine ale_inirbo
