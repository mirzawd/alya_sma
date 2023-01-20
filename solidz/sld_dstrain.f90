!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_dstrain(ielem,pnode,pgaus,pmate,eldis,gpsha,gpcar,gpgi0,gpgdi,deltaE)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_dstrain
  ! NAME
  !    sld_dstrain
  ! DESCRIPTION
  !    Compute strain increment for plasticity models
  ! INPUT
  !    ELDIS ...
  !    GPGI0 ...   Previous deformation gradient tensor ............ F(n)
  ! OUTPUT
  !    GPGDI ...   Updated deformation gradient .....................F(n+1) = grad(U
  !    GPDSTRA ... Incremental strain (logaritmic if problem is nonlinear)
  ! USES
  ! USED BY
  !    sld_elmPRE
  !   end if
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  mnode,ndime
  use def_master, only       :  coupling, ITER_K, TIME_N, dtime
  use mod_sld_stress_model_comput

  implicit none
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: pmate
  real(rp),    intent(in)    :: eldis(ndime,pnode,*)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpgi0(ndime,ndime,pgaus)
  real(rp)                   :: gpig0(ndime,ndime,pgaus)

  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: deltaE(ndime,ndime,pgaus)
  real(rp)                   :: deltaF(ndime,ndime,pgaus)
  real(rp)                   :: dRrot(ndime,ndime,pgaus), dUdef(ndime,ndime,pgaus)
  real(rp)                   :: gpeva(ndime,pgaus), gpeve(ndime,ndime,pgaus)
  real(rp)                   :: gpgre(ndime,ndime,pgaus), gpgre0(ndime,ndime,pgaus)
  real(rp)                   :: bidon, gpdet(pgaus)
  integer(ip)                :: idime,jdime,kdime,mdime,igaus,inode
  integer(ip)                :: optionStrain

  !computation of strain increment as in vumat_bridge
  real(rp)                   :: gpgi1(ndime,ndime,pgaus) !F(n+1)
  real(rp)                   :: rstr0_sld(ndime,ndime,pgaus), rstr1(ndime,ndime,pgaus), udefo(ndime,ndime,pgaus)  ! F = R*U
  real(rp)                   :: temm1(ndime,ndime), temm2(ndime,ndime), temp(ndime,ndime,pgaus)
  real(rp)                   :: gpsl0_sld(ndime,ndime,pgaus), gpsl1(ndime,ndime,pgaus)
  real(rp)                   :: gpfin(ndime,ndime,pgaus)  !B=F F^t
  real(rp)                   :: gpder(ndime,ndime,pgaus)  ! DeltaR
  real(rp)                   :: gpsti(ndime,ndime,pgaus)  ! Delta(log(epsilon))
  real(rp)                   :: gpome(ndime,ndime,pgaus)
  real(rp)                   :: gplep1(ndime,ndime), gplep0(ndime,ndime)
  real(rp)                   :: gpcrt(ndime,ndime) ! F F^T
  real(rp)                   :: eval(ndime), evec(ndime,ndime)
  integer(ip)                :: nrank

  optionStrain = 2_ip

  GAUSS_POINTS_LOOP: do igaus = 1,pgaus

    if (optionStrain==1) then
    !
    ! Calculate deltaF = F_n F_{n-1}^-1
    !
    call invmtx(gpgi0(1,1,igaus),gpig0(1,1,igaus),gpdet(igaus),ndime)

    deltaF(:,:,igaus) = 0.0_rp
    do idime=1,ndime
      do jdime=1,ndime
        deltaF(idime,jdime,igaus) = deltaF(idime,jdime,igaus) + gpgdi(idime,jdime,igaus)*gpig0(jdime,idime,igaus)
      end do
    end do

    !
    ! Polar decomposition of deltaF = deltaR deltaU
    !
    call SM_polar_decomposition(1_ip,ndime,deltaF(1,1,igaus),dRrot(1,1,igaus),dUdef(1,1,igaus))

    !
    ! Calculate the strain increment = ln(delta U) = sum(log(lambda_i) e_i*e_i)
    !
    call spcdec(dUdef(1,1,igaus),gpeva(1,igaus),gpeve(1,1,igaus),bidon,1_ip,'SLD_DSTRAIN')

    do idime=1,ndime
      do jdime=1,ndime
        deltaE(idime,jdime,igaus) = 0.0_rp
      end do
    end do

    do kdime=1,ndime
      do idime=1,ndime
        do jdime=1,ndime

          if (gpeva(kdime,igaus) > 0.0_rp) then
            deltaE(idime,jdime,igaus) = deltaE(idime,jdime,igaus) +&
               log(gpeva(kdime,igaus))*gpeve(idime,kdime,igaus)*gpeve(jdime,kdime,igaus)
          end if
        end do
      end do
    end do

    else if (optionStrain==2) then
      !
      ! Compute Green-Lagrange at iter K
      !
      gpgre(:,:,igaus) = 0.0_rp
      do jdime=1,ndime
        do idime=1,ndime
          do kdime=1,ndime
            gpgre(idime,jdime,igaus) = gpgre(idime,jdime,igaus) +&
                   0.5_rp*(gpgdi(kdime,idime,igaus)*gpgdi(kdime,jdime,igaus))
          end do
        end do
      end do
      gpgre0(:,:,igaus) = 0.0_rp
      do jdime=1,ndime
        do idime=1,ndime
          do kdime=1,ndime
            gpgre0(idime,jdime,igaus) = gpgre0(idime,jdime,igaus) +&
                   0.5_rp*(gpgi0(kdime,idime,igaus)*gpgi0(kdime,jdime,igaus))
          end do
        end do
      end do
      do idime=1,ndime
         gpgre(idime,idime,igaus) = gpgre(idime,idime,igaus)-0.5_rp
         gpgre0(idime,idime,igaus) = gpgre0(idime,idime,igaus)-0.5_rp
      end do

      deltaE(:,:,igaus) = gpgre(:,:,igaus) - gpgre0(:,:,igaus)

    else if (optionStrain==3) then
      do inode=1,pnode
        deltaE(1,1,igaus) = gpcar(1,inode,igaus)*(eldis(1,inode,ITER_K)-eldis(1,inode,TIME_N))
        deltaE(2,2,igaus) = gpcar(2,inode,igaus)*(eldis(2,inode,ITER_K)-eldis(2,inode,TIME_N))
        deltaE(3,3,igaus) = gpcar(3,inode,igaus)*(eldis(3,inode,ITER_K)-eldis(3,inode,TIME_N))

        deltaE(2,3,igaus) = gpcar(2,inode,igaus)*(eldis(1,inode,ITER_K)-eldis(1,inode,TIME_N)) +&
                             gpcar(1,inode,igaus)*(eldis(2,inode,ITER_K)-eldis(2,inode,TIME_N))
        deltaE(1,3,igaus) = gpcar(3,inode,igaus)*(eldis(2,inode,ITER_K)-eldis(2,inode,TIME_N)) +&
                             gpcar(2,inode,igaus)*(eldis(3,inode,ITER_K)-eldis(3,inode,TIME_N))
        deltaE(1,2,igaus) = gpcar(3,inode,igaus)*(eldis(1,inode,ITER_K)-eldis(1,inode,TIME_N)) +&
                             gpcar(1,inode,igaus)*(eldis(3,inode,ITER_K)-eldis(3,inode,TIME_N))
      end do

    else if (optionStrain==4) then   !como en vumat_bridge

      gpgi1 = gpgdi

      !Get R and U (F=R*U : gpgi1 = rstr1*udefo)
      call SM_polar_decomposition(1_ip,ndime,gpgi1(1,1,igaus),rstr1(1,1,igaus),udefo(1,1,igaus))

      !Calculate b = F F^T
      do kdime=1,ndime
        do jdime=1,ndime
          gpfin(jdime,kdime,igaus) = 0.0_rp
          do idime=1,ndime
            gpfin(jdime,kdime,igaus) = gpfin(jdime,kdime,igaus) + gpgi1(jdime,idime,igaus)*gpgi1(kdime,idime,igaus)
          end do
        end do
      end do
      !Find eigenvalues and eigenvectors of b
      call spcdec(gpfin(1,1,igaus),gpeva(1,igaus),gpeve(1,1,igaus),bidon,1_ip,'SLD_DSTRAIN')

      !Calculate updated epsilon_log
      do kdime=1,ndime
        do idime=1,ndime
          do jdime=1,ndime
            gpsl1(idime,jdime,igaus) = gpsl1(idime,jdime,igaus) +&
                                     log(sqrt(gpeva(kdime,igaus)))*gpeve(idime,kdime,igaus)*gpeve(jdime,kdime,igaus)
          end do
        end do
      end do

      !Calculate the angular velocity omega
          do idime=1,ndime
             do jdime=1,ndime
                  temp(idime,jdime,igaus)= (rstr1(idime,jdime,igaus)-rstr0_sld(idime,jdime,igaus))/dtime
             end do
          end do
           !temp*R^T
          do idime=1,ndime
             do jdime=1,ndime
                gpome(idime,jdime,igaus)=0.0_rp
                do kdime=1,ndime
                  gpome(idime,jdime,igaus)= gpome(idime,jdime,igaus)+&
                                            temp(idime,kdime,igaus)*rstr1(jdime,kdime,igaus)
                end do
             end do
          end do

         !Archive R
         do idime=1,ndime
             do jdime=1,ndime
               rstr0_sld(idime,jdime,igaus)=rstr1(idime,jdime,igaus)
             end do
         end do

          !Calculate deltaR
          do idime=1,ndime
             do jdime=1,ndime
                  temm1(idime,jdime)= -0.5_rp*dtime*gpome(idime,jdime,igaus)
             end do
             temm1(idime,idime)= 1.0_rp + temm1(idime,idime)
          end do

          call invmtx(temm1(1,1),temm2(1,1),bidon,ndime)

               do idime=1,ndime
             do jdime=1,ndime
                  temm1(idime,jdime)= 0.5_rp*dtime*gpome(idime,jdime,igaus)
             end do
             temm1(idime,idime)= 1.0_rp + temm1(idime,idime)
          end do

           do idime=1,ndime
             do jdime=1,ndime
                gpder(idime,jdime,igaus)=0.0_rp
                do kdime=1,ndime
                  gpder(idime,jdime,igaus)=gpder(idime,jdime,igaus)+&
                                           temm2(idime,kdime)*temm1(kdime,jdime)
                end do
             end do
           end do

          !Calculate StrainInc (gpsti in Alya and )

          do idime=1,ndime
            do jdime=1,ndime
               temm1(idime,jdime)=0.0_rp
               do kdime=1,ndime
                  do mdime=1,ndime
                     temm1(idime,jdime)=temm1(idime,jdime) &
                          + gpder(idime,kdime,igaus) * gpder(jdime,mdime,igaus)&
                          * gpsl0_sld(kdime,mdime,igaus)
                  end do
                end do
            end do
          end do

          do idime=1,ndime
            do jdime=1,ndime
                  gpsti(idime,jdime,igaus) = 0.0_rp
            end do
          end do

        do idime=1,ndime
            do jdime=1,ndime
                  gpsti(idime,jdime,igaus) = gpsl1(idime,jdime,igaus)-temm1(idime,jdime)
            end do
          end do

        ! Archive strain_L
        do idime=1,ndime
          do jdime=1,ndime
            gpsl0_sld(idime,jdime,igaus) = gpsl1(idime,jdime,igaus)
          end do
        end do

        deltaE(:,:,igaus) = gpsti(:,:,igaus)

    else if (optionStrain == 5_ip) then

        gplep1 = 0.0_rp
        gpcrt = matmul(gpgdi(:,:,igaus),transpose(gpgdi(:,:,igaus)))

        call spcdec(gpcrt,eval,evec,nrank,1_ip,'SLD_DSTRAIN,NDIME=3,COMPUTING LOG STRAIN TENSOR')

        do idime=1,ndime
          do jdime=1,ndime
            do kdime=1,ndime
              gplep1(idime,jdime) = gplep1(idime,jdime) + log(sqrt(eval(kdime)))*evec(idime,kdime)*evec(jdime,kdime)
            end do
          end do
        end do

        gplep0 = 0.0_rp
        gpcrt = matmul(gpgi0(:,:,igaus),transpose(gpgi0(:,:,igaus)))

        call spcdec(gpcrt,eval,evec,nrank,1_ip,'SLD_DSTRAIN,NDIME=3,COMPUTING LOG STRAIN TENSOR')

        do idime=1,ndime
          do jdime=1,ndime
            do kdime=1,ndime
              gplep0(idime,jdime) = gplep0(idime,jdime) + log(sqrt(eval(kdime)))*evec(idime,kdime)*evec(jdime,kdime)
            end do
          end do
        end do

        deltaE(:,:,igaus) = gplep1 - gplep0

    !
    ! logE = sum_i(0.5_rp*log(lambda_i) N_i X N_i
    !
    else if (optionStrain == 6_ip) then

      gplep1 = 0.0_rp
      !Get R and U (F=R*U : gpgi1 = rstr1*udefo)
      gpgi1 = gpgdi
      call SM_polar_decomposition(1_ip,ndime,gpgi1(1,1,igaus),rstr1(1,1,igaus),udefo(1,1,igaus))

      ! Vap's and vep's of udefo
      call spcdec(udefo(:,:,igaus),eval,evec,nrank,1_ip,'SLD_DSTRAIN,NDIME=3,COMPUTING LOG STRAIN TENSOR')

      do idime=1,ndime
        do jdime=1,ndime
          do kdime=1,ndime
            gplep1(idime,jdime) = gplep1(idime,jdime) + log(eval(kdime))*evec(idime,kdime)*evec(jdime,kdime)
          end do
        end do
      end do

      gplep0 = 0.0_rp
      !Get R and U (F=R*U : gpgi1 = rstr1*udefo)
      call SM_polar_decomposition(1_ip,ndime,gpgi0(1,1,igaus),rstr1(1,1,igaus),udefo(1,1,igaus))

      ! Vap's and vep's of udefo
      call spcdec(udefo(:,:,igaus),eval,evec,nrank,1_ip,'SLD_DSTRAIN,NDIME=3,COMPUTING LOG STRAIN TENSOR')

      do idime=1,ndime
        do jdime=1,ndime
          do kdime=1,ndime
            gplep0(idime,jdime) = gplep0(idime,jdime) + log(eval(kdime))*evec(idime,kdime)*evec(jdime,kdime)
          end do
        end do
      end do

      deltaE(:,:,igaus) = gplep1 - gplep0

    end if

  end do GAUSS_POINTS_LOOP

end subroutine
