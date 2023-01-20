!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Exmedi
!> @{
!> @name    Fractional difussion for electrophisiology
!> @file    mod_exm_fractional_diffusion.f90
!> @author  Alfonso Santiago
!> @brief   Fractional diffusion 
!> @details Fractional diffusion for electrophysiology. Based on
!>          Nicole Cusimano's mathematical developement and 
!>          Mathlab code.
!> @{
!
!-----------------------------------------------------------------------
module mod_exm_fractional_diffusion

    use def_master,         only : dtime, momod, modul
    use def_master,         only : amatr
    use def_exmedi,         only : amatr_auxi_exm
    use def_kintyp,         only : ip, rp, lg
    use def_kintyp_solvers, only : soltyp
    use def_parame,         only : pi
    use mod_messages,       only : livinf
    use def_master,         only : intost, retost


    implicit none


    character(300)          :: messa

    real(rp),  parameter :: epsil = epsilon(1.0_rp)

    integer(ip), save                             :: nintp ! Number of INTegration Points
    real(rp), save                                :: coeff

    real(rp), save                                :: MSeigenmax, &
                                                     MSeigenmin
    real(rp), save                                :: ellipK, &
                                                     ellipKp
    real(rp), save                                :: k

   complex(rp), save, dimension(:), allocatable  :: quadr_nodes(:)
   complex(rp), save, dimension(:), allocatable  :: quadr_weigh(:)
   complex(rp), save, dimension(:), allocatable  :: dxidt(:)
 
    interface jacobi_elliptic_function
        module procedure  jacobi_elliptic_scalar, &
                          jacobi_elliptic_matrix
    end interface


   private :: jacobi_elliptic_function
   private :: jacobi_elliptic_scalar
   private :: jacobi_elliptic_matrix

   public :: frac_diff_initialisation
   public :: solve_fractional_diffusion_time_step

   contains
    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-FEB-18
    !> @brief   Initialisation for the fractional diffusion
    !> @details This subroutine initialises the fractional diffusion
    !>          for electrophysiology. The objective is to obtain
    !>                 MSeigenmax  : the largest eigenvalue of M{^1}S
    !>                 MSeigenmin  : the smalles eigenvalue of M{-1}S
    !>                 quadr_nodes : the quadrature nodes
    !>                 quadr_weigh : the quadrature weights
    !>         
    !>         
    !-----------------------------------------------------------------------
    subroutine frac_diff_initialisation()

        use def_exmedi, only: fract_diff_coef_exm 
        use def_exmedi, only: fract_diff_nintp_exm
        use def_domain, only: nmate

        implicit none
        complex(rp), dimension(:), allocatable  :: midpn
        complex(rp), dimension(:), allocatable  :: sn, cn, dn
        real(rp)                                :: errmax, errmin
        integer(ip)                             :: i


        if (nmate.ge.2_ip) then
          call runend ('FRACTIONAL DIFFUSION: FRACTIONAL DIFFUSION NOT PREPARED FOR MORE THAN ONE MATERIAL (but youre invited to code it)')
        else
          nintp=fract_diff_nintp_exm(1)
          coeff= fract_diff_coef_exm(1)
        endif

        allocate(quadr_nodes(nintp))
        allocate(quadr_weigh(nintp))
        allocate(midpn(nintp))
        allocate(sn(nintp))
        allocate(cn(nintp))
        allocate(dn(nintp))
        allocate(dxidt(nintp))


        MSeigenmin=0.0_rp
        Mseigenmax=0.0_rp
        k=0.0_rp

        messa = '     INITIATING FRACTIONAL DIFFUSION... ' 
        call livinf(0_ip,messa,1_ip)
        messa = '       COMPUTING EIGENVALUES BY POWER METHOD ' 
        call livinf(0_ip,messa,1_ip)

        call obtain_eigenvalues(MSeigenmax, MSeigenmin, errmax, errmin)

        messa = '           MAX EIGENVALUE: '//trim(retost(MSeigenmax))//' ERROR: '//trim(retost(errmax))
        call livinf(0_ip,messa,1_ip)
        messa = '           MIN EIGENVALUE: '//trim(retost(MSeigenmin))//' ERROR: '//trim(retost(errmin))
        call livinf(0_ip,messa,1_ip)


        k=(sqrt(MSeigenmax/MSeigenmin)-1.0_rp)/(sqrt(MSeigenmax/MSeigenmin)+1.0_rp)
 

        call elliptic_integral(-log(k)/pi,ellipK,ellipKp)
        ! Elliptic integral is OK


        do i=1,nintp
          midpn(i)=nintp-0.5_rp-i
          midpn(i) = 0.5_rp*cmplx(0.0,1.0,kind=rp)*ellipKp - ellipK + midpn(i)*2.0_rp*K/nintp
        enddo
        
        do  i=1,nintp
           call jacobi_elliptic_function(midpn(i), -log(real(k,kind=rp))/pi, sn(i), cn(i), dn(i)) 
           ! Elliptic function is OK
           quadr_nodes(i) = sqrt(MSeigenmax*MSeigenmin)*(1/k+sn(i))*(1/k-sn(i))         
           dxidt(i)=cn(i)*dn(i)/(1/k-sn(i))**2.0_rp
        enddo

        write(6,*) 'quadr_nodes:', quadr_nodes
        write(6,*) 'dxidt:', dxidt
        write(6,*) '------------------'
        write(6,*) 'END INITIALISATION'
        write(6,*) '------------------'
        write(6,*) '------------------'


    end subroutine frac_diff_initialisation
    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-JULY-17
    !> @brief   
    !> @details 
    !>         
    !>         
    !>         
    !-----------------------------------------------------------------------
    subroutine obtain_eigenvalues(lambda_max_mine, lambda_min_mine, errormax, errormin)

        use mod_memory, only :  memory_alloca, memory_deallo
        use def_exmedi, only: kfl_timet_exm, ndofn_exm
        use def_domain, only: npoin, vmass
        use def_domain, only: r_dom, c_dom ! CSR matrix indexes
        use def_solver, only: solve_sol, memit, SOL_NO_PRECOND
        use mod_solver, only : solver_condition_number
        use mod_solver, only : solver_parallel_scalar_product
        use mod_solver, only : solver_parallel_SpMV
        use mod_solver, only : solver_parallel_double_scalar_product
        use mod_solver, only : solver_parallel_vector_L2norm
        use mod_matrix , only: matrix_scaling_CSR

        implicit none

        real(rp), intent(out)                            :: lambda_min_mine, lambda_max_mine
        real(rp), intent(out)                            :: errormax, errormin
        integer(ip)                                      :: iiter
        integer(ip)                                      :: ii
        integer(ip)                                      :: i, j
        type(soltyp), pointer                            :: solve_local(:)
        integer(ip)                                      :: ndof,nn, nsize, ncols, nrows
        real(rp)                                         :: ynorm, xnorm, xdoty
        real(rp)                                         :: lambda, lambda_new
        integer(ip),  parameter                          :: maxit = 10000
        real(rp), parameter                              :: zeror =  epsilon(1.0_rp)
        real(rp)                                         :: toler, eps
        real(rp),     pointer                            :: xx(:), xx_maxeig(:), xx_maxeig2(:)
        real(rp),     pointer                            :: yy(:)
        real(rp),     pointer                            :: auxi(:)
        real(rp)                                         :: auxsca1, auxsca2

        ! Unity test variables
        logical, parameter                               :: unity_test=.false.
        real(rp), allocatable, target, dimension(:)      :: aa_unity, diag_unity
        integer(ip), allocatable, target, dimension(:)   :: ia_unity, ja_unity
        real(rp), allocatable, dimension(:,:) :: A(:,:)                             

       if (kfl_timet_exm .ne. 2_ip) call runend ("MOD_EXM_FRACTIONAL_DIFFUSION: EIGENVALUES CAN BE ONLY COMPUTED WITH IMPLICIT SCHEME")
        ! Local solver metadata

        allocate(solve_local(1))
        solve_local(1) = solve_sol(1)
        solve_local(1) % kfl_preco = SOL_NO_PRECOND
        nn    = solve_local(1) % nequa
        ndof  = solve_local(1) % ndofn
        nrows = solve_local(1) % nequa * ndof
        ncols = solve_local(1) % ncols * ndof
        nsize = max(1_ip,ncols,nrows)

        nullify(xx)
        nullify(xx_maxeig)
        nullify(xx_maxeig2)
        nullify(yy)
        nullify(auxi)
        call memory_alloca(memit,'XX','fract_diff_eigenv',xx,nsize)
        call memory_alloca(memit,'XX_MAXEIG','frac_diff_eigenv',xx_maxeig,nsize)
        call memory_alloca(memit,'XX_MAXEIG2','frac_dif_eigenv',xx_maxeig2,nsize)
        call memory_alloca(memit,'YY','frac_diff_eigenv',yy,nsize)
        call memory_alloca(memit,'auxi','frac_diff_eigenv',auxi,nsize)
    
        !
        ! how to copy matrix an(ii) in CSR format in the 
        ! dense matrix A.
        !
        !   A=0.0_rp
        !   do i=1,7
        !    do j=solve%ia(i),solve%ia(i+1)-1
        !        A(i,solve%ja(j))=an(j)
        !    enddo
        !   enddo

        ! unity test
        !
        !                   +-          -+
        !                   |  1   2   3 |
        ! Original matrix = |  4   5   6 |
        !                   |  7   8   9 |
        !                   +-          -+
        !
        !                   +-          -+
        !                   |  2   0   0 |
        ! Mass Scaling =    |  0   4   0 |
        !                   |  0   0   8 |
        !                   +-          -+
        !
        !                   +-          -+
        !                   |  2   4   6 |
        ! Left scaling =    | 16  20  24 |
        !                   | 56  64  72 |
        !                   +-          -+
        !
        !                  
        !                  
        ! Eigenvalues = 0, 97.2095594029714, -3.21616018663792
        !                 
        !                 

        if(unity_test) then
            write(6,*) 'WARNING, WARNING. INITIATING UNITY TEST'
            npoin=7_ip
            ndofn_exm=1_ip

            allocate(A(npoin,npoin))
            allocate(ia_unity(npoin+1),ja_unity(npoin*npoin))
            allocate(aa_unity(npoin*npoin))
            allocate(diag_unity(npoin))

            aa_unity   =(/10.842555555555553_rp, 1.6263888888888887_rp, 1.6263888888888887_rp,&
            -0.54213888888888884_rp, -5.9634166666666655_rp, -5.9634166666666655_rp,&
            -1.6263611111111109_rp, 0.40659722222222217_rp, 0.67765972222222204_rp,&
            -0.81319444444444433_rp, -0.27106249999999998_rp, 0.40659722222222217_rp,&
            0.67765972222222204_rp, -0.81319444444444433_rp, -0.27106249999999998_rp,&
            -0.13553472222222221_rp, 0.67765972222222204_rp, -0.27106249999999998_rp,&
            -0.27106249999999998_rp, -4.4725624999999996_rp, -2.4395833333333332_rp,&
            -0.81318749999999995_rp, 8.5385208333333331_rp, 0.40660416666666666_rp,&
            -1.2197916666666666_rp, -4.4725624999999996_rp, -2.4395833333333332_rp,&
            -0.81318749999999995_rp, 0.40660416666666666_rp, 8.5385208333333331_rp,&
            -1.2197916666666666_rp, -1.2197708333333332_rp, -0.81318749999999995_rp,&
            -0.81318749999999995_rp, -1.2197916666666666_rp, -1.2197916666666666_rp,&
            5.2857291666666661_rp/)

            ia_unity = (/1, 8, 12, 16, 20, 26, 32, 38/)
            ja_unity = (/1, 2, 3, 4, 5, 6, 7,1, 2, 5, 7, 1, 3, 6, 7, &
                        1, 4, 5, 6, 1, 2, 4, 5, 6, 7, 1, 3, 4, 5, 6,&
                        7, 1, 2, 3, 5, 6, 7/)

            diag_unity = (/ 1.0_rp,1.0_rp,1.0_rp,1.0_rp,1.0_rp,1.0_rp,1.0_rp/)

            amatr_auxi_exm=>aa_unity
            amatr=>aa_unity
            r_dom=>ia_unity
            c_dom=>ja_unity
            vmass=>diag_unity

            solve_local(1) % nequa = 7_ip
            solve_local(1) % ncols = 7_ip
            solve_local(1) % nequa_own = 7_ip 
            solve_local(1)%ia=> ia_unity
            solve_local(1)%ja=> ja_unity

            write(6,*) 'MATRIX UNITY TEST:'
            A=0.0_rp
            do i=1,npoin
             do j=solve_local(1)%ia(i),solve_local(1)%ia(i+1)-1
                 A(i,solve_local(1)%ja(j))=amatr(j)
             enddo
            enddo

            do i=1,7
                do j=1,7
                   write(6, fmt="(f8.5,a1)", advance="no" ) A(i,j), " "
                enddo
                write(6,*) ' ' 
            enddo

            write(6,*) '-----------------------------------'

        endif
        !
        ! End input unity test
        !

        !! Obtain M^-1 S
        !
        do i=1,size(amatr, KIND=ip)
          amatr(i) = amatr_auxi_exm(i)
        enddo

   
        call matrix_scaling_CSR(1_ip,npoin,ndofn_exm,r_dom,c_dom,amatr,vmass,LEFT_SCALING=.true.,RIGHT_SCALING=.false.)
    



        
        !------------------------------------------------------------------------
        !       COMPUTE EIGENVALUES
        !------------------------------------------------------------------------
        !
        ! Lambda_max
        !
        iiter  = 0
        toler=1.0e-09_rp
        eps    = huge(1.0_rp)
        lambda = -1.0_rp

        xx = 1.0_rp
        call solver_parallel_vector_L2norm(solve_local(1),xx,xnorm,OPENMP=.true.)
        do ii = 1,ncols
           xx(ii) = xx(ii) / xnorm
        end do

        do while( iiter < maxit .and. eps >= toler )
           iiter = iiter + 1
           
           !A*x
           call solver_parallel_SpMV(solve_local(1),amatr,xx,yy,OPENMP=.true.)


           call solver_parallel_double_scalar_product(solve_local(1),yy,yy,xx,ynorm,xdoty,OPENMP=.true.)
           ynorm = max(sqrt(ynorm),zeror)
           
           do ii = 1,ncols
              xx(ii) = yy(ii) / ynorm
              xx_maxeig(ii)=xx(ii)
           end do
           lambda_new = xdoty 
           eps        = abs(lambda_new-lambda) / max(abs(lambda_new),zeror)
           lambda     = lambda_new

           !if(unity_test) write(6,*)'max=', iiter, eps,lambda

        end do


        ! 
        ! confirm A*xx_maxeig = lambda*xx_maxeig
        ! 
        call solver_parallel_SpMV(solve_local(1),amatr,xx_maxeig,xx,OPENMP=.true.)

        auxi=1.0_rp
        call solver_parallel_scalar_product(solve_local(1),xx,auxi,auxsca1)

        call solver_parallel_scalar_product(solve_local(1),lambda*xx_maxeig,auxi,auxsca2)
            
        errormax=abs(auxsca1-auxsca2)

        if(unity_test) write(6,*) 'maximum eigenvalue:', lambda, 'ERROR:', errormax

        if( iiter <= maxit ) then
          lambda_max_mine = lambda
        else
          lambda_max_mine=-888.0_rp ! Failed to converge.
        endif


        if(unity_test) write(6,*) '--------------------'
        !
        ! Lambda_min
        !
        iiter  = 0
        toler=1.0e-05_rp
        eps    = huge(1.0_rp)
        lambda = -100.0_rp
        xx=lambda_max_mine

        !call solver_parallel_vector_L2norm(solve_local(1),xx,xnorm,OPENMP=.true.)
        !do ii = 1,ncols
        !   xx(ii) = xx(ii) / xnorm
        !end do

        do while( iiter < maxit .and. eps >= toler )
           ! do not compute if maximum eigenvalue not converged 
           if(lambda_max_mine.eq.-888.0_rp) then
             lambda=-888.0_rp
             exit
           endif

           iiter = iiter + 1

           !A*x
           call solver_parallel_SpMV(solve_local(1),amatr,xx,yy,OPENMP=.true.)

           !A*x-lambda_max*I
           do ii = 1,ncols
              yy(ii) = yy(ii) - lambda_max_mine * xx(ii)
           end do

           call solver_parallel_double_scalar_product(solve_local(1),yy,yy,xx,ynorm,xdoty,OPENMP=.true.)
           ynorm = max(sqrt(ynorm),zeror)
           do ii = 1,ncols
              xx(ii) = yy(ii) / ynorm
              xx_maxeig2(ii)=xx(ii)
           end do
           lambda_new = xdoty + lambda_max_mine
           eps        = abs(lambda_new-lambda) / max(abs(lambda_new),zeror)
           lambda     = lambda_new

           if(unity_test) write(6,*)'min1=',iiter, eps,lambda
        end do

        ! 
        ! confirm A*xx_maxeig2 = 0
        ! 
        call solver_parallel_SpMV(solve_local(1),amatr,xx_maxeig2,xx,OPENMP=.true.)
        auxi=1.0_rp
        call solver_parallel_scalar_product(solve_local(1),xx,auxi,errormin)

        !write(6,*) 'Zero eigenvector found with error: ', errormin
        !
        ! NOTE TODO NOTE TODO NOTE TODO
        !
        ! Another good way to obtain xx_maxeig2 is to solve
        ! A*xx_maxeig2=0
        !

        if(unity_test) write(6,*) '--------------------'
        !
        ! Smallest non-zero eigenvalue
        !
        if( (eps.eq.0.0_rp) .or. (iiter .ge. maxit) .or. (eps .gt. toler) .or. (lambda .lt. 1e-013_rp) ) then

          iiter  = 0
          eps    = huge(1.0_rp)
          toler=1.0e-05_rp
          lambda = -100.0_rp
          xx=lambda_max_mine


          do while( iiter < maxit .and. eps >= toler )
             ! do not compute if maximum eigenvalue not converged 
             if(lambda_max_mine.eq.-888.0_rp) then
               lambda=-888.0_rp
               exit
             endif
             iiter = iiter + 1

             !A*x
             call solver_parallel_SpMV(solve_local(1),amatr,xx,yy,OPENMP=.true.)


             call solver_parallel_vector_L2norm(solve_local(1),yy,xnorm,OPENMP=.true.)
             !write(6,*) 'yynorm:', xnorm
            

             !w_min*x
             call solver_parallel_scalar_product(solve_local(1),xx_maxeig2,xx,xdoty)

             !write(6,*) 'xdoty: ', xdoty

             ! Eigenvalue shift + Hotelling deflation
             do ii = 1,ncols
               yy(ii) = yy(ii) - lambda_max_mine * xx(ii)
               yy(ii) = yy(ii) + lambda_max_mine * xx_maxeig2(ii)*xdoty
             end do

             !write(6,*) '|after deflation: '
             !write(6,*) 'yy', yy(1),yy(2),yy(3)

             call solver_parallel_vector_L2norm(solve_local(1),yy,xnorm,OPENMP=.true.)
             !write(6,*) 'yynorm:', xnorm

             call solver_parallel_double_scalar_product(solve_local(1),yy,yy,xx,ynorm,xdoty,OPENMP=.true.)
             ynorm = max(sqrt(ynorm),zeror)
             do ii = 1,ncols
                xx(ii) = yy(ii) / ynorm
                xx_maxeig(ii)=xx(ii)
             end do
             lambda_new = xdoty + lambda_max_mine
             eps        = abs(lambda_new-lambda) / max(abs(lambda_new),zeror)
             lambda     = lambda_new

             !if(unity_test) write(6,*)'min2=',iiter, eps,lambda

             !if(iiter.eq.2_ip) exit
          end do



        endif

        ! 
        ! confirm A*xx_maxeig2 = lambda*xx_maxeig2
        ! 
        call solver_parallel_SpMV(solve_local(1),amatr,xx_maxeig,xx,OPENMP=.true.)

        auxi=1.0_rp
        call solver_parallel_scalar_product(solve_local(1),xx,auxi,auxsca1)

        call solver_parallel_scalar_product(solve_local(1),lambda*xx_maxeig,auxi,auxsca2)
            
        errormin=abs(auxsca1-auxsca2)


        if( iiter <= maxit ) then
          lambda_min_mine = lambda
        else
          lambda_min_mine=-888.0_rp ! Failed to converge.
        endif

        if(unity_test) write(6,*) 'minimum eigenvalue: ', lambda_min_mine, 'ERROR:', eps


        call memory_deallo(memit,'XX','fract_diff_eigenv',xx)
        call memory_deallo(memit,'XX_MAXEIG','frac_diff_eigenv',xx_maxeig)
        call memory_deallo(memit,'XX_MAXEIG2','frac_dif_eigenv',xx_maxeig2)
        call memory_deallo(memit,'YY','frac_diff_eigenv',yy)
        call memory_deallo(memit,'auxi','frac_diff_eigenv',auxi)
        !
        !
        !------------------------------------------------------------------------
        !    END COMPUTE EIGENVALUES
        !------------------------------------------------------------------------





        !
        ! OLD CALL
        !

        !write(6,*) '--------------------------------------------------------'
        !write(6,*) '   USING THE OLD SUBROUTINE      '
        !write(6,*) '--------------------------------------------------------'
        !call solver_condition_number(solve_local(1), amatr, condition_number, lambda_min_mine, lambda_max_mine)
        !
        ! END OLD CALL
        !

        if( (lambda_min_mine .eq. -888.0_rp) .or. (lambda_max_mine .eq. -888.0_rp) ) call RUNEND("MOD_EXM_FRACTIONAL_DIFFUSION::initialisation. Eigenvalues power method not converged")

        !
        !Unity test verification
        !
        if(unity_test) then
          if(abs(lambda_max_mine-17.454033265288_rp).gt.1e-4_rp .or. abs(lambda_min_mine-0.565070727370468_rp).gt.1e-4_rp)then
            write(6,*) 'eig_max: ', lambda_max_mine, '| eig_min', lambda_min_mine
            call runend('EXM_MOD_FRACTIONAL_DIFFFUSION: WRONG EIGENVALUES ON UNITY TEST')
          endif
          call runend('MOD_EXM_FRACTIONAL_DIFFUSION: Unity test executed')
        endif
        !
        !End verification unity test
        !


    end subroutine obtain_eigenvalues
    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-JAN-30
    !> @brief   Complete eliptic integral of the first kind, with complement
    !> @details Returns the value of the complete elliptic integral of the 
    !>          first kind, evaluated at M=exp(-2*pi*L), 0<L<inf
    !>          Also returns the result for the complementary parameter 1-M
    !>          which is useful when M<EPS. Even when M<1e-6.
    !>          Recall that the elliptic modulus k is related to the parameter
    !>          M by M=k**2
    !>         
    !>         This subrutine uses the method of the arithmetic-geometric 
    !>         mean described in 17.6 of M. Abramowitz and IA Stegun
    !>         "handbook of mathematical functions", Dover ,1965.
    !>         
    !>         This subroutine was adapted  from Toby Driscoll's
    !>         ellipkkp.m.        
    !>         
    !>         
    !-----------------------------------------------------------------------
    subroutine elliptic_integral(L,K,Kp)

    implicit none
    real(rp), intent(in)       :: L !Must be a scalar
    real(rp), intent(out)      :: K, Kp
    real(rp)                   :: m

      
     K=0.0_rp
     Kp=0.0_rp

    if(L .gt. 10.0_rp) then
    !If m=exp(e,-2*pi*L) is extremely small, use 0(m) approximations
      K=pi/2.0_rp
      Kp=pi*L+log(4.0_rp)

    else
    ! If  not, do the proper computations

      m=exp(-2.0_rp*pi*L)
      call computeK(m,K,'GETKK')

      call computeK(m,Kp,'GETKP')

    endif
    

    contains
     
        subroutine computeK(m_in,K_out,option)
            implicit none
            real(rp), intent(in)  :: m_in
            real(rp), intent(out) :: K_out
            character(len=5)      :: option
            real(rp)              :: a0, b0,s0, a1, b1, c1 , mm, w1, m_comp
            real(rp),parameter    :: eps=1e-10
            integer(ip)           ::i1

            m_comp = 0_rp
            b0 = 0_rp

            !! Initialisation depending on the case
            if(option .eq. 'GETKK') then
              a0=1.0_rp
              b0=sqrt(1.0_rp-m_in)
              s0=m_in
            elseif(option .eq. 'GETKP') then
              a0=1.0_rp
              b0=sqrt(m_in)
              s0=1.0_rp-m_in
            else
              call runend ('MOD_EXMEDI: SUBRU ELLIPTIC: OPTION NO DEFINED')
            endif

            i1=0_ip
            mm=1.0_rp
            K_out=0.0_rp

            !! Loop independent of the case
            !
            do while(mm .gt. eps)
              a1=(a0+b0)/2.0_rp
              b1=sqrt(a0*b0)
              c1=(a0-b0)/2.0_rp
              i1=i1+1_ip
              w1=(2.0_rp**i1)*(c1**2.0_rp)

              !! Code adaptation as w1 is not a matrix/vector
              !
              !mm=maxval(maxval(w1))
              mm=abs(w1)
              !
              !!
              s0=s0+w1
              a0=a1
              b0=b1
            end do
            
            K_out = pi/(2*a1)
            
            !! Ending, dependent on the case
            !
            if(option .eq. 'GETKK') then
              m_comp=1.0_rp
            elseif(option .eq. 'GETKP') then
              m_comp=0.0_rp
            else
              call runend ('MOD_EXMEDI: SUBRU ELLIPTIC: OPTION NO DEFINED')
            endif
            
            !! Code adaptation as m_in is not a vector
            !! And we dont have find and isempty functions
            !
            !im = find(m_in==m_comp)
            !if( the vector IM is not empty ) then
              !K_out(im) = K_out(im)*huge(K_out(im))
            !endif
            if(m_in .eq. m_comp) then
              K_out = K_out*huge(K_out)
            endif
            !
            !!

        end subroutine computeK

    end subroutine elliptic_integral

    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-JAN-30
    !> @brief   Jacobi Elliptic functions for complex arguments
    !> @details Returns the values of the Jacobi elliptic functions evaluated
    !>          at complex argument U, and parameter M=exp(-2*pi*L), 0<L<inf.
    !>          Recall that M=k**2 where k is thee elliptic modulus.
    !>          U may be a matrix, L must be a scalar. The entries of U
    !>          are expected to lie within the rectangle |Re(u)|<K,
    !>          0<Im(U)<Kp, where K and Kp were computed from the
    !>          elliptic integral.
    !>         
    !>         This algorithm is the descending Landen transformation
    !>         described in L. Howell's PhD thesis from MIT. Additional
    !>         formulas from Gradshteyn and Ryzhik, 5th ed., and  
    !>         Abbramowitz and Stegun
    !>         
    !>         This subroutine was adapted  from Toby Driscoll's
    !>         ellipjc.m.        
    !>         
    !>         
    !-----------------------------------------------------------------------

    recursive subroutine jacobi_elliptic_matrix(U_in, L , sn, cn, dn, flag)
      implicit none
      complex(rp), intent(in)            :: U_in(:,:) ! U may be a matrix (of complex values)
      real(rp), intent(in)               :: L      ! Must be a scalar
      logical(lg), optional, intent(in)  :: flag
     
      complex(rp), intent(out)              :: sn(:,:), cn(:,:), dn(:,:) !No allocatable needed

      complex(rp), allocatable           :: v(:,:)
      complex(rp), allocatable           :: u(:,:)
      complex(rp), allocatable           :: sn1(:,:),cn1(:,:),dn1(:,:), denom(:,:)

      real(rp),parameter                 :: eps=1e-6_rp
      real(rp)                           :: m
      real(rp)                           :: K,Kp
      logical(lg), allocatable           :: high(:,:)
      integer(ip)                        :: ud1,ud2,i,j
      real(rp)                           :: kappa, x, mu


     call runend ('MOD_EXM_FRACTIONAL_DIFFUSION_JACOBI_ELLIPTIC_MATRIX. THIS MATRIX FUNCTION HASN BEEN TESTED')

      sn(:,:)=cmplx(0.0,0.0,kind=rp)
      cn(:,:)=cmplx(0.0,0.0,kind=rp)
      dn(:,:)=cmplx(0.0,0.0,kind=rp)

      ud1=size(U,1,KIND=ip)
      ud2=size(U,2,KIND=ip)

      allocate(high (ud1,ud2))
      allocate(sn1 (ud1,ud2))
      allocate(cn1 (ud1,ud2))
      allocate(dn1 (ud1,ud2))
      allocate(v (ud1,ud2))
      allocate(u (ud1,ud2))
      allocate(denom (ud1,ud2))

      u=U_in
      high=.false.
      sn1(:,:)=cmplx(0.0,0.0,kind=rp)
      cn1(:,:)=cmplx(0.0,0.0,kind=rp)
      dn1(:,:)=cmplx(0.0,0.0,kind=rp)
      v(:,:)=cmplx(0.0,0.0,kind=rp)
      denom(:,:)=cmplx(0.0,0.0,kind=rp)

      if(present(flag))then
          high=.false.
          m=L
      else
          call elliptic_integral(L,K,Kp)
             
          do i = 1, ud1
             do j = 1, ud2
               if(aimag(u(i,j)).gt.Kp/2.0_rp) then
                 high(i,j) = .true.
                 u(i,j)=cmplx(0.0,1.0,kind=rp)*Kp-u(i,j)
               endif
             enddo
          enddo
          m=exp(-2.0_rp*pi*L)

      endif


      if (m.lt.4.0_rp*eps) then
          do i = 1, ud1
             do j = 1, ud2
               sn=sin(u(i,j))+m/4.0_rp*(sin(u(i,j))*cos(u(i,j))-u)*cos(u(i,j))
               cn=cos(u(i,j))+m/4.0_rp*(-sin(u(i,j))*cos(u(i,j))+u)*sin(u(i,j))
               dn=1.0_rp+m/4*(cos(u(i,j))**2.0_rp-sin(u(i,j))**2.0_rp-1.0_rp)
             enddo
          enddo

      else
          if (m.gt.1e-3_rp) then
            kappa = (1-sqrt(1-m))/(1+sqrt(1-m))
          else
            x=m/4.0_rp
            kappa = 132.0_rp*x**6.0_rp + 42.0_rp*x**5.0_rp + &  
                     14.0_rp*x**4.0_rp +  5.0_rp*x**3.0_rp + &  
                     2.0_rp*x**2.0_rp +  1.0_rp*x**1.0_rp
          endif
        
          mu=kappa**2.0_rp

          do i = 1, ud1
             do j = 1, ud2
                v(i,j)=u(i,j)/(1.0_rp+kappa)
             enddo
          enddo

          call jacobi_elliptic_matrix(v, mu , sn1, cn1, dn1, flag=.true.)

          do i = 1, ud1
             do j = 1, ud2
                  denom(i,j) = (1.0_rp+kappa*sn1(i,j)**2.0_rp)
                  sn(i,j) = (1.0_rp+kappa)*sn1(i,j)/denom(i,j)
                  cn(i,j) = cn1(i,j)*dn1(i,j)/denom(i,j)
                  dn(i,j) = (1.0_rp-kappa*sn1(i,j)**2)/denom(i,j)
             enddo
          enddo
      endif


      do i = 1, ud1
         do j = 1, ud2
           if(high(i,j))then
             dn(i,j) = cmplx(0.0,1.0,kind=rp)*cn(i,j)/sn(i,j)
             cn(i,j) = cmplx(0.0,1.0,kind=rp)*dn(i,j)/(sqrt(m)*sn(i,j))
             sn(i,j) = -1.0_rp/sqrt(m)*sn(i,j)
           endif
         enddo
      enddo
      
    end subroutine jacobi_elliptic_matrix

    !-----------------------------------------------------------------------
    !> 
    !>
    !>  SCALAR VERSION OF THE SUBROUTINE
    !>
    !> @author  Alfonso Santiago
    !> @date    2019-JAN-30
    !> @brief   Jacobi Elliptic functions for complex arguments
    !> @details Returns the values of the Jacobi elliptic functions evaluated
    !>          at complex argument U, and parameter M=exp(-2*pi*L), 0<L<inf.
    !>          Recall that M=k**2 where k is thee elliptic modulus.
    !>          U may be a matrix, L must be a scalar. The entries of U
    !>          are expected to lie within the rectangle |Re(u)|<K,
    !>          0<Im(U)<Kp, where K and Kp were computed from the
    !>          elliptic integral.
    !>         
    !>         This algorithm is the descending Landen transformation
    !>         described in L. Howell's PhD thesis from MIT. Additional
    !>         formulas from Gradshteyn and Ryzhik, 5th ed., and  
    !>         Abbramowitz and Stegun
    !>         
    !>         This subroutine was adapted  from Toby Driscoll's
    !>         ellipjc.m.        
    !>         
    !>         
    !-----------------------------------------------------------------------

    recursive subroutine jacobi_elliptic_scalar(U_in, L , sn, cn, dn, flag)
      implicit none
      complex(rp), intent(in)            :: U_in ! U is a complex scalar
      real(rp), intent(in)               :: L      ! Must be a scalar
      logical(lg), optional, intent(in)  :: flag
     
      complex(rp), intent(out)              :: sn, cn, dn
      complex(rp)                           :: snh, cnh, dnh
      complex(rp)                           :: u
      complex(rp)                           :: v
      complex(rp)                           :: sn1,cn1,dn1
      complex(rp)                           :: denom

      real(rp),parameter                 :: eps=1e-08_rp
      real(rp)                           :: m
      real(rp)                           :: K,Kp
      logical(lg)                        :: high
      real(rp)                           :: kappa, x, mu

      u=U_in
      sn=cmplx(0.0,0.0,kind=rp)
      cn=cmplx(0.0,0.0,kind=rp)
      dn=cmplx(0.0,0.0,kind=rp)
      high=.false.

      if(present(flag))then
          high=.false.
          m=L
      else
          call elliptic_integral(L,K,Kp)
             
          if(aimag(u).gt.Kp/2.0_rp) then
            high = .true.
            u=cmplx(0.0,1.0,kind=rp)*Kp-u
          endif
          m=exp(-2.0_rp*pi*L)
      endif


      if (m.lt.4.0_rp*eps) then
          sn=sin(u)+m/4.0_rp*(sin(u)*cos(u)-u)*cos(u)
          cn=cos(u)+m/4.0_rp*(-sin(u)*cos(u)+u)*sin(u)
          dn=1.0_rp+m/4.0_rp*(cos(u)**2.0_rp-sin(u)**2.0_rp-1.0_rp)

      else
          if (m.gt.1e-3_rp) then
            kappa = (1.0_rp-sqrt(1.0_rp-m))/(1.0_rp+sqrt(1.0_rp-m))
          else
            x=m/4.0_rp
            kappa = 132.0_rp*x**6.0_rp + 42.0_rp*x**5.0_rp + &  
                     14.0_rp*x**4.0_rp +  5.0_rp*x**3.0_rp + &  
                     2.0_rp*x**2.0_rp +  1.0_rp*x**1.0_rp
          endif
        
          mu=kappa**2.0_rp
          v=u/(1.0_rp+kappa)

          call jacobi_elliptic_scalar(v, mu , sn1, cn1, dn1, .true.)

          denom = (1.0_rp+kappa*sn1**2.0_rp)
          sn = (1.0_rp+kappa)*sn1/denom
          cn = cn1*dn1/denom
          dn = (1.0_rp-kappa*sn1**2.0_rp)/denom
      endif


      if(high)then
        snh=sn
        cnh=cn
        dnh=dn
        sn = -1.0_rp/(sqrt(m)*snh)
        cn = cmplx(0.0_rp,1.0_rp,kind=rp)*dnh/(sqrt(m)*snh) 
        dn = cmplx(0.0_rp,1.0_rp,kind=rp)*cnh/snh  
      endif

    end subroutine jacobi_elliptic_scalar
    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-JULY-17
    !> @brief   
    !> @details 
    !>         
    !>         
    !>         
    !-----------------------------------------------------------------------
    subroutine solve_fractional_diffusion_time_step()
        use def_domain, only: npoin, vmass
        use mod_solver, only : solver_parallel_SpMV
        use mod_matrix, only: matrix_scaling_CSR
        use def_solver, only: solve_sol

        implicit none

        integer(ip)                            :: i
        real(rp), allocatable, dimension(:)    :: rhsid, elmag_real
        !real(rp), allocatable, dimension(:)    :: THE_TRUE_AND_ONLY_ONE_ELMAG   ! This is here so the compilation doesn't blows
        complex(rp), allocatable, dimension(:) :: elmag_auxi_cx, elmag_cx, aux_vec_cx
        type(soltyp), pointer, dimension(:)    :: solve_local


        complex(rp), pointer, dimension(:) :: amatr_cx ! USE POINTER AND SUBROUTINE MEMALLOC
        real(rp), pointer, dimension(:)    :: xx       ! USE POINTER AND SUBROUTINE MEMALLOC

        call runend("solve_fractional_diffusion_time_step has bugs! E.g. elmag_cx is not allocated")
        allocate(elmag_cx(1))      ! TODO: added these to eliminate a warning that the allocation is missing, remove when necessary
        allocate(elmag_auxi_cx(1)) !
        allocate(elmag_real(1))    !
        allocate(rhsid(1))         !


        allocate(solve_local(1))
        solve_local(1) = solve_sol(1)
        solve_local(1) % kfl_cmplx = 1

        do i=1,size(amatr,KIND=ip)
          amatr_cx(i) = amatr_auxi_exm(i)
          amatr(i) = amatr_auxi_exm(i)
        enddo

        !Build rhsid
        do i=1,npoin
            rhsid(i)=vmass(i)*0.0_rp ! * 0.0_rp should be replaced by xioni
        enddo

        elmag_cx=0.0_rp
        do  i=1,nintp
           !Compute quadrature weights for the current time step
           quadr_weigh(i)=f(quadr_nodes(i))*dxidt(i)
            
           ! Build xi(r) * M
           aux_vec_cx=quadr_nodes(i)*vmass

           ! Build xi(r)*M - amatr
           call modify_amatr_diag(-1.0_rp, amatr_cx, aux_vec_cx)

           ! Solve the following system(solve_local(1), amatr, rhsid, elmag_auxi_cx)
           ! WRITE THE LINE OF CODE

           elmag_cx=elmag_cx + elmag_auxi_cx*quadr_weigh(i)/quadr_nodes(i)
        enddo

         elmag_real = -4.0_rp*ellipK*sqrt(Mseigenmax*MSeigenmin)*aimag(elmag_cx)/(k*pi*nintp)




        call solver_parallel_SpMV(solve_sol(1),amatr,elmag_real,xx,OPENMP=.true.)

        
        do i=1,npoin
          elmag_real(i)=xx(i)/vmass(i)
        enddo


       ! Aqui va la siguiente linea (que para mi no tiene sentido):
       !
       ! e=ones(size(elmag_real))
       !
       ! EL_UNICO_Y_VERDADERO_ELMAG = elmag_real+(e'*M*(xioni-elmag_real)/(e'*M*e))*e
       !
       ! Y digo que no tiene sentido porque:
       !     -> E es una matriz densa llena de "1"
       !     -> E' no tiene sentido si tengo en cuenta lo anterior
       !     -> En la division, las dimensiones no dan bien.
       !     -> Lo que hay después del "+" no tiene sentido
       !     -> Quien conyo hace operaciones con una matriz densa llena de "1"
       !     
       !




        contains
            function f(z) result(r)                                             
              implicit none
              complex(rp), intent(in)    :: z
              complex(rp)                :: r
              real(rp)                   :: dt

              dt=dtime
              r = 1.0_rp/(1.0_rp+dt*coeff*z)

            endfunction
            !-----------------------------------------------------------------------
            !> 
            !> @author  Alfonso Santiago
            !> @date    2019-JULY-17
            !> @brief   creates amatr to send to solver by asigning a sign to it and
            !>          multiplying, substracting or dividing a vector if size npoin.
            !> @details creates amatr from amatr_aux that holds the stiffness matrix
            !>          the exmedi problem.
            !>         
            !>         
            !>         
            !-----------------------------------------------------------------------
            subroutine modify_amatr_diag(amatr_sign, amatr_cx, vector_cx)

                use def_domain, only: npoin
                use def_domain, only: r_dom !, c_dom ! CSR matrix indexes
                use def_solver, only: solve_sol

                implicit none
                complex(rp), pointer, dimension(:)        :: amatr_cx
                real(rp), intent(in)                      :: amatr_sign
                complex(rp), dimension(:), intent(inout)  :: vector_cx(:)
                integer(ip)                               :: i, j
                
                if(size(vector_cx, KIND=ip) .ne. npoin) call runend('MOD_FRACTIONAL_DIFFUSION. CREATE_AMATR. INPUT VECTOR NOT NPOIN SIZE')
                !! Obtain sign*S
                !
                do i=1,size(amatr, KIND=ip)
                  amatr_cx(i) = sign(1.0_rp,amatr_sign)*amatr_cx(i)
                enddo

                if( solve_sol(1) % kfl_symme == 1 ) then
                    do i=1_ip,npoin
                      j=r_dom(i+1)-1
                      amatr_cx(i)=amatr_cx(i)+vector_cx(i)
                    enddo
                else
                    call runend('MOD_EXM_FRACTIONAL_DIFFUSION: CODE NOT PREPARED FOR UNSIMMETRIC MATRICES')
                    ! Check matrix_diagonal_CSR in mod_matrix to implement the code.
                endif


            end subroutine modify_amatr_diag

        

    end subroutine solve_fractional_diffusion_time_step
end module mod_exm_fractional_diffusion
!> @}
