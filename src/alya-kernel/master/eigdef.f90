!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine eigdef(itask)
  !-----------------------------------------------------------------------
  !****f* master/eigdef
  ! NAME
  !    eigdef
  ! DESCRIPTION
  !    This subroutine deos the following:
  !    ITASK = 0 ..... Initialize the Eigen solver type
  !    ITASK = 1,2 ... Bridge between modules and parall service
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master 
  use def_solver
  use def_inpout
  use mod_memchk
  implicit none  
  integer(ip), intent(in) :: itask
  integer(ip)             :: ivari,jtask,nvari
  integer(4)              :: istat

  jtask = itask

  if( itask < 0 ) then
     nvari = -itask
     allocate(momod(modul) % eigen(nvari),stat=istat)   ! Allocate memory   
     eigeg     => momod(modul) % eigen
     eigen_sol => momod(modul) % eigen
     jtask = 0
  end if

  if(jtask==0_ip) then
     !
     ! Initialize solver
     !
     do ivari=1,size(eigen_sol)
        eigen_sol(ivari)%wprob     = 'PROBLEM NAME'
        eigen_sol(ivari)%wsolv     = 'SOLVER NAME'
        eigen_sol(ivari)%wprec     = 'PRECONDITIONER NAME'
        eigen_sol(ivari)%ndofn     = 1           ! D.o.f
        eigen_sol(ivari)%ndof2     = 1           ! D.o.f^2
        eigen_sol(ivari)%neiva     = 0           ! 1 system
        eigen_sol(ivari)%neige     = 0           ! number of eigen values requiered
        eigen_sol(ivari)%nzmbt     = 0           ! RHS eigne matrix size
        eigen_sol(ivari)%nzmat     = 0           ! LHS
        eigen_sol(ivari)%nzrhs     = 0           ! RHS 
        eigen_sol(ivari)%kfl_algso = 1           ! Solver: dir
        eigen_sol(ivari)%kfl_facto = 0           ! Factorize matrix
        eigen_sol(ivari)%kfl_massm = 0           ! Diagonal mass matrix
        eigen_sol(ivari)%miter     = 100         ! Max. number of iterations 
        eigen_sol(ivari)%solco     = 1.0d-6      ! Solver tolerance
        eigen_sol(ivari)%itsol(1)  =  huge(1_ip) ! Max. # of solver iterations
        eigen_sol(ivari)%itsol(2)  = -huge(1_ip) ! Min. # of solver iterations
        eigen_sol(ivari)%itsol(3)  = 0           ! Ave. # of solver iterations
        eigen_sol(ivari)%nsolv     = 0           ! # of solves
        eigen_sol(ivari)%lun_cvgei = 0           ! Convergence unit
        eigen_sol(ivari)%lun_solei = 0           ! Output unit
        eigen_sol(ivari)%kfl_preco = 0           ! Preconditioner
        eigen_sol(ivari)%error     = 0           ! error en eigpack
     end do

  else if(jtask==1_ip.or.jtask==2_ip) then
     !
     ! Used for Parall service
     !
     do ivari=1,size(eigen_sol)

        call iexcha(eigen_sol(ivari)%ndofn)
        call iexcha(eigen_sol(ivari)%ndof2)
        call iexcha(eigen_sol(ivari)%neiva)
        call iexcha(eigen_sol(ivari)%neige)
        call iexcha(eigen_sol(ivari)%nzmbt)
        call iexcha(eigen_sol(ivari)%nzmat)
        call iexcha(eigen_sol(ivari)%nzrhs)
        call iexcha(eigen_sol(ivari)%kfl_algso)
        call iexcha(eigen_sol(ivari)%kfl_facto)
        call iexcha(eigen_sol(ivari)%kfl_massm)
        call iexcha(eigen_sol(ivari)%miter)
        call rexcha(eigen_sol(ivari)%solco)
        call iexcha(eigen_sol(ivari)%itsol(1) )
        call iexcha(eigen_sol(ivari)%itsol(2) )
        call iexcha(eigen_sol(ivari)%itsol(3) )
        call iexcha(eigen_sol(ivari)%nsolv    )
        call iexcha(eigen_sol(ivari)%lun_cvgei)
        call iexcha(eigen_sol(ivari)%lun_solei)
        call iexcha(eigen_sol(ivari)%kfl_preco)
        call iexcha(eigen_sol(ivari)%error)

        if(parii==2.and.IMASTER) parch(nparc+1:nparc+30)  = eigen_sol(ivari)%wprob(1:30)
        if(parii==2.and.ISLAVE)  eigen_sol(ivari)%wprob(1:30) = parch(nparc+1:nparc+30)
        nparc=nparc+30
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+30)  = eigen_sol(ivari)%wsolv(1:30)
        if(parii==2.and.ISLAVE)  eigen_sol(ivari)%wsolv(1:30) = parch(nparc+1:nparc+30)
        nparc=nparc+30
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+30)  = eigen_sol(ivari)%wprec(1:30)
        if(parii==2.and.ISLAVE)  eigen_sol(ivari)%wprec(1:30) = parch(nparc+1:nparc+30)
        nparc=nparc+30
        if(nparc>len(parch)) call runend('SOLDEF: TOO MANY CHARACTERS')

     end do

  else if( jtask == 4_ip ) then

     do ivari=1,size(eigen_sol)

        if(eigen_sol(ivari)%kfl_algso==-999) then

           call runend('SOLDEF: NO SOLVER HAS BEEN DEFINED '&
                //'FOR PROBLEM '//trim(eigen_sol(ivari)%wprob)) 

        else

           !eigen_sol(ivari)%nunkn = eigen_sol(ivari)%ndofn*npoin
           eigen_sol(ivari)%ndof2 = eigen_sol(ivari)%ndofn**2
           !eigen_sol(ivari)%nequa = npoin     

           if( eigen_sol(ivari)%kfl_algso == -2 .or. eigen_sol(ivari)%kfl_algso == -3 ) then
              !
              ! Sparce solver symetric!!
              !
              eigen_sol(ivari)%neige = npoin*eigen_sol(ivari)%neiva
              
           else if( eigen_sol(ivari)%kfl_algso >= 1 ) then
              !
              ! Iterative solvers
              !
              eigen_sol(ivari)%neige = npoin*eigen_sol(ivari)%neiva

           else if( eigen_sol(ivari)%kfl_algso == -1 ) then

              eigen_sol(ivari)%neige = npoin*npoin
 
           end if
           !
           ! RHS Eigen matrix
           !
           if( eigen_sol(1)%kfl_massm == 0 ) then
              eigen_sol(ivari)%nzmbt = npoin
           else
              if( eigen_sol(ivari)%kfl_algso == -3 ) then
                 eigen_sol(ivari)%nzmbt = nzsym*eigen_sol(ivari)%ndof2
              else
                 eigen_sol(ivari)%nzmbt = nzdom*eigen_sol(ivari)%ndof2
              end if
           end if

        end if

        eigen_sol(ivari)%nzmat = nzsol
        eigen_sol(ivari)%nzrhs = npoin

        neige = max(neige,eigen_sol(ivari)%neige)
        neiva = max(neiva,eigen_sol(ivari)%neiva)
        nzmbt = max(nzmbt,eigen_sol(ivari)%nzmbt)
        nzmat = max(nzmat,eigen_sol(ivari)%nzmat)
        nzrhs = max(nzrhs,eigen_sol(ivari)%nzrhs)

     end do

  end if

end subroutine eigdef
