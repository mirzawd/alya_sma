!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_precon

  !--------------------------------------------------------------------------------------
  !****f* solite/mod_precon
  ! NAME
  !    mod_precon
  ! DESCRIPTION
  !    This module contains different implementations of the deflated preconditioning
  ! USES
  ! USED BY
  !***
  !--------------------------------------------------------------------------------------

  use def_master, only  :  IMASTER,INOTMASTER,INOTSLAVE,ISEQUEN,ISLAVE,kfl_paral,NPOIN_TYPE,modul
  use def_kintyp, only  :  ip,rp
  use def_solver, only  :  memit
  use mod_communications
  use mod_memchk
  use mod_csrdirx

  implicit none

contains


  subroutine precon1(nbnodes,ngrou,nbvar,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,pp,qq)

  !---------------------------------------------------------------------------------------------------------------------------------------------------
  ! NAME
  !    precon1
  ! DESCRIPTION
  !    This routine performes simple deflation on a complex sparse linear system of algebraic equations A * q = p
  !                                      q = W * Ac^(-1) * W^(T) * p
  ! INPUT ARGUMENTS
  !    NBNODES ..... Number of nodes in the interior of a mesh
  !    NGROU ....... Number of groups for deflation
  !    NBVAR ....... Number of unknowns in each node
  !    INVPR .......
  !    INVPC .......
  !    IL .......... CSR format: index vector for beginning of a row block for lower triangular matrix L derived from the coarse matrix A'
  !    JL .......... CSR format: index vector for column numbers for lower triangular matrix L derived from the coarse matrix A'
  !    LN .......... Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format derived from the coarse matrix A'
  !    IU .......... CSR format: index vector for beginning of a row block for upper triangular matrix U derived from the coarse matrix A'
  !    JU .......... CSR format: index vector for column numbers for upper triangular matrix U derived from the coarse matrix A'
  !    UN .......... Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format derived from the coarse matrix A'
  !    PP .......... RHS complex vector of the system
  ! OUTPUT ARGUMENTS
  !    QQ .......... Vector of complex unknowns of the system
  !---------------------------------------------------------------------------------------------------------------------------------------------------

    implicit none
  !------------------------------- Input Vars ----

    integer(ip), intent(in)  :: nbnodes, ngrou, nbvar

    integer(ip), pointer     :: invpR(:),invpC(:)
    integer(ip), intent(in)  :: iL(*),jL(*)
    complex(rp), intent(in)  :: Ln(*)
    integer(ip), intent(in)  :: iU(*),jU(*)
    complex(rp), intent(in)  :: Un(*)

    complex(rp), intent(in)  :: pp(*)

  !------------------------------- Work Arrays ---

    complex(rp), intent(out) :: pc(*),qc(*)

  !------------------------------- Output Vars ----

    complex(rp), intent(out)  :: qq(*)

    !pc = W^(T) * p
    call wtvectx(nbnodes,ngrou,nbvar,pc,pp)

    !qc = Ac^(-1) * pc
    call CSR_LUsolx(ngrou,nbvar,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc)

    !q = W * qc
    call wvectx(nbnodes,nbvar,qc,qq)

  end subroutine precon1


  subroutine precon2(nbnodes,ngrou,nbvar,an,ja,ia,diag,invdiag,invpR,invpC,iL,jL,Ln,iU,jU,Un,rc,ec,pp,qq)

  !---------------------------------------------------------------------------------------------------------------------------------------------------
  ! NAME
  !    precon2
  ! DESCRIPTION
  !    This routine performes V-cycle (multigrid method) on a complex sparse linear system of algebraic equations A * q = p
  !    First, it performs NITE=5 iterations of some basic iterative method on the system A * q = p, where q_0 = 0 ---> q_5
  !    Then, it performs deflation on the residual equation A * e_5 = r_5, where r_5 = p - A * q_5 ---> e_5 ---> q_5_new = q_5 + e_5
  !    Finally, it performs again NITE=5 iterations of some basic iterative method on the system A * q = p, where q_0 = q_5_new
  ! INPUT ARGUMENTS
  !    NBNODES ..... Number of nodes in the interior of a mesh
  !    NGROU ....... Number of groups for deflation
  !    NBVAR ....... Number of unknowns in each node
  !    AN .......... Sparse complex matrix of the original system in BCSR (Blocked Compressed Sparse Row) format
  !    JA .......... Compressed Sparse format: index vector for column numbers
  !    IA .......... Compressed Sparse format: index vector for beginning of a row block
  !    DIAG ........ Diagonal of the matrix A
  !    INVDIAG ..... Inverse of the diagonal of the matrix A
  !    INVPR .......
  !    INVPC .......
  !    IL .......... CSR format: index vector for beginning of a row block for lower triangular matrix L derived from the coarse matrix A'
  !    JL .......... CSR format: index vector for column numbers for lower triangular matrix L derived from the coarse matrix A'
  !    LN .......... Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format derived from the coarse matrix A'
  !    IU .......... CSR format: index vector for beginning of a row block for upper triangular matrix U derived from the coarse matrix A'
  !    JU .......... CSR format: index vector for column numbers for upper triangular matrix U derived from the coarse matrix A'
  !    UN .......... Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format derived from the coarse matrix A'
  !    PP .......... RHS complex vector of the system
  ! OUTPUT ARGUMENTS
  !    QQ .......... Vector of complex unknowns of the system
  !---------------------------------------------------------------------------------------------------------------------------------------------------

  !-----------------------------------------------------
  !
  !     q = (I + W A'(-1) W^(t)) p
  !
  !-----------------------------------------------------

    implicit none

  !------------------------------- Input Vars ----

    integer(ip), intent(in)  :: nbnodes,ngrou,nbvar

    complex(rp), intent(in)  :: an(*)
    integer(ip), intent(in)  :: ja(*),ia(*)
    integer(ip), pointer     :: invpR(:),invpC(:)
    integer(ip), intent(in)  :: iL(*),jL(*)
    complex(rp), intent(in)  :: Ln(*)
    integer(ip), intent(in)  :: iU(*),jU(*)
    complex(rp), intent(in)  :: Un(*)

    complex(rp), intent(in)  :: diag(*),invdiag(*),pp(*)

  !------------------------------- Work Arrays ---

    complex(rp), intent(out) :: rc(*),ec(*)

  !------------------------------- Output Vars ----

    complex(rp), intent(out) :: qq(*)

  !------------------------------- Local Vars ----

    integer(ip)              :: ii,totnode
    integer(4)               :: istat
!    real(rp)                 :: cpu_smal1,cpu_smal2,cpu_smal
    complex(rp), pointer     :: rr(:),ee(:)

!#ifdef EVENT
! call mpitrace_user_function(1)
!#endif

    totnode = nbnodes * nbvar
    !Allocate memory for working arrays
    allocate(rr(totnode), stat=istat)
    call memchk(0_ip,istat,memit,'RR','precon2',rr)
    allocate(ee(totnode), stat=istat)
    call memchk(0_ip,istat,memit,'EE','precon2',ee)

   !call cputim(cpu_smal1)
   if (INOTMASTER) then


   !Perform NITE=3 iterations of some basic iterative method on A * q = p
   !q_0 = 0
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static)
      do ii = 1,totnode
        qq(ii) = (0.0_rp,0.0_rp)
      enddo
! !$omp end parallel do
      !call jacobi(nbnodes,nbvar,an,ja,ia,invdiag,pp,qq)
      call sor(nbnodes,nbvar,an,ja,ia,diag,invdiag,pp,qq)


    !Compute the residula of the previous system: r_5 = p - A * q_5
    call bcsplx(nbnodes,nbvar,an,ja,ia,qq,rr)
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static)
    do ii = 1,totnode
      rr(ii) = pp(ii) - rr(ii)
    enddo
! !$omp end parallel do


   endif
   call PAR_BARRIER()
   !call cputim(cpu_smal2)
 	 !if (INOTSLAVE) then
 	 !  cpu_smal = cpu_smal2 - cpu_smal1
 	 !  write(*,*)'Time to relax and residuo:',cpu_smal,kfl_paral
   !endif

 !Perform deflation on the residual equation A * e_5 = r_5
 !call cputim(cpu_smal1)


   !rc = W^(T) * r_5
   call wtvectx(nbnodes,ngrou,nbvar,rc,rr)



 call PAR_BARRIER()
 !call cputim(cpu_smal2)
 !if (INOTSLAVE) then
 !  cpu_smal = cpu_smal2 - cpu_smal1
 !	write(*,*)'Time to small:',cpu_smal,kfl_paral
 !endif

 !call cputim(cpu_smal1)
   if (INOTMASTER) then



    !ec = A'^(-1) * rc
    call CSR_LUsolx(ngrou,nbvar,invpR,invpC,iL,jL,Ln,iU,jU,Un,rc,ec)



   endif
 call PAR_BARRIER()
 !call cputim(cpu_smal2)
 !if (INOTSLAVE) then
 !  cpu_smal = cpu_smal2 - cpu_smal1
 !	write(*,*)'Time to solve small:',cpu_smal,kfl_paral
 !endif

 !call cputim(cpu_smal1)
   !e_5 = W * ec
   call wvectx(nbnodes,nbvar,ec,ee)
 call PAR_BARRIER()
 !call cputim(cpu_smal2)
 !if (INOTSLAVE) then
 !  cpu_smal = cpu_smal2 - cpu_smal1
 !	write(*,*)'Time to grow:',cpu_smal,kfl_paral
 !endif

 !call cputim(cpu_smal1)
   if (INOTMASTER) then
    !Make a correction of q_5: q_5_new = q_5 + e_5
! !$omp parallel do                 &
! !$omp            default(shared)  &
! !$omp            private(ii)      &
! !$omp            schedule(static)
    do ii = 1,totnode
      qq(ii) = qq(ii) + ee(ii)
    enddo
! !$omp end parallel do
    !Perform NITE=3 iterations of some basic iterative method on A * q = p



      !q_0 = q_5_new
      !call jacobi(nbnodes,nbvar,an,ja,ia,invdiag,pp,qq)
      call sor(nbnodes,nbvar,an,ja,ia,diag,invdiag,pp,qq)



   endif
  call PAR_BARRIER()
  !call cputim(cpu_smal2)
  !if (INOTSLAVE) then
  !  cpu_smal = cpu_smal2 - cpu_smal1
  !	write(*,*)'Time to finish and relax:',cpu_smal,kfl_paral
  !endif
    !Deallocate memory of working arrays
    call memchk(2_ip,istat,memit,'EE','jacobi',ee)
    deallocate(ee,stat=istat)
    if ( istat /= 0_ip ) call memerr(2_ip,'EE','jacobi',0_ip)
    call memchk(2_ip,istat,memit,'RR','jacobi',rr)
    deallocate(rr,stat=istat)
    if ( istat /= 0_ip ) call memerr(2_ip,'RR','jacobi',0_ip)



  end subroutine precon2

  subroutine jacobi(nbnodes,nbvar,an,ja,ia,invdiag,pp,qq)

  !-----------------------------------------------------------------------------------------------------------
  ! NAME
  !    jacobi
  ! DESCRIPTION
  !    This routine performes NITE=3 Jacobi iterations on a complex sparse linear system of algebraic equations
  !                                                  A * q = p
  ! INPUT ARGUMENTS
  !    NBNODES ..... Number of nodes in the interior of a mesh
  !    NBVAR ....... Number of unknowns in each node
  !    AN .......... Sparse complex matrix of the original system in BCSR (Blocked Compressed Sparse Row) format
  !    JA .......... Compressed Sparse format: index vector for column numbers
  !    IA .......... Compressed Sparse format: index vector for beginning of a row block
  !    INVDIAG ..... Inverse of the diagonal of the matrix A
  !    PP .......... RHS complex vector of the system
  !    QQ .......... Vector of complex unknowns of the system -
  !                  input argument is an initial guess passed to the routine from outside
  ! OUTPUT ARGUMENTS
  !    QQ .......... Vector of complex unknowns of the system
  !-----------------------------------------------------------------------------------------------------------

    implicit none

  !------------------------------- Input Vars ----

    integer(ip), intent(in)     :: nbnodes,nbvar

    complex(rp), intent(in)     :: an(nbvar,nbvar,*)
    integer(ip), intent(in)     :: ja(*),ia(*)

    complex(rp), intent(in)     :: pp(nbvar,*),invdiag(nbvar,*)

  !------------------------------- Output Vars ----

    complex(rp), intent(inout)  :: qq(nbvar,*)

  !------------------------------- Local Vars ----

    integer(ip)                 :: i,ii,jj,col,kk,ll
    integer(4)                  :: istat
!    real(rp)                    :: cpu_smal1,cpu_smal2,cpu_smal
    complex(rp), pointer        :: aux(:,:)

    !Allocate memory for working arrays
    allocate(aux(nbvar,nbnodes), stat=istat)
    call memchk(0_ip,istat,memit,'AUX','jacobi',aux)

  !Jacobi iteration: q_k+1 = D^(-1) * (L + U) * q_k + D^(-1) * p
  !A = D - L - U, where D is the diagonal of A, -L is strictly lower part of A, -U is strictly upper part of A

    !Do NITE=3 Jacobi iterations
    DO i = 1,1
    !!call cputim(cpu_smal1)
      !Put an initial approximation, q_0, or previous approximation, q_k, in a working array
      do ii = 1,nbnodes
        do jj = 1,nbvar
          aux(jj,ii) = qq(jj,ii)
          qq(jj,ii)  = (0.0_rp,0.0_rp)
        enddo
      enddo
      !(-L - U) * q_k
      do ii = 1,nbnodes
        do jj = ia(ii),ia(ii+1)-1
          col = ja(jj)
          if (col == ii) then
            do ll = 1,nbvar
              do kk = 1,ll-1
                qq(kk,ii) = qq(kk,ii) + an(kk,ll,jj) * aux(ll,ii)
              enddo
              do kk = ll+1,nbvar
                qq(kk,ii) = qq(kk,ii) + an(kk,ll,jj) * aux(ll,ii)
              enddo
            enddo
          else
            do ll = 1,nbvar
              do kk = 1,nbvar
                qq(kk,ii) = qq(kk,ii) + an(kk,ll,jj) * aux(ll,col)
              enddo
            enddo
          endif
        enddo
      enddo
    !!call cputim(cpu_smal2)
    !!cpu_smal = cpu_smal2 - cpu_smal1
  	!!write(*,*)'Time to jacobi1:',cpu_smal,kfl_paral
    !!call cputim(cpu_smal1)
      call pararx('SLX',NPOIN_TYPE,nbnodes*nbvar,qq)
    !!call cputim(cpu_smal2)
    !!cpu_smal = cpu_smal2 - cpu_smal1
  	!1write(*,*)'Time to jacobi exch:',cpu_smal,kfl_paral
    !!call cputim(cpu_smal1)
      !D^(-1) * (L + U) * q_k + D^(-1) * p
      do ii = 1,nbnodes
        do jj = 1,nbvar
          qq(jj,ii) = invdiag(jj,ii) * (pp(jj,ii) - qq(jj,ii))
        enddo
      enddo
    !!call cputim(cpu_smal2)
    !!cpu_smal = cpu_smal2 - cpu_smal1
  	!!write(*,*)'Time to jacobi2:',cpu_smal,kfl_paral
    ENDDO

    !Deallocate memory of working arrays
    call memchk(2_ip,istat,memit,'AUX','jacobi',aux)
    deallocate(aux,stat=istat)
    if ( istat /= 0_ip ) call memerr(2_ip,'AUX','jacobi',0_ip)

  end subroutine jacobi

  subroutine sor(nbnodes,nbvar,an,ja,ia,diag,invdiag,pp,qq)

  !-----------------------------------------------------------------------------------------------------------
  ! NAME
  !    sor
  ! DESCRIPTION
  !    This routine performes NITE=3 SOR iterations on a complex sparse linear system of algebraic equations
  !                                                  A * q = p
  ! INPUT ARGUMENTS
  !    NBNODES ..... Number of nodes in the interior of a mesh
  !    NBVAR ....... Number of unknowns in each node
  !    AN .......... Sparse complex matrix of the original system in BCSR (Blocked Compressed Sparse Row) format
  !    JA .......... Compressed Sparse format: index vector for column numbers
  !    IA .......... Compressed Sparse format: index vector for beginning of a row block
  !    DIAG ........ Diagonal of the matrix A
  !    INVDIAG ..... Inverse of the diagonal of the matrix A
  !    PP .......... RHS complex vector of the system
  !    QQ .......... Vector of complex unknowns of the system -
  !                  input argument is an initial guess passed to the routine from outside
  ! OUTPUT ARGUMENTS
  !    QQ .......... Vector of complex unknowns of the system
  !-----------------------------------------------------------------------------------------------------------

    implicit none

  !------------------------------- Input Vars ----

    integer(ip), intent(in)     :: nbnodes,nbvar

    complex(rp), intent(in)     :: an(nbvar,nbvar,*)
    integer(ip), intent(in)     :: ja(*),ia(*)

    complex(rp), intent(in)     :: diag(nbvar,*),invdiag(nbvar,*),pp(nbvar,*)

  !------------------------------- Output Vars ----

    complex(rp), intent(inout)  :: qq(nbvar,*)

  !------------------------------- Local Vars ----

    integer(ip)                 :: i,ii,jj,col,kk,ll
    integer(4)                  :: istat
    complex(rp)                 :: omega
    complex(rp), pointer        :: aux1(:,:),aux2(:,:)

    !Allocate memory for working arrays
    allocate(aux1(nbvar,nbnodes), stat=istat)
    call memchk(0_ip,istat,memit,'AUX1','sor',aux1)
    allocate(aux2(nbvar,nbnodes), stat=istat)
    call memchk(0_ip,istat,memit,'AUX2','sor',aux2)

  !SOR iteration: (D-omega*L) * q_k+1 = [omega*U + (1-omega)*D] * q_k + omega * p
  !A = D - L - U, where D is the diagonal of A, -L is strictly lower part of A, -U is strictly upper part of A

    omega = (0.1_rp,0.0_rp)
    !Do NITE=3 SOR iterations
    DO i = 1,3
      !Put an initial approximation, q_0, or previous approximation, q_k, in a working array
      do ii = 1,nbnodes
        do jj = 1,nbvar
          aux1(jj,ii) = qq(jj,ii)
          aux2(jj,ii) = (0.0_rp,0.0_rp)
        enddo
      enddo
      DO ii = 1,nbnodes
        !omega*U * q_k
        do jj = ia(ii),ia(ii+1)-1
          col = ja(jj)
          if (col == ii) then
            do ll = 1,nbvar
              do kk = 1,ll-1
                aux2(kk,ii) = aux2(kk,ii) - omega * an(kk,ll,jj) * aux1(ll,ii)
              enddo
            enddo
          elseif (col > ii) then
            do ll = 1,nbvar
              do kk = 1,nbvar
                aux2(kk,ii) = aux2(kk,ii) - omega * an(kk,ll,jj) * aux1(ll,col)
              enddo
            enddo
          endif
        enddo
        !Solve [D-omega*L] * q_k+1 = aux2
        do jj = ia(ii),ia(ii+1)-1
          col = ja(jj)
          if (col == ii) then
            do ll = 1,nbvar
              qq(ll,ii) = aux2(ll,ii) * invdiag(ll,ii)
              do kk = ll+1,nbvar
                aux2(kk,ii) = aux2(kk,ii) - omega * an(kk,ll,jj) * qq(ll,ii)
              enddo
            enddo
          elseif (col < ii) then
            do ll = 1,nbvar
              do kk = 1,nbvar
                aux2(kk,ii) = aux2(kk,ii) - omega * an(kk,ll,jj) * qq(ll,col)
              enddo
            enddo
          endif
        enddo
      ENDDO
      call pararx('SLX',NPOIN_TYPE,nbnodes*nbvar,qq)
      !q_k+1 = q_k+1 + omega * p * D^(-1) + (1-omega) * q_k
      DO ii = 1,nbnodes
        do jj = 1,nbvar
          qq(jj,ii) = qq(jj,ii) + omega * pp(jj,ii) * invdiag(jj,ii) + ((1.0_rp,0.0_rp)-omega) * aux1(jj,ii)
        enddo
      ENDDO
    ENDDO

    !Deallocate memory of working arrays
    call memchk(2_ip,istat,memit,'AUX2','sor',aux2)
    deallocate(aux2,stat=istat)
    if ( istat /= 0_ip ) call memerr(2_ip,'AUX2','sor',0_ip)

    call memchk(2_ip,istat,memit,'AUX1','sor',aux1)
    deallocate(aux1,stat=istat)
    if ( istat /= 0_ip ) call memerr(2_ip,'AUX1','sor',0_ip)

  end subroutine sor

end module mod_precon


