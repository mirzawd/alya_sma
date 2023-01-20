!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  module mod_iqnls
    use def_kintyp
    use def_master,         only : mem_modul
    use def_master,         only : modul,kfl_paral
    use mod_memory,         only : MEMORY_ALLOCA
    use mod_memory,         only : MEMORY_SIZE
    use mod_memory,         only : MEMORY_ALLOCA_MIN
    use mod_memory,         only : MEMORY_DEALLO
    use mod_communications, only : PAR_SUM
    use mod_communications, only : PAR_BARRIER
    use def_master,         only : IMASTER, INOTMASTER
    use mod_communications, only : PAR_MAX
    use mod_parall,         only : PAR_MY_CODE_RANK
    use mod_operations,     only : OPERATIONS_PARALLEL_VECTOR_L2NORM

    implicit none

    private :: Q_mat
    private :: obtain_Q_mat
    private :: vQ_times_Qaux
    private :: vQ_times_vQ
    private :: vQ_times_vector
    private :: par_norm

    public :: compute_alpha
    public :: QRfiltering

    contains
    !----------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    15/06/2016
    !> @brief   QRfiltering algorithm
    !> @details  This algorithm filters the matrix V to use in compute alpha
    !>           if the column in V is linearly dependent. To know if the column is linearly
    !>           dependent you make a QR decomposition and see if |R(i,i)|<||R||, if it is,
    !>           you delete the column.
    !>
    !>            So the algorithm do the following steps:
    !>                  1) Build the dense matrix V (and W) from the not dense martices
    !>                     currV, currW, histV, histW into the matrices densefiltV and
    !>                     densefiltW.
    !>                   2) Proceed with the matrix QR decomposition in densefiltV.
    !>                   3) Rebuild the matrices densefiltV and densefiltW after knwoing
    !>                      R_{ij} and only if |R(i,i)|>\epsilon*||R||, this is, not linearly
    !>                      dependent.
    !>                   4) update current_tracking and history_tracking depending on the
    !>                      deleted columns in currV, currW, histV and histW.
    !>
    !----------------------------------------------------
    subroutine QRfiltering(currV, histV, currW, histW, current_tracking, history_tracking, epsi, densefiltV, densefiltW )


        implicit none
        real(rp), pointer, intent(inout)         :: currV(:,:), currW(:,:)
        real(rp), pointer, intent(inout)         :: histV(:,:,:), histW(:,:,:)
        integer(ip), intent(inout)               :: current_tracking
        integer(ip), intent(inout)               :: history_tracking(:)
        real(rp), intent(in)                     :: epsi
        real(rp), pointer, intent(inout)         :: densefiltV(:,:), densefiltW(:,:)

        integer(ip)                              :: par_leader_rank
        real(rp)                                 :: par_maxndof
        logical(lg)                              :: IAMHEAD

        integer(ip)                              :: ndof, nrank, nhist
        integer(ip)                              :: max_rank_dense, upd_tracking
        integer(ip)                              :: idof, irank, ihist, iaux, iaux2, iaux3, counter1, counter2, counter3
        integer(ip), pointer                     :: history_tracking_static(:)

        real(rp), pointer                        :: Vdecomp(:,:)
        real(rp), pointer                        :: denseaux(:,:)
        real(rp), pointer                        :: R(:,:)
        real(rp)                                 :: norm_column, Rnorm

        nullify(history_tracking_static)
        nullify(Vdecomp)
        nullify(denseaux)
        nullify(R)

       ndof =size(currV,1)
       nrank=size(currV,2)


       nhist=size(history_tracking)

        if(.not.associated(history_tracking_static)) then
            call memory_alloca(mem_modul(1:2,modul),'history_tracking_static','mod_iqnls', history_tracking_static, nhist)
        endif
        do ihist=1,nhist
            history_tracking_static(ihist)=history_tracking(ihist)
        enddo

       ! Choose the head of the matrix
       !
       par_maxndof=real(ndof,rp)
       call PAR_MAX(par_maxndof,'IN CURRENT TARGET COLOR',par_leader_rank)
       IAMHEAD=.false.
       if(par_leader_rank.eq.PAR_MY_CODE_RANK) IAMHEAD=.true.

!IF(IAMHEAD) then
!WRITE(6,*)'HEAD IS RANK', PAR_MY_CODE_RANK
!WRITE(6,*) 'NDOF, NRANK, NHIST', ndof, nrank, nhist
!WRITE(6,*) '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)'
!WRITE(6,*) 'C', current_tracking
!WRITE(6,*) currV(5,:)
!WRITE(6,*) 'H:', history_tracking
!if(history_tracking(1) .gt. 0) then
!WRITE(6,*) histV(5,:,:)
!endif
!ENDIF


       !Buid the dense matrices
       !
       densefiltV=0.0_rp
       densefiltW=0.0_rp
       iaux=0_ip

         do irank=1,current_tracking-1
           if(INOTMASTER) then
             do idof=1, ndof
               densefiltV(idof,irank)=currV (idof, irank)
               densefiltW(idof,irank)=currW (idof, irank)
             enddo
           endif
           iaux=irank
           max_rank_dense=iaux
         enddo

         do ihist=1, nhist
            if(history_tracking(ihist).ne. 0_ip) then
              iaux2=iaux
              do irank=iaux+1,iaux+history_tracking(ihist)-1
                if(INOTMASTER) then
                  do idof=1, ndof
                    densefiltV(idof,irank) = histV(idof, irank-iaux2, ihist)
                    densefiltW(idof,irank) = histW(idof, irank-iaux2, ihist)
                  enddo
                endif
                iaux=irank
              enddo
              iaux2=iaux
              max_rank_dense=iaux
            endif
         enddo

!IF(IAMHEAD) then
!WRITE(6,*) '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)'
!WRITE(6,*) 'denseV', max_rank_dense
!WRITE(6,*) densefiltV(5,:)
!ENDIF



      if(epsi.ge.0.0_rp) then
         ! QR decompostion
         !
         if(IAMHEAD .and. (max_rank_dense .gt. ndof)) call runend('CIQN FILTERING. QR DECOMPOSITION. NOT ENOUGH DOFS IN HEAD RANK')

         if(.not.associated(R)) then
            call memory_alloca(mem_modul(1:2,modul),'Vdecomp', 'mod_iqnls', Vdecomp,  ndof,           max_rank_dense)
            call memory_alloca(mem_modul(1:2,modul),'denseaux','mod_iqnls', denseaux, ndof,           max_rank_dense)
            call memory_alloca(mem_modul(1:2,modul),'R',       'mod_iqnls', R       , max_rank_dense, max_rank_dense)
         endif


         Vdecomp=0.0_rp
         denseaux=0.0_rp
         if(INOTMASTER) denseaux=densefiltV(1:ndof,1:max_rank_dense)
         R=0.0_rp

         QRDECOMPOSITION: do irank=1, max_rank_dense
    
           if(IAMHEAD) then
             iaux=irank
             Vdecomp(irank:ndof,irank) =denseaux(irank:ndof,irank)
           else
             iaux=1
             Vdecomp(:,irank) = denseaux(:,irank)
           endif

           norm_column=par_norm(Vdecomp(iaux:ndof,irank))

           if(IAMHEAD) Vdecomp(irank,irank) = Vdecomp(irank,irank) - norm_column
          
           norm_column=par_norm(Vdecomp(iaux:ndof,irank))

           do idof=iaux,ndof
             Vdecomp(idof,irank) = Vdecomp(idof,irank) / norm_column
           enddo

           if(irank>size(vdecomp,2)) print*,'PIPI 1'
           if(irank.eq.max_rank_dense) call vQ_times_Qaux(vdecomp(:,irank),irank,denseaux)

        enddo QRDECOMPOSITION


        ! Obtain matrix R
        !
        denseaux=0.0_rp
        if(INOTMASTER) denseaux=densefiltV(1:ndof,1:max_rank_dense)
        do irank=1,max_rank_dense
          call vQ_times_Qaux(vdecomp(:,irank), irank, denseaux(:,:))
        enddo


        Rnorm=0.0_rp
        R=0.0_rp
        if(IAMHEAD) then
          do irank=1,max_rank_dense
            do idof=1,irank
              R(idof,irank)=denseaux(idof,irank)
              Rnorm=Rnorm+R(idof,irank)*R(idof,irank)
            enddo
          enddo
          Rnorm=sqrt(Rnorm)
        endif

!if(IAMHEAD) then
!WRITE(6,*) '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)'
!WRITE(6,*) 'MATRIX R:'
!do idof=1, max_rank_dense
!do irank=1,max_rank_dense
!WRITE(6,'(F10.2)', advance='no') R(idof,irank)
!WRITE(6,'(1a)', advance='no') '|'
!enddo
!WRITE(6,*)
!enddo
!WRITE(6,*) '(-)(-)'
!WRITE(6,*) 'Rnorm:', Rnorm
!WRITE(6,*) '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)'
!endif

        do idof=1,max_rank_dense
          do irank=1,max_rank_dense
            call PAR_SUM(R(idof,irank), 'IN CURRENT TARGET COLOR')
          enddo
        enddo
        call PAR_SUM(Rnorm, 'IN CURRENT TARGET COLOR')
        


         densefiltV=0.0_rp
         densefiltW=0.0_rp
         iaux=0_ip
         iaux2=0_ip
         iaux3=1_ip
         counter1=1_ip
         upd_tracking=current_tracking
         CHECK_CURRENT_COLUMNS: do irank=1, current_tracking-1

           if((abs(R(irank,irank)) .lt. epsi*Rnorm) .and. (upd_tracking .gt. 3_ip) ) then !IF R(i,i)<e*||R||

                currV(:,irank)=0.0_rp
                currW(:,irank)=0.0_rp

             upd_tracking=upd_tracking-1

           else !ELSE R(i,i)<e*||R||
             if(INOTMASTER) then
               do idof=1, ndof
                 densefiltV(idof,counter1)=currV (idof, irank)
                 densefiltW(idof,counter1)=currW (idof, irank)
               enddo
             endif
             counter1=counter1+1
           endif!ENDIF R(i,i)<e*||R||

         enddo CHECK_CURRENT_COLUMNS

         currV=0.0_rp
         currW=0.0_rp
         !do irank=1, iaux
         do irank=1, counter1
          if(INOTMASTER) then
            do idof=1, ndof
             currV(idof,irank)=densefiltV(idof,irank)
             currW(idof,irank)=densefiltW(idof,irank)
            enddo
          endif

          norm_column=par_norm(currV(:,irank))
          if((norm_column.eq.0.0_rp).and.(irank.lt.counter1)) call runend ('CIQN FILTERING: ZERO NORM in currV')
         enddo

         current_tracking=upd_tracking

         counter2=counter1


         CHECK_HISTORY:do ihist=1, nhist
            if(history_tracking_static(ihist).ne. 0_ip) then ! IF THERE ARE HISTORIES

              upd_tracking=history_tracking_static(ihist)
              IN_HISTORY_COLUMNS: do irank=1, history_tracking_static(ihist)-1

                if((abs(R(counter2,counter2)) .lt. epsi*Rnorm) ) then !IF R(i,i)<e*||R||
                      histV(:, irank, ihist) =0.0_rp
                      histW(:, irank, ihist) =0.0_rp
                  upd_tracking=upd_tracking-1
                else !ELSE R(i,i)<e*||R||
                  if(INOTMASTER) then
                    DOFS: do idof=1, ndof
                      densefiltv(idof,counter2) = histv(idof, irank, ihist)
                      densefiltw(idof,counter2) = histw(idof, irank, ihist)
                    enddo DOFS
                  endif
                  counter2=counter2+1
                endif!ENDIF R(i,i)<e*||R||
                
              enddo IN_HISTORY_COLUMNS

              history_tracking(ihist)=upd_tracking
            endif
         enddo CHECK_HISTORY


         max_rank_dense=counter2


         if(INOTMASTER) then
           histV(:,:,:)=0.0_rp
           histW(:,:,:)=0.0_rp
         endif

         counter3=counter1
         do ihist=1, nhist
            if(history_tracking(ihist).ne. 0_ip) then ! IF THERE ARE HISTORIES
              do irank=1,history_tracking(ihist)-1
                if(INOTMASTER) then; do idof=1, ndof
                    histV(idof,irank,ihist) = densefiltV(idof,counter3)
                    histW(idof,irank,ihist)=densefiltW(idof,counter3)
                  enddo; endif
                counter3=counter3+1

                norm_column=par_norm(histV(:,irank,ihist))
                if((norm_column.eq.0.0_rp).and.(irank.lt.(history_tracking(ihist)-1))) call runend ('CIQN FILTERING: ZERO NORM in histV')
              enddo
            endif
         enddo


       endif !END IF epsi.ge.0

!IF(IAMHEAD)THEN
!WRITE(6,*) '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)'
!write(6,*) 'densefiltv:'
!write(6,*) densefiltV(5,:)
!ENDIF

        call memory_deallo(mem_modul(1:2,modul),'history_tracking_static','mod_iqnls', history_tracking_static)
        call memory_deallo(mem_modul(1:2,modul),'Vdecomp',                'mod_iqnls', Vdecomp)
        call memory_deallo(mem_modul(1:2,modul),'denseaux',               'mod_iqnls', denseaux)
        call memory_deallo(mem_modul(1:2,modul),'R',                      'mod_iqnls', R)

    end subroutine QRfiltering




    !----------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    15/06/2016
    !> @brief   Parallel implementation of a portion of the interface quasi newton algorithm
    !> @details 
    !>
    !> This subrotuine obtains the alpha vector for the interface
    !> quasi newton algoritm.
    !> With this alpha:
    !>
    !>     x^k+1 = x^k + W*alpha
    !>
    !> As alpha =f(V,r) tricks can be made in the inside
    !> to avoid buling some matrices and hereafter generate
    !> a more efficient parallel code.
    !>
    !> This algoritm is based in the SEQUENTIAL
    !> economy QR decomposition.
    !>
    !> The overall algorithm in this funciton reads:
    !>
    !>     1. Decompose matrix A in a set of vectors
    !>        v_i where:
    !>         Q=(I - 2*v_1*v_1^T)*(I - 2*v_2*v_2^T)...
    !>
    !>     2. Obtain upper triangular matrix U where:
    !>         Q*U=A
    !> 
    !>     3. Compute Q^T*r, where r is the residue and
    !>        Q=f(v1,v2,v3...) as a matrix vector loop
    !>
    !>     4. Backsubstitute R*alpha=Q^T*r (rhs previously
    !>        computed)
    !>
    !>
    !>                 n
    !>            +--+--+--+
    !>            |  |  |  |
    !>            +--+--+--+
    !>            |  |  |  |                 +-+-+-+
    !>  A(m,n)=   +--+--+--+ m     alpha(n)= | | | |
    !>            |  |  |  |                 +-+-+-+
    !>            +--+--+--+
    !>            |  |  |  |
    !>            +--+--+--+
    !>
    !>
    !>                n
    !>            +--+--+--+       
    !>            |  |  |  |                   n
    !>            +--+--+--+               +--+--+--+  
    !>            |  |  |  |               |  |  |  |
    !>  Q(m,n)=   +--+--+--+ m             +--+--+--+  
    !>            |  |  |  |      U(m,n)=  |  |  |  | n
    !>            +--+--+--+               +--+--+--+  
    !>            |  |  |  |               |  |  |  |
    !>            +--+--+--+               +--+--+--+ 
    !>             
    !
    !
    !----------------------------------------------------
    subroutine compute_alpha(A, residue, mmax, nmax, alpha)



      use mod_maths,          only : maths_backsu

      implicit none

      ! Interface variables
      real(rp), intent (in)               :: A(:,:)
      real(rp), intent (in)               :: residue(:)
      real(rp), intent (out), target      :: alpha(:)
      integer(ip), intent(in)             :: mmax, nmax !! maximum values of the matrix to decompose

      ! Parallel variables
      integer(ip)                     :: par_leader_rank
      real(rp)                        :: par_maxndof
      logical(lg)                     :: IAMHEAD

!      integer(ip)                         :: leader_rank


      ! Pointers
      real(rp), pointer                   :: palpha(:)
      real(rp), pointer                   :: pR(:,:)

      ! Internal variables
      real(rp), pointer                   :: Q(:,:)
      real(rp), pointer                   :: R(:,:)
      integer(ip)                         :: i_col, i, j, iaux
      real(rp)                            :: norm_column
      real(rp), pointer                   :: v(:,:)
      real(rp), pointer                   :: A_aux(:,:)

      real(rp), pointer                   :: vecaux(:)
      real(rp), pointer                   :: Qt_times_r(:)

!      real(rp)                            :: maxdofs ! declarated as real to enforce compatibility with PAR_MAX

      nullify(Q)
      nullify(R)
      nullify(v)
      nullify(A_aux)
      nullify(vecaux)
      nullify(Qt_times_r)


      call memory_alloca(mem_modul(1:2,modul),'R',       'mod_iqnls', R       , nmax,nmax)

  
      if(INOTMASTER .and. (mmax.gt.0_ip))then
        call memory_alloca(mem_modul(1:2,modul),'Q',          'mod_iqnls', Q           , mmax,nmax)
        call memory_alloca(mem_modul(1:2,modul),'v',          'mod_iqnls', v           , mmax,nmax)
        call memory_alloca(mem_modul(1:2,modul),'A_aux',      'mod_iqnls', A_aux       , mmax,nmax)
        call memory_alloca(mem_modul(1:2,modul),'Qt_times_r', 'mod_iqnls', Qt_times_r  , mmax)
        call memory_alloca(mem_modul(1:2,modul),'vecaux',     'mod_iqnls', vecaux      , mmax)
      elseif(IMASTER .or. (mmax.eq.0_ip))then
        call memory_alloca(mem_modul(1:2,modul),'Q',          'mod_iqnls', Q           , 1_ip,nmax)
        call memory_alloca(mem_modul(1:2,modul),'v',          'mod_iqnls', v           , 1_ip,nmax)
        call memory_alloca(mem_modul(1:2,modul),'A_aux',      'mod_iqnls', A_aux       , 1_ip,nmax)
        call memory_alloca(mem_modul(1:2,modul),'Qt_times_r', 'mod_iqnls', Qt_times_r  , 1_ip)
        call memory_alloca(mem_modul(1:2,modul),'vecaux',     'mod_iqnls', vecaux      , 1_ip)
      else
        call runend('something is really wrong')
      endif

      pR => R
      palpha => alpha
      Q=0.0_rp

      A_aux=0.0_rp
      if(INOTMASTER .and. mmax.gt.0_ip) A_aux=A(1:mmax,1:nmax)

      !! IDEAS PARALL
      !!
      !! Elejir el subdominio con mas nodos mojados como el 
      !! "subdominio lider", el subdominio que tiene la cabeza
      !! de la matriz a descomponer.
      !!

      ! Find leader rank: the subdomain that heads the decomposition
      !
!      maxdofs=real(mmax,rp)
!      call PAR_MAX(maxdofs,'IN CURRENT TARGET COLOR','EXCLUDE MASTER',leader_rank)
       par_maxndof=real(mmax,rp)
       call PAR_MAX(par_maxndof,'IN CURRENT TARGET COLOR',par_leader_rank)
       IAMHEAD=.false.
       if(par_leader_rank.eq.PAR_MY_CODE_RANK) IAMHEAD=.true.

       if(IAMHEAD .and. (nmax.gt.mmax)) call runend('CIQN FILTERING. QR DECOMPOSITION. NOT ENOUGH DOFS IN HEAD RANK')
      !!
      !! END IDEAS PARALL

      columns: do i_col=1, nmax 


        ! a_i=A(i_col,i_col:mmax)
        
        v(:,i_col)=0.0_rp
        !
        ! If I'm the leader subdomain, I crop the vector
        ! If Im not the leader, I use the full vector
        !
!        if(PAR_MY_CODE_RANK .eq. leader_rank) then
        if(IAMHEAD) then
          iaux=i_col
          v(i_col:mmax,i_col) = A_aux(i_col:mmax,i_col)
        else
          iaux=1
          if(INOTMASTER) v(:,i_col) = A_aux(:,i_col)
        endif

        ! ||alpha|| = sqrt(sum(A(i)^2))
        !
        norm_column = 0.0_rp
        !
        ! If Im the leader use the cropped vector, If'm not
        ! use the full vector
        !
!        if(PAR_MY_CODE_RANK .eq. leader_rank) then
!!        if(IAMHEAD) then
!!          do i=i_col,mmax
!!              norm_column =norm_column + v(i,i_col) * v(i,i_col)
!!          enddo
!!        else
!!          do i=1,mmax
!!              norm_column =norm_column + v(i,i_col) * v(i,i_col)
!!          enddo
!!        endif
!!        
!!        call PAR_SUM(norm_column,'IN CURRENT TARGET COLOR')
!!        
         norm_column=par_norm(v(iaux:mmax,i_col))

        ! u=a_i-alpha*e_1
        !
        !! IMPORTANT THING HERE. Some algoritmhs (i.e. matlab)
        !! uses a_i + alpha* e_i. This algorithm also converges
        !! but slightly different. The minus sign has been
        !! left here trivially.

        ! If I'm the leader, for sure I will have the i_col
        ! element. So If Im the leader do the substraction
        !
!        if(PAR_MY_CODE_RANK .eq. leader_rank) then
        if(IAMHEAD) then
          v(i_col,i_col) = v(i_col,i_col) - norm_column
        endif

        ! ||u|| = sqrt(sum(u(i)^2))
        !
        norm_column = 0.0_rp
        !
        ! If Im the leder use the cropped vector, 
        ! if not, use the full vector.
        !
!        if(PAR_MY_CODE_RANK .eq. leader_rank) then
!!        if(IAMHEAD) then
!!          do i=i_col,mmax
!!            norm_column =norm_column + v(i,i_col) * v(i,i_col)
!!          enddo
!!        else
!!          do i=1,mmax
!!            norm_column =norm_column + v(i,i_col) * v(i,i_col)
!!          enddo
!!        endif
!!
!!        call PAR_SUM(norm_column,'IN CURRENT TARGET COLOR')
!!        
!!        norm_column = sqrt(norm_column)

        norm_column=par_norm(v(iaux:mmax,i_col))

        if(norm_column.eq.0.0_rp) call runend('COMPUTE_ALPHA: ZERO COLUMN IN ALGORITHM')
        
        ! v=u/||u||
        !
        !
        ! If Im the leder use the cropped vector, 
        ! if not, use the full vector.
        !
!        if(PAR_MY_CODE_RANK .eq. leader_rank) then
!!        if(IAMHEAD) then
!!          do i=i_col,mmax
!!            v(i,i_col)=v(i,i_col)/norm_column   
!!          enddo
!!        else
!!          do i=1,mmax
!!            v(i,i_col)=v(i,i_col)/norm_column   
!!          enddo
!!        endif
        do i=iaux,mmax
          v(i,i_col)=v(i,i_col)/norm_column   
        enddo

        ! Q_i=I-2*v*v^T
        !
        ! Step not necesary, because Q_i is obtained in function Q_mat(v,i_ini,i,j)

        ! A_i+1=Q_i * A
        !
        if(.not.(i_col.eq.nmax))then
           call vQ_times_Qaux(v(:,i_col),i_col,A_aux)
        endif
        
      enddo columns



     
      ! Now we need the upper triangular matrix R
      !
      !    A=Q*R  -> R=Q^T*A
      !
      ! This means
      ! 
      ! Q^T*A=(Q_1*Q_2*Q_3*....Q_n)^T * A
      !
      !     = Q_n^T * ... Q_3^T * Q_2^T * Q_1^T * A
      !
      ! In that way we can multiply from right to left.

      ! First we do Q_1^T *A
      !A_aux=A
      if(INOTMASTER .and. mmax.gt.0_ip) A_aux=A(1:mmax,1:nmax)
      call vQ_times_Qaux(v(:,1), 1_ip, A_aux(:,:))

      ! Now we keep multipling 
      do i_col=2,nmax
        call vQ_times_Qaux(v(:,i_col), i_col, A_aux(:,:))
      enddo
      
      !! Copy the first nxn square matrix different from zero
      !
      R=0.0_rp
      !! IDEAS PARALL
      !! 
      !! este bucle solo lo realizara el lider porque nmax
      !! (el numero mas grande de rank de iteraciones)
      !! es mas chico que su numero de degrees of freedom.
      !! 
      !! END IDEAS PARALL

      !
      ! If I'm the master, I'll have the head of
      ! the matrix, and only I should have the R
      !
!      if(PAR_MY_CODE_RANK .eq. leader_rank) then

      R=0.0_rp
      if(IAMHEAD) then
        !do i=1,nmax
        do j=1,nmax
          !do j=1,nmax
          do i=1,j
            R(i,j)=A_aux(i,j)
          enddo
        enddo
      endif

      ! Now sum all the results of every subdomain
      !
      call PAR_SUM(pR,'IN CURRENT TARGET COLOR')





      ! Finally in the backsubstitution we are going to
      ! compute:
      !
      !      R*alpha = - Q^T * r
      !
      ! So first we are going to compute
      ! the RHS as:
      !
      ! Q^T*r=(Q_1*Q_2*Q_3*....Q_n)^T * r
      !
      !     = Q_n^T * ... Q_3^T * Q_2^T * Q_1^T * r
      !
      ! In that way we can multiply from right to left.

      ! First we do Q_1^T *r
      Qt_times_r =0.0_rp

      Qt_times_r = vQ_times_vector(v(:,1),1_ip,.true.,residue(:))


      ! Now we keep multipling 
      do i_col=2,nmax
          Qt_times_r = vQ_times_vector(v(:,i_col),i_col,.true.,Qt_times_r(:))
      enddo

      Qt_times_r = - Qt_times_r


      !
      ! Everyvody should have R and Qt_times_r in this point

      ! Now we call backsubstitution to obtain alpha
      !
      ! R*alpha = - Q^T * r

      alpha=0.0_rp
      !
      ! If I'm the master, Ive the good R and I've
      ! the Qt_times_r that is going to be backsu
      ! bstituted, I should do the backsu
      !
!      if(PAR_MY_CODE_RANK .eq. leader_rank) then
      if(IAMHEAD) then
        call maths_backsu(R(:,:),alpha(:), Qt_times_r,nmax)
      endif

      call PAR_SUM(palpha, 'IN CURRENT TARGET COLOR')

      call memory_deallo(mem_modul(1:2,modul),'Q',          'mod_iqnls', Q)
      call memory_deallo(mem_modul(1:2,modul),'R',          'mod_iqnls', R)
      call memory_deallo(mem_modul(1:2,modul),'v',          'mod_iqnls', v)
      call memory_deallo(mem_modul(1:2,modul),'A_aux',      'mod_iqnls', A_aux)
      call memory_deallo(mem_modul(1:2,modul),'Qt_times_r', 'mod_iqnls', Qt_times_r)
      call memory_deallo(mem_modul(1:2,modul),'vecaux',     'mod_iqnls', vecaux)

    end subroutine compute_alpha
    !  end subroutine compute_alpha
  !contains
    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    function Q_mat(v,i_ini,i,j) result(q_ij)   
        implicit none

        real(rp), intent(in)        :: v(:)
        integer(ip), intent(in)     :: i_ini,i,j
        real(rp)                    :: q_ij
       
        if (i .lt. i_ini) then
            if(i .eq. j) then
                q_ij = 1.0_rp
            else
                q_ij = 0.0_rp 
            endif
        else
            if(i .eq. j) then
                q_ij = 1.0_rp - 2.0_rp*v(i)*v(j)
            else  
                q_ij = -2.0_rp*v(i)*v(j)
            endif
        endif

    end function Q_mat
    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    function obtain_Q_mat(v,i_ini,max_col) result(Q)   
        implicit none

        real(rp), intent(in)        :: v(:)
        integer(ip), intent(in)     :: i_ini, max_col
        integer(ip)                 :: i, j, m
        real(rp)                    :: Q(size(v,1),size(v,1))
      
        m=size(v)

        
        if(m.eq.0_ip) then
          Q=0.0_rp
          return
        endif

        if ( (max_col .gt. m)     .or. &
             (max_col .lt. i_ini) ) then
              call runend('mod_maths: problem with obtain_Q_mat. Wrong max column number')
        endif
        do i=1,m
            do j=1,max_col
                Q(i,j) = Q_mat(v,i_ini,i,j) 
            enddo
        enddo

    end function obtain_Q_mat

    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    subroutine vQ_times_Qaux(v,i_ini,Qaux)   
        implicit none

        real(rp), intent(in)          :: v(:)
        integer(ip), intent(in)       :: i_ini
        real(rp), intent(inout)       :: Qaux(:,:)
        real(rp), pointer             :: mat_aux(:,:)
        real(rp), pointer             :: vec_aux(:)
        real(rp), pointer             :: pvec_aux(:)
        real(rp)                      :: re_aux
        integer(ip)                   :: mv,mQ,nQ
        integer(ip)                   :: i, j
      
        mv=size(v)
        mQ=size(Qaux,1)
        nQ=size(Qaux,2)
    
        if ( (mv .ne. mQ ) ) call runend('mod_maths: vQ_times_Qaux: wrong dimensions')

        nullify(mat_aux)
        nullify(vec_aux)
        if(.not.associated(mat_aux)) then
            call memory_alloca(mem_modul(1:2,modul),'mat_aux', 'mod_iqnls', mat_aux, mQ, nQ)
            call memory_alloca(mem_modul(1:2,modul),'vec_aux', 'mod_iqnls', vec_aux, nQ)
        endif

        pvec_aux => vec_aux

        ! vQ_times_Qaux is :
        ! 
        !  (I-2*v*v^T) * A
        !
        !  A - 2*v*v^T * A
        !
        ! So I can start with v^T * A
        !
        do i=1,nQ
          re_aux=0.0_rp
          do j=1,mQ
            !! IDEAS PARALL
            !!
            !! en este if habria que agregar " Y si soy el subdominio lider"
            !! ya que i_ini siempre sera menor que el rank de iteraciones
            !! guardadas por lo tanto unicamente correspondera a el lider.
            !! algo asi:
            !!
            !!  if ((j .lt. i_ini) .and. IAMLEADER) then
            !!
            !! END IDEAS PARALL

            if (j .lt. i_ini) then
              re_aux = re_aux + 0.0_rp
            else
              re_aux = re_aux + v(j)*Qaux(j,i) !TODO aca creo que es el segudno indice
            endif
          enddo
          vec_aux(i)=re_aux ! here v^T * A is stored
        enddo

        call PAR_SUM(pvec_aux,'IN CURRENT TARGET COLOR')
        
        ! Now we can multiply v * (v^T * A)
        ! and add the A - 2*(...) in the same
        ! line
        !
        do i=1,mv
          do j=1,nQ
            !! IDEAS PARALL
            !!
            !! En este if pasa lo mismo que en el anterior
            !! solo se hace lo de i_ini si soy el subdomino
            !! lider
            !! 
            !!  if ((j .lt. i_ini) .and. IAMLEADER) then
            !!
            !! END IDEAS PARALL
            if (j .lt. i_ini) then
              mat_aux(i,j) =  Qaux(i,j) - 0.0_rp
            else
              mat_aux(i,j) =  Qaux(i,j) - 2.0_rp * v(i) * vec_aux(j)
            endif
          enddo
        enddo


        do i=1,mv
          do j=1,nQ
            Qaux(i,j)=mat_aux(i,j)
          enddo
        enddo

        call memory_deallo(mem_modul(1:2,modul),'mat_aux', 'mod_iqnls', mat_aux)
        call memory_deallo(mem_modul(1:2,modul),'vec_aux', 'mod_iqnls', vec_aux)

    end subroutine vQ_times_Qaux


    !--------------------------------------
    ! Function that multiplies two Q
    ! special matrices
    !--------------------------------------
    function vQ_times_vQ(v1, i_ini_1, v2, i_ini_2) result(Q)
        implicit none

        real(rp), intent(in)        :: v1(:), v2(:)
        integer(ip), intent(in)     :: i_ini_1, i_ini_2
        real(rp)                    :: Q(size(v1,1),size(v1,1))
        real(rp)                    :: rea_aux
        integer(ip)                 :: m,i,j,k
        
        m=size(v1,1)
        Q=0.0_rp
        


    do i=1,m !swipe in rows
        do j=1,m
            rea_aux=0.0_rp
            do k=1,m !swipe in columns
                rea_aux = rea_aux + Q_mat(v1, i_ini_1, i, k) * Q_mat(v2, i_ini_2, k, j)!contract columns
            enddo
            Q(i,j)=rea_aux
        enddo
    enddo


    end function vQ_times_vQ
    !--------------------------------------
    ! Function that multiplies a special
    ! matrix Q with a vector
    !--------------------------------------
    function vQ_times_vector(v1, i_ini_1, transp, vector) result(vec_res)
        implicit none

        real(rp), intent(in)          :: v1(:)
        integer(ip), intent(in)       :: i_ini_1
        real(rp), intent(in)          :: vector(:)
        logical, intent(in)           :: transp

        real(rp), allocatable, target :: vec_res(:)
        real(rp), pointer             :: pvec_res(:)

        real(rp)                      :: rea_aux
        integer(ip)                   :: m,n
        integer(ip)                   :: i

        m=size(v1,1)
        n=size(vector,1)

        allocate(vec_res(n))
        pvec_res => vec_res

        vec_res=0.0_rp
      

        if((m .ne. n) .and. (m.ne.0_ip) ) then
          call runend('vq_times_vector: vector and matrix with different dimension')
        endif



        !
        ! This subroutine does:
        !
        !   (I - 2*v*v^T)*r
        !
        ! by computing:
        !
        !    r - 2*v*v^T*r
        !
        ! THe solution for the transpose vQ matrix is trivial as:
        !
        !   (I - 2*v*v^T)^T
        !
        !   I - 2 * (v^T)^T * v^T
        !
        !   I - 2 * v*v^T
        !
        ! So finally:
        !
        !   (I - 2*v*v^T)^T = I - 2*v*v^T
        !
        ! wich is understandable because it is symetric
        !

        !
        ! so we start with v^T * r
        !
        rea_aux=0.0_rp
        do i=1,m
          !! IDEAS PARALL
          !!
          !! En este if pasa lo mismo que en el anterior
          !! solo se hace lo de i_ini si soy el subdomino
          !! lider
          !! 
          !!  if ((j .lt. i_ini) .and. IAMLEADER) then
          !!
          !! END IDEAS PARALL
          if(i .lt. i_ini_1) then
              rea_aux=rea_aux + 0.0_rp
          else
              rea_aux=rea_aux + v1(i)*vector(i)
          endif
        enddo

        call PAR_SUM(rea_aux,'IN CURRENT TARGET COLOR')


        !
        ! Now we can mulyiply v * (v^T * r)
        ! and we can take profit of the loop
        ! and add r - 2*(...)
        !
        do i=1,m
          !! IDEAS PARALL
          !!
          !! En este if pasa lo mismo que en el anterior
          !! solo se hace lo de i_ini si soy el subdomino
          !! lider
          !! 
          !!  if ((j .lt. i_ini) .and. IAMLEADER) then
          !!
          !! END IDEAS PARALL
          if(i .lt. i_ini_1) then
            vec_res(i) = vector(i) - 0.0_rp
          else
            vec_res(i) = vector(i) - 2.0_rp * v1(i) * rea_aux 
          endif
        enddo


    end function vQ_times_vector
    !--------------------------------------
    ! Function that computes a norm for
    ! a vector in parallel
    !--------------------------------------
    function par_norm(v) result(norm)
        implicit none

        real(rp), intent(inout)          :: v(:)
        real(rp)                      :: norm
        integer(ip)                   :: i

       do i=1,size(v)
         if(abs(v(i)).lt.1.0E-50_rp) v(i)=0.0_rp
       enddo

       norm=sum(v(:)*v(:))

       call PAR_SUM(norm,'IN CURRENT TARGET COLOR')
       norm = sqrt(norm)

    endfunction par_norm  
  end module mod_iqnls
