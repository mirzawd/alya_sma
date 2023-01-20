!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_memchk

  use def_parame
  use def_kintyp
  use mod_memory, only : lbytm
  use mod_memory, only : memory_output_info ! For retrocompatibility
  
  !-----------------------------------------------------------------------
  !
  ! Memory: - allocation check
  !         - deallocation check
  !         - reallocation
  !
  !-----------------------------------------------------------------------

  interface memrea
     module procedure merrp1s,&                 ! reals
          &           merip1s,&                 ! integers
          &           merlg1s,&                 ! logical
          &           mercp1s,&                 ! cells
          &           meriip1s,&                ! integers (:,:)
          &           merrrp1s,&                ! reals (:,:)
          &           merrrrp1s                 ! reals (:,:,:)
  end interface

  interface memchk
     module procedure memrp1,memrp2,memrp3,&   ! reals
          &           memrp4,&
          &           memxp1,memxp2,memxp3,&   ! complex
          &           memip1,memip2,memip3,&   ! 4-byte integers 
          &           memi81,memi82,memi83,&   ! 8-byte integers 
          &           mei1p1,mei1p2,&          ! types
          &           mei1pp1,&      
          &           mei2p1,&
          &           mer1p1,mer1p2,&
          &           mer2p1,mer2p2,&
          &           mer3p1,mer3p2,&
          &           mermem,&
          &           memc11,&                   ! characters
          &           memlg1,&                   ! logical
          &           memcel
  end interface

contains

  subroutine merrp1s(nsnew,memor,vanam,vacal,varia) ! reallocate real
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    real(rp),     pointer       :: varia(:),merrp1(:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig
    integer(4)                  :: istat
    nsold=size(varia)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(merrp1(1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,merrp1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,istat)
    do isize=1,nsold
       merrp1(isize)=varia(isize)
    end do
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    if(istat/=0) call memerr(one,vanam,vacal,zero)
    varia=>merrp1
  end subroutine merrp1s

  subroutine merlg1s(nsnew,memor,vanam,vacal,varia) ! reallocate logical
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    logical(lg),  pointer       :: varia(:),merlg1(:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig
    integer(4)                  :: istat
    nsold=size(varia)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(merlg1(1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,merlg1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,zero)
    do isize=1,nsold
       merlg1(isize)=varia(isize)
    end do
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    if(istat/=0) call memerr(one,vanam,vacal,zero)
    varia=>merlg1
  end subroutine merlg1s

  subroutine merip1s(nsnew,memor,vanam,vacal,varia) ! reallocate integer
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  pointer       :: varia(:),merip1(:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig
    integer(4)                  :: istat
    nsold=size(varia)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(merip1(1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,merip1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,zero)
    do isize=1,nsold
       merip1(isize)=varia(isize)
    end do
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    if(istat/=0) call memerr(one,vanam,vacal,zero)
    varia=>merip1
  end subroutine merip1s

  subroutine mercp1s(nsnew,memor,vanam,vacal,varia) ! reallocate cell
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    type(cell),   pointer       :: varia(:),mercp1(:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig
    integer(4)                  :: istat
    nsold=size(varia)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(mercp1(1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,mercp1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,zero)
    do isize=1,nsold
       mercp1(isize)=varia(isize)
    end do
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    if(istat/=0) call memerr(one,vanam,vacal,zero)
    varia=>mercp1
  end subroutine mercp1s

  subroutine meriip1s(nsnew,memor,vanam,vacal,varia) ! reallocate integer(:,:)
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  pointer       :: varia(:,:),meriip1(:,:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig,ndime,idime
    integer(4)                  :: istat
    nsold=size(varia,2)
    ndime=size(varia,1)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(meriip1(1:ndime,1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,meriip1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,zero)
    do isize=1,nsold
       do idime=1,ndime
          meriip1(idime,isize)=varia(idime,isize)
       enddo
    enddo
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    if(istat/=0) call memerr(one,vanam,vacal,zero)
    varia=>meriip1
  end subroutine meriip1s

  subroutine merrrp1s(nsnew,memor,vanam,vacal,varia) ! reallocate reals(:,:)
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    real(rp),  pointer       :: varia(:,:),merrrp1(:,:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig,ndime,idime
    integer(4)                  :: istat
    nsold=size(varia,2)
    ndime=size(varia,1)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(merrrp1(1:ndime,1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,merrrp1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,zero)
    do isize=1,nsold
       do idime=1,ndime
          merrrp1(idime,isize)=varia(idime,isize)
       enddo
    enddo
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    if(istat/=0) call memerr(one,vanam,vacal,zero)
    varia=>merrrp1
  end subroutine merrrp1s

  subroutine merrrrp1s(nsnew,memor,vanam,vacal,varia) ! reallocate reals(:,:,:)
    implicit none
    integer(ip),  intent(in)    :: nsnew
    character(*), intent(in)    :: vanam,vacal
    real(rp),  pointer          :: varia(:,:,:),merrrrp1(:,:,:)
    integer(8),   intent(inout) :: memor(2)
    integer(ip)                 :: nsold,isize,nsbig,ndim1,ndim2,idim1,idim2
    integer(4)                  :: istat
    nsold=size(varia,3)
    ndim2=size(varia,2)
    ndim1=size(varia,1)
    if(nsnew<=nsold)then
       return
    endif
    nsbig=int(1.5d+00*nsnew)
    allocate(merrrrp1(1:ndim1,1:ndim2,1:nsbig),stat=istat)
    call memchk(zero,istat,memor,vanam,vacal,merrrrp1)
    if(.not.associated(varia)) call memerr(one,vanam,vacal,zero)
    do isize=1,nsold
       do idim2=1,ndim2
          do idim1=1,ndim1
             merrrrp1(idim1,idim2,isize)=varia(idim1,idim2,isize)
          enddo
       enddo
    enddo
    call memchk(two,istat,memor,vanam,vacal,varia)
    deallocate(varia,stat=istat)
    !if(istat/=0) call memerrr(one,vanam,vacal,zero)
    varia=>merrrrp1

  end subroutine merrrrp1s

  subroutine memrp1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvari,ivari
    integer(8),   intent(inout) :: memor(2) 
    real(rp)                    :: varia(:)
    if(itask==0) then
       if(istat==0) then
          nvari=int(size(varia),8)
          lbytm=nvari*int(kind(varia),8)
          do ivari=1,nvari
             varia(ivari)=0.0_rp
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'real')
  end subroutine memrp1

  subroutine memrp2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,ivar1,ivar2
    integer(8),   intent(inout) :: memor(2)
    real(rp)                    :: varia(:,:)
    if(itask==0) then
       if(istat==0) then
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          lbytm=int(size(varia),8)*int(kind(varia),8)
          do ivar2=1,nvar2
             do ivar1=1,nvar1
                varia(ivar1,ivar2)=0.0_rp
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'real')
  end subroutine memrp2

  subroutine memrp3(itask,istat,memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,nvar3,ivar1,ivar2,ivar3
    integer(8),   intent(inout) :: memor(2)
    real(rp)                    :: varia(:,:,:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          nvar3=int(size(varia,3),8)
          do ivar3=1,nvar3
             do ivar2=1,nvar2
                do ivar1=1,nvar1
                   varia(ivar1,ivar2,ivar3)=0.0_rp
                end do
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'real')
  end subroutine memrp3

  subroutine memrp4(itask,istat,memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,nvar3,nvar4,ivar1,ivar2,ivar3,ivar4
    integer(8),   intent(inout) :: memor(2)
    real(rp)                    :: varia(:,:,:,:) 
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          nvar3=int(size(varia,3),8)
          nvar4=int(size(varia,4),8)
          do ivar4=1,nvar4
             do ivar3=1,nvar3
                do ivar2=1,nvar2
                   do ivar1=1,nvar1
                      varia(ivar1,ivar2,ivar3,ivar4)=0_rp
                   end do
                end do
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'real')
  end subroutine memrp4

  subroutine memip1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvari,ivari
    integer(8),   intent(inout) :: memor(2)
    integer(4)                  :: varia(:)
    if(itask==0) then
       if(istat==0) then
          nvari=int(size(varia),8)
          lbytm=nvari*int(kind(varia),8)
          do ivari=1,nvari
             varia(ivari)=0_4
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'integer(4)')
  end subroutine memip1

  subroutine memip2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,ivar1,ivar2
    integer(8),   intent(inout) :: memor(2)
    integer(4)                  :: varia(:,:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          lbytm=int(size(varia),8)*int(kind(varia),8)
          do ivar2=1,nvar2
             do ivar1=1,nvar1
                varia(ivar1,ivar2)=0_4
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'integer(4)')
  end subroutine memip2

  subroutine memip3(itask,istat,memor,vanam,vacal,varia)
    !
    ! Integer(4)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,nvar3,ivar1,ivar2,ivar3
    integer(8),   intent(inout) :: memor(2)
    integer(4)                  :: varia(:,:,:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          nvar3=int(size(varia,3),8)
          do ivar3=1,nvar3
             do ivar2=1,nvar2
                do ivar1=1,nvar1
                   varia(ivar1,ivar2,ivar3)=0_4
                end do
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'integer(4)')
  end subroutine memip3

  subroutine mei1p1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(i1p)                   :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%l)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(i1p)')
  end subroutine mei1p1

  subroutine mei1pp1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(i1pp)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(i1pp)                  :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=nsize*ip
          do isize=1,nsize
             varia(isize)%n=0
             nullify(varia(isize)%l)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(i1p)')
  end subroutine mei1pp1

  subroutine mei1p2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(i1p)                   :: varia(:,:)
    integer(ip)                 :: nsize,isiz1,nsiz1,isiz2,nsiz2
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          nsiz1=size(varia,1)
          nsiz2=size(varia,2)
          lbytm=nsize*ip
          do isiz2=1,nsiz2
             do isiz1=1,nsiz1
                nullify(varia(isiz1,isiz2)%l)
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(i1p)')
  end subroutine mei1p2

  subroutine mei2p1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(i2p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(i2p)                   :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%l)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(i2p)')
  end subroutine mei2p1

  subroutine mer1p1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(r1p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(r1p)                   :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%a)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(r1p)')
  end subroutine mer1p1

  subroutine mer1p2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(r1p)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(r1p)                   :: varia(:,:)
    integer(ip)                 :: nsize,isiz1,nsiz1,isiz2,nsiz2
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          nsiz1=size(varia,1)
          nsiz2=size(varia,2)
          lbytm=nsize*ip
          do isiz2=1,nsiz2
             do isiz1=1,nsiz1
                nullify(varia(isiz1,isiz2)%a)
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(r1p)')
  end subroutine mer1p2

  subroutine mer2p1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(r2p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(r2p)                   :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%a)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(r2p)')
  end subroutine mer2p1

  subroutine mer2p2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(r2p)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(r2p)                   :: varia(:,:)
    integer(ip)                 :: nsize,isiz1,nsiz1,isiz2,nsiz2
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          nsiz1=size(varia,1)
          nsiz2=size(varia,2)
          lbytm=nsize*ip
          do isiz2=1,nsiz2
             do isiz1=1,nsiz1
                nullify(varia(isiz1,isiz2)%a)
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(r2p)')
  end subroutine mer2p2

  subroutine mer3p1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(r3p)                   :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%a)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(r3p)')
  end subroutine mer3p1

  subroutine mer3p2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(r3p)                   :: varia(:,:)
    integer(ip)                 :: nsize,isiz1,nsiz1,isiz2,nsiz2
    if(itask==0) then
       if(istat==0) then 
          nsize=size(varia)
          nsiz1=size(varia,1)
          nsiz2=size(varia,2)
          lbytm=nsize*ip
          do isiz2=1,nsiz2
             do isiz1=1,nsiz1
                nullify(varia(isiz1,isiz2)%a)
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(r3p)')
  end subroutine mer3p2

  subroutine mermem(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(elm)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(elm)                   :: varia(:)
    if(itask==0) then
       if(istat==0) then 
          lbytm=int(size(varia),8)*int(ip,8)
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*ip
    end if
    call memory_output_info(memor,vanam,vacal,'type(elm)')
  end subroutine mermem

  subroutine memc11(itask,istat,memor,vanam,vacal,varia)
    !
    ! Character(5)(:)
    !
    implicit none 
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    character(1)                :: varia(:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          varia=""
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'character')
  end subroutine memc11

  subroutine memi81(itask,istat,memor,vanam,vacal,varia)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    integer(8)                  :: nvari,ivari
    integer(8)                  :: varia(:)
    if(itask==0) then
       if(istat==0) then
          nvari=int(size(varia),8)
          lbytm=nvari*int(kind(varia),8)
          do ivari=1,nvari
             varia(ivari)=0_8
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'integer(8)')
  end subroutine memi81

  subroutine memi82(itask,istat,memor,vanam,vacal,varia)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    integer(8)                  :: nvar1,nvar2,ivar1,ivar2
    integer(8)                  :: varia(:,:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          lbytm=int(size(varia),8)*int(kind(varia),8)
          do ivar2=1,nvar2
             do ivar1=1,nvar1
                varia(ivar1,ivar2)=0_8
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'integer(8)')
  end subroutine memi82

  subroutine memi83(itask,istat,memor,vanam,vacal,varia)
    !
    ! Integer(8)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    integer(8)                  :: nvar1,nvar2,nvar3,ivar1,ivar2,ivar3
    integer(8)                  :: varia(:,:,:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          nvar3=int(size(varia,3),8)
          do ivar3=1,nvar3
             do ivar2=1,nvar2
                do ivar1=1,nvar1
                   varia(ivar1,ivar2,ivar3)=0_8
                end do
             end do
          end do         
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'integer(8)')
  end subroutine memi83

  subroutine memlg1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    logical(lg)                 :: varia(:)
    integer(8)                  :: ivari

    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)
          do ivari=1,int(size(varia),8)
             varia(ivari)=.false.
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)
    end if
    call memory_output_info(memor,vanam,vacal,'logical')
  end subroutine memlg1

  subroutine memcel(itask,istat,memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(cell)                  :: varia(:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(ip,8)
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(ip,8)
    end if
    call memory_output_info(memor,vanam,vacal,'logical')
  end subroutine memcel

  subroutine memxp1(itask,istat,memor,vanam,vacal,varia)
    !
    ! Complex(rp)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvari,ivari
    integer(8),   intent(inout) :: memor(2) 
    complex(rp)                 :: varia(:)
    if(itask==0) then
       if(istat==0) then
          nvari=int(size(varia),8)
          lbytm=nvari*int(kind(varia),8)*2
          do ivari=1,nvari
             varia(ivari)=CMPLX(0.0_rp,0.0_rp,kind=rp)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)*2
    end if
    call memory_output_info(memor,vanam,vacal,'cmpx')
  end subroutine memxp1

  subroutine memxp2(itask,istat,memor,vanam,vacal,varia)
    !
    ! Complex(rp)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,ivar1,ivar2
    integer(8),   intent(inout) :: memor(2)
    complex(rp)                 :: varia(:,:)
    if(itask==0) then
       if(istat==0) then
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          lbytm=int(size(varia),8)*int(kind(varia),8)*2
          do ivar2=1,nvar2
             do ivar1=1,nvar1
                varia(ivar1,ivar2)=CMPLX(0.0_rp,0.0_rp,kind=rp)
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)*2
    end if
    call memory_output_info(memor,vanam,vacal,'cmpx')
  end subroutine memxp2

  subroutine memxp3(itask,istat,memor,vanam,vacal,varia)
    !
    ! Complex(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8)                  :: nvar1,nvar2,nvar3,ivar1,ivar2,ivar3
    integer(8),   intent(inout) :: memor(2)
    complex(rp)                 :: varia(:,:,:)
    if(itask==0) then
       if(istat==0) then
          lbytm=int(size(varia),8)*int(kind(varia),8)*2
          nvar1=int(size(varia,1),8)
          nvar2=int(size(varia,2),8)
          nvar3=int(size(varia,3),8)
          do ivar3=1,nvar3
             do ivar2=1,nvar2
                do ivar1=1,nvar1
                   varia(ivar1,ivar2,ivar3)=CMPLX(0.0_rp,0.0_rp,kind=rp)
                end do
             end do
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-int(size(varia),8)*int(kind(varia),8)*2
    end if
    call memory_output_info(memor,vanam,vacal,'cmpx')
  end subroutine memxp3

end module mod_memchk

