!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_partition_sfc_2.f90
!> @author  Ricard Borrell
!> @date    22-10-2020
!> @brief   SFC partitioning
!> @details Paralle partitioning using SFC
!-----------------------------------------------------------------------

module mod_partition_sfc_2

   use def_kintyp,          only : ip,rp,lg,r1p,i1p
   use mod_maths_basic,     only : maths_sfc_ExaIndex

   use mod_memory,          only : memory_alloca
   use mod_memory,          only : memory_deallo
   use mod_parall,          only : par_memor
   use mod_communications,  only : PAR_SUM
   use mod_communications,  only : PAR_MAX
   use mod_communications,  only : PAR_COMM_RANK_AND_SIZE
   use def_mpi
#include "def_mpi.inc"
   implicit none

   public :: partition_sfc_2
   
   private   
    
   character(100), parameter                      :: vacal = "mod_partition_sfc_2"
   real(rp),    parameter  :: zeror               = epsilon(1.0_rp) ! Almost zero...

contains  
   !----------------------------------------------------------------------
   !
   !> @author  Ricard Borrell
   !> @date    22/10/2020
   !> @brief   Partition of the mesh using Space Filling Curves (SFC)
   !> @details Main subroutine
   !
   !----------------------------------------------------------------------

   subroutine partition_sfc_2(lepar,npart_sfc,lenti,nenti,sdim,lweig,lcorr,PAR_COMM)

      integer(ip), pointer,           intent(inout) :: lepar(:)         !< partition rank of each entity 
      integer(ip),                    intent(inout) :: npart_sfc        !< number of partitions
      real(rp),    pointer,           intent(in)    :: lenti(:,:)       !< coordinates of the entities LENTI_(SDIM_,NENTI_)
      integer(ip),                    intent(in)    :: nenti            !< number of entities 
      integer(ip),                    intent(in)    :: sdim             !< space dimension (2D/3D)
      integer(ip), pointer, optional, intent(in)    :: lweig(:)         !< weight of each entitiy LWEIG_(NENTI_)
      real(rp),    pointer, optional, intent(in)    :: lcorr(:)         !< partitition correction coeficients           
      MY_MPI_COMM,          optional, intent(in)    :: PAR_COMM         !< communicator (empty if sequential)


      integer(ip)                                   :: ii, idime, ind,jj,kk
      MY_MPI_COMM                                   :: PAR_COMM_SFC2
      integer(4)                                    :: PAR_MY_RANK_SFC2
      integer(4)                                    :: PAR_MY_SIZE_SFC2
      real(rp)                                      :: min_coord(3)      
      real(rp)                                      :: max_coord(3)      
      real(rp),     pointer                         :: aux_coord(:)
      real(rp)                                      :: totalw, maxw, avgw
      real(rp)                                      :: obj_lb, lb    
      real(rp)                                      :: auxr, objw
      real(rp),     pointer                         :: lweigr(:) 
      real(rp),     pointer                         :: lsegw(:)
      real(rp),     pointer                         :: lparw(:)
      integer(8),   pointer                         :: lsfcid(:)
      integer(8)                                    :: maxIdx
      integer(8),   pointer                         :: lsegm(:)
      integer(8),   pointer                         :: lpart(:)
      integer(8)                                    :: auxi8, ind8
      integer(ip)                                   :: zz

      !0) Inicialitzacions
      if(present(PAR_COMM)) then
         PAR_COMM_SFC2=PAR_COMM
         call PAR_COMM_RANK_AND_SIZE(PAR_COMM_SFC2,PAR_MY_RANK_SFC2,PAR_MY_SIZE_SFC2)
      else   
         PAR_MY_SIZE_SFC2=1_4
      endif
      nullify(lweigr,lsegw,lpart,lparw,lsfcid,lsegm)
      call memory_alloca(par_memor,'LWEIGR',vacal,lweigr,nenti)
      call memory_alloca(par_memor,'LSEGW' ,vacal,lsegw,2*(npart_sfc+1))
      call memory_alloca(par_memor,'LPARW' ,vacal,lparw,npart_sfc+1)
      call memory_alloca(par_memor,'lsfcid',vacal,lsfcid,nenti)
      call memory_alloca(par_memor,'lsegm' ,vacal,lsegm,2*(npart_sfc+1))
      call memory_alloca(par_memor,'lpart' ,vacal,lpart,npart_sfc+1)
      if( .not. associated(lepar) ) then
         call memory_alloca(par_memor,'LEPAR',vacal,lepar,nenti)
      end if

      !1) Calculate bounding box
      nullify(aux_coord)
      allocate(aux_coord(6))
      aux_coord(1:6)=-huge(1.0_rp)*0.1_rp
      if(nenti>0) then
         do idime = 1_ip,sdim
            aux_coord(idime)        =  maxval(lenti(idime,1:nenti))
            aux_coord(3_ip+idime)   = -minval(lenti(idime,1:nenti))
         end do
      endif
      call PAR_MAX(aux_coord,PAR_COMM_SFC2,INCLUDE_ROOT=.true.)

      max_coord(1:3) =  aux_coord(1:3)
      min_coord(1:3) = -aux_coord(4:6)
      max_coord(1:3) =  max_coord(1:3) + (max_coord(1:3)-min_coord(1:3))*0.001_rp + zeror
      min_coord(1:3) =  min_coord(1:3) - (max_coord(1:3)-min_coord(1:3))*0.001_rp - zeror
      deallocate(aux_coord)

      !2) Calcular max weight, average weight and total weight
      totalw=0.0_rp
      maxw=0.0_rp
      do ii = 1_ip,nenti
         
         if(present(lweig)) then
            lweigr(ii) = lweig(ii)
         else
            lweigr(ii) = 1.0_rp
         endif
         totalw = totalw + lweigr(ii)
         if(lweigr(ii) > maxw) maxw = lweigr(ii) 

      enddo
      call PAR_SUM(totalw, PAR_COMM_SFC2,INCLUDE_ROOT=.true.)
      call PAR_MAX(maxw,PAR_COMM_SFC2,INCLUDE_ROOT=.true.)
      avgw = totalw/npart_sfc

      !3) Calcular minim entre 0.99 i granularitat:
      obj_lb = min(0.99_rp,avgw/(avgw+maxw))

      !5) Calcular index per cada element
      zz= 0_ip 
      maxIdx = 0
      do ii = 1,nenti
         call maths_sfc_ExaIndex(lenti(1:sdim,ii),min_coord,max_coord,sdim,lsfcid(ii),maxIdx)
         if(lsfcid(ii)<=0 .or. lsfcid(ii)>maxIdx) then
            print *, "lsfcid: ",lsfcid(ii), maxIdx
            call runend("somthing wrong with lsfcid")
         endif
         if(lsfcid(ii)==5) zz = zz+1;
         auxi8=maxIdx/npart_sfc
         lepar(ii)=ceiling((real(lsfcid(ii),rp)/real(maxIdx,rp))*npart_sfc)
         if(lepar(ii)<0 .or. lepar(ii)> npart_sfc) then
            call runend("somthing wrong with lepar")
         endif
      enddo
      
      !4) Inicialitzar segments
      do ii=1,npart_sfc+1
         ind=2*(ii-1)+1
         lsegm(ind)=0_8
         lsegw(ind)=0.0_rp
         ind = ind+1
         lsegm(ind)=maxIdx
         lsegw(ind)=totalw
         lparw(ii)=0.0_rp
      enddo
      lpart(1)=0_8
      lpart(npart_sfc+1)=maxIdx
      objw = real(totalw,rp)/real(npart_sfc,rp)

      !6) Iterative loop
      jjloop : do jj=1, 100

            lparw=0.0
            !6.1) Calcular nous spliting points e inicialitzar pesos
            do ii=2,npart_sfc
               ind=2*(ii-1)+1
               call linear_inverse(ii,lsegm(ind),  lsegw(ind),   &
                  &                lsegm(ind+1),lsegw(ind+1), & 
                  &                (ii-1)*objw ,lpart(ii))
            enddo

            !6.2) Calcular contribució local en cada particio
            iiloop : do ii=1,nenti
               ind8 = lsfcid(ii)
               if(ind8 == maxIdx ) then
                  lepar(ii)=npart_sfc
                  cycle iiloop
               endif
               kkloop : do kk= lepar(ii)+1,2,-1
                  if(ind8>=lpart(kk)) then
                     exit kkloop
                  else
                     if(ind8>=lpart(kk-1)) then
                        lparw(kk)=lparw(kk) + lweigr(ii)
                        lepar(ii)= kk-1
                        cycle iiloop
                     endif
                  endif
               end do kkloop
               kkloop2 : do kk= lepar(ii)+1,npart_sfc
                  if(ind8<lpart(kk)) then
                     exit kkloop2
                  else
                     if(ind8<lpart(kk+1)) then
                        lparw(kk+1)=lparw(kk+1) + lweigr(ii)
                        lepar(ii)= kk
                        cycle iiloop
                     endif
                  endif
               end do kkloop2
               print *, "lsfcid: ",ind8,maxIdx, lepar(ii),"lpart",lpart
               call runend("partition_sfc_2: something wrong")
            end do iiloop

         !6.3) Reducció
         call PAR_SUM(lparw,PAR_COMM_SFC2,INCLUDE_ROOT=.true.)

         !6.4) calcular lb

         lb = avgw/maxval(lparw)
         if(lb >= obj_lb) exit jjloop

         !6.5) actualitzar borders, sempre mirart-los tots
         
         !6.5.1) Evaluation of accumulated weight
         do ii = 2,npart_sfc+1
            lparw(ii) = lparw(ii) + lparw(ii-1)
         enddo
         !6.5.2) Calculation of new segment
         do ii = 2,npart_sfc
            auxr = objw*(ii-1)
            ind  = 2*(ii-1)+1
            !this is P² -> to be optimized
            do kk = 2, npart_sfc 
               if(lparw(kk) >= lsegw(ind) .and. &
                  &   lparw(kk) <= lsegw(ind+1)) then

                  if(lparw(kk)> auxr .and. lpart(kk)<lsegm(ind+1)) then
                     lsegm(ind+1)=lpart(kk)
                     lsegw(ind+1)=lparw(kk)
                  else if(lpart(kk) > lsegm(ind) .and. lparw(kk)<= auxr) then
                     lsegm(ind)=lpart(kk)
                     lsegw(ind)=lparw(kk)
                  endif

               endif
            enddo
         enddo 

   end do jjloop

   !7/ Deallocate
   call memory_deallo(par_memor,'LWEIGR',vacal,lweigr)
   call memory_deallo(par_memor,'LSEGW' ,vacal,lsegw)
   call memory_deallo(par_memor,'LPARW ',vacal,lparw)
   call memory_deallo(par_memor,'LSFCID',vacal,lsfcid)
   call memory_deallo(par_memor,'LPART' ,vacal,lpart)
   call memory_deallo(par_memor,'LSEGM' ,vacal,lsegm)

end subroutine partition_sfc_2


subroutine linear_inverse(ii,x1,y1,x2,y2,y3,x3)

   integer(ip),intent(in)  :: ii
   integer(8), intent(in)  :: x1,x2
   real(rp),   intent(in)  :: y1,y2,y3
   integer(8), intent(out) :: x3

   real(rp)                :: x1f,x2f,x3f
   x1f=real(x1,rp)
   x2f=real(x2,rp)

   if(y2-y1<1) then
      x3f = 0.5_rp*(x1f+x2f)
   else
      x3f=(x2-x1)*((y3-y1)/(y2-y1))
   endif
   x3 = ceiling(x3f,8)
   x3=x3+x1

end subroutine

! optimitzacio:
!
!    0. Fer bounding boxes tancades
!    1. Afegir nous punts a les communications
!    2. Eliminar punts subdomins tancats de les communicacions 
!    3. Posar-ho en objectes
!    4. Afgir partmeter profunditat i iteracions
!    5. Poar els pesos amb floats 

end module mod_partition_sfc_2


!> @}
