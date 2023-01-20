!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @addtogroup Electromechanical coupling
!> @{
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @author  Constantine Butakoff      : 
!> @date    Febreuary, 2021
!> @}
!------------------------------------------------------------------------------
module mod_biofibers
   use def_kintyp,   only                   :  ip, rp, lg
   use def_domain,   only                   :  ndime
   ! ----------------------------------------
   implicit none

   ! ----------------------------------------
   type bio_fibers
      logical(lg)                          :: initialised
      integer(ip)                          :: location
      integer(ip)                          :: field_lng_X
      integer(ip)                          :: field_sht_X
      integer(ip)                          :: field_nrm_X
      real(rp)                             :: vector_lng(3)
      real(rp)                             :: vector_nrm(3)
      real(rp)                             :: vector_sht(3)
      real(rp)                             :: basis(3,3)
      real(rp),   pointer                  :: lng_X(:,:), nod_lng(:,:,:)
      real(rp),   pointer                  :: sht_X(:,:), nod_sht(:,:,:)
      real(rp),   pointer                  :: nrm_X(:,:), nod_nrm(:,:,:)
      real(rp),   pointer                  :: nod_F(:,:,:)
      real(rp),   pointer                  :: dummy(:,:)
      contains
         procedure,   pass                 :: init => init_biofibers
         procedure,   pass                 :: destroy => delete_biofibers
         procedure,   pass                 :: manage_arrays => biofib_extra_arrays_work
         procedure,   pass                 :: allocate_memory => biofib_allocate_memory
         procedure,   pass                 :: deallocate_memory => biofib_deallocate_memory
         procedure,   pass                 :: set_data => biofib_init_from_parameters
         procedure,   pass                 :: read_data => biofib_read_data
         procedure,   pass                 :: send_data => biofib_send_data
         procedure,   pass                 :: update_fibers_at_nodes => biofib_update_fibers_at_nodes
         procedure,   pass                 :: get_reference_fibers_at_gp => biofib_get_reference_fibers_at_gp
         procedure,   pass                 :: get_reference_basis_at_gp => biofib_get_reference_basis_at_gp
         procedure,   pass                 :: get_current_fibers_at_gp => biofib_get_current_fibers_at_gp
         procedure,   pass                 :: get_reference_fibers_at_nodes => biofib_compute_reference_fibers_at_nodes
         procedure,   pass                 :: get_current_fibers_at_nodes => biofib_compute_current_fibers_at_nodes
         procedure,   pass                 :: assign_fields => biofib_point_fields
         procedure,   pass                 :: compute_deformation_gradient_at_GP => biofib_compute_deformation_gradient_at_GP
   end type  
   ! ----------------------------------------
   integer(ip),  parameter                 :: NOT_DEFINED     = 0_ip
   integer(ip),  parameter                 :: FIBER_NODAL     = 1_ip
   integer(ip),  parameter                 :: FIBER_ELEMENTAL = 2_ip
   integer(ip),  parameter                 :: FIBER_VECTOR    = 3_ip
   integer(ip),  parameter                 :: FIBER_BASIS     = 4_ip
   ! ----------------------------------------
   logical(lg)                             :: kfl_biofibers = .false.
   logical(lg),            protected       :: kfl_reference_fibers_at_nodes = .false.
   ! ----------------------------------------
   type(bio_fibers), public                :: biofibers
   ! ----------------------------------------
   interface biofib_point_nodal_fibers
      module procedure biofib_point_nodal_fibers_all, &
           &           biofib_point_nodal_fibers_wtime
   end interface biofib_point_nodal_fibers
   ! ----------------------------------------
   public  :: biofib_display_warning
   private :: biofib_read_data, biofib_send_data 
   ! ----------------------------------------
   contains

      ! ----------------------------------------
      !
      ! Public methods
      !
      ! ----------------------------------------
      subroutine biofib_init_from_parameters( self, wfiber, kfl )
         use mod_messages,         only :  messages_live
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         character(len=*),  intent(in)    :: wfiber
         integer(ip),       intent(in)    :: kfl
         ! -------------------------------
         call messages_live('   SETTING BIO-FIBERS MODEL') 

         ! Check if it is a 3D case
         call check_consistency(0_ip)

         ! Initialize data
         call self % init()
         kfl_biofibers = .true.

         ! Set data according to the one read
         if( wfiber == 'LONGI' )then
            if( kfl < 0_ip )then
               if ( self % location == FIBER_VECTOR ) then
                  call runend("kernel: Mixing fiber definitions is not allowed. Either all have to be field or XALIG/YALIG/ZALIG")
               end if
               self % location = FIBER_NODAL
               self % field_lng_X = -kfl
            else
               self % location = FIBER_VECTOR
               self % vector_lng(1:3) = get_vector( kfl )
            endif
         elseif( wfiber == 'SHEET' )then
            if( kfl < 0_ip )then
               if ( self % location == FIBER_VECTOR ) then
                  call runend("kernel: Mixing fiber definitions is not allowed. Either all have to be field or XALIG/YALIG/ZALIG")
               end if
               self % location = FIBER_NODAL
               self % field_sht_X = -kfl
            else
               self % location = FIBER_VECTOR
               self % vector_sht(1:3) = get_vector( kfl )
            endif
         elseif( wfiber == 'NORMA' )then
            if( kfl < 0_ip )then
               if ( self % location == FIBER_VECTOR ) then
                  call runend("kernel: Mixing fiber definitions is not allowed. Either all have to be field or XALIG/YALIG/ZALIG")
               end if
               self % location = FIBER_NODAL
               self % field_nrm_X = -kfl  
            else
               self % location = FIBER_VECTOR
               self % vector_nrm(1:3) = get_vector( kfl )
            endif
         endif
         ! -------------------------------
         contains
             function get_vector( kfl ) result(v)
                 integer(ip), intent(in) :: kfl
                 real(rp)                :: v(3)
                 select case(kfl)
                     case(1); v(1:3) = (/ 1.0_rp, 0.0_rp, 0.0_rp /)
                     case(2); v(1:3) = (/ 0.0_rp, 1.0_rp, 0.0_rp /)
                     case(3); v(1:3) = (/ 0.0_rp, 0.0_rp, 1.0_rp /)
                 end select
             end function get_vector
         ! -------------------------------
      end subroutine biofib_init_from_parameters


      subroutine biofib_read_data( self )
         use def_domain,           only :  ndime
         use def_master,           only :  intost, coupling
         use def_inpout,           only :  words, param, exists
         use mod_ecoute,           only :  ecoute
         use mod_messages,         only :  messages_live
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
         call messages_live('   READING BIO-FIBERS MODEL')

         ! Check if it is a 3D case
         call check_consistency(0_ip)

         ! Initialize data
         call self % init()
         kfl_biofibers = .true.
   
         ! This function is called without going to the next line. Forcing read of the next line here
         call ecoute('biofib_read_data')

         ! Read the fiber model options
         do while(words(1)/='ENDBI')

            if( words(1) == 'LOCAT' )then
               if( words(2) == 'NODAL' )then
                     self % location = FIBER_NODAL
                     call messages_live('      FIBERS DEFINED AT NODES')

               elseif( words(2) == 'ELEME' )then                  
                     call check_consistency(1_ip)
                     self % location = FIBER_ELEMENTAL 
                     call messages_live('      FIBERS DEFINED AT ELEMENTS')

               endif
               
            elseif( words(1) == 'VECTO' )then
               self % location = FIBER_VECTOR
               self % vector_lng(1:3) = normalize(ndime,param(1:3))
               call messages_live('      FIBERS DEFINED THROUGH A NORMALIZED VECTOR')
            
            elseif( words(1) == 'BASIS' )then
               self % location = FIBER_BASIS
               self % basis(1:3,1) = normalize(ndime,param(1:3))
               self % basis(1:3,2) = normalize(ndime,param(4:6))
               self % basis(1:3,3) = normalize(ndime,param(7:9))
               call messages_live('      FIBERS DEFINED THROUGH A NORMALIZED BASIS')   

            elseif( words(1) == 'FIBER' .or. words(1) == 'LONGI' )then
               self % field_lng_X = int(param(1)) 
               call messages_live('      LONGITUDINAL VECTORS READ FROM FIELD '//intost(self % field_sht_X))
                           
            elseif( words(1) == 'SHEET' )then
               self % field_sht_X = int(param(1))
               call messages_live('      SHEET VECTORS READ FROM FIELD '//intost(self % field_sht_X))

            elseif( words(1) == 'NORMA' )then
               self % field_nrm_X = int(param(1)) 
               call messages_live('      NORMAL VECTORS READ FROM FIELD '//intost(self % field_nrm_X))

            end if

            call ecoute('biofib_read_data')

         enddo

         ! Perform some comprovations
         if( self % initialised )then
            if( self % location == 0_ip )then
               call runend('BIOFIB_READ_DATA: LOCATION ONLY CAN BE NODAL OR ELEMENTAL OR BASIS OR VECTOR')
            elseif( self % location == FIBER_VECTOR )then
               if( maxval(abs(self % vector_lng(:))) < epsilon(0.0_rp) )then
                  call runend('BIOFIB_READ_DATA: THE LONGITUDINAL FIBER VECTOR MUST BE DEFINED') 
               endif
            elseif( self % location == FIBER_BASIS )then
               ! TODO : DEFINE THE BASIS CHECKING
            elseif( self % location == FIBER_ELEMENTAL .or. self % location == FIBER_NODAL  )then
               if( self % field_lng_X == 0_ip )then
                  call runend('BIOFIB_READ_DATA: AT LEAST A FIELD FOR THE FIBERS HAS TO BE DEFINED') 
               endif
            endif
         endif            
         ! -------------------------------
         contains 

            pure function normalize(ndime,xx) result(rr)
               ! ------------------------
               implicit none
               ! ------------------------
               integer(ip), intent(in) :: ndime
               real(rp),    intent(in) :: xx(:)
               real(rp)                :: norm
               real(rp)                :: rr(ndime)
               integer(ip)             :: ii
               ! ------------------------
               norm = 0.0_rp
               do ii = 1, ndime
                  norm = norm + xx(ii)*xx(ii)
               enddo
               if( abs(norm) > 0.0_rp )then
                  do ii = 1, ndime
                     rr(ii) = xx(ii)/sqrt(norm)
                  enddo
               endif
               ! ------------------------
            end function normalize

      end subroutine biofib_read_data


      subroutine biofib_send_data( self )
         ! -------------------------------
         use def_master,            only : ISEQUEN, INOTMASTER
         use mod_exchange,          only : exchange_end, exchange_add, exchange_init
         ! -------------------------------
         implicit none
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
         ! Only continue if it is parall
         if( ISEQUEN ) return

         ! Interchange initialisation 
         call exchange_init()
         call exchange_add( kfl_biofibers )
         call exchange_end()

         ! Only continue if self are defined
         if( kfl_biofibers )then
            call self % init()
            call exchange_init()
            call exchange_add( self % location    )
            call exchange_add( self % field_lng_X )
            call exchange_add( self % field_sht_X )
            call exchange_add( self % field_nrm_X ) 
            call exchange_add( self % vector_lng  )
            call exchange_add( self % vector_sht  )
            call exchange_add( self % vector_nrm  )
            call exchange_add( self % basis       )
            call exchange_end()
         endif
         ! -------------------------------
      end subroutine biofib_send_data

      
      subroutine biofib_allocate_memory( self )
         use def_domain,           only  : nelem, npoin
         use def_master,           only  : gdepo, modul, mem_modul
         use mod_memory,           only  : memory_alloca
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip)                      :: nsize
         ! -------------------------------
         ! Allocate dummy pointer according to fibers time
         if( .not. associated(self%dummy)) then
            if( self % location == FIBER_ELEMENTAL )then
               nsize = max(1_ip,nelem)
            elseif( self % location == FIBER_NODAL )then
               nsize = max(1_ip,npoin)
            else
               nsize = 1_ip
            endif
            call memory_alloca(mem_modul(1:2,modul),'DUMMY','biofibers',self%dummy,ndime,nsize)
         endif

         ! Allocate memory for deformation gradient 
         if( .not. associated(gdepo) )then
            call memory_alloca(mem_modul(1:2,modul),'GDEPO','biofibers',gdepo,ndime,ndime,npoin)
         endif
         ! -------------------------------
      end subroutine biofib_allocate_memory


      subroutine biofib_deallocate_memory( self )
         use def_master,           only  : gdepo, modul, mem_modul
         use mod_memory,            only : memory_deallo
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
          if( associated( self%dummy ) ) call memory_deallo(mem_modul(1:2,modul),'DUMMY','biofibers',self%dummy)
          if( associated( gdepo )      ) call memory_deallo(mem_modul(1:2,modul),'GDEPO','biofibers',gdepo)
         ! -------------------------------
      end subroutine biofib_deallocate_memory


      subroutine biofib_extra_arrays_work( self, wtask )
         use def_master,            only : gdepo
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         character(len=*),  intent(in)    :: wtask
         ! -------------------------------
         ! Allocate extra variables
         if(     wtask == 'ALLOCATE' )then
             call self % allocate_memory()
         elseif( wtask == 'DEALLOCATE' )then
             call self % deallocate_memory()    
         elseif( wtask == 'REDISTRIBUTE' )then
             call self % deallocate_memory()    
             call self % allocate_memory()
         endif
 
         ! Point to variables allocated by Kernel
         self % nod_F => gdepo 

         ! Point to read fields or the dummy variable 
         call biofib_point_fields( self )
         ! -------------------------------
      end subroutine biofib_extra_arrays_work


      subroutine biofib_point_fields( self )
         use def_master,           only  : INOTEMPTY
         use def_domain,           only  : xfiel
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
         if( INOTEMPTY )then
    
            ! Associate fiber direction to the corresponding field
            if( self % field_lng_X > 0_ip )then
               self % lng_X => xfiel(self % field_lng_X) % a(:,:,1)
               call normalize(self % lng_X)
            else
               self % lng_X => self % dummy
            endif
            
            ! Associate sheet direction to the corresponding field
            if( self % field_sht_X > 0_ip )then
               self % sht_X => xfiel(self % field_sht_X) % a(:,:,1)
               call normalize(self % sht_X)
            else
               self % sht_X => self % dummy
            endif
      
            ! Associate normal direction to the corresponding field
            if( self % field_nrm_X > 0_ip )then
               self % nrm_X => xfiel(self % field_nrm_X) % a(:,:,1)
               call normalize(self % nrm_X)
            else
               self % nrm_X => self % dummy
            endif

         end if   
         
         ! -------------------------------
         contains 

            subroutine normalize(array)
               ! ------------------------
               implicit none
               ! ------------------------
               real(rp), intent(inout) :: array(:,:) 
               real(rp)                :: norm
               integer(ip)             :: ii, jj
               ! ------------------------
               do jj = 1, size(array,dim=2)
                  norm = 0.0_rp
                  do ii = 1, size(array,dim=1)
                        norm = norm + array(ii,jj)**2
                  enddo
                  if( abs(norm) > 0.0_rp )then
                        array(:,jj) = array(:,jj)/sqrt(norm)
                  endif
               enddo
               ! ------------------------
            end subroutine normalize      

         ! -------------------------------
      end subroutine biofib_point_fields


      subroutine biofib_get_reference_fibers_at_gp( self, ielem, pnode, pgaus, gp_N, gp_lng, gp_sht, gp_nrm )
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip),       intent(in)    :: ielem, pnode, pgaus
         real(rp),          intent(in)    :: gp_N(pnode,pgaus)
         real(rp),          intent(out)   :: gp_lng(3,pgaus)
         real(rp),          intent(out)   :: gp_sht(3,pgaus)
         real(rp),          intent(out)   :: gp_nrm(3,pgaus)
         ! -------------------------------
         integer(ip)                      :: ii
         real(rp)                         :: el_fiber(3,pnode)
         real(rp)                         :: el_sheet(3,pnode)
         real(rp)                         :: el_norma(3,pnode)
         ! -------------------------------
         ! Initialise
         gp_lng = 0.0_rp 
         gp_sht = 0.0_rp 
         gp_nrm = 0.0_rp

         ! Get fibers at the GP according to the model 
         if( self % initialised )then

            if(     self % location == FIBER_VECTOR )then 
               do ii = 1, pgaus
                   gp_lng(1:3,ii) = self % vector_lng(1:3)
               enddo
               do ii = 1, pgaus
                   gp_sht(1:3,ii) = self % vector_sht(1:3)
               enddo
               do ii = 1, pgaus
                   gp_nrm(1:3,ii) = self % vector_nrm(1:3)
               enddo

            elseif( self % location == FIBER_BASIS )then
               do ii = 1, pgaus
                  gp_lng(1:3,ii) = self % basis(1:3,1)
               enddo
               do ii = 1, pgaus
                  gp_sht(1:3,ii) = self % basis(1:3,2)
               enddo
               do ii = 1, pgaus
                  gp_nrm(1:3,ii) = self % basis(1:3,3)
               enddo

            elseif( self % location == FIBER_ELEMENTAL )then
               if( self % field_lng_X > 0_ip)then
                  do ii = 1, pgaus
                     gp_lng(1:3,ii) = self % lng_X(1:3,ielem)
                  enddo
               else
                  do ii = 1, pgaus
                     gp_lng(1:3,ii) = self % vector_lng(1:3)
                  enddo
               endif
               if( self % field_sht_X > 0_ip)then
                  do ii = 1, pgaus
                     gp_sht(1:3,ii) = self % sht_X(1:3,ielem)
                  enddo
               else
                  do ii = 1, pgaus
                     gp_sht(1:3,ii) = self % vector_sht(1:3)
                  enddo
               endif
               if( self % field_nrm_X > 0_ip)then
                  do ii = 1, pgaus
                     gp_nrm(1:3,ii) = self % nrm_X(1:3,ielem)
                  enddo
               else
                  do ii = 1, pgaus
                     gp_nrm(1:3,ii) = self % vector_nrm(1:3)
                  enddo
               endif
               
            elseif( self % location == FIBER_NODAL )then
               if( self % field_lng_X > 0_ip )then
                  el_fiber = gather_from_alya_to_array(ielem,pnode,self % lng_X)
                  gp_lng = interpolate_and_normalize(pnode,pgaus,gp_N,el_fiber)
               else
                  do ii = 1, pgaus
                     gp_lng(1:3,ii) = self % vector_lng(1:3)
                  enddo
               endif
               if( self % field_sht_X > 0_ip )then
                  el_sheet = gather_from_alya_to_array(ielem,pnode,self % sht_X)
                  gp_sht = interpolate_and_normalize(pnode,pgaus,gp_N,el_sheet)
               else
                  do ii = 1, pgaus
                     gp_sht(1:3,ii) = self % vector_sht(1:3)
                  enddo
               endif 
               if( self % field_nrm_X > 0_ip )then
                  el_norma = gather_from_alya_to_array(ielem,pnode,self % nrm_X)
                  gp_nrm = interpolate_and_normalize(pnode,pgaus,gp_N,el_norma)
               else
                  do ii = 1, pgaus
                     gp_nrm(1:3,ii) = self % vector_nrm(1:3)
                  enddo
               endif

            endif

         else
            ! Default is oriented witht he CSYS_GLOBAL
            do ii = 1, pgaus
               gp_lng(1,ii) = 1.0_rp 
               gp_sht(2,ii) = 1.0_rp 
               gp_nrm(3,ii) = 1.0_rp 
            enddo
         endif
         ! -------------------------------
         contains

            function gather_from_alya_to_array(ielem,nnode,alyavar) result(elemarr)
               use def_domain, only : lnods
               ! ------------------------
               integer(ip), intent(in) :: ielem
               integer(ip), intent(in) :: nnode
               real(rp),    intent(in) :: alyavar(:,:)
               real(rp)                :: elemarr(3,nnode)
               integer(ip)             :: inode, ipoin
               ! ------------------------
               do inode = 1,nnode
                  ipoin = lnods(inode,ielem)
                  elemarr(1:3,inode) = alyavar(1:3,ipoin)
               end do
               ! ------------------------
            end function gather_from_alya_to_array


            function interpolate_and_normalize(pnode,pgaus,gpsha,var_el) result(var_gp)
               ! ------------------------
               integer(ip), intent(in) :: pnode, pgaus
               real(rp),    intent(in) :: gpsha(:,:)
               real(rp),    intent(in) :: var_el(:,:)
               real(rp)                :: var_gp(3,pgaus)
               real(rp)                :: norm
               integer(ip)             :: idof, inode, igaus
               ! ------------------------
               var_gp(:,:) = 0.0_rp
               do igaus = 1, pgaus
                  ! Interpolate
                  do inode = 1, pnode
                        do idof = 1, 3
                           var_gp(idof,igaus) = var_gp(idof,igaus) + gpsha(inode,igaus)*var_el(idof,inode)
                        enddo
                  end do
                  ! Normalise
                  norm = 0.0_rp
                  do idof = 1, 3
                        norm = norm + var_gp(idof,igaus)**2
                  enddo
                  if( abs(norm) > 0.0_rp )then
                        var_gp(:,igaus) = var_gp(:,igaus)/sqrt(norm)
                  endif
               enddo
               ! ------------------------
            end function interpolate_and_normalize

         ! -------------------------------
      end subroutine biofib_get_reference_fibers_at_gp


      subroutine biofib_get_reference_basis_at_gp( self, ielem, pnode, pgaus, gp_N, gp_basis )
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip),       intent(in)    :: ielem, pnode, pgaus
         real(rp),          intent(in)    :: gp_N(pnode,pgaus)
         real(rp),          intent(out)   :: gp_basis(3,3,pgaus)
         ! -------------------------------
         integer(ip)                      :: gp
         real(rp)                         :: gp_lng(3,pgaus)
         real(rp)                         :: gp_sht(3,pgaus)
         real(rp)                         :: gp_nrm(3,pgaus)
         ! -------------------------------
         ! Get each direction
         call biofib_get_reference_fibers_at_gp( self, ielem, pnode, pgaus, gp_N, gp_lng, gp_sht, gp_nrm )
         ! Pack as a basis
         do gp = 1, pgaus
            gp_basis(1:3,1,gp) = gp_lng(1:3,gp)
            gp_basis(1:3,2,gp) = gp_sht(1:3,gp)
            gp_basis(1:3,3,gp) = gp_nrm(1:3,gp)
         enddo
         ! -------------------------------
      end subroutine biofib_get_reference_basis_at_gp


      subroutine biofib_get_current_fibers_at_gp( self, ielem, pnode, pgaus, gp_N, gp_lng_y, gp_sht_y, gp_nrm_y )
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
         integer(ip),       intent(in)    :: ielem, pnode, pgaus
         real(rp),          intent(in)    :: gp_N(pnode,pgaus)
         real(rp),          intent(out)   :: gp_lng_y(3,pgaus)
         real(rp),          intent(out)   :: gp_sht_y(3,pgaus)
         real(rp),          intent(out)   :: gp_nrm_y(3,pgaus)
         real(rp)                         :: gp_lng_X(3,pgaus)
         real(rp)                         :: gp_sht_X(3,pgaus)
         real(rp)                         :: gp_nrm_X(3,pgaus)
         real(rp)                         :: gp_F(3,3,pgaus)
         ! -------------------------------
         ! Get reference fibers at GP
         call biofib_get_reference_fibers_at_gp( self, ielem, pnode, pgaus, gp_N, gp_lng_X, gp_sht_X, gp_nrm_X )

         ! Get deformation gradient at GP
         call biofib_compute_deformation_gradient_at_GP( self, ielem, gp_F) 

         ! Push forward of f=F'*f0 ; s=F'*s0 ; n=F'*n0
         gp_lng_y = MxV(pgaus,gp_F,gp_lng_X)
         gp_sht_y = MxV(pgaus,gp_F,gp_sht_X)
         gp_nrm_y = MxV(pgaus,gp_F,gp_nrm_X)
         ! -------------------------------
         contains

            pure function MxV(pgaus,M,v) result(r)
               ! -------------------------------
               integer(ip),        intent(in) :: pgaus
               real(rp),           intent(in) :: M(3,3,pgaus)
               real(rp),           intent(in) :: v(3,pgaus)
               real(rp)                       :: r(3,pgaus)
               integer(ip)                    :: ii, jj, igaus
               ! -------------------------------
               do igaus = 1, pgaus
                  r(:,igaus) = 0.0_rp
                  do jj = 1, 3
                     do ii = 1, 3
                        r(ii,igaus) = r(ii,igaus) + M(ii,jj,igaus)*v(jj,igaus)
                     enddo
                  enddo 
               enddo
               ! -------------------------------
         end function MxV

      end subroutine biofib_get_current_fibers_at_gp


      subroutine biofib_compute_reference_fibers_at_nodes( self )
         use def_master,  only             : IEMPTY
         use def_domain,  only             : ltype, lnods, elmar, coord, vmass 
         use def_domain,  only             : mgaus, mnode, nnode, ngaus, ndime, nelem, npoin
         use mod_elmgeo,  only             : elmgeo_jacobian_matrix
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip)                      :: ielem, inode, igaus, ipoin, idime
         integer(ip)                      :: pgaus, pnode, pelty
         real(rp)                         :: xmean
         real(rp)                         :: el_coord(ndime,mnode)
         real(rp)                         :: gp_factor, gp_W(mgaus), gp_N(mnode,mgaus), gp_dNde(ndime,mnode,mgaus)
         real(rp)                         :: gp_J_det
         ! -------------------------------
         ! Check dime == 3
         call check_consistency( 0_ip )

         ! Empty MPI?
         if( IEMPTY ) return

         ! Previous computed?
         if( kfl_reference_fibers_at_nodes ) return

         if(     self % location == FIBER_ELEMENTAL )then
            ! Check no eccoupling
            call check_consistency( 1_ip )

            ! Initialise
            call init_array(self % nod_lng(:,:,3))
            call init_array(self % nod_sht(:,:,3))
            call init_array(self % nod_nrm(:,:,3))

            ! Average the fibers at nodes from integration points
            do ielem = 1, nelem
               pelty = ltype(ielem)

               if( pelty > 0_ip )then
                  pnode = nnode(pelty)
                  pgaus = ngaus(pelty)

                  ! Gather coordinates
                  do inode = 1, pnode
                     ipoin = lnods(inode,ielem)
                     do idime = 1, ndime
                        el_coord(idime,inode) = coord(idime,ipoin)
                     end do
                  end do

                  ! Get interpolation
                  gp_W(1:pgaus) = elmar(pelty) % weigp(1:pgaus)
                  gp_N(1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
                  gp_dNde(1:ndime,1:pnode,1:pgaus) = elmar(pelty) % deriv(1:ndime,1:pnode,1:pgaus)

                  ! From GP to nodes
                  do igaus = 1, pgaus
                     call elmgeo_jacobian_matrix(ndime,pnode,el_coord,gp_dNde(1:ndime,1:pnode,igaus),gp_J_det)
                     xmean = gp_W(igaus) * gp_J_det
                     do inode = 1, pnode    
                        ipoin = lnods(inode,ielem)
                        gp_factor = gp_N(inode,igaus) * xmean
                        self % nod_lng(1:3,ipoin,3) = self % nod_lng(1:3,ipoin,3) + gp_factor * self % lng_X(1:3,ielem)
                        self % nod_sht(1:3,ipoin,3) = self % nod_sht(1:3,ipoin,3) + gp_factor * self % sht_X(1:3,ielem)
                        self % nod_nrm(1:3,ipoin,3) = self % nod_nrm(1:3,ipoin,3) + gp_factor * self % nrm_X(1:3,ielem)
                     end do
                  end do

               endif
            end do

            ! Interchange with other workers
            call rhsmod(ndime,self % nod_lng(:,:,3))
            call rhsmod(ndime,self % nod_sht(:,:,3))
            call rhsmod(ndime,self % nod_nrm(:,:,3))

            ! Average
            do ipoin=1, npoin
               self % nod_lng(1:3,ipoin,3) = self % nod_lng(1:3,ipoin,3) / vmass(ipoin)
               self % nod_sht(1:3,ipoin,3) = self % nod_sht(1:3,ipoin,3) / vmass(ipoin)
               self % nod_nrm(1:3,ipoin,3) = self % nod_nrm(1:3,ipoin,3) / vmass(ipoin)
            end do

         elseif( self % location == FIBER_NODAL )then
            do ipoin = 1, npoin
               self % nod_lng(1:3,ipoin,3) = self % lng_X(1:3,ipoin)
               self % nod_sht(1:3,ipoin,3) = self % sht_X(1:3,ipoin)
               self % nod_nrm(1:3,ipoin,3) = self % nrm_X(1:3,ipoin)
            enddo
            
         elseif( self % location == FIBER_VECTOR )then
            do ipoin = 1, npoin
               self % nod_lng(1:3,ipoin,3) = self % vector_lng(1:3)
               self % nod_sht(1:3,ipoin,3) = self % vector_sht(1:3)
               self % nod_nrm(1:3,ipoin,3) = self % vector_nrm(1:3)
            enddo

         elseif( self % location == FIBER_BASIS )then
            do ipoin = 1, npoin
               self % nod_lng(1:3,ipoin,3) = self % basis(1:3,1)
               self % nod_sht(1:3,ipoin,3) = self % basis(1:3,2)
               self % nod_nrm(1:3,ipoin,3) = self % basis(1:3,3)
            enddo
         
         endif

         ! Set fibers 
         kfl_reference_fibers_at_nodes = .true.

      end subroutine biofib_compute_reference_fibers_at_nodes


      subroutine biofib_compute_current_fibers_at_nodes( self )
         use def_domain,  only             : npoin
         use def_master,  only             : displ, IEMPTY
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip)                      :: ipoin
         ! -------------------------------
         ! Check dime == 3
         call check_consistency( 0_ip )
  
         if( IEMPTY ) return

         ! If displ is not associated, F = I, which means that F_t+1 = F_t = F_0
         if( associated(displ) )then
            ! Get current deformation gradient at nodes
            call biofib_compute_deformation_gradient_at_nodes( self )

            ! Push-forward the fibers 
            do ipoin = 1, npoin
               self % nod_lng(1:3,ipoin,1) = MxV(self % nod_F(1:3,1:3,ipoin),self % nod_lng(1:3,ipoin,3))
               self % nod_sht(1:3,ipoin,1) = MxV(self % nod_F(1:3,1:3,ipoin),self % nod_sht(1:3,ipoin,3))
               self % nod_nrm(1:3,ipoin,1) = MxV(self % nod_F(1:3,1:3,ipoin),self % nod_nrm(1:3,ipoin,3))
            enddo

         else
              self % nod_lng(1:3,1:npoin,1) = self % nod_lng(1:3,1:npoin,3)
              self % nod_sht(1:3,1:npoin,1) = self % nod_sht(1:3,1:npoin,3)
              self % nod_nrm(1:3,1:npoin,1) = self % nod_nrm(1:3,1:npoin,3)
         endif

         contains

            pure function MxV(M,v) result(r)
               ! -------------------------------
               real(rp),           intent(in) :: M(3,3)
               real(rp),           intent(in) :: v(3)
               real(rp)                       :: r(3)
               integer(ip)                    :: ii, jj
               ! -------------------------------
               r(:) = 0.0_rp
               do jj = 1, 3
                  do ii = 1, 3
                     r(ii) =  r(ii) + M(ii,jj)*v(jj)
                  enddo
               enddo
               ! -------------------------------
            end function MxV

      end subroutine biofib_compute_current_fibers_at_nodes


      subroutine biofib_compute_deformation_gradient_at_nodes( self )
         use def_master,  only             : INOTEMPTY, displ
         use def_domain,  only             : ltype, lnods, elmar, coord, vmass
         use def_domain,  only             : mgaus, mnode, nnode, ngaus, ndime, nelem, npoin
         use mod_elmgeo,  only             : elmgeo_cartesian_derivatives_jacobian
         use mod_array_operations,  only   : array_operations_initialization
         use mod_matrix,  only             : matrix_assemble_element_RHS 
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip)                      :: ielem, ipoin, gp, ii, jj, aa
         integer(ip)                      :: pgaus, pnode, pelty
         real(rp)                         :: el_X(ndime,mnode), el_U(ndime,mnode), el_F(ndime,ndime,mnode)
         real(rp)                         :: gp_factor, gp_W(mgaus), gp_N(mnode,mgaus)
         real(rp)                         :: gp_dNde(ndime,mnode,mgaus), gp_dNdx(ndime,mnode,mgaus)
         real(rp)                         :: gp_inv_J(ndime,ndime,mgaus), gp_det_J(mgaus), gp_F(ndime,ndime,mgaus)
         ! -------------------------------
         if( INOTEMPTY ) then
            ! Initialise GDEPO
            call array_operations_initialization(self % nod_F) 
         
            ! Assemble the Deformation gradient
            if( associated(displ) )then
               do ielem = 1, nelem
                  pelty = ltype(ielem)
                  if( pelty > 0_ip )then
                      pnode = nnode(pelty)
                      pgaus = ngaus(pelty)
            
                      ! Gather coordinates
                      do aa = 1, pnode
                         ipoin = lnods(aa,ielem)
                         do ii = 1, ndime
                              el_X(ii,aa) = coord(ii,ipoin)
                              el_U(ii,aa) = displ(ii,ipoin,1)
                         end do
                      end do
            
                      ! Get interpolation information
                      gp_W(1:pgaus) = elmar(pelty) % weigp(1:pgaus)
                      gp_N(1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
                      gp_dNde(1:ndime,1:pnode,1:pgaus) = elmar(pelty) % deriv(1:ndime,1:pnode,1:pgaus)
            
                      ! Geometric variables at GP
                      call elmgeo_cartesian_derivatives_jacobian(ndime,pnode,pnode,pgaus,el_X,gp_dNde,gp_inv_J,gp_dNdx,gp_det_J)
                      
                      ! Deformation gradient at GP
                      gp_F(:,:,:) = 0.0_rp
                      do gp = 1, pgaus
                         do ii = 1, ndime
                            do jj = 1, ndime
                               do aa = 1, pnode
                                  gp_F(ii,jj,gp) = gp_F(ii,jj,gp) + el_U(ii,aa) * gp_dNdx(jj,aa,gp)
                               enddo
                            enddo
                         enddo
                         do ii = 1, ndime
                            gp_F(ii,ii,gp) = gp_F(ii,ii,gp) + 1.0_rp
                         enddo
                      enddo
                      
                      ! Deformation gradient at NODES
                      el_F(:,:,:) = 0.0_rp
                      do gp = 1, pgaus
                          gp_factor = gp_det_J(gp) * gp_W(gp)
                          do aa = 1,pnode
                             ipoin = lnods(aa,ielem)
                             do ii = 1,ndime
                                do jj = 1,ndime
                                   el_F(ii,jj,aa) = el_F(ii,jj,aa) + gp_factor * gp_N(aa,gp) * gp_F(ii,jj,gp)
                                end do
                             end do
                          end do
                      enddo
                     
                      ! Assemble
                      call matrix_assemble_element_RHS(ndime*ndime,ndime*ndime,pnode,lnods(:,ielem),el_F,self % nod_F)
            
                   endif
            
               enddo
            
               ! Interchange with other workers
               call rhsmod(ndime*ndime,self % nod_F)
            
               ! Average
               do ipoin=1, npoin
                  self % nod_F(1:ndime,1:ndime,ipoin) = self % nod_F(1:ndime,1:ndime,ipoin) / vmass(ipoin)
               end do
            
            else
         
              ! When not DISPL is associated F = I
               do ipoin=1, npoin
                  do ii =1, ndime
                     self % nod_F(ii,ii,ipoin) = 1.0_rp 
                  end do
               end do
         
            end if
         end if
         ! -------------------------------
      end subroutine biofib_compute_deformation_gradient_at_nodes

      subroutine biofib_compute_deformation_gradient_at_GP( self, ielem, gp_F )

         use def_master,  only             : displ
         use def_domain,  only             : ltype, lnods, elmar, coord
         use def_domain,  only             : mgaus, mnode, nnode, ngaus, ndime
         use mod_elmgeo,  only             : elmgeo_cartesian_derivatives_jacobian
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip),       intent(in)    :: ielem
         real(rp),          intent(out)   :: gp_F(:,:,:)
         integer(ip)                      :: ipoin, gp, ii, jj, aa
         integer(ip)                      :: pgaus, pnode, pelty
         real(rp)                         :: el_X(ndime,mnode), el_U(ndime,mnode)
         real(rp)                         :: gp_W(mgaus), gp_N(mnode,mgaus)
         real(rp)                         :: gp_dNde(ndime,mnode,mgaus), gp_dNdx(ndime,mnode,mgaus)
         real(rp)                         :: gp_inv_J(ndime,ndime,mgaus), gp_det_J(mgaus)
         ! -------------------------------
         ! Get some dimensions
         pelty = ltype(ielem)
         pnode = nnode(pelty)
         pgaus = ngaus(pelty)
 
         ! Gather coordinates
         do aa = 1, pnode
            ipoin = lnods(aa,ielem)
            do ii = 1, ndime
               el_X(ii,aa) = coord(ii,ipoin)
               el_U(ii,aa) = displ(ii,ipoin,1)
            end do
         end do

         ! Get interpolation information
         gp_W(1:pgaus) = elmar(pelty) % weigp(1:pgaus)
         gp_N(1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
         gp_dNde(1:ndime,1:pnode,1:pgaus) = elmar(pelty) % deriv(1:ndime,1:pnode,1:pgaus)

         ! Geometric variables at GP
         call elmgeo_cartesian_derivatives_jacobian(&
             ndime,pnode,pnode,pgaus,el_X,gp_dNde,gp_inv_J,gp_dNdx,gp_det_J)

         ! Deformation gradient at GP
         gp_F(:,:,:) = 0.0_rp
         do gp = 1, pgaus
            do ii = 1, ndime
               do jj = 1, ndime
                  do aa = 1, pnode
                     gp_F(ii,jj,gp) = gp_F(ii,jj,gp) + el_U(ii,aa) * gp_dNdx(jj,aa,gp)
                  enddo
               enddo
            enddo
            do ii = 1, ndime
               gp_F(ii,ii,gp) = gp_F(ii,ii,gp) + 1.0_rp
            enddo
         enddo
         ! -------------------------------
      end subroutine biofib_compute_deformation_gradient_at_GP


      subroutine biofib_point_nodal_fibers_wtime( vecptr, wfiber, wtime )
         ! ------------------------------------------
         real(rp),         intent(inout), pointer  :: vecptr(:,:)
         character(len=*), intent(in)              :: wfiber
         character(len=*), intent(in),    optional :: wtime
         integer(ip)                               :: itime
         ! ------------------------------------------
         ! Reference or current time?
         if( present(wtime) )then
            if( wtime == 'REFERENCE' )then
               itime = 3_ip
            elseif( wtime == 'PAST' )then
               itime = 2_ip
            elseif( wtime == 'CURRENT' )then
               itime = 1_ip
            else
               call runend('biofib_point_nodal_fibers: wtime option not valid')
            endif
         else
            ! Assuming current
            itime = 1_ip
         endif
         ! Update fibers if is necessary
         !if( itime == 2_ip )then
         !   call biofibers % get_reference_fibers_at_nodes()
         !else
         !   call biofibers % get_reference_fibers_at_nodes()
         !   call biofibers % get_current_fibers_at_nodes()
         !endif
         ! Pointing the fiber LONGITUDINAl/SHEET/NORMAL
         if(     wfiber == 'LONGI' .or. wfiber == 'LONGITUDINAL' .or. wfiber == 'FIBER' )then
            vecptr => biofibers % nod_lng(:,:,itime)
         elseif( wfiber == 'SHEET' )then
            vecptr => biofibers % nod_sht(:,:,itime)
         elseif( wfiber == 'NORMA' .or.  wfiber == 'NORMAL' )then
            vecptr => biofibers % nod_nrm(:,:,itime)
         else
            call runend('biofib_point_nodal_fibers: wfiber option not valid')
         endif
         ! -------------------------------
      end subroutine biofib_point_nodal_fibers_wtime


      subroutine biofib_update_fibers_at_nodes( self, itask )
         use def_master,             only  : INOTEMPTY
         use def_master,             only  : ITASK_INIUNK, ITASK_BEGSTE, ITASK_ENDSTE
         use def_domain,             only  : npoin
         ! ---------------------------------
         class(bio_fibers), intent(inout) :: self
         integer(ip),      intent(in)     :: itask
         ! ---------------------------------
         if( INOTEMPTY )then

            select case(itask)

               case(ITASK_INIUNK)
                  ! Compute (:,:,2)
                   call biofib_compute_reference_fibers_at_nodes(self)
                  ! (:,:,2) <= (:,:,3) 
                  self % nod_lng(1:3,1:npoin,2) = self % nod_lng(1:3,1:npoin,3)
                  self % nod_sht(1:3,1:npoin,2) = self % nod_sht(1:3,1:npoin,3)
                  self % nod_nrm(1:3,1:npoin,2) = self % nod_nrm(1:3,1:npoin,3)
                  ! (:,:,1) <= (:,:,2) 
                  self % nod_lng(1:3,1:npoin,1) = self % nod_lng(1:3,1:npoin,2)
                  self % nod_sht(1:3,1:npoin,1) = self % nod_sht(1:3,1:npoin,2)
                  self % nod_nrm(1:3,1:npoin,1) = self % nod_nrm(1:3,1:npoin,2)

               case(ITASK_BEGSTE)
                  ! (:,:,1) <= (:,:,2) 
                  self % nod_lng(1:3,1:npoin,2) = self % nod_lng(1:3,1:npoin,1)
                  self % nod_sht(1:3,1:npoin,2) = self % nod_sht(1:3,1:npoin,1)
                  self % nod_nrm(1:3,1:npoin,2) = self % nod_nrm(1:3,1:npoin,1)

               case(ITASK_ENDSTE)
                  ! Compute (:,:,1)
                  call biofib_compute_current_fibers_at_nodes(self)

            end select

         endif
         ! --------------------------------
      end subroutine biofib_update_fibers_at_nodes


      subroutine biofib_point_nodal_fibers_all( vecptr, wfiber )
         ! ------------------------------------------
         real(rp),         intent(inout), pointer  :: vecptr(:,:,:)
         character(len=*), intent(in)              :: wfiber
         ! ------------------------------------------
         if(     wfiber == 'LONGI' .or. wfiber == 'LONGITUDINAL' .or. wfiber == 'FIBER' )then
            vecptr => biofibers % nod_lng
         elseif( wfiber == 'SHEET' )then
            vecptr => biofibers % nod_sht
         elseif( wfiber == 'NORMA' .or.  wfiber == 'NORMAL' )then
            vecptr => biofibers % nod_nrm
         else
            call runend('biofib_point_nodal_fibers: wfiber option not valid')
         endif
         ! -------------------------------
      end subroutine biofib_point_nodal_fibers_all

 
      subroutine biofib_display_warning( )
         use mod_messages,          only :  messages_live
         ! -------------------------------
         call messages_live("-------------------------------------------------------","WARNING")
         call messages_live(" A deprecated syntax is used to define the FIBER_MODEL:",'WARNING')
         call messages_live("    FIBER_MODEL FIELD=1                                ",'WARNING')
         call messages_live("    ORTHOTROPIC_MODEL FIELD=1,2                        ",'WARNING')
         call messages_live("    FIBER_MODEL                                        ",'WARNING')
         call messages_live("       ...                                             ",'WARNING')
         call messages_live("    END_FIBER_MODEL                                    ",'WARNING')
         call messages_live(" Instead of the old syntax, use the in KER file:       ","WARNING")
         call messages_live("    PHYSICAL_PROBLEM                                   ","WARNING")
         call messages_live("       BIO_FIBER_MODEL                                 ","WARNING")
         call messages_live("         LOCATION: NODAL | ELEMENTAL                   ","WARNING")
         call messages_live("         FIBER:    1     $ (mandatory ID_FIELD)        ","WARNING")
         call messages_live("         SHEET:    2     $ (optional ID_FIELD)         ","WARNING")
         call messages_live("         NORMA:    3     $ (optional ID_FIELD)         ","WARNING")
         call messages_live("       END_BIO_FIBER_MODEL                             ","WARNING")
         call messages_live("    END_PHYSICAL_PROBLEM                               ","WARNING")
         call messages_live("-------------------------------------------------------","WARNING")
         call runend('MOD_BIOFIBERS: (L495)')
         ! -------------------------------
      end subroutine biofib_display_warning


      ! ----------------------------------------
      !
      ! Private methods
      !
      ! ----------------------------------------

      subroutine init_biofibers( self )
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
         if( self % initialised ) return
         self % initialised = .true.
         self % location    = NOT_DEFINED
         self % field_lng_X = 0_ip
         self % field_sht_X = 0_ip
         self % field_nrm_X = 0_ip
         self % vector_lng  = 0.0_rp
         self % vector_sht  = 0.0_rp
         self % vector_nrm  = 0.0_rp
         self % basis       = 0.0_rp
         self % lng_X => null()
         self % sht_X => null()
         self % nrm_X => null()
         self % nod_lng => null()
         self % nod_sht => null()
         self % nod_nrm => null() 
         ! -------------------------------
      end subroutine init_biofibers


      subroutine delete_biofibers( self )
         ! -------------------------------
         class(bio_fibers), intent(inout) :: self
         ! -------------------------------
         self % initialised = .false.
         self % location    = NOT_DEFINED
         self % field_lng_X = 0_ip
         self % field_sht_X = 0_ip
         self % field_nrm_X = 0_ip
         call self % deallocate_memory()
         self % lng_X => null()
         self % sht_X => null()
         self % nrm_X => null()
         self % nod_lng => null() 
         self % nod_sht => null() 
         self % nod_nrm => null() 
         ! -------------------------------
      end subroutine delete_biofibers

      subroutine check_consistency( icase )
         use def_master, only            : coupling
         use def_domain, only            : ndime
         use mod_messages, only          :  messages_live
         ! -------------------------------
         integer(ip), intent(in)        :: icase
         ! -------------------------------
         select case(icase)

            case(0_ip)

               if( ndime /= 3_ip )then
                  call messages_live("-------------------------------------------------------","WARNING")
                  call messages_live(" BIO-FIBERS NOT ENABLED FOR 2D PROBLEMS                ",'WARNING')
                  call messages_live("-------------------------------------------------------","WARNING")
                  call runend('BIOFIB: check_consistency(0)')
               endif

            case(1_ip)

               if( coupling('SOLIDZ','EXMEDI') >= 1_ip .or. coupling('EXMEDI','SOLIDZ') >= 1_ip )then
                  call messages_live("-------------------------------------------------------","WARNING")
                  call messages_live(" ELEMENTAL FIBERS NOT ENABLED WHEN ECCOUPLING          ",'WARNING')
                  call messages_live("-------------------------------------------------------","WARNING")
                  call runend('BIOFIB_READ_DATA: line 100')
               endif

         end select
         ! -------------------------------
      end subroutine check_consistency

      subroutine init_array( array )
         ! -------------------------------
         real(rp),   intent(inout) :: array(:,:)
         ! -------------------------------
         integer(ip)               :: ii, jj
         ! -------------------------------
         do jj = 1, size(array,dim=2)
            do ii = 1, size(array,dim=1)
               array(ii,jj) = 0.0_rp
            enddo
         enddo
         ! -------------------------------
      end subroutine init_array

end module mod_biofibers
