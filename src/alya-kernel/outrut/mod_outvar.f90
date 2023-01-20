!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup IO
!> @{
!> @file    mod_outvar.f90
!> @author  houzeaux
!> @date    2020-08-27
!> @brief   Output variables 
!> @details Output variables for time step jttim and time value dutim
!-----------------------------------------------------------------------

module mod_outvar

  use def_kintyp_basic,      only : ip,rp,lg,r3p
  use def_domain,            only : ndime
  use def_domain,            only : ntens
  use def_domain,            only : npoin
  use def_domain,            only : nelem
  use def_domain,            only : nboun
  use def_domain,            only : memor_dom
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_size
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_parame,            only : zero,one
  use def_master,            only : ivapo,ITASK_ENDRUN,INOTMASTER,iasca,iavec,iar3p
  use def_master,            only : IMASTER,mitim,ittyp,intost,kfl_paral
  use def_master,            only : gescx,gesca,gevec,ger3p,gevex,postp
  use def_kermod,            only : witness_mesh
  use mod_postpr,            only : postpr
  use def_postpr,            only : kfl_ivari
  use mod_communications,    only : PAR_MAX
  use mod_arrays,            only : arrays
  use mod_interp_fe,         only : interp_fe
  use mod_interp_fe,         only : interp_fe_deallocate
  implicit none

  private

  public :: outvar
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-08-27
  !> @brief   Output
  !> @details Main output subroutine
  !> 
  !-----------------------------------------------------------------------


  subroutine outvar(jvari,jttim,dutim,wopos,MESH,MESH_ID)

    integer(ip),                      intent(in) :: jvari
    integer(ip),                      intent(in) :: jttim
    real(rp),                         intent(in) :: dutim
    character(5),                     intent(in) :: wopos(*)
    type(mesh_type_basic),  optional, intent(in) :: MESH
    integer(ip),            optional, intent(in) :: MESH_ID    
    integer(ip)                                  :: itens,ipoin
    integer(ip)                                  :: idime,ivari,kdime,dim_max,imesh
    integer(ip)                                  :: icoun,jdime
    character(5)                                 :: wauxi(3,6)
    character(1)                                 :: waux1(9)
    character(2)                                 :: waux2(9)
    character(5)                                 :: waux3(3)
    character(20)                                :: wtype
    real(rp),               pointer              :: ggggg(:)
    real(rp),               pointer              :: vvvvv(:,:)
    real(rp),               pointer              :: gesca_loc(:)
    real(rp),               pointer              :: gevec_loc(:,:)
    logical(lg)                                  :: if_interpolated

    nullify(vvvvv)
    nullify(ggggg)
    nullify(gesca_loc)
    nullify(gevec_loc)

    ivari = abs(jvari)
    ivapo = ivari
    if( mitim == 0 .and. ittyp == ITASK_ENDRUN ) goto 10
    !
    ! A mesh is present
    !
    if(      present(MESH_ID) ) then
       imesh = MESH_ID
    else if( present(MESH)    ) then
       imesh = MESH % id
    else
       imesh = 0
    end if
    !
    ! Array should be interpolated 
    !
    if_interpolated = .false.
    if( imesh /= 0 ) then
       if( associated(witness_mesh(imesh) % inte % lelem) ) then
          if_interpolated = .true.
       end if
    end if

    if( kfl_ivari(1) > 0 ) then

       select case ( wopos(2) )

       case ( 'SCALA' ) 

          !-----------------------------------------------------------------
          !
          ! SCALA
          !
          !-----------------------------------------------------------------

          if( imesh == 0 ) then
             call postpr(gesca,wopos,jttim,dutim)
          else
             if( if_interpolated .and. INOTMASTER) then
                call interp_fe(gesca,witness_mesh(imesh) % inte,gesca_loc,witness_mesh(imesh) % mesh % npoin)
                call postpr(gesca_loc,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
             else
                call postpr(gesca,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
             end if
          end if

       case ( 'VECTO' , 'MATRI' ) 

          !-----------------------------------------------------------------
          !
          ! VECTOR and MATRIX
          !
          !-----------------------------------------------------------------

          if( wopos(2) == 'VECTO' ) then
             kdime = ndime
          else
             kdime = memory_size(gevec,1_ip)
             call PAR_MAX(kdime)
          end if

          if( imesh == 0 ) then
             call postpr(gevec,wopos,jttim,dutim,kdime)
          else
             if( if_interpolated .and. INOTMASTER ) then
                call interp_fe(gevec,witness_mesh(imesh) % inte,gevec_loc,witness_mesh(imesh) % mesh % npoin)
                call postpr   (gevec_loc,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
             else
                call postpr   (gevec,wopos,jttim,dutim,kdime,MESH=witness_mesh(imesh) % mesh)
             end if
          end if

       case ( 'R3P  ' , 'R3PVE' ) 

          !-----------------------------------------------------------------
          !
          ! R3P
          !
          !-----------------------------------------------------------------

          if( kfl_ivari(1) > 0 ) then 
             if( imesh == 0 ) then
                call postpr(ger3p,wopos,jttim,dutim)
             else
                call runend('OUTVAR: CANNOT POSTPROCESS R3P VARIABLES ON WTNESS MESHES')
             end if
          end if

       case ( 'MULTI' ) 

          !-----------------------------------------------------------------
          !
          ! Multi-dimensional
          !
          !-----------------------------------------------------------------

          waux3(1:3) = wopos(1:3)
          waux3(2)   = 'SCALA'
          kdime      = memory_size(gevec,1_ip)

          call PAR_MAX(kdime)

          select case ( wopos(3) )
          case ( 'NPOIN' ) ; dim_max = npoin
          case ( 'NELEM' ) ; dim_max = nelem
          case ( 'NBOUN' ) ; dim_max = nboun
          end select

          call memory_alloca(memor_dom,'GGGGG','outvar',ggggg,dim_max)

          do idime = 1,kdime
             wtype = intost(idime)
             if( idime < 10 ) then
                waux3(1) = wopos(1)(1:3)//'0'//trim(wtype(1:1))
             else if( idime < 100 ) then
                waux3(1) = wopos(1)(1:3)//trim(wtype(1:2))
             else if( idime < 1000 ) then               
                waux3(1) = wopos(1)(1:2)//trim(wtype(1:3))
             else              
                waux3(1) = trim(wtype(1:5))           
             end if 
             do ipoin = 1,dim_max
                ggggg(ipoin) = gevec(idime,ipoin)
             end do
             call postpr(ggggg,waux3,jttim,dutim) 
          end do
          call memory_deallo(memor_dom,'GGGGG','outvar',ggggg)

       case ( 'SCALX' )         

          !-----------------------------------------------------------------
          !
          ! SCALX
          !
          !-----------------------------------------------------------------

          call memory_alloca(memor_dom,'GGGGG','outvar',ggggg,npoin)

          wauxi(2,1) = 'SCALA'
          wauxi(3,1) = wopos(3)           

          wauxi(1,1) = wopos(1)(1:4)//'r'
          do ipoin = 1,npoin
             ggggg(ipoin) = real(gescx(ipoin))
          end do
          call postpr(ggggg,wauxi(:,1),jttim,dutim) 
           
          wauxi(1,1) = wopos(1)(1:4)//'i'           
          do ipoin = 1,npoin
             ggggg(ipoin) = aimag(gescx(ipoin))
          end do
          call postpr(ggggg,wauxi(:,1),jttim,dutim) 

          call memory_deallo(memor_dom,'GGGGG','outvar',ggggg)

       case ('VECTX' )        

          !-----------------------------------------------------------------
          !
          ! VECTX
          !
          !-----------------------------------------------------------------

          call memory_alloca(memor_dom,'VVVVV','outvar',vvvvv,ndime,npoin)

          wauxi(2,1) = 'VECTO'
          wauxi(3,1) = wopos(3)    

          wauxi(1,1) = wopos(1)(1:4)//'r'
          do ipoin = 1,npoin
             do idime = 1,ndime
                vvvvv(idime,ipoin) = real(gevex(idime,ipoin))
             end do
          end do
          call postpr(vvvvv,wauxi(:,1),jttim,dutim) 
          
          wauxi(1,1) = wopos(1)(1:4)//'i'  
          do ipoin = 1,npoin
             do idime = 1,ndime
                vvvvv(idime,ipoin) = aimag(gevex(idime,ipoin))
             end do
          end do
          call postpr(vvvvv,wauxi(:,1),jttim,dutim) 
          
          call memory_deallo(memor_dom,'VVVVV','outvar',vvvvv)

       case ( 'TENSO' ) 

          !-----------------------------------------------------------------
          !
          ! TENSO
          !
          !-----------------------------------------------------------------

          call memory_alloca(memor_dom,'GGGGG','outvar',ggggg,npoin)

          select case ( ndime )
          case ( 1_ip ) ; call runend('CANNOT POSTPROCESS A TENSOR IN 1D')
          case ( 2_ip ) ; waux2(1:6) = ['XX','YY','XY','TH','NU','ER']
          case ( 3_ip ) ; waux2(1:6) = ['XX','YY','ZZ','YZ','XZ','XY']
          end select

          do itens = 1,ntens
             wauxi(1,itens) = wopos(1)(1:3)//waux2(itens)
             wauxi(2,itens) = 'SCALA'
             wauxi(3,itens) = wopos(3)           
             do ipoin = 1,npoin
                ggggg(ipoin) = gevec(itens,ipoin)
             end do
             if( imesh == 0 ) then
                call postpr(ggggg,wauxi(:,itens),jttim,dutim)   
             else
                call postpr(ggggg,wauxi(:,itens),jttim,dutim,MESH=witness_mesh(imesh) % mesh)              
             end if
          end do

          call memory_deallo(memor_dom,'GGGGG','outvar',ggggg)

       case ( 'VECT2' ) 

          !-------------------------------------------------------------
          !
          ! VECT2: vector of vectors
          !
          !-------------------------------------------------------------

          select case ( wopos(3) )
          case ( 'NPOIN' ) ; dim_max = npoin
          case ( 'NELEM' ) ; dim_max = nelem
          case ( 'NBOUN' ) ; dim_max = nboun
          end select

          call memory_alloca(memor_dom,'VVVVV','outvar',vvvvv,ndime,dim_max)

          if ( ndime == 1 ) then

             call runend('CANNOT POSTPROCESS A TENSOR IN 1D')

          else

             waux1(1:3) = ['1','2','3']
             do idime = 1,ndime
                wauxi(1,idime) = wopos(1)(1:4)//waux1(idime)
                wauxi(2,idime) = 'VECTO'
                wauxi(3,idime) = wopos(3) 
                do ipoin = 1,dim_max
                   do jdime = 1,ndime
                      icoun              = (idime-1)*ndime+jdime
                      vvvvv(jdime,ipoin) = gevec(icoun,ipoin)
                   end do
                end do
                if( imesh == 0 ) then
                   call postpr(vvvvv,wauxi(:,idime),jttim,dutim)
                else
                   call postpr(vvvvv,wauxi(:,idime),jttim,dutim,MESH=witness_mesh(imesh) % mesh)
                end if

             end do

          end if
          
          call memory_deallo(memor_dom,'VVVVV','outvar',vvvvv)

       end select
    end if

10  continue

    if( iasca == 1 ) call memgen(2_ip,one,zero)
    if( iavec == 1 ) call memgen(2_ip,one,one)
    nullify(gesca)
    nullify(gevec)
    if( if_interpolated ) then
       call interp_fe_deallocate(gevec_loc)
       call interp_fe_deallocate(gesca_loc)
    end if
    ivapo = 0

  end subroutine outvar

end module mod_outvar
!> @}
