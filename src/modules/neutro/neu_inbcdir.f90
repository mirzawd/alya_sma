!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 !----------------------------------------------------------------------- 
!> @addtogroup NeutroTurnon 
!> @ingroup Neutro !> @{ 
!> @file neu_inbcDir.f90 
!> @author A.S 
!> @brief Impose boundary conditions 
!> @details Impose boundary conditions 
!> @} !----------------------------------------------------------------------- 
subroutine neu_inbcdir()
  use def_parame 
  use def_inpout 
  use def_master 
  use def_kermod 
  use def_domain 
  use def_neutro 
  use mod_bouder, only : bouder
  implicit none 
  integer(ip) :: ienergy,idirection,iboun,pnodb,inode,ipoin,pblty,ielem,pelty,pgaub,igaus,pnode,pgaus,igaub!,iunkn,ii 
  real(rp)    :: xjaci(9),xjacm(9),DOTPR_RETORNA,eucta,gpdet
  real(rp)    :: gpcar(ndime,mnode,mgaus) 
  real(rp)    :: baloc(ndime,ndime,mgaub) !> Local directions, normal to boundary element (baloc(1:ndime, ndime)) and 
                                          !> two tangencial (baloc(1:ndime,1) 
!  and baloc(1:ndime,2) for 3D 
  real(rp)    :: elcod(ndime,mnode) !> Element code
  real(rp)    :: bocod(ndime,mnodb) !> Boundary code
  external :: elmder, chenor


  if( INOTMASTER ) then

    do idirection = 1,num_directions_neu 
      do ienergy = 1,num_energies_neu 
        do iboun = 1,nboun 
          if(kfl_fixbo_neu(ienergy,idirection) % l(iboun) == 1 ) then ! es vacuum
                 
            pblty = ltypb(iboun) !> We put the boundary type of the current boundary element in pblty 
            pnodb = nnode(pblty) !> From the boundary type we get the number of nodes in the boundary 
            ielem = lelbo(iboun) !> We determine the identifier for the current element touching the boundary

            pelty = ltype(ielem) 
            pnode = nnode(pelty) !> Number of nodes in current element 
            pgaus = ngaus(pelty) !> Number of Gauss points in current element 
            pgaub =     ngaus(pblty) !> Number of Gauss points in the boundary
            ! Gather operations ! 
            bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun)) !> Boundary coord determination 
            elcod(1:ndime,1:pnode) =    coord(1:ndime,lnods(1:pnode,ielem)) !> Element coord determination 
            do igaus = 1,pgaus
              call elmder(pnode,ndime,elmar(pelty) % deriv(1:ndime,1:pnode,igaus),& ! Cartesian derivative 
                          elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci) ! and Jacobian (elmar is     the element data base matrix)
            end do 
            do igaub = 1,pgaub !> This routine calculates the baloc (normal and tangencial directions in the element) and eucta 
              call bouder(pnodb,ndime,ndimb,elmar(pblty) % deriv(1:ndimb,1:pnodb,igaub), bocod,baloc(:,:,igaub),eucta) 
              call chenor(pnode,baloc(:,:,igaub),bocod,elcod) ! Check normal
            end do 
            dotpr_retorna = dot_product(direc_neu(1:ndime,idirection),baloc(1:ndime,ndime,1)) 

            if(dotpr_retorna<=0.0_rp .and. kfl_fixbo_neu(ienergy,idirection) % l(iboun)   == 1) then
                    
              do inode=1,pnodb 
                  ipoin=lnodb(inode,iboun) 
                  kfl_fixno_neu(ienergy,idirection) % l(1,ipoin)=1 

              !  write(1356,*) idirection,ipoin,dotpr_retorna 
              enddo 
            endif 
          endif 
        enddo
      enddo 
    enddo 
  end if 

end subroutine neu_inbcdir
