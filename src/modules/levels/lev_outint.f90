!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_outint
  !****f* Levels/lev_outint
  ! NAME 
  !    lev_outint
  ! DESCRIPTION
  !    Postprocess Interface at it = ittim
  ! USES
  !    
  ! USED BY
  !    lev_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  implicit none
  integer(ip)             :: icomp,idime,jdime,ninterf,inter,j,p1,p2,p3
  integer(ip)             :: compp(ndime)
  real(rp)                :: xx(ndime,3),xinter,yinter,xmin
  real(rp)                :: ni,nj,nk,vec(3,3),nn

  if( kfl_alloc_lev == 0 ) call lev_conint()

  if( INOTSLAVE ) then

     if( ndime == 3 .and. 2 == 1 ) then
        write(lun_inter_lev,1)
     else
        write(lun_inter_lev,'(a)') 'solid levels'
     end if

     if(ndime==2) then

        ninterf = nboun_lev*2_ip

        do icomp=1,ninterf 

           xmin = 1000000000.0_rp

           do j=icomp,ninterf
              if(coord_lev(1,j)<xmin) then
                 xmin = coord_lev(1,j)
                 inter = j
              endif   
           enddo 

           xinter = coord_lev(1,icomp)
           yinter = coord_lev(2,icomp)
           coord_lev(1,icomp) = coord_lev(1,inter)
           coord_lev(2,icomp) = coord_lev(2,inter)
           coord_lev(1,inter) = xinter         
           coord_lev(2,inter) = yinter         

        enddo   

     endif   

     do icomp=1,nboun_lev

        do idime=1,ndime
           compp(idime)=lnodb_lev(idime,icomp)
        end do

        do idime = 1,ndime
           do jdime = 1,ndime
              xx(idime,jdime) = coord_lev(idime,compp(jdime))
           end do
        end do

!!$           write(lun_inter_lev,2) x1, y1, z1, x2, y2, z2, x3, y3, z3
        if( ndime == 3 .and. 2 == 1 ) then
           write(lun_inter_lev,7) ( ( xx(idime,jdime), idime = 1,ndime ) , jdime = 1,ndime ) 
           write(lun_inter_lev,4)
        else if( ndime == 3 .and. 1 == 1 ) then
           p1 = lnodb_lev(1,icomp)
           p2 = lnodb_lev(2,icomp)
           p3 = lnodb_lev(3,icomp)
           call nortri(p1,p2,p3,coord_lev,vec,ndime)
           ni = vec(1,3)
           nj = vec(2,3)
           nk = vec(3,3)
           nn = 1.0_rp / sqrt( ni*ni + nj*nj + nk*nk )
           ni = ni * nn
           nj = nj * nn
           nk = nk * nn
           write(lun_inter_lev,'(a,3(1x,e12.6))') 'facet normal ',ni,nj,nk
           write(lun_inter_lev,'(a)')             'outer loop'
           write(lun_inter_lev,'(a,3(1x,e12.6))') 'vertex ',coord_lev(1,p1),coord_lev(2,p1),coord_lev(3,p1)
           write(lun_inter_lev,'(a,3(1x,e12.6))') 'vertex ',coord_lev(1,p2),coord_lev(2,p2),coord_lev(3,p2)
           write(lun_inter_lev,'(a,3(1x,e12.6))') 'vertex ',coord_lev(1,p3),coord_lev(2,p3),coord_lev(3,p3)
           write(lun_inter_lev,'(a)')             'endloop'
           write(lun_inter_lev,'(a)')             'endfacet'
        else
           write(lun_inter_lev,*) xx(1,1),xx(2,1)
           write(lun_inter_lev,*) xx(1,2),xx(2,2)
        end if
!!$           write(lun_inter_lev,3) x1, y1, z1, normt(1), normt(2), normt(3), x2, y2, z2, normt(1), &
!!$           normt(2), normt(3), x3, y3, z3 normt(1), normt(2), normt(3)
!!$           write(lun_inter_lev,5) x1, y1, z1, x2, y2, z2, x3, y3, z3

     enddo
     if( ndime == 3 .and. 1 == 1 ) then
        write(lun_inter_lev,'(a)') 'endsolid levels'
     end if
  endif

  !! POVRay Format 
1 format(' union { ')
  !! Triangle
2 format(' triangle { <',F10.8,' ,',F10.8,' ,',F10.8,'>, <',F10.8,' ,',F10.8,&
       ' ,',F10.8,'>, <',F10.8,' ,',F10.8,' ,',F10.8,'> }')
  !! Smooth Triangle
3 format(' smooth_triangle { <',F10.8,' ,',F10.8,' ,',F10.8,'>, <',F10.8,' ,',&
       F10.8,' ,',F10.8,'>, <',F10.8,' ,',F10.8,' ,',F10.8,'>, <',F10.8,' ,',&
       F10.8,' ,',F10.8,'>, <',F10.8,' ,',F10.8,' ,',F10.8,'>, <',F10.8,' ,',&
       F10.8,' ,',F10.8,'> }')
4 format(' } ')
  !! Simple Formate to write description of interface triangles
5 format(F10.8,' ',F10.8,' ',F10.8,' ',F10.8,' ',F10.8,' ',F10.8,' ',F10.8,' ',&
       F10.8,' ',F10.8)
6 format(' triangle { <',E24.16,' ,',E24.16,' ,',E24.16,'>, <',E24.16,' ,',&
       E24.16,' ,',E24.16,'>, <',E24.16,' ,',E24.16,' ,',E24.16,'> }')
7 format(' triangle { <',F12.8,' ,',F12.8,' ,',F12.8,'>, <',F12.8,' ,',F12.8,' ,',&
       F12.8,'>, <',F12.8,' ,',F12.8,' ,',F12.8,'> }')
  !! Smooth Triangle

end subroutine lev_outint
