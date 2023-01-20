!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    bounor.f90
!> @author  Guillaume Houzeaux
!> @date    18/10/2012
!> @brief   Computes the boundary normals
!> @details Comppute the boundary normals. Only valid for 
!!          TRI03, TRI06, QUA04, QUA08, QUA09 boundaries.
!!          For QUA04 and QUA09, thye exterior normal is averaged
!!          over the two triangles forming the boundary element. 
!!          According to the mesher used, we have two options:\n
!!          - Clock wise numbering: normal points inwards
!!          - Counterclock wise numbering: normal points outwards
!!            - GiD: normal points inwards
!!            - Fensap: normal points outwards
!!            - Windmesh: normal points outwards
!!          \verbatim
!!          - BOUNO(1:NDIME,IBOUN) .... Exterior normal to IBOUN
!!          - NINVE ................... Number of inverted normals
!!                                      that were pointing inwards
!!          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine bounor(kboun,lnodb,ltypb,lelbo,ninve,bouno)
  use def_kintyp, only     :  ip,rp
  use def_elmtyp, only     :  TRI03,QUA04,TRI06,QUA09
  use def_domain, only     :  ndime,mnodb,mnode,coord
  use def_domain, only     :  nnode,ltype,lnods
  implicit none
  integer(ip), intent(in)  :: kboun                          !< Number of boundaries
  integer(ip), intent(in)  :: lnodb(mnodb,kboun)             !< Boundary connectivity
  integer(ip), intent(in)  :: ltypb(kboun)                   !< Boundary type
  integer(ip), intent(in)  :: lelbo(kboun)                   !< Boundary/element connectvity
  integer(ip), intent(out) :: ninve                          !< Number of inverted boundaries
  real(rp),    intent(out) :: bouno(ndime,kboun)             !< Boundary exterior normals
  integer(ip)              :: iboun,p1,p2,p3,idime
  integer(ip)              :: ipoin,inode,inodb,ielem
  integer(ip)              :: pnodb,pblty,pnode
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: bocod(ndime,mnodb)
  real(rp)                 :: vec(3,3),rmod

  if( ndime == 2 ) then

     !-------------------------------------------------------------------
     !
     ! 2D case
     !
     !-------------------------------------------------------------------

     do iboun = 1,kboun
        p1             =  lnodb(1,iboun)
        p2             =  lnodb(2,iboun)
        vec(1,1)       =  coord(1,p2) - coord(1,p1)
        vec(2,1)       =  coord(2,p2) - coord(2,p1)
        bouno(1,iboun) =  vec(2,1) 
        bouno(2,iboun) = -vec(1,1) 
     end do

  else if( ndime == 3 ) then

     !-------------------------------------------------------------------
     !
     ! 3D case
     !
     !-------------------------------------------------------------------

     do iboun = 1,kboun
        pblty = abs(ltypb(iboun))

        if( pblty == TRI03 .or.  pblty == TRI06 ) then

           p1                 = lnodb(1,iboun)
           p2                 = lnodb(2,iboun)
           p3                 = lnodb(3,iboun)
           call nortri(p1,p2,p3,coord,vec,ndime)
           bouno(    1,iboun) = 0.5_rp*vec(1,3)
           bouno(    2,iboun) = 0.5_rp*vec(2,3)
           bouno(ndime,iboun) = 0.5_rp*vec(3,3)

        else if( pblty == QUA04 .or. pblty == QUA09 ) then

           p1             = lnodb(1,iboun)
           p2             = lnodb(2,iboun)
           p3             = lnodb(3,iboun)
           call nortri(p1,p2,p3,coord,vec,ndime)
           bouno(1,iboun) = vec(1,3)
           bouno(2,iboun) = vec(2,3)
           bouno(3,iboun) = vec(3,3)

           p1             = lnodb(1,iboun)
           p2             = lnodb(3,iboun)
           p3             = lnodb(4,iboun)
           call nortri(p1,p2,p3,coord,vec,ndime)
           bouno(1,iboun) = bouno(1,iboun) + vec(1,3)
           bouno(2,iboun) = bouno(2,iboun) + vec(2,3)
           bouno(3,iboun) = bouno(3,iboun) + vec(3,3)

           p1             = lnodb(1,iboun)
           p2             = lnodb(2,iboun)
           p3             = lnodb(4,iboun)
           call nortri(p1,p2,p3,coord,vec,ndime)
           bouno(1,iboun) = bouno(1,iboun) + vec(1,3)
           bouno(2,iboun) = bouno(2,iboun) + vec(2,3)
           bouno(3,iboun) = bouno(3,iboun) + vec(3,3)
           
           p1             = lnodb(4,iboun)
           p2             = lnodb(2,iboun)
           p3             = lnodb(3,iboun)
           call nortri(p1,p2,p3,coord,vec,ndime)
           bouno(1,iboun) = bouno(1,iboun) + vec(1,3)
           bouno(2,iboun) = bouno(2,iboun) + vec(2,3)
           bouno(3,iboun) = bouno(3,iboun) + vec(3,3)          
           !
           ! HEX: We assume that ordering points towards inside
           !
           rmod = sqrt(   bouno(1,iboun) * bouno(1,iboun) &
                &      + bouno(2,iboun) * bouno(2,iboun) &
                &      + bouno(3,iboun) * bouno(3,iboun) )
           bouno(1,iboun) = bouno(1,iboun) / rmod
           bouno(2,iboun) = bouno(2,iboun) / rmod
           bouno(3,iboun) = bouno(3,iboun) / rmod

        end if

     end do
  end if
  !
  ! Ensure normal points outwards and normalize it
  !
  do iboun = 1,kboun
     pblty = abs(ltypb(iboun))
     pnodb = nnode(pblty)
     ielem = lelbo(iboun)
     pnode = nnode(ltype(ielem))
     do inodb = 1,pnodb
        ipoin = lnodb(inodb,iboun)
        do idime = 1,ndime
           bocod(idime,inodb) = coord(idime,ipoin)
        end do
     end do
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        do idime = 1,ndime
           elcod(idime,inode) = coord(idime,ipoin)
        end do
     end do
     call chebou(pnode,pnodb,bouno(1,iboun),bocod,elcod,ninve)
     call vecuni(ndime,bouno(1,iboun),rmod)
  end do

end subroutine bounor
      
