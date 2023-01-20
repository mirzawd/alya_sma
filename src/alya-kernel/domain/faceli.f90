!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine faceli(&
     sabox,blink,ltypb,lnodb,coord,icoor,jcoor,icoo2,mnodb,nboun,proje,iboun,dista)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_faliib
  ! DESCRIPTION
  !    This routines find the intersection between a line segment and a particle.
  !    Also, returns the normalized distance measure from the first point 
  !    INPUT 
  !       poin1, poin2: points of the line segment 
  !       imbou: particle identification
  !       nboun: number of ibouns of the particle
  !    OUTPUT
  !       inter: intersection point
  !       iboun: intersection iboun
  !       dista: normalized distance measure from the first point. 
  !              = 0 if the intersection point is the first point
  !              > 0 and < 1 if the intersection point is inside the segment
  !              = 1 if the intersection point is the second point
  ! USED BY
  !    fielib 
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,nnode
  implicit none  
  
  integer(ip), intent(in)  :: mnodb    
  integer(ip), intent(in)  :: nboun
  real(rp),    intent(in)  :: sabox(2,ndime,2*nboun-1)
  integer(ip), intent(in)  :: blink(2*nboun-1)
  integer(ip), intent(in)  :: ltypb(nboun)
  integer(ip), intent(in)  :: lnodb(mnodb,nboun)
  real(rp),    intent(in)  :: coord(ndime,*)
  real(rp),    intent(in)  :: icoor(ndime),jcoor(ndime)  
  real(rp),    intent(in)  :: icoo2(ndime)
  real(rp),    intent(out) :: proje(ndime),dista
  integer(ip), intent(out) :: iboun

  integer(ip), pointer     :: canfa(:) => null()
  integer(ip)              :: nlist,inodb,ipoin
  integer(ip)              :: idime,ilist,pblty,pnodb,ifoun
  real(rp)                 :: bobox(3,2),bouno(3),numer,denom,toler
  real(rp)                 :: bocod(ndime,mnodb),dummr
  real(rp)                 :: plapo(3),temdi,dist2,dumma(ndime)

  allocate( canfa(nboun))
  toler =  1.0e-1_rp
  iboun = -100_ip
  dista =  1.0e9_rp
  bobox(3,1) = 0.0_rp
  bobox(3,2) = 0.0_rp

  !
  ! Maybe the segment doesn't intersect the particle (i.e.: Lohner method)
  !
  do idime = 1,ndime
     bobox(idime,1) = min(icoor(idime),jcoor(idime))  
     bobox(idime,2) = max(icoor(idime),jcoor(idime)) 
     temdi          = bobox(idime,2) - bobox(idime,1)
     bobox(idime,1) = bobox(idime,1) - abs(temdi*toler) 
     bobox(idime,2) = bobox(idime,2) + abs(temdi*toler)
  end do
  call facbox(sabox,blink,bobox,nboun,canfa,nlist)
  ilist = 0

  do while ( ilist < nlist ) 
     ilist = ilist + 1
     !
     ! Element properties and dimensions
     !
     iboun = canfa(ilist)
     pblty = ltypb(iboun)
     pnodb = nnode(pblty)
     !
     ! bouno: Exterior normal           
     !
     call extbou(1_ip,pnodb,lnodb(1,iboun),coord,bouno)
     !
     ! inode: Number of the first node in the actual particle iboun 
     ! 
     do inodb = 1,pnodb
        ipoin = lnodb(inodb,iboun)                  
        do idime = 1,ndime
           bocod(idime,inodb) = coord(idime,ipoin)
        end do
     end do
   
     numer = 0.0_rp
     denom = 0.0_rp
     do idime = 1,ndime
        numer = numer + ( bouno(idime) * ( bocod(idime,1) - icoor(idime) )  )
        denom = denom + ( bouno(idime) * ( jcoor(idime)   - icoor(idime) )  )
     end do
       
     if ( denom /= 0.0_rp ) then       
        !
        ! temdi: Normalized Distance measure from the first node
        !
        temdi = numer / denom 
        !
        ! The segment intersects the plane that contains the iboun
        !
        if ( temdi >= -0.01_rp .and. temdi <= 1.01_rp ) then           
           do idime = 1,ndime
              plapo(idime) = icoor(idime) + temdi * ( jcoor(idime) - icoor(idime) )
           end do
           !
           ! Determine if the projection point on plane is inside the iboun
           !
           if( ndime == 3 ) then
              call instr2(pnodb,plapo,bocod,ifoun,dummr,dummr)
           else if( ndime == 2 ) then
              call dposeg(plapo,bocod(1,1),bocod(1,2),dummr,dumma,ifoun)
           end if

           dist2 = 0.0_rp
           do idime = 1,ndime
              dist2 = dist2 +  (plapo(idime) - icoo2(idime))*(plapo(idime) - icoo2(idime)) 
           end do
           !
           ! The projected point is in the iboun
           !
           if( ifoun == 1 .and. abs(dist2) < abs(dista)) then
              do idime = 1,ndime
                 proje(idime) = plapo(idime) 
              end do
              dista = dist2
              iboun = canfa(ilist)
           end if
        end if
     else
        ifoun = 0
        call runend('FACELI: ZERO DENOM')
     end if
  end do

  deallocate( canfa )

end subroutine faceli
