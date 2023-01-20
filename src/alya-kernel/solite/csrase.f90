!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine csrase(elmat,amatr,ndofn,pnode,mevat,lnode,iprob)
  !-----------------------------------------------------------------------
  !****f* solite/csrase
  ! NAME
  !    csrase
  ! DESCRIPTION
  !    This routine assemble the matrix using the bcsr format
  ! INPUT
  !    ELMAT
  !    NDOFN
  !    PNODE
  !    MEVAT
  !    LNODE
  !    IPROB
  ! OUTPUT
  !    AMATR
  ! USED BY
  !    ***_assmat
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,r_sol,c_sol,r_sym,c_sym
  implicit none 
  integer(ip), intent(in)    :: iprob
  integer(ip), intent(in)    :: ndofn,pnode,mevat
  integer(ip), intent(in)    :: lnode(pnode)
  real(rp),    intent(in)    :: elmat(mevat,mevat)
  real(rp),    intent(inout) :: amatr(ndofn*ndofn,*)
  integer(ip)                :: ievat,jevat,idofn,jdofn,idof2,idof3
  integer(ip)                :: inode,jnode,ndofi
  integer(ip)                :: ipoin,jpoin,izsol,jcolu,idime,jdime
  integer(ip)                :: iposi,jposi

  if(iprob==1) then
     !
     ! Old way of assembling NASTIN
     !     
     do inode=1,pnode
        ipoin = lnode(inode)
        do jnode=1,pnode
           jpoin = lnode(jnode)
           izsol = r_sol(ipoin)
           jcolu = c_sol(izsol)
           do while(jcolu/=jpoin)
              izsol = izsol + 1
              jcolu = c_sol(izsol)
           end do
           
           ievat=(inode-1)*ndime
           idof3=0

           do idime=1,ndime
              ievat=ievat+1
              jevat=(jnode-1)*ndime
              idof2=(idime-1)*ndofn
              do jdime=1,ndime
                 jevat=jevat+1
                 idof2=idof2+1
                 !$OMP ATOMIC
                 amatr(idof2,izsol) = amatr(idof2,izsol) + elmat(ievat,jevat)  ! Top left                    
              end do
              idof3=idof3+ndofn
              jevat=ndime*pnode+jnode
              !$OMP ATOMIC
              amatr(idof3,izsol) = amatr(idof3,izsol) + elmat(ievat,jevat)     ! Top right

           end do
           ievat=ndime*pnode+inode
           jevat=(jnode-1)*ndime
           idof2=ndime*ndofn
           do jdime=1,ndime
              jevat=jevat+1
              idof2=idof2+1
              !$OMP ATOMIC
              amatr(idof2,izsol) = amatr(idof2,izsol) + elmat(ievat,jevat)    ! Bot. left

           end do
           idof2=ndofn*ndofn
           jevat=ndime*pnode+jnode
           !$OMP ATOMIC
           amatr(idof2,izsol) = amatr(idof2,izsol) + elmat(ievat,jevat)       ! Bot. right

        end do
     end do

  else if(iprob==2) then

     if(ndofn==1) then
        !
        ! General case: 1 unknown
        !   
        do inode = 1,pnode
           ipoin = lnode(inode)
           do jnode = 1,pnode
              jpoin = lnode(jnode)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1)
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jcolu == jpoin ) then

                 !$OMP ATOMIC
                 amatr(1,izsol) = amatr(1,izsol) + elmat(inode,jnode)

              end if
           end do
        end do

     else if(ndofn==2) then
        !
        ! General case: 2 unknowns
        !
        do inode = 1,pnode
           ipoin = lnode(inode)
           iposi = 2*inode-1
           do jnode = 1,pnode
              jpoin = lnode(jnode)
              jposi = 2*jnode-1
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while(jcolu/=jpoin .and. izsol < r_sol(ipoin+1)-1 )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jpoin == jcolu ) then

                 ievat          = iposi
                 jevat          = jposi
                 !$OMP ATOMIC
                 amatr(1,izsol) = amatr(1,izsol) + elmat(ievat,jevat)
                 jevat          = jposi+1
                 !$OMP ATOMIC
                 amatr(2,izsol) = amatr(2,izsol) + elmat(ievat,jevat)

                 ievat          = iposi+1
                 jevat          = jposi
                 !$OMP ATOMIC
                 amatr(3,izsol) = amatr(3,izsol) + elmat(ievat,jevat)
                 jevat          = jposi+1
                 !$OMP ATOMIC
                 amatr(4,izsol) = amatr(4,izsol) + elmat(ievat,jevat)   

              end if
           end do
        end do

     else if(ndofn==3) then
        !
        ! General case: 3 unknowns
        !
        do inode = 1,pnode
           ipoin = lnode(inode)
           iposi = 3*inode-2
           do jnode = 1,pnode
              jpoin = lnode(jnode)
              jposi = 3*jnode-2
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while( jcolu /= jpoin .and. izsol < r_sol(ipoin+1)-1 )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jcolu == jpoin ) then

                 ievat          = iposi
                 jevat          = jposi
                 !$OMP ATOMIC
                 amatr(1,izsol) = amatr(1,izsol) + elmat(ievat,jevat)
                 jevat          = jposi+1
                 !$OMP ATOMIC
                 amatr(2,izsol) = amatr(2,izsol) + elmat(ievat,jevat)
                 jevat          = jposi+2
                 !$OMP ATOMIC
                 amatr(3,izsol) = amatr(3,izsol) + elmat(ievat,jevat)

                 ievat          = iposi+1
                 jevat          = jposi
                 !$OMP ATOMIC
                 amatr(4,izsol) = amatr(4,izsol) + elmat(ievat,jevat)
                 jevat          = jposi+1
                 !$OMP ATOMIC
                 amatr(5,izsol) = amatr(5,izsol) + elmat(ievat,jevat)   
                 jevat          = jposi+2
                 !$OMP ATOMIC
                 amatr(6,izsol) = amatr(6,izsol) + elmat(ievat,jevat)  

                 ievat          = iposi+2
                 jevat          = jposi
                 !$OMP ATOMIC
                 amatr(7,izsol) = amatr(7,izsol) + elmat(ievat,jevat)
                 jevat          = jposi+1
                 !$OMP ATOMIC
                 amatr(8,izsol) = amatr(8,izsol) + elmat(ievat,jevat)   
                 jevat          = jposi+2
                 !$OMP ATOMIC
                 amatr(9,izsol) = amatr(9,izsol) + elmat(ievat,jevat)   

              end if
           end do
        end do

     else
        do inode=1,pnode
           ipoin = lnode(inode)
           do jnode=1,pnode
              jpoin = lnode(jnode)
              izsol = r_sol(ipoin)
              jcolu = c_sol(izsol)
              do while(jcolu/=jpoin .and. izsol < r_sol(ipoin+1)-1 )
                 izsol = izsol + 1
                 jcolu = c_sol(izsol)
              end do
              if( jpoin == jcolu ) then

                 do idofn=1,ndofn
                    ievat=(inode-1)*ndofn+idofn
                    do jdofn=1,ndofn
                       jevat=(jnode-1)*ndofn+jdofn
                       idof2=(idofn-1)*ndofn+jdofn
                       !$OMP ATOMIC
                       amatr(idof2,izsol) = amatr(idof2,izsol) + elmat(ievat,jevat)
                    end do
                 end do

              end if
           end do
        end do
     end if

  else if(iprob<0) then
     !
     ! Assemble only equation for dof NDOFI
     !
     ndofi=-iprob
     do inode=1,pnode
        ipoin = lnode(inode)
        do jnode=1,pnode
           jpoin = lnode(jnode)
           izsol = r_sol(ipoin)
           jcolu = c_sol(izsol)
           do while(jcolu/=jpoin)
              izsol = izsol + 1
              jcolu = c_sol(izsol)
           end do
           ievat=(inode-1)*ndofn+ndofi
           jevat=(jnode-1)*ndofn+ndofi

           !$OMP ATOMIC
           amatr(1,izsol) = amatr(1,izsol) + elmat(ievat,jevat) 

        end do
     end do

  else if(iprob==3) then
     !
     ! Assemble symmetric matrix
     !
     if(ndofn==1) then
        do inode=1,pnode
           ipoin = lnode(inode)
           do jnode=1,pnode
              jpoin = lnode(jnode)
              if(jpoin<=ipoin) then
                 izsol = r_sym(ipoin)
                 jcolu = c_sym(izsol)
                 do while(jcolu/=jpoin .and. izsol < r_sym(ipoin+1)-1 )
                    izsol = izsol + 1
                    jcolu = c_sym(izsol)
                 end do
                 if( jcolu == jpoin ) then

                    !$OMP ATOMIC
                    amatr(1,izsol) = amatr(1,izsol)+elmat(inode,jnode)

                 end if
              end if
           end do
        end do
     else
        do inode=1,pnode
           ipoin = lnode(inode)
           do jnode=1,pnode
              jpoin = lnode(jnode)
              if(jpoin<=ipoin) then
                 izsol = r_sym(ipoin)
                 jcolu = c_sym(izsol)
                 do while( jcolu /= jpoin .and. izsol < r_sym(ipoin+1)-1 )
                    izsol = izsol + 1
                    jcolu = c_sym(izsol)
                 end do
                 if( jcolu == jpoin ) then

                    do idofn=1,ndofn
                       ievat=(inode-1)*ndofn+idofn
                       do jdofn=1,ndofn
                          jevat=(jnode-1)*ndofn+jdofn
                          idof2=(idofn-1)*ndofn+jdofn
                          !$OMP ATOMIC
                          amatr(idof2,izsol) = amatr(idof2,izsol) + elmat(ievat,jevat)
                       end do
                    end do

                 end if
              end if
           end do
        end do
     end if
  end if

end subroutine csrase
