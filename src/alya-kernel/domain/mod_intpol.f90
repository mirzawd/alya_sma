!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_intpol

  use def_kintyp, only       :  ip,rp
  implicit none

contains


  subroutine krigin(coord,uncoo,limit,lnode,shapl,incoo)
    !--------------------------------------------------------------------------
    !****f* ibm_fielib
    ! NAME
    !    ibm_fielib
    ! DESCRIPTION
    !    Find the best element formed by a given node and their neiborgs that 
    !    contains the point of intersection with the particle surface    
    !
    !    Returns the point of prejection, a list of elemental nodes 
    !    and the value of the shape functions
    ! USED BY
    !    nastin/nsi_embedd
    !--------------------------------------------------------------------------
    use def_domain,  only    :  ndime

    real(rp),    intent(in)             :: coord(ndime,*)
    real(rp),    intent(in)             :: uncoo(ndime)
    integer(ip), intent(in)             :: limit
    integer(ip), intent(in)             :: lnode(limit)
    real(rp),    intent(out)            :: shapl(limit+1)
    real(rp),    intent(in),  optional  :: incoo(ndime)


    integer(ip)                         :: idime,iinte,jinte
!    integer(ip)                         :: dummi,infor
    integer(ip)                         :: nukno,ndata
    integer(ip), pointer                :: linde(:)   => null()
    real(rp)                            :: coori(ndime),coorj(ndime)
    real(rp),    pointer                :: covma(:,:) => null()
    real(rp),    pointer                :: covve(:)   => null()

    ndata = limit
    if( present(incoo) ) then
       ndata = ndata + 1
    end if

    ! Use linear mean for kriging
    nukno = ndata+ndime+1
    ! Use constant mean for kriging
    if (ndata < ndime+1) nukno = ndata+1

    allocate(covma(nukno,nukno))
    allocate(covve(nukno))
    allocate(linde(nukno))

    do iinte = 1,nukno
       !do iinte = 1,limit+1
       covve(iinte) = 0.0_rp
       do jinte = 1,nukno
          !do jinte = 1,limit+1
          covma(iinte,jinte) = 0.0_rp
       end do
    end do
    !
    ! Assembly the covariance matrix
    !
    do iinte = 1,ndata
       do jinte = iinte,ndata
          do idime = 1,ndime
             if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
             if (iinte >  limit) coori(idime) = incoo(idime)
             if (jinte <= limit) coorj(idime) = coord(idime,lnode(jinte))
             if (jinte >  limit) coorj(idime) = incoo(idime)
          end do
          covma(iinte,jinte) = 0.0_rp
          do idime = 1,ndime
             covma(iinte,jinte) = covma(iinte,jinte) + (coori(idime) - coorj(idime))**2.0_rp           
          end do
          covma(iinte,jinte) = covma(iinte,jinte)**1.5_rp
       end do
    end do

    do iinte = 1,ndata-1
       do jinte = iinte+1,ndata
          covma(jinte,iinte) = covma(iinte,jinte)
       end do
    end do
    !
    ! Nugget effect
    !
    !do iinte = 1,limit
    !   covma(iinte,iinte) =  covma(iinte,iinte) + (1.0_rp)**2.0_rp
    !end do
    !
    ! Assembly the mean value basis
    !
    do iinte = 1,ndata
       covma(iinte,ndata+1) = 1.0_rp
       if (nukno == ndata+ndime+1) then
          do idime = 1,ndime
             if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
             if (iinte >  limit) coori(idime) = incoo(idime)
          end do
          do idime = 1,ndime
             covma(iinte,ndata+idime+1)       = coori(idime)
          end do
       end if
    end do
    do iinte = 1,ndata
       covma(ndata+1,iinte) = covma(iinte,ndata+1)
       if (nukno == ndata+ndime+1) then
          do idime = 1,ndime     
             covma(ndata+idime+1,iinte)       = covma(iinte,ndata+idime+1)
          end do
       end if
    end do
    !
    ! Assembly the right hand side
    !
    do iinte = 1,ndata
       do idime = 1,ndime
          if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
          if (iinte >  limit) coori(idime) = incoo(idime)
       end do
       covve(iinte) = 0.0_rp
       do idime = 1,ndime
          covve(iinte) =  covve(iinte) + ( uncoo(idime)  - coori(idime) )**2.0_rp
       end do
       covve(iinte) = covve(iinte)**1.5_rp
    end do

    covve(ndata+1) = 1.0_rp
    if (nukno == ndata+ndime+1) then
       do idime = 1,ndime
          covve(ndata+idime+1)       = uncoo(idime)
       end do
    end if
    !
    ! Obtain the krigging interpolation coefficinets
    !  
    do iinte = 1,ndata
       shapl(iinte) = covve(iinte)
    end do

    deallocate(covma)
    deallocate(covve)
    deallocate(linde)

  end subroutine krigin

end module mod_intpol
