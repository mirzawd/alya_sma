!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_cshder(itask)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_cshder
  ! NAME
  !    ale_cshder
  ! DESCRIPTION
  !    This routine calculates some element arrays
  ! OUTPUT
  !    - For each element type, using user integration rule: 
  !      WEIGP(ngaus)
  !      SHAPE(nnode,ngaus)
  !      DERIV(ndime,nnode,ngaus)
  !      HESLO(ntens,nnode,ngaus)
  !    - For each element type, using a closed integration rule:
  !      WEIGC(nnode)
  !      SHAPC(nnode,nnode)
  !      DERIC(ndime,nnode,nnode)
  !      HESLC(ntens,nnode,nnode)
  !    - Center of gravity:
  !      SHACG(nnode)
  !      DERCG(ndime,nnode)
  !      WEICG
  !    - Element Gauss points to nodes:
  !      SHAGA(ngaus,nnode) 
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use mod_memchk 
  use def_alefor
  use mod_rulepw,     only : rulepw
  use def_quadrature, only : GAUSS_LEGENDRE_RULE
  use mod_elmgeo,     only : elmgeo_shapf_deriv_heslo
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielty,mdime,iboun,iimbo,pblty
  integer(ip)             :: pnode,pgaus,pdime,prule,ptopo,pquad
  real(rp)                :: hesib(ntens,max(1_ip,mnodi,mnoib),max(1_ip,mgaib))
  real(rp)                :: posib(ndime,mnodi+mgaib)
  integer(ip)             :: ierro
  integer(4)              :: istat
  !
  ! Element shape function and derivatives: SHAIB,DERIB,HESIB,WEIIB
  !
  if( nrbod > 0 ) then

     do ielty = 1,nelty

        if( lexib(ielty) == itask ) then

           pnode = nnode(ielty)
           pgaus = ngaib(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)
           prule = lruib(ielty)
           ptopo = ltopo(ielty)
           pquad = GAUSS_LEGENDRE_RULE
           
           allocate(elmar(ielty)%shaib(pnode,pgaus),stat=istat)
           call memchk(zero,istat,memor_dom,'SHAIB','ale_cshder',elmar(ielty)%shaib)
           allocate(elmar(ielty)%derib(mdime,pnode,pgaus),stat=istat)
           call memchk(zero,istat,memor_dom,'DERIB','ale_cshder',elmar(ielty)%derib)
           allocate(elmar(ielty)%weiib(pgaus),stat=istat) 
           call memchk(zero,istat,memor_dom,'WEIIB','ale_cshder',elmar(ielty)%weiib)

           call rulepw(pdime,pgaus,ptopo,pquad,posib,elmar(ielty)%weiib,ierro)
           
           call elmgeo_shapf_deriv_heslo(&
                pdime,pnode,pgaus,posib,elmar(ielty)%shaib,&
                elmar(ielty)%derib,hesib,ierro)

        end if
     end do

  end if
  !
  ! LNUIB
  !
  do iimbo = 1,nrbod
     do iboun = 1,rbbou(iimbo) % nboib
        pblty = rbbou(iimbo) % ltyib(iboun)
        lnuib(pblty) = lnuib(pblty) + 1
     end do
  end do


end subroutine ale_cshder
