!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine smooth(velem,vpoin)
  !------------------------------------------------------------------------
  !****f* Master/smooth
  ! NAME 
  !    smooth
  ! DESCRIPTION
  !    Smooth a variable defined on elements
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only              :  ip,rp,r1p
  use def_master, only              :  INOTMASTER,mem_modul,modul
  use def_domain, only              :  ndime,npoin,nelem,nnode,mnode
  use def_domain, only              :  lnods,ltype,coord,elmar
  use def_domain, only              :  ngaus,kfl_naxis,lelch
  use def_elmtyp  
  use mod_memchk
  implicit none
  type(r1p),intent(in)              :: velem(nelem)
  real(rp),intent(out)              :: vpoin(npoin) 
  integer(ip)                       :: ipoin,idime,inode,ielem,igaus
  integer(ip)                       :: pnode,pelty,pgaus
  integer(ip)                       :: jnode,knode
  integer(4)                        :: istat
  real(rp)                          :: detjm,gpvol,gpcar(ndime,mnode)
  real(rp)                          :: elcod(ndime,mnode)
  real(rp)                          :: xjaci(9),xjacm(9),xfact
  real(rp),     pointer             :: vmass_xxx(:)

  if( INOTMASTER ) then
     !
     ! Smooth property 
     !
     allocate( vmass_xxx(npoin) , stat = istat )
     call memchk(0_ip,istat,mem_modul(1:2,modul),'VMASS_XXX','smooth',vmass_xxx)
     !
     ! Initialization
     !
     do ipoin = 1,npoin
        vmass_xxx(ipoin) = 0.0_rp
        vpoin(ipoin)     = 0.0_rp
     end do
     !
     ! Loop over elements
     !
     elements: do ielem = 1,nelem
        pelty = ltype(ielem)

        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! Gather vectors
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           !
           ! Loop over Gauss points 
           !
           gauss_points: do igaus = 1,pgaus
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                   elcod,gpcar,detjm,xjacm,xjaci)
              gpvol = elmar(pelty) % weigp(igaus) * detjm
              if( kfl_naxis == 1 ) then
                 call runend('MOD_GRADIE: NOT CODED')
              end if
              !
              ! Extension
              !
              if( lelch(ielem) == ELEXT ) then
                 knode = 1
              else
                 knode = pnode
              end if
              !
              ! Assemble
              !
              do inode = 1,knode
                 ipoin        = lnods(inode,ielem)
                 xfact        = gpvol * elmar(pelty) % shape(inode,igaus)
                 vpoin(ipoin) = vpoin(ipoin) + xfact * velem(ielem) % a(igaus)
                 do jnode = 1,pnode
                    vmass_xxx(ipoin) = vmass_xxx(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                 end do
              end do

           end do gauss_points
        end if
     end do elements
     !
     ! Parallelization
     !
     call rhsmod(1_ip,vpoin)
     call rhsmod(1_ip,vmass_xxx)
     !
     ! Solve diagonal system
     !
     do ipoin = 1,npoin
        if( vmass_xxx(ipoin) /= 0.0_rp ) &
             vpoin(ipoin) = vpoin(ipoin) / vmass_xxx(ipoin)
     end do
     !
     ! Deallocate memory
     !
     call memchk(2_ip,istat,mem_modul(1:2,modul),'VMASS_XXX','smooth',vmass_xxx)
     deallocate(vmass_xxx,stat=istat)
     if( istat /= 0 ) call memerr(2_ip,'VMASS_XXX','smooth',0_ip)

  end if

end subroutine smooth
