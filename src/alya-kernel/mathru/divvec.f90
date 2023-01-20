!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine divvec(unkno,diunk)
  !------------------------------------------------------------------------
  !****f* Mathru/divvec
  ! NAME
  !    divvec
  ! DESCRIPTION
  !    This routine computes the divergence of a vector
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
    use def_kintyp
    use def_master, only             : INOTMASTER
    use def_domain, only             : ndime,npoin,nelem,nnode,mnode,&
         &                             lnods,ltype,coord,vmasc,elmar
    implicit none
    real(rp),   intent(in)          :: unkno(ndime,npoin)
    real(rp),   intent(out), target :: diunk(npoin)
    integer(ip)                     :: ipoin,idime,inode,ielem
    integer(ip)                     :: pnode,pelty,jnode
    real(rp)                        :: detjm,gpvol,gpcar(ndime,mnode) 
    real(rp)                        :: eldiv
    real(rp)                        :: elunk(ndime,mnode),elcod(ndime,mnode)
    real(rp)                        :: xjaci(9),xjacm(9)

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ipoin=1,npoin
          diunk(ipoin)=0.0_rp
       end do
       !
       ! Loop over elements
       !
       if( ndime==3 ) then
          do ielem=1,nelem
             pelty = ltype(ielem) 
             pnode = nnode(pelty)
             !
             ! Gather vectors
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1,inode) = coord(1,ipoin)
                elcod(2,inode) = coord(2,ipoin)
                elcod(3,inode) = coord(3,ipoin)
                elunk(1,inode) = unkno(1,ipoin)
                elunk(2,inode) = unkno(2,ipoin)
                elunk(3,inode) = unkno(3,ipoin)
             end do
             !
             ! Loop over Gauss points (which are nodes)
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                call elmder(&
                     pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty)%weigc(inode)*detjm
                !
                ! Divergence
                !                
                eldiv = 0.0_rp
                do jnode = 1,pnode
                   eldiv = eldiv + gpcar(1,jnode) * elunk(1,jnode)
                   eldiv = eldiv + gpcar(2,jnode) * elunk(2,jnode)
                   eldiv = eldiv + gpcar(3,jnode) * elunk(3,jnode)
                end do
                diunk(ipoin) = diunk(ipoin) + gpvol * eldiv
             end do
          end do
       else    ! ndime/=3
          do ielem=1,nelem
             pelty=ltype(ielem) 
             pnode=nnode(pelty)
             !
             ! Gather vectors
             !
             do inode=1,pnode
                ipoin=lnods(inode,ielem)
                elcod(1:ndime,inode)=coord(1:ndime,ipoin)
                elunk(1:ndime,inode)=unkno(1:ndime,ipoin)
             end do
             !
             ! Loop over Gauss points (which are nodes)
             !
             do inode=1,pnode
                ipoin=lnods(inode,ielem)
                call elmder(&
                     pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol=elmar(pelty)%weigc(inode)*detjm
                !
                ! Velocity strain rates
                !
                eldiv=0.0_rp
                do idime=1,ndime
                   eldiv = eldiv + &
                           dot_product(gpcar(idime,1:pnode),elunk(idime,1:pnode))
                end do
                diunk(ipoin) = diunk(ipoin) + gpvol * eldiv
             end do
          end do
       end if
       !
       ! Periodicity
       !
       call rhsmod(1_ip,diunk)
       !
       ! Solve diagonal system
       !
       do ipoin=1,npoin
          diunk(ipoin) = diunk(ipoin) / vmasc(ipoin)
       end do       
    end if

end subroutine divvec


