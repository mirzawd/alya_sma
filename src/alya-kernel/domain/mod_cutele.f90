!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_cutele

  use def_kintyp, only     :  ip,rp,lg
  use def_master, only     :  lntib,imbou,netyp
  use def_elmtyp, only     :  ELCUT,TRI03,QUA04,TET04,HEX08
  use def_elmgeo, only     :  element_type
  use def_domain
  use mod_kdtree
  use mod_elmgeo, only     :  elmgeo_shapf_deriv_heslo
  implicit none
  integer(ip)              :: pgaus_sub
  integer(ip)              :: pelty_sub
  integer(ip)              :: pnode_sub

contains

  subroutine buicut (itask)
    use def_domain, only                :  nelem

    integer(ip), intent(in)             :: itask
    integer(ip)                         :: ielem,ielem_sub

    if( itask == 1 ) then

       if( .not. associated(cutel) ) then
          allocate( cutel(nelem) )
          do ielem = 1,nelem
             nullify(cutel(ielem) % l)
             nullify(cutel(ielem) % lb)
             nullify(cutel(ielem) % linou)
             cutel(ielem) % nelem = 0_ip
             cutel(ielem) % iimbo = 0_ip
             cutel(ielem) % nboun = 0_ip
          end do
       end if

       if( ndime == 2 ) then
          pelty_sub = TRI03
       else
          pelty_sub = TET04
       end if
       pnode_sub = nnode(pelty_sub)
       pgaus_sub = ngaus(pelty_sub)


    elseif( itask == 2 ) then

       do ielem = 1,nelem
          if ( associated(cutel(ielem) % l) ) then
             do ielem_sub = 1,cutel(ielem) % nelem
                deallocate(cutel(ielem)%l(ielem_sub)%elcod)             
             end do
             deallocate(cutel(ielem)%l)
          end if
       end do
       deallocate(cutel)

    end if

  end subroutine buicut


  subroutine onecut(iimbo,ielem,norma,pcoor,messa,lcrkf_sld)
    !-----------------------------------------------------------------------
    !***
    ! NAME
    !    onecut
    ! DESCRIPTION
    !    This routines divide an element cut by a particle or a plane
    !    IIMBO > 0  : The element is cut by a particle
    !    IIMBO <= 0 : The element is cut by a a plane
    !           any other value: It doesn't delete any subelements    
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------
    implicit none
    integer(ip),  intent(in)            :: iimbo,ielem
    real(rp),     intent(inout)         :: norma(ndime),pcoor(ndime)
    character(*), intent(in), optional  :: messa
    integer(ip),  intent(in), optional  :: lcrkf_sld(*)
    integer(ip)                         :: pelty,inode,idime
    integer(ip)                         :: ielem_sub
    integer(ip)                         :: ipoin_sub,npoin_sub,ipoib_sub
    integer(ip)                         :: npoib_sub,nelem_sub
    integer(ip)                         :: lelem_sub(ndime+1,8),iboun_tot
    integer(ip)                         :: npoib_tot,npoin_tot,ielem_tot
    integer(ip)                         :: nelem_tot,lelem_tot(ndime+1,50)
    integer(ip)                         :: nboun_tot,lboun_tot(ndime+1,50)
    integer(ip)                         :: ielem_div,nelem_div
    integer(ip)                         :: eleme_div(ndime+1,6)
    integer(ip)                         :: lntib_sub(ndime+1),lntib_tot(50)
    real(rp)                            :: lcoor_sub(ndime,40)
    real(rp)                            :: lcoor_tot(ndime,50)
    real(rp)                            :: lcoob_sub(ndime,4)
    real(rp)                            :: lcoob_tot(ndime,50)
    real(rp)                            :: elcod(ndime,ndime+1),sumdi
    real(rp)                            :: pladi   
!    real(rp)                            :: pcoo2(ndime)
    logical(lg)                         :: gmove

    gmove = .true.
    if( present(messa) ) then
       if( trim(messa) == 'DO_NOT_MOVE_NODES' ) gmove = .false.
    end if
    !
    ! Subdivide the element.
    !           
    pelty = ltype(ielem)
    if (ndime == 2) then 
       if( pelty == TRI03 ) then
          nelem_div  = 1

          eleme_div(1,1) = lnods(1,ielem)
          eleme_div(2,1) = lnods(2,ielem)
          eleme_div(3,1) = lnods(3,ielem)
       elseif( pelty ==  QUA04 ) then
          nelem_div  = 2

          eleme_div(1,1) = lnods(1,ielem)
          eleme_div(2,1) = lnods(2,ielem)
          eleme_div(3,1) = lnods(3,ielem)

          eleme_div(1,2) = lnods(1,ielem)
          eleme_div(2,2) = lnods(3,ielem)
          eleme_div(3,2) = lnods(4,ielem)
       end if

    elseif (ndime == 3) then
       if( pelty == TET04 ) then
          nelem_div  = 1

          eleme_div(1,1) = lnods(1,ielem)
          eleme_div(2,1) = lnods(2,ielem)
          eleme_div(3,1) = lnods(3,ielem)
          eleme_div(4,1) = lnods(4,ielem)                
       elseif( pelty == HEX08 ) then
          nelem_div  = 6

          eleme_div(1,1) = lnods(1,ielem)
          eleme_div(2,1) = lnods(3,ielem)
          eleme_div(3,1) = lnods(4,ielem)
          eleme_div(4,1) = lnods(5,ielem)

          eleme_div(1,2) = lnods(3,ielem)
          eleme_div(2,2) = lnods(4,ielem)
          eleme_div(3,2) = lnods(5,ielem)
          eleme_div(4,2) = lnods(8,ielem)

          eleme_div(1,3) = lnods(1,ielem)
          eleme_div(2,3) = lnods(2,ielem)
          eleme_div(3,3) = lnods(3,ielem)
          eleme_div(4,3) = lnods(5,ielem)

          eleme_div(1,4) = lnods(2,ielem)
          eleme_div(2,4) = lnods(3,ielem)
          eleme_div(3,4) = lnods(5,ielem)
          eleme_div(4,4) = lnods(6,ielem)

          eleme_div(1,5) = lnods(3,ielem)
          eleme_div(2,5) = lnods(5,ielem)
          eleme_div(3,5) = lnods(6,ielem)
          eleme_div(4,5) = lnods(7,ielem)

          eleme_div(1,6) = lnods(3,ielem)
          eleme_div(2,6) = lnods(5,ielem)
          eleme_div(3,6) = lnods(7,ielem)
          eleme_div(4,6) = lnods(8,ielem)
       end if
    end if

    nelem_tot = 0
    npoin_tot = 0
    npoib_tot = 0
    nboun_tot = 0


    !do idime = 1,ndime
    !   pcoo2(idime) = pcoor(idime)
    !end do
    !if( gmove ) call moveme(ielem,norma,pcoor,lcrkf_sld)   

    do ielem_div = 1,nelem_div       
       do inode = 1,ndime+1
          do idime = 1,ndime
             elcod(idime,inode) = coord(idime,eleme_div(inode,ielem_div)) 
          end do
          if (iimbo > 0) lntib_sub(inode) = lntib(eleme_div(inode,ielem_div))
       end do

       call subdiv(iimbo,lntib_sub,elcod,norma,pcoor,lelem_sub,nelem_sub,lcoor_sub,npoin_sub,lcoob_sub,npoib_sub)       

       if (npoib_sub >= ndime) then       
          !
          ! Build the boundary
          !
          do ipoib_sub = 1,npoib_sub
             do idime = 1,ndime
                lcoob_tot(idime,ipoib_sub + npoib_tot) = lcoob_sub(idime,ipoib_sub)
             end do
          end do
          nboun_tot = nboun_tot + 1_ip
          do inode = 1,ndime
             lboun_tot(inode,nboun_tot) = inode + npoib_tot       
          end do
          !
          ! There are two boundary elements (only in 3D)
          !      
          if (npoin_sub == 4) then
             nboun_tot = nboun_tot + 1_ip
             do inode = 1,ndime
                lboun_tot(inode,nboun_tot) = inode + 1 + npoib_tot
             end do
          end if
          npoib_tot = npoib_tot + npoib_sub
       end if
       !
       ! Build the cut elements
       !
       do ielem_sub = 1,nelem_sub
          do inode = 1,ndime+1
             lelem_tot(inode,ielem_sub + nelem_tot) = lelem_sub(inode,ielem_sub) + npoin_tot
          end do
       end do
       !do inode = 1,ndime+1
       !   do idime = 1,ndime
       !      lcoor_tot(idime,inode + npoin_tot) = elcod(idime,inode)
       !   end do
       !   if (iimbo > 0) lntib_tot(inode + npoin_tot) = lntib_sub(inode)
       !end do
       do ipoin_sub = 1,npoin_sub
          do idime = 1,ndime
             lcoor_tot(idime,ipoin_sub + npoin_tot) = lcoor_sub(idime,ipoin_sub)
          end do
          if (iimbo > 0) lntib_tot(ipoin_sub +npoin_tot) = 0_rp
       end do

       nelem_tot = nelem_tot + nelem_sub
       npoin_tot = npoin_tot + npoin_sub !+ ndime+1       
    end do
    !
    ! Allocate arrays
    !                  
    if ( associated(cutel(ielem) % lb) ) then
       do iboun_tot = 1,cutel(ielem) % nboun
          deallocate(cutel(ielem)%lb(iboun_tot)%bocod)             
       end do
       deallocate(cutel(ielem)%lb)
    end if

    cutel(ielem) % nboun = nboun_tot
    allocate(cutel(ielem)%lb(nboun_tot))
    do iboun_tot = 1,nboun_tot
       allocate(cutel(ielem)%lb(iboun_tot)%bocod(ndime,ndime))
    end do

    if ( associated(cutel(ielem) % l) ) then                   
       do ielem_tot = 1,cutel(ielem) % nelem
          deallocate(cutel(ielem)%l(ielem_tot)%elcod)             
       end do
       deallocate(cutel(ielem)%l)
    end if

    cutel(ielem) % nelem = nelem_tot
    allocate(cutel(ielem)%l(nelem_tot))
    do ielem_tot = 1,nelem_tot
       allocate(cutel(ielem)%l(ielem_tot)%elcod(ndime,ndime+1))
    end do

    do iboun_tot = 1,nboun_tot
       do inode = 1,ndime
          do idime = 1,ndime
             cutel(ielem) % lb(iboun_tot) % bocod(idime,inode) = lcoob_tot(idime,lboun_tot(inode,iboun_tot))
          end do
       end do
    end do

    do ielem_tot = 1,nelem_tot
       sumdi = 1.0_rp
       do inode = 1,ndime+1
          do idime = 1,ndime
             cutel(ielem) % l(ielem_tot) % elcod(idime,inode) = lcoor_tot(idime,lelem_tot(inode,ielem_tot))
          end do
       end do
       if (iimbo > 0) then
          do inode = 1,ndime+1
             if (lntib_tot(lelem_tot(inode,ielem_tot)) < 0) sumdi = -1.0_rp
          end do
       else
          sumdi = 0.0_rp          
          do inode = 1,ndime+1
             pladi = 0.0_rp
             do idime = 1,ndime
                pladi = pladi + norma(idime) * ( cutel(ielem) % l(ielem_tot) % elcod(idime,inode) - pcoor(idime) )
             end do
             sumdi = sumdi + pladi
          end do
       end if
       if (sumdi > 0.0_rp ) then
          cutel(ielem) % l(ielem_tot) % inout = 1_ip          
       else
          cutel(ielem) % l(ielem_tot) % inout = -1_ip
       end if
    end do

  end subroutine onecut


!!$  subroutine moveme(ielem,norma,pcoor)    
!!$    integer(ip), intent(in)     :: ielem
!!$    real(rp),    intent(in)     :: norma(ndime)
!!$    real(rp),    intent(inout)  :: pcoor(ndime)
!!$    real(rp)                    :: facto,dista,plapo(ndime)
!!$    integer(ip)                 :: pelty,pnode,inode,jnode
!!$    integer(ip)                 :: idime,ifoun,inter,itera
!!$
!!$
!!$    pelty = ltype(ielem)
!!$    pnode = nnode(pelty)
!!$
!!$    dista = 0.0_rp
!!$    do inode = 1,pnode-1     
!!$       do jnode = inode,pnode
!!$          !
!!$          ! Unit normal factor to reach the current node position 
!!$          ! from the point of projection of this node in the plane
!!$          !
!!$          facto = 0.0_rp
!!$          do idime = 1,ndime
!!$             facto = facto + norma(idime)*( pcoor(idime) - coord(idime,lnods(inode,ielem)) )
!!$          end do
!!$          !do idime = 1,ndime
!!$          !   ppoin(idime) = coord(idime,lnods(inode,ielem)) -  facto*norma(idime)
!!$          !end do
!!$          !tempo = 0.0_rp
!!$          !do idime = 1,ndime
!!$          !   tempo = tempo + norma(idime) * (coord(idime,lnods(inode,ielem)) - ppoin(idime))
!!$          !end do
!!$          dista = dista + facto
!!$       end do
!!$    end do
!!$
!!$    inter = 1
!!$    itera = 0
!!$    do while (inter == 1 .and. itera < 1000)
!!$       itera = itera + 1
!!$       inter = 0
!!$       do inode = 1,pnode-1
!!$          do jnode = inode,pnode
!!$             call plaseg(norma,pcoor,coord(:,lnods(inode,ielem)),coord(:,lnods(jnode,ielem)),plapo,ifoun,dista)
!!$             if (ifoun == -1 .or. ifoun == -2) then
!!$                inter = 1
!!$                do idime = 1,ndime
!!$                   pcoor(idime) = plapo(idime)
!!$                end do
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine moveme






  subroutine subdiv(iimbo,lntib,elcod,norma,pcoor,lelem,nelem_cut,lcoor,npoin_cut,lboun,nboun_cut)    
    !-----------------------------------------------------------------------
    !***
    ! NAME
    !    subdiv
    ! DESCRIPTION
    !    This routines divide a triangular or an tetrahedron element cut by a particle or a plane
    !    IIMBO > 0  : The element is cut by a particle
    !    IIMBO <= 0 : The element is cut by a a plane
    !           any other value: It doesn't delete any subelements    
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)     :: iimbo
    integer(ip), intent(in)     :: lntib(ndime+1)
    real(rp),    intent(in)     :: elcod(ndime,ndime+1)
    real(rp),    intent(in)     :: norma(ndime),pcoor(ndime)
    integer(ip), intent(out)    :: lelem(ndime+1,8)
    real(rp),    intent(out)    :: lcoor(ndime  ,40)
    real(rp),    intent(out)    :: lboun(ndime  ,4)
    integer(ip), intent(out)    :: nelem_cut,npoin_cut,nboun_cut

    integer(ip)                :: inode,jnode,knode,idime,ifoun,ielem,nfoun,dummi
    real(rp)                   :: plapo(ndime),dummr

    nelem_cut = 1
    npoin_cut = 0
    nboun_cut = 0

    do inode = 1,ndime+1
       npoin_cut = npoin_cut + 1
       lcoor(1:ndime,npoin_cut) = elcod(1:ndime,inode)
    end do

    do inode = 1,ndime+1      
       lelem(inode,1) = inode       
    end do
    !
    ! Loop over edges INODE-JNODE
    !
    do inode = 1,ndime       
       do jnode = inode,ndime+1
          nfoun = 0
          if (iimbo <= 0) then

             call plaseg(norma,pcoor,elcod(1,inode),elcod(1,jnode),plapo,ifoun)
             !
             ! Capture the boundary
             !
             !!!OJOO if (ifoun /= 0) then             
             if (ifoun > 0) then             
                nboun_cut = nboun_cut + 1_ip
                do idime = 1,ndime
                   lboun(idime,nboun_cut) = plapo(idime)
                end do
             end if
          else
             if ((lntib(inode) < 0 .and. lntib(jnode) == 0)  .or. (lntib(jnode) < 0 .and. lntib(inode) == 0)) then
                ifoun = 1
                call faceli(imbou(iimbo) % sabox,imbou(iimbo) % blink,imbou(iimbo)%ltyib,imbou(iimbo)%lnoib,imbou(iimbo)%cooib, & 
                     elcod(1,inode),elcod(1,jnode),elcod(1,inode),mnoib,imbou(iimbo)%nboib,plapo,dummi,dummr)
             end if
          end if
          !
          ! Divide the element
          !          
          if (ifoun /= 0) then             
             npoin_cut = npoin_cut + 1
             lcoor(1:ndime,npoin_cut) = plapo(1:ndime)

             do ielem = 1,nelem_cut
                if ( lelem(inode,ielem) == inode .and. lelem(jnode,ielem) == jnode ) then
                   nfoun = nfoun + 1 
                   do knode = 1,ndime+1
                      lelem(knode,nelem_cut+nfoun) = lelem(knode,ielem) 
                   end do

                   lelem(inode,ielem)           = npoin_cut                  
                   lelem(jnode,nelem_cut+nfoun) = npoin_cut                  
                end if
             end do
             nelem_cut = nelem_cut + nfoun
          end if
       end do
    end do

  end subroutine subdiv






  subroutine moveme(ielem,norma,pcoor,lcrkf_sld)    
    !-----------------------------------------------------------------------
    !***
    ! NAME
    !    subdiv
    ! DESCRIPTION
    !    This routines divide a triangular or an tetrahedron element cut by a particle or a plane
    !    IIMBO > 0  : The element is cut by a particle
    !    IIMBO <= 0 : The element is cut by a a plane
    !           any other value: It doesn't delete any subelements    
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------
    use def_domain
    use def_master
    use def_elmtyp, only         :  ELCUT,TRI03,QUA04,TET04,HEX08

    implicit none

    integer(ip), intent(in)     :: ielem
    real(rp),    intent(inout)  :: norma(ndime),pcoor(ndime)
    integer(ip), intent(in)     :: lcrkf_sld(*)

    integer(ip)                 :: pelty,nedge_cut,ledge(2,16) 
    integer(ip)                 :: idime,iedge,ipoin,jpoin,kpoin
    integer(ip)                 :: ifoun,nnode_cut
    integer(ip)                 :: lfoun(16),inter
    integer(ip)                 :: ielty,iface,pface,inodb,pnodb        
    integer(ip)                 :: ifacg,inode,lcrac(npoin)

    real(rp)                    :: rmod,vec(ndime),coors(ndime,10)
    real(rp)                    :: proje(ndime),a(ndime),b(ndime)
    real(rp)                    :: c(ndime),xcoor(ndime),dotp1,dotp2
    real(rp)                    :: itera


    !
    ! Edge information
    !
    pelty = ltype(ielem)
    if (ndime == 2) then 
       if( pelty == TRI03 ) then
          nedge_cut = 3_ip
          ledge(1,1) = lnods(1,ielem);ledge(2,1) = lnods(2,ielem); 
          ledge(1,2) = lnods(1,ielem);ledge(2,2) = lnods(3,ielem);

          ledge(1,3) = lnods(2,ielem);ledge(2,3) = lnods(3,ielem);          
       elseif( pelty ==  QUA04 ) then
          nedge_cut = 4_ip
          ledge(1,1) = lnods(1,ielem);ledge(2,1) = lnods(2,ielem); 
          ledge(1,2) = lnods(2,ielem);ledge(2,2) = lnods(3,ielem); 
          ledge(1,3) = lnods(3,ielem);ledge(2,3) = lnods(4,ielem); 
          ledge(1,4) = lnods(4,ielem);ledge(2,4) = lnods(1,ielem);
       end if
    elseif (ndime == 3) then 

       if( pelty == TET04 ) then
          nedge_cut = 6_ip
          ledge(1,1) = lnods(1,ielem);ledge(2,1) = lnods(2,ielem); 
          ledge(1,2) = lnods(1,ielem);ledge(2,2) = lnods(3,ielem); 
          ledge(1,3) = lnods(1,ielem);ledge(2,3) = lnods(4,ielem);

          ledge(1,4) = lnods(2,ielem);ledge(2,4) = lnods(3,ielem);
          ledge(1,5) = lnods(2,ielem);ledge(2,5) = lnods(4,ielem);

          ledge(1,6) = lnods(3,ielem);ledge(2,6) = lnods(4,ielem);
       elseif( pelty ==  HEX08 ) then
          nedge_cut = 12_ip
          ledge(1,1) = lnods(1,ielem);ledge(2,1) = lnods(2,ielem); 
          ledge(1,2) = lnods(2,ielem);ledge(2,2) = lnods(3,ielem); 
          ledge(1,3) = lnods(3,ielem);ledge(2,3) = lnods(4,ielem); 
          ledge(1,4) = lnods(4,ielem);ledge(2,4) = lnods(1,ielem);

          ledge(1,5) = lnods(1,ielem);ledge(2,5) = lnods(5,ielem); 
          ledge(1,6) = lnods(2,ielem);ledge(2,6) = lnods(6,ielem); 
          ledge(1,7) = lnods(3,ielem);ledge(2,7) = lnods(7,ielem); 
          ledge(1,8) = lnods(4,ielem);ledge(2,8) = lnods(8,ielem); 

          ledge(1,9)  = lnods(5,ielem);ledge(2,9)  = lnods(6,ielem); 
          ledge(1,10) = lnods(6,ielem);ledge(2,10) = lnods(7,ielem); 
          ledge(1,11) = lnods(7,ielem);ledge(2,11) = lnods(8,ielem); 
          ledge(1,12) = lnods(8,ielem);ledge(2,12) = lnods(5,ielem); 
       end if
    end if


    do ipoin = 1,npoin
       lcrac(ipoin) = 0_ip
    end do
    !
    ! Loop over faces
    !          
    ielty = ltype(ielem)
    pface = element_type(ielty) % number_faces               


    do iface=1,pface
       ifacg = lelfa(ielem) % l(iface)
       !
       !
       if( lcrkf_sld(ifacg) > 0 ) then

          pnodb = element_type(ielty) % node_faces(iface)
          do inodb = 1,pnodb
             ! Local face node                
             inode = element_type(ielty) % list_faces(inodb,iface)
             ! Global face node (into domain)                
             ipoin  = lnods(inode,ielem)
             lcrac(ipoin) = 1_ip
          end do
       end if
    end do


    rmod   = -1.0_rp
    kpoin = 0_ip
    nnode_cut = 0_ip
    !
    ! Loop over edges
    !          
    do iedge = 1,nedge_cut
       lfoun(iedge) = 0
       ipoin = ledge(1,iedge)
       jpoin = ledge(2,iedge)
       call plaseg(norma,pcoor,coord(:,ipoin),coord(:,jpoin),proje,ifoun)

       if (ifoun /= 0) then
          !
          ! The edge contains a previous crack point
          !
          if (lcrac(ipoin) == 1_ip .and. lcrac(jpoin) == 1_ip) then
             lfoun(iedge) = 1
             ! Store the previous crack point
             nnode_cut = nnode_cut + 1
             do idime = 1,ndime
                coors(idime,nnode_cut) = proje(idime)
             end do
          elseif (ifoun == -1) then
             kpoin = ipoin
             do idime = 1,ndime
                vec(idime) = coord(idime,jpoin) - coord(idime,ipoin)
             end do
             call vecuni(ndime,vec,rmod)
          elseif (ifoun == -2) then
             kpoin = jpoin
             do idime = 1,ndime
                vec(idime) = coord(idime,ipoin) - coord(idime,jpoin)
             end do
             call vecuni(ndime,vec,rmod)
          end if
       end if
    end do
    !
    ! Move the node kpoin
    !
    if (kpoin > 0) then
       inter = 1
       itera = 0.0_rp
       do while (inter == 1 .and. itera < 10.0_rp)
          inter = 0
          itera = itera + 1.0_rp
          if (nnode_cut == 0) then
             do idime = 1,ndime
                pcoor(idime) = pcoor(idime) + itera*1.0e-2_rp*rmod*vec(idime) 
             end do
          else
             if (ndime == 2) then
                do idime = 1,ndime
                   xcoor(idime) = coord(idime,kpoin) + itera*1.0e-2_rp*rmod*vec(idime) 
                   pcoor(idime) = (xcoor(idime)+coors(idime,1))/2.0_rp
                end do
                c(1) = xcoor(2) - coors(2,1)
                c(2) = -(xcoor(1) - coors(1,1))
             elseif (ndime == 3) then
                do idime = 1,ndime
                   xcoor(idime) = coord(idime,kpoin) + itera*1.0e-2_rp*rmod*vec(idime) 
                   pcoor(idime) = (xcoor(idime)+coors(idime,1)+coors(idime,2))/3.0_rp
                end do

                do idime = 1,ndime
                   a(idime)     = xcoor(idime) - coors(idime,1)
                   b(idime)     = xcoor(idime) - coors(idime,2)
                end do
                call vecpro(c,a,b,3_ip)          
             end if
          end if

          call vecuni(ndime,c,rmod)
          do idime = 1,ndime                   
             dotp1 =  c(idime)*norma(idime)
             dotp2 = -c(idime)*norma(idime)
          end do
          if (dotp1 > dotp2) then
             do idime = 1,ndime
                norma(idime) = c(idime)
             end do
          else
             do idime = 1,ndime
                norma(idime) = -c(idime)
             end do
          end if

          do iedge = 1,nedge_cut            
             if (lfoun(iedge) == 0) then
                ipoin = ledge(1,iedge)                
                jpoin = ledge(2,iedge)             
                call plaseg(norma,pcoor,coord(:,ipoin),coord(:,jpoin),proje,ifoun)
                if (ifoun /= 0) inter = 1
             end if
          end do

       end do
    end if

  end subroutine moveme

  subroutine cutele(itask,xnorm,xpcor,lcrkf_sld)
    !-----------------------------------------------------------------------
    !***
    ! NAME
    !    cutele
    ! DESCRIPTION
    !    This routines divide the elements cut by the surface mesh defined by an skd-tree
    !    ITASK. 1:               Use a given plane normal and point
    !           2:               Use least squares to determine the plane normal and point (with imbou)
    !           3:               Use the intersetion with the particle
    ! USED BY
    !    domain
    !***
    !-----------------------------------------------------------------------

    use def_domain, only                :  ndime,nelem,lelch
    use def_master, only                :  lntib

    integer(ip), intent(in)             :: itask
    real(rp),    intent(inout)          :: xnorm(ndime,*)
    real(rp),    intent(inout)          :: xpcor(ndime,*)
    integer(ip),    intent(in),optional    :: lcrkf_sld(*)

    integer(ip)                :: idime,ielem,inode,pnode
    integer(ip)                :: ipoin,dummi,iimbo

    real(rp)                   :: norma(ndime),pcoor(ndime)
    real(rp), pointer          :: posit(:,:) => null()
    real(rp), pointer          :: value(:)   => null()
    real(rp), pointer          :: temp2(:,:) => null()
    real(rp)                   :: deter,tempo(ndime,ndime),invte(ndime,ndime)
    real(rp)                   :: coeff(ndime),prcod(ndime,mnode)
    real(rp)                   :: dista,proje(ndime),mindi,maxdi,incoo(ndime),excoo(ndime)
    real(rp)                   :: prnor(ndime)
    real(rp)                   :: pladi


    do ielem = 1,nelem          
       if ( lelch(ielem)  == ELCUT ) then

          if (itask == 1) then
             !
             ! Use a given plane normal and point
             !
             do idime = 1,ndime
                norma(idime) = xnorm(idime,ielem)
                pcoor(idime) = xpcor(idime,ielem)
             end do
             call onecut(0_ip,ielem,xnorm(:,ielem),xpcor(:,ielem),'',lcrkf_sld)
             do idime = 1,ndime
                xnorm(idime,ielem) = norma(idime)
                xpcor(idime,ielem) = pcoor(idime)
             end do
          elseif (itask == 2) then
             !
             ! Use the intersection with the particle
             !             
             do inode = 1,pnode
                ipoin  = lnods(inode,ielem)
                if (lntib(ipoin) < 0) then
                   iimbo =  abs(lntib(ipoin))
                end if
             end do
             do idime = 1,ndime
                norma(idime) = 0.0_rp
                pcoor(idime) = 0.0_rp
             end do
             call onecut(iimbo,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
          elseif (itask == 3) then
             !
             ! Least square: assembly
             !
             mindi = 0.0_rp
             maxdi = 0.0_rp
             do idime = 1,ndime
                prnor(idime) = 0.0_rp
             end do
             do inode = 1,pnode
                ipoin  = lnods(inode,ielem)
                if (lntib(ipoin) < 0) then
                   iimbo =  abs(lntib(ipoin))
                end if
             end do

             do inode = 1,pnode
                !
                ! POSITIONS: COORDINATES
                !     
                call dpopar(                                                                          &
                     ndime, imbou(iimbo) % cooib(1:ndime,ipoin) , imbou(iimbo) % npoib , mnoib ,       &
                     imbou(iimbo) % nboib , 1.0e10_rp , imbou(iimbo) % ltyib , imbou(iimbo) % lnoib , &
                     imbou(iimbo) % cooib , dista , norma , proje , dummi , imbou(iimbo) % fabox ,    &
                     imbou(iimbo) % sabox , imbou(iimbo) % blink , imbou(iimbo) % struc ,             &
                     imbou(iimbo) % ldist , imbou(iimbo) % lnele )
                do idime = 1,ndime
                   prnor(idime)       = prnor(idime) + norma(idime)
                   prcod(idime,inode) = proje(idime)
                end do

                if ( dista < mindi  ) then
                   do idime = 1,ndime
                      incoo(idime) = coord(idime,ipoin)
                   end do
                   mindi = dista          
                elseif ( dista > maxdi ) then
                   do idime = 1,ndime
                      excoo(idime) = coord(idime,ipoin)
                   end do
                   maxdi = dista
                end if


             end do

             allocate(posit(pnode,ndime))
             allocate(value(pnode))
             allocate(temp2(ndime,pnode))
             do inode = 1,pnode
                posit(inode,1) = 1.0_rp
                if     (abs(prnor(1)) >= max(abs(prnor(2)),abs(prnor(ndime)))) then
                   value(inode)       = prcod(1,inode)                   
                   posit(inode,ndime) = prcod(ndime,inode)
                   posit(inode,2)     = prcod(2,inode)
                elseif (abs(prnor(2)) >= max(abs(prnor(1)),abs(prnor(ndime)))) then
                   value(inode)       = prcod(2,inode)
                   posit(inode,ndime) = prcod(ndime,inode)                   
                   posit(inode,2)     = prcod(1,inode)
                elseif (abs(prnor(ndime)) >= max(abs(prnor(1)),abs(prnor(2)))) then
                   value(inode)       = prcod(ndime,inode)                   
                   posit(inode,ndime) = prcod(2,inode)
                   posit(inode,2)     = prcod(1,inode)
                end if
             end do
             !
             ! Least square: solution
             !                         
             call mbmatb(tempo,posit,posit,ndime,ndime,pnode)
             call invmtx(tempo,invte,deter,ndime)           
             call mbmabt(temp2,invte,posit,ndime,pnode,ndime)
             call mbvab0(coeff,temp2,value,ndime,pnode)   

             pcoor(1)     = 0.0_rp
             pcoor(2)     = 0.0_rp
             pcoor(ndime) = 0.0_rp             
             if     (abs(prnor(1)) >= max(abs(prnor(2)),abs(prnor(ndime)))) then
                norma(1)     = -1.0_rp 
                norma(2)     = coeff(2)
                norma(ndime) = coeff(ndime)
                pcoor(1)     = coeff(1)
             elseif (abs(prnor(2)) >= max(abs(prnor(1)),abs(prnor(ndime)))) then
                norma(1)     = coeff(2)
                norma(ndime) = coeff(ndime)
                norma(2)     = -1.0_rp
                pcoor(2)     = coeff(1)
             elseif (abs(prnor(ndime)) >= max(abs(prnor(1)),abs(prnor(2)))) then
                norma(1)     = coeff(2) 
                norma(2)     = coeff(ndime)
                norma(ndime) = -1.0_rp
                pcoor(ndime) = coeff(1)
             end if

             if (abs(mindi) > abs(maxdi)) then
                pladi = 0.0_rp
                do idime = 1,ndime
                   pladi = pladi + norma(idime) * ( incoo(idime) - pcoor(idime) )
                end do
                !
                ! Check if the normal is exterior
                !
                if (pladi > 0.0_rp) then
                   norma(1) = -norma(1)
                   norma(2) = -norma(2)
                   if (ndime == 3) norma(3) = -norma(3)
                end if
             elseif (abs(maxdi) > abs(mindi)) then
                pladi = 0.0_rp
                do idime = 1,ndime
                   pladi = pladi + norma(idime) * ( excoo(idime) - pcoor(idime) )
                end do
                !
                ! Check if the normal is exterior
                !
                if (pladi < 0.0_rp) then
                   norma(1) = -norma(1)
                   norma(2) = -norma(2)
                   if (ndime == 3) norma(3) = -norma(3)
                end if
             end if

             deallocate(posit)
             deallocate(value)
             deallocate(temp2)
             call onecut(0_ip,ielem,norma,pcoor,'DO_NOT_MOVE_NODES')
          end if
       end if
    end do
  end subroutine cutele

  subroutine plaseg(norma,pcoor,icoor,jcoor,plapo,ifoun,facto)
    use def_domain, only              :  ndime
    use def_master, only              :  zeror
    real(rp),    intent(in)           :: norma(ndime),pcoor(ndime)
    real(rp),    intent(in)           :: icoor(ndime),jcoor(ndime)
    real(rp),    intent(out)          :: plapo(ndime)
    integer(ip), intent(out)          :: ifoun
    real(rp),    intent(in), optional :: facto


    integer(ip)                      :: idime
    real(rp)                         :: numer,denom,dista,dist2,toler
    !
    ! Intersection between the plane and the line
    !  
    toler = 0.02_rp
    numer = 0.0_rp
    denom = 0.0_rp
    do idime = 1,ndime
       numer = numer + ( norma(idime) * ( pcoor(idime) - icoor(idime) )  )
       denom = denom + ( norma(idime) * ( jcoor(idime) - icoor(idime) )  )
    end do
    if( abs(denom) < 1.0e-5_rp ) then
       ifoun = 0
    else
       dista = numer / denom 
       do idime = 1,ndime
          plapo(idime) = icoor(idime) + dista * ( jcoor(idime) - icoor(idime) )
       end do
       ifoun = 0_ip
       if( dista >= toler .and. dista <= 1.0_rp-toler ) then
          ifoun =  1_ip
       else if( dista >= 0.0_rp .and. dista < toler ) then
          if( present(facto) ) then
             dist2 = 0.0_rp
             do idime = 1,ndime
                dist2 = dist2 + sign(1.0_rp,facto)*norma(idime)*(jcoor(idime) - icoor(idime))
             end do
             if (dist2 > 0) then
                do idime = 1,ndime
                   plapo(idime) = icoor(idime) + 2.0_rp*toler * ( jcoor(idime) - icoor(idime) )
                end do          
             else
                do idime = 1,ndime
                   plapo(idime) = icoor(idime)
                end do          
             end if             
          end if
          ifoun = -1_ip
       else if( dista > 1.0_rp-toler .and. dista <= 1.0_rp ) then
          if( present(facto) ) then
             dist2 = 0.0_rp
             do idime = 1,ndime
                dist2 = dist2 + sign(1.0_rp,facto)*norma(idime)*(jcoor(idime) - icoor(idime))
             end do
             if (dist2 > 0) then
                do idime = 1,ndime
                   plapo(idime) = icoor(idime) + (1.0_rp-2.0_rp*toler) * ( jcoor(idime) - icoor(idime) )                  
                end do          
             else
                do idime = 1,ndime
                   plapo(idime) = jcoor(idime)
                end do          
             end if
          end if
          ifoun = -2_ip
       end if
    end if

  end subroutine plaseg

  subroutine elmcar_cut(&
       pelty,pnode,plapl,ielem,elcod,pgaus,gpvol,gpsha,gpder,gpcar,gphes)

    !-----------------------------------------------------------------------
    !****f* domain/elmcar_cut
    ! NAME 
    !    elmcar_cut
    ! DESCRIPTION
    ! USES
    ! USED BY
    !    nsi_matrix
    !***
    !-----------------------------------------------------------------------
    use def_kintyp, only       :  ip,rp
    use def_elmtyp, only       :  TRI03,TET04
    use def_domain, only       :  ltopo,mgaus,ndime,ntens
    use def_domain, only       :  mnode,elmar
    use mod_elmgeo, only       :  elmgeo_natural_coordinates 
!    use def_domain, only       :  lnods
!    use def_master, only       :  fleve,kfl_paral
    use def_isoparametric, only : 
    implicit none
    integer(ip), intent(in)    :: pelty
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: plapl
    integer(ip), intent(in)    :: ielem 
    real(rp),    intent(in)    :: elcod(ndime,pnode)

    integer(ip), intent(inout) :: pgaus

    real(rp),    intent(out)   :: gpvol(mgaus)
    real(rp),    intent(out)   :: gpsha(pnode,mgaus)
    real(rp),    intent(out)   :: gpder(ndime,pnode,*)
    real(rp),    intent(out)   :: gpcar(ndime,mnode,mgaus)
    real(rp),    intent(out)   :: gphes(ntens,mnode,mgaus)

    integer(ip)                :: ptopo,ifoun,idime,inode
    integer(ip)                :: ielem_sub,igaus_sub
    integer(ip)                :: igaus!,ipoin
    real(rp)                   :: lmini,lmaxi,xfact
    real(rp)                   :: coloc(3),coglo(3)
    real(rp)                   :: deri1(ndime,mnode)
    real(rp)                   :: shap1(mnode)

    real(rp)                   :: gpwei(mgaus)
    real(rp)                   :: heslo(ntens,pnode,mgaus)

    real(rp)                   :: sumwe,dummr(9)
    !
    ! Sum of weight of master element
    !
    sumwe = 0.0_rp
    do igaus = 1,pgaus
       sumwe = sumwe + elmar(pelty) % weigp(igaus)
    end do

    ptopo     = ltopo(pelty) 
    lmini     = 0.0_rp
    lmaxi     = 1.0_rp
    pgaus     = 0

    do ielem_sub = 1,cutel(ielem) % nelem
       !
       ! XFACT: volume of sub element
       !
       call jacdet(&
            ndime,pnode_sub,cutel(ielem) % l(ielem_sub) % elcod,&
            elmar(pelty_sub)%dercg,dummr,xfact)
       xfact = xfact * elmar(pelty_sub) % weicg

       !call elmgeo_natural_coordinates(    &
       !     ndime,pelty,pnode,elcod,shap1, &
       !     deri1,cutel(ielem) % l(ielem_sub) % elcod(1:ndime,1),coloc,ifoun)
       !if( kfl_paral==1) print*,'TRI NODE1=',ielem_sub,coloc(1:2)
       !call elmgeo_natural_coordinates(    &
       !     ndime,pelty,pnode,elcod,shap1, &
       !     deri1,cutel(ielem) % l(ielem_sub) % elcod(1:ndime,2),coloc,ifoun)
       !if( kfl_paral==1) print*,'TRI NODE2=',ielem_sub,coloc(1:2) 
       !call elmgeo_natural_coordinates(    &
       !     ndime,pelty,pnode,elcod,shap1, &
       !     deri1,cutel(ielem) % l(ielem_sub) % elcod(1:ndime,3),coloc,ifoun)
       !if( kfl_paral==1) print*,'TRI NODE3=',ielem_sub,coloc(1:2) 
       
       do igaus_sub = 1,pgaus_sub
          !
          ! COGLO: Global coordinate of sub-element Gauss point
          !
          coglo = 0.0_rp
          do inode = 1,pnode_sub
             do idime = 1,ndime
                coglo(idime) = coglo(idime) &
                     + cutel(ielem) % l(ielem_sub) % elcod(idime,inode) &
                     * elmar(pelty_sub) % shape(inode,igaus_sub)
             end do
          end do
          !
          ! COLOC: Local coordinate in master-element
          !
          call elmgeo_natural_coordinates(    &
               ndime,pelty,pnode,elcod,shap1, &
               deri1,coglo,coloc,ifoun)
          
          if( ifoun == 0 ) then
             print*,'error 1=',ielem,lelch(ielem),pelty_sub,pnode_sub
             print*,'error 2=',coglo(1:ndime)
             call runend('OUPS mod_cutele')
          end if
          !
          ! Compute shape function and derivatives of master-element at COLOC
          !        
          pgaus        = pgaus + 1
          gpwei(pgaus) = elmar(pelty_sub) % weigp(igaus_sub) * xfact

          call elmgeo_shapf_deriv_heslo(&
               ndime,pnode,coloc,gpsha(:,pgaus),&
               gpder(:,:,pgaus),heslo(:,:,pgaus))
          
       end do

    end do
    if( associated(cutel(ielem) % linou) ) then                   
       deallocate(cutel(ielem)%linou)
    end if
    allocate(cutel(ielem) % linou(pgaus))
    igaus = 0
    do ielem_sub = 1,cutel(ielem) % nelem       
       do igaus_sub = 1,pgaus_sub
          igaus = igaus + 1
          cutel(ielem) % linou(igaus) = cutel(ielem) % l(ielem_sub) % inout        
       end do
    end do
    
    xfact = 0.0_rp
    do igaus = 1,pgaus
       xfact = xfact + gpwei(igaus)
    end do
    xfact = sumwe/xfact
    do igaus = 1,pgaus
       gpwei(igaus) = gpwei(igaus) * xfact
    end do    

    call elmcar(&
         pnode,pgaus,plapl,gpwei,gpsha,gpder,heslo,elcod,gpvol,gpcar,&
         gphes,ielem)

!!$    if( 1 == 1 ) then ! TESTEO
!!$       do igaus = 1,pgaus
!!$          xfact = 0.0_rp
!!$          dummr = 0.0_rp
!!$          do inode = 1,pnode
!!$             ipoin = lnods(inode,ielem)
!!$             xfact = xfact + fleve(ipoin,1) * gpsha(inode,igaus)
!!$             dummr(1:ndime) = dummr(1:ndime) + coord(1:ndime,ipoin) * gpsha(inode,igaus)
!!$          end do
!!$          if( xfact < 0.0_rp ) then
!!$          !   gpvol(igaus) = 0.0_rp
!!$          !   print*,dummr(1:ndime)
!!$          end if
!!$       end do
!!$    end if

  end subroutine elmcar_cut

 subroutine elmder_cut(&
       pelty,pnode,plapl,ielem,elcod,pgaus,gpvol,gpsha,gpcar,gphes)

    !-----------------------------------------------------------------------
    !****f* domain/elmder_cut
    ! NAME 
    !    elmder_cut
    ! DESCRIPTION
    ! USES
    ! USED BY
    !    elmder
    !***
    !-----------------------------------------------------------------------
    use def_kintyp,        only : ip,rp
    use def_elmtyp,        only : TRI03,TET04
    use def_domain,        only : ltopo,mgaus,ndime,ntens
    use def_domain,        only : mnode,elmar
    use mod_elmgeo,        only : elmgeo_natural_coordinates
    use def_isoparametric, only : shape2
    use def_isoparametric, only : shape3
    implicit none
    integer(ip), intent(in)    :: pelty
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: plapl
    integer(ip), intent(in)    :: ielem 
    real(rp),    intent(in)    :: elcod(ndime,pnode)

    integer(ip), intent(out)   :: pgaus
    real(rp),    intent(out)   :: gpvol(mgaus)
    real(rp),    intent(out)   :: gpsha(pnode,mgaus)
    real(rp),    intent(out)   :: gpcar(ndime,mnode,mgaus)
    real(rp),    intent(out)   :: gphes(ntens,mnode,mgaus)

    integer(ip)                :: ptopo,ifoun,idime,inode
    integer(ip)                :: ielem_sub,igaus_sub
    integer(ip)                :: ierro
    real(rp)                   :: lmini,lmaxi,xfact
    real(rp)                   :: coloc(3),coglo(3)
    real(rp)                   :: deri1(ndime,mnode)
    real(rp)                   :: shap1(mnode)

    real(rp)                   :: gpwei(mgaus)
    real(rp)                   :: deriv(ndime,pnode,mgaus)
    real(rp)                   :: heslo(ntens,pnode,mgaus)

    ptopo     = ltopo(pelty) 
    lmini     = 0.0_rp
    lmaxi     = 1.0_rp

    !xfact     = 1.0_rp / real(cutel(ielem) % nelem,rp)
    xfact     = 1.0_rp 
    pgaus     = 0

    do ielem_sub = 1,cutel(ielem) % nelem

       do igaus_sub = 1,pgaus_sub
          !
          ! COGLO: Global coordinate of sub-element Gauss point
          !
          do idime = 1,ndime
             coglo(idime) = 0.0_rp
          end do
          do inode = 1,pnode_sub
             do idime = 1,ndime
                coglo(idime) = coglo(idime) + cutel(ielem)%l(ielem_sub) % elcod(idime,inode) &
                     * elmar(pelty_sub) % shape(inode,igaus_sub)
             end do
          end do
          !
          ! COLOC: Local coordinate in master-element
          !
          call elmgeo_natural_coordinates(    &
               ndime,pelty,pnode,elcod,shap1, &
               deri1,coglo,coloc,ifoun)
          !call elsest_chkelm(&
          !     ndime,ptopo,pnode,elcod,shap1,deri1,&
          !     coglo,coloc,ifoun,lmini,lmaxi)

          if( ifoun == 0 ) call runend('OUPS mod_cutele 2')
          !
          ! Compute shpae function and derivatives of master-element at COLOC
          !        
          pgaus        = pgaus + 1
          gpwei(pgaus) = elmar(pelty_sub) % weigp(igaus_sub) * xfact
          ierro        = 0
          
          call elmgeo_shapf_deriv_heslo(&
               ndime,pnode,coloc,gpsha(:,pgaus),&
               deriv(:,:,pgaus),heslo(:,:,pgaus),ierro)

       end do

    end do 
 
    !do igaus = 1,pgaus
    !   call elmder(&
    !        pnode,ndime,deriv(1,1,igaus),&
    !        elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
    !end do

  end subroutine elmder_cut

end module mod_cutele
