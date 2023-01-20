!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_evalfa( x ,k,M,neigh,gradiF,coord_new)

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  implicit none  
  integer(ip),intent(in)         :: x ,k,M
  integer(ip),intent(in)         :: neigh(mnode,M) 
  real(rp) ,intent(in)           :: gradiF(ndime)
  real(rp) ,intent(inout)        :: coord_new(ndime) !!!solo out!
  integer(ip)                    :: lcara(3),ipoin,jpoin,qpoin,idime,jdime,flag
  integer(ip)                    :: inode,tetraop,jzdom,kelem,pnode
  integer(ip)                    :: ii,kk,kk_min,paso 
  real(rp)                       :: bocod(ndime,3),codel(ndime,4),cosa,cosa_max
  real(rp)                       :: numer,denom,vec12(ndime),vec13(ndime),alfa
  real(rp)                       :: coord_ini(ndime),coord_aux(ndime),coord_aux2(ndime),gradiF_local(ndime)
  real(rp)                       :: normal(ndime,ndime),norca(ndime),norma,t,numer1,numer2,baric(ndime)
  real(rp)                       :: targe_new,targe_ini,targe_min,dummr

  !testear = 0
  tetraop = 0_ip
  ipoin = x 
  paso = 10_ip
  !bucle tetras conectados a ipoin
  cosa_max  = -10.0_rp
  numer = 0.0_rp
  numer1=0.0_rp
  numer2=0.0_rp
  denom = 0.0_rp
  do jdime=1,ndime
     do idime=1,ndime
        normal(idime,jdime) = 0.0_rp
     end do
     norca(jdime) = -10.0_rp
  end do
  targe_new = 0.0_rp
  targe_ini = 0.0_rp
  targe_min = 10000.0_rp

  do jzdom=pelpo(ipoin), pelpo(ipoin+1)-1
     kelem = lelpo(jzdom)
     if(ltype(kelem)>0)then
        pnode = nnode(ltype(kelem))
        do inode=1,pnode
           do idime=1,ndime
              codel(idime,inode) = coord_ale(idime,lnods(inode,kelem),1)
           end do
        end do

        ii = 0
        do inode = 1,pnode
           jpoin = lnods(inode,kelem)
           if(jpoin/=ipoin)then
              ii = ii + 1
              lcara(ii) = jpoin
           end if
        end do
        do idime = 1,ndime
           bocod(idime,1)  = coord_ale(idime,lcara(1),1)
        end do
        do idime = 1,ndime
           bocod(idime,2)  = coord_ale(idime,lcara(2),1)
        end do
        do idime = 1,ndime
           vec12(idime) = bocod(idime,2)-bocod(idime,1)
        end do

        if( ndime == 3 ) then
           do idime = 1,ndime
              bocod(idime,3)  = coord_ale(idime,lcara(3),1)
           end do
           do idime = 1,ndime
              vec13(idime) = bocod(idime,3)-bocod(idime,1)
           end do
           normal(1,3) = vec12(2)*vec13(3) - vec12(3)*vec13(2)
           normal(2,3) = vec12(3)*vec13(1) - vec12(1)*vec13(3)
           normal(3,3) = vec12(1)*vec13(2) - vec12(2)*vec13(1)          
           norma       = sqrt(normal(1,3)*normal(1,3) + normal(2,3)*normal(2,3) + normal(3,3)*normal(3,3))
           norma       = 1.0_rp/norma
           normal(1,3) = normal(1,3)*norma
           normal(2,3) = normal(2,3)*norma
           normal(3,3) = normal(3,3)*norma
        else
           normal(1,2) = -vec12(2)
           normal(2,2) =  vec12(1)
           norma = sqrt( normal(1,2)*normal(1,2) + normal(2,2)*normal(2,2) )
           norma = 1.0_rp/norma
           normal(1,2) = normal(1,2)*norma
           normal(2,2) = normal(2,2)*norma
        end if

        call chenor(pnode,normal,bocod(1:ndime,1:pnode-1),codel)
        !print*,coord_ale(1,ipoin)+normal(1,3),coord_ale(2,ipoin)+normal(2,3),coord_ale(3,ipoin)+normal(3,3)

        !
        !producto escalar entre el vector p y la normal cara opuesta a ipoin
        !
        cosa = 0.0_rp
        do idime = 1,ndime
           cosa = cosa + gradiF(idime)*normal(idime,ndime)
        end do

        if(cosa_max < cosa) then
           tetraop = kelem 
           cosa_max = cosa
           !!centro la normal en el baricentro de la cara!!!no es necesario...norca unitario
           baric = 0.0_rp
           do idime=1,ndime
              do inode=1,pnode-1
                 baric(idime) = baric(idime) + coord_ale(idime,lcara(inode),1) 
              end do
           end do
           dummr = 1.0_rp/(pnode-1)

           do idime=1,ndime
              !baric(idime) = baric(idime)*dummr
              !norca(idime) = baric(idime)+normal(idime,3) 
              !norca(idime) = coord_ale(idime,x)+normal(idime,3)
              norca(idime) = normal(idime,ndime)
           end do
        end if

     end if
  end do
  !
  ! Distancia entre ipoin y la cara:alpha
  !
  ii = 1
  do while( ii <= pnode )
     if(lnods(ii,tetraop) /= ipoin) then
        qpoin = lnods(ii,tetraop)
        ii = 4
     end if
     ii = ii + 1
  end do
  
  do idime = 1,ndime
     numer  = numer + ( norca(idime)*coord_ale(idime,qpoin,1) - &
          norca(idime)*coord_ale(idime,ipoin,1))
     denom = denom + norca(idime)*gradiF(idime)
  end do
  alfa = numer/denom 
  alfa = alfa/2.0_rp !!!!!!  

  do idime = 1,ndime
     coord_ini(idime) = coord_ale(idime,ipoin,1)
  end do
  !print*,'dentro de evalf llamo a evafun'
  call ale_evafun( x ,coord_ini,k , M , neigh ,targe_ini)
  !print*, 'targe_ini',targe_ini
  targe_min = targe_ini

  kk_min = -1
  flag   = 0
  kk     = 0
  !maxiter = 3
  !iter = 0
  !do while (iter<maxiter)
  !   iter = iter+1
  t =  alfa/real(paso,rp) 
  do idime = 1,ndime
     gradiF_local(idime) = t * gradiF(idime)
  end do

  do while( kk < paso )
     kk = kk + 1
     !print*,'CAMBIO COORDSSSSSSSSS', kk , gradiF_local
     do idime = 1,ndime
        coord_aux(idime) = coord_ale(idime,ipoin,1) + real(kk,rp) * gradiF_local(idime) 
        !coord_aux(idime) =  kk * gradiF_local(idime) 
     end do
     !print*,'coords aux para el paso',kk,'son',coord_aux(:)

     !print*,'EL NUEVO COORDS COORDSSSSSSSSS', coord_aux

     !print*,targe_new
     !print*,'dentro del BUCLE-kk llamo a evafun'
     call ale_evafun( x , coord_aux ,k , M , neigh ,targe_new) 

     if(targe_new<targe_min)then
        flag = 1
        targe_min = targe_new
        kk_min    = kk
        do idime = 1,ndime
           coord_aux2(idime) = coord_aux(idime)
        end do
        !!!!!kk = paso
     end if
     !print*,'targe y kk',targe_min,kk_min
     if(flag == 0 .and. kk >=10)then
        kk = paso!!!SALGO SI NO MEJORO EN 10 PASOS
        !print*,'ME PIROOOOOOOOOOOOOOOOOOOOOOOOO'
     end if
  end do
  !print*,'FLAGGGGGGGGGGGGGGGGG',flag,coord_new
  if (flag == 1) then
     do idime = 1,ndime
        coord_new(idime) = coord_aux2(idime) 
     end do
     !print*,targe_min
     !   iter = maxiter
     !else
     !   paso = paso*10
  end if
  !end do
  ! print*,'coord_ini',coord_ini
  ! print*,'coord_new',coord_new
end subroutine ale_evalfa
