!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_evafun( x , coori ,k , M , neigh ,targe)

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  implicit none  
  integer(ip),intent(in)         :: x , k , M
  integer(ip),intent(in)         :: neigh(mnode,M)
  real(rp) ,intent(in)           :: coori(ndime)
  real(rp) ,intent(inout)        :: targe
  integer(ip)                    :: idime,nodo,tetra,ll,maxiter,iter
  real(rp)                       :: invW(ndime,k),tinvW(k,ndime),invWt(k,ndime)
  real(rp)                       :: A(ndime,k,M),S(ndime,k,M),invS(ndime,k,M),normgradif
  real(rp)                       :: detS(M),detSt(M), normS(M),norminvS(M),StS(ndime,k,M)
  real(rp)                       :: kappaS(M)
  real(rp)                       :: aux(ndime,k),aux2(ndime,k),dkdS(ndime,k,M)
  real(rp)                       :: detinvWt,eps
 

  if (ndime==3)then
     invW(1,1) = 1.0_rp
     invW(1,2) = -1.0_rp/3.0_rp*sqrt(3.0_rp)            !!!-5.77350269189626e-01!!!
     invW(1,3) = -1.0_rp/6.0_rp*sqrt(6.0_rp)             !!-4.08248290463863e-01 
     invW(2,1) = 0.0_rp
     invW(2,2) = 2.0_rp/3.0_rp*sqrt(3.0_rp)           !!!1.15470053837925e+00 
     invW(2,3) = -1.0_rp/6.0_rp*sqrt(6.0_rp)         !!!-4.08248290463863e-01
     invW(3,1) = 0.0_rp
     invW(3,2) = 0.0_rp
     invW(3,3) = 1.0_rp/2.0_rp*sqrt(6.0_rp)           !!!!1.22474487139159e+00
  else
     invW(1,1) =  1.0_rp
     invW(1,2) =  1.0_rp/2.0_rp              
     invW(2,1) =  0.0_rp
     invW(2,2) =  sqrt(3.0_rp)/2.0_rp           !!!1.15470053837925e+00 
  end if

  detinvWt = 0.0_rp
  normgradif = 10000.0_rp
  maxiter = 100_ip
  iter = 0_ip
  eps  = 0.001_rp
  targe = 0.0_rp
  !
  !Initialize
  !
  do tetra=1,M
     do nodo=1,k
        do idime=1,ndime
           A(idime,nodo,tetra) = 0.0_rp
           S(idime,nodo,tetra) = 0.0_rp
           StS (idime,nodo,tetra) = 0.0_rp
           invS(idime,nodo,tetra) = 0.0_rp
           dkdS (idime,nodo,tetra) = 0.0_rp
        end do
     end do
     detSt(tetra) = 0.0_rp
     normS(tetra) = 0.0_rp
     norminvS(tetra) = 0.0_rp
  end do

  do nodo=1,k
     do idime = 1,ndime
        tinvW (idime,nodo)= invW(nodo,idime)
     end do
  end do
  call invmtx(tinvW,invWt,detinvWt,ndime)

  do tetra=1,M
     do nodo=1,k
        do idime=1,ndime
           A(idime,nodo,tetra) = coord_ale(idime,neigh(nodo,tetra),1) - coori(idime)
        end do
     end do
     do nodo=1,k
        do idime=1,ndime
           do ll = 1,k
              S(idime,nodo,tetra) = S(idime,nodo,tetra) + A(idime,ll,tetra) * invW(ll,nodo)
           end do
        end do
     end do
     do nodo=1,k
        do idime=1,ndime
           normS(tetra) = normS(tetra) + S(idime,nodo,tetra) * S(idime,nodo,tetra)  
           do ll = 1,k
              StS (idime,nodo,tetra) = StS(idime,nodo,tetra) + S(ll,idime,tetra)* S(ll,nodo,tetra)
           end do
        end do
     end do
     normS(tetra) = sqrt(normS(tetra))
  end do

  ! targe_aux=0.0_rp
  do tetra=1,M
     do nodo=1,k
        do idime=1,ndime
           aux(idime,nodo) = S(idime,nodo,tetra)
        end do
     end do
     call invmtx(aux,aux2,detS(tetra),ndime)
     do nodo=1,k
        do idime=1,ndime
           invS(idime,nodo,tetra) = aux2(idime,nodo) 
        end do
     end do

     do nodo=1,k
        do idime=1,ndime
           norminvS(tetra) = norminvS(tetra) + invS(idime,nodo,tetra) * invS(idime,nodo,tetra)  
        end do
     end do

     norminvS(tetra) = sqrt(norminvS(tetra))
     kappaS(tetra) = normS(tetra) * norminvS(tetra)
     !if(kappaS(tetra)>1000.0_rp)then
      !  print*,'el tetra',tetra,neigh(:,tetra)
      !  print*,'antes targe:',targe
      !  print*,'S y Sinv', normS(tetra), norminvS(tetra)
     !end if
     targe = targe + kappaS(tetra)* kappaS(tetra)
     !if(kappaS(tetra)>1000)print*,'al final el targe',targe
     !targe_aux = targe_aux + kappaS(tetra)* kappaS(tetra)
     !print*,'targe despues',targe,'por el kappas',kappaS(tetra)
  end do
  !targe=targe_aux
  !print*,'targeeee',targe
 ! print*,targe

end subroutine ale_evafun
