!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_evagra ( x , k , M , neigh,gradiF,normgradif_out,kappaS_out )

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  implicit none  
  integer(ip), intent(in)        :: x , k , M
  integer(ip), intent(in)        :: neigh(mnode,M)
  real(rp),    intent(out)       :: gradiF(ndime),normgradif_out,kappaS_out(M)
  integer(ip)                    :: idime,nodo,tetra,ll,inodo,jnodo,maxiter,iter
  real(rp)                       :: invW(ndime,k),tinvW(k,ndime),invWt(k,ndime)
  real(rp)                       :: A(ndime,k,M),S(ndime,k,M),invS(ndime,k,M),normgradif
  real(rp)                       :: detS(M),detSt(M), normS(M),norminvS(M),StS(ndime,k,M)
  real(rp)                       :: invSt(ndime,k,M),kappaS(M),Sprima(ndime,k,M)
  real(rp)                       :: aux(ndime,k),aux2(ndime,k),dkdS(ndime,k,M),gradik(ndime,M)
  real(rp)                       :: detinvWt,caca(k),eps
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
  maxiter = 1000_ip
  iter = 0_ip
  eps  = 0.001_rp
  !
  !Initialize
  !
  do idime = 1,ndime
     gradiF(idime) = 0.0_rp
  end do
  do tetra=1,M
     do nodo=1,k
        do idime=1,ndime
           A(idime,nodo,tetra)  = 0.0_rp
           S(idime,nodo,tetra) = 0.0_rp
           StS (idime,nodo,tetra) = 0.0_rp
           dkdS (idime,nodo,tetra) = 0.0_rp
           caca(idime) = 0.0_rp
        end do
     end do
     do idime = 1,ndime
        gradik(idime,tetra)=0.0_rp 
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
           A(idime,nodo,tetra) = coord_ale(idime,neigh(nodo,tetra),1) - coord_ale(idime,x,1)
           !if(abs(A(idime,nodo,tetra))>=1.0e8 .or. abs(A(idime,nodo,tetra))<=1.0e-8  )print*,'las coords...',coord_ale(idime,neigh(nodo,tetra)),coord_ale(idime,x),'idime',idime,'nodo',nodo,'tetra',tetra
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
     !print*,'para cada tetra el kappa', kappaS(tetra) 
     kappaS_out(tetra) = kappaS(tetra)
     do nodo=1,k
        do idime=1,ndime
           aux(idime,nodo) = S(nodo,idime,tetra)
        end do
     end do
     call invmtx(aux,aux2,detSt(tetra),ndime)
     do nodo=1,k
        do idime=1,ndime
           invSt(idime,nodo,tetra) = aux2(idime,nodo) 
           if(idime==nodo) then
              Sprima(idime,nodo,tetra) = normS(tetra)*normS(tetra)-StS(idime,nodo,tetra)
           else
              Sprima(idime,nodo,tetra) = -StS(idime,nodo,tetra)
           end if
        end do
     end do



     do nodo=1,k
        do idime=1,ndime
           do ll =1,k
              dkdS(idime,nodo,tetra) = (normS(tetra) * normS(tetra) /&
                   (detS(tetra)*detS(tetra)*kappaS(tetra)))* S(idime,ll,tetra)*&
                   Sprima(ll,nodo,tetra)
           end do
           dkdS(idime,nodo,tetra) = dkdS(idime,nodo,tetra)- kappaS(tetra)*invSt(idime,nodo,tetra)&
                + norminvS(tetra)*norminvS(tetra)/kappaS(tetra)*&
                S(idime,nodo,tetra)

        end do
     end do
  end do

  do tetra=1,M
     do idime=1,ndime
        do inodo=1,k
           do jnodo =1,k
              gradik(idime,tetra) = gradik(idime,tetra) - dkdS(idime,inodo,tetra) * invWt(inodo,jnodo) * 1.0_rp
           end do
        end do
     end do
  end do
  !
  !Calculamos el gradiente del funcional norma 2 del condition number: gradiF
  !
  do idime = 1,ndime
     do tetra = 1,M
        gradiF(idime) = gradiF(idime) + 2.0_rp * kappaS(tetra) * gradik(idime,tetra)
     end do
  end do

  normgradif = 0.0_rp
  do idime = 1,ndime
     normgradif = normgradif + gradiF(idime) * gradiF(idime)
  end do
  normgradif     = sqrt(normgradif)
  normgradif_out = normgradif 
 
  if( normgradif /= 0.0_rp ) then
     normgradif = 1.0_rp/normgradif
     do idime = 1,ndime
        gradiF(idime) = - gradiF(idime) * normgradif
     end do
  end if

end subroutine ale_evagra
