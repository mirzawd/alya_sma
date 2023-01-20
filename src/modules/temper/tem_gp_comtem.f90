!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_gp_comtem
  !-----------------------------------------------------------------------
  !****f* Temper/solve_tem
  ! NAME 
  !    solve_tem
  ! DESCRIPTION
  !    This routine computes the temperature from the enthalpy using the 
  !    polynomial coefficients for the specific heat
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp,r1p
  use def_master, only          : sphec_gp,sphec_gp_ht,sphec_gp_lt,conce,&
       therm,conce, tempe_gp, mem_modul,modul,&
       tempe,INOTEMPTY

  use def_domain, only          : nnode,nelem,ltype,&
       ngaus,elmar,lnods,&
       mgaus,mnode,npoin

  use def_temper, only          : cfi_hmax_tem, cfi_hmin_tem,&
       kfl_adiab_tem,bvtem_tem,kfl_fixno_tem
  use mod_memory, only          : memory_alloca, memory_deallo
  use mod_physics, only         : physics_H_2_TCp

  !use def_parame
  !use def_master
  !use def_domain
  !use mod_ker_proper

  implicit none
  integer(ip)             :: inode,igaus,ielem,ipoin
  real(rp)                :: cploc(6,2),ent_loc,acval
  integer(ip)             :: lnods_loc(mnode)
  real(rp)                :: elcon(mnode)
  real(rp)                :: gpcon(mgaus)
  real(rp)                :: elent(mnode)
  real(rp)                :: gpent(mgaus)
  integer(ip)             :: pgaus,pelty,pnode
  type(r1p),pointer       :: aux_r1p(:)
  !
  ! total enthalpy equation, source term implicit in cp coefficients
  !
  elements: do ielem = 1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)
     if( pelty > 0 ) then
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)

        sphec_gp(ielem) % a = 0.0_rp
        !
        ! Gather all
        !
        lnods_loc(1:pnode) = lnods(1:pnode,ielem)
        do inode=1,pnode
           if (kfl_adiab_tem /= 0) then
              elcon(inode) = min(1.0_rp,max(0.0_rp,conce( lnods_loc(inode), 3, 1 ) ))
           endif
           elent(inode) = therm( lnods_loc(inode), 1 ) 
        enddo

        !
        ! Enthalpy on GP 
        !
        do igaus = 1,pgaus
           gpcon(igaus)         = 0.0_rp
           gpent(igaus)         = 0.0_rp
           do inode = 1,pnode
              if (kfl_adiab_tem /= 0) then
                 gpcon(igaus) = gpcon(igaus)&
                      + elmar(pelty)%shape(inode,igaus) * elcon(inode)
              endif
              gpent(igaus)         = gpent(igaus)&
                   + elmar(pelty)%shape(inode,igaus) * elent(inode)
           end do
        end do

        do igaus = 1,pgaus
           !
           ! Check if adiabatic calculation
           !
           if (kfl_adiab_tem /= 0) then           ! Enthalpy is given by the mixture fraction
              acval   = min(1.0_rp,max(0.0_rp,(gpcon(igaus))))
              ent_loc = acval * cfi_hmax_tem + (1.0_rp - acval) * cfi_hmin_tem
           else
              ent_loc = gpent(igaus)              ! Enthalpy is recomputed each time-step
           endif

           cploc(1:6,1) = sphec_gp_lt(ielem) % a(igaus,1:6,1)
           cploc(1:6,2) = sphec_gp_ht(ielem) % a(igaus,1:6,1)


           call physics_H_2_TCp(ent_loc, cploc, tempe_gp(ielem) % a(igaus,1,1), sphec_gp(ielem) % a(igaus,1,1)) 

        enddo
     endif



  end do elements
  !
  ! Put temperature on nodes
  !
  if (INOTEMPTY) then
     nullify(aux_r1p)
     call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','tem_gp_comtem',aux_r1p,max(1_ip,nelem))
     do ielem = 1,nelem
        pelty = ltype(ielem)
        pgaus = ngaus(pelty)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','tem_gp_comtem',aux_r1p(ielem)%a,pgaus)
        aux_r1p(ielem) % a = tempe_gp(ielem) % a(:,1,1)
     end do
     do ipoin=1,npoin
        tempe(ipoin,1) = 0.0_rp
     end do
     call smooth(aux_r1p, tempe(:,1))
     call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','tem_gp_comtem',aux_r1p)
  endif
  !
  ! Correct boundary values
  ! 
  do ipoin = 1,npoin
     if(kfl_fixno_tem(1,ipoin)==1) then
        tempe(ipoin,1) = bvtem_tem(1,ipoin,1)
     end if
  end do


end subroutine tem_gp_comtem
