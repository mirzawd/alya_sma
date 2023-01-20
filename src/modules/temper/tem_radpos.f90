!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_radpos

!------------------------------------------------------------------------
!
! This routine dumps the radiation geometry and results
!      
!------------------------------------------------------------------------
  use      def_domain
  use      def_master
  use      def_temper
  use      def_inpout
  use      mod_iofile
  use      mod_memchk
  implicit none
  integer(ip)    :: idime,ipoin,inodb,iboun,ielty
  integer(ip)    :: iview,jboun,lboun
  character(13)  :: elemt
  character(150) :: dumml
  character(8)   :: state
  character(20)  :: wopos

  if(kfl_viewf_tem(1)==0) return
  kfl_radia_tem=0                  ! Do not postprocess again
!
! Geometry: loop over element types
!      
  call iofile(zero,lun_ramsh_tem,fil_ramsh_tem,'POST-PROCESS')
  do ielty=1,nelty
     if(ndime==2) then
        elemt='Linear' 
     else
        if(nnode(ielty)==3) then
           elemt='Triangle' 
        else
           elemt='Quadrilateral'
        end if
     end if
     !
     ! Header
     !
     write(lun_ramsh_tem,1)&
          adjustl(trim(title)//'_'//intost(ielty)),ndime,&
          adjustl(trim(elemt)),nnode(ielty)
     !    
     ! Coordinates
     !
     if(ielty==1) then
        write(lun_ramsh_tem,2)'coordinates'           
        do ipoin=1,npoin
           if(lpoty(ipoin)>0) then
              write(lun_ramsh_tem,3) ipoin,(coord(idime,ipoin),idime=1,ndime)
           end if
        end do
        write(lun_ramsh_tem,2)'end coordinates'
     end if
     !         
     ! Connectivity
     !
     write(lun_ramsh_tem,2)'elements'
     do iboun=1,nboun
        if(ltypb(iboun)==ielty) then
           write(lun_ramsh_tem,4) iboun,&
                (lnodb(inodb,iboun),inodb=1,nnode(ielty)),ltypb(iboun)
        end if
     end do
     write(lun_ramsh_tem,2)'end elements'
  end do 
  call iofile(two,lun_ramsh_tem,dumml,'MESH')  ! GiD: close mesh file
!
! Results
!
  call iofile(zero,lun_rares_tem,fil_rares_tem,'POST-PROCESS')
  write(lun_rares_tem,'(a)')'GiD Post Results File 1.0'
  write(lun_rares_tem,'(a)')' '

  do ielty=1,nelty
     if(ndime==2) then
        elemt='Linear' 
     else
        if(nnode(ielty)==3) then
           elemt='Triangle' 
        else
           elemt='Quadrilateral'
        end if
     end if
     write(lun_rares_tem,5)&
          adjustl(trim(title))//'_Gauss_'//adjustl(trim(intost(ielty))),&
          adjustl(trim(elemt)),&
          adjustl(trim(title))//'_'//adjustl(trim(intost(ielty)))
     state='ANALYSIS'
     do iview=1,10
        iboun=kfl_viewf_tem(iview)
        if(iboun>=1) then
           wopos='VIEW_FACTOR_'//adjustl(trim(intost(kfl_viewf_tem(iview))))
           write(lun_rares_tem,72) &
                adjustl(trim(wopos)),state,0,'Scalar',&
                adjustl(trim(title))//'_Gauss_'//adjustl(trim(intost(ielty)))
           write(lun_rares_tem,71) 'Values'
           do jboun=1,nboun
              if(jboun>iboun) then
                 if(ltypb(iboun)==ielty) then
                    lboun=jboun-iboun
                    write(lun_rares_tem,74) jboun,viewf_tem(iboun)%a(lboun)
                 end if
              else if(jboun<iboun) then
                  if(ltypb(iboun)==ielty) then
                     lboun=iboun-jboun
                    write(lun_rares_tem,74) jboun,viewf_tem(jboun)%a(lboun)
                 end if
              else
                 if(ltypb(iboun)==ielty) then
                    write(lun_rares_tem,74) jboun,0.0_rp
                 end if
              end if
           end do
           write(lun_rares_tem,71) 'End Values'
           flush(lun_rares_tem)
        end if
     end do
  end do
   
1 format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2 format(a)
3 format(i7, 3(1x,e16.8e3))
4 format(i7,27(1x,i7))
5 format('GaussPoints ',a,' ElemType ',a,' ',a,/,&
       & 'Number Of Gauss Points: 1',/,&
       & 'Natural Coordinates: internal',/,&
       & 'End GaussPoints',/)
!
! GiD formats
!
71   format(a)
72   format('Result ',a,' ',a,' ',e12.6,' ',a,' OnGaussPoints ',a)
73   format('ComponentNames ',a)
74   format(i7, 3(1x,e16.8E3))
!
! Femview formats
!
710  format(1x,i4,a1,a6,e12.5,32x,i2,i5)
720  format(1x,i2,2x,a5,3x,3i5)
730  format(1x,i2,2x,a5,3x,2i5)
740  format(1x,i2,i5,e12.5)
750  format(1x,i2)
   
end subroutine tem_radpos

