!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_livinf(itask,message,inume)
!-----------------------------------------------------------------------
!****f* outrut/par_livinf
! NAME
!    par_livinf
! DESCRIPTION
!    This routine write live information on run.
! USES
! USED BY
!    Reapro
!***
!-----------------------------------------------------------------------    
  use def_master
  use def_domain
  use def_parall
  use mod_parall
  use mod_messages, only : livinf
  use mod_messages, only : messages_live
  use mod_live_info_config, only : live_info_config
  implicit none
  integer(ip),   intent(in) :: itask,inume
  character(*)              :: message
  character(300)            :: messa,messb
  character(20)             :: mess1,mess2

  if( INOTSLAVE .and. live_info_config%lun_livei /= 0 ) then     

     if(itask==0) then
        if(inume==nproc_par-1) then
           messa='ALL SLAVES HAVE INITIATED'
        else
           mess1=intost(nproc_par-inume-1_ip)
           call runend(trim(mess1)//'PROCESSES HAVE NOT BEEN INITIATED')
        end if

     else if(itask==3) then
        if( kfl_partition_par == PAR_METIS4 ) then
           if(kfl_parti_par==1) then
              messa='MASTER PARTITIONS ELEMENT GRAPH WITH METIS USING NODES CONNECTIVITY'
           else
              messa='MASTER PARTITIONS ELEMENT GRAPH WITH METIS USING FACE CONNECTIVITY'
           end if
        else if( kfl_partition_par == PAR_SFC ) then
           messa='MASTER PARTITIONS MESH USING A HILBERT SPACE FILLING CURVE'
        else if( kfl_partition_par == PAR_ZOLTAN ) then
           messa='MASTER PARTITIONS MESH USING ZOLTAN A HSFC'
        else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
           messa='MASTER PARTITIONS MESH USING AN ORIENTED BIN'
        else if( kfl_partition_par < 0 ) then
           messa='MASTER PARTITIONS MESH USING A FIELD'
        end if

     else if(itask==4) then
        messa='MASTER COMPUTES COMMUNICATION STRATEGY'

     else if(itask==5) then
        messa='MASTER COMPUTES PERMUTATION ARRAYS'

     else if(itask==6) then
        if(kfl_ptask==0) then
           messa='MASTER WRITES MESH AND PARTITION DATA IN RESTART FILES'
        else if(kfl_ptask==1) then
           messa='MASTER/SLAVES EXCHANGE MESH AND PARTITION DATA'
        else if(kfl_ptask==2) then
           messa='MASTER/SLAVES READ MESH AND PARTITION DATA FROM RESTART FILE'
        end if
        if(kfl_fileh_par==1.and.kfl_ptask/=1) messa=trim(messa)//' WITH HIERARCHY'

     else if(itask==7) then
        continue

     else if(itask==8) then
        messa='MASTER/SLAVES READ MESH AND PARTITION DATA FROM RESTART FILE'

     else if(itask==9) then
        messa='MASTER OUTPUT MESH'

     else if(itask==10) then
        if(inume/=0) then
           call runend('SOME SLAVES COULD NOT READ THEIR RESTART FILES')
        else
           messa='MASTER/SLAVES: ALL SLAVES HAVE READ THEIR RESTART FILES'
        end if

     else if(itask==11) then
        mess1=intost(nproc_par-1_ip)
        messa='CHECK MPI. 1 MASTER + '//trim(mess1)//' SUBDOMAINS'

     else if(itask==12) then
        messa='MPI IS WORKING WELL'

     else if(itask==13) then
        messa='MPI IS WORKING WELL'

     else if(itask==14) then
        messa=trim(namod(modul))//': MASTER SENDS TO SLAVES GROUPS OF DEFLATED CG'

     else if(itask==15) then
        if(kfl_ptask==0) then
           messa=trim(namod(modul))//': MASTER WRITES BOUNDARY CONDITIONS IN RESTART FILES'
        else if(kfl_ptask==1) then
           messa=trim(namod(modul))//': MASTER SENDS BOUNDARY CONDITIONS TO SLAVES'
        else if(kfl_ptask==2) then
           messa=trim(namod(modul))//': SLAVES READ BOUNDARY CONDITIONS FROM RESTART FILES'
        end if

     else if(itask==16) then
        if(kfl_ptask==0) then
           messa=trim(namod(modul))//': MASTER WRITES PHYSICAL AND NUMERICAL DATA IN RESTART FILES'
        else if(kfl_ptask==1) then
           messa=trim(namod(modul))//': MASTER SENDS PHYSICAL AND NUMERICAL DATA TO SLAVES'
        else if(kfl_ptask==2) then
           messa=trim(namod(modul))//': SLAVES READ PHYSICAL AND NUMERICAL DATA FROM RESTART FILES'
        end if

     else if(itask==17) then
        if(kfl_ptask==0) then
           messa=trim(namod(modul))//': MASTER WRITES PHYSICAL AND NUMERICAL ARRAYS IN RESTART FILES'
        else if(kfl_ptask==1) then
           messa=trim(namod(modul))//': MASTER SENDS PHYSICAL AND NUMERICAL ARRAYS TO SLAVES'
        else if(kfl_ptask==2) then
           messa=trim(namod(modul))//': SLAVES READ PHYSICAL AND NUMERICAL ARRAYS FROM RESTART FILES'
        end if

     else if(itask==18) then
        messa='MASTER ORDERS INTERIOR AND BOUNDARY NODES'

     else if (itask == 20) then
        messa='PARALL PREPROCESS -->  '//trim(message)

     else if (itask == 21) then
        messa='SLAVE START PARTITION PREPROCESS -->  '//trim(message)

     else if (itask == 22) then
        messa='SLAVE ENDS PARTITION PREPROCESS -->  '//trim(message)

     end if
     
     if(itask==2) then
        
        messa = 'MASTER COMPUTES ELEMENT GRAPH'
        call messages_live(trim(messa))
        
     else if(itask==1) then
        mess1 = intost(npart_par)
        messb = 'MESH PARTITION (# OF SUBDOMAINS= '//trim(mess1)//')'
        call messages_live(trim(messb),'START SECTION')

     else if(itask==7) then
        mess1 = intost(nedge)
        mess2 = intost(medge)
        messb = 'MASTER FOUND '//trim(mess1)//' EDGES AND '//trim(mess2)//' MAX # EDGES/ELEMENT'
        call messages_live(trim(messb))
        
     else if(itask==1000) then
        
        call messages_live('MESH PARTITION','END SECTION')
        
     else
        
        call messages_live(trim(messa))

     end if

  end if

  !if( inews == 1 ) isect = isect + 3

1 format(' (# EDGES= ',a,', MAX # EDGES/ELEMENT= ',a,')')

end subroutine par_livinf
