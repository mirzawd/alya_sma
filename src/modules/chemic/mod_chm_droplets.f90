!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup chemic
!> @{
!> @file    chm_droplet_id.f90
!> @author  margarida moragues
!> @date    2020-03-30
!> @brief   Eulerian droplets
!> @details Toolbox for Eulerian droplets identification
!-----------------------------------------------------------------------

module mod_chm_droplets

  use def_kintyp,    only : ip, rp
  use def_master,    only : INOTSLAVE, INOTEMPTY, cutim, itcou, ittim, modul, itinn, mem_modul, modul, conce
  use def_chemic,    only : diameter_drop_chm, lun_droplet_chm, ndrop_chm, centroid_drop_chm, droplet_compactness_limit_chm,&
                            droplet_max_diameter_chm, levelSet_threshold_chm, ndrop_chm, volume_cluster_chm, volume_drop_chm,&
                            compactness2_drop_chm, diameter_drop_chm

  implicit none

  private

  public :: chm_droplet_id
  public :: chm_droplet_output

contains

!------------------------------------------------------------------------
!>
!> @author  margarida moragues
!> @date    2020-03-30
!> @brief   Identification of Eulerian droplets on a level set field
!> @details Computation of its volume, diameter, centroid, compactness, etc.
!>
!------------------------------------------------------------------------

  subroutine chm_droplet_id()

    use def_domain,            only : ndime, nelem, meshe
    use def_parame,            only : pi
    use mod_clusters,          only : clusters, mask_from_field, clusters_volume
    use def_kermod,            only : ndivi
    use mod_memory,            only : memory_alloca, memory_deallo, memory_resize
    use mod_shape_descriptors, only : compute_compactness_2

    implicit none
    integer(ip),       pointer     :: lmask(:),legro(:)
    integer(ip)                    :: idime,iclus,nclus
    real(rp),          pointer     :: xx(:)

    external                       :: runend

    nullify(lmask)
    nullify(legro)
    nullify(xx)

    call memory_alloca(mem_modul(1:2,modul),'LMASK','chm_droplet_id',lmask,max(1_ip,nelem))

    if( INOTEMPTY ) xx => conce(:,3,1)

    call mask_from_field(lmask,xx,levelSet_threshold_chm,ON_ELEMENTS=.true.)

    call clusters(lmask,legro,nclus,meshe(ndivi),ON_ELEMENTS=.true.)

    call memory_deallo(mem_modul(1:2,modul),'VOLUME_DROP_CHM',      'chm_droplet_id',volume_drop_chm)
    call memory_deallo(mem_modul(1:2,modul),'DIAMETER_DROP_CHM',    'chm_droplet_id',diameter_drop_chm)
    call memory_deallo(mem_modul(1:2,modul),'CENTROID_DROP_CHM',    'chm_droplet_id',centroid_drop_chm)
    call memory_deallo(mem_modul(1:2,modul),'COMPACTNESS2_DROP_CHM','chm_droplet_id',compactness2_drop_chm)
    call memory_deallo(mem_modul(1:2,modul),'VOLUME_CLUSTER_CHM',   'chm_droplet_id',volume_cluster_chm)

    call memory_alloca(mem_modul(1:2,modul),'VOLUME_DROP_CHM',      'chm_droplet_id',volume_drop_chm,      nclus)
    call memory_alloca(mem_modul(1:2,modul),'DIAMETER_DROP_CHM',    'chm_droplet_id',diameter_drop_chm,    nclus)
    call memory_alloca(mem_modul(1:2,modul),'CENTROID_DROP_CHM',    'chm_droplet_id',centroid_drop_chm,    ndime,nclus)
    call memory_alloca(mem_modul(1:2,modul),'COMPACTNESS2_DROP_CHM','chm_droplet_id',compactness2_drop_chm,nclus)
    call memory_alloca(mem_modul(1:2,modul),'VOLUME_CLUSTER_CHM',   'chm_droplet_id',volume_cluster_chm,   nclus)

    call clusters_volume(legro,nclus,volume_cluster_chm)

    call compute_compactness_2(legro,nclus,xx,volume_drop_chm,centroid_drop_chm,compactness2_drop_chm)

    call memory_deallo(mem_modul(1:2,modul),'LMASK','chm_droplet_id',lmask)
    call memory_deallo(mem_modul(1:2,modul),'LEGRO','chm_droplet_id',legro)

    ndrop_chm = 0
    do iclus=1,nclus
       if ( ndime == 2 ) then
          diameter_drop_chm(iclus) = sqrt(4.0_rp * volume_drop_chm(iclus) / pi)
       else if (ndime == 3) then
          diameter_drop_chm(iclus) = (6.0_rp * volume_drop_chm(iclus) / pi) ** (1.0_rp/3.0_rp)
       else
        call runend('CHM_INIUNK: WRONG DIMENSION NDIME')
       end if

       if ( compactness2_drop_chm(iclus) >= droplet_compactness_limit_chm&
           .and. diameter_drop_chm(iclus) <= droplet_max_diameter_chm ) then
          ndrop_chm = ndrop_chm + 1
          volume_drop_chm(ndrop_chm) = volume_drop_chm(iclus)
          diameter_drop_chm(ndrop_chm) = diameter_drop_chm(iclus)
          volume_cluster_chm(ndrop_chm) = volume_cluster_chm(iclus)
          do idime = 1, ndime
             centroid_drop_chm(idime,ndrop_chm) = centroid_drop_chm(idime,iclus)
          end do
          compactness2_drop_chm(ndrop_chm) = compactness2_drop_chm(iclus)
       end if
    end do

    call memory_resize(mem_modul(1:2,modul),'VOLUME_DROP_CHM',       'chm_droplet_id',volume_drop_chm,      ndrop_chm)
    call memory_resize(mem_modul(1:2,modul),'DIAMETER_DROP_CHM',     'chm_droplet_id',diameter_drop_chm,    ndrop_chm)
    call memory_resize(mem_modul(1:2,modul),'CENTROID_DROP_CHM',     'chm_droplet_id',centroid_drop_chm,    ndime,ndrop_chm)
    call memory_resize(mem_modul(1:2,modul),'COMPACTNESS2_DROP_CHM', 'chm_droplet_id',compactness2_drop_chm,ndrop_chm)
    call memory_resize(mem_modul(1:2,modul),'VOLUME_CLUSTER_CHM',    'chm_droplet_id',volume_cluster_chm,   ndrop_chm)

  end subroutine chm_droplet_id



!------------------------------------------------------------------------
!>
!> @author  margarida moragues
!> @date    2020-04-20
!> @brief   Eulerian droplets output results
!> @details For a given time step, it writes on a file the list of
!>          droplets and its volume, diameter, centroid, compactness, etc.
!>
!------------------------------------------------------------------------

  subroutine chm_droplet_output()

    use def_domain,      only :  ndime

    implicit none
    integer(ip)                 :: idime,idrop
    integer(ip), save           :: ipass=0

    if( INOTSLAVE ) then

       if( ipass == 0 ) write(lun_droplet_chm,100)
       ipass = ipass + 1
       write(lun_droplet_chm,110) ittim,itcou,itinn(modul),cutim,ndrop_chm
       write(lun_droplet_chm,120)
       do idrop=1,ndrop_chm
          write(lun_droplet_chm,130) idrop, volume_cluster_chm(idrop), volume_drop_chm(idrop), compactness2_drop_chm(idrop),&
              diameter_drop_chm(idrop), (centroid_drop_chm(idime,idrop), idime=1,ndime)
       end do
       write(lun_droplet_chm,*) ' '
    end if

    !
    ! Formats
    !
  100 format('# --| ALYA Eulerian Droplet Data                ',/,&
         &   '# --| Columns displayed:                        ',/,&
         &   '# --| 1. Id number                              ',/,&
         &   '# --| 2. Cluster Volume                         ',/,&
         &   '# --| 3. Spherical droplet Volume (*)           ',/,&
         &   '# --| 4. Compactness_2 (*)                      ',/,&
         &   '# --| 5. Spherical droplet Diameter (*)         ',/,&
         &   '# --| 6. Centroid_x (*)                         ',/,&
         &   '# --| 7. Centroid_y (*)                         ',/,&
         &   '# --| 8. Centroid_z (*)                         ',/,&
         &   '# --| (*) it takes into account the fluid concentration (level set function) ',/,&
         &   ' ')
  110 format('# Time step =', i9, '   Global Iteration =', i9, '   Inner Iteration =', i9, '   Current time =', e12.5, '   Total&
      & number of droplets =', i9)

  120 format('#           1','                 2','                 3','                 4',&
         &   '                 5','                 6','                 7','                 8')
  130 format(4x,i9,100(2x,e16.8e3))

  end subroutine chm_droplet_output

end module mod_chm_droplets
!> @}

