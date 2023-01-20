!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_enrich(itask)
  !-----------------------------------------------------------------------
  !****f* solidz/sld_enrich
  ! NAME
  !    sld_enrich
  ! DESCRIPTION
  !    This routines computes things related to enrichement
  ! OUTPUT
  !    VMASS_SLD(NPOIN) : Diagonal mass matrix
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_solidz
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ielem,idofn,idime,inode
  
  logical                 :: debugging
  debugging = .false.

  if( kfl_xfeme_sld == 0 ) return

  select case ( itask )

  case ( 1_ip )
     !
     ! Identify enriched nodes => LNENR_SLD(IPOIN) = 1
     !
     if( INOTMASTER ) then

        do ipoin = 1,npoin
           lnenr_sld(ipoin) = 0
        end do

        do ielem = 1,nelem
           if (leenr_sld(ielem) /= 0) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 lnenr_sld(ipoin) = leenr_sld(ielem)
              end do
           end if
        end do

        call parari('SLX',NPOIN_TYPE,npoin,lnenr_sld)

        do ipoin = 1,npoin
           lnenr_sld(ipoin) = min(lnenr_sld(ipoin),1_ip)
        end do
        
        if (debugging) then
           write(*,*)'ipon, lnenr_sld = '
           do ipoin = 1,npoin
              if (lnenr_sld(ipoin) > 0) then
                 write(*,*)ipoin,lnenr_sld(ipoin)
              end if
           enddo
        end if

     end if


  case ( 2_ip )
     !
     ! Put VMASS = 1 and RHSID = 0 for non-enriched nodes
     !
     if( INOTMASTER ) then

        do ipoin = 1,npoin
           if( lnenr_sld(ipoin) == 0 ) then
              idofn = ( ipoin - 1 ) * ndofn_sld + ndime
              vmasx_sld(ipoin) = 1.0_rp
              do idime = 1,ndime
                 idofn = idofn + 1
                 rhsid(idofn) = 0.0_rp
              end do
           end if
        end do

     end if

  end select

end subroutine sld_enrich
