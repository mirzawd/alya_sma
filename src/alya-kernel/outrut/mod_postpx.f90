!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  module mod_postpx
  !-----------------------------------------------------------------------
  !****f* outrut/mod_postpx
  ! NAME
  !   mod_postpr
  ! DESCRIPTION
  !   This routine manages the postprocess
  ! USES
  ! USED BY
  !   output_***
  !   outvar_***
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_postpr
  
  interface postpx
     module procedure &
          posscx,posvex
  end interface
  
contains
  
  subroutine posscx(bridge,wopos,itste,ttime,kpoin)

    !-----------------------------------------------------------------------
    !
    ! Write a scalar in postprocess file
    !
    !-----------------------------------------------------------------------
    implicit none
    character(*), intent(in)         :: wopos(*)
    complex(rp),  intent(in), target :: bridge(:)
    integer(ip) , intent(in)         :: itste
    real(rp)    , intent(in)         :: ttime
    integer(ip),  optional           :: kpoin
    integer(ip)                      :: ipoin,iallo,ii,ierror
    character(8)                     :: state
    character(25)                    :: filename
    character(5)                     :: wopex(2)
        
    if( IMASTER .and. kfl_postp_par == 0 ) then
       return
    else
       wopex(1) = '_REAL'
       wopex(2) = '_IMAG'
       !
       ! If postprocess is carried out by Master, scatter vector
       !
       iallo = 0
       if( ISEQUEN .or. ( ISLAVE .and. kfl_postp_par == 0 ) ) then
         parx1 => bridge
       else if( kfl_postp_par == 1 ) then
         if( INOTMASTER ) then         
           parx1 => bridge
         else           
           iallo =  1
         end if       
           party =  3
           pardi =  1
           parki =  4
           pard1 =  ndime
           call par_scatte()
           if( INOTMASTER ) return      
       end if
    
       state='ANALYSIS'
       !
       ! GiD format
       !
       filename = 'potential_scalar.out'
       open(unit=1,file=filename,iostat=ierror)
       if (ierror == 0_ip) then
         write(1,'(a)')'GiD Post Results File 1.0'
         write(1,'(a)')' '              
         do ii = 1,2
           write(1,2) wopos(1)//wopex(ii),state,ttime,'Scalar'
           write(1,3) wopos(1)//wopex(ii)
           write(1,1) 'Values'              
           if( ii == 1 ) then
             do ipoin=1,npoin
               write(1,4)ipoin,real(parx1(ipoin))
             enddo
           else
             do ipoin=1,npoin
               write(1,4)ipoin,aimag(parx1(ipoin))
             enddo
           endif
           write(1,1) 'End Values'
           !flush(1)
         enddo
         close (unit=1)
       else
  	     write(*,*) 'ERROR OPENING FILE: ',filename
       endif         
       !
       ! Parall: deallocate memory if necessary
       !
    end if    
    !
    ! GiD formats
    !
1   format(a)
2   format('Result ',a,' ',a,' ',e14.8,' ',a,' OnNodes')
3   format('ComponentNames ',a)
4   format(i9, 3(1x,e16.8E3))

  end subroutine posscx

  subroutine posvex(bridge,wopos,itste,ttime,kpoin)

    !-----------------------------------------------------------------------
    !
    ! Write a vector in postprocess file
    !
    !-----------------------------------------------------------------------
    implicit none
    character(*), intent(in)         :: wopos(*)
    complex(rp),  intent(in), target :: bridge(:,:)
    integer(ip),  intent(in)         :: itste
    real(rp),     intent(in)         :: ttime
    integer(ip),  optional           :: kpoin
    integer(ip)                      :: ipoin,idime,ierror
    integer(ip)                      :: iallo,ii
    character(8)                     :: state
    character(20)                    :: wopo2(3)
    character(25)                    :: filename
    character(5)                     :: wopex(2)

    if( IMASTER .and. kfl_postp_par == 0 ) then
       return
    else   
       wopex(1) = '_REAL'
       wopex(2) = '_IMAG'
       !
       ! If postprocess is carried out by Master, scatter vector
       !
       iallo = 0
       if( ISEQUEN .or. ( ISLAVE .and. kfl_postp_par == 0 ) ) then
         parx2 => bridge
       else if( kfl_postp_par == 1 ) then
         if( INOTMASTER ) then         
           parx2 => bridge
         else           
           iallo =  1
         end if       
           party =  3
           pardi =  2
           parki =  4
           pard1 =  ndime
           call par_scatte()
           if( INOTMASTER ) return                 
       end if
      
       wopo2(1) = trim(wopos(1)) // '_X'
       wopo2(2) = trim(wopos(1)) // '_Y'
       wopo2(3) = trim(wopos(1)) // '_Z'
       state='ANALYSIS'
       !
       ! GiD format
       ! 
       filename = 'potential_vector.out'
       open(unit=2,file=filename,iostat=ierror)
       if (ierror == 0_ip) then
         write(2,'(a)')'GiD Post Results File 1.0'
         write(2,'(a)')' '       
         do ii = 1,2     
           write(2,2) wopos(1)//wopex(ii),state,ttime,'Vector'
           write(2,3) trim(wopo2(1))//wopex(ii), trim(wopo2(2))//wopex(ii), trim(wopo2(3))//wopex(ii)
           write(2,1) 'Values'                
           if( ii == 1 ) then 
             do ipoin=1,npoin
               write(2,4)ipoin,(real(parx2(idime,ipoin)),idime=1,ndime)
             enddo
           else
             do ipoin=1,npoin
               write(2,4)ipoin,(aimag(parx2(idime,ipoin)),idime=1,ndime)
             enddo
           endif
           write(2,1) 'End Values'
           !flush(2)
         enddo
         close (unit=2)
       else
  	     write(*,*) 'ERROR OPENING FILE: ',filename
       endif                      
       !
       ! Parall: deallocate memory if necessary
       !
    end if     
    !
    ! GiD formats
    !
1   format(a)
2   format('Result ',a,' ',a,' ',e14.8,' ',a,' OnNodes')
3   format('ComponentNames ',a,',',a,',',a)
4   format(i9, 3(1x,e16.8E3))

  end subroutine posvex
  
end module mod_postpx
