!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine poscod(itask)
  !------------------------------------------------------------------------
  !****f* domain/poscod
  ! NAME 
  !    poscod
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    reacod
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ifint,krang(50,5),irang,kfixi,kfibo,ifoun
  integer(ip)             :: kcode,nrang,idime,kcods(10),nfoun,ifaul
  integer(ip)             :: kfoun,kfaul,ipoin,iboun,pnodb,pblty,inodb
  real(rp)                :: vrang(50,6,3)

  !----------------------------------------------------------------------
  !
  ! Reading the conditions
  !
  !----------------------------------------------------------------------

  ifint= 0
  ifaul= 0
  kfaul= 0
  krang= 0
  vrang= 0.0_rp
  irang= 0

  do while( words(1)/='ENDON' )

     if( words(1) =='INTER' ) then

        if(exists('ADDIN')) ifint = 1 

     else if( words(1) =='DEFAU' ) then

        ifaul = 1
        kfaul = int(param(1))

     else if( words(1) =='BCBOU' ) then      

        irang = irang + 1
        do while(words(1)/='ENDBC')
           !
           ! boundary conditions inside a bounding box
           !
           if( words(1) =='SMALL' ) then 
              vrang(irang,1,1) = param(1)
              vrang(irang,1,2) = param(2)
              vrang(irang,1,3) = param(3)
           else if(words(1)=='LARGE' ) then                             
              vrang(irang,2,1) = param(1)
              vrang(irang,2,2) = param(2)
              vrang(irang,2,3) = param(3)
           else if(words(1)=='CODEN' ) then                             
              krang(irang,1) = 1
              krang(irang,2) = int(param(1))
           end if

           call ecoute('poscod')
        end do

     else if( words(1) =='RANGE' ) then         

        irang = irang + 1

        do while(words(1)/='ENDRA')
           !
           ! limiting hyper-surface
           !
           if( words(1) =='CODEN' ) then                             
              kcode = getint('CODEN',1_ip,'#CODE NUMBER FOR THIS BODY')
              krang(irang,2) = kcode
           else if( words(1) =='COORD' ) then
              if(words(2)=='XCOOR') krang(irang,3) = 1
              if(words(2)=='YCOOR') krang(irang,3) = 2
              if(words(2)=='ZCOOR') krang(irang,3) = 3
           else if(words(1)=='LOWER' ) then                             
              krang(irang,1)   =  10
              vrang(irang,1,1) = param(1)
           else if(words(1)=='GREAT' ) then
              ! krang(irang,1)   =  krang(irang,1) + 1
              krang(irang,1)   =  20
              vrang(irang,1,1) =  param(1)
           end if

           call ecoute('poscod')
        end do

     end if

     call ecoute('poscod')

  end do

  nrang = irang

  !----------------------------------------------------------------------
  !
  ! Conditions on nodes
  !
  !----------------------------------------------------------------------

  if( itask == 1  ) then

     points: do ipoin=1,npoin

        kfixi = -1
        kfibo =  0
        ifoun =  0

        nodal_conditions: do irang= 1,nrang

           kfoun= 0
           if(krang(irang,1) == 1 ) then   ! bounding box
              if(coord(1,ipoin) > vrang(irang,1,1) ) then
                 if(coord(1,ipoin) < vrang(irang,2,1) ) then
                    kfoun= kfoun+1
                 end if
              end if
              if((kfoun == 1) .and. (ndime>=2) ) then
                 if(coord(2,ipoin) > vrang(irang,1,2) ) then
                    if(coord(2,ipoin) < vrang(irang,2,2) ) then
                       kfoun= kfoun+1
                    end if
                 end if
              end if
              if((kfoun == 2) .and. (ndime==3) ) then
                 if(coord(3,ipoin) > vrang(irang,1,3) ) then
                    if(coord(3,ipoin) < vrang(irang,2,3) ) then
                       kfoun= kfoun+1
                    end if
                 end if
              end if
              if(kfoun == ndime ) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)                 
              end if
           else if(krang(irang,1) == 10  ) then   ! lower than the coordinate value
              idime= krang(irang,3)
              if(coord(idime,ipoin) < vrang(irang,1,1) ) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)
              end if
           else if(krang(irang,1) == 20  ) then   ! greater than the coordinate value
              idime= krang(irang,3)
              if(coord(idime,ipoin) > vrang(irang,1,1) ) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)
              end if

           end if

        end do nodal_conditions
 
        nfoun = ifoun
        if( ( ifint == 0 ) .and. ( nfoun > 1 ) ) nfoun = 1  !Truncates to the first condition found

        if( ( nfoun == 0 ) .and. ( ifaul == 1 ) ) then     ! a default code is given
           nfoun = 1      !Make as if one code was found
           kcods(nfoun) = kfaul
        end if

        do ifoun = 1,min(nfoun,mcono) 
           kfl_codno(ifoun,ipoin) = kcods(ifoun)
        end do

     end do points


  !----------------------------------------------------------------------
  !
  ! Conditions on boundaries
  !
  !----------------------------------------------------------------------

  else if( itask == 2  ) then

     boundaries: do iboun=1,nboun

        kfixi = -1
        kfibo =  0
        ifoun =  0
        pblty = abs(ltypb(iboun))
        pnodb = nnode(pblty)

        bound_conditions: do irang= 1,nrang

           kfoun= 0
           if(krang(irang,1) == 1 ) then   ! bounding box
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 if(coord(1,ipoin) > vrang(irang,1,1) ) then
                    if(coord(1,ipoin) < vrang(irang,2,1) ) then
                       kfoun= kfoun+1
                    end if
                 end if
                 if((kfoun == 1) .and. (ndime>=2) ) then
                    if(coord(2,ipoin) > vrang(irang,1,2) ) then
                       if(coord(2,ipoin) < vrang(irang,2,2) ) then
                          kfoun= kfoun+1
                       end if
                    end if
                 end if
                 if((kfoun == 2) .and. (ndime==3) ) then
                    if(coord(3,ipoin) > vrang(irang,1,3) ) then
                       if(coord(3,ipoin) < vrang(irang,2,3) ) then
                          kfoun= kfoun+1
                       end if
                    end if
                 end if
              end do
              if( kfoun == ndime*pnodb ) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)                 
              end if
           else if(krang(irang,1) == 10  ) then   ! lower than the coordinate value
              idime= krang(irang,3)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun) 
                 if(coord(idime,ipoin) < vrang(irang,1,1) ) then
                    kfoun= kfoun + 1
                 end if
              end do
              if( kfoun == pnodb ) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)                 
              end if
           else if(krang(irang,1) == 20  ) then   ! greater than the coordinate value
              idime= krang(irang,3)
              do inodb = 1,pnodb 
                 ipoin = lnodb(inodb,iboun)
                 if(coord(idime,ipoin) > vrang(irang,1,1) ) then
                    kfoun= kfoun + 1
                 end if
              end do
              if( kfoun == pnodb ) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)                 
              end if
           end if

        end do bound_conditions
 
        nfoun = ifoun
        if( ( ifint == 0 ) .and. ( nfoun > 1 ) ) nfoun = 1  ! Truncates to the first condition found

        if( ( nfoun == 0 ) .and. ( ifaul == 1 ) ) then      ! a default code is given
           nfoun = 1      !Make as if one code was found
           kcods(nfoun) = kfaul
        end if

        do ifoun = 1,min(nfoun,mcono) 
           kfl_codbo(iboun) = kcods(ifoun)  ! IT KEEPS THE LAST CONDITION FOUND, UNLIKE NODES BOUNDARIES CAN ONLY HAVE 1 COND 

        end do

     end do boundaries

  end if

end subroutine poscod
