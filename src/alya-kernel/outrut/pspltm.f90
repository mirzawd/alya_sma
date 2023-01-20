!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pspltm(&
     nrow2,ncol2,ndof,mode,ja,ia,amatr,title,ptitle,dsize,&
     munt,nlines,lines,izero,iunt)
  !-----------------------------------------------------------------------
  ! PSPLTM - PostScript PLoTer of a (sparse) Matrix
  ! This version by loris renggli (renggli@masg1.epfl.ch), Dec 1991
  ! and Youcef Saad 
  !------
  ! Loris RENGGLI, Swiss Federal Institute of Technology, Math. Dept
  ! CH-1015 Lausanne (Switzerland)  -- e-mail:  renggli@masg1.epfl.ch
  ! Modified by Youcef Saad -- June 24, 1992 to add a few features:
  ! separation lines + acceptance of MSR format.
  !-----------------------------------------------------------------------
  ! input arguments description :
  ! 
  ! nrow   = number of rows in matrix
  ! 
  ! ncol   = number of columns in matrix 
  ! 
  ! izero  = 0,1,2: null coef. do not appear/always appear/are clearer
  !
  ! mode   = integer indicating whether the matrix is stored in 
  ! CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2) 
  ! 
  ! ja     = column indices of nonzero elements when matrix is
  ! stored rowise. Row indices if stores column-wise.
  ! ia     = integer array of containing the pointers to the 
  ! beginning of the columns in arrays a, ja.
  ! 
  ! title  = character*(*). a title of arbitrary length to be printed 
  ! as a caption to the figure. Can be a blank character if no
  ! caption is desired.
  ! 
  ! ptitle = position of title; 0 under the drawing, else above
  ! 
  ! dsize   = dsize of the drawing  
  ! 
  ! munt   = units used for dsize : 'cm' or 'in'
  ! 
  ! nlines = number of separation lines to draw for showing a
  ! partionning
  ! of the matrix. enter zero if no partition lines are wanted.
  ! 
  ! lines  = integer array of length nlines containing the coordinates
  ! of 
  ! the desired partition lines . The partitioning is symmetric: 
  !      a horizontal line across the matrix will be drawn in 
  !      between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
  !      an a vertical line will be similarly drawn between columns
  !      lines(i) and lines(i)+1 for i=1,2,...,nlines 
  !
  ! iunt   = logical unit number where to write the matrix into.
  !----------------------------------------------------------------------- 
  ! additional note: use of 'cm' assumes european format for paper dsize
  ! (21cm wide) and use of 'in' assumes american format (8.5in wide).
  ! The correct centering of the figure depends on the proper choice. Y.S.
  !-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip),   intent(in) :: nrow2,ncol2,ndof,ptitle,mode,iunt,nlines
  integer(ip),   intent(in) :: lines(*),ja(*),ia(nrow2+1)
  integer(ip),   intent(in) :: izero
  real(rp),      intent(in) :: amatr(ndof,ndof,*),dsize
  character*(*), intent(in) :: title*(*)
  character(2),  intent(in) :: munt
  character(50)             :: messa
  integer(ip)               :: n,nr,nc,maxdim,istart,ilast,nrow,ncol
  integer(ip)               :: ii,k,ltit,m,kol,isep,idof,jdof
  real(rp)                  :: lrmrgn,botmrgn,xtit,ytit,ytitof
  real(rp)                  :: fnstit,siz,xl,xr,yb,yt,scfct,u2dot
  real(rp)                  :: frlw,delt,paperx,xx,yy,sif,xf,yf
  real(rp)                  :: haf=0.5_rp,zero=0.0_rp,conv=2.54_rp
  !
  ! square = .true. to draw a square frame around rectangular matrix
  ! number = .true. to put node number on diagonal
  ! box    = .true. to put a box around the drawing
  ! tit    = .true. to write title
  !
  logical(lg)               :: square = .false.
  logical(lg)               :: number = .false.
  logical(lg)               :: box    = .true.
  logical(lg)               :: tit    = .true.

  nrow = nrow2*ndof
  ncol = ncol2*ndof

  siz = dsize
  nr  = nrow
  nc  = ncol
  n   = nc
  if(mode==0) n = nr
  maxdim = max(nrow, ncol)
  m   = 1 + maxdim
  nc  = nc+1
  nr  = nr+1
  sif = 0.50_rp*28.0_rp/real(nrow)
  xf  = 0.0_rp
  yf  = 0.15_rp*28.0_rp/real(nrow)
  !
  ! units (cm or in) to dot conversion factor and paper size
  ! 
  if (munt=='cm'.or.munt=='CM') then
     u2dot  = 72.0_rp/conv
     paperx = 21.0_rp
  else
     u2dot  = 72.0_rp
     paperx = 8.5_rp*conv
     siz    = siz*conv
  end if
  !
  ! left and right margins (drawing is centered)
  ! 
  lrmrgn = (paperx-siz)/2.0_rp
  !
  ! bottom margin : 2 cm
  !
  botmrgn = 2.0_rp
  !
  ! scaling factor
  !
  scfct = siz*u2dot/m
  !
  ! matrix frame line witdh
  !
  frlw = 0.25_rp
  !
  ! font size for title (cm)
  !
  fnstit = 0.5_rp
  ltit = len(trim(title))
  !
  ! position of title : centered horizontally
  ! at 1.0 cm vertically over the drawing
  !
  ytitof = 1.0_rp
  xtit   = paperx/2.0_rp
  ytit   = botmrgn+siz*nr/m + ytitof
  !
  ! almost exact bounding box
  !
  xl = lrmrgn*u2dot - scfct*frlw/2
  xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
  yb = botmrgn*u2dot - scfct*frlw/2
  yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2
  if(ltit>0) then
     yt = yt + (ytitof+fnstit*0.70)*u2dot
  end if
  !
  ! add some room to bounding box
  !
  delt = 10.0_rp
  xl   = xl-delt
  xr   = xr+delt
  yb   = yb-delt
  yt   = yt+delt
  !
  ! correction for title under the drawing
  !
  if(ptitle==0.and.ltit>0) then
     ytit = botmrgn + fnstit*0.3_rp
     botmrgn = botmrgn + ytitof + fnstit*0.7_rp
  end if
  !
  ! begin of output
  !
  write(iunt,10) '%!'
  write(iunt,10) '%%Creator: PSPLTM routine'
  write(iunt,12) '%%BoundingBox:',xl,yb,xr,yt
  write(iunt,10) '%%EndComments'
  write(iunt,10) '/cm {72 mul 2.54 div} def'
  write(iunt,10) '/mc {72 div 2.54 mul} def'
  write(iunt,10) '/pnum { 72 div 2.54 mul 20 string'
  write(iunt,10) 'cvs print ( ) print} def'
  write(iunt,10) '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
  !
  ! Define colors:
  !
  write(iunt,10) '%% ----------------------------------------------------- COLORS '
  write(iunt,10) '%% To change color, modify the following constants:'
  write(iunt,10) '%% cnum: numbering of diagonal'
  write(iunt,10) '%% czer: zero coefficients'
  write(iunt,10) '%% cnze: non-zero coefficients'
  write(iunt,10) '%% ctex: text and border' 
  write(iunt,10) '/cnum { 1.0 0.0 0.0 setrgbcolor } bind def' ! numbering of diagonal
  !write(iunt,10) '/czer { 0.7 0.7 0.7 setrgbcolor } bind def' ! zero coefficients
  write(iunt,10) '/czer { 0.0 0.0 0.0 setrgbcolor } bind def' ! zero coefficients
  write(iunt,10) '/cnze { 0.0 0.0 0.0 setrgbcolor } bind def' ! non-zero coefficients
  write(iunt,10) '/ctex { 0.0 0.0 0.0 setrgbcolor } bind def' ! text and border
  write(iunt,10) '%% ----------------------------------------------------- END COLORS '
  !
  ! Define font
  !
  write(iunt,10) '%% ----------------------------------------------------- FONT '
  write(iunt,10) '%% To change font of diagonal number, modify the following:'
  write(iunt,10) '/fnum {'
  !write(iunt,10) '/Helvetica findfont'
  !write(iunt,15) sif,' scalefont'
  !write(iunt,10) 'setfont'
  write(iunt,10) '      } bind def'
  write(iunt,10) '%% ----------------------------------------------------- END FONT '
  !
  ! we leave margins etc. in cm so it is easy to modify them if
  ! needed by editing the output file
  !
  write(iunt,10) 'ctex'
  write(iunt,10) 'gsave'
  if(tit) then
     if (ltit>0) then
        write(iunt,10) '%% ----------------------------------------------------- TITLE '
        write(iunt,*) '/Helvetica findfont ',fnstit,' cm scalefont setfont '
        write(iunt,*) xtit,' cm ',ytit,' cm moveto '
        write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
        write(iunt,10) '%% ----------------------------------------------------- END TITLE '
     end if
  end if
  write(iunt,*) lrmrgn,' cm ',botmrgn,' cm translate'
  write(iunt,*) siz,' cm ',m,' div dup scale '
  !
  ! draw a frame around the matrix
  !
  if(box) then
     write(iunt,10) '%% ----------------------------------------------------- BOX '
     write(iunt,*)  frlw,' setlinewidth'
     write(iunt,10) 'newpath'
     write(iunt,11) 0, 0, ' moveto'
     if (square) then
        write(iunt,11) m, 0, ' lineto'
        write(iunt,11) m, m, ' lineto'
        write(iunt,11) 0, m, ' lineto'
     else
        write(iunt,11) nc, 0,' lineto'
        write(iunt,11) nc,nr,' lineto'
        write(iunt,11) 0, nr,' lineto'
     end if
     write(iunt,10) 'closepath stroke'
     write(iunt,10) '%% ----------------------------------------------------- END BOX '
  end if
  !
  ! drawing the separation lines 
  !
  write(iunt,*)  ' 0.2 setlinewidth'
  do kol=1, nlines 
     isep = lines(kol) 
     !
     ! horizontal lines 
     !
     yy = real(nrow-isep) + haf 
     xx = real(ncol+1) 
     write(iunt,13) zero, yy, ' moveto '
     write(iunt,13) xx,   yy, ' lineto stroke '
     !
     ! vertical lines 
     !
     xx = real(isep) + haf 
     yy = real(nrow+1)  
     write(iunt,13) xx, zero,' moveto '
     write(iunt,13) xx, yy,  ' lineto stroke '             
  end do
  ! 
  !----------- plotting loop ---------------------------------------------
  !
  write(iunt,10) '1 1 translate'
  write(iunt,10) '0.8 setlinewidth'
  write(iunt,10) '/p {moveto 0 -.40 rmoveto '
  write(iunt,10) '           0  .80 rlineto stroke} def'
  write(iunt,10) 'cnze'

  n=n/ndof
  do ii=1,n
     istart = ia(ii)
     ilast  = ia(ii+1)-1 
     if(mode==1) then
        do k=istart,ilast
           write(iunt,11) ii-1, nrow-ja(k), ' p'
        end do
     else
        do k=istart,ilast
           do idof=1,ndof
              do jdof=1,ndof
                 if( izero == 0 ) then
                    !
                    ! Null coefficients do not appear
                    !
                    if(abs(amatr(jdof,idof,k))<=1.0e-12_rp) then
                       continue
                    else
                       write(iunt,11) (ja(k)-1)*ndof+jdof-1, nrow-((ii-1)*ndof+idof), ' p'
                    end if
                 else if( izero == 1 ) then
                    !
                    ! Null coefficients are like others
                    !
                    write(iunt,11) (ja(k)-1)*ndof+jdof-1, nrow-((ii-1)*ndof+idof), ' p'
                 else if( izero == 2 ) then
                    !
                    ! Null coefficients are clearer
                    !
                    if(abs(amatr(jdof,idof,k))<=1.0e-12_rp) then
                       write(iunt,10) 'czer'
                       write(iunt,11) (ja(k)-1)*ndof+jdof-1, nrow-((ii-1)*ndof+idof), ' p'
                       write(iunt,10) 'cnze'
                    else
                       write(iunt,11) (ja(k)-1)*ndof+jdof-1, nrow-((ii-1)*ndof+idof), ' p'
                    end if
                 end if
                 if(number) then
                    if(ja(k)==ii) then
                       write(iunt,10) 'cnum'
                       write(iunt,10) 'fnum'
                       write(iunt,13) real((ja(k)-1)*ndof+jdof-1)-xf,real(nrow-((ii-1)*ndof+idof))-yf,' moveto'
                       write(messa,*) ii
                       messa=adjustl(messa)
                       messa='( '//adjustl(trim(messa))//')  Cshow'
                       write(iunt,10) trim(messa)
                       write(iunt,10) 'cnze'
                    end if
                 end if
              end do
           end do
        end do
        !
        ! add diagonal element if MSR mode
        !
        if(mode==2) write(iunt,11) ii-1, nrow-ii, ' p'
     endif
  end do
  write(iunt,10) 'showpage'

10 format (a)
11 format (2(i10,1x),a)
12 format (a,4(1x,f9.2))
13 format (2(f9.2,1x),a)
14 format (a,i2,a)
15 format (f9.2,1x,a)

end subroutine pspltm
