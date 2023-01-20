!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_sld_shell
  use def_kintyp, only :ip,rp
  implicit none
  private

  public :: STIFCSTDKT
  public :: STIFPOR3D

contains

  !*************************************************************************
  !     [STIFCSTDKT.F]
  !*************************************************************************
  !
  !     OBTIENE MATRIZ DE RIGIDEZ [KE] - ELEMENTO CST+DKT EN EL ESPACIO 3D
  !     INPUT: T(ESPESOR), A(AREA),E-MODULO DE YOUNG, POI-POISSON
  !            XC(1:NDIME,3) - COORDENADAS DE LOS NODOS (X-Y-Z)
  !            EST=1, ESTADO DE TENSION PLANA / EST=2-DEF. PLANA
  !     OUTPUT: KE(18,18) - MATRIZ DE RIGIDEZ
  !*************************************************************************

  SUBROUTINE STIFCSTDKT(XC,E,POI,T,EST,KE)
    IMPLICIT NONE
    !
    !*** GLOBAL
    !
    integer(ip), intent(in)  :: EST         !< State of tension
    real(rp),    intent(in)  :: T           !< Thickness
    real(rp),    intent(in)  :: E           !< Young modulus
    real(rp),    intent(in)  :: POI         !< Poisson coefficient
    real(rp),    intent(in)  :: XC(3,3)     !< Node coordinates
    real(rp),    intent(out) :: KE(18,18)   !< Stifness matrix
    !
    !*** LOCAL
    !
    integer(ip) :: I,J,K,CON(3),NODO,VET1(6),VET2(9),COL
    real(rp)    :: XCL(3,2),KECST(6,6),KEDKT(9,9),KEL(15,15)
    real(rp)    :: X(3),Y(3),Z(3),AREA,L,M,N,XL(2),YL(2),ZL(2),LONG
    real(rp)    :: VEC(3),X_X1,X_Y2,X_Z3,Y_X1,Y_Y2,Y_Z3,MATA(15,18)
    real(rp)    :: Z_X1,Z_Y2,Z_Z3
    real(rp)    :: VERX(3),VERY(3),VERZ(3)
    !
    !*** OBTENCION DE LOS VERSORES DEL SISTEMA LOCAL
    !
    DO I = 1,3
       X(I) = XC(I,1)
       Y(I) = XC(I,2)
       Z(I) = XC(I,3)
       IF (I /=3 ) THEN
          XL(I) = XC(I,1)
          YL(I) = XC(I,2)
          ZL(I) = XC(I,3)
       ENDIF
    ENDDO
    ! VERSOR SEGUN DIRECCION 12- VERSOR x
    CALL LONGB (XL,YL,ZL,LONG,L,M,N)
    VERX(1) = L
    VERX(2) = M
    VERX(3) = N
    X_X1    = L
    X_Y2    = M
    X_Z3    = N
    ! CALCULO DEL AREA Y DEL VERSOR NORMAL AL ELEMENTO - VERSOR z=VEC
    CALL AREATR3D (X,Y,Z,VERZ,AREA)
    CALL AREATR3D (X,Y,Z,VEC,AREA)

    Z_X1    = VEC(1)
    Z_Y2    = VEC(2)
    Z_Z3    = VEC(3)
    ! OBTENCIÓN DEL VERSOR y - PRODUCTO VECTORIAL
    VERY(1) =  VEC(2)*N-VEC(3)*M
    VERY(2) = -(VEC(1)*N-VEC(3)*L)
    VERY(3) =  VEC(1)*M-VEC(2)*L
    Y_X1    =  VEC(2)*N-VEC(3)*M
    Y_Y2    = -(VEC(1)*N-VEC(3)*L)
    Y_Z3    =  VEC(1)*M-VEC(2)*L
    !
    !*** MATRIZ DE TRANSFORMACIÓN-[MATA(15,18)] XYZ=>xyz
    !
    MATA = 0.0_rp
    ! NODO 1
    DO I = 1,3 !NODOS
       COL = 6*I-6
       K   = 5*I
       DO J=1,3
          MATA(K-4,COL+J)   = VERX(J)
          MATA(K-3,COL+J)   = VERY(J)
          MATA(K-2,COL+J)   = VERZ(J)
          MATA(K-1,COL+3+J) = VERX(J)
          MATA(K,COL+3+J)   = VERY(J)
       ENDDO
    ENDDO

    GOTO 20
    MATA(1,1)   = X_X1
    MATA(1,2)   = X_Y2
    MATA(1,3)   = X_Z3
    MATA(2,1)   = Y_X1
    MATA(2,2)   = Y_Y2
    MATA(2,3)   = Y_Z3
    MATA(3,1)   = Z_X1
    MATA(3,2)   = Z_Y2
    MATA(3,3)   = Z_Z3
    MATA(4,4)   = X_X1
    MATA(4,5)   = X_Y2
    MATA(4,6)   = X_Z3
    MATA(5,4)   = Y_X1
    MATA(5,5)   = Y_Y2
    MATA(5,6)   = Y_Z3
    ! NODO 2
    MATA(6,7)   = X_X1
    MATA(6,8)   = X_Y2
    MATA(6,9)   = X_Z3
    MATA(7,7)   = Y_X1
    MATA(7,8)   = Y_Y2
    MATA(7,9)   = Y_Z3
    MATA(8,7)   = Z_X1
    MATA(8,8)   = Z_Y2
    MATA(8,9)   = Z_Z3
    MATA(9,10)  = X_X1
    MATA(9,11)  = X_Y2
    MATA(9,12)  = X_Z3
    MATA(10,10) = Y_X1
    MATA(10,11) = Y_Y2
    MATA(10,12) = Y_Z3
    ! NODO 3
    MATA(11,13) = X_X1
    MATA(11,14) = X_Y2
    MATA(11,15) = X_Z3
    MATA(12,13) = Y_X1
    MATA(12,14) = Y_Y2
    MATA(12,15) = Y_Z3
    MATA(13,13) = Z_X1
    MATA(13,14) = Z_Y2
    MATA(13,15) = Z_Z3
    MATA(14,16) = X_X1
    MATA(14,17) = X_Y2
    MATA(14,18) = X_Z3
    MATA(15,16) = Y_X1
    MATA(15,17) = Y_Y2
    MATA(15,18) = Y_Z3
    !
    !*** CONSTRUCCION DE LA MATRIZ DE RIGIDEZ [KeCST(6,6)] LOCAL
    !
    !  COORDENADAS LOCALES DE LOS NODOS
20  CONTINUE
    XCL(1,1) = 0.0_rp
    XCL(1,2) = 0.0_rp
    XCL(2,1) = LONG
    XCL(2,2) = 0.0_rp
    XCL(3,1) = (X(3)-X(1))*L+(Y(3)-Y(1))*M+(Z(3)-Z(1))*N
    XCL(3,2) = (X(3)-X(1))*Y_X1+(Y(3)-Y(1))*Y_Y2+(Z(3)-Z(1))*Y_Z3

    CALL STIFCST(XCL,E,POI,T,EST,KECST)
    !
    !*** CONSTRUCCION DE LA MATRIZ DE RIGIDEZ [KeDKT(9,9)] LOCAL
    !
    CALL STIFDKT(XCL,E,POI,T,KEDKT)

    !     =================================================================
    !	   CONSTRUCTRUCCION [KEL(15,15)]=[KECST]+[KEDKT]
    !     =================================================================

    CON(1) =  1
    CON(2) =  6
    CON(3) = 11
    !*** CONSTRUCCION DEL VECTOR AUXILIAR VET1(6), VET2(9) PARA ENSAMBLAJE
    K = 1
    DO I = 1,3 !NODOS DEL ELEMENTO
       NODO = CON(I)
       DO J = 1,2 !Num. grados de libertad del elemento
          VET1(K) = NODO+(J-1)
          K = K+1
       ENDDO
    ENDDO
    K = 1
    DO I = 1,3 !NODOS DEL ELEMENTO
       NODO = CON(I)
       DO J = 1,3 !Num. grados de libertad del elemento
          VET2(K) = NODO+(J+1)
          K = K+1
       ENDDO
    ENDDO

    KEL = 0.0_rp

    DO I = 1,9
       DO J = 1,9
          IF ( I<=6 .AND. J<=6 ) THEN !CONTRIBUCION CST
             KEL(VET1(I),VET1(J)) = KEL(VET1(I),VET1(J))+KECST(I,J)
          ENDIF
          KEL(VET2(I),VET2(J)) = KEL(VET2(I),VET2(J))+KEDKT(I,J) !CONTRIBUCION DKT
       ENDDO
    ENDDO

    GOTO 10
    ! CST-1 FILA

    KEL(1,1)=KECST(1,1)
    KEL(1,2)=KECST(1,2)
    KEL(1,6)=KECST(1,3)
    KEL(1,7)=KECST(1,4)
    KEL(1,11)=KECST(1,5)
    KEL(1,12)=KECST(1,6)
    ! CST-2 FILA
    KEL(2,2)=KECST(2,2)
    KEL(2,6)=KECST(2,3)
    KEL(2,7)=KECST(2,4)
    KEL(2,11)=KECST(2,5)
    KEL(2,12)=KECST(2,6)
    ! CST-3FILA
    KEL(6,6)=KECST(3,3)
    KEL(6,7)=KECST(3,4)
    KEL(6,11)=KECST(3,5)
    KEL(6,12)=KECST(3,6)
    ! CST-4FILA
    KEL(7,7)=KECST(4,4)
    KEL(7,11)=KECST(4,5)
    KEL(7,12)=KECST(4,6)
    ! CST -5FILA
    KEL(11,11)=KECST(5,5)
    KEL(11,12)=KECST(5,6)
    ! CST -6FILA
    KEL(12,12)=KECST(6,6)
    !
    !** CONTRIBUCION DEL DKT
    !
    ! DKT-1FILA
    KEL(3,3)=KEDKT(1,1)
    KEL(3,4)=KEDKT(1,2)
    KEL(3,5)=KEDKT(1,3)
    KEL(3,8)=KEDKT(1,4)
    KEL(3,9)=KEDKT(1,5)
    KEL(3,10)=KEDKT(1,6)
    KEL(3,13)=KEDKT(1,7)
    KEL(3,14)=KEDKT(1,8)
    KEL(3,15)=KEDKT(1,9)
    ! DKT-2FILA
    KEL(4,4)=KEDKT(2,2)
    KEL(4,5)=KEDKT(2,3)
    KEL(4,8)=KEDKT(2,4)
    KEL(4,9)=KEDKT(2,5)
    KEL(4,10)=KEDKT(2,6)
    KEL(4,13)=KEDKT(2,7)
    KEL(4,14)=KEDKT(2,8)
    KEL(4,15)=KEDKT(2,9)
    ! DKT-3FILA
    KEL(5,5)=KEDKT(3,3)
    KEL(5,8)=KEDKT(3,4)
    KEL(5,9)=KEDKT(3,5)
    KEL(5,10)=KEDKT(3,6)
    KEL(5,13)=KEDKT(3,7)
    KEL(5,14)=KEDKT(3,8)
    KEL(5,15)=KEDKT(3,9)
    ! DKT-4FILA
    KEL(8,8)=KEDKT(4,4)
    KEL(8,9)=KEDKT(4,5)
    KEL(8,10)=KEDKT(4,6)
    KEL(8,13)=KEDKT(4,7)
    KEL(8,14)=KEDKT(4,8)
    KEL(8,15)=KEDKT(4,9)
    ! DKT-5FILA
    KEL(9,9)=KEDKT(5,5)
    KEL(9,10)=KEDKT(5,6)
    KEL(9,13)=KEDKT(5,7)
    KEL(9,14)=KEDKT(5,8)
    KEL(9,15)=KEDKT(5,9)
    ! DKT-6FILA
    KEL(10,10)=KEDKT(6,6)
    KEL(10,13)=KEDKT(6,7)
    KEL(10,14)=KEDKT(6,8)
    KEL(10,15)=KEDKT(6,9)
    ! DKT-7FILA
    KEL(13,13)=KEDKT(7,7)
    KEL(13,14)=KEDKT(7,8)
    KEL(13,15)=KEDKT(7,9)
    ! DKT-8FILA
    KEL(14,14)=KEDKT(8,8)
    KEL(14,15)=KEDKT(8,9)
    ! DKT-9FILA
    KEL(15,15)=KEDKT(9,9)
    !
    !*** PARTE SIMETRICA
    !
    DO I=1,15
       DO J=1,(I-1)
          KEL(I,J)=KEL(J,I)
       ENDDO
    ENDDO

10  CONTINUE
    !
    !     =================================================================
    !	   CALCULO DE [KE]=[MATA]t*[KeL(15,15)]*[MATA(15,18)] - GLOBAL
    !     =================================================================
    !
    CALL BDBCO1(MATA,KEL,KE,15_ip,18_ip,1_ip)
    !

    RETURN
  END SUBROUTINE STIFCSTDKT

  !*************************************************************************
  !     [ESCSTDKT.F]
  !*************************************************************************
  !
  !     OBTIENE DEFORMACIONxTENSION - ELEMENTO CST+DKT EN EL ESPACIO 3D (SISTEMA GLOBAL)
  !     INPUT: T(ESPESOR), A(AREA),E-MODULO DE YOUNG, POI-POISSON
  !            XC(3,3) - COORDENADAS DE LOS NODOS (X-Y-Z)
  !            EST=1, ESTADO DE TENSION PLANA / EST=2-DEF. PLANA
  !     OUTPUT: STRAIN_A, STRES_A
  !*************************************************************************

  SUBROUTINE ESCSTDKT(XC,NPAR,P_MAT,UE,STRAIN_A,STRES_A)
    IMPLICIT NONE
    !
    !*** GLOBAL
    !
    integer(ip) :: NPAR
    real(rp)    :: XC(3,6),UE(18),STRAIN_A(3,6),STRES_A(3,6)
    real(rp)    :: P_MAT(NPAR)
    !
    !*** LOCAL
    !
    integer(ip) :: I,J,INODE
    real(rp)    :: XCL(3,2)
    real(rp)    :: X(3),Y(3),Z(3),AREA,L,M,N,XL(2),YL(2),ZL(2),LONG
    real(rp)    :: VEC(3),X_X1,X_Y2,X_Z3,Y_X1,Y_Y2,Y_Z3,MATA(15,18)
    real(rp)    :: Z_X1,Z_Y2,Z_Z3
    real(rp)    :: UEL(15),UELCST(6),MA(3,3)
    real(rp)    :: EPL(3,3),EPG(3,3)
    real(rp)    :: STRAIN_AL(3,4),STRES_AL(3,4),YOUNG,POI,DMATX(4,4)
    !
    !*** OBTENCIÓN DE LOS  VERSORES DEL SISTEMA LOCAL
    !
    DO I=1,3
       X(I)=XC(I,1)
       Y(I)=XC(I,2)
       Z(I)=XC(I,3)
       IF (I/=3) THEN
          XL(I)=XC(I,1)
          YL(I)=XC(I,2)
          ZL(I)=XC(I,3)
       ENDIF
    ENDDO
    ! VERSOR SEGUN DIRECCION 12- VERSOR x
    CALL LONGB (XL,YL,ZL,LONG,L,M,N)
    X_X1=L
    X_Y2=M
    X_Z3=N
    ! CALCULO DEL AREA Y DEL VERSOR NORMAL AL ELEMENTO - VERSOR z=VEC
    CALL AREATR3D (X,Y,Z,VEC,AREA)
    Z_X1=VEC(1)
    Z_Y2=VEC(2)
    Z_Z3=VEC(3)
    ! OBTENCIÓN DEL VERSOR y - PRODUCTO VECTORIAL
    Y_X1 =  VEC(2)*N-VEC(3)*M
    Y_Y2 = -(VEC(1)*N-VEC(3)*L)
    Y_Z3 =  VEC(1)*M-VEC(2)*L
    !
    !*** MATRIZ DE TRANSFORMACIÓN-[MATA(6,9)] XYZ=>xyz
    !
    MATA=0.0_rp
    ! NODO 1
    MATA(1,1)=X_X1
    MATA(1,2)=X_Y2
    MATA(1,3)=X_Z3
    MATA(2,1)=Y_X1
    MATA(2,2)=Y_Y2
    MATA(2,3)=Y_Z3
    MATA(3,1)=Z_X1
    MATA(3,2)=Z_Y2
    MATA(3,3)=Z_Z3
    MATA(4,4)=X_X1
    MATA(4,5)=X_Y2
    MATA(4,6)=X_Z3
    MATA(5,4)=Y_X1
    MATA(5,5)=Y_Y2
    MATA(5,6)=Y_Z3
    ! NODO 2
    MATA(6,7)=X_X1
    MATA(6,8)=X_Y2
    MATA(6,9)=X_Z3
    MATA(7,7)=Y_X1
    MATA(7,8)=Y_Y2
    MATA(7,9)=Y_Z3
    MATA(8,7)=Z_X1
    MATA(8,8)=Z_Y2
    MATA(8,9)=Z_Z3
    MATA(9,10)=X_X1
    MATA(9,11)=X_Y2
    MATA(9,12)=X_Z3
    MATA(10,10)=Y_X1
    MATA(10,11)=Y_Y2
    MATA(10,12)=Y_Z3
    ! NODO 3
    MATA(11,13)=X_X1
    MATA(11,14)=X_Y2
    MATA(11,15)=X_Z3
    MATA(12,13)=Y_X1
    MATA(12,14)=Y_Y2
    MATA(12,15)=Y_Z3
    MATA(13,13)=Z_X1
    MATA(13,14)=Z_Y2
    MATA(13,15)=Z_Z3
    MATA(14,16)=X_X1
    MATA(14,17)=X_Y2
    MATA(14,18)=X_Z3
    MATA(15,16)=Y_X1
    MATA(15,17)=Y_Y2
    MATA(15,18)=Y_Z3

    !
    !*** CONSTRUCCION DE LA MATRIZ DE RIGIDEZ [KeCST(6,6)] LOCAL
    !
    !  COORDENADAS LOCALES DE LOS NODOS
    XCL(1,1) = 0.0_rp
    XCL(1,2) = 0.0_rp
    XCL(2,1) = LONG
    XCL(2,2) = 0.0_rp
    XCL(3,1) = (X(3)-X(1))*L+(Y(3)-Y(1))*M+(Z(3)-Z(1))*N
    XCL(3,2) = (X(3)-X(1))*Y_X1+(Y(3)-Y(1))*Y_Y2+(Z(3)-Z(1))*Y_Z3

    !
    !     =================================================================
    !	   CALCULO DE LOS DESPLAZAMIENTOS LOCALES [UEL]=[MATA(15,18)]*[UE(18)] - GLOBAL
    !     =================================================================
    !
    DO I=1,15
       UEL(I)=0.0_rp
       DO J=1,18
          UEL(I)=UEL(I)+MATA(I,J)*UE(J)
       ENDDO
    ENDDO

    UELCST(1) = UEL(1)
    UELCST(2) = UEL(2)
    UELCST(3) = UEL(6)
    UELCST(4) = UEL(7)
    UELCST(5) = UEL(11)
    UELCST(6) = UEL(12)
    !
    !*** DEFORMACION EN EL SISTEMA LOCAL
    !
    CALL STRACST(XCL,UELCST,STRAIN_AL)
    !
    !*** STRESS
    !
    YOUNG = P_MAT(1)
    POI   = P_MAT(2)

    CALL MATD_2D(1_ip,YOUNG,POI,DMATX)

    DO INODE=1,3
       DO I=1,3
          STRES_AL(INODE,I) = 0.0_rp
          DO J=1,3
             STRES_AL(INODE,I)=STRES_AL(INODE,I)+DMATX(I,J)*STRAIN_AL(INODE,J)
          ENDDO
       ENDDO
       STRAIN_AL(INODE,4) = -(POI/YOUNG)*(STRES_AL(INODE,1)+STRES_AL(INODE,2))
       STRES_AL(INODE,4)  = 0.0_rp
    ENDDO

    !
    !*** TRANSFORMACION DE LAS COMPONETNES DEL TENSOR DE DEFORMACION P/ SISTEMA GLOBAL
    !       [E]=[A]t[E'][A]
    !

    MA(1,1)=X_X1
    MA(1,2)=X_Y2
    MA(1,3)=X_Z3
    MA(2,1)=Y_X1
    MA(2,2)=Y_Y2
    MA(2,3)=Y_Z3
    MA(3,1)=Z_X1
    MA(3,2)=Z_Y2
    MA(3,3)=Z_Z3
    !
    !*** DEFORMACION
    !
    EPL=0.0_rp
    EPL(1,1)=STRAIN_AL(1,1)
    EPL(2,2)=STRAIN_AL(1,2)
    EPL(1,2)=STRAIN_AL(1,3)
    EPL(2,1)=EPL(1,2)
    EPL(3,3)=STRAIN_AL(1,4)

    CALL BDBCO1(MA,EPL,EPG,3_ip,3_ip,1_ip)
    !
    !** organizacion: Exx,Eyy,Ezz,Exy,Eyz,Exz
    !
    DO I=1,3
       STRAIN_A(I,1)=EPG(1,1)
       STRAIN_A(I,2)=EPG(2,2)
       STRAIN_A(I,3)=EPG(3,3)
       STRAIN_A(I,4)=EPG(1,2)
       STRAIN_A(I,5)=EPG(2,3)
       STRAIN_A(I,6)=EPG(1,3)
    ENDDO
    !
    !*** STRESS
    !
    EPL=0.0_rp
    EPL(1,1)=STRES_AL(1,1)
    EPL(2,2)=STRES_AL(1,2)
    EPL(1,2)=STRES_AL(1,3)
    EPL(2,1)=EPL(1,2)
    CALL BDBCO1(MA,EPL,EPG,3_ip,3_ip,1_ip)
    !
    !** organizacion: Sxx,Syy,Szz,Sxy,Syz,Sxz
    !
    DO I=1,3
       STRES_A(I,1)=EPG(1,1)
       STRES_A(I,2)=EPG(2,2)
       STRES_A(I,3)=EPG(3,3)
       STRES_A(I,4)=EPG(1,2)
       STRES_A(I,5)=EPG(2,3)
       STRES_A(I,6)=EPG(1,3)
    ENDDO

  END SUBROUTINE ESCSTDKT

  !*************************************************************************
  !    [MATB_CST.F]
  !*************************************************************************
  !*** MATRIZ [B] - CST
  !    INPUT:
  !          XC: COORDENADAS DE LOS VÉRTICES
  !    OUTPUT:
  !          AREA-Area del elmento; BMAT-matriz [B] del elemento CST
  !*************************************************************************

  SUBROUTINE MATB_CST(XC,AREA,BMAT)
    IMPLICIT NONE
    real(rp) :: XC(3,2),AREA,BMAT(3,6)
    integer(ip) :: I
    real(rp) :: B(3),C(3),X(3),Y(3)

    DO I=1,3
       X(I)=XC(I,1)
       Y(I)=XC(I,2)
    ENDDO
    !
    !*** CALCULO DEL AREA
    !

    CALL AREATR2D (X,Y,AREA)

    B(1)=Y(2)-Y(3)
    B(2)=Y(3)-Y(1)
    B(3)=Y(1)-Y(2)
    C(1)=X(3)-X(2)
    C(2)=X(1)-X(3)
    C(3)=X(2)-X(1)

    !
    BMAT=0.0_rp
    !
    BMAT(1,1)=B(1)/(2.0_rp*AREA)
    BMAT(2,2)=C(1)/(2.0_rp*AREA)
    BMAT(3,1)=BMAT(2,2)
    BMAT(3,2)=BMAT(1,1)
    !
    BMAT(1,3)=B(2)/(2.0_rp*AREA)
    BMAT(2,4)=C(2)/(2.0_rp*AREA)
    BMAT(3,3)=BMAT(2,4)
    BMAT(3,4)=BMAT(1,3)
    !
    BMAT(1,5)=B(3)/(2.0_rp*AREA)
    BMAT(2,6)=C(3)/(2.0_rp*AREA)
    BMAT(3,5)=BMAT(2,6)
    BMAT(3,6)=BMAT(1,5)
    !
  END SUBROUTINE MATB_CST

  !*****************************************************************
  !    [STIFDKT.F]
  !*****************************************************************
  !     MATRIZ DE RIGIDEZ DEL ELEMENTO DE PLACA DKT- COORDENADAS GLOBALES
  !    Input: T:espesor; E: Módulo de Young; NU: Poisson
  !           X(3),Y(3): coordenadas de los nodos
  !    Output: KE(12,12): matriz de rigidez del elemento DKT
  !*****************************************************************
  SUBROUTINE STIFDKT(XC,E,NU,T,KE)
    IMPLICIT NONE
    !
    !*** GLOBAL
    !
    real(rp) :: XC(3,2),T,KE(9,9),E,NU
    !
    !*** LOCAL
    !
    integer(ip) :: I,J,K1,I1,II,K2,JJ,KOD(2,9)
    real(rp)    :: X(3),Y(3),DET
    real(rp)    :: D(3,3),GG(10,9),DD(9,9),SS(9,9),QQ(9,9),PP(3,3)
    real(rp)    :: PT(2,3),RS(2,3),Q(3),BB(3),CC(3),ALS(3)

    DATA KOD/1,1,2,3,3,2,4,4,5,6,6,5,7,7,8,9,9,8/
    DATA PP/12.0_rp,3*4.0_rp,2.0_rp,1.0_rp,4.0_rp,1.0_rp,2.0_rp/

    DO I=1,3
       X(I)=XC(I,1)
       Y(I)=XC(I,2)
    ENDDO

    CALL MATD_PLACA(T,E,NU,D)

    BB(1)=Y(2)-Y(3)
    BB(2)=Y(3)-Y(1)
    BB(3)=Y(1)-Y(2)
    CC(1)=X(3)-X(2)
    CC(2)=X(1)-X(3)
    CC(3)=X(2)-X(1)
    DET=(BB(1)*CC(2)-BB(2)*CC(1))*24
    DO I=1,3
       DO J=1,3
          DO K1=1,3
             II=(I-1)*3+K1
             DO K2=1,3
                JJ=(J-1)*3+K2
                DD(II,JJ)=D(I,J)*PP(K1,K2)/DET
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO I=1,3
       ALS(I)=BB(I)*BB(I)+CC(I)*CC(I)
       PT(1,I)=6*CC(I)/ALS(I)
       PT(2,I)=6*BB(I)/ALS(I)
       RS(1,I)=3*CC(I)**2/ALS(I)
       RS(2,I)=3*BB(I)**2/ALS(I)
       Q(I)=3*BB(I)*CC(I)/ALS(I)
    ENDDO

    DO I=1,2
       II=(I-1)*5
       GG(II+1,KOD(I,1))=PT(I,3)
       GG(II+2,KOD(I,1))=-PT(I,2)
       GG(II+3,KOD(I,1))=-PT(I,3)
       GG(II+4,KOD(I,1))=PT(I,2)-PT(I,3)
       GG(II+5,KOD(I,1))=PT(I,2)
       GG(II+1,KOD(I,2))=-Q(3)
       GG(II+2,KOD(I,2))=-Q(2)
       GG(II+3,KOD(I,2))=Q(3)
       GG(II+4,KOD(I,2))=Q(2)+Q(3)
       GG(II+5,KOD(I,2))=Q(2)
       GG(II+1,KOD(I,3))=-1-RS(I,3)
       GG(II+2,KOD(I,3))=-1-RS(I,2)
       GG(II+3,KOD(I,3))=RS(I,3)
       GG(II+4,KOD(I,3))=RS(I,2)+RS(I,3)
       GG(II+5,KOD(I,3))=RS(I,2)
       GG(II+1,KOD(I,3))=-1-RS(I,3)
       GG(II+1,KOD(I,4))=-PT(I,3)
       GG(II+3,KOD(I,4))=PT(I,3)
       GG(II+4,KOD(I,4))=PT(I,1)+PT(I,3)
       GG(II+1,KOD(I,5))=-Q(3)
       GG(II+3,KOD(I,5))=Q(3)
       GG(II+4,KOD(I,5))=Q(3)-Q(1)
       GG(II+1,KOD(I,6))=1-RS(I,3)
       GG(II+3,KOD(I,6))=RS(I,3)
       GG(II+4,KOD(I,6))=RS(I,3)-RS(I,1)
       GG(II+2,KOD(I,7))=PT(I,2)
       GG(II+4,KOD(I,7))=-PT(I,1)-PT(I,2)
       GG(II+5,KOD(I,7))=-PT(I,2)
       GG(II+2,KOD(I,8))=-Q(2)
       GG(II+4,KOD(I,8))=Q(2)-Q(1)
       GG(II+5,KOD(I,8))=Q(2)
       GG(II+2,KOD(I,9))=1-RS(I,2)
       GG(II+4,KOD(I,9))=RS(I,2)-RS(I,1)
       GG(II+5,KOD(I,9))=RS(I,2)
    ENDDO
    DO I=1,9
       QQ(1,I)=BB(2)*GG(1,I)+BB(3)*GG(2,I)
       QQ(2,I)=2*BB(2)*GG(3,I)+BB(3)*GG(4,I)
       QQ(3,I)=BB(2)*GG(4,I)+2*BB(3)*GG(5,I)
       QQ(4,I)=-CC(2)*GG(6,I)-CC(3)*GG(7,I)
       QQ(5,I)=-2*CC(2)*GG(8,I)-CC(3)*GG(9,I)
       QQ(6,I)=-CC(2)*GG(9,I)-2*CC(3)*GG(10,I)
       QQ(7,I)=CC(2)*GG(1,I)+CC(3)*GG(2,I)-BB(2)*GG(6,I)-BB(3)*GG(7,I)
       QQ(8,I)=2*CC(2)*GG(3,I)+CC(3)*GG(4,I)-2*BB(2)*GG(8,I)-BB(3)*GG(9,I)
       QQ(9,I)=CC(2)*GG(4,I)+2*CC(3)*GG(5,I)-BB(2)*GG(9,I)-2*BB(3)*GG(10,I)
    ENDDO

    DO I=1,9
       DO J=1,9
          SS(I,J)=0.0_rp
          DO I1=1,9
             SS(I,J)=SS(I,J)+DD(I,I1)*QQ(I1,J)
          ENDDO
       ENDDO
    ENDDO
    DO I=1,9
       DO J=1,9
          KE(I,J)=0.0_rp
          DO I1=1,9
             KE(I,J)=KE(I,J)+QQ(I1,I)*SS(I1,J)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE STIFDKT

  !************************************************************
  !    [STRACST.F]
  !************************************************************
  ! CALCULO DE LA DEFORMACION EN EL ELEMENTO CST
  !  DEFORMACION CONSTANTES
  !************************************************************
  SUBROUTINE STRACST(XC,UE,STRAIN_A)
    IMPLICIT NONE
    !
    !*** GLOBALES
    !
    real(rp)    :: XC(3,2),UE(6),STRAIN_A(3,4)
    !
    !*** LOCALES
    !
    integer(ip) :: I,J
    real(rp)    :: BMAT(3,6),VEC(4),AREA

    !
    !** CALCULO DE AREA Y DE LOS COEFICIENTES c,b, DE LA MATRIZ [B]
    !

    CALL MATB_CST(XC,AREA,BMAT)

    !
    !** EFECTUANDO LA MULTIPLICACION DE MATRICIESO [B][UE]
    !
    DO I=1,3
       VEC(I)=0.0_rp
       DO J=1,6
          VEC(I)=VEC(I)+BMAT(I,J)*UE(J)
       ENDDO
    ENDDO
    !
    !*** LA TENSION ES CONSTANTE EN EL ELEMENTO
    !
    DO I=1,3
       DO J=1,3
          STRAIN_A(I,J)=VEC(J)
       ENDDO
    ENDDO

    !
    RETURN
  END SUBROUTINE STRACST
  !*****************************************************
  !       SUBROUTINE AREATR2D.F
  !*****************************************************
  ! Esta subroutina obtiene el área de un triángulo definido por tres puntos
  ! INPUT: X(3),Y(3): Coordenadas de los nodos
  ! OUTPUT: AREA (area del elemento triangular)
  !
  !*****************************************************

  SUBROUTINE AREATR2D (X,Y,AREA)
    IMPLICIT NONE
    real(rp) :: X(3),Y(3),AREA

    AREA = X(1)*(Y(2)-Y(3))-X(2)*(Y(1)-Y(3))+X(3)*(Y(1)-Y(2))
    AREA = AREA/2.0_rp

    IF (AREA<=0.0_rp) THEN
       WRITE(*,*) 'AREA NEGATIVA, VERIFICAR NUMERACIÓN DEL ELEMENTO'
    ENDIF

    RETURN
  END SUBROUTINE AREATR2D

  !*****************************************************
  !       [AREATR3D.F]
  !*****************************************************
  ! Esta subroutina obtiene el área de un triángulo definido por tres puntos
  ! Obtiene tambien el versor(vector unitario) normal al elemento de área
  ! INPUT: X(3),Y(3),Z(3): Coordenadas de los nodos
  ! OUTPUT: AREA (area del elemento triangular)
  !         VEC(3) - Componentes del versor normal al elemento de área
  !
  !*****************************************************

  SUBROUTINE AREATR3D (X,Y,Z,VEC,AREA)
    IMPLICIT NONE
    real(rp) :: X(3),Y(3),Z(3)
    real(rp) :: VEC(3),AREA
    !
    VEC(1) = ((Y(2)-Y(1))*(Z(3)-Z(1)))-((Z(2)-Z(1))*(Y(3)-Y(1)))
    VEC(2) = -(((X(2)-X(1))*(Z(3)-Z(1)))-((Z(2)-Z(1))*(X(3)-X(1))))
    VEC(3) = ((X(2)-X(1))*(Y(3)-Y(1)))-((Y(2)-Y(1))*(X(3)-X(1)))
    ! modulo del vector VEC
    AREA = SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2)
    !
    !*** VERSOR NORMAL AL ELEMENTO DE AREA
    !
    VEC(1) = VEC(1)/AREA
    VEC(2) = VEC(2)/AREA
    VEC(3) = VEC(3)/AREA
    !
    !*** AREA DEL ELEMENTO TRIANGULAR
    !
    AREA = AREA/2.0_rp

  END SUBROUTINE AREATR3D

  !**********************************************************************
  !     [BDBCO1.F]
  !**********************************************************************
  !                               T
  !****THIS ROUTINE PERFORMS E = B D B OF TWO GENERAL MATRICES
  !                   ( upper triangle )
  !
  !                  BMATX = NRAWS * NCOLU
  !                  DMATX = NRAWS * NRAWS
  !                  EMATX = NCOLU * NCOLU
  !
  !**********************************************************************

  SUBROUTINE BDBCO1(BMATX,DMATX,EMATX,NRAWS,NCOLU,KSYMM)
    IMPLICIT NONE
    !
    integer(ip) :: KSYMM
    integer(ip) :: NCOLU,NRAWS
    real(rp) ::  BMATX(NRAWS,NCOLU), DMATX(NRAWS,NRAWS), EMATX(NCOLU,NCOLU)
    !
    integer(ip) :: ICOLU,JCOLU,KRAWS,LRAWS,INDEX
    !
    !
    DO ICOLU=1,NCOLU
       INDEX=ICOLU            ! 'SYMMETRIC CASE'
       IF(KSYMM==1) INDEX=1   ! 'UNSYMMETRIC CASE'
       DO JCOLU=INDEX,NCOLU
          EMATX(ICOLU,JCOLU)=0.0_rp
          DO KRAWS=1,NRAWS
             DO LRAWS=1,NRAWS
                EMATX(ICOLU,JCOLU) = EMATX(ICOLU,JCOLU)&
                     & + BMATX(KRAWS,ICOLU) * DMATX(KRAWS,LRAWS) * BMATX(LRAWS,JCOLU)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE BDBCO1

  !*********
  !  ESTA SUBROUTINA OBTIENE LA LONGITUD ENTRE DOS PUNTOS
  ! Input: Coordenadas X(2), Y(2), Z(2)
  ! Output: Longitud (LONG), cosenos directores L,M,N
  !*******
  SUBROUTINE LONGB (X,Y,Z,LONG,L,M,N)
    IMPLICIT NONE
    real(rp) :: X(2),Y(2),Z(2),LONG,L,M,N
    L    = (X(2)-X(1))**2+(Y(2)-Y(1))**2+(Z(2)-Z(1))**2
    LONG = SQRT(L)
    L    = (X(2)-X(1))/LONG
    M    = (Y(2)-Y(1))/LONG
    N    = (Z(2)-Z(1))/LONG

  END SUBROUTINE LONGB

  !*************************************************************************
  !     SUBROUTINE MATD_2D.F
  !*************************************************************************
  !
  !     OBTIENE MATRIZ CONSTITUTIVA [D(4,4)] - PARA CASO 2D
  !     INPUT: E-MODULO DE YOUNG, POI-POISSON
  !            EST=1, ESTADO DE TENSION PLANA / EST=2-DEF. PLANA
  !     OUTPUT: D(4,4) - MATRIZ CONSTITUTIVA-2D
  !  organizacion: Sxx,Syy,Sxy,Szz
  !*************************************************************************

  SUBROUTINE MATD_2D(EST,E,POI,D)
    IMPLICIT NONE

    integer(ip) :: EST
    real(rp)    :: E,POI,D(4,4),CT
    !
    !** DEFINE MATRIZ CONSTITUTIVA
    !
    D = 0.0_rp
    IF ( EST == 1 ) THEN
       !
       ! TENSION PLANA
       !
       D(1,1) = E/(1.0_rp-POI**2)
       D(1,2) = D(1,1)*POI
       D(3,3) = E/(2.0_rp*(1.0_rp+POI))
    ELSEIF ( EST == 2 ) THEN
       !
       ! DEFORMACION PLANA
       !
       CT     = E/((1.0_rp+POI)*(1.0_rp-2.0_rp*POI))
       D(1,1) = CT*(1.0_rp-POI)
       D(1,2) = CT*POI
       D(3,3) = E/(2.0_rp*(1.0_rp+POI))
       D(4,1) = D(1,2)
       D(4,2) = D(1,2)
    ENDIF

    D(2,1)=D(1,2)
    D(2,2)=D(1,1)

  END SUBROUTINE MATD_2D

  !     *****************************************************************
  !     SUB-ROUTINA PARA CALCULO DE LA MATRIZ CONSTITUTIVA ELASTICA
  !
  !     *****************************************************************
  SUBROUTINE MATD_PLACA(T,YOUNG,NU,D)
    IMPLICIT NONE
    real(rp), intent(in)  :: T
    real(rp), intent(in)  :: YOUNG
    real(rp), intent(in)  :: NU
    real(rp), intent(out) :: D(3,3)
    real(rp)              :: AUX

    AUX    = (YOUNG*T**3)/(12.0_rp*(1.0_rp-NU**2))
    D(1,1) = AUX
    D(1,2) = AUX*NU
    D(2,1) = D(1,2)
    D(2,2) = D(1,1)
    D(3,3) = AUX*(1-NU)/2.0_rp

  END SUBROUTINE MATD_PLACA

  !*************************************************************************
  !     [STIFCST.F]
  !*************************************************************************
  !
  !     OBTIENE MATRIZ DE RIGIDEZ [KE] - ELEMENTO CST
  !     INPUT: T(ESPESOR), A(AREA),E-MODULO DE YOUNG, POI-POISSON
  !            XC(3,2) - COORDENADAS DE LOS NODOS
  !            EST=1, ESTADO DE TENSION PLANA / EST=2-DEF. PLANA
  !     OUTPUT: KE(6,6) - MATRIZ DE RIGIDEZ DEL ELEMENTO CST
  !*************************************************************************

  SUBROUTINE STIFCST(XC,E,POI,T,EST,KE)
    IMPLICIT NONE
    !
    !*** GLOBAL
    !
    real(rp)    :: T,E,POI,XC(3,2),KE(6,6)
    real(rp)    :: D1,D2,D3,B(3),C(3),D(4,4)
    integer(ip) :: I,J,EST
    !
    !*** LOCAL
    !
    real(rp) :: X(3),Y(3),A

    DO I=1,3
       X(I) = XC(I,1)
       Y(I) = XC(I,2)
    ENDDO
    !
    !** MATRIZ CONSTITUTIVA [D]
    !
    CALL MATD_2D(EST,E,POI,D)
    D1 = D(1,1)
    D2 = D(1,2)
    D3 = D(3,3)
    !
    !*** CALCULO DEL AREA
    !

    CALL AREATR2D (X,Y,A)

    B(1) = Y(2)-Y(3)
    B(2) = Y(3)-Y(1)
    B(3) = Y(1)-Y(2)
    C(1) = X(3)-X(2)
    C(2) = X(1)-X(3)
    C(3) = X(2)-X(1)
    !
    !     =================================================================
    !	   CALCULO DE [KE] -triangular superior
    !     =================================================================
    !
    ! primera fila
    KE(1,1) = C(1)**2*D3+B(1)**2*D1
    KE(1,2) = B(1)*D2*C(1)+C(1)*D3*B(1)
    KE(1,3) = B(1)*D1*B(2)+C(1)*D3*C(2)
    KE(1,4) = B(1)*D2*C(2)+C(1)*D3*B(2)
    KE(1,5) = B(1)*D1*B(3)+C(1)*D3*C(3)
    KE(1,6) = B(1)*D2*C(3)+C(1)*D3*B(3)
    ! segunda fila
    KE(2,2) = C(1)**2*D1+B(1)**2*D3
    KE(2,3) = C(1)*D2*B(2)+B(1)*D3*C(2)
    KE(2,4) = C(1)*D1*C(2)+B(1)*D3*B(2)
    KE(2,5) = C(1)*D2*B(3)+B(1)*D3*C(3)
    KE(2,6) = C(1)*D1*C(3)+B(1)*D3*B(3)
    ! tercera fila
    KE(3,3) = C(2)**2*D3+B(2)**2*D1
    KE(3,4) = B(2)*D2*C(2)+C(2)*D3*B(2)
    KE(3,5) = B(2)*D1*B(3)+C(2)*D3*C(3)
    KE(3,6) = B(2)*D2*C(3)+C(2)*D3*B(3)
    ! cuarta fila
    KE(4,4) = C(2)**2*D1+B(2)**2*D3
    KE(4,5) = C(2)*D2*B(3)+B(2)*D3*C(3)
    KE(4,6) = C(2)*D1*C(3)+B(2)*D3*B(3)
    ! quinta fila
    KE(5,5) = C(3)**2*D3+B(3)**2*D1
    KE(5,6) = B(3)*D2*C(3)+C(3)*D3*B(3)
    ! sexta fila
    KE(6,6) = C(3)**2*D1+B(3)**2*D3

    !
    !*** PARTE SIMETRICA
    !
    DO I=1,6
       DO J=1,(I-1)
          KE(I,J)=KE(J,I)
       ENDDO
    ENDDO
    !
    KE=T/(4.0_rp*A)*KE
    !
  END SUBROUTINE STIFCST

!    *****************************************************************
!       [STIFPOR3D.F]
!     *****************************************************************
!     SUB-ROUTINA QUE PROPORCIONA LA MATRIZ DE RIGIDEZ DEL ELEMENTO PORTICO ESPACIAL
!             EN LAS COORDENADAS LOCALES (x-y-z)
!   Input:  L-longitud de la barra;
!           YOUNG: Módulo de Young;
!           NU: coeficiente de Poisson;
!           AREA: Area de la seccion transversal;
!           JT: Momento de inercia a torsión efectivo;
!           IY: Momento de inercia a flexion según eje y
!           IZ: Momento de inercia a flexion según eje z
!   Output: KE: matriz de rigediez en coordenadas locales
!
!     *****************************************************************
  SUBROUTINE STIFPOR3D(L,YOUNG,NU,AREA,JT,IY,IZ,KE)
    !
    !*** globales
    !
    REAL(rp)    :: L,JT,IY,IZ,KE(12,12),YOUNG,NU,AREA
    !
    !*** locales
    !
    INTEGER(ip) :: I,J
    REAL(rp)    :: EA,EIZ,EIY,GIX,G

    G=YOUNG/(2.0D0*(1.0D0+NU))

    KE=0.0D0

    EA=YOUNG*AREA/L
    EIz=YOUNG*IZ
    EIy=YOUNG*IY
    GIx=G*JT/L
    !
    !** TRIANGULAR SUPERIOR
    !
    ! FILA 1
    KE(1,1)=EA
    KE(1,7)=-EA
    ! FILA2
    KE(2,2)=12.0D0*EIZ/(L**3)
    KE(2,6)=6.0D0*EIZ/(L**2)
    KE(2,8)=-12.0D0*EIZ/(L**3)
    KE(2,12)=6.0D0*EIZ/(L**2)
    ! FILA 3
    KE(3,3)=12.0D0*EIY/(L**3)
    KE(3,5)=-6.0D0*EIY/(L**2)
    KE(3,9)=-12.0D0*EIY/(L**3)
    KE(3,11)=-6.0D0*EIY/(L**2)
    ! FILA 4
    KE(4,4)=GIX
    KE(4,10)=-GIX
    ! FILA 5
    KE(5,5)=4.0D0*EIY/L
    KE(5,9)=6.0D0*EIY/(L**2)
    KE(5,11)=2.0D0*EIY/L
    ! FILA 6
    KE(6,6)=4.0D0*EIZ/L
    KE(6,8)=-6.0D0*EIZ/(L**2)
    KE(6,12)=2.0D0*EIZ/L
    ! FILA 7
    KE(7,7)=EA
    ! FILA 8
    KE(8,8)=12.0D0*EIZ/(L**3)
    KE(8,12)=-6.0D0*EIZ/(L**2)
    ! FILA 9
    KE(9,9)=12.0D0*EIY/(L**3)
    KE(9,11)=6.0D0*EIY/(L**2)
    ! FILA 10
    KE(10,10)=GIX
    ! FILA 11
    KE(11,11)=4.0D0*EIY/L
    ! FILA 12
    KE(12,12)=4.0D0*EIZ/L

    !
    !*** PARTE SIMETRICA
    !
    DO I=1,12
       DO J=1,(I-1)
          KE(I,J)=KE(J,I)
       ENDDO
    ENDDO

  END SUBROUTINE STIFPOR3D

end module mod_sld_shell
