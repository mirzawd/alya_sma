!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



#ifdef VECTOR_SIZE
#define DVS VECTOR_SIZE
#else
#ifdef VECTOR_SIZE_CPU
#define DVS VECTOR_SIZE_CPU
#else
#define DVS 1
#endif
#endif

#if defined PNODE_VALUE && PGAUS_VALUE
#define ZELTY 37
#define ZDIME 3
#define ZDOFN 3
#define ZNODE PNODE_VALUE
#define ZGAUS PGAUS_VALUE
#define ZVOIG 6
#define ZDOFU 24
#else
#define ZELTY pelty
#define ZDIME ndime
#define ZDOFN ndofn_sld
#define ZGAUS pgaus
#define ZNODE pnode
#define ZVOIG (ndime-1)*3
#define ZDOFU pnode*ndofn_sld
#endif