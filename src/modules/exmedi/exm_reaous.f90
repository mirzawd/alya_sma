!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup exmedi
!> @{
!> @file    exm_reaous.f90
!> @author  multiple
!> @brief   Read postprocess data
!> @details Read postprocess data
!> @} 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!.md<module>exmedi
!.md<input>case.exm.dat
!.md<pos>2
!.md<sec>
      
subroutine exm_reaous
   !-----------------------------------------------------------------------
   !****f* Exmedi/exm_reaous
   ! NAME 
   !    exm_reaous
   ! DESCRIPTION
   !    This routine reads the output strategy 
   ! USES
   !    listen
   ! USED BY
   !    exm_turnon
   !***
   !-----------------------------------------------------------------------
   use      def_parame
   use      def_inpout
   use      def_master
   use      def_domain
   use      mod_memchk

   use      def_exmedi
   use mod_ecoute, only :  ecoute
   use mod_output_postprocess, only : output_postprocess_read

   implicit none
   
   if( INOTSLAVE ) then
      
      !
      ! Reach the section
      !
      rewind(lisda)
      call ecoute('exm_reaous')
      do while(words(1)/='OUTPU')
         call ecoute('exm_reaous')
      end do
      !
      ! Begin to read data
      !
      do while(words(1)/='ENDOU')
         call ecoute('exm_reaous')

         call output_postprocess_read()
         
         if(words(1)=='PARAM') then
         ! Postprocess general parameters
            call ecoute('exm_reaous')
            do while(words(1)/='ENDPA')
               if(words(1)=='ISOCH') then
                  call runend('EXM_REAOUS: ISOCHRONES THRESHOLD MUST BE WRITTEN IN PROPERTIES SECTION')
               end if
               call ecoute('exm_reaous')
            end do
         
         else if (words(1) == "SAVEC") then
            if (exists("CELLM")) then
               ! SAVE_CONVERGENCE CELL_MODEL -- saves ohara.mXXcXX files with curves for different 
               !concentrations during the ODE iterations for each material and cell type
               kfl_save_convergence_cellmodel = 1_ip 
            end if

         else if (words(1) == "DUMPI") then
            ! DUMP_INITIAL_CONDITIONS -- saves initial conditions to be hardcoded for a cell model 
            !concentrations during the ODE iterations for each material and cell type
            kfl_save_init_cellmodel = 1_ip 

         end if
      end do
   end if

end subroutine exm_reaous
       
!.md# Output and Post-process
!.md<code>
!.md<0>OUTPUT_&_POST_PROCESS
!
!.md<1>POSTPROCESS XXXXX, STEPS=int1 [, FILTER=int2]                                                   $ Variables computed on all nodes or filtered
!.md<field>POSTPROCESS
!.md<com>Postprocess variables on nodes at each int1 time steps. the name of the file
!.md<com>where the variables is stored at time step nnn is: <tt> sample-XXXXX-00000nnn.post.alyabin</tt>.
!.md<com>A fliter can be applied to the variable. Filters are defined in kermod. Variables:
!.md<com>
!.md<com>    - INTRA: intracellular potential
!.md<com>    - CALCI: calcium
!.md<com>    - EXTRA: extracellular potential. <b>Probably not implemented</b>
!.md<com>    - POTAS: potassium
!.md<com>    - EXTRA: extracellular potential
!.md<com>    - RECOV: ???
!.md<com>    - FIBER: fibers
!.md<com>    - SHEET: Sheet fiber direction
!.md<com>    - NORMA: Normal fiber direction
!.md<com>    - FREE4: ???
!.md<com>    - FREE5: ???
!.md<com>    - ISOCH: Isochrones (depolarisation time)
!.md<com>    - REPOL: Repolarisation time
!.md<com>    - CECO1:
!.md<com>    - CECO2:
!.md<com>    - CECO3:
!.md<com>    - CURRE: Currents ???
!.md<com>    - QNETO: integral of patassium currents (by FDA)
!.md<com>    - VCO0x: ???
!.md<com>    - VCOx : ???
!.md<com>    - VAU0x: ???
!.md<com>    - VAUx : ??? 
!.md<com>    - GRAFI: <b>TODO: DESCIBE</b>
!.md<com>    - IOCON: Save ion concentrations, a ,multidimensional field on nodes
!.md<com>
!
!.md<1>ELEMENT_SET
!.md<1>  XXXXX                                                                                         $ Variables computed in volume integrals
!.md<1>END_ELEMENT_SET
!.md<field>ELEMENT_SET
!.md<com>Variables computed as volume integrals on element sets. The results are stored
!.md<com>in file <tt>sample-element.nsi.set</tt>. Let V be the size of the bounday set. 
!.md<com>The following variables can be postprocessed:
!.md<com>
!.md<com>    - INTRA: intracellular potential
!.md<com>    - CALCI: Calcium
!.md<com>    - SODIU: Sodium 
!.md<com>    - POTAS: Potassium  
!.md<com>
!
!.md<1>NODE_SET
!.md<1>  XXXXX                                                                                         $ Variables computed on nodes
!.md<1>END_NODE_SET
!.md<field>NODE_SET
!.md<com>Variables computed as boundary integrals on boundary sets. The results are stored
!.md<com>in file <tt>sample-node.nsi.set</tt>. 
!.md<com>The following variables can be postprocessed:
!.md<com>
!.md<com>    - INTRA: intracellular potential
!.md<com>    - CALCI: Calcium
!.md<com>    - SODIU: Sodium 
!.md<com>    - POTAS: Potassium  
!.md<com>
!.md<com> 
!.md<1>WITNESS_POINTS
!.md<1>  XXXXX                                                                                         $ Variables computed witness points
!.md<1>END_WITNESS_POINTS
!.md<field>WITNESS_POINTS
!.md<com>
!.md<com>    - INTRA 
!.md<com>    - COORX
!.md<com>    - CAI    : The next 11 variables are: VCONC(1:11)
!.md<com>    - NASS 
!.md<com>    - KI   
!.md<com>    - KISS 
!.md<com>    - NAI  
!.md<com>    - CASS 
!.md<com>    - CANSR
!.md<com>    - CAJSR
!.md<com>    - JRLNP  : Save Jrelnp
!.md<com>    - JRLP   : Save Jrelp
!.md<com>    - CAMKT
!.md<com>    - INA    : The next variables are VICEL(1:26) 
!.md<com>    - INAL 
!.md<com>    - ITO  
!.md<com>    - ICAL 
!.md<com>    - IKR  
!.md<com>    - IKS  
!.md<com>    - IK1  
!.md<com>    - INCI   : Save INaCa_i
!.md<com>    - INCSS  : Save INaCa_ss
!.md<com>    - INAK 
!.md<com>    - IKB  
!.md<com>    - INAB 
!.md<com>    - ICAB 
!.md<com>    - IPCA 
!.md<com>    - JDI    : Save Jdiff
!.md<com>    - JDINA  : Save JdiffNa
!.md<com>    - JDIK 
!.md<com>    - JUP 
!.md<com>    - JLEAK
!.md<com>    - JTR  
!.md<com>    - JREL 
!.md<com>    - CAMKA
!.md<com>    - STIM 
!.md<com>    - ICANA
!.md<com>    - ICAK 
!.md<com>    - CAMKB
!.md<com>
!.md<com>
!.md<1>SAVE_CONVERGENCE CELL_MODEL $ saves ohara.mXXcXX files with curves for different concentrations during the ODE iterations for each material and cell type
!.md<1>DUMP_INITIAL_CONDITIONS     $ saves initial conditions to be hardcoded for a cell model concentrations during the ODE iterations for each material and cell type
!.md<0>END_OUTPUT_&_POST_PROCESS
!.md</code>
