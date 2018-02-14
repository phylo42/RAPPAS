/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import java.io.File;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public class TEST_Main_v2 {
    public static void main(String[] args) {

            ///////////////////////////////////////////////////////////////////
            //TEST ZONE, forces arguments
            String HOME = System.getenv("HOME");
            
            //DATASET 4BRANCHES TESTS UNROOTED(!) --PAML--
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/4leaves_tree_benchmark/viromeplacer_unrooted";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/4leaves_tree_benchmark";
//            String a=inputsPath+File.separator+"basic.aln";
//            String t=inputsPath+File.separator+"RAxML_bestTree.basic_tree";
            
            //DATASET 4BRANCHES TESTS FORCE ROOTING(!) --PAML--
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/4leaves_tree_benchmark/viromeplacer_force_rooting";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/4leaves_tree_benchmark";
//            String a=inputsPath+File.separator+"basic.aln";
//            String t=inputsPath+File.separator+"RAxML_bestTree.basic_tree";
            
            //DATASET BASIC RAPID TESTS  --PAML--
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD2";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
//            String a=inputsPath+"mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
//            String t=inputsPath+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";
            
            //DATASET LARGER SET --PAML--
            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML";
            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL";
            String a=inputsPath+File.separator+"bv_refs_aln_stripped_99.5.fasta";
            String t=inputsPath+File.separator+"RAxML_result.bv_refs_aln";
            String arDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/AR";
            String exTree=HOME+"/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/extended_trees";

            //--------------------------------

            //DATASET BASIC RAPID TESTS --PHYML--
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD_SMALL_PHYML";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
//            String a=inputsPath+"mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
//            String t=inputsPath+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";


            //DATASET LARGER SET --PHYML--
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD_LARGE_PHYML";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL";
//            String a=inputsPath+File.separator+"bv_refs_aln_stripped_99.5.fasta";
//            String t=inputsPath+File.separator+"RAxML_result.bv_refs_aln";
            

            //QUERIES::
            
            //DATASET 4BRANCHES TESTS
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/4leaves_tree_benchmark/queries.fasta";
//            String db=workDir+File.separator+"DB_session_k5_a1.0_t9.765625E-4.full";
            
            
            //DATASET BASIC RAPID TESTS:
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/alphaTest1";
//            String db=HOME+"/Dropbox/viromeplacer/test_datasets/WD2/DB_session_k8_a1.5_t3.9106607E-4.full";

            //pplacer benchmark queries 
            String q=inputsPath+File.separator+"mod_p4z1r36_query_only2.fasta";
//          String q=inputsPath+"mod_p4z1r36_query_1st_seq_expanded.fasta";
//          String q=inputsPath+"mod_p4z1r36_query_ancestrals.fasta";
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/mod_p4z1r36_query_only2.fasta";
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/mod_2VGB.qc.fasta";
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/mod_2VGB.qc.fasta_1000000_reads_only.fasta";
            //String q=HOME+"/Dropbox/viromeplacer/test_datasets/mod_2VGB_100000.qc.fasta";
            //String db=workDir+File.separator+"DB_session_k8_a1.5_t3.9106607E-4.medium";
            //String db=workDir+File.separator+"DB_session_k8_a1.5_t3.9106607E-4.small";
            //String db=workDir+File.separator+"DB_session_k8_a1.5_f1_t-3.40775.union";
            
            //String db=workDir+File.separator+"DB_session_k10_a0.75_f1_t-7.269987.union";
            String db=workDir+File.separator+"DB_session_k10_a1.0_f1_t-6.0206.union";
            //String db=workDir+File.separator+"DB_session_k8_a1.0_f1_t-4.81648.union";
            
//            String q="/home/ben/Downloads/R5_nx648_la_r150.fasta";
//            String db=workDir+File.separator+"DB_session_k5_a1.0_t9.765625E-4.medium";

            //db build launch
            String arguments=
                              "-m B "
                            + "-w "+workDir+" "
                            + "-r "+a+" "
                            + "-t "+t+" "
                            + "-k "+String.valueOf(10)+" "
                            + "-a "+String.valueOf(1.0)+" "
                            + "-v 1 "
                            + "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/phyml-AR/bin/phyml "
                            + "--ardir "+arDir+" "
                            //+ "--extree "+exTree+" "
                            //+ "--forceroot"
                            //+ "--dbinram "
//                            + "-q "+q+" "
                            + "--nsbound -10000000.0 "
                            //+ "--do-gap-jumps "
                            //+ "--only-1-jump"
                            ;
            
            
            
            
            // placement launch
//            arguments=
//                              "-m p "
//                            + "-w "+workDir+" "
//                            + "-q "+q+" "
//                            + "-d "+db+" "
//                            + "-v 0 "
//                            //+ "--nsbound -1000.0 "  -Infinity if not given
//                            + "--keep-at-most 7 "
//                            + "--keep-factor 0.01"
//                            ;            
                            
                            
/////////FOR TESTS CORRECTION: HCV, A34, k8_a1.0, R500bp  --> case where our algo cannot fails compared to EPA/PPlacer          
//            arguments=
//                              "-m B "
//                            + "-w /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0 "
//                            + "-r /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Ax/A34_nx34_la.align "
//                            + "-t /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Tx/T34_nx34_la.tree "
//                            + "-k "+String.valueOf(8)+" "
//                            + "-a "+String.valueOf(1.0)+" "
//                            + "-v 1 "
//                            + "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/baseml "
//                            + "--ardir /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/AR "
//                            //+ "--extree "+exTree+" "
//                            //+ "--builddbfull "
//                            + "--froot "
//                            + "--dbinram "
//                            + "-q /home/ben/Desktop/k8_a1.1/R33_nx348_la_r300.fasta "
//                            + "--nsbound -100000.0 "
//                            + "--nocalib "
//                            + "--unihash"
//                            ;     
//            arguments=
//                              "-m p "
//                            + "-w /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0 "
//                            + "-q /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Rx/R34_nx34_la_r500.fasta "
//                            //+ "-d /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0/DB_session_k8_a1.0_f1_t-4.81648.medium "
//                            //+ "-d /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0/DB_session_k8_a1.0_f1_t-4.81648.union "
//                            //+ "-d /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0/DB_session_k11_a1.0_f1_t-6.6226597.medium "
//                            + "-d /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0/DB_session_k11_a0.5_f1_t-6.6226597.union "
//                            //+ "-d /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV/Dx/A34_nx34_la/k8_a1.0/DB_session_k11_a1.0_f1_t-6.6226597.union "
//                            + "-v 1 "
//                            + "--nsbound -100000.0"
//                            ;    
/////////FOR DB SMALL CORRECTION: HCV, A33, k8_a1.1, R300bp              
                            

/////////FOR TESTS OF DB BUILD SPEED: HCV, A34, k12_a0.75 --> case where our RAPPAS seems to slow down with number of nodes explored       
//            arguments=
//                              "-m B "
//                            + "-w /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV_sunion/Dx/A34_nx34_la/k12_a0.75 "
//                            + "-r /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV_sunion/Ax/A34_nx34_la.align "
//                            + "-t /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV_sunion/Tx/T34_nx34_la.tree "
//                            + "-k "+String.valueOf(11)+" "
//                            + "-a "+String.valueOf(0.75)+" "
//                            + "-v 1 "
//                            + "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/baseml "
//                            + "--ardir /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/HCV_sunion/Dx/A34_nx34_la/AR "
//                            //+ "--extree "+exTree+" "
//                            //+ "--builddbfull "
//                            + "--froot "
//                            //+ "--dbinram "
//                            //+ "-q /home/ben/Desktop/k8_a1.1/R33_nx348_la_r300.fasta "
//                            + "--nsbound -100000.0 "
//                            + "--nocalib "
//                            + "--unihash"
//                            ;     
/////////FOR TESTS OF DB BUILD SPEED: HCV, A34, k12_a0.75 --> case where our RAPPAS seems to slow down with number of nodes explored       


/////////FOR TESTS OF DB BUILD ERROR: BOLD_chordata, A30, k6_a1.50 --> case where our RAPPAS produces "null" value for likelihood_weight_ratio
//            arguments=
//                              "-m B "
//                            + "-w /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Dx/A30_nx759_la/k6_a1.50 "
//                            + "-r /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Ax/A30_nx759_la.align "
//                            + "-t /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Tx/T30_nx759_la.tree "
//                            + "-k "+String.valueOf(6)+" "
//                            + "-a "+String.valueOf(1.50)+" "
//                            + "-v 1 "
//                            + "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/phyml-AR/bin/phyml "
//                            + "--ardir /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Dx/A30_nx759_la/AR "
//                            + "--dbinram "
//                            + "--no-reduction "
//                            + "--nsbound -10000000.0 "
//                            + "-q /media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Rx/R30_nx759_la_r150.fasta"
                            //+    ",/media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Rx/R30_nx759_la_r300.fasta"
                            //+    ",/media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Rx/R30_nx759_la_r600.fasta"
                            //+    ",/media/ben/STOCK/DATA/viromeplacer/accu_tests/imports/BOLD_chordata/Rx/R30_nx759_la_r1200.fasta"
                            ;     
/////////FOR TESTS OF DB BUILD ERROR: BOLD_chordata, A30, k6_a1.50 --> case where our RAPPAS produces "null" value for likelihood_weight_ratio




                            
//  FOR PROTEIN ANALYSIS TESTS
/////////////////////////////////
//
            workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD_PROT_PAML";
            inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml_prot";
            a=inputsPath+File.separator+"cox1_euka_oma.mfa";
            t=inputsPath+File.separator+"RAxML_bestTree.test";
            arDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD_PROT_PAML/AR";
            //String exTree=HOME+"/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/extended_trees";

            q=inputsPath+File.separator+"queries.fasta";
            db=HOME+"/Dropbox/viromeplacer/test_datasets/WD_PROT_PAML/DB_session_k6_a1.0_f1_t-3.61236.union";

            //db build launch
            arguments=
                              "-m B "
                            + "-s amino "
                            + "-w "+workDir+" "
                            + "-r "+a+" "
                            + "-t "+t+" "
                            + "-k "+String.valueOf(3)+" "
                            + "-a "+String.valueOf(1.0)+" "
                            + "-v 1 "
                            //+ "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/codeml "
                            + "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/phyml-AR/bin/phyml "
                            + "--ardir "+arDir+" "
                            //+ "--extree "+exTree+" "
                            //+ "--dbfull "
                            //+ "--froot"
                            //+ "--dbinram "
//                            + "-q "+q+" "
//                            + "--nsbound -100000.0 "
//                            + "--nocalib"
                            ;
//
////            // placement launch
//            arguments=
//                              "-m p "
//                            + "-s amino "
//                            + "-w "+workDir+" "
//                            + "-q "+q+" "
//                            + "-d "+db+" "
//                            + "-v 0 "
//                            + "--nsbound -1000.0"
//                            ;   
//
//

//////////////////////////////////


//  FOR PROTEIN ANALYSIS TESTS, large dataset
///////////////////////////////////////////////
//
            workDir=HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/D140";
            inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/D140";
            a=inputsPath+File.separator+"mod_140.phy.fasta";
            t=inputsPath+File.separator+"RAxML_bipartitions.140.BEST.WITH";
            arDir=HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/D140/AR";
            //String exTree=HOME+"/Dropbox/viromeplacer/test_datasets/WD_LARGE_PAML/extended_trees";

            q=inputsPath+File.separator+"R2_nx5_la20_r300.fasta";
            db=HOME+"/Dropbox/viromeplacer/test_datasets/accuracy_tests/DATA/D140/DB_session_k3_a1.0_f1_t-3.90309.union";

            //db build launch
            arguments=
                              "-m B "
                            + "-s amino "
                            + "-w "+workDir+" "
                            + "-r "+a+" "
                            + "-t "+t+" "
                            + "-k "+String.valueOf(6)+" "
                            + "-a "+String.valueOf(5.0)+" "
                            + "-v 1 "
                            //+ "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/codeml "
                            + "--arbinary "+HOME+"/Dropbox/viromeplacer/test_datasets/software/phyml-AR/bin/phyml "
                            + "--ardir "+arDir+" "
                            //+ "--extree "+exTree+" "
                            //+ "--dbfull "
                            //+ "--froot"
                            //+ "--dbinram "
//                            + "-q "+q+" "
//                            + "--nsbound -100000.0 "
//                            + "--nocalib"
                            ;
//
////            // placement launch
            arguments=
                              "-m p "
                            + "-s amino "
                            + "-w "+workDir+" "
                            + "-q "+q+" "
                            + "-d "+db+" "
                            + "-v 0 "
                            + "--nsbound -100000.0"
                            ;   



//////////////////////////////////


    
           
            //TEST ZONE//
            System.out.println(Arrays.toString(args));
            Main_v2.main(arguments.split(" "));
            
            
    }
}
