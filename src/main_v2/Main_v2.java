/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main_v2;

import core.AAStates;
import core.DNAStates;
import core.States;
import java.io.File;
import static main.Main_DBBUILD.TYPE_DNA;


/**
 *
 * @author linard
 */
public class Main_v2 {

    private final static String consoleVersion="0.7";

    public static void main (String[] args) {
        try {
            long startTime=System.currentTimeMillis();
            System.out.println("#############################################################");
            System.out.println("## RApid Phylogenetic Placement via Ancestral Sequences");
            System.out.println("## (RAPPAS) v"+consoleVersion);
            System.out.println("## Author: benjamin.linard (LIRMM-CNRS, Montpellier, France)");
            System.out.println("#############################################################");
            //System.out.println(VM.current().details());
            System.setProperty("viromeplacer_version", consoleVersion);
            
            
            
            
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
//            String q=inputsPath+File.separator+"mod_p4z1r36_query_only2.fasta";
//          String q=inputsPath+"mod_p4z1r36_query_1st_seq_expanded.fasta";
//          String q=inputsPath+"mod_p4z1r36_query_ancestrals.fasta";
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/mod_p4z1r36_query_only2.fasta";

            String q=HOME+"/Dropbox/viromeplacer/test_datasets/mod_2VGB.qc.fasta";
            String db=workDir+File.separator+"DB_session_k8_a1.5_t3.9106607E-4.full";
            
//            String q="/home/ben/Downloads/R5_nx648_la_r150.fasta";
//            String db=workDir+File.separator+"DB_session_k5_a1.0_t9.765625E-4.medium";

            //db build launch
//            String arguments=
//                              "-m B "
//                            + "-w "+workDir+" "
//                            + "-i "+a+" "
//                            + "-t "+t+" "
//                            + "-k "+String.valueOf(5)+" "
//                            + "-a "+String.valueOf(1.0)+" "
//                            + "-v 1 "
//                            + "--ardir "+arDir+" "
//                            + "--extree "+exTree+" "
//                            + "--skipdbfull "
//                            + "--froot";
            
//            // placement launch
            String arguments=
                              "-m p "
                            + "-w "+workDir+" "
                            + "-q "+q+" "
                            + "-d "+db+" "
                            + "-s medium "
                            + "-v 0"
                            ;            
            
            //force args
            //args=arguments.split(" ");
            
            
            
            
            //System.out.println(Arrays.toString(args));
            
           
            //type of Analysis, currently not in command line parameters
            States s=null; 
            int analysisType=TYPE_DNA;
            //States: DNA or AA
            if (analysisType==TYPE_DNA)
                s=new DNAStates();   
            else if (analysisType==TYPE_DNA)
                s=new AAStates();
           
            //TEST ZONE//

            
            
            //parse program arguments
            ArgumentsParser_v2 argsParser = new ArgumentsParser_v2(args);
            //argsParser.ARBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/paml4.9b_hacked/bin/baseml");
            //argsParser.ARBinary=new File(HOME+"/Dropbox/viromeplacer/test_datasets/software/phyml/src/phyml");
            
            //HACK FOR CURRENT DEBUGING, avoids check if it exists or not (done by ArgumentsParser)
            argsParser.ARBinary=new File("baseml");
            
            
            if (argsParser.mode==ArgumentsParser_v2.DBBUILD_MODE) {
                System.out.println("Starting db_build pipeline...");

                
                Main_DBBUILD_2.DBGeneration(null,
                                            argsParser.k,
                                            argsParser.alpha,
                                            argsParser.fakeBranchAmount,
                                            s,
                                            argsParser.alignmentFile,
                                            argsParser.treeFile,
                                            argsParser.workingDir,
                                            argsParser.ARBinary,
                                            argsParser.ARDirToUse,
                                            argsParser.exTreeDir,
                                            argsParser.skipdbfull,
                                            argsParser.forceRooting
                                            );
                
            } else if (argsParser.mode==ArgumentsParser_v2.PLACEMENT_MODE) {
                System.out.println("Starting placement pipeline...");
//                int placed=Main_PLACEMENT_V05_align_scoring_separated_for_time_eval.Main_PLACEMENT_V05_align_scoreallnodes(
//                                            argsParser.queriesFile,
//                                            argsParser.databaseFile,
//                                            argsParser.workingDir
//                                            );
                int placed=Main_PLACEMENT_FIN.doPlacements(
                                            argsParser.queriesFile,
                                            argsParser.databaseFile,
                                            argsParser.workingDir,
                                            argsParser.callString
                                            );
                
                System.out.println("# query reads placed: "+placed);
            }
            
            
            long endTime=System.currentTimeMillis();
            System.out.println("Total execution: "+(endTime-startTime)+" ms");
            System.exit(0);
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

}
