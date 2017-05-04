/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main;

import core.AAStates;
import core.DNAStates;
import core.States;
import etc.Infos;
import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import static main.Main_DBBUILD.TYPE_DNA;



/**
 *
 * @author linard
 */
public class Main {

    private final static String consoleVersion="0.2";

    public static void main (String[] args) {
        try {
            long startTime=System.currentTimeMillis();
            System.out.println("#####################################");
            System.out.println("## Viromeplacer v"+consoleVersion);
            System.out.println("#####################################");
            
            System.setProperty("viromeplacer_version", consoleVersion);
            
            
            
            
            ///////////////////////////////////////////////////////////////////
            //TEST ZONE, forces arguments
            String HOME = System.getenv("HOME");

            //DATASET BASIC RAPID TESTS:
//            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD2";
//            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
//            String a=inputsPath+"mod_mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
//            String t=inputsPath+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";
//            
            //DATASET LARGER SET:
            String workDir=HOME+"/Dropbox/viromeplacer/test_datasets/WD";
            String inputsPath=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL";
            String a=inputsPath+"bv_refs_aln_stripped_99.5.fasta";
            String t=inputsPath+"RAxML_result.bv_refs_aln";


            //QUERIES::
            
            //DATASET BASIC RAPID TESTS:
//            String q=HOME+"/Dropbox/viromeplacer/test_datasets/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/alphaTest1";
//            String db=HOME+"/Dropbox/viromeplacer/test_datasets/WD2/PAML_session_params_k8_mk8_f1.5_t3.9106607E-4";

            //pplacer benchmark queries 
            //String q=inputsPath+File.separator+"mod_p4z1r36_query_only2.fasta";
            String q="/home/benclaff/mod_2VGB.qc.fasta";
//          String q=inputsPath+"mod_p4z1r36_query_1st_seq_expanded.fasta";
//          String q=inputsPath+"mod_p4z1r36_query_ancestrals.fasta";
            String db=workDir+File.separator+"PAML_session_params_k8_mk8_f1.5_t3.9106607E-4";



            //db build launch
//            String arguments=
//                              "-m b "
//                            + "-w "+workDir+" "
//                            + "-i "+a+" "
//                            + "-t "+t+" "
//                            + "-k "+String.valueOf(8)+" "
//                            + "-a "+String.valueOf(1.5)+" "
//                            + "-v 1"
//                            ;
            
            // placement launch
            String arguments=
                              "-m p "
                            + "-w "+workDir+" "
                            + "-q "+q+" "
                            + "-d "+db+" "
                            + "-v 0"
                            ;            
            
            
            
          
            args=arguments.split(" ");
            
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
            ArgumentsParser argsParser = new ArgumentsParser(args);
            argsParser.pamlPath=new File("/media/ben/STOCK/SOFTWARE/paml4.9b_hacked/bin/baseml");
            
            
            if (argsParser.mode==ArgumentsParser.DBBUILD_MODE) {
                System.out.println("Starting db_build pipeline...");
                Main_DBBUILD.DBGeneration(  null,
                                            argsParser.k,
                                            argsParser.alpha,
                                            argsParser.fakeBranchAmount,
                                            s,
                                            argsParser.alignmentFile,
                                            argsParser.treeFile,
                                            argsParser.workingDir,
                                            argsParser.pamlPath
                                            );
            } else if (argsParser.mode==ArgumentsParser.PLACEMENT_MODE) {
                System.out.println("Starting placement pipeline...");
                int placed=Main_PLACEMENT_V04_align_scoreallnodes_diagsumeltsremoved.Main_PLACEMENT_V03_align_scoreallnodes(
                                            argsParser.queriesFile,
                                            argsParser.databaseFile,
                                            argsParser.workingDir
                                            );
                System.out.println("# query reads placed: "+placed);
            }
            
            
            long endTime=System.currentTimeMillis();
            System.out.println("Total execution: "+(endTime-startTime)+" ms");
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

}
