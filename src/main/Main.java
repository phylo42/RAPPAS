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
import static main.Main_DBBUILD.TYPE_DNA;



/**
 *
 * @author linard
 */
public class Main {

    private final static String consoleVersion="0.1";

    public static void main (String[] args) {
        try {
            long startTime=System.currentTimeMillis();
            System.out.println("#####################################");
            System.out.println("## Viromeplacer v"+consoleVersion);
            System.out.println("#####################################");
            
            System.setProperty("viromeplacer_version", consoleVersion);
            
            
            
            
            ///////////////////////////////////////////////////////////////////
            //TEST ZONE, forces arguments
            
            String workDir="/media/ben/STOCK/DATA/viromeplacer/WD";
            
            //here,pplacer benchmark to build DB
            String inputsPath="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/";
//            String a=wd+"bv_refs_aln.fasta";
            String a=inputsPath+"bv_refs_aln_stripped_99.5.fasta";
            String t=inputsPath+"RAxML_result.bv_refs_aln";

//            a=inputsPath+"bv_refs_aln_stripped_99.5_SMALL_SUBSET.fasta";
//            t=inputsPath+"RAxML_result.bv_refs_aln_SMALL_SUBSET.tree";

            //pplacer benchmark queries 
            String q=inputsPath+"mod_p4z1r36_query_only2.fasta";
            //String q=inputsPath+"mod_p4z1r36_query_1st_seq_expanded.fasta";
            //String q=inputsPath+"mod_p4z1r36_query_ancestrals.fasta";
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
//            String arguments=
//                              "-m p "
//                            + "-w "+workDir+" "
//                            + "-q "+q+" "
//                            + "-d "+db+" "
//                            + "-v 0"
//                            ;            
//            
//            
//            
//            
//            args=arguments.split(" ");
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
            //System.out.println(argsParser);
            
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
                int placed=Main_PLACEMENT_V03_align_scoreallnodes.Main_PLACEMENT_V03_align_scoreallnodes(
                                            argsParser.queriesFile,
                                            argsParser.databaseFile,
                                            argsParser.workingDir
                                            );
                System.out.println("# queries placed: "+placed);
            }
            
            
            long endTime=System.currentTimeMillis();
            System.out.println("TOTAL EXECUTION TIME: "+(endTime-startTime)+" ms");
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

}
