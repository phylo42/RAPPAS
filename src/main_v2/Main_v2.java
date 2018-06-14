/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main_v2;

import core.AAStates;
import core.DNAStatesShifted;
import core.DNAStates;
import core.States;
import etc.NullPrintStream;
import java.io.File;
import javax.swing.UIManager;


/**
 *
 * @author linard
 */
public class Main_v2 {

    private final static String consoleVersion="1.02";

    public static void main (String[] args) {
        try {
            long startTime=System.currentTimeMillis();
            System.out.println("################################################");
            System.out.println("## Rapid Alignment-free Phylogenetic Placement ");
            System.out.println("## via Ancestral Sequences");
            System.out.println("## RAPPAS v"+consoleVersion);
            System.out.println("## benjamin.linard, fabio.pardi ([at].lirmm.fr)");
            System.out.println("## LIRMM, Univ. of Montpellier, CNRS");
            System.out.println("################################################");
            //System.out.println(VM.current().details());
            System.setProperty("viromeplacer_version", consoleVersion);
            
            //hack related to Problems under MAC OS implementation of
            //the Aqua (mac Look and feel)
            //in some implementations, (Mac implementation if Java... not Oracle or open JDK)
            //for unknown reason, the use of Jtree prompts the virtual machine to 
            //use the class com.apple.laf.AquaTreeUI
            //which is not Serializable and causes crashes when Jtree is serialized
            //(i.e. when the database is saved on the disk)
            // Set cross-platform Java L&F (also called "Metal")
            UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            
            //parse program arguments
            ArgumentsParser_v2 argsParser = new ArgumentsParser_v2(args);
            
            
            //set verbosity to null, if required
            if (argsParser.verbose<0) {
                //PrintStream original = System.out;
                System.setOut(new NullPrintStream());
                //System.out.println("Message not shown.");
                //System.setOut(original); //to get back output stream
            }
            
            
            //type of Analysis, DNA or AA
            States s=null; 
            if (argsParser.analysisType==ArgumentsParser_v2.TYPE_DNA) {
                //s=new DNAStatesShifted();
                s=new DNAStates();
                System.out.println("Set analysis for DNA");
            } else if (argsParser.analysisType==ArgumentsParser_v2.TYPE_PROTEIN) {
                s=new AAStates();
                System.out.println("Set analysis for PROTEIN");
            }

            
            
            //////////////////////
            //DB_BUILD MODE
            
            if (argsParser.mode==ArgumentsParser_v2.DBBUILD_MODE) {
                System.out.println("Starting db_build pipeline...");
                



                Main_DBBUILD_3.DBGeneration(argsParser.analysisType,
                //Main_DBBUILD_HashTriplet.DBGeneration(argsParser.analysisType,
                                            null,
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
                                            argsParser.builddbfull,
                                            argsParser.forceRooting,
                                            argsParser.dbInRAM,
                                            argsParser.queriesFiles,
                                            argsParser.callString,
                                            argsParser.nsBound,
                                            argsParser.noCalibration,
                                            argsParser.unionHash,
                                            argsParser.reduction,
                                            argsParser.reducedAlignFile,
                                            argsParser.reductionRatio,
                                            argsParser.onlyFakeNodes,
                                            argsParser.keepAtMost,
                                            argsParser.keepFactor,
                                            argsParser.doGapJumps,
                                            argsParser.limitTo1Jump,
                                            argsParser.gapJumpThreshold
                        
                                            );
                
            //////////////////////
            //PLACEMENT MODE
                
            } else if (argsParser.mode==ArgumentsParser_v2.PLACEMENT_MODE) {
                System.out.println("Starting placement pipeline...");
                //load session itself (i.e the DB)
                System.out.println("Loading ancestral words DB... ("+argsParser.databaseFile.getName()+")");
                long startLoadTime=System.currentTimeMillis();
                SessionNext_v2 session= SessionNext_v2.load(argsParser.databaseFile,true);
                long endLoadTime=System.currentTimeMillis();
                System.out.println("Loading the database took "+(endLoadTime-startLoadTime)+" ms");
            
                Main_PLACEMENT_v07 placer=new Main_PLACEMENT_v07(session, false);
                for (int i = 0; i < argsParser.queriesFiles.size(); i++) {
                    File query = argsParser.queriesFiles.get(i);
                    int placed=placer.doPlacements(query,
                                                argsParser.databaseFile,
                                                argsParser.workingDir,
                                                argsParser.callString,
                                                argsParser.nsBound,
                                                argsParser.keepAtMost,
                                                argsParser.keepFactor
                                                );
                }
                
            }
            
            
            long endTime=System.currentTimeMillis();
            System.out.println("Total execution time: "+(endTime-startTime)+" ms ("+((endTime-startTime)/60000)+" min)");
            System.out.println("Have a coffee, you \"placed\" your world.");
            //System.exit(0);
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

}
