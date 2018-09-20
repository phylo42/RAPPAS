/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main_v2;

import core.AAStates;
import core.DNAStatesShifted;
import core.States;
import etc.NullPrintStream;
import java.io.File;
import javax.swing.UIManager;
import models.EvolModel;


/**
 *
 * @author linard
 */
public class Main_v2 {

    private final static String consoleVersion="1.05";

    public static void main (String[] args) {
        try {
            long startTime=System.currentTimeMillis();
            
	    //System.out.println(VM.current().details());
            System.setProperty("viromeplacer_version", consoleVersion);
            
            //hack related to Problems under MAC OS implementation of
            //the Aqua (mac Look and feel) in some implementations,
            //(Mac implementation of Java is screwed... not Oracle or openJDK)
            //for unknown reason, the use of Jtree prompts the virtual machine to 
            //use the class com.apple.laf.AquaTreeUI
            //which is not Serializable and causes crashes when Jtree is serialized
            //(i.e. when the database is saved on the disk)
            // Set cross-platform Java L&F as defaults (also called "Metal")
            UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            
            //parse program arguments
            ArgumentsParser_v2 argsParser = new ArgumentsParser_v2(args,consoleVersion);
            
            
            //set verbosity to null, if required
            if (argsParser.verbose<0) {
                //PrintStream original = System.out;
                System.setOut(new NullPrintStream());
                //System.out.println("Message not shown.");
                //System.setOut(original); //to get back output stream
            }
            
            System.out.println("################################################");
            System.out.println("## RAPPAS v"+consoleVersion);
            System.out.println("## ---------------------------------------------");
            System.out.println("## Rapid Alignment-free Phylogenetic Placement ");
            System.out.println("## via Ancestral Sequences");
            System.out.println("## Linard B, Swenson KM, Pardi F");
            System.out.println("## LIRMM, Univ. of Montpellier, CNRS, France");
            System.out.println("## https://doi.org/10.1101/328740");
            System.out.println("## benjamin/dot/linard/at/lirmm/dot/fr");
            System.out.println("################################################");
            
            
            
            //type of Analysis, DNA or AA
            States s=null; 
            if (argsParser.states==ArgumentsParser_v2.STATES_DNA) {
                s=new DNAStatesShifted();
                System.out.println("Set analysis for DNA");
            } else if (argsParser.states==ArgumentsParser_v2.STATES_PROTEIN) {
                s=new AAStates(argsParser.convertUOX);
                System.out.println("Set analysis for PROTEIN");
            }
            
            

            
            
            //////////////////////
            //DB_BUILD MODE
            
            if (argsParser.phase==ArgumentsParser_v2.DBBUILD_PHASE) {
                System.out.println("Starting db_build pipeline...");
                
                //set default model for ASR if user set nothing
                EvolModel model=null;
                if ( argsParser.modelString==null ) {
                    model=EvolModel.getDefault(argsParser.states);
                    System.out.println("User did not set model parameters, using default: "+model.toString());
                } else {
                    model=new EvolModel(
                        argsParser.states,
                        argsParser.modelString, 
                        argsParser.alpha, 
                        argsParser.categories
                    );
                    System.out.println("Model parameters: "+model.toString());
                }
                //test if model compatible to nucl/prot, test not done is user provides --arparameters
                if ( !model.isProteinModel() && (s instanceof AAStates) && (argsParser.arparameters==null)) {
                    System.out.println("You ask for a nucleic model but use amino acid states, please correct (see -h).");
                    System.exit(1);
                } else if (model.isProteinModel() && !(s instanceof AAStates) && (argsParser.arparameters==null)) {
                    System.out.println("You ask for a proteic model but use nucleic states, please correct (see -h).");
                    System.exit(1);
                }


                Main_DBBUILD_3.DBGeneration(
                                            null,
                                            argsParser.k,
                                            argsParser.omega,
                                            argsParser.ghostsAmount,
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
                                            argsParser.gapJumpThreshold,
                                            model,
                                            argsParser.arparameters
                        
                                            );
                System.out.println("Have a coffee, you \"built\" your world.");

                
            //////////////////////
            //PLACEMENT MODE
            } else if (argsParser.phase==ArgumentsParser_v2.PLACEMENT_PHASE) {
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
                                                argsParser.keepFactor,
                                                argsParser.guppyCompatible
                                                );
                }
                System.out.println("Have a coffee, you \"placed\" your world.");

            }
            
            
            long endTime=System.currentTimeMillis();
            System.out.println("Total execution time: "+(endTime-startTime)+" ms ("+((endTime-startTime)/60000)+" min)");
            //System.exit(0);
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

}
