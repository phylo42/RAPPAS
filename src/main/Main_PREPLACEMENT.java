/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import au.com.bytecode.opencsv.CSVWriter;
import core.DNAStates;
import core.PProbas;
import core.QueryWord;
import core.States;
import core.algos.ReadKnife;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PAMLWrapper;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.NewickParser;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.RelaxedTree;

/**
 *
 * @author ben
 */
public class Main_PREPLACEMENT {
    
    public static void main(String[] args) {

        try {
            System.setProperty("debug.verbose","1");
            
            
            // INPUT FILES//////////////////////////////////////////////////////
            
            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/";
//            String a=wd+"bv_refs_aln.fasta";
            String a=wd+"bv_refs_aln_stripped_99.5.fasta";
            String t=wd+"RAxML_result.bv_refs_aln";
            String q=wd+"mod_p4z1r36_query_only.fasta";
            String pp=wd+"rst";
            
//            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
//            String a=wd+"mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
//            String t=wd+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";
//            String q=wd+"alphaTest1";
//            String pp=wd+"rst";
            
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            //build of relaxed tree/////////////////////////////////////////////
            double minBranchLength=0.001;
            int numberOfFakeBranchesPerEdge=1;
            
            //posterior probas analysis/////////////////////////////////////////
            //mers size
            int k=8;
            int min_k=8;
            //site and word posterior probas thresholds
            double sitePPThreshold=1e-6;
            double wordPPThreshold=1e-6;
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            
            //pre-treatment parameters//////////////////////////////////////////
            //for a number of nodes =total_nodes/nodeShift , build diagsum vectors
            //think to keep sum of the diagsums to do mean at the end and highlight positions > to mean
            int nodeShift=200;
            //peek detection parameters
            int w=31; //sliding window size
            double r_win=0.50; //peek/windows_mean ratio
            double r_nodes=0.50; //node_where_peek_detected/scanned_nodes ratio
            int maxAllowedPeeks=50; //maximum number of localizations allowed for a read, if greater read discarded because clearly a repeat
            //States: DNA or AA
            States s=new DNAStates();    
            
            //debug/////////////////////////////////////////////////////////////
            int queryLimit=10000;
            boolean logDiagsums=false;
            boolean logDetailedDiagsums=false;
            
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            
            //////////////////////
            //LOAD ORIGINAL ALIGNMENT
            FASTAPointer fp=new FASTAPointer(new File(a), false);
            Fasta f=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((f=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(f);
            }
            Alignment align=new Alignment(fastas);
            Infos.println(align.describeAlignment(false));
            fp.closePointer();
            
            /////////////////////
            //PARSE ORIGINAL TREE
            NewickParser np=new NewickParser();
            PhyloTree tree = np.parseNewickTree(new File(t));
            
            /////////////////////
            //BUILD RELAXED TREE
//            RelaxedTree relaxedTreeOnBranches=new RelaxedTree(tree,minBranchLength,fakeBranchPerEdge);
//            ArrayList<PhyloNode> listOfNewFakeLeaves = relaxedTreeOnBranches.getListOfNewFakeLeaves();
//            //add new leaves to alignment
//            for (int i = 0; i < listOfNewFakeLeaves.size(); i++) {
//                PhyloNode node = listOfNewFakeLeaves.get(i);
//                char[] gapSeq=new char[align.getLength()];
//                Arrays.fill(gapSeq, '-');
//                align.addSequence(node.getLabel(), gapSeq);
//            }
//            //write alignment and tree for BrB
//            align.writeAlignmentAsFasta(new File(wd+"relaxed_align_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge.fasta"));
//            align.writeAlignmentAsPhylip(new File(wd+"relaxed_align_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge.phylip"));
//            np.writeNewickTree(relaxedTreeOnBranches, new File(wd+"relaxed_tree_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge.tree"),true,false);
//            np.writeNewickTree(relaxedTreeOnBranches, new File(wd+"relaxed_tree_BrB_minbl"+minBranchLength+"_"+fakeBranchPerEdge+"peredge_withBL.tree"),true,true);

            
            //////////////////////////////////////
            // HERE LAUNCH BASEML ON RELAXED TREE
            //todo
            //launch paml from java, can it be done whithout modifying
            //the ctl file but through ocmmand parameters ?
            
            
            
            //////////////////////////////////////////////
            //LOAD THE POSTERIOR PROBAS AND PAML TREE IDS
            Infos.println("Loading PAML tree ids and POsterior Probas...");
            PAMLWrapper pw=new PAMLWrapper(align, s);
            FileInputStream input = new FileInputStream(new File(pp));
            tree= pw.parseTree(input);
            input.close();
            input = new FileInputStream(new File(pp));
            PProbas pprobas = pw.parseProbas(input);
            input.close();

            
            
            /////////////////////
            //PLACEMENT OF THE QUERIES
            fp=new FASTAPointer(new File(q), false);
            FileWriter fw =new FileWriter(new File("Queries.fasta"));
            
            int queryCounter=0;
            while ((f=fp.nextSequenceAsFastaObject())!=null) {
                if (queryCounter>=queryLimit)
                    break;
                
                Infos.println("#######################################################################");
                Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+f.getHeader());
                Infos.println("#######################################################################");
                fw.append(f.getFormatedFasta()+"\n");
                long startScanTime=System.currentTimeMillis();
                
                int scannedNodesCounter=0;
                Infos.println("Number of internal nodes: "+tree.getInternalNodesByDFS().size());
                int scannedNodes= tree.getInternalNodesByDFS().size()/nodeShift ;
                Infos.println("Number of internal nodes to scan: "+scannedNodes);
                
                //query length
                int queryLength=f.getSequence().length();
                Infos.println("Query length: "+queryLength);
                //diagSums shift on the left
                int diagsumShift=queryLength-minOverlap;
                Infos.println("Diagsum shift: "+diagsumShift);
                //diagSums table
                int[] diagsumIndex = new int[scannedNodes]; //for instance if nodeShift=20 [20,40,60...]
                Arrays.fill(diagsumIndex,-1);
                //the diagSums vector has the length of the ref alignment for now
                //this is FOR NOW taking into account partial overlaps on left and right
                //but now when the read is a complementary sequence, these diagsums are ignored... for now.
                double[][] allDiagSums = new double[scannedNodes][align.getLength()+diagsumShift-minOverlap];
                Infos.println("Diagsum vector size: "+(align.getLength()+diagsumShift-minOverlap));
                
                
                /////////////////////////////
                // LOG OUTPUT 1:words detais
                CSVWriter writerDetails=null;
                if (logDiagsums) {
                    writerDetails=new CSVWriter(new FileWriter(new File("diagSums_details.tsv")), '\t');
                    String[] headerData=new String[2+align.getLength()];
                    headerData[0]="node";
                    headerData[1]="mer";
                    for (int i = 0; i < align.getLength(); i++) {
                        headerData[i+2]="pos="+i;
                    }
                    writerDetails.writeNext(headerData);
                }
                /////////////////
                
                
                /////////////////////////////////////
                //FOR EACH TESTED NODE
                for (int i=0;i<tree.getInternalNodesByDFS().size();i++) {
                    if (i%nodeShift!=0) {continue;}
                    if (scannedNodesCounter>scannedNodes-1) {break;} //all nodes done

                    int nodeId=tree.getInternalNodesByDFS().get(i);
                    PhyloNode node = tree.getById(nodeId);
                    //Infos.println("Building diagSum vector #"+scannedNodesCounter+" for node: "+node);
                    diagsumIndex[scannedNodesCounter]=nodeId; //nodes ids are from 0 to N
                    
                    ///////////////////////////////////
                    // LOOP ON QUERY WORDS
                    QueryWord qw=null;
                    ReadKnife rk=new ReadKnife(f, k, min_k, s, ReadKnife.SAMPLING_NON_OVERLAPPING);
                    //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                    int scannedMersCounter=0;
                    while ((qw=rk.getNextWord())!=null) {
                        //Infos.println("Query mer: "+qw.toString());
                        
                        /////////////////
                        // LOG OUTPUT 1
                        String[] data=new String[2+align.getLength()];
                        data[0]=String.valueOf(nodeId);
                        data[1]=String.valueOf(qw.getOriginalPosition());
                        /////////////////
                        
                        
                        //////////////////////////////////////
                        //BUILD DIAGSUMS, ONE PER NODE
                        for (int refPos=0;refPos<(align.getLength()-minOverlap);refPos++) {  
                            double PPStar=1.0;
                            for (int queryPos=0;queryPos<qw.getWord().length;queryPos++) {
                                PPStar*=pprobas.getPP(nodeId, refPos+queryPos, qw.getWord()[queryPos]);
                            }
//                            System.out.println("refPos:"+refPos);
//                            System.out.println("readLength-minOverlap:"+(readLength-minOverlap) );
//                            System.out.println("(refPos-qw.getOriginalPosition()):"+(refPos-qw.getOriginalPosition()));
//                            System.out.println("(readLength-minOverlap)+(refPos-qw.getOriginalPosition()):"+((readLength-minOverlap)+(refPos-qw.getOriginalPosition())));
                            
                            //the condition below is to avoid the last mers of 
                            //the query which would contradict the minOverlap condition:
                            //in the first bases of the alignment (pos<readlength)
                            //do not scan this word if it would represent an overlap 
                            //with the ref align < minOverlap
                            //this could be optimized to avoid some iterations ???
                            if (diagsumShift+(refPos-qw.getOriginalPosition())>-1) 
                                allDiagSums[scannedNodesCounter][diagsumShift+(refPos-qw.getOriginalPosition())]+=Math.log10(PPStar);
                            if (logDiagsums) {
                                data[2+refPos]=String.valueOf(Math.log10(PPStar)); 
                            }
                        }
                        if(logDiagsums)
                            writerDetails.writeNext(data);
                        
                        scannedMersCounter++;

                    }
                    
                    //Infos.println("# Mer scanned: "+scannedMersCounter);
                    scannedNodesCounter++;
                }
                if(logDiagsums) {
                    writerDetails.flush();
                    writerDetails.close();
                }
                
                
                long endScanTime=System.currentTimeMillis();
                Infos.println("Scan on "+(scannedNodesCounter+1)+" nodes took "+(endScanTime-startScanTime)+" ms");
                
                
                
                /////////////////
                // LOG OUTPUT 2
                CSVWriter writerDiagsum=null;
                if (logDetailedDiagsums) {
                    writerDiagsum=new CSVWriter(new FileWriter(new File("diagSums.tsv")), '\t');
                    String[] headerData=new String[2+scannedNodes];
                    headerData[0]="diagsum_position";
                    headerData[1]="align_site";
                    for (int i = 0; i < allDiagSums.length; i++) {
                        headerData[i+2]="nodeid="+diagsumIndex[i];
                    }
                    writerDiagsum.writeNext(headerData);

                    for (int i = 0; i < allDiagSums[0].length; i++) { //for each site
                        String[] data=new String[2+scannedNodes];
                        data[0]=String.valueOf(i);
                        data[1]=String.valueOf(i-diagsumShift);
                        for (int j = 0; j < scannedNodes; j++) {
                            data[j+2]=String.valueOf(allDiagSums[j][i]);
                        }
                        writerDiagsum.writeNext(data);
                    }
                    writerDiagsum.flush();
                    writerDiagsum.close();
                }
                /////////////////
                
                
                
                
                ///////////////////////////////////////////////
                //TO DO HERE: DETECTION OF PEEKS VIA RATIOS

                //loop to build PP* on PPstats table from 0 to align length-k.
                startScanTime=System.currentTimeMillis();
                double[][] ratios=new double[allDiagSums.length][allDiagSums[0].length];
                ArrayList<ArrayList<Integer>> detectedPeeks=new ArrayList<>(allDiagSums[0].length); //list.get(allDiagsumsPisiotn)=list(nodeIds_that_validated_it)
                for (int i = 0; i < allDiagSums[0].length; i++) {
                    detectedPeeks.add(i,new ArrayList<>());
                }
                
                int peekCounter=0;
                for (int node = 0; node < allDiagSums.length; node++) {
                    int nodeId=diagsumIndex[node];
                    //Infos.println("Peek detection in nodeid="+nodeId);
                    //build the inital sum from the w first positions
                    double PPStar0=allDiagSums[node][0];
                    double PPStarSum=PPStar0;
                    for (int j = 1; j < w-1; j++) {
                        PPStarSum+=allDiagSums[node][j];
                    }
                    //compare central value to the mean of the window
                    //if superiior to set ratio, it's a hit
                    double ratio=allDiagSums[node][0+(w/2)]/(PPStarSum/w);
                    ratios[node][0]=ratio;
                    if(ratio<=r_win) {
                        detectedPeeks.get(0).add(nodeId);
                        ratios[node][0]=ratio;
                        //Infos.println("  PEEK DETECTED AT "+0);
                    }
                    //now shift to right 
                    for (int pos=1;pos<align.getLength()-minOverlap-(w/2);pos++) {
                        //substract previous 1st position  
                        PPStarSum-=PPStar0;
                        //add next w-th position
                        PPStarSum+=allDiagSums[node][pos+w-1];
                        ratio=allDiagSums[node][pos+(w/2)]/(PPStarSum/w);
                        ratios[node][pos]=ratio;
                        if(allDiagSums[node][pos+(w/2)]/(PPStarSum/w)<=r_win) {
                            detectedPeeks.get(pos).add(nodeId);
                            ratios[node][pos]=ratio;
                            //Infos.println("  PEEK DETECTED AT "+pos);
                        }
                        PPStar0=allDiagSums[node][pos];
                    }
                }
                
                //check r_nodes, the ratio of nodes that vlaidated a position
                double[] ratios2=new double[allDiagSums[0].length];
                ArrayList<Integer> positionsToAnalyse=new ArrayList<>();
                for (int pos = 0; pos < detectedPeeks.size(); pos++) {
                    double ratio=detectedPeeks.get(pos).size()/allDiagSums.length;
                    ratios2[pos]=ratio;
                    if ( ratio >= r_nodes) {                            
                        Infos.println("Position "+pos+" passes r_nodes test !");
                        positionsToAnalyse.add(pos);
                    }
                }
                
                endScanTime=System.currentTimeMillis();
                Infos.println("Peek search done on "+(endScanTime-startScanTime)+" ms");
                
                
                //will not analyse this query further if it appears to be a repeat
                if ( positionsToAnalyse.size()>maxAllowedPeeks) {
                    Infos.println("This query is probably a repeat, more than "+maxAllowedPeeks+" possible localizations detected !" );
                    Infos.println("QUERY DISCARDED FROM THE ANALYSIS.");
                    queryCounter++;
                    continue;
                } else if (positionsToAnalyse.size()<1) {
                    Infos.println("No peek detected from this query." );
                    Infos.println("QUERY DISCARDED FROM THE ANALYSIS.");
                    queryCounter++;
                    continue;
                }
                Infos.println(positionsToAnalyse.size()+" localizations will be analysed further, for precise placement.");

                
                
                
                ///////////////////////////////////////////////////////////////////////////
                //DEEP ANALYSIS AROUND EACH SELECTED PEEK --> final phylogenetic placement
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                queryCounter++;
            }
            
            
            
            fw.close();
            fp.closePointer();
            
            
            

            
            
            
            
            
            
            
            System.exit(0);
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PREPLACEMENT.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    

    
    
    
    
    
    
    
    
}
