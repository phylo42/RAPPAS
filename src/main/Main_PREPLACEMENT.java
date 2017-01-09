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
            
            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/";
//            String a=wd+"bv_refs_aln.fasta";
            String a=wd+"bv_refs_aln.fasta";
            String t=wd+"RAxML_result.bv_refs_aln";
            String q=wd+"mod_p4z1r36_query_only.fasta";
            String pp=wd+"rst";
            
//            String wd="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/alpha_RNApol/model_GTRnuc/";
//            String a=wd+"mafft_centroids.derep_prefix.Coronovirinae_alpha_RNApol_all_VIPR_20-07-2016_CdsFastaResults_CORRECTED.fasta";
//            String t=wd+"RAxML_bipartitionsBranchLabels.result_alpha_RNApol_REROOTED.tree";
//            String q=wd+"alphaTest1";
//            String pp=wd+"rst";
            
            
            
            
            ///////////////////////////////////////////////////////
            //PARAMETERS
            
            //build of relaxed tree
            double minBranchLength=0.001;
            int numberOfFakeBranchesPerEdge=1;
            //mers size
            int k=4;
            int min_k=4;
            //site and word posterior probas thresholds
            double sitePPThreshold=1e-6;
            double wordPPThreshold=1e-6;
            //for a number of nodes =total_nodes/nodeShift , build diagsum vectors
            //think to keep sum of the diagsums to do mean at the end and highlight positions > to mean
            int nodeShift=50;
            
            //States: DNA or AA
            States s=new DNAStates();    
            
            //debug
            int queryLimit=2;
            
            ///////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////
            
            
            //////////////////////
            //LOAD ORIGINAL ALIGNEMENT
            FASTAPointer fp=new FASTAPointer(new File(a), false);
            Fasta f=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((f=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(f);
            }
            Alignment align=new Alignment(fastas);
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
                if (queryCounter>queryLimit)
                    break;
                Infos.println("### PLACEMENT FOR QUERY: "+f.getHeader());
                fw.append(f.getFormatedFasta()+"\n");
                long startScanTime=System.currentTimeMillis();
                
                int scannedNodesCounter=0;
                Infos.println("Number of internal nodes: "+tree.getInternalNodesByDFS().size());
                int scannedNodes= tree.getInternalNodesByDFS().size()/nodeShift ;
                Infos.println("Number of internal nodes to scan: "+scannedNodes);
                
                //diagSums table
                int[] diagSumIndex = new int[scannedNodes]; //for instance if nodeShift=20 [20,40,60...]
                Arrays.fill(diagSumIndex,-1);
                //the diagSums vector has the length of the ref alignment for now
                //this is NOT taking into account partial overlaps on left and right, these diagsums are ignored...
                double[][] allDiagSums = new double[scannedNodes][align.getLength()];
                
                
                /////////////////
                // LOG OUTPUT 1
                CSVWriter writerDetails=new CSVWriter(new FileWriter(new File("diagSums_details.tsv")), '\t');
                String[] headerData=new String[2+align.getLength()];
                headerData[0]="node";
                headerData[1]="mer";
                for (int i = 0; i < align.getLength(); i++) {
                    headerData[i+2]="pos="+i;
                }
                writerDetails.writeNext(headerData);
                /////////////////
                
                
                /////////////////////////////////////
                //FOR EACH TESTED NODE
                for (int i=0;i<tree.getInternalNodesByDFS().size();i++) {
                    if (i%nodeShift!=0) {continue;}
                    if (scannedNodesCounter>scannedNodes-1) {break;} //all nodes done

                    int nodeId=tree.getInternalNodesByDFS().get(i);
                    PhyloNode node = tree.getById(nodeId);
                    Infos.println("Building diagSum vector #"+scannedNodesCounter+" for node: "+node);
                    diagSumIndex[scannedNodesCounter]=nodeId; //nodes ids are from 0 to N
                    
                    ///////////////////////////////////
                    // GENERATE ONE DIAGSUM PER NODE
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
                        
                        //loop to build PP* on PPstats table from 0 to align length-k.
                        for (int refPos=f.getSequence().length();refPos<(align.getLength()-f.getSequence().length());refPos++) {  //init=qw.getOriginalPosition() pr que rpos-qw.getOriginalPosition() > -1
                            double PPStar=1.0;
                            for (int queryPos=0;queryPos<qw.getWord().length;queryPos++) {
                                PPStar*=pprobas.getPP(nodeId, refPos+queryPos, qw.getWord()[queryPos]);
                            }
                            allDiagSums[scannedNodesCounter][refPos-qw.getOriginalPosition()]+=Math.log10(PPStar);
                            data[2+refPos]=String.valueOf(Math.log10(PPStar)); 
                        }
                        writerDetails.writeNext(data);
                        scannedMersCounter++;
                    }
                    Infos.println("# Mer scanned: "+scannedMersCounter);
                    scannedNodesCounter++;
                }
                
                writerDetails.flush();
                writerDetails.close();
                
                
                long endScanTime=System.currentTimeMillis();
                Infos.println("Scan on "+(scannedNodesCounter+1)+" nodes took "+(endScanTime-startScanTime)+" ms");
                
                /////////////////
                // LOG OUTPUT 2
                CSVWriter writerDiagsum=new CSVWriter(new FileWriter(new File("diagSums.tsv")), '\t');
                headerData=new String[1+scannedNodes];
                headerData[0]="site";
                for (int i = 0; i < allDiagSums.length; i++) {
                    headerData[i+1]="node:id="+diagSumIndex[i];
                }
                writerDiagsum.writeNext(headerData);
                
                for (int i = 0; i < allDiagSums[0].length; i++) { //for each site
                    String[] data=new String[1+scannedNodes];
                    data[0]=String.valueOf(i);
                    for (int j = 0; j < scannedNodes; j++) {
                        data[j+1]=String.valueOf(allDiagSums[j][i]);
                    }
                    writerDiagsum.writeNext(data);
                }
                writerDiagsum.flush();
                writerDiagsum.close();
                /////////////////
                
                
                ///////////////////////////////////////////////
                //TO DO HERE: DETECTION OF PEEKS VIA RATIOS
                int w=41;
                
                
                
                
                ///////////////////////////////////////////////////////////////////////////
                //DEEP ANALYSIS AROUND EACH SELECTED PEEK --> final phylogenetic placement
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                queryCounter++;
            }
            
            
            
            fw.close();
            fp.closePointer();
            
            
            
                        /////////////////
                        // CODE TO SCAN OVER. Finally, not very adapted here but may be used later.
//                        //loop to build PP* on PPstats table from 0 to align length-k.
//                        ArrayList<Double> allPPStar=new ArrayList<>();
//                        //intial product
//                        double old_pos0=pprobas.getPP(nodeId, 0, qw.getWord()[0]);
//                        double PPStar=old_pos0;
//                        for (int i = 1; i < k; i++) {
//                            PPStar*=pprobas.getPP(nodeId, i, qw.getWord()[i]);
//                        }
//                        allPPStar.add(PPStar);
//                        //now shift to right 
//                        for (int pos=1;pos<align.getLength()-k;pos++) {
//                            //divide by previous 1st position  
//                            PPStar/=old_pos0;
//                            //multiply by next k-th position
//                            PPStar*=pprobas.getPP(nodeId, pos+k-1, qw.getWord()[k-1]);   
//                            old_pos0=pprobas.getPP(nodeId, pos, qw.getWord()[0]);
//                            allPPStar.add(PPStar);
//                        }
            
            
            
            
            
            
            
            System.exit(0);
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PREPLACEMENT.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    

    
    
    
    
    
    
    
    
}
