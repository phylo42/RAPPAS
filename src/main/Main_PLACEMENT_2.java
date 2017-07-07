/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import au.com.bytecode.opencsv.CSVWriter;
import core.AAStates;
import core.DNAStates;
import core.PProbasSorted;
import core.QueryWord;
import core.hash.SimpleHash;
import core.States;
import core.algos.SequenceKnife;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class Main_PLACEMENT_2 {
    
    //sequence type
    public static int TYPE_DNA=1;
    public static int TYPE_PROTEIN=2;
    //memory mode
    public static int MEMORY_LOW=1;
    public static int MEMORY_LARGE=2;
    
    public static void main(String[] args) {

        try {
            System.setProperty("debug.verbose","1");
            
            ////////////////////////////////////////////////////////////////////
            // INPUT QUERY
            //pplacer benchmark, stat not based on relaxed tree
            String inputsPath="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/";
            String q=inputsPath+"mod_p4z1r36_query_only2.fasta";

            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //Search mode (memory_low = no hash, PPStar calculated directly from the PPProbas);
            int memory_mode=MEMORY_LARGE;
            
            
            //preplacement parameters//////////////////////////////////////////
            //for a number of nodes =total_nodes/nodeShift , build diagsum vectors
            //think to keep sum of the diagsums to do mean at the end and highlight positions > to mean
            int nodeShift=100;
            // pour stocker ou non dans hash
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            //word sampling method
            int preplacementWordSampling=SequenceKnife.SAMPLING_LINEAR;
            //peek detection parameters
            int w=31; //sliding window size
            int w_min=17;//minimum window size (for diagSums borders, basically, the w_win/2 1st position will not show a ratio, i.e. output value is -1.0)
            if (w_min>w) {w_min=w;}//!!! must be w_min<=w. if not, the program force w_min=w
            double r_win=0.8; //peek/windows_mean ratio
            double r_nodes=0.20; //node_where_peek_detected/scanned_nodes ratio
            int maxAllowedPeeks=50; //maximum number of localizations allowed for a read, if greater read discarded because clearly a repeat
            //fine placement parameters//////////////////////////////////////////
            //word sampling method
            int wordSampling=SequenceKnife.SAMPLING_NON_OVERLAPPING;
            //searchMethod
            boolean deepSearch=false;
            //when doing deepsearch, bp on left and right around the (peek;peek+query_length) interval
            int peekBoundary=10; 
            //the 'top' nodes (of highest score) are reported for a peek
            //int topNodesSize=5;
            
            //debug/////////////////////////////////////////////////////////////
            //max number of queries treated 
            
            int queryLimit=1000000;
            //which log to write, !!!
            //more logs= much slower placement because of disk access latency
            boolean logDetailedDiagsums=true;
            boolean logDiagsums=true;
            boolean logPeekRatios=true;
            boolean logPreplacementDiagsums=true;
            boolean logPreplacementDetailedDiagsums=true;

;
            

            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOADING ALL DATA FROM THE INPUT DATABASE
            long startLoadTime=System.currentTimeMillis();
            
            //LOAD SESSION//////////////////////////////////////////////////////
            //base path for outputs
            String path="/media/ben/STOCK/DATA/viromeplacer/WD/";
            //logs
            String logPath=path+"logs/";
            //trees
            String relaxedTreePath=path+"relaxed_trees/";
            //ancestral reconstruciton
            String ARPath=path+"AR/";
            //session itself
            boolean loadHash=true;
            if (memory_mode==MEMORY_LOW) {loadHash=false;}
            SessionNext session= SessionNext.load(new File(path+"PAML_session_params_k10_mk10_f2.0_t9.765625E-4"),loadHash);
            
            //type of Analysis//////////////////////////////////////////////////
            States s=session.states; 
            int analysisType=-1;
            //States: DNA or AA
            if (s instanceof DNAStates)
                analysisType=Main_PLACEMENT_2.TYPE_DNA;
            else if (s instanceof AAStates)
                analysisType=Main_PLACEMENT_2.TYPE_PROTEIN;
            
            //posterior probas analysis/////////////////////////////////////////
            //mers size
            int k=session.k;
            int min_k=session.minK;
            Infos.println("k="+k);
            Infos.println("min_k="+min_k);
            //site and word posterior probas thresholds
            float sitePPThreshold=session.stateThreshold;
            float wordPPStarThreshold=session.wordThreshold;
            float thresholdFactor=session.factor;
            double thresholdAsLog=Math.log10(wordPPStarThreshold);
            double wordAbsent=Math.pow(0.25,k); // not used in hash, but for scoring queries

            Infos.println("factor="+thresholdFactor);
            Infos.println("sitePPThreshold="+sitePPThreshold+" (for info, unused)");
            Infos.println("wordPPStarThreshold="+wordPPStarThreshold);
            Infos.println("wordPPStarThreshold(log10)="+thresholdAsLog);
            Infos.println("wordAbsent="+wordAbsent);
            



            //PREPARE DIRECTORIES///////////////////////////////////////////////
            if (!new File(path).exists()) {new File(path).mkdir();}
            if (!new File(logPath).exists()) {new File(logPath).mkdir();}
            if (!new File(relaxedTreePath).exists()) {new File(relaxedTreePath).mkdir();}
            if (!new File(ARPath).exists()) {new File(ARPath).mkdir();}
            
            
            //LOAD ORIGINAL ALIGNMENT///////////////////////////////////////////
            Alignment align=session.align;
            Infos.println(align.describeAlignment(false));
            
            //LOAD TREES////////////////////////////////////////////////////////
            NewickReader np=new NewickReader();
            PhyloTree tree = session.tree;
            NewickWriter nw=new NewickWriter(new File("null"));
            String relaxedTreeForJplace=nw.getNewickTree(tree, true, true, true);
            //Infos.println(relaxedTreeForJplace);
            nw.close();
            
            //LOAD THE POSTERIOR PROBAS/////////////////////////////////////////
            PProbasSorted pprobas = session.parsedProbas;
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(pprobas.getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateIndexSet(0, 0, 5)));
            
            SimpleHash hash=null;
            if (memory_mode==MEMORY_LARGE) {
                hash=session.hash;
                Infos.println("Hash test: "+hash.getRegisteredWords().size()+" keys");
            }
            
            Infos.println(Environement.getMemoryUsage());
            long endLoadTime=System.currentTimeMillis();
            Infos.println("Loading the database took "+(endLoadTime-startLoadTime)+" ms");

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOAD FINISHED, TIME TO START THE PLACEMENT PROCESS
            


            ////////////////////////////////////////////////////////////////////
            //MAP TO REGISTER PLACEMENTS: map(query)=(map(nodeId)=positions_in_align)  !position in align, not diagsum
            HashMap<Fasta,HashMap<Integer,ArrayList<Integer>>> placementsPerQuery=new HashMap();
            LinkedHashMap<Fasta,ArrayList<Fasta>> identicalQueries=new LinkedHashMap();

            /////////////////////
            //PLACEMENT OF THE QUERIES
            FASTAPointer fp=new FASTAPointer(new File(q), false);
            FileWriter fw =new FileWriter(new File("Queries.fasta"));
            
            int queryCounter=0;
            int queryMatchingRefCounter=0;
            Fasta fasta=null;
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                if (queryCounter>=queryLimit)
                    break;
                
                Infos.println("#######################################################################");
                Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+fasta.getHeader());
                Infos.println("#######################################################################");
                fw.append(fasta.getFormatedFasta()+"\n");
                
                //check if another query didn't had exactly the same sequence
                //fasta is a Comparable based on the sequence
                if (identicalQueries.containsKey(fasta)) {
                    identicalQueries.get(fasta).add(fasta);
                    Infos.println("Sequence identicial to "+identicalQueries.get(fasta).get(0).getHeader()+". Placement skipped.");
                    queryCounter++;
                    continue;
                }
                identicalQueries.put(fasta, new ArrayList<>());
                identicalQueries.get(fasta).add(fasta);
                
                
                
                
                if (!placementsPerQuery.containsKey(fasta))
                    placementsPerQuery.put(fasta, new HashMap<>());
                
                long startScanTime=System.currentTimeMillis();
                
                int scannedNodesCounter=0;
                Infos.println("Number of internal nodes: "+tree.getInternalNodesByDFS().size());
                int scannedNodes= tree.getInternalNodesByDFS().size()/nodeShift ;
                Infos.println("Number of internal nodes to scan: "+scannedNodes);
                
                //query length
                int queryLength=fasta.getSequence(false).length();
                Infos.println("Query length: "+queryLength);
                //diagSums shift on the left
                int diagsumShift=queryLength-minOverlap;
                Infos.println("Diagsum shift: "+diagsumShift);
                //diagSums table
                int[] preplacementDiagsumIndex = new int[scannedNodes]; //f index[0] nodeId23 ; index[1]=nodeID112 ...
                Arrays.fill(preplacementDiagsumIndex,-1);
                //the diagSums vector has the length of the ref alignment for now
                //this is FOR NOW taking into account partial overlaps on left and right
                //but now when the read is a complementary sequence, these diagsums are ignored... for now.
                float[][] preplacementAllDiagsums = new float[scannedNodes][align.getLength()+diagsumShift-minOverlap];
                Infos.println("Diagsum vector size: "+(align.getLength()+diagsumShift-minOverlap));
                
                
                /////////////////////////////
                // LOG OUTPUT 1:words detais
                CSVWriter writerDetails=null;
                if (logPreplacementDetailedDiagsums) {
                    Infos.println("Write logs: preplacement_diagSums_details.tsv");
                    writerDetails=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_preplacement_diagSums_details.tsv")), '\t');
                    String[] headerData=new String[2+align.getLength()];
                    headerData[0]="node";
                    headerData[1]="mer";
                    for (int i = 0; i < align.getLength(); i++) {
                        headerData[i+2]="pos="+i;
                    }
                    writerDetails.writeNext(headerData);
                }
                /////////////////
                
                
                ////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                //LOW_MEMORY PRE-PLACEMENT, ON LIMITED NODE NUMBER
                // --> WORD SCANNED OVER THE REFERENCE
                ////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                if (memory_mode==MEMORY_LOW) {
                    ArrayList<Integer> preplacementNodeIds=new ArrayList<>();
                    for (int i=0;i<tree.getInternalNodesByDFS().size();i++) {
                        if (i%nodeShift!=0) {continue;}
                        if (scannedNodesCounter>scannedNodes-1) {break;} //all nodes done

                        int nodeId=tree.getInternalNodesByDFS().get(i);
                        preplacementNodeIds.add(nodeId);
                        //Infos.println("Building diagSum vector #"+scannedNodesCounter+" for node: "+nodeId);
                        preplacementDiagsumIndex[scannedNodesCounter]=nodeId; //nodes ids are from 0 to N

                        ///////////////////////////////////
                        // LOOP ON QUERY WORDS
                        QueryWord qw=null;
                        SequenceKnife sk=new SequenceKnife(fasta, k, min_k, s, preplacementWordSampling);
                        //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                        int queryWordCounter=0;
                        while ((qw=sk.getNextWord())!=null) {
                            //Infos.println("Query mer: "+qw.toString());

                            /////////////////
                            // LOG OUTPUT 1
                            String[] data=null;
                            if (logPreplacementDetailedDiagsums) {
                                data=new String[2+align.getLength()];
                                data[0]=String.valueOf(nodeId);
                                data[1]=String.valueOf(qw.getOriginalPosition());
                            }
                            /////////////////


                            //////////////////////////////////////
                            //BUILD DIAGSUMS, ONE PER NODE
                            for (int refPos=0;refPos<align.getLength()-k;refPos++) {  
                                double PPStar=0.0;
                                for (int queryPos=0;queryPos<qw.getWord().length;queryPos++) {
                                    PPStar+=pprobas.getPP(nodeId, refPos+queryPos, pprobas.getStateIndex(nodeId, refPos+queryPos, qw.getWord()[queryPos]) );
//                                    //DEBUG WITH GLKT0ZE01C2HN1
//                                    if ( (refPos>(622+queryWordCounter*k)) && (refPos<(626+queryWordCounter*k)) ) {
//                                        System.out.println("Proba at "+(refPos+queryPos)+" : "+pprobas.getPP(nodeId, refPos+queryPos, index));
//                                        System.out.println("Current Product: PP* = "+PPStar);
//                                    }
                                }
                                //normalization to PPStar
                                //if (PPStar<thresholdAsLog)
                                //    PPStar=thresholdAsLog;

                                //Infos.println("refPos: "+refPos+"  PP* ="+PPStar);
                                //the conditions below are to avoid the last mers of 
                                //the query which would contradict the minOverlap condition:
                                //in the first bases of the alignment (pos<readlength)
                                //do not scan this word if it would represent an overlap 
                                //with the ref align < minOverlap
                                //this could be optimized to avoid some iterations ???
                                int currentDiagSumPos=diagsumShift+refPos-qw.getOriginalPosition();
                                if (currentDiagSumPos>-1 && currentDiagSumPos<preplacementAllDiagsums[0].length) {
                                    preplacementAllDiagsums[scannedNodesCounter][currentDiagSumPos]+=PPStar;
//                                    //DEBUG WITH GLKT0ZE01C2HN1
//                                    if ( (refPos>(622+queryWordCounter*k)) && (refPos<(626+queryWordCounter*k)) ) {
//                                        System.out.println("X= "+(diagsumShift+(refPos-qw.getOriginalPosition())));
//                                        System.out.println("preplacementAllDiagsums[scannedNodesCounter][X]+="+PPStar);
//                                        System.out.println("preplacementAllDiagsums[scannedNodesCounter][X]="+preplacementAllDiagsums[scannedNodesCounter][diagsumShift+(refPos-qw.getOriginalPosition())]);
//                                    }

                                }


                                if (logPreplacementDetailedDiagsums) {
                                    data[2+refPos]=String.valueOf(PPStar); 
                                }
                            }
                            if(logPreplacementDetailedDiagsums)
                                writerDetails.writeNext(data);


                            queryWordCounter++;

                        }

                        //Normalize the first and last positions of the diagsums.
                        //indeed, in the "duagSumShift" first positions
                        //and in the "minoverlap" last positions, not all the words
                        //of the query were summed, (minoverlap condition).
                        //we compensate here by calculating the mean word score in the
                        //overlap and use this mean for all wordw which don't overlap...
                        //the value itself is increased by (wordsInOverlap/totalWordsInTheQuery)
                        //System.out.println(queryWordCounter);
                        for (int diagSumPos=0;diagSumPos<diagsumShift;diagSumPos++) {
                            int wordsInQueryRefOverlap=(minOverlap+diagSumPos-(queryLength%k))/k;
                            //System.out.println("wordsCountInMinOverlap:"+wordsInQueryRefOverlap);
                            preplacementAllDiagsums[scannedNodesCounter][diagSumPos]*=(0.0+queryWordCounter)/wordsInQueryRefOverlap;

                        }
                        for (int diagSumPos=preplacementAllDiagsums[0].length-diagsumShift;diagSumPos<preplacementAllDiagsums[0].length;diagSumPos++) {
                            int wordsInQueryRefOverlap=(preplacementAllDiagsums[0].length-diagSumPos+minOverlap)/k;
                            //System.out.println("diagsumPos:"+diagSumPos+" wordsCountInMinOverlap:"+wordsInQueryRefOverlap);
                            preplacementAllDiagsums[scannedNodesCounter][diagSumPos]*=(0.0+queryWordCounter)/wordsInQueryRefOverlap;
                        }

                        //Infos.println("# Mer scanned: "+scannedMersCounter);

                        scannedNodesCounter++;
                    }


                    if(logPreplacementDetailedDiagsums) {
                        writerDetails.flush();
                        writerDetails.close();
                    }


                ////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                //LARGE_MEMORY PRE-PLACEMENT, HASH-BASED
                // --> HIGHEST SCORING NODES PER WORD ARE PICKED-UP
                // --> words absent from hash have proba: factor.(1/4^k)
                ////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                } else if (memory_mode==MEMORY_LARGE) {
                    
                    HashMap<Integer,Integer> nodeVisited=new HashMap<>();
                    preplacementDiagsumIndex = new int[session.tree.getNodeCount()]; //f index[0] nodeId23 ; index[1]=nodeID112 ...
                    Arrays.fill(preplacementDiagsumIndex,-1);
                    preplacementAllDiagsums = new float[session.tree.getNodeCount()][align.getLength()+diagsumShift-minOverlap];
                    Infos.println("Diagsum vector size: "+(align.getLength()+diagsumShift-minOverlap));                    
                    ///////////////////////////////////
                    // LOOP ON QUERY WORDS
                    QueryWord qw=null;
                    SequenceKnife sk=new SequenceKnife(fasta, k, min_k, s, preplacementWordSampling);
                    //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                    int queryWordCounter=0;
                    int queryWordFoundCounter=0;
                    int nodeVisitedCounter=0;
                    
                    float bottomVal=(float)(sk.getMerOrder().length*Math.log10(wordAbsent));
                    Infos.println("bottomVal="+bottomVal);
                    while ((qw=sk.getNextWord())!=null) {
                        //Infos.println("Query mer: "+qw.toString());
                        
                        SimpleHash.Tuple topTuple = hash.getTopTuple(qw);
                        //System.out.println("Match: "+topTuple);
                        //if this word is not registered in the hash
                        if (topTuple==null) {
                            
                            
                            
                            continue;
                        }
                        queryWordFoundCounter++;

                        
                        if (!nodeVisited.containsKey(k)) {
                            //System.out.println("NEW NODE:"+topTuple.getNodeId());
                            nodeVisited.put(topTuple.getNodeId(), nodeVisitedCounter);
                            
                            preplacementDiagsumIndex[nodeVisitedCounter]=topTuple.getNodeId(); //nodes ids are from 0 to N
                            Arrays.fill(preplacementAllDiagsums[nodeVisitedCounter],bottomVal);
                            //System.out.println("Set table:"+Arrays.toString(Arrays.copyOfRange(preplacementAllDiagsums[nodeVisitedCounter],0,5)));
                            nodeVisitedCounter++;
                        }
                        int diagSumPos=topTuple.getRefPos()-qw.getOriginalPosition()+diagsumShift;
                        //System.out.println("Match diagsumPos:"+diagSumPos);
                        if (diagSumPos>-1 && diagSumPos<preplacementAllDiagsums[0].length)
                            preplacementAllDiagsums[nodeVisitedCounter][diagSumPos]+=(-wordAbsent+topTuple.getPPStar());
                        
                        queryWordCounter++;

                    }
                    System.out.println("proportion of words found: "+queryWordFoundCounter+"/"+queryWordCounter);
                    System.out.println("#nodes hit: "+nodeVisitedCounter);

      
                    
                }
                
                
                
                
                

                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                long endScanTime=System.currentTimeMillis();
                Infos.println("Scan on "+(scannedNodesCounter)+" nodes took "+(endScanTime-startScanTime)+" ms");
                
                /////////////////
                // LOG OUTPUT 2
                Infos.println("Write logs: preplacement_diagSums.tsv");
                CSVWriter writerDiagsum=null;
                if (logPreplacementDiagsums) {
                    writerDiagsum=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_preplacement_diagSums.tsv")), '\t');
                    String[] headerData=new String[2+preplacementAllDiagsums.length];
                    headerData[0]="diagsum_position";
                    headerData[1]="reference_align_site";
                    for (int i = 0; i < preplacementAllDiagsums.length; i++) {
                        headerData[i+2]="nodeid="+preplacementDiagsumIndex[i];
                    }
                    writerDiagsum.writeNext(headerData);

                    for (int i = 0; i < preplacementAllDiagsums[0].length; i++) { //for each site
                        String[] data=new String[2+scannedNodes];
                        data[0]=String.valueOf(i);
                        data[1]=String.valueOf(i-diagsumShift);
                        for (int j = 0; j < scannedNodes; j++) {
                            data[j+2]=String.valueOf(preplacementAllDiagsums[j][i]);
                        }
                        writerDiagsum.writeNext(data);
                    }
                    writerDiagsum.flush();
                    writerDiagsum.close();
                }
                /////////////////
                
                                
                
                
                
                
                
                ////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                //DETECTION OF PEEKS VIA RATIOS
                ////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                
                startScanTime=System.currentTimeMillis();
                //the ratios table will be of length = preplacementDiagsum.length
                //pour les w première/dernières positions, la moyenne est caculée
                //sur un nombre de positions <w et égale à i+/-(w/2-1)
                double[][] ratios=new double[preplacementAllDiagsums.length][preplacementAllDiagsums[0].length];
                //list.get(allDiagsumsPosition)=list(nodeIds_that_validated_it)
                ArrayList<ArrayList<Integer>> preplacementPeeks=new ArrayList<>(preplacementAllDiagsums[0].length); 
                for (int i = 0; i < preplacementAllDiagsums[0].length; i++) {
                    preplacementPeeks.add(i,new ArrayList<>());
                }
                
                int peekCounter=0;
                for (int node = 0; node < preplacementAllDiagsums.length; node++) {
                    int nodeId=preplacementDiagsumIndex[node];
                    //Infos.println("Peek detection in nodeid="+nodeId);
                    
                    //build the inital ratios for the w-1 first positions of diagsums
                    //these are all w_min<windows and 1st windows of size w
                    //in this case, w_min is taken into account because the full
                    double ratio=-1.0;
                    double PPStarSum=0.0;
                    for (int diagsumPos = 0; diagsumPos < w; diagsumPos++) {
                        PPStarSum+=preplacementAllDiagsums[node][diagsumPos];
                        //System.out.print("diagsumPos="+diagsumPos);
                        //we build ratios only for non pair windows
                        if (diagsumPos%2==0) {
                            if ((diagsumPos+1)>=w_min) {
                                //System.out.print(" *PPStarSum="+PPStarSum+" ratio["+(diagsumPos)/2+"]="+preplacementAllDiagsums[node][diagsumPos/2]+"/("+PPStarSum+"/"+(diagsumPos+1)+")");
                                ratio=preplacementAllDiagsums[node][diagsumPos/2]/(PPStarSum/(diagsumPos+1));
                                ratios[node][diagsumPos/2]=ratio;
                                //compare central value to the mean of the window
                                //if superiior to set ratio, it's a hit
                                if(ratio<=r_win) {
                                    preplacementPeeks.get(diagsumPos/2).add(nodeId);
                                    //System.out.println("  PEEK DETECTED AT "+(diagsumPos/2));
                                }
                            }
                        }
                        //System.out.println();
                        if (diagsumPos<(w_min/2)) {
                            //System.out.println(" *PPStarSum="+PPStarSum+" -->skip diagsumPos="+diagsumPos+" because of w_min. ");
                            ratios[node][diagsumPos]=-1.0;
                        }
                    }                    
                    //now shift to right, windows are all of size w
                    //2nd to last windows of size w
                    for (int diagsumPos=w;diagsumPos<preplacementAllDiagsums[0].length;diagsumPos++) {
                        //System.out.print("diagsumPos="+diagsumPos+" ");
                        //substract previous 1st position  
                        PPStarSum-=preplacementAllDiagsums[node][diagsumPos-w];
                        //and add next position
                        PPStarSum+=preplacementAllDiagsums[node][diagsumPos];
                        //System.out.print(" PPStarSum="+PPStarSum+" ratio["+(diagsumPos-w/2)+"]="+preplacementAllDiagsums[node][diagsumPos-w/2]+"/("+PPStarSum+"/"+w+")");
                        ratio=preplacementAllDiagsums[node][diagsumPos-w/2]/(PPStarSum/w);
                        ratios[node][diagsumPos-w/2]=ratio;
                        if(ratio<=r_win) {
                            preplacementPeeks.get(diagsumPos-w/2).add(nodeId);
                            //System.out.println("  PEEK DETECTED AT "+(diagsumPos-w/2));
                        }
                        //System.out.println();
                    }
                    //build theratios ratios for the last w-1 positions
                    for (int diagsumPos = preplacementAllDiagsums[0].length-(w-1); diagsumPos < preplacementAllDiagsums[0].length; diagsumPos++) {
                        PPStarSum-=preplacementAllDiagsums[node][diagsumPos-1];
                        //System.out.print("diagsumPos="+diagsumPos);
                        int w_size=preplacementAllDiagsums[0].length-diagsumPos;
                        //we build ratios only for non pair windows
                        if (diagsumPos%2==0) {
                            if (w_size>=w_min) {
                                //System.out.print(" *PPStarSum="+PPStarSum+" ratio["+(diagsumPos+w_size/2)+"]="+preplacementAllDiagsums[node][diagsumPos+w_size/2]+"/("+PPStarSum+"/"+(w_size)+")");
                                ratio=preplacementAllDiagsums[node][diagsumPos+w_size/2]/(PPStarSum/w_size);
                                ratios[node][diagsumPos+w_size/2]=ratio;
                                //compare central value to the mean of the window
                                //if superiior to set ratio, it's a hit
                                if(ratio<=r_win) {
                                    preplacementPeeks.get(diagsumPos+w_size/2).add(nodeId);
                                    //System.out.println("  PEEK DETECTED AT "+(diagsumPos+w_size/2));
                                }
                            } 
                        }
                        if ((diagsumPos+w_min/2)>preplacementAllDiagsums[0].length-1){
                            //System.out.print(" *PPStarSum="+PPStarSum+" -->skip diagsumPos="+diagsumPos+" because of w_min. ");
                            ratios[node][diagsumPos]=-1.0;
                        }
                        //System.out.println();
                    } 
                }
                
                //check r_nodes, the ratio of nodes that validated a position
                double[] ratios2=new double[preplacementAllDiagsums[0].length];
                ArrayList<Integer> diagsumPosToAnalyseFurther=new ArrayList<>();
                int positionWithPeeks=0;
                for (int diagsumPos = 0; diagsumPos < preplacementPeeks.size(); diagsumPos++) {
                    //System.out.println(diagsumPos+" "+preplacementPeeks.get(diagsumPos));
                    double ratio=(0.0+preplacementPeeks.get(diagsumPos).size())/preplacementAllDiagsums.length;
                    ratios2[diagsumPos]=ratio;
                    if (preplacementPeeks.get(diagsumPos).size()>0)
                        positionWithPeeks++;
                    if ( ratio >= r_nodes) {                            
                        Infos.println("Position "+(diagsumPos-diagsumShift)+" (diagsum["+diagsumPos+"]) passes test r_nodes>="+ratio+" !");
                        diagsumPosToAnalyseFurther.add(diagsumPos);
                    }
                }
                Infos.println(diagsumPosToAnalyseFurther.size()+"/"+positionWithPeeks+" positions pass r_nodes test.");
                endScanTime=System.currentTimeMillis();
                Infos.println("Peek search took "+(endScanTime-startScanTime)+" ms");
                                
                
                ///////////////////////////////////////////////////////////////////////////
                // LOG RATIOS
                if (logPeekRatios) {
                    Infos.println("Write logs: ratio_node.tsv");
                    CSVWriter writerRatios=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_ratio_window.tsv")), '\t');
                    CSVWriter writerRatios2=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_ratio_node.tsv")), '\t');
                    
                    //ratios2[pos]=ratio;
                    String[] data=new String[3];
                    data[0]="diasum_pos";
                    data[1]="reference_align_site";
                    data[2]="ratio_node";
                    writerRatios2.writeNext(data);
                    for (int i = 0; i < ratios2.length; i++) {
                        data=new String[3];
                        data[0]=String.valueOf(i);
                        data[1]=String.valueOf(i-diagsumShift);
                        data[2]=String.valueOf(ratios2[i]);
                        writerRatios2.writeNext(data);
                    }
                    writerRatios2.close();
                    
                    Infos.println("Write logs: ratio_windows.tsv");
                    //double[][] ratios=new double[preTreatmentAllDiagsums.length][preTreatmentAllDiagsums[0].length];
                    String[] headerData=new String[2+scannedNodes];
                    headerData[0]="diagsum_position";
                    headerData[1]="reference_align_site";
                    for (int i = 0; i < preplacementAllDiagsums.length; i++) {
                        headerData[i+2]="nodeid="+preplacementDiagsumIndex[i];
                    }
                    writerRatios.writeNext(headerData);
                    for (int i = 0; i < ratios[0].length; i++) { //for each site
                        data=new String[2+scannedNodes];
                        data[0]=String.valueOf(i);
                        data[1]=String.valueOf(i-diagsumShift);
                        for (int j = 0; j < scannedNodes; j++) {
                            data[j+2]=String.valueOf(ratios[j][i]);
                        }
                        writerRatios.writeNext(data);
                    }
                    writerRatios.close();
                }
                ///////////////////////////////////////////////////////////////////////////
                
                
                ///////////////////////////////////////////////////////////////////////////
                //DISCARD POTENTIAL REPEATS

                //will not analyse this query further if it appears to be a repeat
                if ( diagsumPosToAnalyseFurther.size()>maxAllowedPeeks) {
                    Infos.println("This query is probably a repeat, more than "+maxAllowedPeeks+" possible localizations detected !" );
                    Infos.println("QUERY DISCARDED FROM THE ANALYSIS.");
                    queryCounter++;
                    continue;
                } else if (diagsumPosToAnalyseFurther.size()<1) {
                    Infos.println("No peek detected from this query." );
                    Infos.println("QUERY DISCARDED FROM THE ANALYSIS.");
                    queryCounter++;
                    continue;
                }
                Infos.println(diagsumPosToAnalyseFurther.size()+" localization(s) will undergo precise placement.");
                queryMatchingRefCounter++;
                
                
                
                
                ///////////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////
                //REMAING NODES ANALYSIS, AROUND EACH SELECTED PEEK 
                ///////////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////

                ////////////////////////////////////////////////////////////////
                //we can avoid the nodes for which diagSums where allReady built
                //copy the corresponding regions of pretreatmentAllDiagsums
                //
                //IMPORTANT NOTE, HERE THE SELECTION OF NODES COULD BE DONE
                //AROUND THE BEST NODES (best-scoring diagsums), 
                ////////////////////////////////////////////////////////////////
                //Infos.println(scannedNodesCounter+" nodes were analysed for the preplacement.");
                int[] nodesToAnalyse=new int[tree.getInternalNodesByDFS().size()-scannedNodesCounter];
                int counter=0;
                for (int nodeId: tree.getInternalNodesByDFS()) {
                    boolean alreadyDone=false;
                    for (int i = 0; i < preplacementDiagsumIndex.length; i++) {
                        if (preplacementDiagsumIndex[i]==nodeId)
                            alreadyDone=true;
                    }
                    if (!alreadyDone) {
                        nodesToAnalyse[counter]=nodeId;
                        counter++;
                    }
                }
                Infos.println("Remaining #nodes for localized analysis:"+nodesToAnalyse.length);
                //copy preplacement diagsums to final diagsums
                int[] diagsumIndex=new int[preplacementDiagsumIndex.length+nodesToAnalyse.length];
                float[][] selectedDiagsums=new float[preplacementAllDiagsums.length+nodesToAnalyse.length][preplacementAllDiagsums[0].length];
                for (int node = 0; node < preplacementDiagsumIndex.length; node++) {
                    diagsumIndex[node]=preplacementDiagsumIndex[node];
                    selectedDiagsums[node]=preplacementAllDiagsums[node];
                }
                for (int column = 0; column < nodesToAnalyse.length; column++) {
                    diagsumIndex[preplacementDiagsumIndex.length+column]=nodesToAnalyse[column];
                }
                //System.out.println(Arrays.toString(diagsumIndex));
                
                ////////////////////////////
                // LOG OUTPUT 1:words details
                writerDetails=null;
                
                if (logDetailedDiagsums) {
                    Infos.println("Write logs: diagSums_details.tsv");
                    writerDetails=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_diagSums_details.tsv")), '\t');
                    String[] headerData=new String[2+align.getLength()];
                    headerData[0]="node";
                    headerData[1]="mer";
                    for (int i = 0; i < align.getLength(); i++) {
                        headerData[i+2]="pos="+i;
                    }
                    writerDetails.writeNext(headerData);
                }
                /////////////////
                
                
                /////////////////
                // LOG OUTPUT 2: diagsums
                writerDiagsum=null;
                if (logDiagsums) {
                    Infos.println("Write logs: diagSums.tsv");
                    writerDiagsum=new CSVWriter(new FileWriter(new File(logPath+fasta.getHeader()+"_diagSums.tsv")), '\t');
                    String[] headerData=new String[2+diagsumIndex.length];
                    headerData[0]="diagsum_position";
                    headerData[1]="reference_align_site";
                    for (int i = 0; i < diagsumIndex.length; i++) {
                        headerData[i+2]="nodeid="+diagsumIndex[i];
                    }
                    writerDiagsum.writeNext(headerData);
                }
                /////////////////


                
                ////////////////////////////////////////////////////////////////
                // LOOP ON PEEKS, AND SELECT HIGHEST SCORE FOR EACH PEEK
                for (int diagsumPosToAnalyse:diagsumPosToAnalyseFurther) {
                    Infos.println("-->Analysing remaining nodes for peek="+(diagsumPosToAnalyse-diagsumShift)+" (diagsum["+diagsumPosToAnalyse+"])");
                    


                    startScanTime=System.currentTimeMillis();
                    
                    ///////////////////////////////////
                    // LOOP ON REMAINING NODES
                    int scannedNodesCounter2=preplacementAllDiagsums.length;
                    for (int i = 0; i < nodesToAnalyse.length; i++) {
                        int nodeId=nodesToAnalyse[i];
                        diagsumIndex[scannedNodesCounter2]=nodeId;        
                        //Infos.println(" NodeId="+nodeId+" ("+scannedNodesCounter2+")");
                        
                        ///////////////////////////////////
                        // LOOP ON QUERY WORDS
                        QueryWord qw=null;
                        SequenceKnife rk=new SequenceKnife(fasta, k, min_k, s, wordSampling);
                        while ((qw=rk.getNextWord())!=null) {
                            //Infos.println("Query mer: "+qw.toString());

                            /////////////////
                            // LOG OUTPUT 1
                            String[] data= null;
                            if (logDetailedDiagsums) {
                                data=new String[2+selectedDiagsums[0].length];
                                data[0]=String.valueOf(nodeId);
                                data[1]=String.valueOf(qw.getOriginalPosition());
                            }
                            /////////////////

                            if (deepSearch) {
                                ///////////////////////////////////////////////////////////////
                                //BUILD DIAGSUMS, ONE PER NODE, AROUND THE READ + PEEKBOUNDARY
                                //NOTE: NECESSARY FOR READS WITH LOTS OF INDELS ? (RNAs ?)
                                for (int diagsumPos=diagsumPosToAnalyse-peekBoundary;diagsumPos<diagsumPosToAnalyse+queryLength+peekBoundary+1;diagsumPos++) {  
                                    //System.out.println("Pos: "+diagsumPos);
                                    //to avoid boundary issues
                                    if (diagsumPos<0 || diagsumPos>selectedDiagsums[0].length) {
                                        //System.out.println("Avoid diagsum overflow (due to peekBoundary):"+diagsumPos);
                                        continue;
                                    }

                                    double PPStar=0.0;
                                    for (int queryPos=0;queryPos<qw.getWord().length;queryPos++) {
                                        PPStar+=pprobas.getPP(nodeId, diagsumPos+queryPos-diagsumShift, pprobas.getStateIndex(nodeId, diagsumPos+queryPos-diagsumShift, qw.getWord()[queryPos]) );
                                    }

                                    //the condition below is to avoid the last mers of 
                                    //the query which would contradict the minOverlap condition:
                                    if (diagsumShift+(diagsumPos-qw.getOriginalPosition())>-1) 
                                        selectedDiagsums[scannedNodesCounter2][diagsumPos-qw.getOriginalPosition()]+=PPStar; 
                                    if (logDetailedDiagsums) {
                                        data[2+diagsumPos]=String.valueOf(PPStar); 
                                    }
                                }
                            } else {
                                ///////////////////////////////////////////////////////////////
                                //READS WITHOUT INDELS, SHOULD REQUIRE A  SINGLE
                                //DIAGSUM POSITION CALCULATION
                                                            
                                double PPStar=0.0;
                                for (int queryPos=0;queryPos<qw.getWord().length;queryPos++) {
                                    try {
                                        int alignPos=diagsumPosToAnalyse+qw.getOriginalPosition()+queryPos-diagsumShift;
                                        if (alignPos<(selectedDiagsums[0].length-diagsumShift) && alignPos>-1) {
                                            PPStar+=pprobas.getPP(nodeId,alignPos, pprobas.getStateIndex(nodeId, alignPos, qw.getWord()[queryPos]) );
                                        }
                                    } catch (Exception ex) {
                                        System.out.println("nodeId="+nodeId+"   diagsumPosToAnalyse="+diagsumPosToAnalyse+"+qw.getOriginalPosition()="+qw.getOriginalPosition()+"+queryPos="+queryPos+"-diagSumShift="+diagsumShift+"  == "+(diagsumPosToAnalyse+qw.getOriginalPosition()+queryPos-diagsumShift));
                                        ex.printStackTrace();
                                        System.exit(1);
                                    }
                                }
                                //if (diagsumShift+(diagsumPosToAnalyse-qw.getOriginalPosition())>-1) 
                                    selectedDiagsums[scannedNodesCounter2][diagsumPosToAnalyse]+=PPStar;
                                if (logDetailedDiagsums) {
                                    data[2+diagsumPosToAnalyse]=String.valueOf(PPStar); 
                                }
                                
                                
                            }
                            
                            
                            
                            
                            if(logDetailedDiagsums)
                                writerDetails.writeNext(data);


                        }
                        //System.out.println("Final diagsum in ["+diagsumPosToAnalyse+"]="+selectedDiagsums[scannedNodesCounter][diagsumPosToAnalyse]);
                        
                        


                        scannedNodesCounter2++;

                    }
                    
                    if (logDetailedDiagsums)
                        writerDetails.close();
                    
                    if (logDiagsums) {
                        for (int i = 0; i < selectedDiagsums[0].length; i++) { //for each site
                            String[] data=new String[2+scannedNodesCounter2];
                            data[0]=String.valueOf(i);
                            data[1]=String.valueOf(i-diagsumShift);
                            for (int j = 0; j < scannedNodesCounter2; j++) {
                                data[j+2]=String.valueOf(selectedDiagsums[j][i]);
                            }
                            writerDiagsum.writeNext(data);
                        }
                    }
 
                    
                    ////////////////////////////////////////////////////////////////
                    // REPORT BEST NODE PLACEMENT FOR CURRENT PEEK
                    int bestNode=-1;
                    double associatedPPStar=Double.NEGATIVE_INFINITY;
                    for (int i = 0; i < selectedDiagsums.length; i++) {
                        double PPStar=selectedDiagsums[i][diagsumPosToAnalyse];
                        if (PPStar>associatedPPStar) {
                            //System.out.println("PPStar>associatedPPStar:"+PPStar+">"+associatedPPStar+"    nodeOdAssociated="+diagsumIndex[i]);
                            associatedPPStar=PPStar;
                            bestNode=diagsumIndex[i];
                        }
                    }
                    if (!placementsPerQuery.get(fasta).containsKey(bestNode))
                        placementsPerQuery.get(fasta).put(bestNode, new ArrayList<Integer>());
                    placementsPerQuery.get(fasta).get(bestNode).add(diagsumPosToAnalyse-diagsumShift);
                    
                    Infos.println("Best placement is: "+tree.getById(bestNode));

                    
                    endScanTime=System.currentTimeMillis();
                    Infos.println("-->Localized analysis around peek="+diagsumPosToAnalyse+" took "+(endScanTime-startScanTime)+" ms");
                    ////////////////////////////////////////////////////////////
                }
                

                
                if (logDiagsums) {
                    writerDiagsum.flush();
                    writerDiagsum.close();
                }
                
                
                
                
                
                queryCounter++;
            }
            
            fw.close();
            fp.closePointer();
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OUTPUT THE PLACEMENTS IN TSV FORMAT
            ////////////////////////////////////////////////////////////////////
            
            CSVWriter fwPlacement=new CSVWriter(new FileWriter(new File(logPath+"placements.tsv")), '\t');
            String[] header=new String[4];
            header[0]="Query";
            header[1]="NodeId";
            header[2]="ExtNodeId";
            header[3]="Loc_start";
            fwPlacement.writeNext(header);
            String[] data=new String[4];
            for (Iterator<Fasta> iterator = placementsPerQuery.keySet().iterator(); iterator.hasNext();) {
                Fasta fastaMaster = iterator.next();
                for (Iterator<Fasta> iteratorSlave = identicalQueries.get(fastaMaster).iterator(); iteratorSlave.hasNext();) {
                    Fasta fastaSlave = iteratorSlave.next();
                    data[0]=fastaSlave.getHeader();
                    for (Iterator<Integer> iterator1 = placementsPerQuery.get(fastaSlave).keySet().iterator(); iterator1.hasNext();) {
                        Integer nextNode = iterator1.next();
                        data[1]=String.valueOf(nextNode);
                        data[2]=String.valueOf("none");//tree.getById(nextNode).getExternalId());
                        for (Iterator<Integer> iterator2 = placementsPerQuery.get(fastaSlave).get(nextNode).iterator(); iterator2.hasNext();) {
                            Integer nextPosition = iterator2.next();
                            data[3]=String.valueOf(nextPosition);
                            fwPlacement.writeNext(data);
                        }
                    }
                }
            }
            fwPlacement.close();
            
            
            
            
            ////////////////////////////////////////////////////////////////////
            //OUTPUT THE PLACEMENTS IN JPLACE (JSON) FORMAT
            ////////////////////////////////////////////////////////////////////
            //with library json.simple
            int level=0;
            JSONObject top=new JSONObject(); 
            LinkedHashMap topMap=new LinkedHashMap();
            //object tree (mandatory)
            topMap.put("tree",relaxedTreeForJplace);
            //object placements (mandatory)
            topMap.put("placements", new ArrayList());
            for (Iterator<Fasta> iterator = identicalQueries.keySet().iterator(); iterator.hasNext();) {
                Fasta currentFasta = iterator.next();
                JSONObject p=new JSONObject();
                //will add section "p" and section "nm" in the object p
                
                
                for (Iterator<Integer> iterator1 = placementsPerQuery.get(currentFasta).keySet().iterator(); iterator1.hasNext();) {
                    Integer next1 = iterator1.next();

                    for (Iterator<Integer> iterator2 = placementsPerQuery.get(currentFasta).get(next1).iterator(); iterator2.hasNext();) {
                        Integer next2 = iterator2.next();

                        
                    }
                }
                //topMap.put("p", );
                //topMap.put("nm", );
                
            }
            top.put("placements","XXX");
            
            //object version (mandatory
            top.put("version",3);
            //object metadata (mandatory): sub-objet "version is mandatory"
            JSONObject invoc=new JSONObject();
            invoc.put("invocation", "viromeplacer xxxxxxxxxxxxxxx");
            top.put("metadata", invoc);
            //object fields
            JSONArray fList=new JSONArray();
            fList.add("edge_id");
            fList.add("son_node_id");
            fList.add("son_node_externa_id");
            fList.add("reference_localization");
            top.put("fields", fList);
            
            //put all the elements in the top JSON object
            top.putAll(topMap);
            String out=top.toJSONString();
            System.out.println(out);
            
            
            

            System.out.println("Result: "+queryMatchingRefCounter+"/"+queryCounter+" needed to be localized.");
            
            
            
            
            
            
            System.exit(0);
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PLACEMENT_2.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        
        
    }
    

    
    
    public static void inputStreamToOutputStream(final InputStream inputStream, final OutputStream out) {
        Thread t = new Thread(new Runnable() {

            @Override
            public void run() {
                try {
                    int d;
                    while ((d = inputStream.read()) != -1) {
                        out.write(d);
                    }
                } catch (IOException ex) {
                    Infos.println(ex.getCause());
                }
            }
        });
        t.setDaemon(true);
        t.start();
    }
    
    
    
    
    
}
