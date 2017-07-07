/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import au.com.bytecode.opencsv.CSVWriter;
import charts.ChartsForNodes;
import core.AAStates;
import core.DNAStates;
import core.DiagSum;
import core.PProbasSorted;
import core.QueryWord;
import core.Score;
import core.hash.SimpleHash;
import core.States;
import core.algos.SequenceKnife;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.awt.GridLayout;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.ui.RefineryUtilities;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloTree;

/**
 * Algo in this version:
 * 
 * 1. all words of the query are searched in the hash
 *    --> we store the corresponding nodes in a list
 *    --> positions are used to define the alignment
 * 2. using the alignment, for all nodes stored in the previous list, 
 *    we search query words associated to the positions 
 *    --> we build diagsums for this node selection nodes (1 per node)
 *    --> we keep as peeks only diagsum position holding more than X words
 *    --> the sum of these peeks define the final score and the placement
 * 
 * 
 * 
 * @author ben
 */
public class Main_PLACEMENT_V04_align_scoreallnodes_diagsumeltsremoved {
    
    //sequence type
    public static int TYPE_DNA=1;
    public static int TYPE_PROTEIN=2;
    //memory mode
    public static int MEMORY_LOW=1;
    public static int MEMORY_LARGE=2;

    
    public static int Main_PLACEMENT_V03_align_scoreallnodes(File q, File db, File workDir) {

        try {
            
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            
            //preplacement parameters//////////////////////////////////////////
            //for a number of nodes =total_nodes/nodeShift , build diagsum vectors
            //think to keep sum of the diagsums to do mean at the end and highlight positions > to mean
            int nodeShift=100; // carefull, brings an error when nodeShift<2
            if (nodeShift<2) {nodeShift=2;}
            // pour stocker ou non dans hash
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            //word sampling method
            int queryWordSampling=SequenceKnife.SAMPLING_LINEAR;

            
            //debug/////////////////////////////////////////////////////////////
            //max number of queries treated 
            int queryLimit=100000000;
            //which log to write, !!!
            //more logs= much slower placement because of disk access latency

            //graph of words alignment
            boolean graphAlignment=false;
            

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOADING ALL DATA FROM THE INPUT DATABASE
            long startLoadTime=System.currentTimeMillis();
            
            //LOAD SESSION//////////////////////////////////////////////////////
            //logs
            String logPath=workDir+File.separator+"logs"+File.separator;
            //trees
            String relaxedTreePath=workDir+File.separator+"relaxed_trees"+File.separator;
            //ancestral reconstruciton
            String ARPath=workDir+File.separator+"AR"+File.separator;
            //session itself
            boolean loadHash=true;
            System.out.println("Loading ancestral words DB...");
            SessionNext session= SessionNext.load(db,loadHash);
            
            //type of Analysis//////////////////////////////////////////////////
            States s=session.states; 
            int analysisType=-1;
            //States: DNA or AA
            if (s instanceof DNAStates)
                analysisType=Main_PLACEMENT_V04_align_scoreallnodes_diagsumeltsremoved.TYPE_DNA;
            else if (s instanceof AAStates)
                analysisType=Main_PLACEMENT_V04_align_scoreallnodes_diagsumeltsremoved.TYPE_PROTEIN;
            
            //posterior probas parameters///////////////////////////////////////
            //mers size
            int k=session.k;
            int min_k=session.minK;
            Infos.println("k="+k);
            Infos.println("min_k="+min_k);
            //site and word posterior probas thresholds
            float sitePPThreshold=session.stateThreshold;
            float wordPPStarThreshold=session.wordThreshold;
            float factor=session.factor;
            float thresholdAsLog=(float)Math.log10(wordPPStarThreshold);

            Infos.println("factor="+factor);
            Infos.println("sitePPThreshold="+sitePPThreshold+" (for info, unused)");
            Infos.println("wordPPStarThreshold="+wordPPStarThreshold);
            Infos.println("wordPPStarThreshold(log10)="+thresholdAsLog);

            //PREPARE DIRECTORIES///////////////////////////////////////////////
            if (!workDir.exists()) {workDir.mkdir();}
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
            Infos.println("# nodes in the tree: "+tree.getNodeCount());
            Infos.println("# leaves in the tree: "+tree.getLeavesCount());
            Infos.println("# internal nodes in the tree: "+tree.getInternalNodesByDFS().size());
            String relaxedTreeForJplace=nw.getNewickTree(tree, true, true, true);
            //Infos.println(relaxedTreeForJplace);
            nw.close();
            
            //LOAD THE POSTERIOR PROBAS/////////////////////////////////////////
            PProbasSorted pprobas = session.parsedProbas;
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(pprobas.getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateIndexSet(0, 0, 5)));
            
            SimpleHash hash=session.hash;
            Infos.println(Environement.getMemoryUsage());
            long endLoadTime=System.currentTimeMillis();
            System.out.println("Loading the database took "+(endLoadTime-startLoadTime)+" ms");
            Environement.printMemoryUsageDescription();

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOAD FINISHED, TIME TO START THE PLACEMENT PROCESS
            

            long startTotalTime=System.currentTimeMillis();


            ////////////////////////////////////////////////////////////////////
            //MAP TO REGISTER PLACEMENTS: map(query)=(map(score)=positions_in_align)  !position in align, not diagsum
            HashMap<String,LinkedList<Score>> placementCompleteScore=new HashMap<>(); //fasta-> score=(score;PP*)

           ////////////////////////////////////////////////////////////////////
           //DIAGSUMS OUTPUT THE PLACEMENTS IN TSV FORMAT (completeScores)
           ////////////////////////////////////////////////////////////////////
           CSVWriter fwPlacement=new CSVWriter(new FileWriter(new File(logPath+"placements.tsv")), '\t');
           String[] header=new String[4];
           header[0]="Query";
           header[1]="NodeId";
           header[2]="ExtNodeId";
           header[2]="OriginalNodeId";
           header[3]="PP*";
           fwPlacement.writeNext(header);
            
            
            /////////////////////
            //PLACEMENT OF THE QUERIES
            FASTAPointer fp=new FASTAPointer(q, false);
            int totalQueries=fp.getContentSize();
            Infos.println("Input fasta contains "+totalQueries+" sequences");
            fp.resetPointer();
            //FileWriter fw =new FileWriter(new File(logPath+"queries.fasta"));
            
            int queryCounter=0;
            Fasta fasta=null;
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
               
                
                if (queryCounter>=queryLimit)
                    break;
                if ((queryCounter%10000)==0) {
                    System.out.println(queryCounter+"/"+totalQueries+" queries placed ("+(((0.0+queryCounter)/totalQueries)*100)+"%)");
                }
                
                Infos.println("#######################################################################");
                Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+fasta.getHeader());
                Infos.println("#######################################################################");
                //fw.append(fasta.getFormatedFasta()+"\n");
                long startScanTime=System.currentTimeMillis();
                int queryLength=fasta.getSequence(false).length();
                Infos.println("Query length: "+queryLength);
                

                
                



                
                ///////////////////////////////////
                // LOOP ON QUERY WORDS
                QueryWord qw=null;
                SequenceKnife sk=new SequenceKnife(fasta, k, min_k, s, queryWordSampling);
                //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                int queryWordCounter=0;
                int queryWordFoundCounter=0;
                
               
                //a diagsum based only on a the single top best PP* associated to a query word.
                //accumulating all nodes, just taing higher PP* and corresponding positions
                DiagSum bestPPStarsDiagsum=new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep());
                

                //preparing preplacement graph data :
                double[][] graphDataForTopTuples =null;
                int xSize=bestPPStarsDiagsum.getSize();
                int ySize=queryLength;
                
                if (graphAlignment) {
                    graphDataForTopTuples =new double[3][xSize*ySize];
                    //init the Z axis (PP*) to very small values for all possible (X,Y)
                        for (int y = 0; y < ySize; y++) {
                            for (int x = 0; x < xSize; x++) {
                                //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
                                graphDataForTopTuples[0][x+(y*xSize)]=x;
                                graphDataForTopTuples[1][x+(y*xSize)]=y;
                                graphDataForTopTuples[2][x+(y*xSize)]=thresholdAsLog;
                            }
                        }
                }
                    
                
                
                

                // MATRIX DESCRIBING ALIGNMENT
                //////////////////////////////////////////
                
                //[x][y] = [queryPosition][n_ieme_hit_in_alignment]
                ArrayList<Integer>[] alignmentMatrix=new ArrayList[queryLength-k+1];
                for (int i = 0; i < alignmentMatrix.length; i++) {
                    alignmentMatrix[i]=new ArrayList<>(5);
                }
                //Arrays.fill(alignmentMatrix,new ArrayList<Integer>()); //not working, ask why in internet
                
                // LIST OF NODES ENCOUNTERED DURING ALIGNMENT, USED FOR LATER PLACEMENT
                ///////////////////////////////////////////////////////
                boolean[] nodesEncountered=new boolean[tree.getNodeCount()]; //tab[score]=true/false
                int encounteredCount=0;
                
                
                
                //ON TOP OF DIAGSUM APPROACH, BUILD THE REFERENCE COMPLETE SCORE
                //////////////////////////////////////////////////////////////////
                //when all tuples of all nodes are parsed, to compare these placements.
                //here, all alignmnet built from topTuples is used to define the ref positions
                //then ALL tuples corresponding to this position are read
                //if present sum the allTuples PP* t, if not, sum PP* threshold 
                Float[] completeScore=new Float[tree.getNodeCount()];
                float baseline=thresholdAsLog*sk.getMerOrder().length;
                Arrays.fill(completeScore,baseline);
                //completeScore[score]=baseline, for now do not take into account 
                //correction do to partial overlap should be later
                
                
                
                // ITERTION ON ALL QUERY WORDS, TO BUILD THE ALIGNMENT
                ///////////////////////////////////////////////////////
                
                while ((qw=sk.getNextWord())!=null) {
                    //Infos.println("Query mer: "+qw.toString());


                    SimpleHash.Tuple topTuple = hash.getTopTuple(qw);
                    //if this word is not registered in the hash
                    if (topTuple==null) {
                        queryWordCounter++;
                        continue;
                    }
                    //word is in the hash
                    queryWordFoundCounter++;
                    
  
                    
                    //get best PP* of each query word (whatever node/position)
                    //////////////////////////////////////////////////////////
                    
                    //a diagsum based only on a the single top best PP* associate dto a word.
                    int topDiagSumPos=topTuple.getRefPos()-qw.getOriginalPosition()+(queryLength-minOverlap);
                    if (topDiagSumPos>-1 && topDiagSumPos<bestPPStarsDiagsum.getSize()) {
                        if(graphAlignment) {
                            graphDataForTopTuples[2][queryWordCounter*xSize+topTuple.getRefPos()]= topTuple.getPPStar();                        
                        }
                        alignmentMatrix[qw.getOriginalPosition()].add(topTuple.getRefPos());
                        //System.out.println("nodesEncountered[topTuple.getNodeId()]"+nodesEncountered[topTuple.getNodeId()]);
                        if (nodesEncountered[topTuple.getNodeId()]==false)
                            encounteredCount++;
                        nodesEncountered[topTuple.getNodeId()]=true;
                    }
                    
                    
                    queryWordCounter++;
                    
                    //DEBUG
                    //if (queryWordCounter>1000)
                    //        break;
                    //DEBUG
                    
                }
                
                //check alignment vs diagsum stats
//                for (int i = 0; i < alignmentMatrix.length; i++) {
//                    System.out.println("alignment i="+i+": "+alignmentMatrix[i]);
//                }
                //check diagsum stats
//                for (int i = 0; i < bestPPStarsPerPositionDiagsum.getSize(); i++) {
//                    System.out.println("pre-diagsum i="+i+": #words="+bestPPStarsPerPositionDiagsum.getEffectiveWordCount(i)+" PP*="+bestPPStarsPerPositionDiagsum.getSum(i));
//                }

                
                if (graphAlignment) {
                    //graph for topTuplePerPos
                    DefaultXYZDataset datasetForGraph=new DefaultXYZDataset();
                    //datasetForGraph.addSeries(0, graphDataForTopTuples);
                    datasetForGraph.addSeries(0, graphDataForTopTuples);
                    JFrame infos=new JFrame();
                    infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                    infos.setLayout(new GridLayout(1, 1));
                    infos.add(ChartsForNodes.buildReadMatchForANode4("Top_tuple_per_pos of read="+fasta.getHeader(),bestPPStarsDiagsum.getSize(),queryLength, datasetForGraph,thresholdAsLog,0));
                    infos.setSize(1024, 250);
                    infos.pack();
                    RefineryUtilities.centerFrameOnScreen(infos);
                    infos.setVisible(true);
                }
               

                Infos.println("Proportion of query words retrieved in the hash: "+queryWordFoundCounter+"/"+queryWordCounter);
                long endScanTime=System.currentTimeMillis();
                Infos.println("Alignment ("+encounteredCount+" candidate nodes registered) took "+(endScanTime-startScanTime)+" ms");
           
                

                // HERE reuse the alignment vector to determine which positions 
                //should be tested for all nodes
                
                //load the PP* of these positions nodes per nodes
                //build the diagsums per nodes 
                //select peeks including #words > threshold
                //select the node with the best peeks as the best placement
                
                startScanTime=System.currentTimeMillis();
                //to register which node is assined to which line of preplacementAllDiagsum
                Infos.println("Launching search on candidate nodes...");

                
                
                //now, pass again on all algnment positions, but search them in 
                //all internal nodes.
                //only alignment positions whose coordinates matches a
                //diagsumPos associated to more than (20-k+1)/s words
                //(20 consecutive represented by consecutive words)
                //will be tested
                
                

                
                //V2, pickup nodes directly from tuples corrsponding to position, avoiding to check all nodes
                for (int alignPos = 0; alignPos < alignmentMatrix.length; alignPos++) {
                    boolean firstOfResPos=true;
                    QueryWord word=sk.getWordAt(alignPos);
                    if (alignmentMatrix[alignPos].size()>0) {
                        Integer bestRefPos=alignmentMatrix[alignPos].get(0);
                        List<SimpleHash.Tuple> allTuples = hash.getTuplePerRefPosition(word, bestRefPos);
                        //System.out.println("allTuples:"+allTuples);
                        for (int i = 0; i < allTuples.size(); i++) {
                            SimpleHash.Tuple tuple = allTuples.get(i);
                            //completeScore, based only on 1st refPos (eq. topTuple)
                            if (firstOfResPos) {
                                //System.out.println("completeScore: refPos="+refPos+" score:"+tuple.getNodeId()+" "+tuple.getPPStar());
                                completeScore[tuple.getNodeId()]+=-thresholdAsLog+tuple.getPPStar();
                                break;
                            }
                        }
                        firstOfResPos=false;
                    }
                }            
                
                
                endScanTime=System.currentTimeMillis();
                Infos.println("Targeted scan on "+tree.getInternalNodesByDFS().size()+" nodes took "+(endScanTime-startScanTime)+" ms");
                
                
                
                //check the completeScore content
                placementCompleteScore.put(fasta.getHeader(), new LinkedList<>());
                for (int i = 0; i < completeScore.length; i++) {
                    if (completeScore[i]!=baseline) {
                        placementCompleteScore.get(fasta.getHeader()).add(new Score(i, completeScore[i]));
                    }
                }
                Collections.sort(placementCompleteScore.get(fasta.getHeader()));
                
                String[] data=new String[4];

               data[0]=fasta.getHeader();
               int n=0;
               for (Iterator<Score> iterator1= placementCompleteScore.get(fasta.getHeader()).iterator(); iterator1.hasNext();) {
                   Score nextNodePlacement = iterator1.next();
                   data[1]=String.valueOf(nextNodePlacement.getNodeId());
                   data[2]=String.valueOf("none");//tree.getById(nextNode).getExternalId());
                   data[3]=String.valueOf(nextNodePlacement.getPPStar());
                   fwPlacement.writeNext(data);
                   n++;
                   if (n>5)
                       break;
               }



                
                
                queryCounter++;
            }
            
            fwPlacement.close();
            fp.closePointer();
            
            
            

            long endTotalTime=System.currentTimeMillis();
            System.out.println("Placement took "+(endTotalTime-startTotalTime)+" ms");
            Infos.println("#######################################################################");
            Infos.println("### DONE, placement execution took (excluding DB load): "+(endTotalTime-startTotalTime)+" ms");
            Infos.println("#######################################################################");
            Infos.println("### "+queryCounter+" READS WERE ANALYZED");
            Infos.println("#######################################################################");
            
            
//            for (int i = 0; i < tree.getNodeIdsByDFS().size(); i++) {
//                System.out.println("nodeId: "+tree.getNodeIdsByDFS().get(i)+
//                        "  -->  FakeToOriginal: "+((ExtendedTree)tree).getFakeToOriginalId(tree.getNodeIdsByDFS().get(i)));
//                
//            }
//            
            
            


            
            

            

            
            //tree.displayTree();
            
            
            //System.exit(0);
            return queryCounter;
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PLACEMENT_V04_align_scoreallnodes_diagsumeltsremoved.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
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
