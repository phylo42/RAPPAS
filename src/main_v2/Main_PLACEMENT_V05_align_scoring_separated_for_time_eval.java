/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import alignement.Alignment;
import au.com.bytecode.opencsv.CSVWriter;
import charts.ChartsForNodes;
import core.AAStates;
import core.DNAStates;
import core.DiagSum;
import core.PProbasSorted;
import core.QueryWord;
import core.Score;
import core.States;
import core.algos.SequenceKnife;
import core.hash.Pair;
import core.hash.CustomHash_v2;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.awt.GridLayout;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import jonelo.jacksum.JacksumAPI;
import jonelo.jacksum.algorithm.AbstractChecksum;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.ui.RefineryUtilities;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import tree.ExtendedTree;
import tree.NewickReader;
import tree.NewickWriter;
import tree.PhyloTree;

/**
 * Algo in this version:
 * 
 * 1. all words of the query are searched in the hash
 *    --> we storeFullHash the corresponding nodes in a list
    --> positions are used to define the alignment
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
public class Main_PLACEMENT_V05_align_scoring_separated_for_time_eval {
    
    //sequence type
    public static int TYPE_DNA=1;
    public static int TYPE_PROTEIN=2;
    //memory mode
    public static int MEMORY_LOW=1;
    public static int MEMORY_LARGE=2;
    
    
    public static int Main_PLACEMENT_V05_align_scoreallnodes(File q, File db, File workDir) {

        try {
                        
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
              
            // pour stocker ou non dans hash
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            //word sampling method
            int queryWordSampling=SequenceKnife.SAMPLING_LINEAR;

            
            //debug/////////////////////////////////////////////////////////////
            //max number of queries treated 
            int queryLimit=100000;
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
            SessionNext_v2 session= SessionNext_v2.load(db,loadHash);
            
            //type of Analysis//////////////////////////////////////////////////
            States s=session.states; 
            int analysisType=-1;
            //States: DNA or AA
            if (s instanceof DNAStates)
                analysisType=Main_PLACEMENT_V05_align_scoring_separated_for_time_eval.TYPE_DNA;
            else if (s instanceof AAStates)
                analysisType=Main_PLACEMENT_V05_align_scoring_separated_for_time_eval.TYPE_PROTEIN;
            
            //posterior probas parameters///////////////////////////////////////
            //mers size
            int k=session.k;
            int min_k=session.minK;
            Infos.println("k="+k);
            Infos.println("min_k="+min_k);
            //site and word posterior probas thresholds
            float sitePPThreshold=session.stateThreshold;
            float wordPPStarThreshold=session.PPStarThreshold;
            float factor=session.alpha;
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
            PhyloTree originalTree=session.originalTree;
            //write a newick version of original tree, setting the nodeId
            //as internal node labels 
            NewickWriter nw=new NewickWriter(new File(logPath+"tree_with_edge_labels.nwk"));
            nw.writeNewickTree(originalTree, true, true, true);
            PhyloTree ARtree = session.ARTree;
            ExtendedTree extendedtree = session.extendedTree;

            
            Infos.println("# nodes in the tree: "+ARtree.getNodeCount());
            Infos.println("# leaves in the tree: "+ARtree.getLeavesCount());
            Infos.println("# internal nodes in the tree: "+ARtree.getInternalNodesByDFS().size());
            String relaxedTreeForJplace=nw.getNewickTree(originalTree, true, true, true);
            //Infos.println(relaxedTreeForJplace);
            nw.close();
            
            //LOAD THE POSTERIOR PROBAS/////////////////////////////////////////
            PProbasSorted pprobas = session.parsedProbas;
            //to raidly check that sorted probas are OK
            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(pprobas.getPPSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateSet(0, 0, 5)));
            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateIndexSet(0, 0, 5)));
            
            CustomHash_v2 hash=session.hash;
            Infos.println(Environement.getMemoryUsage());
            long endLoadTime=System.currentTimeMillis();
            System.out.println("Loading the database took "+(endLoadTime-startLoadTime)+" ms");
            Environement.printMemoryUsageDescription();

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOAD FINISHED, TIME TO START THE PLACEMENT PROCESS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            //global timer for align + scoring
            long startTotalTime=System.currentTimeMillis();

            ////////////////////////////////////////////////////////////////////
            //MAP TO REGISTER PLACEMENTS: map(query)=(map(score)=positions_in_align)  !position in align, not diagsum
            HashMap<String,LinkedList<Score>> placementCompleteScore=new HashMap<>(); //fasta-> score=(score;PP*)

            ////////////////////////////////////////////////////////////////////
            //TO OUTPUT THE PLACEMENTS IN TSV FORMAT (completeScores)
            ////////////////////////////////////////////////////////////////////
            int bufferSize=2097152; // buffer of 2mo
            BufferedWriter fwPlacement=new BufferedWriter(new FileWriter(new File(logPath+"placements.tsv")),bufferSize);
            StringBuffer sb=new StringBuffer("Query\tARTree_NodeId\tARTree_NodeName\tExtendedTree_NodeId\tARTree_NodeName\tOriginal_NodeId\tOriginal_NodeName\tPP*\n");
            
            
            /////////////////////
            //PLACEMENT OF THE QUERIES
            FASTAPointer fp=new FASTAPointer(q, false);
            int totalQueries=fp.getContentSize();
            Infos.println("Input fasta contains "+totalQueries+" sequences");
            fp.resetPointer();
            //FileWriter fw =new FileWriter(new File(logPath+"queries.fasta"));
            
            
            ///////////////////////////////////////////////////////            
            // VECTORS USED TO ALIGN AND SCORE NODES
            
            //list of variable size, reset at each read
            HashMap<Integer,Boolean> selectedNodes=new HashMap<>(); 
            //instanciated once
            int[] nodeOccurences=new int[ARtree.getNodeCount()]; // tab[#times_encoutered] --> index=nodeId
            float[] nodeScores=new float[ARtree.getNodeCount()]; // tab[score] --> index=nodeId
            //Arrays.fill(nodeOccurences, 0);
            //Arrays.fill(nodeScores, 0.0f);     
            
            
            ///////////////////////////////////////////////////////////////////       
            // PREPARE BUFFER FOR CSV FILE WRITER
            
            ///////////////////////////////////////////////////////////////////       
            // PREPARE CHECKSUM FOR IDENTICAL READS REGISTER
            AbstractChecksum checksumGenerator=null;
            try {
                checksumGenerator = JacksumAPI.getChecksumInstance("sha256");
                checksumGenerator.setEncoding(AbstractChecksum.HEX);
            } catch (NoSuchAlgorithmException ex) {
                Logger.getLogger(Main_PLACEMENT_V05_align_scoring_separated_for_time_eval.class.getName()).log(Level.SEVERE, null, ex);
            }
            HashMap<String,ArrayList<String>> identicalSeqsRegistry=new HashMap<>();
            
            ////////////////////////////////////////////////////////////////////
            //PREPARE STRUCTURE (JSON OBJECT) FOR JPLACE OUTPUT 
            //this need to be done outside the alignment/scoring loop
            //because identical sequences (same score) will be in the same 
            //json "placement" (p) object
            //with library json.simple
            
            //map to associate checksums to JSONObject
            HashMap<String,JSONObject> checksumToJSONObject=new HashMap<>();
            
            //top level of the json tree
            JSONObject top=new JSONObject(); 
            LinkedHashMap topMap=new LinkedHashMap();
            //object tree (mandatory)
            topMap.put("tree",relaxedTreeForJplace);  
            //we do an array of placement object
            //all identical reads with be injected in the same object
            JSONArray placements=new JSONArray();
            
            
            ////////////////////////////////////////////////////////////////////
            // MAIN ALGORITHM BELOW !!!
            ////////////////////////////////////////////////////////////////////
            
            int queryCounter=0;
            long totalAlignTime=0;
            long totalScoringTime=0;
            long totalWritingTime=0;
            long totalResetTime=0;
            long totalChecksumTime=0;
            
            Fasta fasta=null;
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
               
                //debug
                if (queryCounter>=queryLimit)
                    break;
                //console display to follow the process
                if ((queryCounter%10000)==0) {
                    System.out.println(queryCounter+"/"+totalQueries+" queries placed ("+(((0.0+queryCounter)/totalQueries)*100)+"%)");
                }
                //checksum built, if already present do not compute placement
                //again. this might need compression if fasta sequences headers
                //are too heavy in memory
                long startChecksumTime=System.currentTimeMillis();
                checksumGenerator.reset(); //make it ready before next checksum computation
                checksumGenerator.update(fasta.getSequence(false).getBytes());
                String checksum = checksumGenerator.getFormattedValue();
                int cutIndex=fasta.getHeader().indexOf(" ");
                if (cutIndex<0) //basically, space not found
                    cutIndex=fasta.getHeader().length();
                String subHeader=fasta.getHeader().substring(0,cutIndex);
                if (identicalSeqsRegistry.containsKey(checksum)) {
                    identicalSeqsRegistry.get(checksum).add(subHeader);
                    Infos.println("! SKIPPED BECAUSE DUPLICATE: "+fasta.getHeader());
                    queryCounter++;
                    //before skipping, update the outputs
                    //for now, only jplace is done here...
                    //note that detailed comments about the jplace are
                    //in the jplace block on the bottom of the loop
                    //get back the placement out from the JSONObject
                    JSONObject placement = checksumToJSONObject.get(checksum);
                    //the "p" object values are the same, no changes
                    JSONArray pMetadata=(JSONArray)placement.get("p");
                    //the "nm" object has to be extended with the identifier
                    //and multiplicity of this read
                    JSONArray allIdentifiers=(JSONArray)placement.get("nm");
                    JSONArray readMultiplicity=new JSONArray();
                    readMultiplicity.add(fasta.getHeader());
                    readMultiplicity.add(1);
                    allIdentifiers.add(readMultiplicity);                    
                    
                    continue;
                } else {
                    ArrayList<String> a=new ArrayList<>();
                    a.add(subHeader);
                    identicalSeqsRegistry.put(checksum,a);
                }
                long endChecksumTime=System.currentTimeMillis();
                totalChecksumTime+=endChecksumTime-startChecksumTime;
                
                
                Infos.println("#######################################################################");
                Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+fasta.getHeader());
                Infos.println("#######################################################################");
                //fw.append(fasta.getFormatedFasta()+"\n");
                long startAlignTime=System.currentTimeMillis();
                int queryLength=fasta.getSequence(false).length();
                Infos.println("Query length: "+queryLength);
                
                
                ///////////////////////////////////
                // PREPARE QUERY K-MERS
                QueryWord qw=null;
                SequenceKnife sk=new SequenceKnife(fasta, k, min_k, s, queryWordSampling);
                //Infos.println("Mer order: "+Arrays.toString(rk.getMerOrder()));
                int queryWordCounter=0;
                int queryWordFoundCounter=0;
                
  
                
                
                //if alignment graph are required
                DiagSum bestPPStarsDiagsum=null;
                double[][] graphDataForTopTuples =null;
                int xSize=-1;
                if (graphAlignment) {
                    //a diagsum based only on a the single top best PP* associated to a query word.
                    //accumulating all nodes, just taing higher PP* and corresponding positions
                    bestPPStarsDiagsum=new DiagSum(queryLength, align.getLength(), minOverlap, k, sk.getStep());
                    //preparing preplacement graph data :
                    xSize=bestPPStarsDiagsum.getSize();
                    int ySize=queryLength;
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
                
                
                // PREPARE MATRIX DESCRIBING ALIGNMENT
                //////////////////////////////////////////
                
                //[query_position]=reference_position
                int[] alignmentMatrix=new int[queryLength-k+1];
                //filled with -1 ; a -1 val means the word was not found in hash
                Arrays.fill(alignmentMatrix,-1);
                
                
                // BUILD THE ALIGNMENT
                ///////////////////////////////////////////////////////

                while ((qw=sk.getNextWord())!=null) {
                    //Infos.println("Query mer: "+qw.toString());

                    //get top Pair_16_32_bit associated to this word
                    Pair topPair = hash.getTopPair(qw.getWord());
                    //if this word is not registered in the hash
                    if (topPair==null) {
                        queryWordCounter++;
                        continue;
                    }
                    //word is in the hash
                    queryWordCounter++;
                    queryWordFoundCounter++;
                    //top position of this word
                    int topPosition=hash.getPositions(qw.getWord())[0];
                    //graph of alignment, if asked
                    if(graphAlignment) {
                        int topDiagSumPos=topPosition-qw.getOriginalPosition()+(queryLength-minOverlap);
                        if (topDiagSumPos>-1 && topDiagSumPos<bestPPStarsDiagsum.getSize()) {
                                graphDataForTopTuples[2][queryWordCounter*xSize+topPosition]= topPair.getPPStar();                        
                        }
                    }
                    //register the alignment itself
                    alignmentMatrix[qw.getOriginalPosition()]=topPosition;
                    //register nodes that were hit during alignment (only those will be scored)
                    selectedNodes.put(topPair.getNodeId(), true);
                    
                    //DEBUG
                    //if (queryWordCounter>1000)
                    //        break;
                    //DEBUG
                    
                }
                
                //check alignment vs diagsum stats
//                for (int i = 0; i < alignmentMatrix.length; i++) {
//                    System.out.println("alignment i="+i+": "+alignmentMatrix[i]);
//                }

                //display alignment result of requested
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
                long endAlignTime=System.currentTimeMillis();
                totalAlignTime+=(endAlignTime-startAlignTime);
                Infos.println("Alignment ("+selectedNodes.keySet().size()+" candidate nodes registered) took "+(endAlignTime-startAlignTime)+" ms");
           
                

                // ALIGNMENT DONE, NOW SCORING SELECTED NODES
                ///////////////////////////////////////////////////////

                long startScoringTime=System.currentTimeMillis();
                //to register which node is assined to which line of preplacementAllDiagsum
                Infos.println("Launching scoring on candidate nodes...");
                Infos.println("Candidate nodes: ("+selectedNodes.keySet().size()+") "+selectedNodes.keySet());
                //V3, pickup PP* associated to position defined by the alignment
                //and only for the nodes encountered in the alignment phase
                for (int queryPos = 0; queryPos < alignmentMatrix.length; queryPos++) {
                    //if this position was indeed aligned
                    int refPos;
                    if ((refPos=alignmentMatrix[queryPos])>-1) {
                        QueryWord word=sk.getWordAt(queryPos);
                        //pairs of (nodeId,PP*)
                        List<Pair> allPairs = hash.getPairs(word.getWord(), refPos);
                        for (int i = 0; i < allPairs.size(); i++) {
                            Pair p = allPairs.get(i);
                            //score only encountered nodes
                            if (selectedNodes.get(p.getNodeId())!=null) {
                                nodeOccurences[p.getNodeId()]+=1;
                                nodeScores[p.getNodeId()]+=p.getPPStar();
                            }
                        }
                    //position not aligned, 
                    } else {
                        
                    }
                
                }            
                
                //now add the score corresponding to the words not found,
                // i.e. threshold*#_words_not_scored (because no in hash)
                int maxWords=sk.getMaxMerCount();

                float bestScore=Float.NEGATIVE_INFINITY;
                int bestNode=-1;
                for (Iterator<Integer> iterator = selectedNodes.keySet().iterator(); iterator.hasNext();) {
                    Integer nodeId = iterator.next();
                    nodeScores[nodeId]+=thresholdAsLog*(maxWords-nodeOccurences[nodeId]);
                    if(nodeScores[nodeId]>bestScore) {
                        bestScore=nodeScores[nodeId];
                        bestNode=nodeId;
                    }
                }
                Infos.println("Best node (ARTree) is : "+bestNode+" (score="+bestScore+")");                
                
                long endScoringTime=System.currentTimeMillis();
                totalScoringTime+=endScoringTime-startScoringTime;
                
                
                //TO DEBUG
//                System.out.println(fasta.getHeader());
//                System.out.print("bestNode:"+bestNode);
//                System.out.print("\t\tbestScore:"+bestScore);
//                System.out.println("\t\tbestSCore (normalized by #mers):"+(bestScore/sk.getMaxMerCount()));
//                System.out.println("nodeScores:"+Arrays.toString(nodeScores));
//                System.out.println("nodeOccurences:"+Arrays.toString(nodeOccurences));
//                System.out.println("selectedNodes:"+selectedNodes.keySet().toString());
                
                
                //write result in file if a node was hit
                if (bestNode>-1) {
                    long startWritingTime=System.currentTimeMillis();
                    sb.append(fasta.getHeader().split(" ")[0]+"\t");
                    sb.append(String.valueOf(bestNode)+"\t"); //ARTree nodeID
                    sb.append(String.valueOf(session.ARTree.getById(bestNode).getLabel())+"\t"); //ARTree nodeName
                    if (session.nodeMapping.get(bestNode)==null)
                        System.out.println("bestNode not found: "+bestNode+" "+String.valueOf(session.ARTree.getById(bestNode).getLabel()));
                    int extendedTreeId=session.nodeMapping.get(bestNode);
                    sb.append(String.valueOf(extendedTreeId)+"\t"); //extended Tree nodeID
                    sb.append(String.valueOf(session.extendedTree.getById(extendedTreeId).getLabel())+"\t"); //extended Tree nodeName
                    //check if this was a fake node or not
                    Integer originalNodeId = extendedtree.getFakeToOriginalId(extendedTreeId);//will return null if this nodeId was not a Fake node
                    if (originalNodeId==null) {originalNodeId=extendedTreeId;}
                    sb.append(String.valueOf(originalNodeId)+"\t"); //edge of original tree (original nodeId)
                    sb.append(String.valueOf(session.originalTree.getById(originalNodeId).getLabel())+"\t"); //edge of original tree (original nodeName
                    sb.append(String.valueOf(nodeScores[bestNode])+"\n");
                    long endWritingTime=System.currentTimeMillis();
                    totalWritingTime+=endWritingTime-startWritingTime;
                }
                
                //push the stringbuffer to the bufferedwriter every 10000 sequences
                if ((queryCounter%10000)==0) {
                    int size=sb.length();
                    fwPlacement.append(sb);
                    fwPlacement.flush();
                    sb=null;
                    sb=new StringBuffer(size);
                }                
                

                //fill the JSON placement object
                //2 possibilities:
                //-either do a new placement object and add ot to the list
                // of placements (block above, at checksum test)
                //-or add it to previous placement object if this is a duplicate
                // sequence
                
                //write result in file if a node was hit
                if (bestNode>-1) {
                    //object representing placements (mandatory), it contains 2 key/value couples
                    JSONObject placement=new JSONObject();
                    //first we build the "p" array, containing the position/scores of all reads
                    JSONArray pMetadata=new JSONArray();
                    //in pplacer/EPA several placements can be associated to a query
                    //we input only the best one, but that can be changed in the future
                    JSONArray placeColumns=new JSONArray();

                    //fake fields for compatibility with current tools (guppy, archeopteryx)
                    //should be provided as an option
                    placeColumns.add(0.1);
                    placeColumns.add(0.1);
                    placeColumns.add(0.1);


                    //placeColumns.add(session.ARTree.getById(bestNode).getLabel()); // 1. ARTree nodeName
                    int extendedTreeId=session.nodeMapping.get(bestNode);
                    //placeColumns.add(session.extendedTree.getById(extendedTreeId).getLabel()); // 2. extended tree nodeName
                    //check if this was a fake node or not
                    Integer originalNodeId = extendedtree.getFakeToOriginalId(extendedTreeId);//will return null if this nodeId was not a Fake node
                    if (originalNodeId==null) {originalNodeId=extendedTreeId;}
                    placeColumns.add(originalNodeId); // 3. edge of original tree (original nodeId=edgeID)
                    //placeColumns.add(session.originalTree.getById(originalNodeId).getLabel()); // 4. edge of original tree (original nodeName)
                    placeColumns.add(nodeScores[bestNode]); // 4. PP*
                    pMetadata.add(placeColumns);

                    placement.put("p", pMetadata);

                    //second we build the "nm" array, containing the reads identifiers
                    //And their read multiplicity
                    JSONArray allIdentifiers=new JSONArray();
                    //these are simply all the identical sequences
                    JSONArray readMultiplicity=new JSONArray();
                    readMultiplicity.add(fasta.getHeader());
                    readMultiplicity.add(1);
                    allIdentifiers.add(readMultiplicity);
                    placement.put("nm", allIdentifiers);

                    //store the placement in the list of placements
                    placements.add(placement);
                    //JSON for this read DONE, if duplicates are found later,
                    //will be added to the corresponding "p" and "nm" array
                    //using the checksumToJSONObject map
                    //for now, just register the reference
                    checksumToJSONObject.put(checksum, placement);                
                }

                //reset the scoring vectors
                long startResetTime=System.currentTimeMillis();
                for (Iterator<Integer> iterator = selectedNodes.keySet().iterator(); iterator.hasNext();) {
                    Integer nodeId = iterator.next();
                    nodeScores[nodeId]=0.0f;
                    nodeOccurences[nodeId]=0;
                }
                selectedNodes.clear();
                long endResetTime=System.currentTimeMillis();
                totalResetTime+=endResetTime-startResetTime;
                
                queryCounter++;
            }
            
            //just for coherent output, close the percentage
            System.out.println(queryCounter+"/"+totalQueries+" queries placed ("+(((0.0+queryCounter)/totalQueries)*100)+"%)");
            
            //flush to disk the last CSV buffer
            fwPlacement.append(sb);
            fwPlacement.close();
            fp.closePointer();
            
            //FINISH the json strucutre and output it to a file
            //associate the list of placements to the top level
            top.put("placements",placements);
            //object version (mandatory
            top.put("version",3);
            //object metadata (mandatory): sub-objet "version is mandatory"
            JSONObject invoc=new JSONObject();
            invoc.put("invocation", "viromeplacer xxxxxxxxxxxxxxx");
            top.put("metadata", invoc);
            //object fields
            //for info:
            //- in pplacer: "distal_length", "edge_num", "like_weight_ratio", "likelihood", "pendant_length"
            //- in EPA: "edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"
            JSONArray fList=new JSONArray();
            //add fake fields to be compatible with current visualisation tools
            fList.add("distal_length");
            fList.add("like_weight_ratio");
            fList.add("pendant_length");
            
            //fList.add("ARTree_nodeName");
            //fList.add("ExtendedTree_nodeName");
            fList.add("edge_num"); //i.e equal to the id of the son node
            //fList.add("edge_label"); //i.e equal to the id of the son node
            fList.add("likelihood"); //rename to likelihood even if it is not, but for compatibility with other programs
            top.put("fields", fList);
            
            //put all the elements in the top JSON object
            top.putAll(topMap);
            String out=top.toJSONString();

            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\],\n"); //]}]
            //out=out.replace("]},", "]},"); //]}
            
            FileWriter fwJSON =new FileWriter(new File(logPath+File.separator+"placements.jplace"));
            fwJSON.append(out);
            fwJSON.close();
            
            
            
  

            long endTotalTime=System.currentTimeMillis();
            System.out.println("############################################################");
            System.out.println("Checksum registry took in total: "+totalChecksumTime+" ms");
            System.out.println("Alignments took in total: "+totalAlignTime+" ms");
            System.out.println("Scoring took in total: "+totalScoringTime+" ms");
            System.out.println("Writing CSV took in total: "+totalWritingTime+" ms");
            System.out.println("Reset took in total: "+totalResetTime+" ms");
            System.out.println("------------------------------------------------------------");
            System.out.println("Checksum registry took on average: "+(totalChecksumTime/queryCounter)+" ms");
            System.out.println("Alignments took on average: "+(totalAlignTime/queryCounter)+" ms");
            System.out.println("Scoring took on average: "+(totalScoringTime/queryCounter)+" ms");
            System.out.println("Writing CSV took on average: "+(totalWritingTime/queryCounter)+" ms");
            System.out.println("Reset took on average: "+(totalResetTime/queryCounter)+" ms");
            System.out.println("------------------------------------------------------------");
            System.out.println("Process (without DB load) took: "+(endTotalTime-startTotalTime)+" ms");
            System.out.println("############################################################");
            Infos.println("#######################################################################");
            Infos.println("### DONE, placement execution took (excluding DB load): "+(endTotalTime-startTotalTime)+" ms");
            Infos.println("#######################################################################");
            Infos.println("### "+queryCounter+" READS WERE ANALYZED");
            Infos.println("#######################################################################");
            
            
//            for (int i = 0; i < ARTree.getNodeIdsByDFS().size(); i++) {
//                System.out.println("nodeId: "+ARTree.getNodeIdsByDFS().get(i)+
//                        "  -->  FakeToOriginal: "+((ExtendedTree)ARTree).getFakeToOriginalId(ARTree.getNodeIdsByDFS().get(i)));
//                
//            }
//            
            



            
            

            

            
            //tree.displayTree();
            
            
            //System.exit(0);
            return queryCounter;
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PLACEMENT_V05_align_scoring_separated_for_time_eval.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
        }
        
        
        

    }
    

    
    

    
    
    
    
    
}
