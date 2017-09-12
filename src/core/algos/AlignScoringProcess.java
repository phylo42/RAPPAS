/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import charts.ChartsForNodes;
import com.google.common.math.Quantiles;
import core.DiagSum;
import core.QueryWord;
import core.hash.Pair;
import etc.Infos;
import inputs.Fasta;
import inputs.SequencePointer;
import java.awt.GridLayout;
import java.io.BufferedWriter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import jonelo.jacksum.JacksumAPI;
import jonelo.jacksum.algorithm.AbstractChecksum;
import main_v2.Main_PLACEMENT_v07;
import main_v2.SessionNext_v2;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.ui.RefineryUtilities;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import tree.PhyloNode;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class AlignScoringProcess {
    
    
    //debug/////////////////////////////////////////////////////////////
    //max number of queries treated 
    int queryLimit=10000;

    //graph of words alignment
    boolean graphAlignment=false;
    
    //parameters asked when program launched
    SessionNext_v2 session=null;
    Float nsBound=Float.NEGATIVE_INFINITY;
    
    //elements related to fil outputs
    //filled while the placement process is running
    JSONArray placements=null;
    BufferedWriter bwTSV=null;
    StringBuffer sb=null;
    
    /**
     *
     * @param session
     * @param nsBound
     * @param queryLimit the value of queryLimit
     */
    public AlignScoringProcess(SessionNext_v2 session, Float nsBound, int queryLimit) {
        this.session=session;
        this.nsBound=nsBound;
        this.queryLimit=queryLimit;
    }
    
    /**
     * 
     * @param rs
     * @param samplingAmount
     * @param bwTSV
     * @param queryWordSampling
     * @param minOverlap
     * @param q_quantile
     * @param n_quantile
     * @return 
     * @throws java.io.IOException 
     */
    public float processCalibration(RandomSeqGenerator rs, int samplingAmount,BufferedWriter bwTSV, int queryWordSampling, int minOverlap, int q_quantile, int n_quantile) throws IOException {

        //list of normalizedScore
        ArrayList<Float> normalizedScores=new ArrayList<>(samplingAmount);

        ///////////////////////////////////////////////////////            
        // PREPARE VECTORS USED TO ALIGN AND SCORE NODES

        //list of nodes encountered during matches search in the hash
        //variable size, expanded only at 1st encounter with the node
        //when nodeOccurences[nodeId]==0
        //,reset at each read
        ArrayList<Integer> selectedNodes=new ArrayList<>(session.ARTree.getNodeCount());
        //instanciated once
        int[] nodeOccurences=new int[session.ARTree.getNodeCount()]; // tab[#times_encoutered] --> index=nodeId
        float[] nodeScores=new float[session.ARTree.getNodeCount()]; // tab[score] --> index=nodeId
        //Arrays.fill(nodeOccurences, 0);
        //Arrays.fill(nodeScores, 0.0f);     

        
        
        //////////////////////////////////////////////////////////////////
        // PREPARE TSV OUTPUT
        //header of CSV output
        if (bwTSV!=null) {
            sb=new StringBuffer("Query\tARTree_NodeId\tARTree_NodeName\tExtendedTree_NodeId\tARTree_NodeName\tOriginal_NodeId\tOriginal_NodeName\tPP*\n");
        }
   
        
        ////////////////////////////////////////////////////////////////////
        // MAIN ALGORITHM LOOP BELOW !!!
        // DO KMERS ALIGNMENT AND SCORING FOR ALL SEQUENCE QUERIES
        ////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////

        int queryCounter=0;
        int queryPlacedCounter=0;
        for (int query = 0; query < samplingAmount; query++) { //<-- MAIN LOOP: QUERY PER QUERY, TODO PARALLELIZED VERSION

            queryCounter++;

            //console display to follow the process
            if ((queryCounter%200000)==0) {
                System.out.println(queryCounter+"/"+samplingAmount+" queries placed ("+(((0.0+queryCounter)/samplingAmount)*100)+"%)");
            }
            
            //debug
            if (queryCounter>queryLimit)
                break;
            
            
            ///////////////////////////////////
            //GENERATE RANDOM SEQ
            Fasta fasta=rs.generateSequence();

//            Infos.println("#######################################################################");
//            Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+fasta.getHeader());
//            Infos.println("#######################################################################");
            //fw.append(fasta.getFormatedFasta()+"\n");

            ///////////////////////////////////
            // PREPARE QUERY K-MERS
            QueryWord qw=null;
            SequenceKnife sk=new SequenceKnife(fasta, session.k, session.minK, session.states, queryWordSampling);

            
            ////////////////////////////////////////////////////////////////
            // BUILD THE ALIGNMENT AND SCORE IN A SIGNLE LOOP ON QUERY WORDS
            ////////////////////////////////////////////////////////////////
//            Infos.println("Launching scoring on candidate nodes...");
            int queryWordFoundCounter=0;
            int queryWordCounter=0;

            //loop on words
            while ((qw=sk.getNextWord())!=null) {
                queryWordCounter++;
                //Infos.println("Query mer: "+qw.toString());
                //position of this word
                int position=session.hash.getTopPosition(qw);
                //if this word is not registered in the hash
                if (position<0) {
                    continue;
                }
                queryWordFoundCounter++;
                //get Pairs associated to this word
                List<Pair> allPairs = session.hash.getPairsOfTopPosition(qw);
                //System.out.println("Pairs: "+allPairs);
                for (int i = 0; i < allPairs.size(); i++) {
                    Pair p = allPairs.get(i);
                    //we will score only encountered nodes, originalNode registered
                    //at 1st encouter
                    if (nodeOccurences[p.getNodeId()]==0) {
                        selectedNodes.add(p.getNodeId());
                    }
                    //count # times originalNode encountered
                    nodeOccurences[p.getNodeId()]+=1;
                    //score associated to originalNode x for current read
                    nodeScores[p.getNodeId()]+=p.getPPStar();
                }
                //System.out.println("  nodeOccurences:"+Arrays.toString(nodeOccurences));
                //System.out.println("  nodeScores:"+Arrays.toString(nodeScores));

            }

//            Infos.println("Proportion of query words retrieved in the hash: "+queryWordFoundCounter+"/"+queryWordCounter);
//            Infos.println("Candidate nodes: ("+selectedNodes.size()+") ");      

            //if selectedNodes is empty (no node could be associated)
            //for instance when no query words could be found in the hash
            if (selectedNodes.size()<1) {
                //Infos.println("Read cannot be placed.");
                //TODO: currently this query will not be in output csv 
                //and jplace... put them in special output ?
                continue; //to next query
            }
            
            queryPlacedCounter++;


            // NOW CORRECTING SCORING BY UNMATCHED WORDS
            ///////////////////////////////////////////////////////
            //now add the score corresponding to the words not found,
            // query.e. threshold*#_words_not_scored (because no in hash)
            int maxWords=sk.getMerCount();
            int bestNodeId=-1;
            float bestScore=Float.NEGATIVE_INFINITY;
            int secondBest=-1;
            float secondScore=Float.NEGATIVE_INFINITY;
            for (Integer nodeId:selectedNodes) {
                //System.out.println("Scoring originalNode:"+nodeId);
                //System.out.println("nodeMapping:"+session.nodeMapping.get(nodeId));
                //int extendedTreeId=session.nodeMapping.get(nodeId);
                //System.out.println("extendedTreeId:"+extendedTreeId);
                //Integer originalNodeId = extendedTree.getFakeToOriginalId(extendedTreeId);
                //System.out.println("originalNodeId:"+originalNodeId);
                //System.out.println("  scoring originalNode: ARTree="+session.ARTree.getById(nodeId)+" ExtendedTree="+session.extendedTree.getById(extendedTreeId)+" OriginalTree="+session.originalTree.getById(originalNodeId));
                nodeScores[nodeId]+=session.PPStarThresholdAsLog10*(maxWords-nodeOccurences[nodeId]);
                if (nodeScores[nodeId]>bestScore) {
                    secondBest=bestNodeId;
                    secondScore=bestScore;
                    bestNodeId=nodeId;
                    bestScore=nodeScores[nodeId];
                } else if (nodeScores[nodeId]>secondScore) {
                    secondBest=nodeId;
                    secondScore=nodeScores[nodeId];
                }
            }

//            Infos.println("Best node (ARTree) is : "+bestNodeId+" (score="+bestScore+")");
//            Infos.println("mapping: ARTree="+session.ARTree.getById(bestNodeId)+" ExtendedTree="+session.extendedTree.getById(session.nodeMapping.get(bestNodeId))+" OriginalTree="+session.originalTree.getById(session.extendedTree.getFakeToOriginalId(session.nodeMapping.get(bestNodeId))));

            //simple debug test
            if (session.nodeMapping.get(bestNodeId)==null) { //simple test
                System.out.println("bestNode not found: "+bestNodeId+" "+String.valueOf(session.ARTree.getById(bestNodeId).getLabel()));
                System.exit(1);
            }

            // SELECT BEST NEIGHBOOR FAKE NODE IF BEST NODE IS ORIGINAL NODE
            ////////////////////////////////////////////////////////////////
            //check if this was a fake originalNode or not
            //to do that, retromapping from ARTree to extended tree 
            int extendedTreeId=session.nodeMapping.get(bestNodeId);
            int originalNodeId = session.extendedTree.getFakeToOriginalId(extendedTreeId);
            PhyloNode nodeToTest = session.extendedTree.getById(extendedTreeId);
            //if this is an original originalNode, select adjacent branch 
            //leading to 2nd best PP*, if 2nd best PP* is original, do the same
            //for 3rd and so on...
            if (!nodeToTest.isFakeNode()) {
                //System.out.println("############### change best node to neighboors !");                    
//                Infos.println("Current best node is an original node...");
                PhyloNode firstNode = null;
                PhyloNode secondNode = null;
                //if there was no other scored nodes ? this case happened when a single query mer was found in the DB
                if (secondBest<0) {
                    //arbitrary choice, take FAKE node which is left son.
                    firstNode = session.ARTree.getById(bestNodeId);
                    secondNode = firstNode.getChildAt(0);
                    bestNodeId =secondNode.getId();

                } else {
                    firstNode = session.ARTree.getById(bestNodeId);
                    //select node of 2nd best score
                    secondNode = session.ARTree.getById(secondBest);
                    //get path from this node to bestNodeId
                    //System.out.println("1st node: "+ARTree.getById(bestNodeId)+" 2nd node:"+ARTree.getById(secondNodeId));
                    PhyloTree.Path shortestPath = session.ARTree.shortestPath(session.ARTree.getRoot(), firstNode, secondNode);
                    //System.out.println("Path to second: "+shortestPath.path);
                    //path will be as:
                    //firstNode-X0-...-secondNode
                    //or nodeToTest-secondNode(X0) if immediate neighboor
                    //in all case the 2nd elt of the path is the X0 chosen 
                    //for the placement
                    bestNodeId=shortestPath.path.get(1).getId();
                  
                }
                extendedTreeId=session.nodeMapping.get(bestNodeId);
                originalNodeId = session.extendedTree.getFakeToOriginalId(extendedTreeId);         
                //System.out.println("NEW Selected node (ARTree) is : "+bestNodeId+" (score="+bestScore+")");
                //System.out.println("mapping: ARTree="+ARTree.getById(bestNodeId)+" ExtendedTree="+extendedTree.getById(extendedTreeId)+" OriginalTree="+session.originalTree.getById(originalNodeId));

                if (!session.extendedTree.getById(extendedTreeId).isFakeNode()) {
                    System.out.println("Something went wrong in neighboor node search !!!!");
                    System.exit(1);
                }

            } 


            //basic normalization, divide score by number of words present
            //in the query
            float normalizedScore=bestScore/maxWords;
            normalizedScores.add(normalizedScore);
//            Infos.println("Normalized score: "+normalizedScore);




            //TO DEBUG
//                System.out.println(fasta.getHeader());
//                System.out.print("bestNodeId:"+bestNodeId);
//                System.out.print("\t\tbestScore:"+bestScore);
//                System.out.println("\t\tbestSCore (normalized by #mers):"+(bestScore/sk.getMaxMerCount()));
//                System.out.println("nodeScores:"+Arrays.toString(nodeScores));
//                System.out.println("nodeOccurences:"+Arrays.toString(nodeOccurences));
//                System.out.println("selectedNodes:"+selectedNodes.keySet().toString());


            //write result in file if a originalNode was hit
            if (bwTSV!=null) {
                //OUTPUT n°1: the CSV report of placement
                //allow in particular to check that nodes were correclty
                //mappes at every step (original tree, extended tree,
                //AR modified tree)
                sb.append(fasta.getHeader().split(" ")[0]).append("\t");
                sb.append(String.valueOf(bestNodeId)).append("\t"); //ARTree nodeID
                sb.append(String.valueOf(session.ARTree.getById(bestNodeId).getLabel())).append("\t"); //ARTree nodeName
                sb.append(String.valueOf(extendedTreeId)).append("\t"); //extended Tree nodeID
                sb.append(String.valueOf(session.extendedTree.getById(extendedTreeId).getLabel())).append("\t"); //extended Tree nodeName
                sb.append(String.valueOf(originalNodeId)).append("\t"); //edge of original tree (original nodeId)
                sb.append(String.valueOf(session.originalTree.getById(originalNodeId).getLabel())).append("\t"); //edge of original tree (original nodeName
                sb.append(String.valueOf(normalizedScore)).append("\n");
                
                //push the stringbuffer to the CSV bufferedwriter every 10000 sequences
                if ((queryCounter%10000)==0) {
                    int size=sb.length();
                    bwTSV.append(sb);
                    bwTSV.flush();
                    sb=null;
                    sb=new StringBuffer(size);
                }
            }
            //reset the scoring vectors
            for(Integer nodeId:selectedNodes) {
                nodeScores[nodeId]=0.0f;
                nodeOccurences[nodeId]=0;
            }
            selectedNodes.clear();

        }
        
        //flush to disk the last CSV buffer
        if (bwTSV!=null) {
            bwTSV.append(sb);
        }
        System.out.println("Queries actually placed:"+queryPlacedCounter);

        //use the normalized scores to define quantile
        return new Double(Quantiles.scale(q_quantile).index(n_quantile).compute(normalizedScores)).floatValue();        
        
    }
    
    
    
    
    
    
    
    
    
    /**
     * 
     * @param fp
     * @param placements
     * @param bwTSV
     * @param queryWordSampling
     * @param minOverlap
     * @return the number of queries effectively placed (kmers were found in DB)
     * @throws java.io.IOException
     */
    public int processQueries(SequencePointer fp, JSONArray placements, BufferedWriter bwTSV, int queryWordSampling, int minOverlap) throws IOException {
        
        ///////////////////////////////////////////////////////            
        // PREPARE VECTORS USED TO ALIGN AND SCORE NODES

        //list of nodes encountered during matches search in the hash
        //variable size, expanded only at 1st encounter with the node
        //when nodeOccurences[nodeId]==0
        //,reset at each read
        ArrayList<Integer> selectedNodes=new ArrayList<>();
        //instanciated once
        int[] nodeOccurences=new int[session.ARTree.getNodeCount()]; // tab[#times_encoutered] --> index=nodeId
        float[] nodeScores=new float[session.ARTree.getNodeCount()]; // tab[score] --> index=nodeId
        //Arrays.fill(nodeOccurences, 0);
        //Arrays.fill(nodeScores, 0.0f);     


        ///////////////////////////////////////////////////////////////////       
        // PREPARE CHECKSUM FOR IDENTICAL READS REGISTER
        AbstractChecksum checksumGenerator=null;
        try {
            checksumGenerator = JacksumAPI.getChecksumInstance("sha256");
            checksumGenerator.setEncoding(AbstractChecksum.HEX);
        } catch (NoSuchAlgorithmException ex) {
            Logger.getLogger(Main_PLACEMENT_v07.class.getName()).log(Level.SEVERE, null, ex);
        }
        //map to associate seqeunce checksums to Fasta headers
        //map(checksum)=List of headers corresponding to identical sequences
        HashMap<String,ArrayList<String>> identicalSeqsRegistry=new HashMap<>();
        //map to associate sequence checksums to JSONObject
        //map(checksum)=JSON placement object to which identical reads
        //              are associated (same score and placement)
        HashMap<String,JSONObject> checksumToJSONObject=new HashMap<>();

        //////////////////////////////////////////////////////////////////
        // PREPARE TSV OUTPUT
        //header of CSV output
        if (bwTSV!=null) {
            sb=new StringBuffer("Query\tARTree_NodeId\tARTree_NodeName\tExtendedTree_NodeId\tExtendedTree_NodeName\tOriginal_NodeId\tOriginal_NodeName\tPP*\n");
        }
   
        
        


        ////////////////////////////////////////////////////////////////////
        // MAIN ALGORITHM LOOP BELOW !!!
        // DO KMERS ALIGNMENT AND SCORING FOR ALL SEQUENCE QUERIES
        ////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////
        long startTotalPlacementTime=System.currentTimeMillis();


        int queryCounter=0;
        long totalAlignTime=0;
        long totalScoringTime=0;
        long totalWritingTime=0;
        long totalResetTime=0;
        long totalChecksumTime=0;

        int totalQueries=fp.getContentSize();
        Fasta fasta=null;
        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {  //<-- MAIN LOOP: QUERY PER QUERY, TODO PARALLELIZED VERSION

            queryCounter++;

            //console display to follow the process
            if ((queryCounter%200000)==0) {
                System.out.println(queryCounter+"/"+totalQueries+" queries placed ("+(((0.0+queryCounter)/totalQueries)*100)+"%)");
            }
            
            //debug
            if (queryCounter>queryLimit)
                break;
            
            ///////////////////////////////////
            // CHECKSUM BUILD
            //if already present do not compute placement
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
            //if this query sequence was already encountered
            if (identicalSeqsRegistry.containsKey(checksum)) {
                identicalSeqsRegistry.get(checksum).add(subHeader);
                Infos.println("! SKIPPED BECAUSE DUPLICATE: "+fasta.getHeader());
                queryCounter++;
                //before skipping, update the outputs
                //for now, only jplace is done here...
                //note that detailed comments about the jplace are
                //in the jplace block on the bottom of the loop
                //get back the placement out from the JSONObject
                //if it passed the --nsbound debug option (if not, do not exists)
                JSONObject placement = checksumToJSONObject.get(checksum);
                if (placement!=null) {
                    //the "p" object values are the same, no changes
                    //JSONArray pMetadata=(JSONArray)placement.get("p");
                    //the "nm" object has to be extended with the identifier
                    //and multiplicity of this read
                    JSONArray allIdentifiers=(JSONArray)placement.get("nm");
                    JSONArray readMultiplicity=new JSONArray();
                    readMultiplicity.add(fasta.getHeader());
                    readMultiplicity.add(1);
                    allIdentifiers.add(readMultiplicity);   
                }
                //go to next query
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
            SequenceKnife sk=new SequenceKnife(fasta, session.k, session.minK, session.states, queryWordSampling);
            int queryWordCounter=0;
            int queryWordFoundCounter=0;




            //if alignment graph are required
            DiagSum bestPPStarsDiagsum=null;
            double[][] graphDataForTopTuples =null;
            int xSize=-1;
            if (graphAlignment) {
                //a diagsum based only on a the single top best PP* associated to a query word.
                //accumulating all nodes, just taing higher PP* and corresponding positions
                bestPPStarsDiagsum=new DiagSum(queryLength, session.align.getLength(), minOverlap, session.k, sk.getStep());
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
                            graphDataForTopTuples[2][x+(y*xSize)]=session.PPStarThresholdAsLog10;
                        }
                    }
            }

            ////////////////////////////////////////////////////////////////
            // BUILD THE ALIGNMENT AND SCORE IN A SINGLE LOOP ON QUERY WORDS
            ////////////////////////////////////////////////////////////////
            Infos.println("Launching scoring on candidate nodes...");
            //loop on words
            while ((qw=sk.getNextWord())!=null) {
                //Infos.println("Query mer: "+qw.toString());
                queryWordCounter++;
                //position of this word
                int topPosition=session.hash.getTopPosition(qw);
                //if this word is not registered in the hash
                if (topPosition<0) {
                    continue;
                }
                queryWordFoundCounter++;
                //get Pairs associated to this word
                List<Pair> allPairs = session.hash.getPairsOfTopPosition(qw);
                //System.out.println("Pairs: "+allPairs);
                for (int i = 0; i < allPairs.size(); i++) {
                    Pair p = allPairs.get(i);
                    //we will score only encountered nodes, originalNode registered
                    //at 1st encouter
                    if (nodeOccurences[p.getNodeId()]==0) {
                        selectedNodes.add(p.getNodeId());
                    }
                    //count # times originalNode encountered
                    nodeOccurences[p.getNodeId()]+=1;
                    //score associated to originalNode x for current read
                    nodeScores[p.getNodeId()]+=p.getPPStar();
                }
                //System.out.println("  nodeOccurences:"+Arrays.toString(nodeOccurences));
                //System.out.println("  nodeScores:"+Arrays.toString(nodeScores));

                //graph of alignment, if asked
                if(graphAlignment) {
                    int topDiagSumPos=topPosition-qw.getOriginalPosition()+(queryLength-minOverlap);
                    if (topDiagSumPos>-1 && topDiagSumPos<bestPPStarsDiagsum.getSize()) {
                            graphDataForTopTuples[2][queryWordCounter*xSize+topPosition]= allPairs.get(0).getPPStar();                        
                    }
                }

                //DEBUG
                //if (queryWordCounter>1000)
                //        break;
                //DEBUG

            }

            //display alignment result of requested
            if (graphAlignment) {
                //graph for topTuplePerPos
                DefaultXYZDataset datasetForGraph=new DefaultXYZDataset();
                //datasetForGraph.addSeries(0, graphDataForTopTuples);
                datasetForGraph.addSeries(0, graphDataForTopTuples);
                JFrame infos=new JFrame();
                infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                infos.setLayout(new GridLayout(1, 1));
                infos.add(ChartsForNodes.buildReadMatchForANode4("Top_tuple_per_pos of read="+fasta.getHeader(),bestPPStarsDiagsum.getSize(),queryLength, datasetForGraph,session.PPStarThresholdAsLog10,0));
                infos.setSize(1024, 250);
                infos.pack();
                RefineryUtilities.centerFrameOnScreen(infos);
                infos.setVisible(true);
            }


            Infos.println("Proportion of query words retrieved in the hash: "+queryWordFoundCounter+"/"+queryWordCounter);
            long endAlignTime=System.currentTimeMillis();
            totalAlignTime+=(endAlignTime-startAlignTime);
            Infos.println("Candidate nodes: ("+selectedNodes.size()+") ");      

            //if selectedNodes is empty (no node could be associated)
            //for instance when no query words could be found in the hash
            if (selectedNodes.size()<1) {
                //Infos.println("Read cannot be placed.");
                //TODO: currently this query will not be in output csv 
                //and jplace... put them in special output ?
                continue; //to next query
            }


            // NOW CORRECTING SCORING BY UNMATCHED WORDS
            ///////////////////////////////////////////////////////
            long startScoringTime=System.currentTimeMillis();
            //now add the score corresponding to the words not found,
            // query.e. threshold*#_words_not_scored (because no in hash)
            int maxWords=sk.getMerCount();
            int bestNodeId=-1;
            float bestScore=Float.NEGATIVE_INFINITY;
            int secondBest=-1;
            float secondScore=Float.NEGATIVE_INFINITY;
            for (Integer nodeId:selectedNodes) {
                //System.out.println("Scoring originalNode:"+nodeId);
                //System.out.println("nodeMapping:"+session.nodeMapping.get(nodeId));
                //int extendedTreeId=session.nodeMapping.get(nodeId);
                //System.out.println("extendedTreeId:"+extendedTreeId);
                //Integer originalNodeId = extendedTree.getFakeToOriginalId(extendedTreeId);
                //System.out.println("originalNodeId:"+originalNodeId);
                //System.out.println("  scoring originalNode: ARTree="+session.ARTree.getById(nodeId)+" ExtendedTree="+session.extendedTree.getById(extendedTreeId)+" OriginalTree="+session.originalTree.getById(originalNodeId));
                nodeScores[nodeId]+=session.PPStarThresholdAsLog10*(maxWords-nodeOccurences[nodeId]);
                if (nodeScores[nodeId]>bestScore) {
                    secondBest=bestNodeId;
                    secondScore=bestScore;
                    bestNodeId=nodeId;
                    bestScore=nodeScores[nodeId];
                } else if (nodeScores[nodeId]>secondScore) {
                    secondBest=nodeId;
                    secondScore=nodeScores[nodeId];
                }
            }

            Infos.println("Best node (ARTree) is : "+bestNodeId+" (score="+bestScore+")");
            Infos.println("mapping: ARTree="+session.ARTree.getById(bestNodeId)+" ExtendedTree="+session.extendedTree.getById(session.nodeMapping.get(bestNodeId))+" OriginalTree="+session.originalTree.getById(session.extendedTree.getFakeToOriginalId(session.nodeMapping.get(bestNodeId))));

            //simple debug test
            if (session.nodeMapping.get(bestNodeId)==null) { //simple test
                System.out.println("bestNode not found: "+bestNodeId+" "+String.valueOf(session.ARTree.getById(bestNodeId).getLabel()));
                System.exit(1);
            }

            // SELECT BEST NEIGHBOOR FAKE NODE IF BEST NODE IS ORIGINAL NODE
            ////////////////////////////////////////////////////////////////
            //check if this was a fake originalNode or not
            //to do that, retromapping from ARTree to extended tree 
            int extendedTreeId=session.nodeMapping.get(bestNodeId);
            int originalNodeId = session.extendedTree.getFakeToOriginalId(extendedTreeId);
            PhyloNode nodeToTest = session.extendedTree.getById(extendedTreeId);
            //if this is an original originalNode, select adjacent branch 
            //leading to 2nd best PP*
            if (!nodeToTest.isFakeNode()) {
                //System.out.println("############### change best node to neighboors !");                    
                Infos.println("Current best node is an original node...");

                PhyloNode firstNode = null;
                PhyloNode secondNode = null;
                //if there was no other scored nodes ? this case happened when a single query mer was found in the DB
                if (secondBest<0) {
                    //arbitrary choice, take FAKE node which is left son.
                    firstNode = session.ARTree.getById(bestNodeId);
                    secondNode = firstNode.getChildAt(0);
                    bestNodeId =secondNode.getId();

                } else {
                    firstNode = session.ARTree.getById(bestNodeId);
                    //select node of 2nd best score
                    secondNode = session.ARTree.getById(secondBest);
                    //get path from this node to bestNodeId
                    //System.out.println("1st node: "+ARTree.getById(bestNodeId)+" 2nd node:"+ARTree.getById(secondNodeId));
                    PhyloTree.Path shortestPath = session.ARTree.shortestPath(session.ARTree.getRoot(), firstNode, secondNode);
                    //System.out.println("Path to second: "+shortestPath.path);
                    //path will be as:
                    //firstNode-X0-...-secondNode
                    //or nodeToTest-secondNode(X0) if immediate neighboor
                    //in all case the 2nd elt of the path is the X0 chosen 
                    //for the placement
                    bestNodeId=shortestPath.path.get(1).getId();
                  
                }
                extendedTreeId=session.nodeMapping.get(bestNodeId);
                originalNodeId = session.extendedTree.getFakeToOriginalId(extendedTreeId);        
                
                
                //System.out.println("NEW Selected node (ARTree) is : "+bestNodeId+" (score="+bestScore+")");
                //System.out.println("mapping: ARTree="+ARTree.getById(bestNodeId)+" ExtendedTree="+extendedTree.getById(extendedTreeId)+" OriginalTree="+session.originalTree.getById(originalNodeId));

                if (!session.extendedTree.getById(extendedTreeId).isFakeNode()) {
                    System.out.println("Something went wrong in neighboor node search !!!!");
                    System.exit(1);
                }

            } 


            //basic normalization, divide score by number of words present
            //in the query
            float normalizedScore=bestScore/maxWords;
            Infos.println("Normalized score: "+normalizedScore);


            long endScoringTime=System.currentTimeMillis();
            totalScoringTime+=endScoringTime-startScoringTime;


            //TO DEBUG
//                System.out.println(fasta.getHeader());
//                System.out.print("bestNodeId:"+bestNodeId);
//                System.out.print("\t\tbestScore:"+bestScore);
//                System.out.println("\t\tbestSCore (normalized by #mers):"+(bestScore/sk.getMaxMerCount()));
//                System.out.println("nodeScores:"+Arrays.toString(nodeScores));
//                System.out.println("nodeOccurences:"+Arrays.toString(nodeOccurences));
//                System.out.println("selectedNodes:"+selectedNodes.keySet().toString());

            long startWritingTime=System.currentTimeMillis();
            //write result in file if a originalNode was hit
            if (bwTSV!=null) {
                //only score passing --nsbound if this debug option is set
                if (normalizedScore>=nsBound) {
                    //OUTPUT n°1: the CSV report of placement
                    //allow in particular to check that nodes were correclty
                    //mappes at every step (original tree, extended tree,
                    //AR modified tree)
                    sb.append(fasta.getHeader().split(" ")[0]).append("\t");
                    sb.append(String.valueOf(bestNodeId)).append("\t"); //ARTree nodeID
                    sb.append(String.valueOf(session.ARTree.getById(bestNodeId).getLabel())).append("\t"); //ARTree nodeName
                    sb.append(String.valueOf(extendedTreeId)).append("\t"); //extended Tree nodeID
                    sb.append(String.valueOf(session.extendedTree.getById(extendedTreeId).getLabel())).append("\t"); //extended Tree nodeName
                    sb.append(String.valueOf(originalNodeId)).append("\t"); //edge of original tree (original nodeId)
                    sb.append(String.valueOf(session.originalTree.getById(originalNodeId).getLabel())).append("\t"); //edge of original tree (original nodeName
                    sb.append(String.valueOf(normalizedScore)).append("\n");
                }
            }
            //OUTPUT n°2: the JSON placement object (jplace file)
            //2 possibilities:
            //-either do a new placement object and add ot to the list
            // of placements (block below is executed)
            //-or a previous placement object corresponding to an
            //identical sequence exists, then we don't the block below
            //as all the alignmnet/placement algo was skipped.

            //only score passing --nsbound if this debug option is set
            if (normalizedScore>=nsBound) {
                //System.out.println("Score pass threshold: normalizedScore="+normalizedScore+" nsBound="+nsBound);
                //object representing placements (mandatory), it contains 2 key/value couples
                JSONObject placement=new JSONObject();
                //first we build the "p" array, containing the position/scores of all reads
                JSONArray pMetadata=new JSONArray();
                //in pplacer/EPA several placements can be associated to a query
                //we input only the best one, but that can be changed in the future
                JSONArray placeColumns=new JSONArray();

                //fake fields for compatibility with current tools (guppy, archeopteryx)
                //should be provided as an option
                placeColumns.add(0.1); //distal_length
                placeColumns.add(normalizedScore); //we put also PP* in the like_weight_ratio column to allow sorting in iTol
                placeColumns.add(0.1); //pendant_length

                //placeColumns.add(session.ARTree.getById(bestNodeId).getLabel()); // 1. ARTree nodeName
                //placeColumns.add(session.extendedTree.getById(extendedTreeId).getLabel()); // 2. extended tree nodeName
                placeColumns.add(originalNodeId); // 3. edge of original tree (original nodeId=edgeID)
                //placeColumns.add(session.originalTree.getById(originalNodeId).getLabel()); // 4. edge of original tree (original nodeName)
                placeColumns.add(normalizedScore); // 4. PP*
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

            long endWritingTime=System.currentTimeMillis();
            totalWritingTime+=endWritingTime-startWritingTime;


            //push the stringbuffer to the CSV bufferedwriter every 50000 sequences
            if (bwTSV!=null) {
                if ((queryCounter%50000)==0) {
                    int size=sb.length();
                    bwTSV.append(sb);
                    sb=null;
                    sb=new StringBuffer(size);
                }   
            }

            //reset the scoring vectors
            long startResetTime=System.currentTimeMillis();
            for (Integer nodeId : selectedNodes) {
                nodeScores[nodeId]=0.0f;
                nodeOccurences[nodeId]=0;
            }
            selectedNodes.clear();
            long endResetTime=System.currentTimeMillis();
            totalResetTime+=endResetTime-startResetTime;

        }
        
        //flush to disk the last CSV buffer
        bwTSV.append(sb);
            
        long endTotalPlacementTime=System.currentTimeMillis();
        System.out.println("############################################################");
        System.out.println("Checksum registry took in total: "+totalChecksumTime+" ms");
        System.out.println("Alignments and pre-scoring took in total: "+totalAlignTime+" ms");
        System.out.println("Scoring took in total: "+totalScoringTime+" ms");
        System.out.println("Writing .csv and .jplace took in total: "+totalWritingTime+" ms");
        System.out.println("Reset took in total: "+totalResetTime+" ms");
        System.out.println("------------------------------------------------------------");
        System.out.println("(per read average)");
        System.out.println("Checksum registry took on average: "+((0.0+totalChecksumTime)/queryCounter)+" ms");
        System.out.println("Alignments and pre-scoring took on average: "+((0.0+totalAlignTime)/queryCounter)+" ms");
        System.out.println("Scoring took on average: "+((0.0+totalScoringTime)/queryCounter)+" ms");
        System.out.println("Writing .csv and .jplace took on average: "+((0.0+totalWritingTime)/queryCounter)+" ms");
        System.out.println("Reset took on average: "+((0.0+totalResetTime)/queryCounter)+" ms");
        System.out.println("------------------------------------------------------------");
        System.out.println("Placement Process (without DB load) took: "+(endTotalPlacementTime-startTotalPlacementTime)+" ms");
        System.out.println("############################################################");
        
        return queryCounter;
    }
    
    
    public boolean parallelProcessReads() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    
    
}
