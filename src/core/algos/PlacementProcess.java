/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import core.DNAStatesShifted;
import core.hash.Pair;
import etc.Infos;
import inputs.Fasta;
import inputs.SequencePointer;
import it.unimi.dsi.fastutil.Hash;
import it.unimi.dsi.fastutil.chars.Char2FloatMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenCustomHashMap;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.security.NoSuchAlgorithmException;
import java.text.NumberFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import jonelo.jacksum.JacksumAPI;
import jonelo.jacksum.algorithm.AbstractChecksum;
import main_v2.Main_PLACEMENT_v07;
import main_v2.SessionNext_v2;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import tree.PhyloNode;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class PlacementProcess {
    
    
    //debug/////////////////////////////////////////////////////////////
    boolean debug=false;
    //csv log
    boolean csvLog=false;
    //max number of queries treated 
    int queryLimit=Integer.MAX_VALUE;
    //graph of words alignment
    boolean graphAlignment=false; //NOTE: This will work only if hash is based on PositionNodes
    boolean merStats=false; //log outputing stats associated with mers, CAUTION produces big files
    //test different way to select best score
    boolean useSelectionAlgo=true; // score using Hoare's selection algorithm
    //debug/////////////////////////////////////////////////////////////
    
    
    //parameters asked when program launched
    SessionNext_v2 session=null;
    Float nsBound=Float.NEGATIVE_INFINITY;
    
    //elements related to fil outputs
    //filled while the placement process is running
    JSONArray placements=null;
    BufferedWriter bwTSV=null;
    StringBuffer sb=null;
    NumberFormat nf = NumberFormat.getInstance(Locale.UK);
    
    /**
     *
     * @param session
     * @param nsBound
     * @param queryLimit the value of queryLimit
     */
    public PlacementProcess(SessionNext_v2 session, Float nsBound, int queryLimit) {
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

        //list of bestNormalizedScore
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
        
        //init kmer extraction
        SequenceKnife sk=new SequenceKnife(session.k, session.minK, session.states, queryWordSampling);

        int queryCounter=0;
        int queryPlacedCounter=0;
        for (int query = 0; query < samplingAmount; query++) { //<-- MAIN LOOP: QUERY PER QUERY, TODO PARALLELIZED VERSION

            queryCounter++;

            //console display to follow the process
            if ((queryCounter%200000)==0) {
                System.out.println(queryCounter+"/"+samplingAmount+" queries placed ("+(((0.0+queryCounter)/samplingAmount)*100)+"%)");
            }
            
            //debug
            if (debug && (queryCounter>queryLimit) )
                break;
            
            ///////////////////////////////////
            //GENERATE RANDOM SEQ
            Fasta fasta=rs.generateSequence();

            ///////////////////////////////////
            // PREPARE QUERY K-MERS
            byte[] qw=null;
            sk.init(fasta);
            
            ////////////////////////////////////////////////////////////////
            // BUILD THE ALIGNMENT AND SCORE IN A SIGNLE LOOP ON QUERY WORDS
            ////////////////////////////////////////////////////////////////
            int queryWordFoundCounter=0;
            int queryWordCounter=0;

            //loop on words
            while ((qw=sk.getNextByteWord())!=null) {
                queryWordCounter++;
                //Infos.println("Query mer: "+qw.toString());
                
                //get Pairs associated to this word
                List<Pair> allPairs = session.hash.getPairsOfTopPosition(qw);
                //if this word is not registered in the hash
                if (allPairs==null) {
                    continue;
                }
                queryWordFoundCounter++;
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

//            Infos.println("Proportion of query words retrieved in the hash: "+queryKmerMatchingDB+"/"+queryKmerCount);
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
//            Infos.println("Normalized score: "+bestNormalizedScore);




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
        //return new Double(Quantiles.scale(q_quantile).index(n_quantile).compute(normalizedScores)).floatValue();

        throw new UnsupportedOperationException();
        
    }
    
    
    public static void merStatsHeader(BufferedWriter bw) throws IOException {
    	bw.append("Query\tARNodeId\tARNodeName\tExtendedNodeId\tExtendedNodeName\tOriginalNodeId\tOriginalNodeName\tMerPos\tPPStar\n");
    }
    
    public static void merStats(SessionNext_v2 session, ArrayList<float[]> merFound, Fasta fasta, BufferedWriter bw) throws IOException {
        for (int nodeId = 0; nodeId < merFound.size(); nodeId++) {
            for (int merPos = 0; merPos < merFound.get(nodeId).length; merPos++) {
                if (!session.ARTree.getById(nodeId).isLeaf()) {
                    int extendedTreeId=session.nodeMapping.get(nodeId);
                    int originalNodeId = session.extendedTree.getFakeToOriginalId(extendedTreeId);
                    PhyloNode extNode = session.extendedTree.getById(extendedTreeId);
                    PhyloNode origNode = session.originalTree.getById(originalNodeId);
                    bw.append(fasta.getHeader()+"\t");
                    bw.append(nodeId+"\t"+session.ARTree.getById(nodeId).getLabel()+"\t");
                    bw.append(extendedTreeId+"\t"+extNode.getLabel()+"\t");
                    bw.append(originalNodeId+"\t"+origNode.getLabel()+"\t");
                    bw.append(merPos+"\t");
                    bw.append(String.valueOf(merFound.get(nodeId)[merPos]));
                    bw.append("\n");
                }
            }
        }
    }
    

    public static float computeWeightRatioShift(float lowest, float best) {
    	float weightRatioShift = 0.0f;
    	if ( -308f >= lowest ) { // this is Double.MIN_NORMAL=2.2250738585072014E-308 as reported by javadoc
            weightRatioShift=best;
    	}
    	return weightRatioShift;
    }
    
    public static double computeWeightRatio(Score s, double weightRatioShift, double allLikelihoodSums) {
    	return Math.pow(10.0, (double)(s.score-weightRatioShift))/allLikelihoodSums;
    }
    
    public double fillBestScoreList(float[] nodeScores, List<Integer> selectedNodes, Score[] bestScoreList, int numBest) {
    	double allLikelihoodSums = 0.0;

        // Support k-best max heap, where k is the number of placements to keep.
        // Visit only 'selected nodes'
        PriorityQueue<Score> heap = new PriorityQueue<Score>();
        for (int nodeId : selectedNodes)
        {
            heap.add(new Score(nodeId, nodeScores[nodeId]));

            if (heap.size() > numBest) {
                heap.remove();
            }
        }

        // Find the max and the min score among selected,
        // find the total score sum
        float lowest = 0.0f;
        float best = -Float.MAX_VALUE;
        for (Score scorePair: heap) {
            // calculate the total likelihood sum (for the likelihood weight ratio)
            // before the power of ten, divide by #words, to use the normalized score
            allLikelihoodSums += Math.pow(10.0, scorePair.score);

            // remind lower power of 10, if goes below double limit,
            // will need to shift the values
            if (scorePair.score < lowest) {
                lowest = scorePair.score;
            }
            if (scorePair.score > best) {
                best = scorePair.score;
            }
        }

        // copy to the output array
        int i = 0;
        for (Score scorePair : heap) {
            bestScoreList[i] = scorePair;
            i++;
        }
        Arrays.sort(bestScoreList);

        // if necessary
        // prepare variable weightRatioShift for allowing
        // weight-ratio computations when proba is < to "double" primitive boundary
        float weightRatioShift = computeWeightRatioShift(lowest, best);
        if (weightRatioShift != 0.0f) {
            // rebuild allLikelihoodSums using the shift
            allLikelihoodSums = 0.0;
            for (int ii = bestScoreList.length-numBest; ii < bestScoreList.length; ii++) {
                allLikelihoodSums += Math.pow(10.0, (double)(bestScoreList[ii].score-weightRatioShift));
            }
        }
        
        return allLikelihoodSums;
    }

    
    /**
     * 
     * @param fp
     * @param placements
     * @param bwTSV
     * @param bwNotPLaced
     * @param sk
     * @param minOverlap
     * @param logDir
     * @param keepAtMost
     * @param guppyCompatible
     * @param keepFactor
     * @param treatAmbiguities
     * @param treatAmbiguitiesWithMax
     * @return the number of queries effectively placed (kmers were found in DB)
     * @throws java.io.IOException
     */
    public int processQueries(  SequencePointer fp,
                                JSONArray placements,
                                BufferedWriter bwTSV,
                                BufferedWriter bwNotPLaced,
                                ISequenceKnife sk,
                                int minOverlap,
                                File logDir,
                                int keepAtMost,
                                float keepFactor,
                                boolean guppyCompatible,
                                boolean treatAmbiguities,
                                boolean treatAmbiguitiesWithMax
                            ) throws IOException {
        
        
        ///////////////////////////////////////////////////////            
        // PREPARE VECTORS USED TO ALIGN AND SCORE NODES

        //list of nodes encountered during matches search in the hash
        //variable size, expanded only at 1st encounter with the node
        //when nodeOccurences[nodeId]==0
        //,reset at each read
        ArrayList<Integer> L=new ArrayList<>(10);
        //instanciated once
        int[] C=new int[session.originalTree.getNodeCount()]; // tab[#times_encoutered] --> index=nodeId
        float[] S=new float[session.originalTree.getNodeCount()]; // tab[score] --> index=nodeId
        int numberOfBestScoreToConsiderForOutput=-1;
        Score[] bestScoreList=new Score[keepAtMost]; //in ascending order
        for (int i = 0; i < bestScoreList.length; i++) {
            bestScoreList[i]=new Score(-1, Float.NEGATIVE_INFINITY);
        }

        ///////////////////////////////////////////////////////////////////       
        // PREPARE CHECKSUM FOR IDENTICAL READS REGISTER
        AbstractChecksum checksumGenerator=null;
        try {
            //128bits checksum should be enough to avoid collision over
            //millions of queries (with 32bits, 4 collisions over 10^6 reads)
            checksumGenerator = JacksumAPI.getChecksumInstance("md5"); 
            checksumGenerator.setEncoding(AbstractChecksum.BIN);
        } catch (NoSuchAlgorithmException ex) {
            Logger.getLogger(Main_PLACEMENT_v07.class.getName()).log(Level.SEVERE, null, ex);
        }
        //map to associate sequence checksums to JSONObject
        //map(checksum)=JSON placement object to which identical reads
        //              are associated (same score and placement)
        Object2ObjectOpenCustomHashMap<byte[],JSONObject> checksumToJSONObject=new Object2ObjectOpenCustomHashMap(new Hash.Strategy<byte[]>() {
            @Override
            public int hashCode(byte[] o) {
                return Arrays.hashCode(o);
            }

            @Override
            public boolean equals(byte[] a, byte[] b) {
                return Arrays.equals(a, b);
            }
        });
        

        //////////////////////////////////////////////////////////////////
        // PREPARE TSV OUTPUT
        //header of CSV output
        if (bwTSV!=null) {
            sb=new StringBuffer("Query\tARTree_NodeId\tARTree_NodeName\tExtendedTree_NodeId\tExtendedTree_NodeName\tOriginal_NodeId\tOriginal_NodeName\tPP*\n");
        }
        //debug file of mers stats
        BufferedWriter bwMerStats =null;
        if (merStats) {
            File fileMerStats=new File(logDir+File.separator+"merStats.csv");
            bwMerStats = Files.newBufferedWriter(fileMerStats.toPath());
            merStatsHeader(bwMerStats);
        }
        


        ////////////////////////////////////////////////////////////////////
        // MAIN ALGORITHM LOOP BELOW !!!
        // DO KMERS ALIGNMENT AND SCORING FOR ALL SEQUENCE QUERIES
        ////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////
        long startTotalPlacementTime=System.currentTimeMillis();
        System.out.println("Starting to place queries...");

        int queryCounter=0;
        int placedQueryCounter=0;
        int selectedNodeSize=0;
//        long totalKnifeTime=0;
//        long totalAlignTime=0;
//        long totalScoringTime=0;
//        long totalWritingTime=0;
//        long totalResetTime=0;
//        long totalChecksumTime=0;
//        long totalT1Time=0;
//        long totalT2Time=0;
        
        int totalQueries=fp.getContentSize();
        Fasta fasta=null;
        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {  //<-- MAIN LOOP: QUERY PER QUERY, TODO PARALLELIZED VERSION

            queryCounter++;
            
            //console display to follow the process
            if ((queryCounter%10000)==0) {
                System.out.println(queryCounter+"/"+totalQueries+
                " queries placed ("+
                (((0.0+queryCounter)/totalQueries)*100)+
                "%)  --  Time elapsed: "+
                ((0.0+(System.currentTimeMillis()-startTotalPlacementTime))/1000)+" s");
            }
            
            //debug
            if (debug && queryCounter>queryLimit)
                break;
            
            ///////////////////////////////////
            // CHECKSUM BUILD
            //if already present do not compute placement
            //again. this might need compression if fasta sequences headers
            //are too heavy in memory
//            long startChecksumTime=System.currentTimeMillis();
            checksumGenerator.reset(); //make it ready before next checksum computation
            checksumGenerator.update(fasta.getSequence(true).getBytes());
            byte[] checksum = checksumGenerator.getByteArray();
            int cutIndex=fasta.getHeader().indexOf(" ");
            if (cutIndex<0) { //basically, space not found
                cutIndex=fasta.getHeader().length();
            }
            String subHeader=fasta.getHeader().substring(0,cutIndex);
            //if this query sequence was already encountered
            JSONObject placement =null;
            if ((placement=checksumToJSONObject.get(checksum))!=null) {
                //ArrayList<String> array = identicalSeqsRegistry.get(checksum);
                //array.add(subHeader);
                Infos.println("! SKIPPED BECAUSE DUPLICATE: "+fasta.getHeader());
                //+" MULTIPLICITY:"+identicalSeqsRegistry.get(checksum));
                
                //System.out.println(array);
                //before skipping, update the outputs
                //in the jplace block on the bottom of the main placement loop
                //get back the placement out from the JSONObject
                //if it passed the --nsbound debug option (if not, do not exists)
                if (placement!=null) {
                    //the "p" object values stay unchanged
                    //the "nm" object has to be extended with the identifier
                    //and multiplicity of this read
                    JSONArray allQueryIdentifiers=(JSONArray)placement.get("nm");
                    JSONArray queryMultiplicity=new JSONArray();
                    queryMultiplicity.add(subHeader);
                    queryMultiplicity.add(1);
                    allQueryIdentifiers.add(queryMultiplicity);   
                }
                //go to next query, as detected as duplicate and jplace file now updated
                continue;
            } else {
                placement=new JSONObject();
                //ArrayList<String> a=new ArrayList<>();
                //a.add(subHeader);
                //identicalSeqsRegistry.put(checksum,1);
            }
//            long endChecksumTime=System.currentTimeMillis();
//            totalChecksumTime+=endChecksumTime-startChecksumTime;

            placedQueryCounter++;

            Infos.println("#######################################################################");
            Infos.println("### PLACEMENT FOR QUERY #"+queryCounter+" : "+fasta.getHeader());
            Infos.println("#######################################################################");
            //fw.append(fasta.getFormatedFasta()+"\n");
            int queryLength=fasta.getSequence(false).length();
            Infos.println("Query length: "+queryLength+ " (Q="+(queryLength-session.k+1)+")");

            ///////////////////////////////////
            // PREPARE QUERY K-MERS
//            long startKnifeTime=System.currentTimeMillis();
            sk.init(fasta);
            int queryKmerCount=0;
            int ambiguousMerTreated=0;
            int queryKmerMatchingDB=0;
            int skippedKmer=0;
//            long endKnifeTime=System.currentTimeMillis();
//            totalKnifeTime+=(endKnifeTime-startKnifeTime);



            //if alignment graph are required
//            double[][] graphDataForTopTuples =null;
//            int xSize=-1;
//            if (graphAlignment) {
//                //preparing preplacement graph data :
//                xSize=session.align.getLength();
//                int ySize=queryLength;
//                graphDataForTopTuples =new double[3][xSize*ySize];
//                //init the Z axis (PP*) to very small values for all possible (X,Y)
//                    for (int y = 0; y < ySize; y++) {
//                        for (int x = 0; x < xSize; x++) {
//                            //System.out.println(col+" "+line+" "+(col+(line*session.align.getLength())));
//                            graphDataForTopTuples[0][x+(y*xSize)]=x;
//                            graphDataForTopTuples[1][x+(y*xSize)]=y;
//                            graphDataForTopTuples[2][x+(y*xSize)]=session.PPStarThresholdAsLog10;
//                        }
//                    }
//            }

            ////////////////////////////////////////////////////////////////
            // BUILD THE ALIGNMENT AND SCORE IN A SINGLE LOOP ON QUERY WORDS
            ////////////////////////////////////////////////////////////////
//            long startAlignTime=System.currentTimeMillis();
            ArrayList<float[]> merDetails=new ArrayList<>(session.ARTree.getNodeCount());
            for (int i = 0; i < session.ARTree.getNodeCount(); i++) {
                float[] t=new float[sk.getMerCount()];
                Arrays.fill(t, session.PPStarThresholdAsLog10);
                merDetails.add(t);
                
            }
            //loop on words
            byte[] qw=null;
            while ((qw=sk.getNextByteWord())!=null) {
                //System.out.println("MER "+queryKmerCount+" "+Arrays.toString(qw));

                //skip kmers with too many ambiguities
                if (qw.length==1) {
                    //System.out.println("SKIPPED "+queryKmerCount);
                    queryKmerCount++;
                    skippedKmer++;
                    continue;
                }
                //treat simple kmers
                if (qw.length==session.k) {
                
                    //get Pairs associated to this word
                    //long startT1Time=System.currentTimeMillis();
                    Char2FloatMap.FastEntrySet allPairs =null;
                    //List<Pair> allPairs=null;
                    if (session.states instanceof DNAStatesShifted) {
                        allPairs = session.hash.getPairsOfTopPosition2(session.states.compressMer(qw));
                    } else {
                        allPairs = session.hash.getPairsOfTopPosition2(qw);
                    }
                    //long endT1Time=System.currentTimeMillis();
                    //totalT1Time+=(endT1Time-startT1Time);

                    //word is not present in hash
                    if (allPairs==null) {
                        queryKmerCount++;
                        continue;
                    }
                    queryKmerMatchingDB++;

                    for (ObjectIterator<Char2FloatMap.Entry> iterator = allPairs.iterator(); iterator.hasNext();) {
                        Char2FloatMap.Entry entry = iterator.next();
                        
                        int x=(int)entry.getCharKey();
                        //int nodeId=entry.getNodeId();
                        //we will score only encountered nodes, originalNode registered
                        //at 1st encouter
                        if (C[x]==0) {
                            L.add(x);
                            S[x]+= sk.getMerCount()*session.PPStarThresholdAsLog10;
                        }
                        //count # times originalNode encountered
                        C[x]+=1;
                        //score associated to originalNode x for current read
                        S[x]+= entry.getFloatValue()-session.PPStarThresholdAsLog10;
                        merDetails.get(x)[queryKmerCount]+= entry.getFloatValue()-session.PPStarThresholdAsLog10;
                    }
                
                } else {
                    if (treatAmbiguities) {
                        ambiguousMerTreated++;
                        if (treatAmbiguitiesWithMax) {
                            treatAmbiguitiesWithMax(qw,sk.getMerCount(),L,C,S,queryKmerMatchingDB,queryKmerCount,merDetails);
                        } else {
                            treatAmbiguitiesWithMean(qw,sk.getMerCount(),L,C,S,queryKmerMatchingDB,queryKmerCount,merDetails);
                        }
                    } else {
                        queryKmerCount++;
                        skippedKmer++;
                        continue;
                    }
                }
                
                
                //long endT2Time=System.currentTimeMillis();
                //totalT2Time+=(endT2Time-startT2Time);
                
                //graph of alignment, if asked
//                if(graphAlignment) {
//                   int topPosition =session.hash.getTopPosition(qw);
//                   //System.out.println((queryKmerCount*xSize+topPosition)+"="+allPairs.get(0).getPPStar());
//                   graphDataForTopTuples[2][queryKmerCount*xSize+topPosition]= new Double((float)allPairs.toArray()[0]);                        
//                }
                queryKmerCount++;

            }
//            System.out.println("  AFTER ALIGNMENT");
//            System.out.println("  nodeOccurences:"+Arrays.toString(Arrays.copyOfRange(nodeOccurences,0,nodeOccurences.length)));
//            System.out.println("  nodeScores:"+Arrays.toString(Arrays.copyOfRange(nodeScores,0,nodeOccurences.length)));
            if (merStats) {
            	merStats(session, merDetails, fasta, bwMerStats);
            }
            

            //display alignment result of requested
//            if (graphAlignment) {
//                //graph for topTuplePerPos
//                DefaultXYZDataset datasetForGraph=new DefaultXYZDataset();
//                //datasetForGraph.addSeries(0, graphDataForTopTuples);
//                datasetForGraph.addSeries(0, graphDataForTopTuples);
//                JFrame infos=new JFrame();
//                infos.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
//                infos.setLayout(new GridLayout(1, 1));
//                infos.add(ChartsForNodes.buildReadMatchForANode4("Top_tuple_per_pos of read="+fasta.getHeader(),session.align.getLength(),queryLength, datasetForGraph,session.PPStarThresholdAsLog10,0));
//                infos.setSize(1024, 250);
//                infos.pack();
//                RefineryUtilities.centerFrameOnScreen(infos);
//                infos.setVisible(true);
//            }


            Infos.println("matching_kmers/ambiguous_kmers_treated/skipped_kmers: "+queryKmerMatchingDB+"/"+ambiguousMerTreated+"/"+skippedKmer);
            //long endAlignTime=System.currentTimeMillis();
            //totalAlignTime+=(endAlignTime-startAlignTime);
            //Infos.println("Candidate nodes: "+selectedNodes.size()+" ");      

            //if selectedNodes is empty (no node could be associated)
            //for instance when no query words could be found in the hash
            if ( (L.size()<1) ) {
                Infos.println("Read cannot be placed.");
                //report queries that could not be placed because
                //none of its kmers found in DB
                if (bwNotPLaced != null) {
	            bwNotPLaced.append(fasta.getHeader());
	            bwNotPLaced.newLine();
                }
                continue; //to next query
            }
            
            
            //to check what is the real ordering
//            float[] copyOf = new float[selectedNodes.size()];
//            int c=0;
//            for (int nodeId:selectedNodes) {
//                copyOf[c]=nodeScores[nodeId];
//                c++;
//            }
//            Arrays.sort(copyOf);
//            System.out.println(Arrays.toString(Arrays.copyOfRange(copyOf, 0, 10)));
//            System.out.println(Arrays.toString(Arrays.copyOfRange(copyOf, copyOf.length-10, copyOf.length)));


            //here keep track of the nth best node scores using selection algorithm
            //this should be on average O(k.log(k) + n)
            double allLikelihoodSums = 0.0;
            float bestScore=-1.0f;
            int bestNodeId =-1;
            int secondBestNodeId =-1;
            if (useSelectionAlgo) {
                numberOfBestScoreToConsiderForOutput=keepAtMost;
                //level down keepAtMost is less nodes were selected
                if (L.size()<keepAtMost) {
                    numberOfBestScoreToConsiderForOutput=L.size();
                }
            	
                allLikelihoodSums = fillBestScoreList(S, L, bestScoreList, numberOfBestScoreToConsiderForOutput);

                bestScore=bestScoreList[bestScoreList.length-1].score;
                bestNodeId=bestScoreList[bestScoreList.length-1].nodeId;
            }
            //System.out.println("bestScoreList:"+Arrays.toString(bestScoreList));
            //System.out.println("numberOfBestScoreToConsiderForOutput: "+numberOfBestScoreToConsiderForOutput);
            
            
            
            
            //System.out.println("  AFTER SCORING");
            //System.out.println("  nodeOccurences:"+Arrays.toString(Arrays.copyOfRange(nodeOccurences,150,250)));
            //System.out.println("  nodeScores:"+Arrays.toString(Arrays.copyOfRange(nodeScores,150,250)));
            Infos.println("Best : "+bestNodeId+" (score="+bestScore+")");
            //Infos.println("mapping: ARTree="+session.ARTree.getById(bestNodeId)+" ExtendedTree="+session.extendedTree.getById(session.nodeMapping.get(bestNodeId))+" OriginalTree="+session.originalTree.getById(session.extendedTree.getFakeToOriginalId(session.nodeMapping.get(bestNodeId))));
            //System.out.println("allLikelihoodSums: "+allLikelihoodSums);
//            if (useSelectionAlgo) {
//                System.out.println("selection algo + sort: "+Arrays.toString(bestScoreList));
//            }
            

            // FOLLOWING BLOCK ONLY IF DEBUG OPTION --orinodes IS CALLED
            // AT DB_BUILD. THEN, THIS BLOCK IS EXECUTED AS DB CONTAIN ARTREE
            // NODEIDS. IF DB WAS BUILT WITHOUT THIS OPTION, IT CONTAINS
            // ORIGINAL_TREE NODEIDS, NOT AR_TREE NODEIDS !!!
            // SEARCH BEST NEIGHBOOR GHOST IF BEST NODE IS ORIGINAL NODE
            ////////////////////////////////////////////////////////////////
            
            int extendedTreeId=-1;
            int originalNodeId = -1;
            PhyloNode nodeToTest = null;
            if (session.onlyFakes==false) {       
                //check if this was a fake originalNode or not
                //to do that, retromapping from ARTree to extended tree 
                extendedTreeId=session.nodeMapping.get(bestNodeId);
                originalNodeId = session.extendedTree.getFakeToOriginalId(extendedTreeId);
                nodeToTest = session.extendedTree.getById(extendedTreeId);
                //if this is an original originalNode, select adjacent branch 
                //leading to 2nd best PP*
                if (!nodeToTest.isFakeNode() ) {
                    //System.out.println("############### change best node to neighboors !");                    
                    Infos.println("Current best node is an original node...");
                    PhyloNode firstNode = null;
                    PhyloNode secondNode = null;
                    //if there was no other scored nodes ? this case happened when a single query mer was found in the DB
                    if (secondBestNodeId<0) {
                        //arbitrary choice, take FAKE node which is left son.
                        firstNode = session.ARTree.getById(bestNodeId);
                        secondNode = firstNode.getChildAt(0);
                        bestNodeId =secondNode.getId();
                    } else {
                        firstNode = session.ARTree.getById(bestNodeId);
                        //select node of 2nd best score
                        secondNode = session.ARTree.getById(secondBestNodeId);
                        //get path from this node to bestNodeId
                        //System.out.println("1st node: "+ARTree.getById(bestNodeId)+" 2nd node:"+ARTree.getById(secondNodeId));
                        PhyloTree.Path shortestPath = session.ARTree.shortestPath(session.ARTree.getRoot(), firstNode, secondNode);
                        //System.out.println("Path to second: "+shortestPath.path.size());
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
                
                
                
                
                
            }


//            long endScoringTime=System.currentTimeMillis();
//            totalScoringTime+=endScoringTime-startScoringTime;


            //TO DEBUG
//                System.out.println(fasta.getHeader());
//                System.out.print("bestNodeId:"+bestNodeId);
//                System.out.print("\t\tbestScore:"+bestScore);
//                System.out.println("\t\tbestSCore (normalized by #mers):"+(bestScore/sk.getMaxMerCount()));
//                System.out.println("nodeScores:"+Arrays.toString(nodeScores));
//                System.out.println("nodeOccurences:"+Arrays.toString(nodeOccurences));
//                System.out.println("selectedNodes:"+selectedNodes.keySet().toString());

//            long startWritingTime=System.currentTimeMillis();
            
            //write result in file if a originalNode was hit
            if (csvLog && (bwTSV!=null) ) {
                //only score passing --nsbound if this debug option is set
                if (bestScoreList[bestScoreList.length-1].score>=nsBound) {
                    //OUTPUT n°1: the CSV report of placement
                    //allow in particular to check that nodes were correclty
                    //mappes at every step (original tree, extended tree,
                    //AR modified tree)
                    sb.append(fasta.getHeader().split(" ")[0]).append("\t");
                    if (session.onlyFakes==false) {
                        sb.append(String.valueOf(bestNodeId)).append("\t"); //ARTree nodeID
                        sb.append(String.valueOf(session.ARTree.getById(bestNodeId).getLabel())).append("\t"); //ARTree nodeName
                        sb.append(String.valueOf(extendedTreeId)).append("\t"); //extended Tree nodeID
                        sb.append(String.valueOf(session.extendedTree.getById(extendedTreeId).getLabel())).append("\t"); //extended Tree nodeName
                        sb.append(String.valueOf(originalNodeId)).append("\t"); //edge of original tree (original nodeId)
                        sb.append(String.valueOf(session.originalTree.getById(originalNodeId).getLabel())).append("\t"); //edge of original tree (original nodeName
                        sb.append(String.valueOf(bestScoreList[bestScoreList.length-1].score)).append("\n");
                    } else {
                        sb.append("").append("\t"); //ARTree nodeID
                        sb.append("").append("\t"); //ARTree nodeName
                        sb.append("").append("\t"); //extended Tree nodeID
                        sb.append("").append("\t"); //extended Tree nodeName
                        sb.append(String.valueOf(bestNodeId)).append("\t"); //edge of original tree (original nodeId)
                        sb.append(String.valueOf(session.originalTree.getById(bestNodeId).getLabel())).append("\t"); //edge of original tree (original nodeName
                        sb.append(String.valueOf(bestScoreList[bestScoreList.length-1].score)).append("\n");
                    }
                }
            }
            
            
            //OUTPUT n°2: the JSON placement object (jplace file)
            //2 possibilities in the process of adding a placement:
            //-either do a new "placement" object and add ot to the list
            // of placements (block below is executed)
            //-or a previous placement object corresponding to an
            //identical sequence exists, then we don't the block below
            //as all the alignmnet/placement algo was skipped.
            //(we didn't reached this block)
            
            //only score passing --nsbound if this debug option is set
            if (bestScoreList[bestScoreList.length-1].score>=nsBound) {
                //first we build the "p" array, containing the position/scores of all reads
                JSONArray pMetadata=new JSONArray();
                
                float best = bestScoreList[bestScoreList.length-1].score;
                float lowest = bestScoreList[bestScoreList.length-numberOfBestScoreToConsiderForOutput].score;
                float weightRatioShift = computeWeightRatioShift(lowest, best);
                
                //we create as many lines in "p" block as asked by --keep-at-most and --keep-ratio
                double bestRatio=-1;
                for (int i = bestScoreList.length-1; i>bestScoreList.length-numberOfBestScoreToConsiderForOutput-1; i--) {
                    //System.out.println("i: "+i);
                    //calculate weight_ratio
                    double weigth_ratio=-1;
                    //System.out.println("bestScoreList[i].score-weightRatioShift : "+(bestScoreList[i].score-weightRatioShift));
                    //System.out.println("Math.pow(10.0, (double)(bestScoreList[i].score-weightRatioShift)):"+ Math.pow(10.0, (double)(bestScoreList[i].score-weightRatioShift)));
                    weigth_ratio=computeWeightRatio(bestScoreList[i], weightRatioShift, allLikelihoodSums);
                    
                    //if best score, memorize this ratio
                    if (i==bestScoreList.length-1) {
                        bestRatio=weigth_ratio;
                    }
                    //System.out.println("Likelihood weight ratio: "+weigth_ratio);
                    //take into account option --keep-factor
                    if (i<bestScoreList.length-1 && weigth_ratio<(bestRatio*keepFactor)) {
                        break;
                    }
                    //in pplacer/EPA several placements can be associated to a query
                    //we input only the best one, but that can be changed in the future
                    //"distal_length","like_weight_ratio","pendant_length","edge_num","likelihood"
                    JSONArray placeColumns=new JSONArray();
                    if (guppyCompatible) {
//                        System.out.println("bestScoreList[i].score : "+bestScoreList[i].score);
//                        System.out.println("bestScoreList[i].nodeId : "+bestScoreList[i].nodeId);
//                        System.out.println(session.originalTree.getById(bestScoreList[i].nodeId));
//                        System.out.println(session.originalTree.getById(bestScoreList[i].nodeId).getJplaceEdgeId());
                        placeColumns.add(session.originalTree.getById(bestScoreList[i].nodeId).getBranchLengthToAncestor()/2); //distal_length
                        placeColumns.add(session.originalTree.getById(bestScoreList[i].nodeId).getJplaceEdgeId()); // 1. edge of original tree (original nodeId=edgeID)
                        placeColumns.add(weigth_ratio); // 3. like_weight_ratio column of ML-based methods
                        placeColumns.add(bestScoreList[i].score); // 2. PP*
                        placeColumns.add(0.0); //pendant_length
                    } else {
                        placeColumns.add(session.originalTree.getById(bestScoreList[i].nodeId).getJplaceEdgeId()); // 1. edge of original tree (original nodeId=edgeID)
                        placeColumns.add(bestScoreList[i].score); // 2. PP*
                        placeColumns.add(weigth_ratio); // 3. like_weight_ratio column of ML-based methods
                        //fake fields for compatibility with current tools (guppy, archeopteryx)
                        //should be provided as an option
                        placeColumns.add(session.originalTree.getById(bestScoreList[i].nodeId).getBranchLengthToAncestor()/2f); //distal_length
                        placeColumns.add(0.0); //pendant_length
                    }
                    pMetadata.add(placeColumns);
                }
                placement.put("p", pMetadata);

                //second we build the "nm" array, containing the read identifier
                //and its multiplicity. 
                //Only one element as this is the 1st read,
                //identical reads will be detected before the placement computation
                //and will simply be added to this list (see checksum block before 
                //the alignment/scoring blocks).
                JSONArray allIdentifiers=new JSONArray();
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

//            long endWritingTime=System.currentTimeMillis();
//            totalWritingTime+=endWritingTime-startWritingTime;


            //push the stringbuffer to the CSV bufferedwriter every 25000 sequences
            if (bwTSV!=null) {
                if ((queryCounter%25000)==0) {
                    int size=sb.length();
                    bwTSV.append(sb);
                    sb=null;
                    sb=new StringBuffer(size);
                }   
            }



            //reset the scoring vectors
            long startResetTime=System.currentTimeMillis();
            for (Integer nodeId : L) {
                S[nodeId]=0.0f;
                C[nodeId]=0;
            }
            for (int i = 0; i < bestScoreList.length; i++) {
                bestScoreList[i]=new Score(-1, Float.NEGATIVE_INFINITY);
            }
            selectedNodeSize+=L.size();
            L.clear();
//            long endResetTime=System.currentTimeMillis();
//            totalResetTime+=endResetTime-startResetTime;
            

        }
        
        if (merStats) {
            bwMerStats.close();
        }
        
        //flush to disk the last CSV buffer
        if (bwTSV != null) {
            bwTSV.append(sb);
        }

            
//        long endTotalPlacementTime=System.currentTimeMillis();
//        System.out.println("############################################################");
//        System.out.println("Checksum registry took in total: "+totalChecksumTime+" ms");
//        System.out.println("Knife took in total: "+totalKnifeTime+" ms");
//        System.out.println("kmer matching took in total: "+totalAlignTime+" ms");
//        System.out.println(" T1 took in total: "+totalT1Time+" ms");
//        System.out.println(" T2 took in total: "+totalT2Time+" ms");
//        System.out.println("Scoring took in total: "+totalScoringTime+" ms");
//        System.out.println("Writing .csv and .jplace took in total: "+totalWritingTime+" ms");
//        System.out.println("Reset took in total: "+totalResetTime+" ms");
//        System.out.println("------------------------------------------------------------");
//        System.out.println("(per placement average)");
//        System.out.println("Checksum registry took on average: "+((0.0+totalChecksumTime)/placedQueryCounter)+" ms");
//        System.out.println("Knife took on average: "+((0.0+totalKnifeTime)/placedQueryCounter)+" ms");
//        System.out.println("k-mer match/score took on average: "+((0.0+totalAlignTime)/placedQueryCounter)+" ms");
//        System.out.println(" T1 took on average: "+((0.0+totalT1Time)/placedQueryCounter)+" ms");
//        System.out.println(" T2 took on average: "+((0.0+totalT2Time)/placedQueryCounter)+" ms");
//        System.out.println("Score correction/ratio took on average: "+((0.0+totalScoringTime)/placedQueryCounter)+" ms");
//        System.out.println("Writing .jplace and/or .csv took on average: "+((0.0+totalWritingTime)/placedQueryCounter)+" ms");
//        System.out.println("Reset took on average: "+((0.0+totalResetTime)/placedQueryCounter)+" ms");
//        System.out.println("L average: "+((0.0+selectedNodeSize)/placedQueryCounter)+" nodes");
//        System.out.println("------------------------------------------------------------");
//        System.out.println("Placement Process (without DB load) took: "+(endTotalPlacementTime-startTotalPlacementTime)+" ms");
//        System.out.println("############################################################");
        
        return queryCounter;
    }
    
    /**
     * 
     * @param w concatenate of kmers of size k
     * @param Q #mer in query
     * @param L list of node encountered in main placement loop
     * @param C occurences per nodeId
     * @param S score per nodeId
     * @param queryKmerMatchingDB debug counter
     */
    private void treatAmbiguitiesWithMean(byte[] w, int Q, ArrayList<Integer> L, int[] C, float[] S,int queryKmerMatchingDB, int queryKmerCount, ArrayList<float[]> merDetails) {
        float[] S_amb=new float[C.length];
        int[] C_amb=new int[C.length];
        ArrayList<Integer> L_amb=new ArrayList<>();
        int W_size=w.length/session.k;
        boolean matched=false;
        //extract probas for each w_prime
        byte[] w_prime=null;
        for (int i=0; i<W_size;i++) {
            w_prime=Arrays.copyOfRange(w, i*session.k, (i+1)*session.k);
            Char2FloatMap.FastEntrySet allPairs =null;
            if (session.states instanceof DNAStatesShifted) {
                allPairs = session.hash.getPairsOfTopPosition2(session.states.compressMer(w_prime));
            } else {
                allPairs = session.hash.getPairsOfTopPosition2(w_prime);
            }
            if (allPairs==null) {
                continue;
            }
            matched=true;
            allPairs.stream().forEach( (entry) -> {
                int x=(int)entry.getCharKey();
                if (C_amb[x]==0) {
                    L_amb.add(x);
                }
                C_amb[x]+=1;
                S_amb[x]+=Math.pow(10,entry.getFloatValue());
            });
        }
        //if one of the w_prime matched increase debug counter
        if (matched) {queryKmerMatchingDB++;}
        //compute mean and add it to scoring/counting tables
        for (int i = 0; i < L_amb.size(); i++) {
            int x = L_amb.get(i);
                if (C[x]==0) {
                    L.add(x);
                    S[x]=Q*session.PPStarThresholdAsLog10;
                }
                C[x]+=1;
                float avgProba=(S_amb[x] + (W_size-C_amb[x])*session.PPStarThreshold) / W_size;
                S[x]+=Math.log10(avgProba)-session.PPStarThresholdAsLog10;
                C_amb[x]=0;
                merDetails.get(x)[queryKmerCount]+= Math.log10(avgProba)-session.PPStarThresholdAsLog10;
        }
        
    }
    
    /**
     * 
     * @param w concatenate of kmers of size k
     * @param Q #mer in query
     * @param L list of node encountered in main placement loop
     * @param C occurences per nodeId
     * @param S score per nodeId
     * @param queryKmerMatchingDB debug counter
     */
    private void treatAmbiguitiesWithMax(byte[] w, int Q, ArrayList<Integer> L, int[] C, float[] S,int queryKmerMatchingDB, int queryKmerCount, ArrayList<float[]> merDetails) {
        float[] S_amb=new float[C.length];
        int[] C_amb=new int[C.length];
        
        ArrayList<Integer> L_amb=new ArrayList<>();
        
        int W_size=w.length/session.k;
        
        //extract probas for each w_prime and retain max
        byte[] w_prime=null;
        boolean matched=false;
        for (int i=0; i<W_size;i++) {
            w_prime=Arrays.copyOfRange(w, i*session.k, (i+1)*session.k);
            Char2FloatMap.FastEntrySet allPairs =null;
            if (session.states instanceof DNAStatesShifted) {
                allPairs = session.hash.getPairsOfTopPosition2(session.states.compressMer(w_prime));
            } else {
                allPairs = session.hash.getPairsOfTopPosition2(w_prime);
            }
            if (allPairs==null) {
                continue;
            }
            matched=true;
            allPairs.stream().forEach( (entry) -> {
                int x=entry.getCharKey();
                if (C_amb[x]==0) {
                    L_amb.add(x);
                    S_amb[x]=entry.getFloatValue();
                }
                C_amb[x]+=1;
                if (entry.getFloatValue() > S_amb[x]) {
                    S_amb[x]=entry.getFloatValue();
                }
            });
        }
        //if one of the w_prime matched increase debug counter
        if (matched) {queryKmerMatchingDB++;}
        //add max to scoring/counting tables
        for (int i = 0; i < L_amb.size(); i++) {
            int x = L_amb.get(i);
                if (C[x]==0) {
                    L.add(x);
                    S[x]=Q*session.PPStarThresholdAsLog10;
                }
                C[x]+=1;
                S[x]+=S_amb[x]-session.PPStarThresholdAsLog10;
                C_amb[x]=0;
                merDetails.get(x)[queryKmerCount]+= S_amb[x]-session.PPStarThresholdAsLog10;

        }
        
    }
    
    
    public boolean parallelProcessReads() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }


    /**
     * simple class to wrap the score and likelihood_weight_ratio associated to the n best nodes
     */
    public static class Score implements Comparable<Score>{
        public int nodeId;
        public float score=Float.NEGATIVE_INFINITY;

        public Score(int nodeId, float score) {
            this.nodeId = nodeId;
            this.score = score;
        }
        
        
        
        @Override
        public int compareTo(Score o) {
            return Float.compare(this.score, o.score);
        }  

        @Override
        public String toString() {
            return nodeId+":"+score;
        }
        
        
    }
    
    private class ScoreComparator implements Comparator<Score> {

//        float[] nodeScores=null;
//        
//        public ScoreComparator(float[] nodeScores) {
//            
//        
//        }
        @Override
        public int compare(Score o1, Score o2) {
            return Float.compare(o1.score, o2.score);
        }
    }
    
    
    
    
}
