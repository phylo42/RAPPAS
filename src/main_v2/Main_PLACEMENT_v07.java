/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import core.AAStates;
import core.DNAStates;
import core.States;
import core.algos.AlignScoringProcess;
import core.algos.SequenceKnife;
import etc.Environement;
import etc.Infos;
import inputs.FASTAPointer;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import tree.NewickWriter;

/**
 * Algo in this version:
 * 
 * 1. all words of the query are searched in the hash
 *    --> we storeFullHash the corresponding nodes in a list
    --> positions are used to define the alignment
 * 2. using the alignment, for all nodes stored in the previous list, 
 *    we search query words associated to the positions 
 *    --> we build diagsums for this originalNode selection nodes (1 per originalNode)
    --> we keep as peeks only diagsum position holding more than X words
 *    --> the sum of these peeks define the final score and the placement
 * 
 * 
 * 
 * @author ben
 */
public class Main_PLACEMENT_v07 {
    

    
    //related to debug operation --dbinram
    SessionNext_v2 session=null;
    boolean dbInRAM=false;
    
    /**
     * call this constructor when using normal placement operation
     */
    public Main_PLACEMENT_v07() { }
    

    /**
     * call this constructor when doing placement via debug operation --dbinram
     * @param session 
     * @param dbInRAM 
     */
    public Main_PLACEMENT_v07(SessionNext_v2 session, boolean dbInRAM) {
        this.session=session;
        this.dbInRAM=dbInRAM;
        
    }
    
    
    

    
    /**
     * 
     * @param q
     * @param db if null, will attempt debug operation (--dbinram)
     * @param workDir
     * @param callString
     * @param nsBound
     * @return 
     */
    public int doPlacements(File q, File db, File workDir, String callString,Float nsBound) {

        try {
                        
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //DEBUG PARAMETERS
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
              
            //minimum read/ref overlap,in bp. When not respected, read not reported
            int minOverlap=100;
            //word sampling method
            int queryWordSampling=SequenceKnife.SAMPLING_LINEAR;

            //debug 
            int queryLimit=100000000;
            

            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            //LOADING ALL DATA FROM THE INPUT DATABASE
            long startLoadTime=System.currentTimeMillis();
            
            //LOAD SESSION//////////////////////////////////////////////////////
            //logs
            String logPath=workDir+File.separator+"logs"+File.separator;
            //ancestral reconstruciton
            String ARPath=workDir+File.separator+"AR"+File.separator;
            //session itself
            boolean loadHash=true;
            if (!dbInRAM) {
                System.out.println("Loading ancestral words DB... ("+db.getName()+")");
                session= SessionNext_v2.load(db,loadHash);
            } 
            //mers size
            Infos.println("k="+session.k);
            Infos.println("min_k="+session.minK);
            //site and word posterior probas thresholds
            Infos.println("factor="+session.alpha);
            Infos.println("sitePPThreshold="+session.stateThreshold+" (for info, unused)");
            Infos.println("PPStarThreshold="+session.PPStarThreshold);
            Infos.println("PPStarThreshold(log10)="+session.PPStarThresholdAsLog10);
            
            
            //type of Analysis AA or DNA////////////////////////////////////////
            States s=session.states; 

            //PREPARE DIRECTORIES///////////////////////////////////////////////
            if (!workDir.exists()) {workDir.mkdir();}
            if (!new File(logPath).exists()) {new File(logPath).mkdir();}
            if (!new File(ARPath).exists()) {new File(ARPath).mkdir();}
            
            
            //SOME BASIC DISPLAY TO CONTROL SESSION LOAD////////////////////////
            Infos.println(session.align.describeAlignment(false));
            Infos.println("# nodes in the tree: "+session.ARTree.getNodeCount());
            Infos.println("# leaves in the tree: "+session.ARTree.getLeavesCount());
            Infos.println("# internal nodes in the tree: "+session.ARTree.getInternalNodesByDFS().size());
            
            //LOAD THE POSTERIOR PROBAS/////////////////////////////////////////
//            PProbasSorted pprobas = session.parsedProbas;
//            //to raidly check that sorted probas are OK
//            Infos.println("NodeId=0, 5 first PP:"+Arrays.deepToString(pprobas.getPPSet(0, 0, 5)));
//            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateSet(0, 0, 5)));
//            Infos.println("NodeId=0, 5 first states:"+ Arrays.deepToString(pprobas.getStateIndexSet(0, 0, 5)));
            
            String[] elts=db.getName().split("\\.");
            String dbSize=elts[elts.length-1];
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


            ////////////////////////////////////////////////
            //LOADING THE QUERIES PROVIDED BY USER
            
            System.out.println("Analyzing query sequences...");
            FASTAPointer fp=new FASTAPointer(q, false);
            int totalQueries=fp.getContentSize();
            Infos.println("Input fasta contains "+totalQueries+" sequences");
            Infos.println("Mean sequence size: "+fp.getContentMean());
            fp.resetPointer();
            //FileWriter fw =new FileWriter(new File(logPath+"queries.fasta"));
            
            
            ////////////////////////////////////////////////////////////////////
            //PREPARE THE WRITER FOR OUTPUT IN TSV FORMAT
            ////////////////////////////////////////////////////////////////////
            int bufferSize=2097152; // buffer of 2mo
            //placement results
            BufferedWriter bwTSVPlacement=new BufferedWriter(new FileWriter(new File(logPath+"placements_"+q.getName()+"_"+dbSize+".tsv")),bufferSize);
            

            ////////////////////////////////////////////////////////////////////
            //PREPARE STRUCTURE (JSON OBJECT) FOR JPLACE OUTPUT 
            //this need to be done outside the alignment/scoring loop
            //because identical sequences (same score) will be in the same 
            //json "placement" (p) object
            //with library json.simple
            
            //top level of the json tree
            JSONObject top=new JSONObject(); 
            LinkedHashMap topMap=new LinkedHashMap();
            //object tree (mandatory)
            //write a newick version of the loaded original tree for control
            NewickWriter nw=new NewickWriter();
            String relaxedTreeForJplace=nw.getNewickTree(session.originalTree, true, true, true);
            nw.close();
            topMap.put("tree",relaxedTreeForJplace);  
            //we do an array of placement object
            //all identical reads with be injected in the same object
            JSONArray placements=new JSONArray();
            
            

                        
            
            ////////////////////////////////////////////////////////////////////
            //PLACEMENT ALGO ITSELF
            ////////////////////////////////////////////////////////////////////
            //NORMALIZED SCORE BELOW THE CALIBRATION RESULT WILL NOT BE OUTPUT
            //IN THE JPLACE OUPUT
            //TODO: HERE EXTEND WITH PARALLELISM
            AlignScoringProcess asp=null;
            if (nsBound!=null) {  //norm score bound was set manually via command line
                System.out.println("User provided nsBound !");
                asp=new AlignScoringProcess(session,nsBound, queryLimit);
            } else {
                asp=new AlignScoringProcess(session,session.calibrationNormScore, queryLimit);
            }
            int queryCounter=asp.processQueries(fp,placements,bwTSVPlacement,queryWordSampling,minOverlap);
            //close TSV logs
            bwTSVPlacement.close();
            fp.closePointer();

            
            ////////////////////////////////////////////////////////////////////
            //JSON OUTPUT
            ////////////////////////////////////////////////////////////////////
            //finish the json strucutre and output it to a file
            //associate the list of placements to the top level
            top.put("placements",placements);
            //object version (mandatory
            top.put("version",3);
            //object metadata (mandatory): sub-objet "version is mandatory"
            JSONObject invoc=new JSONObject();
            invoc.put("invocation", "viromeplacer"+callString);
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
            fList.add("edge_num"); //i.e equal to the id of the son originalNode
            //fList.add("edge_label"); //i.e equal to the id of the son originalNode
            fList.add("likelihood"); //rename to likelihood even if it is not, but for compatibility with other programs
            top.put("fields", fList);
            
            //put all the elements in the top JSON object
            top.putAll(topMap);
            String out=top.toJSONString();

            //just some basic formatting
            out=out.replaceAll("\\},\\{", "\n\\},\\{\n\t"); //},{
            out=out.replaceAll("\\],\"","\\],\n\t\"");   //],"
            out=out.replaceAll("\\]\\}\\],", "\\]\n\\}\n\\],\n"); //]}]
            //out=out.replace("]},", "]},"); //]}
            
            FileWriter fwJSON =new FileWriter(new File(logPath+File.separator+"placements_"+q.getName()+"_"+dbSize+".jplace"));
            fwJSON.append(out);
            fwJSON.close();
            
            long endTotalTime=System.currentTimeMillis();
            
            ////////////////////////////////////////////////////////////////////
            //MAIN ALGO END !!!!
            ////////////////////////////////////////////////////////////////////
            

            Infos.println("#######################################################################");
            Infos.println("### DONE, placement execution took (excluding DB load): "+(endTotalTime-startTotalTime)+" ms");
            Infos.println("#######################################################################");
            //just for coherent output, close the percentage
            System.out.println(queryCounter+"/"+totalQueries+" queries placed ("+(((0.0+queryCounter)/totalQueries)*100)+"%)");
            //just for coherent output, close the percentage
            System.out.println(placements.size()+" significative placements reported in JPlace output ("+((0.0+placements.size()/totalQueries)*100)+"%)");            
            Infos.println("#######################################################################");


            return queryCounter;
            
            
            
        } catch (IOException ex) {
            Logger.getLogger(Main_PLACEMENT_v07.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
        }
        
        
        

    }
    

    /**
     * retrieve the indexes of the nummax highest values in an array 
     * @param orig
     * @param nummax
     * @return 
     */
    public static int[] indexesOfTopElements(float[] orig, int nummax) {
        float[] copy = Arrays.copyOf(orig,orig.length);
        Arrays.sort(copy); //will put all scores=0 at the end...
        System.out.println(Arrays.toString(copy));
        float[] honey = Arrays.copyOfRange(copy,copy.length-nummax,copy.length);
        System.out.println(Arrays.toString(honey));
        int[] result = new int[nummax];
        int resultPos = 0;
        for(int i = 0; i < honey.length; i++) {
            float onTrial = honey[i];
            System.out.println(onTrial);
            int index = Arrays.binarySearch(orig,onTrial);
            if (index<0) {
                System.out.println("ERROR, not found!");
                System.exit(1);
            }
            result[resultPos++] = index;
        }
        return result;
    }
    

    
    
    
    
    
}