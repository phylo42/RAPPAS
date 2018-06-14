/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import etc.Infos;
import java.util.ArrayList;
import java.util.Arrays;
//import main_v2.SessionNext_v2;
import main_v2.SessionNext_Triplet;

/**
 * third version of word explorer
 * (branch and bound algo to build ancestral k-mers)
 * 
 * which takes into account empty columns to generate k-mers
 * the option 'onlyOneJump', restrain the kmer generation to a single jump
 * when several jumps are possible inside the kmer
 * 
 * AND
 * 
 * do hash insertions straight from the branch and bound algorithm to avoid the
 * transfer to placement loop through a useless arraylist which can require
 * too much space when many k-mers are propable
 * 
 * AND 
 * 
 * session passed to constructor to simplify code
 * 
 * @author ben
 */
public class WordExplorer_Triplet {

     
    //external data used in the recursion
    //SessionNext_v2 session=null;
    SessionNext_Triplet session=null;
    int extTreeId=-1;
    int originalId=-1;    
    
    //table to register proba proportion that is left on 
    float[] probaLeftPerCol=null;
    
    //variables for the recursive algorithm
    float currentLogSum=0.0f;
    byte[] word=null;
    boolean boundReached=false;
    boolean wordCompression=false;
    int refPosition=-1;
    int current_k=-1;
    int nodeId=-1;
    int idxOfFirstJump=-1;
    int generateTupleCount=-1;
    
    //represents  gap intervals
    boolean doGapJumps=false;
    ArrayList<Integer>[] gapIntervals =null;
    boolean limitTo1Jump=false;

    /**
     *
     * @param session
     * @param refPosition
     * @param nodeId
     * @param wordCompression
     * @param doGapJumps
     * @param limitTo1Jump
     */
//    public WordExplorer_v3( SessionNext_v2 session,
//                            int refPosition,
//                            int nodeId,
//                            boolean wordCompression,
//                            boolean doGapJumps,
//                            boolean limitTo1Jump) {
    public WordExplorer_Triplet( SessionNext_Triplet session,
                            int refPosition,
                            int nodeId,
                            boolean wordCompression,
                            boolean doGapJumps,
                            boolean limitTo1Jump) {
        //data loaded from session
        this.session=session;
        this.gapIntervals = session.align.getGapIntervals();
        //exploration parameters 
        this.refPosition=refPosition;
        this.nodeId=nodeId;
        this.wordCompression=wordCompression;
        this.doGapJumps=doGapJumps;
        this.limitTo1Jump=limitTo1Jump;
        //accessory variables used in recursion
        this.word=new byte[session.k]; //generated words are of size k
        this.current_k=0;//defines to k-th position where is currently position recursion
        //ARTree id mapped to originals ids, as has will contain these
        this.extTreeId=session.nodeMapping.get(nodeId);
        this.originalId=session.extendedTree.getFakeToOriginalId(extTreeId);

    }
    
    /**
     * 
     * @param i
     * @param j
     */
    public void exploreWords(int i, int j) {
        //for the algorithm below, consider ppSet as a (i,j) matrix
        //- with i the position on the ref alignment, with refPosition=i at start
        //  and refPosition+k being the maximum where i=k-1.
        //- with j the current state, in descending ordered by their PP, from j=0
        //  to a maxumim j=(#states)
        

        //Infos.println("ENTERING : "+i+" "+j+" , current_k:"+current_k);
        //if after the alignment limit (can happen because of gap jumps)
        //then return, assigning a value to this kmer is impossible.
        if (i>session.parsedProbas.getSiteCount()-1) {
            return;
        }
        //reset the jump counter if we come back to 1st mer position
        if (current_k==0) {
            idxOfFirstJump=-1;
        }
        
        word[current_k]=session.parsedProbas.getState(nodeId, i, j);
        //Infos.println("sumCurrentWord="+currentLogSum+"+"+ppSet.getPP(nodeId, i, j));
        currentLogSum+=session.parsedProbas.getPP(nodeId, i, j);
        boundReached = currentLogSum<session.PPStarThresholdAsLog10;
        
        //register word if k-th position
        if (current_k==session.k-1) {
            //register word
            if (!boundReached) {
                byte[] w=null;
                if (wordCompression) {
                    w=session.states.compressMer(Arrays.copyOf(word, word.length));
                } else {
                    w=Arrays.copyOf(word, word.length);
                }
                //System.out.println("REGISTER: "+Arrays.toString(word)+" log10(PP*)="+currentLogSum+" nodeId="+originalId+" refPosition="+refPosition);
                session.hash.addTuple(w, currentLogSum, originalId, refPosition);
                generateTupleCount+=1;
            }
            //decrease before return
            //Infos.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
            currentLogSum-=session.parsedProbas.getPP(nodeId, i, j);
            //force return to parent
            return;
        } else {
            
            //go down recursively, testing each state at position i+1
            for (int j2 = 0; j2 < session.parsedProbas.getStateCount(); j2++) {
                if (boundReached) {break;}
                
                //do this exploration, do not consider gap interval
                //Infos.println("EXPORATION NORMAL");
                //we move to next kmer position
                current_k++;
                exploreWords(i+1, j2);
                current_k--;
                //Infos.println("OUT from : "+(i+1)+" "+j2+" , current_k:"+current_k+" , currently HERE: "+i+" "+j);
                //redo exploration, by jumping over the gap intervals if any
                //do condition on i because i+1 cannot be outside the alignment
                if (doGapJumps && i<session.parsedProbas.getSiteCount()-1) {
                    //verify if next i is associated to a gap interval
                    if (gapIntervals[i+1]!=null) { 
                        //do all jumps combinations 
                        if (!limitTo1Jump) {
                            for (int i_interval = 0; i_interval < gapIntervals[i+1].size(); i_interval++) {
                                //we move to next kmer position and explore
                                current_k++;
                                exploreWords((i+1)+gapIntervals[i+1].get(i_interval), j2); //jump from i to i+1+gapInterval_length
                                current_k--;
                            }
                        //or allow only 1 jump
                        } else {
                            if ((idxOfFirstJump==-1) ) { //do a jump only if not previously done
                                idxOfFirstJump=i;
                                for (int i_interval = 0; i_interval < gapIntervals[i+1].size(); i_interval++) {
                                    //Infos.println("EXPLORATION WITH GAP JUMP: gap_i+1="+(i+1)+" length="+gapIntervals[i+1].get(i_interval));
                                    //Infos.println("NEXT i jump to --> "+((i+1)+gapIntervals[i+1].get(i_interval)));
                                    //we move to next kmer position
                                    current_k++;
                                    exploreWords((i+1)+gapIntervals[i+1].get(i_interval), j2); //jump from i to i+1+gapInterval_length
                                    current_k--;
                                    //Infos.println("OUT from : "+(i+gapIntervals[i+1].get(i_interval))+" "+j2+" , current_k:"+current_k+" , currently HERE: "+i+" "+j);
                                }
                            } else {
                                //Infos.println("DO NOT JUMP, already jumped at idxOfFirstJump="+idxOfFirstJump);
                            }
                        }
                    }
                }
            }
            
            
            
        }
        //decrease before natural return
        //Infos.println("sumCurrentWord="+currentLogSum+"-"+ppSet.getPP(nodeId, i, j));
        currentLogSum-=session.parsedProbas.getPP(nodeId, i, j);
    }
    
    public int getGeneratedTupleCount() {
        return this.generateTupleCount;
    }
    
    
}
