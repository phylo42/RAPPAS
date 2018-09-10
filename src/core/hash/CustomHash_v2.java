/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import core.States;
import etc.Infos;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * In this second hash version, the Nodes composing the buckets contain a table
 * of positions pointing to a list of Pairs representing a
 * (nodeId, PP*)
 * @author ben
 */
public class CustomHash_v2 implements Serializable{
    
    private static final long serialVersionUID = 7000L;
    
    public static final int NODES_POSITION=1;
    public static final int NODES_UNION=2;
    
    int nodeType=NODES_UNION;
    
    HashMap<ByteArrayWrapper,HashPointer> hash=null;

    ByteArrayWrapper wBuf=null;
    
    int maxPointerSize=-1;
    long addTupleCalls=0l;
    long totalContainsTime=0l;
    long totalRegisterTime=0l;
    /**
     *
     * @param k
     * @param s
     * @param nodeType one of CustomHash_v2.NODES_UNION or CustomHash_v2.POSITION_UNION
     */
    public CustomHash_v2(int k, States s, int nodeType) {        
        this.maxPointerSize=new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)).intValue();
        this.nodeType=nodeType;
        this.wBuf=new ByteArrayWrapper(new byte[k]);
        //internal tests showed that with k<14 we generally get at least 75% of the possible k-mers
        hash=new HashMap<>(new Double(maxPointerSize*0.75).intValue(),0.8f);
    }
    
    /**
     * only entry point to fill the hash (used at DB generation)
     * @param w
     * @param PPStar
     * @param nodeId
     * @param refPos 
     */
    public void addTuple(byte[] word, float PPStar,int nodeId,int refPos) {

       long containsStart=System.nanoTime();
       wBuf.setArray(word);
       HashPointer cn=null;
       if ((cn=hash.get(wBuf))==null) {
            switch (nodeType) {
                case NODES_UNION:
                    //System.out.println("new word attached to UnionPointer");
                    hash.put(wBuf, cn=new UnionPointer());
                    break;
                case NODES_POSITION:
                    //System.out.println("new word attached to PositionPointer");
                    hash.put(wBuf, cn=new PositionPointer());
                    break;
            }
            
        }
        long containsEnd=System.nanoTime();
        totalContainsTime+=(containsEnd-containsStart);
        
        long registerStart=System.nanoTime();
        //System.out.println("update HashPointer of type: "+hash.get(w).getClass().toString());
        cn.registerTuple(nodeId, refPos, PPStar);
        long registerEnd=System.nanoTime();
        totalRegisterTime+=(registerEnd-registerStart);
        
        addTupleCalls++;
    }

    /**
     * for debug purposes only (access to CustomHashMap hashtable)
     * @return 
     */
    public HashMap<ByteArrayWrapper, HashPointer> getHash() {
        return hash;
    }
    
    /**
     * to know which type of nodes is backing this hash, 
     * one of CustomHash_v2.NODES_UNION or CustomHash_v2.NODES_POSITION
     * @return 
     */
    public int getHashType() {
        return this.nodeType;
    }

    /**
     * best (nodeId,PP*) associated to this word
     * @param w
     * @return 
     */
    public Pair getTopPair(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(wBuf.setArray(w)))!=null)
            return cn.getBestPair();
        else
            return null;
    }    
    
    /**
     * retrieves only (nodeId,PP*) pairs stored under the position
     * associated to the best PP*
     * @param w
     * @return null if word not in present in hash
     */
    public List<Pair> getPairsOfTopPosition(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(wBuf.setArray(w)))!=null) {
            return cn.getPairList(cn.getBestPosition());
        } else {
            return null;
        }
    }  
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return 
     */
    public int[] getPositions(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(wBuf.setArray(w)))!=null) {
            return cn.getPositions();
        } else {
            return null;
        }
    }
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return -1 if word not in hash
     */
    public int getTopPosition(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(wBuf.setArray(w)))!=null) {
            return cn.getBestPosition();
        } else {
            return -1;
        }
    }
    
    
    
    public void sortData() {
        double startTime=System.currentTimeMillis();
        AtomicInteger pairCount=new AtomicInteger(0);
        hash.values().stream().forEach( l -> {l.sort();pairCount.addAndGet(l.getPairCountInTopPosition());} );
        double endTime=System.currentTimeMillis();
        Infos.println("# Pairs sorted in top positions: "+pairCount.get());
        Infos.println("Pair sorting took "+(endTime-startTime)+" ms");
    }
    
    
    /**
     * pairs whatever the associated position
     * @param w
     * @return 
     */
    public List<Pair> getPairs(byte[] w) {
        ArrayList<Pair> l =new ArrayList();
        for (int p: hash.get(wBuf.setArray(w)).getPositions()) {
            l.addAll(hash.get(wBuf.setArray(w)).getPairList(p));
        }
        return l;
    }
    
    /**
     * pairs selected by associated position
     * @param w
     * @param position
     * @return 
     */
    public List<Pair> getPairs(byte[] w,int position) {
        if (hash.containsKey(wBuf.setArray(w))) {
            return hash.get(wBuf.setArray(w)).getPairList(position);
        } else {
            return null;
        }
    }

    public Set<ByteArrayWrapper> keySet() {
        return hash.keySet();
    }
    
    
    
    
    /**
     * empty all the positions which where not associated to the best PP*
     */
    public void reduceToMediumHash() {
        
        hash.keySet().stream().forEach((next) -> {
            ((PositionPointer)hash.get(next)).clearPairsOfWorsePositions();
        });
        
    }

    /**
     * in top position, empties all the positions which where not associated to the best PP*,
     * and retains only X pairs at the best position
     * @param X
     */
    @Deprecated
    public void reducetoSmallHash(int X) {
        List<ByteArrayWrapper> collect = hash.keySet()  .stream() 
                                            //.peek((w)->System.out.println("REDUCING:"+w))
                                            .filter((w) -> ((PositionPointer)hash.get(w)).limitToXPairsPerPosition(X))
                                            //.peek((w)->System.out.println("TRASHED!:"+w))
                                            .collect(Collectors.toList());
        collect.stream().forEach((w)-> {hash.remove(w);});
        collect=null;
    }
    
    /**
     * discards all words where top positions is associated to more than X nodes
     * @param X
     */
    public void reducetoSmallHash_v2(int X) {
        assert X>0;
        List<ByteArrayWrapper> collect = hash.keySet()  .stream() 
                                            //.peek((w)->System.out.println("REDUCING:"+w))
                                            .filter((w) -> hash.get(w).getPairCountInTopPosition()>X)
                                            //.peek((w)->System.out.println("TRASHED!:"+w))
                                            .collect(Collectors.toList());
        collect.stream().forEach((w)-> {hash.remove(w);});
        collect=null;
    }    
    
    
    public long getTotalContainsTime() {
        return this.totalContainsTime;
    }

    public long getTotalRegisterTime() {
        return this.totalRegisterTime;
    }
    
    public long getTotalAddTupleCalls() {
        return this.addTupleCalls;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
 
    
    
    
    
    
    

    

    

    
    
    
    
}
