/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import core.States;
import etc.Infos;
import it.unimi.dsi.fastutil.chars.Char2FloatMap;
import it.unimi.dsi.fastutil.chars.Char2FloatOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenCustomHashMap;
import java.io.Serializable;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

/**
 * In this second hash version, the Nodes composing the buckets contain a table
 * of positions pointing to a list of Pairs representing a
 * (nodeId, PP*)
 * @author ben
 */
public class CustomHash_v5_FastUtil811 implements Serializable{
    
    private static final long serialVersionUID = 7000L;
    
    public static final int NODES_POSITION=1;
    public static final int NODES_UNION=2;
    
    int nodeType=NODES_UNION;
    
    Char2FloatOpenHashMap preparedNovelMap=null;
    Object2ObjectOpenCustomHashMap<byte[],Char2FloatOpenHashMap> hash;
    
    int maxCapacitySize=-1;

    /**
     *
     * @param k
     * @param s
     * @param nodeType one of NODES_UNION or POSITION_UNION
     */
    public CustomHash_v5_FastUtil811(int k, States s, int nodeType) {        
        this.maxCapacitySize=new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)).intValue();
        this.nodeType=nodeType;
        //internal tests showed that with k<14 we generally get at least 75% of the possible k-mers
        this.hash=new Object2ObjectOpenCustomHashMap<>(
                            new Double(maxCapacitySize/8).intValue(),  //intial capacity
                            0.8f, //inital load factor                
                            new HashStrategy()
                        );
        this.hash.defaultReturnValue(null);
    }
    
    /**
     * only entry point to fill the hash (used at DB generation)
     * @param word
     * @param PPStar
     * @param nodeId
     * @param refPos 
     */
    public void addTuple(byte[] word, float PPStar,int nodeId,int refPos) {
        
        
        Char2FloatOpenHashMap existingMap = hash.putIfAbsent(word, preparedNovelMap=new Char2FloatOpenHashMap());        
        //returnMap is not null when this this kmer was already associated to a Char2FloatOpenHashMap 
        if (existingMap!=null) { 
            float returnVal = existingMap.putIfAbsent((char)nodeId, PPStar);
            if (PPStar > returnVal) {
                existingMap.put((char)nodeId, PPStar);
            }
        } else {
            //set the default value return by puIfAbsent to 10.0
            //which is a value >0 (PP* cannot be >0)
            preparedNovelMap.defaultReturnValue(10f);
            preparedNovelMap.put((char)nodeId, PPStar);
            hash.put(word, preparedNovelMap);
        }
        
            
            
//        if (!hash.containsKey(word)) {            
//            UnionPointerWithMap hp=new UnionPointerWithMap();
//            hp.registerTuple(nodeId, refPos, PPStar);
//            hash.put(word, hp);
//        } else {
//            //update list of pairs
//            hash.get(word).registerTuple(nodeId, refPos, PPStar);
//        }
        
    }

    /**
     * for debug purposes only (access to CustomHashMap hashtable)
     * @return 
     */
    public Object2ObjectOpenCustomHashMap<byte[],Char2FloatOpenHashMap> getHash() {
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
        return null;
    }    
    
    /**
     * retrieves only (nodeId,PP*) pairs stored under the position
     * associated to the best PP*
     * @param w
     * @return null if word not in present in hash
     */
    @Deprecated
    public List<Pair> getPairsOfTopPosition(byte[] w) {
        return null;
    }  
    
    /**
     * 
     * @param w
     * @return 
     */
    public Char2FloatMap.FastEntrySet getPairsOfTopPosition2(byte[] w) {
        Char2FloatOpenHashMap map=null;
        if ((map=hash.get(w))!=null) {
            return map.char2FloatEntrySet();
        } else {
            return null;
        }
    }
//    public List<Pair> getPairsOfTopPosition2(byte[] w) {
//        UnionPointer up=null;
//        if ((up=hash.get(w))!=null) {
//            return up.getPairList(0);
//        } else {
//            return null;
//        }
//    }
    
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return 
     */
    @Deprecated
    public int[] getPositions(byte[] w) {
        return null;
    }
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return -1 if word not in hash
     */
    @Deprecated
    public int getTopPosition(byte[] w) {
        return -1;
    }
    
    
    @Deprecated
    public void sortData() {
        //trim hash to max size
        this.hash.trim(maxCapacitySize);
        //do sorting
        double startTime=System.currentTimeMillis();
        AtomicLong pairCount=new AtomicLong(0);
        hash.values().stream().forEach( l -> {pairCount.addAndGet(l.size());} );
        double endTime=System.currentTimeMillis();
        Infos.println("# Pairs sorted in top positions: "+pairCount.get());
        Infos.println("Pair sorting took "+(endTime-startTime)+" ms");
    }
    
    

    
    /**
     * pairs selected by associated position
     * @param w
     * @param position
     * @return 
     */
    @Deprecated
    public List<Pair> getPairs(byte[] w,int position) {
        return null;
    }    
    
    public Set<byte[]> keySet() {
        return hash.keySet();
    }
    
    /**
     * empty all the positions which where not associated to the best PP*
     */
    @Deprecated
    public void reduceToMediumHash() {
        return;
    }

    /**
     * in top position, empties all the positions which where not associated to the best PP*,
     * and retains only X pairs at the best position
     * @param X
     */
    @Deprecated
    public void reducetoSmallHash(int X) {
        return;
    }
    
    /**
     * discards all words where top positions is associated to more than X nodes
     * @param X
     */
    @Deprecated
    public void reducetoSmallHash_v2(int X) {
        return;
    }    
    
    
    ////////////////////////////////////////////////////////////////////////////

    

    

    

    
    
    
    
}
