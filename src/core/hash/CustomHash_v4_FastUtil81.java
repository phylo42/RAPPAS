/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import core.States;
import etc.Infos;
import it.unimi.dsi.fastutil.chars.Char2FloatMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenCustomHashMap;
import java.io.Serializable;
import java.util.ArrayList;
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
public class CustomHash_v4_FastUtil81 implements Serializable,CustomHash {
    
    private static final long serialVersionUID = 7000L;
    
    
    int nodeType=CustomHash.NODES_UNION;
    
    Object2ObjectOpenCustomHashMap<byte[],UnionPointerWithMap> hash;

    
    ArrayList<Pair> pairsBuffer =new ArrayList(1024);
    
    int maxCapacitySize=-1;

    /**
     *
     * @param k
     * @param s
     * @param nodeType one of NODES_UNION or POSITION_UNION
     */
    public CustomHash_v4_FastUtil81(int k, States s, int nodeType) {        
        this.maxCapacitySize=new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)).intValue();
        this.nodeType=nodeType;
        //internal tests showed that with k<14 we generally get at least 75% of the possible k-mers
        this.hash=new Object2ObjectOpenCustomHashMap<>(
                            new Double(maxCapacitySize/4).intValue(),  //intial capacity
                            0.8f, //inital load factor                
                            new HashStrategy()
                        );
    }
    
    /**
     * only entry point to fill the hash (used at DB generation)
     * @param word
     * @param PPStar
     * @param nodeId
     * @param refPos 
     */
    public void addTuple(byte[] word, float PPStar,int nodeId,int refPos) {
        
        if (!hash.containsKey(word)) {            
            UnionPointerWithMap hp=new UnionPointerWithMap();
            hp.registerTuple(nodeId, refPos, PPStar);
            hash.put(word, hp);
        } else {
            //update list of pairs
            hash.get(word).registerTuple(nodeId, refPos, PPStar);
        }
        
    }

    /**
     * for debug purposes only (access to CustomHashMap hashtable)
     * @return 
     */
    public Object2ObjectOpenCustomHashMap<byte[],UnionPointerWithMap> getHash() {
        return hash;
    }
    
    /**
     * to know which type of nodes is backing this hash, 
     * one of CustomHash.NODES_XXX
     * @return 
     */
    @Override
    public int getHashType() {
        return CustomHash.NODES_UNION;
    }

    /**
     * best (nodeId,PP*) associated to this word
     * @param w
     * @return 
     */
    public Pair getTopPair(byte[] w) {
        HashPointer cn=null;
        if ((cn=hash.get(w))!=null) {
            return cn.getBestPair();
        } else {
            return null;
        }
    }    
    
    /**
     * retrieves only (nodeId,PP*) pairs stored under the position
     * associated to the best PP*
     * @param w
     * @return null if word not in present in hash
     */
    @Deprecated
    public List<Pair> getPairsOfTopPosition(byte[] w) {
        
        HashPointer cn=null;
        if ((cn=hash.get(w))!=null) {
            return cn.getPairList(cn.getBestPosition());
        } else {
            return null;
        }
    }  
    
    @Override
    public Char2FloatMap.FastEntrySet getTuples(byte[] w) {
        UnionPointerWithMap up=null;
        if ((up=hash.get(w))!=null) {
            return hash.get(w).getPairs();
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
        if ((cn=hash.get(w))!=null) {
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
        if ((cn=hash.get(w))!=null) {
            return cn.getBestPosition();
        } else {
            return -1;
        }
    }
    
    
    
    public void sortData() {
        //trim hash to max size
        this.hash.trim(maxCapacitySize);
        //do sorting
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
        pairsBuffer.clear();
        for (int p: hash.get(w).getPositions()) {
            pairsBuffer.addAll(hash.get(w).getPairList(p));
        }
        return pairsBuffer;
    }
    
    /**
     * pairs selected by associated position
     * @param w
     * @param position
     * @return 
     */
    public List<Pair> getPairs(byte[] w,int position) {
        if (hash.containsKey(w)) {
            return hash.get(w).getPairList(position);
        } else {
            return null;
        }
    }    
    
    public Set<byte[]> keySet() {
        return hash.keySet();
    }
    
    /**
     * empty all the positions which where not associated to the best PP*
     */
    public void reduceToMediumHash() {
        return;
    }
    
    /**
     * discards all words where top positions is associated to more than X nodes
     * @param X
     */
    public void reducetoSmallHash_v2(int X) {
        assert X>0;
        List<byte[]> collect = hash.keySet().stream() 
                                            //.peek((w)->System.out.println("REDUCING:"+w))
                                            .filter((w) -> hash.get(w).getPairs().size()>X)
                                            //.peek((w)->System.out.println("TRASHED!:"+w))
                                            .collect(Collectors.toList());
        collect.stream().forEach((w)-> {hash.remove(w);});
        collect=null;
    }    
    
    
    ////////////////////////////////////////////////////////////////////////////

    

    

    

    
    
    
    
}
