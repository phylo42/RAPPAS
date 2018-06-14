/*
 * Copyright (C) 2018 yann
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package core.hash;

import core.States;
import etc.Infos;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenCustomHashMap;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * In this hash version, the Nodes composing the buckets contain 
 * a list of Triplets representing a
 * (nodeId, PP*,refPosition)
 * @author yann
 */
public class CustomHash_Triplet implements Serializable {
    
    public static final int NODES_POSITION=1;
    public static final int NODES_UNION=2;
    
    int nodeType=NODES_UNION;
    
    Object2ObjectOpenCustomHashMap<byte[],Triplet_16_32_32_bit> hash;
    
    ArrayList<Triplet> tripletsBuffer =new ArrayList(1024);
    
    //TripletList list=new TripletList(0);
    
    int maxCapacitySize=-1;

    /**
     *
     * @param k
     * @param s
     */
    public CustomHash_Triplet(int k, States s) {        
        this.maxCapacitySize=new Double(Math.pow(s.getNonAmbiguousStatesCount(), k)).intValue();
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
        Triplet_16_32_32_bit t = null;
        if (!hash.containsKey(word)) {
            t = new Triplet_16_32_32_bit(nodeId, PPStar, refPos);
            t.registerTuple(nodeId, refPos, PPStar);
            hash.put(word, t);
        } else {
            //update list of triplets
            hash.get(word).registerTuple(nodeId, refPos, PPStar);
        }
        
    }

    /**
     * for debug purposes only (access to CustomHashMap hashtable)
     * @return 
     */
    public Object2ObjectOpenCustomHashMap<byte[],Triplet_16_32_32_bit> getHash() {
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
     * best (nodeId,PP*,refPos) associated to this word
     * @param w
     * @return 
     */
    public Triplet getTopTriplet(byte[] w) {
        Triplet tp=null;
        if ((tp=hash.get(w))!=null) {
            hash.trim();
            return tp.getBestTriplet();
        } else {
            return null;
        }
    }
    
    public int getTopPosition(byte[] w) {
        Triplet tp=null;
        if ((tp=hash.get(w))!=null) {
            return tp.getBestPosition();
        } else {
            return -1;
        }
    }
    
    public List<Triplet> getTriplets(byte[] w) {
        tripletsBuffer.clear();
        tripletsBuffer.addAll(hash.get(w).getTripletList(w));
        return tripletsBuffer;
    }
    
    public List<Triplet> getBestTriplets(byte[] w) {
        tripletsBuffer.clear();
        for (Triplet p: hash.get(w).getTripletList(w)) {
            if (p.getPPStar()>hash.get(w).getPPStar()) {
                tripletsBuffer.add(p);
            }           
        }
        return tripletsBuffer;
    }
    
    public void sort() {
        //No need for sorting in Union DB
        //but will use this step to trim the hash
        hash.trim();
    }
    
    public void sortData() {
        //trim hash to max size
        this.hash.trim(maxCapacitySize);
        //do sorting
        double startTime=System.currentTimeMillis();
        AtomicInteger tripletCount=new AtomicInteger(0);
        hash.values().stream().forEach( l -> {this.sort();tripletCount.addAndGet(hash.size());} );
        double endTime=System.currentTimeMillis();
        Infos.println("# Triplets sorted in top positions: "+tripletCount.get());
        Infos.println("Triplet sorting took "+(endTime-startTime)+" ms");
    }
    
    public Set<byte[]> keySet() {
        return hash.keySet();
    }
    
    /**
     * reference alignment positions associated to a word
     * @param w
     * @return 
     */
    public int[] getPositions(byte[] w) {
        Triplet_16_32_32_bit tp=null;
        if ((tp=hash.get(w))!=null) {
            return tp.getPositions();
        } else {
            return null;
        }
    }

}
