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
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import java.io.Serializable;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * In this hash version, the Nodes composing the buckets contain 
 * a list of Triplets representing a
 * (nodeId, PP*,refPosition)
 * @author yann
 */
public class CustomHash_Triplet implements Serializable,CustomHash {
    
    public static final int nodeType=CustomHash.NODES_TRIPLET;
    
    Object2ObjectOpenCustomHashMap<byte[],ObjectOpenHashSet<Triplet_16_32_16_bit>> hash;
        
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
    @Override
    public void addTuple(byte[] word, float PPStar,int nodeId,int refPos) {
        ObjectOpenHashSet<Triplet_16_32_16_bit> set = hash.get(word);
        if (set==null) {
            set =new ObjectOpenHashSet<>();
            set.add(new Triplet_16_32_16_bit(nodeId, PPStar, refPos));
            hash.put(word, set);  
        } else {
            boolean found=false;
            for (Triplet t:set) {
                //nodeid already registered
                if (t.getNodeId()==nodeId) {
                    found=true;
                    //replace previous triplet if better PP*
                    if (PPStar > t.getPPStar()) {
                        t.setPPStar(PPStar);
                        // replace position associated to better PP*
                        t.setRefPosition(refPos);
                    }
                    break;
                }
            }
            //node not yet registered, do it
            if (!found) {
                set.add(new Triplet_16_32_16_bit(nodeId, PPStar, refPos));
            }
        }
        
    }

    /**
     * for debug purposes only (access to CustomHashMap hashtable)
     * @return 
     */
    @Deprecated
    public Object2ObjectOpenCustomHashMap<byte[],ObjectOpenHashSet<Triplet_16_32_16_bit>> getHash() {
        return hash;

    }
    
    /**
     * to know which type of nodes is backing this hash, 
     * one of CustomHash.NODES_XXX
     * @return 
     */
    @Override
    public int getHashType() {
        return CustomHash_Triplet.nodeType;
    }

    
    @Override
    public ObjectOpenHashSet<Triplet_16_32_16_bit> getTuples(byte[] w) {
        return hash.get(w);
    }
    
 
    
    public void sort() {
        //No need for sorting
        //but will use this step to trim the hash
        hash.trim();
    }
    
    @Override
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
    
    

}
