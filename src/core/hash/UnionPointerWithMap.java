/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import it.unimi.dsi.fastutil.chars.Char2FloatMap;
import it.unimi.dsi.fastutil.chars.Char2FloatOpenHashMap;
import java.io.Serializable;
import java.util.List;
import java.util.stream.Collectors;


/**
 * Intermediate table used to select the (nodeId,PP*) pairs through 
 * associated reference position.
 * @author ben
 */
public class UnionPointerWithMap implements HashPointer,Serializable {
    
    private static final long serialVersionUID = 7300L;

    
    //as java as no unsigned int, we simply trick it by using char,
    //this allows a max of  2^16-1 node ids, using 16bit keys instead of integer 32bits
    Char2FloatOpenHashMap map=null;

    /**
     *
     */
    public UnionPointerWithMap() {
        //only one fake position map, with fake position set to -10
        map=new Char2FloatOpenHashMap(20, 0.8f);
        //use positive value (impossible log10(PP*)) as default value
        map.defaultReturnValue(10);
    }

    @Override
    public void registerTuple(int nodeId, int refPosition, float PPStar) {
        float result = map.putIfAbsent((char)nodeId, PPStar);
        if (PPStar > result) {
            map.put((char)nodeId, PPStar);
        }
    }
    
    /**
     * ref alignment positions associated to this node, order by underlying max(PP*)
     * @return 
     */
    @Override
    public int[] getPositions() {
        int[] t={0};
        return t;
    }

    /**
     * get all Pairs associated to a particular reference position.
     * @param refPosition
     * @return 
     */
    @Override
    @Deprecated
    public List<Pair> getPairList(int refPosition) {
        return map.char2FloatEntrySet().stream().map((e)->new Pair_16_32_bit((int)e.getCharKey(),e.getFloatValue())).collect(Collectors.toList());
    }
    
    public Char2FloatMap.FastEntrySet getPairs() {
        return map.char2FloatEntrySet();
    }
    
    
    
    
    /**
     * get number of nodes associated to best position
     * @return 
     */
    @Override
    public int getPairCountInTopPosition() {
        
        return map.size();
    }

    /**
     * get best (nodeId;PP*).
     * @return 
     */
    @Override
    public Pair getBestPair() {
        return map.char2FloatEntrySet().stream().findFirst().map((e)->new Pair_16_32_bit((int)e.getCharKey(),e.getFloatValue())).get();
    }

    /**
     * get best reference position.
     * @return 
     */
    @Override
    public int getBestPosition() {
        return -1;
    }

    @Override
    public void sort() {
        //No need for sorting in Union DB
        //but will use this step to trim the hash
        map.trim();
    }  



    
    
}
