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

import core.DNAStates;
import core.States;
import java.util.Arrays;

/**
 *
 * @author yann
 */
public class TEST_CustomHash_Triplet {
    
    //to test the class
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        States s=new DNAStates();
        
        //ArrayList<Triplet> triplet=new ArrayList<>();
        
//        Triplet t1 = new Triplet_16_32_32_bit(1, 0.00004f, 16);
//        Triplet t2 = new Triplet_16_32_32_bit(2, 0.000000007f, 5);
//        Triplet t3 = new Triplet_16_32_32_bit(3, 0.18f, 25);
        
//        triplet.add(t1);
//        triplet.add(t2);
//        triplet.add(t3);
        
        CustomHash_Triplet ht = new CustomHash_Triplet(8, s);
        
        byte[] word={1, 3, 3, 0, 2, 0, 0, 0}; //byte[] cword=s.compressMer(word);
        byte[] word2={1, 0, 0, 0, 3, 1, 2, 0};
        byte[] word3={1, 0, 0, 0, 3, 1, 2, 1};
        byte[] word4={1, 1, 3, 0, 1, 3, 0, 1};
        byte[] word5={1, 0, 0, 0, 3, 1, 2, 0};
        
//        System.out.println("word:"+word.length);
//        System.out.println("cword:"+cword.length);
        
        
        //Triplet t4 = new Triplet_16_32_32_bit();
        
        //t4.registerTuple(4, 1800, 0.0000000038f);
        //System.out.println("t4: "+t4.getTripletList(word));
        
        //ht.hash.put(word, );
        //triplet.add(t4);
        
        ht.addTuple(word, 0.00004f, 1, 16);
        ht.addTuple(word2, 0.00018f, 2, 25);
        ht.addTuple(word3, 0.000000007f, 2, 5);
        ht.addTuple(word4, 0.8f, 100, 250);
        ht.addTuple(word5, 0.0005f, 3, 20);
        ht.addTuple(word2, 0.001817f, 3, 25);
        ht.addTuple(word2, 0.00015f, 3, 200);
        ht.addTuple(word2, 0.00000798f, 3, 25);
        ht.addTuple(word2, 0.0008f, 3, 5);
        ht.addTuple(word2, 0.0000005f, 7, 25);
        ht.addTuple(word4, 0.00000798f, 3, 39);
        //ht.tripletsBuffer.addAll(triplet);
       
        //System.out.println("word2:"+ht.getTriplets(word2)+" --> "+ht.getTriplets(word2).size());
        // word2:[nodeId=2 PPStar=1.8E-4 refPostion=25, nodeId=3 PPStar=0.001817 refPostion=20, nodeId=7 PPStar=5.0E-7 refPostion=45] --> 3
        //System.out.println("hash"+ht.getHash());
        //System.out.println("tripletsList"+ht.hash.size());
        
        //ht.getHash().entrySet().stream().forEach(e->{System.out.println(Arrays.toString(e.getKey())+" --> "+e.getValue().getNode().length+" --> "+e.getValue().list);});
        
        /*
        [1, 3, 3, 0, 2, 0, 0, 0] --> 1
        [1, 0, 0, 0, 3, 1, 2, 0] --> 2
        [1, 1, 3, 0, 1, 3, 0, 1] --> 1
        [1, 0, 0, 0, 3, 1, 2, 1] --> 1    
        */
        
        System.out.println("test: "+ht.getTuples(word2));
    }    
}
