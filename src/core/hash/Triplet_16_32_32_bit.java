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

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Pair_32_32_32_bit representing a (nodeId,PP*,refPosition) association as (unsigned int,float,int).
 * These are used in the LinkedLists
 * @author yann
 */
public class Triplet_16_32_32_bit implements Triplet,Serializable {
    
    private static final long serialVersionUID = 7102L;

    private char nodeId='\u0000'; //init to 0
    private float PPStar=Float.NEGATIVE_INFINITY;
    private int refPosition=0; //init to 0
    
    TripletList list=new TripletList(0);

    public Triplet_16_32_32_bit(int nodeId, float PPStar, int refPosition) {
        this.nodeId=(char)nodeId;
        this.PPStar=PPStar;
        this.refPosition=refPosition;
    }

    @Override
    public int getNodeId() {
        return (int)nodeId;
    }

    @Override
    public float getPPStar() {
        return PPStar;
    }
    
    @Override
    public int getRefPosition() {
        return refPosition;
    }

    @Override
    public void setNodeId(int nodeId) {
        this.nodeId = (char)nodeId;
    }

    @Override
    public void setPPStar(float PPStar) {
        this.PPStar = PPStar;
    }
    
    @Override
    public void setRefPosition(int refPosition) {
        this.refPosition = refPosition;
    }
    
    @Override
    public void registerTuple(int nodeId, int refPosition, float PPStar) {
        boolean found=false;
        for (Triplet t:list) {
            //nodeid already registered
            //if (t.getNodeId()==nodeId && t.getRefPosition()==refPosition) {
            if (t.getNodeId()==nodeId) {
                found=true;
                //replace previous triplet if better PP*
                if (PPStar > t.getPPStar()) {
                    t.setPPStar(PPStar);
                    // replace position associated to better PP*
                    t.setRefPosition(refPosition);
                }
                break;
            }
        }
        //node not yet registered, do it
        if (!found) {
            list.add(new Triplet_16_32_32_bit(nodeId, PPStar, refPosition));
        }

    }
    
    @Override
    public ArrayList<Triplet> getTripletList(byte[] w) {
        return list;
    }
    
    @Override
    public int[] getNode() {
        int[] node=new int[list.size()];
        for (int i=0;i<list.size();i++) {
            node[i]=list.get(i).getNodeId();
        }
        return node;
    }
    
    @Override
    public int[] getPositions() {
        int[] pos=new int[list.size()];
        for (int i=0;i<list.size();i++) {
            pos[i]=list.get(i).getRefPosition();
        }
        return pos;
    }
    
    /**
     * get best reference position.
     * @return 
     */
    @Override
    public int getBestPosition() {
        return list.get(0).getRefPosition();
    }
    
    @Override
    public Triplet getBestTriplet() {
        return list.get(0);
    }

    /**
     * the comparator is inversed to put highest values first
     * @param o
     * @return 
     */
    public int compareTo(Triplet_16_32_32_bit o) {
        return -Float.compare(this.PPStar,o.getPPStar());
    }

    @Override
    public String toString() {
        return "nodeId="+(int)nodeId+" PPStar="+PPStar+" refPosition="+(int)refPosition;
    }
    
//    private class TripletList extends ArrayList<Triplet> extends Comparable<TripletList>,Serializable {
//        private static final long serialVersionUID = 7210L;
//        
//        public TripletList() {
//            super(1);
//        }
//    }
    
    // mettre dans une classe 
    private class TripletList extends ArrayList<Triplet> implements Comparable<Triplet>, Serializable {
        private int refPosition=-1;
        
        public TripletList(int refPosition) {
            super(1);
            this.refPosition=refPosition;
            
        }
        
        public int getRefPosition() {
            return refPosition;
        }

        @Override
        public int compareTo(Triplet o) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }
  
}
