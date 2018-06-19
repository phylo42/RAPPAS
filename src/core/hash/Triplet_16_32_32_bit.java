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
    

    public Triplet_16_32_32_bit(int nodeId, float PPStar, int refPosition) {
        if (nodeId>='\uFFFF') {
            System.out.println("Extended tree node number reached unsigned short limit (2^16-1). RAPPAS currently not designed to handle such big trees.");
            System.exit(1);
        }
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
    

  
}
