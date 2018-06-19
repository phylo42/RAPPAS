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

import core.Word;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Triplet_16_32_16_bit representing a (nodeId,PP*,refPosition) association as (unsigned int, float, unsigned int) from a size of alignment less than 32,000kb.
 * These are used in the LinkedLists
 * @author yann
 */
public class Triplet_16_32_16_bit implements Triplet,Serializable {
    
    private static final long serialVersionUID = 7101L;

    //as java as no unsigned int, we simply trick it by using char,
    //this allows a max of  2^16-1 node ids, using 16bit instead of the 32bits of integer
    private char nodeId='\u0000'; //init to 0
    private float PPStar=Float.NEGATIVE_INFINITY;
    private char refPosition='\u0000'; //init to 0
        
    public Triplet_16_32_16_bit(int nodeId, float PPStar, int refPosition) {
        if (nodeId>='\uFFFF') {
            System.out.println("Extended tree node number reached unsigned short limit (2^16-1). RAPPAS currently not designed to handle such big trees.");
            System.exit(1);
        }
        if (refPosition>='\uFFFF') {
            System.out.println("Extended reference position number reached unsigned short limit (2^16-1). Triplet_16_32_16_bit class not designed to handle such big sequences.");
            System.exit(1);
        }
        this.nodeId=(char)nodeId;
        this.PPStar=PPStar;
        this.refPosition=(char)refPosition;
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
        return (int)refPosition;
    }

    @Override
    public void setNodeId(int nodeId) {
        if (nodeId>='\uFFFF') {
            System.out.println("Extended tree node number reached unsigned short limit (2^16-1). RAPPAS currently not designed to handle such big trees.");
            System.exit(1);
        }
        this.nodeId = (char)nodeId;
    }

    @Override
    public void setPPStar(float PPStar) {
        this.PPStar = PPStar;
    }
    
    @Override
    public void setRefPosition(int refPosition) {
        this.refPosition = (char)refPosition;
    }
    

    @Override
    public String toString() {
        return "nodeId="+(int)nodeId+" PPStar="+PPStar+" refPostion="+(int)refPosition;
    }
   
}
