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

import java.util.ArrayList;

/**
 * A (nodeId,PP*,refPostion) association
 * @author yann
 */
public interface Triplet extends Comparable<Triplet> {
    public int getNodeId();
    public float getPPStar();
    public int getRefPosition();
    public void setNodeId(int nodeId);
    public void setPPStar(float PPStar);
    public void setRefPosition(int refPostition);
    public void registerTuple(int nodeId, int refPosition, float PPStar);
    public ArrayList<Triplet> getTripletList(byte[] w);
    //public ArrayList<Triplet> getTripletList(int refPosition);
    public Triplet getBestTriplet();
    public int[] getNode();
    public int[] getPositions();
    public int getBestPosition();
    
    
    
    @Override
    public default int compareTo(Triplet o) {
        return -Float.compare(this.getPPStar(),o.getPPStar());
    }
}
