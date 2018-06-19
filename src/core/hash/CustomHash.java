/*
 * Copyright (C) 2018 benclaff
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

import java.util.Set;

/**
 *
 * @author benclaff
 */
public interface CustomHash {
    
    public static final int NODES_POSITION=1;
    public static final int NODES_UNION=2;
    public static final int NODES_TRIPLET=3;
    
    /**
     * only entry point to fill the hash (used at DB generation)
     * @param word
     * @param PPStar
     * @param nodeId
     * @param refPos 
     */
    public void addTuple(byte[] word, float PPStar,int nodeId,int refPos);
    
    /**
     * only entry point to query the hash (used at placement)
     * @param word
     * @return 
     */
    public Set getTuples(byte[] word);

    public void sortData();
    
    public Set<byte[]> keySet();
    
    public int getHashType();

    
}
