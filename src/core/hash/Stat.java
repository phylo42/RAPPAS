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

/**
 * A (PPStar, refPositon) association
 * @author yann
 */
public class Stat {
    
    //as java as no unsigned int, we simply trick it by using char,
    //this allows a max of  2^16-1 node ids, using 16bit instead of the 32bits of integer
    private float PPStar=Float.NEGATIVE_INFINITY;
    private char refPosition='\u0000';
    
    public Stat(float PPStar, int refPosition) {
        this.PPStar=PPStar;
        this.refPosition=(char)refPosition;
    }
    
    public float getPPStar() {
        return PPStar;
    }

    public int getRefPosition() {
        return refPosition;
    }

    public void setPPStar(float PPStar) {
        this.PPStar = PPStar;
    }
    
    public void setRefPositon(int refPostion) {
        this.refPosition = (char)refPosition;
    }

    @Override
    public String toString() {
        return "PPStar="+PPStar+" refPostion="+refPosition;
    }
}
