/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import core.older.ColmerSet;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.DoubleStream;

/**
 *
 * @author ben
 */
public class Locality implements Serializable {

    private static final long serialVersionUID = 4001L;

    private ArrayList<Tuple> tuples=null;
    private ColmerSet cs=null;

    public Locality(ColmerSet cs) {
        //arbitray, the initial list size will be 1/20th of the posible colmers
        this.cs=cs;
        tuples=new ArrayList<>(cs.getColmerCount()/20);
    }

    public ArrayList<Tuple> getTuples() {
        return tuples;
    }
    
    public void addTuple(int colmerId,double ppStar) {
        tuples.add(new Tuple(colmerId,ppStar));
    }
    
    
    @Override
    public String toString() {
        return "Locality[tuples="+tuples.size()+", top5="+Arrays.toString(getBestLocalities(5))+", top5proba="+Arrays.toString(getBestProbas(5))+"]";
    }

    
    
    public int[] getBestLocalities(int limit) {
        Collections.sort(tuples);
        return Arrays.copyOfRange(tuples.stream().mapToInt(t->t.colmerId).toArray(),0,limit);
    }
    
    public double[] getBestProbas(int limit) {
        Collections.sort(tuples);
        return Arrays.copyOfRange(tuples.stream().mapToDouble(t->t.ppStar).toArray(),0,limit);
    }
    
    public DoubleStream getAllProbas() {
        return tuples.stream().mapToDouble(t->t.ppStar);
    }
    
    public List<Tuple> getTuplesSortedByProbas() {
        Collections.sort(tuples);
        return tuples;
    }
    
    
    /**
     * simple structure to associate colmer and PP*\n
     * Is comparable, but return order opposite to natural order to sort such as [5,4,3,2,1]
     */
    public class Tuple implements Comparable<Tuple>,Serializable {
        
        private static final long serialVersionUID = 4011L;

        int colmerId=-1;
        double ppStar=-1.0;

        public Tuple(int colmerId,double ppStar) {
            this.colmerId=colmerId;
            this.ppStar=ppStar;
        }
        
        @Override
        public int compareTo(Tuple o) {
            if (this.ppStar<o.ppStar) {
                return 1;
            } else if (this.ppStar>o.ppStar) {
                return -1;
            } else {
                return 0;
            }
        }

        public int getColmerId() {
            return colmerId;
        }

        public double getPpStar() {
            return ppStar;
        }

        @Override
        public String toString() {
            return "T="+colmerId+":"+ppStar;
        }
        
        
        
    }
    
}
