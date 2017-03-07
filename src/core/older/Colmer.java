/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.older;

import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import tree.PhyloNode;
import tree.PhyloTree;

/**
 * structure representing the words that are probable at different trees levels for a specific mer
 * @author ben
 */
public class Colmer implements Comparable<Colmer>,Serializable {
    
    private static final long serialVersionUID = 4010L;
    
    private int id=-1;
    private int k=-1; //colmer length
    private int startSite=-1; //site of the orginial alignment from which started the colmer
    private int nodeCount=-1; //number of nodes in the associated tree

    
    /**
     * 
     * @param id
     * @param k colmer size
     * @param startSite site position on the original full alignment (colmer 1st position)
     * @param nodeCount # nodes in the reference tree
     * @param states # states (4=DNA,20=AA...)
     */
    public Colmer(int id,int k, int startSite, int nodeCount, int states) {
        this.id=id;
        this.k=k;
        this.startSite=startSite;
        this.nodeCount=nodeCount;

    }
    
    public int getStartSite() {
        return startSite;
    }
    
    public int getEndSite() {
        return startSite+k-1;
    }

    public int getId() {
        return id;
    }
    
    @Override
    public int compareTo(Colmer o) {
        if (this.id>o.getId()) {
            return 1;
        } else if (this.id<o.getId()) {
            return -1;
        } else {
            return 0;
        }
    }
    
}
