/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.older;

import core.older.Colmer;
import alignement.Alignment;
import core.Locality;
import core.SimpleWord;
import core.States;
import core.Word;
import static core.algos.SequenceKnife.SAMPLING_LINEAR;
import static core.algos.SequenceKnife.SAMPLING_SEQUENTIAL;
import static core.algos.SequenceKnife.SAMPLING_STOCHASTIC;
import etc.Infos;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import javax.swing.tree.TreeNode;
import tree.PhyloNode;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class ColmerSet implements Serializable {
    
    private static final long serialVersionUID = 4000L;

    //done through the shuffling of the table merOrder
    final public static int SAMPLING_STOCHASTIC=1; 
    //ie, [0,k],[1,k+1],[2,k+2], ... [n,k+n]
    final public static int SAMPLING_LINEAR=2;
    //ie, [0,k],[k+1,2k],[2k+1,3k], ... ,[1,k+1],[k+2,2K+1], ... ,[2,k+2],[k+3,2k+2]
    final public static int SAMPLING_SEQUENTIAL=3;

    private int k;
    private int minK;
    private States states=null;
    private ConcurrentHashMap<Word,ConcurrentHashMap<Integer,Locality>> register=null; //map(word)=(map(nodeId)=Locality,i.e colmerId,PP*))
    private Alignment align=null;
    private PhyloTree tree=null;
    ArrayList<Colmer> allColmers=null;
    
    //for the mer sampling
    private Long seed=null;
    private int iterator=0; //last returned mer, as index of the merOrder table
    private int[] merOrder=null; //to define the order in which the mer are returned

    private int colmerCount=-1;
    private int currentFill=0;//register how many positions in the alignments have been assigned to colmers.
    
    //to check if an added colmer respects the last autorized size.
    private boolean full=false;
    private int lastPosition=-1;
    private int lastAllowedK=-1;
    
    /**
     * 
     * @param align
     * @param tree
     * @param states
     * @param k
     * @param samplingMode
     */
    public ColmerSet(Alignment align, PhyloTree tree, States states, int k, int minK, int samplingMode) {
        this.k=k;
        this.align=align;
        this.tree=tree;
        this.states=states;
        this.allColmers=new ArrayList<>(1000);
        this.register=new ConcurrentHashMap<>();
        merOrder=new int[align.getLength()];
        Arrays.fill(merOrder, 0);
        switch (samplingMode) {
            case SAMPLING_LINEAR:
                for (int i = 0; i < merOrder.length; i++) {
                    merOrder[i]=i;
                }
                break;
            case SAMPLING_STOCHASTIC:
                shuffledMerOrder();
                break;
            case SAMPLING_SEQUENTIAL:
                sequencialMerOrder();
                break;
            default:
                Infos.println("Sampling mode not recognized, use one of SAMPLING_LINEAR, SAMPLING_STOCHASTIC, SAMPLING_SEQUENTIAL");
                break;
        }
    }
     
    /**
     * must be called before adding a colmer with @addColmer, ensure colmers k consistency with the alignment length\n
     * will return null when no more Colmer can be formed 
     * @param k
     * @param minK
     * @return 
     */
    public Colmer getNextColmer() {
        
        if (iterator>merOrder.length-1) {
            return null;
        }
        int currentPosition=merOrder[iterator];
        int charactersLeft=merOrder.length-currentPosition;
        if (charactersLeft>=minK) {
            Colmer c=null;
            if (charactersLeft<k) {
                c=new Colmer(++colmerCount,charactersLeft, currentPosition, tree.getNodeCount(), states.getStateCount());
            } else {
                c=new Colmer(++colmerCount,k, currentPosition, tree.getNodeCount(), states.getStateCount());
            }
            iterator++;
            return c;
            
        } else {
            //this allow to skipped words that are too short but in the middle
            //of the mer ordering, we just skip them and go to the next one.
            iterator++;
            getNextColmer();
        }
        return null;
        
        
//        if (lastPosition+k>=align.getLength() && !full) { //last colmer smaller than k
//            int sitesLeft=align.getLength()-1-lastPosition;
//            if(sitesLeft<minColmerLength) {return null;}
//            Colmer c=new Colmer(++colmerCount,sitesLeft, lastPosition+1, tree.getNodeCount(), states.getStateCount());
//            lastAllowedK=sitesLeft;
//            lastPosition=lastPosition+lastAllowedK;
//            full=true;
//            return c;
//        } else if (!full){
//            Colmer c=new Colmer(++colmerCount,k, lastPosition+1, tree.getNodeCount(), states.getStateCount());
//            lastAllowedK=k;
//            lastPosition=lastPosition+lastAllowedK;
//            return c;
//        } else {        
//            return null;
//        }
    }
    
    /**
     * a table representing the order in which mers are returned \n
     * each value is the 1st position of the mer which will sample the \n
     * sequence from this position to position+k (excepted if < to the  min k)
     * @return 
     */
    public int[] getMerOrder() {
       return merOrder; 
    }
    
    /**
     * add a colmer to the set, note that it had to be return by @getNextColmer
     * @param c 
     */
    public void addColmer(Colmer c) {
        //check if last colmer of the alignment.
        if (!full) {
            allColmers.add(c);
        }
    }
    
    public Colmer getColmerBySite(int site) {
        System.out.println("Not supported yet");
        return null;
    }
    
    public Colmer getColmerById(int colmerId) {
        return allColmers.get(colmerId);
    }
    
    public int getSetSize() {
        return register.size();
    }
    
    public int getColmerCount() {
        return colmerCount;
    }

    /**
     * returns the structure used for querying the ancestral states.
     * i.e. map(word)=(map(nodeId)=Locality,i.e colmerId,PP*))
     * @return 
     */
    public ConcurrentHashMap<Word,ConcurrentHashMap<Integer,Locality>> getRegister() {
        return register;
    }
    
    
    /**
     * get the full status if all sites of the original alignment have been associated to colmers
     * @return 
     */
    public boolean isFull() {
        return full;
    }

    
    public void registerWord(Word w) {
        if (!register.containsKey(w)) {
            register.put(new SimpleWord(w), new ConcurrentHashMap<>(new Double(tree.getNodeCount()*0.33).intValue()));
        }
    }
    
    public void associateLocality(Word w, PhyloNode n,Colmer c, double proba) {
        if (!register.get(w).containsKey(n.getId())) {
            register.get(w).put(n.getId(),new Locality(this));
        }
        register.get(w).get(n.getId()).addTuple(c.getId(), proba);
    }
    
    
    /**
     * must be called after instantiation if one is interested to retrieve \n
     * the same shuffled mer order
     * @param seed 
     */
    public void forceSeed(long seed) {
        this.seed=seed;
        shuffledMerOrder();
    }
    
    private void shuffledMerOrder() {
        for (int i = 0; i < merOrder.length; i++) {
           merOrder[i]=i;
        }
        Random generator = null;
        if (seed!=null) {
            generator=new Random(seed);
        } else {
            generator=new Random(System.nanoTime());
        }
        for (int i = 0; i < merOrder.length - 1; i++) {
          int j = i + generator.nextInt(merOrder.length - i);
          int t = merOrder[j];
          merOrder[j] = merOrder[i];
          merOrder[i] = t;
        }
        
    }
    
    private void sequencialMerOrder() {
        //ie, [0,k],[k,2k],[2k,3k], ... ,[1,k+1],[k+1,2K+1], ... ,[2,k+2],[k+2,2k+2]
        int counter=0;
        int shift=0; //consumed on the left
        for (int i = 0; i < merOrder.length ; i++) {
            if (shift==k) {break;}
            for (int j = 0; j < (merOrder.length/k)+1; j++) { //+1 to get the last incomplete mer (length<k)
                if ((shift+j*k)<merOrder.length) {
                    merOrder[counter]=shift+j*k;
                    counter++;
                }
            }
            shift++;
        }
    }

    
}
