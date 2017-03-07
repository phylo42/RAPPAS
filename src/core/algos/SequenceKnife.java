/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import core.DNAStates;
import core.QueryWord;
import core.States;
import core.Word;
import etc.Infos;
import inputs.Fasta;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * to easily derive mers from a query sequence,\n
 * the class in intended to retrieve mer in a ordered or stochastic way\n
 * this should allow to abandon the read comparison after a certain number\n
 * of matches and mismatches.
 * @author ben
 */
public class SequenceKnife {
    
    /**
     * done through the shuffling of the table merOrder
     */
    final public static int SAMPLING_STOCHASTIC=1; 
    /**
     * ie, [0,k],[1,k+1],[2,k+2], ... [n,k+n]
     */
    final public static int SAMPLING_LINEAR=2;
    /**
     * ie, [0,k],[k,2k],[2k,3k], ... ,[1,k+1],[k+1,2K+1], ... ,[2,k+2],[k+2,2k+2]
     */
    final public static int SAMPLING_SEQUENTIAL=3;
    /**
     * ie, [0,k],[k+1,2k],[2k+1,3k],...
     */
    final public static int SAMPLING_NON_OVERLAPPING=4;
    
    private Long seed=null;
    
    private int k=-1;
    private int minK=-1;
    private int iterator=0; //last returned mer, as index of the merOrder table
    private byte[] sequence=null; //the inital sequence itself
    private int[] merOrder=null; //to define the order in which the mer are returned
    private States s=null;
    private int step=-1;
    
    /**
     * Basic constructor, will return mers in linear order
     * @param f
     * @param k
     * @param minK
     * @param s 
     */
    public SequenceKnife(Fasta f, int k, int minK, States s) {
        this.k=k;
        this.minK=minK;
        this.s=s;
        String seq = f.getSequence();
        initTables(seq, SAMPLING_LINEAR);
    }
   
    /**
     * Basic constructor, will return mers in linear order
     * @param seq
     * @param k
     * @param minK
     * @param s 
     */
    public SequenceKnife(String seq, int k, int minK, States s) {
        this.k=k;
        this.minK=minK;
        this.s=s;
        initTables(seq, SAMPLING_LINEAR);
    }
    
    /**
     * constructor setting the mer order through SAMPLING_* static variables
     * @param seq
     * @param k
     * @param minK
     * @param s
     * @param samplingMode 
     */
    public SequenceKnife(String seq, int k, int minK, States s, int samplingMode) {
        this.k=k;
        this.minK=minK;
        this.s=s;
        initTables(seq, samplingMode);
    }
    
    /**
     * constructor setting the mer order through SAMPLING_* static variables
     * @param seq
     * @param k
     * @param minK
     * @param s
     * @param samplingMode 
     */
    public SequenceKnife(Fasta f, int k, int minK, States s, int samplingMode) {
        this.k=k;
        this.minK=minK;
        this.s=s;
        String seq = f.getSequence();
        initTables(seq, samplingMode);
    }
    
    private void initTables(String seq, int samplingMode) {
        sequence=new byte[seq.length()];
        merOrder=new int[seq.length()];
        for (int i = 0; i < seq.length(); i++) {
            sequence[i]=s.stateToByte(seq.charAt(i));
        }
        //Infos.println("Binary seq: "+Arrays.toString(sequence));
        Arrays.fill(merOrder, 0);
        switch (samplingMode) {
            case SAMPLING_LINEAR:
                merOrder=new int[seq.length()];
                for (int i = 0; i < merOrder.length; i++) {
                    merOrder[i]=i;
                }
                this.step=1;
                break;
            case SAMPLING_NON_OVERLAPPING:
                merOrder=new int[(seq.length()/k)+1];
                for (int i = 0; i < seq.length(); i++) {
                    if (i%k==0) {
                        merOrder[i/k]=i;
                    }
                }
                this.step=k;
                break;
            case SAMPLING_STOCHASTIC:
                merOrder=new int[seq.length()];
                shuffledMerOrder();
                this.step=1;
                break;
            case SAMPLING_SEQUENTIAL:
                merOrder=new int[seq.length()];
                sequencialMerOrder();
                this.step=1;
                break;
            default:
                Infos.println("Sampling mode not_recognized !");
                break;
        }
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
     * must be called to retireve mers one by one
     * @return the next mer as a @Word, null is no more mers to return
     */
    public QueryWord getNextWord() {
        if (iterator>merOrder.length-1) {
            return null;
        }
        int currentPosition=merOrder[iterator];
        int charactersLeft=sequence.length-currentPosition;
        if (charactersLeft>=minK) {
            byte[] word=null;
            if (charactersLeft<k) {
                word=Arrays.copyOfRange(sequence, currentPosition, currentPosition+charactersLeft);
            } else {
                word=Arrays.copyOfRange(sequence, currentPosition, currentPosition+k);
            }
            iterator++;
            return new QueryWord(word, currentPosition);
            
        } else {
            //Infos.println("Skip word on position "+currentPosition+": length < minK !");
            //this allow to skip words that are too short but in the middle
            //of the shuffled mer order, we just skip them and go to the next one.
            iterator++;
            return getNextWord();
        }
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
    
    public int getStep() {
        return this.step;
    }
    
    public static void main(String[] args) {
        SequenceKnife knife=new SequenceKnife("ATCGCTGATCGATCGA", 7, 4, new DNAStates(), SequenceKnife.SAMPLING_SEQUENTIAL);
        System.out.println(Arrays.toString(knife.getMerOrder()));
        knife=new SequenceKnife("ATCGCTGATCGATCGA", 7, 4, new DNAStates(), SequenceKnife.SAMPLING_STOCHASTIC);
        knife.forceSeed(12345);
        System.out.println(Arrays.toString(knife.getMerOrder()));
        knife=new SequenceKnife("ATCGCTGATCGATCGA", 7, 4, new DNAStates(), SequenceKnife.SAMPLING_LINEAR);
        System.out.println(Arrays.toString(knife.getMerOrder()));
        
        QueryWord w=null;
        while ((w=knife.getNextWord())!=null) {
            System.out.println(w);
            List<QueryWord> mutatedWords = w.getMutatedWords(new DNAStates());
            for (int i = 0; i < mutatedWords.size(); i++) {
                System.out.println("   "+ mutatedWords.get(i));
                
            }
        }
    }
    
}
