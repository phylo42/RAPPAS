/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import core.QueryWord;
import core.States;
import etc.Infos;
import etc.exceptions.NonSupportedStateException;
import inputs.Fasta;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Given a sequence, returns mers as byte[]
 * 
 * !!! THIS NEW VERSION HAS A FUNDAMENTAL CHANGE
 * It doesn't return a byte[] of fixed size (size=k), like SequenceKnife
 * but a byte[] of size: k*n , n being the alternative kmers that can
 * be generated regarding the ambiguities found sequence ;
 * In the placement phase, the number of kmers to treat will be
 * length(byte[])/k
 * @author ben
 */
public class AmbigSequenceKnife implements ISequenceKnife {
    
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
    private int merIterator=0; //last returned mer, as index of the merOrder table
    private byte[] sequence=null; //the inital sequence itself as bytes; value -1 is reserved for ambiguities
    private int[] ambiguityCountPerMer; //register the nmber of ambiguities included in the mer starting at this position
    //relative position of ambiguities in a kmer and alternatives
    //map(sequence_pos)=(map(offset)=[stat_1,...state_n])
    private HashMap<Integer,HashMap<Integer,byte[]>> ambiguityOffsets; 
    private int[] merOrder=null; //to define the order in which the mer are returned
    private States s=null;
    private int step=-1;
    private int samplingMode = -1;
    private int maxAmbigPerMer=0;
    
    //the kmer sent at each getNextByteWord() )call
    byte[] word=null;
        
    public void init(String seq) {
        try {
            initTables(seq, samplingMode);
        } catch (NonSupportedStateException ex) {
            Logger.getLogger(AmbigSequenceKnife.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void init(Fasta f) {
    	init(f.getSequence(false));
    }
    
    /**
     * Basic constructor, will return mers in linear order
     * @param k
     * @param minK
     * @param s 
     * @param samplingMode
     */
    public AmbigSequenceKnife(int k, int minK, States s, int samplingMode) {
        this.k=k;
        this.minK=minK;
        this.s=s;
        this.samplingMode=samplingMode;
        //WARNING: this is currently hard-coded, but should become a parameter
        //this formula allows 2 ambiguity for k>16 in DNA and 1 in amino acids
        this.maxAmbigPerMer=(int)Math.floor(Math.pow(k, 1.0/s.getNonAmbiguousStatesCount()));
    }
    
    private void initTables(String seq, int samplingMode) throws NonSupportedStateException {
        merIterator=0;
        sequence=new byte[seq.length()];
        ambiguityCountPerMer=new int[sequence.length];
        ambiguityOffsets=new HashMap<>();
        for (int i = 0; i < seq.length(); i++) {
            char c=seq.charAt(i);
            //test all characters of query
            if (s.isAmbiguous(c)) { 
                //counts number of ambiguities contained in kmers
                for (int j=i-k+1;j<i+1;j++) {
                    if (j>-1 && j<seq.length()) {
                        ambiguityCountPerMer[j]++;
                        //register ambiguity offset
                        if (!ambiguityOffsets.containsKey(j)) {
                            ambiguityOffsets.put(j,new HashMap<>());
                        }
                        ambiguityOffsets.get(j).put(i-j,s.ambiguityEquivalence(c));
                    }
                }
                //register this position as ambiguity
                sequence[i]=-1;
            } else { 
                //state either non ambiguous or stop process
                try {
                    sequence[i]=s.stateToByte(c);
                } catch (NonSupportedStateException ex) {
                    ex.printStackTrace(System.err);
                    System.out.println("Query contains a non supported state: '"+c+"'");
                    System.exit(1);
                }
            }
        }
//        Infos.println("Binary seq: "+Arrays.toString(sequence));
//        Infos.println("#ambiguity: "+Arrays.toString(ambiguityCountPerMer));
//        Infos.println("ambiguityOffsets: ");
//        for (Iterator<Integer> iterator1 = ambiguityOffsets.keySet().iterator(); iterator1.hasNext();) {
//            Integer next = iterator1.next();
//            for (Iterator<Integer> iterator2 = ambiguityOffsets.get(next).keySet().iterator(); iterator2.hasNext();) {
//                Integer next1 = iterator2.next();
//                Infos.println(next+": "+next1+" :"+Arrays.toString(ambiguityOffsets.get(next).get(next1)));
//            }
//        }
        
        
        switch (samplingMode) {
            case SAMPLING_LINEAR:
                merOrder=new int[seq.length()-k+1];
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
                merOrder=new int[seq.length()-k+1];
                shuffledMerOrder();
                this.step=1;
                break;
            case SAMPLING_SEQUENTIAL:
                merOrder=new int[seq.length()-k+1];
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
     * number of mers that will be provided by this knife (sampling method dependant)
     * @return 
     */
    public int getMerCount() {
        return merOrder.length;
    }
    
    /**
     * max number of words that can be built from this sequence (length-k+1)/s
     * @return 
     */
    public int getMaxMerCount() {
        return (this.sequence.length-this.k+1)/this.step;
    }
    
    /**
     * must be called to retireve mers one by one, when ambiguities are
     * generating alternative kmer, the byte[] will be length%k>1 and its
     * length will be a multiple of k
     * @return the next mer(s) as a byte[], null when no more mers to return
     */
    @Override
    public byte[] getNextByteWord() {
        
        if (merIterator>merOrder.length-1) {
            return null;
        }
        int currentPosition=merOrder[merIterator];
        int charactersLeft=sequence.length-currentPosition;
        if (charactersLeft>=minK) {
            if (charactersLeft<k) {
                word=Arrays.copyOfRange(sequence, currentPosition, currentPosition+charactersLeft);
            } else {
                word=Arrays.copyOfRange(sequence, currentPosition, currentPosition+k);
            }
            //this kmer is registered as containing no ambiguities
            if (ambiguityCountPerMer[merIterator]<1) {
                merIterator++;
                return word;
            //this mer presents ambiguities
            } else {
                //skips it if too many ambiguities
                if (ambiguityCountPerMer[merIterator]>maxAmbigPerMer) {
                    merIterator++;
                    return getNextByteWord();
                } else {           
                    //build and concat all alternatives
                    int altProduct=1;
                    for (Iterator<Integer> iterator = ambiguityOffsets.get(currentPosition).keySet().iterator(); iterator.hasNext();) {
                        Integer offSet = iterator.next();
                        altProduct*=ambiguityOffsets.get(currentPosition).get(offSet).length;
                    }
                    byte[] words=new byte[altProduct*k];
                    //build alternative words
                    //first copy non ambigous bases
                    for (int i = 0; i < k; i++) {
                        if (sequence[currentPosition+i]!=-1) {
                            for (int j = 0; j < altProduct; j++) {
                                words[i+j*k]=sequence[currentPosition+i];
                            }                            
                        } else {
                            //repeat alternatives times other ambiguities
                            int jump=0;
                            for (int step = 0; step < altProduct/ambiguityOffsets.get(currentPosition).get(i).length; step++) {
                                for (int j = 0; j < ambiguityOffsets.get(currentPosition).get(i).length; j++) {
                                    words[i+jump*k]=ambiguityOffsets.get(currentPosition).get(i)[j];
                                    jump++;
                                }
                            }
                        }
                    }
                    merIterator++;
                    return words;
                }
                
            }
            
        } else {
            //Infos.println("Skip word on position "+currentPosition+": length < minK !");
            //this allow to skip words that are too short but in the middle
            //of the shuffled mer order, we just skip them and go to the next one.
            merIterator++;
            return getNextByteWord();
        }
    }

    @Override
    public QueryWord getNextWord() {
        return new QueryWord(getNextByteWord(), merIterator++);
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

    
}
