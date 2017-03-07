/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import core.older.Colmer;
import core.older.ColmerSet;
import core.ProbabilisticWord;
import core.SimpleWord;
import etc.Environement;
import java.util.ArrayList;
import java.util.Arrays;
import tree.PhyloNode;

/**
 *
 * @author ben
 */
public class WordGenerator {
    
    float[] wordProbas=null;
    ArrayList<Float> wordProbas2=null;

    public WordGenerator() {
        
    }
    
    public static void main(String[] args) {
        

        byte[] alphabet= {1, 2, 3, 4};
        Environement.printMemoryUsageDescription();
        long startTime = System.currentTimeMillis();
        byte[][] words=WordGenerator.getAllWords(3, 4);
        long endTime = System.currentTimeMillis();
        System.out.println("ALL word generation used " + (endTime - startTime) + " ms");
        Environement.printMemoryUsageDescription();
        startTime = System.currentTimeMillis();
        
        //display words
        for (int i = 0; i < words.length; i++) {
            byte[] word = words[i];
            for (int j = 0; j < word.length; j++) {
                byte b = word[j];
                System.out.print(b);
            }
            System.out.println("");
        }
  
        
    }
    @Deprecated
    public float[] getAssociatedProbas() {
        return wordProbas;
    }
    @Deprecated
    public ArrayList<Float> getAssociatedProbas2() {
        return wordProbas2;
    }
    
    /**
     * generate all possible words of a defined length\n using bytes[],
     * the alphabet allows max 127 different letters, more than enough for DNA/prot
     * @param wordLength  the length of the word (i.e k)
     * @param alphabetSize will be 4 for DNA, 20 for AA
     * @return 
     */
    @Deprecated
    public static byte[][] getAllWords(int wordLength, int alphabetSize) {
        
        byte[] alphabet=new byte[alphabetSize];
        for (byte i = 0; i < alphabetSize; i++) {
            alphabet[i]=i;
        }
        
        int[] index=new int[wordLength];
        Arrays.fill(index, 0);
        int wordsCounter=0;
        byte[][] words=new byte[new Double(Math.pow(alphabet.length, wordLength)).intValue()][alphabet.length];

        while(true) {
            byte[] word=new byte[wordLength];
            for (int i = 0; i < wordLength; ++i) {
                word[i] = alphabet[index[i]];
            }
            words[wordsCounter]=word;
            wordsCounter++;
            for (int i = wordLength-1; ; --i) { 
                if (i < 0) {return words;}
                    index[i]++;
                if (index[i] == alphabet.length) {
                    index[i] = 0;
                } else {
                    break;
                }
            }
        }
    }
    
    /**
     * generate words for this alphabet if the corresponding probabilities are meaningful\n
     * for a particular node, it may be useless to generate a whole batch of words if the probability of a particular\n
     * character, at a particular position is ~0
     * @param ppSet table of posterior probabilities  ppSet[site][stateProba], ex: ppSet[7][3]= proba of "C" at 8th site
     * @param stateTreshold minimum pp value for a state being considered
     * @param generateProbas generate the pp product (pp*), used only for testing how much costs the product
     * @param wordTreshold minimum pp* for returning a word.
     * @return the list of generated words, which survived the tresholds
     */
    @Deprecated
    public byte[][] generateProbableWords(float[][] ppSet, double stateTreshold, boolean generateProbas, double wordTreshold) {

        int wordsCounter=0;
        
        //set alphabet directly from ppSet, for dna = {1,2,3,4}
        byte[] alphabet=new byte[ppSet[0].length];
        for (byte i = 0; i < alphabet.length; i++) {
            alphabet[i]=i;
        }
        //set wordLength from ppSet
        int wordLength=ppSet.length;
        //index used for the incremental generation of the words
        int[] index=new int[wordLength];
        Arrays.fill(index, 0);

        //table to settle once for all the states that should be skipped
        //true if the state at this position is below the given treshold
        //i.e. the state will be considered in word generation
        //boolean[site][state]
        boolean[][] skippedStates=new boolean[wordLength][alphabet.length];
        int testedWordCount=1;
        for (int i = 0; i < wordLength; i++) {
            int testedStates=alphabet.length; //a priori, all states will be tested
            for (int j = 0; j < alphabet.length; j++) {
                if (ppSet[i][j]<stateTreshold) {
                    skippedStates[i][j]=true;
                    testedStates--; //decrease the number of states that will be tested at this site
                }
            }
            testedWordCount=testedWordCount*testedStates;
        }
        //set the array to the number of words that will be tested
        byte[][] words=new byte[testedWordCount][alphabet.length];
        if (generateProbas) {
            wordProbas=new float[words.length];
        }
        
        //launch algo
        // note that we generat the words from right to left,
        //i.e. for k=4 0000,0001,0002,0003,0010,0011,0012,0013,0020,0021 ...
        //if for this example, skippedStates[2][1]=ture, the loop is shortened such as
        //0000,0001,0002,0003,0020,0021 ... 0100,0101;0102,0103,0120,0121,0122,0123,...
        //all xx1x are skipped
        
        //int loop=0;
        while(true) {
            byte[] word=new byte[wordLength];
            float proba=1.0f;
            //build the word and its prob product from the current index
            for (int i = 0; i<wordLength; ++i) { 
                word[i] = alphabet[index[i]];
                if (generateProbas) {
                    proba=proba*ppSet[i][index[i]];
                }
            }
            //register the word and associated proba if > to wordTreshold
            if (proba>=wordTreshold) {
                words[wordsCounter]=word;
                if (generateProbas) {
                    wordProbas[wordsCounter]=proba;
                }
                wordsCounter++;
            }
            for (int pos = wordLength-1; ; --pos) { 
                //we have process all sites from to to first, so we return the list of words
                if (pos < 0) {                    
                    return words;
                }
                //if the next position/state must be skipped we increase the index by 2
                //at this stage, index is 1 incrementation less than the value that may be skipped, so +1 when getting (index[pos]+1)
                if (index[pos]+1<alphabet.length && skippedStates[pos][index[pos]+1]) {
                    index[pos]=index[pos]+2;
                } else {
                    index[pos]++;
                }
                if (index[pos] == alphabet.length) {
                    index[pos] = 0;
                } else {
                    break;
                }
            }
        }
    }
    
    /**
     * generate words for this alphabet if the corresponding probabilities are meaningful\n
     * for a particular node, it may be useless to generate a whole batch of words if the probability of a particular\n
     * character, at a particular position is ~0
     * @param ppSet table of posterior probabilities  ppSet[nodeId][stateProba], ex: ppSet[23][3], pp of state c for node with id 23
     * @param stateTreshold minimum pp value for a state being considered
     * @param generateProbas generate or note the pp product (pp*), used only for testing
     * @param wordTreshold minimum pp* for returning a word.
     * @return the list of generated words, which survived the tresholds
     */
    @Deprecated
    public ArrayList<byte[]> generateProbableWords2(float[][] ppSet, double stateTreshold, boolean generateProbas, double wordTreshold) {
        
        //set alphabet directly from ppSet, for dna = {1,2,3,4}
        byte[] alphabet=new byte[ppSet[0].length];
        for (byte i = 0; i < alphabet.length; i++) {
            alphabet[i]=i;
        }
        //set wordLength from ppSet
        int wordLength=ppSet.length;
        //index used for the incremental generation of the words
        int[] index=new int[wordLength];
        Arrays.fill(index, 0);

        //table to settle once for all the states that should be skipped
        //true if the state at this position is below the given treshold
        //i.e. the state will be considered in word generation
        //boolean[site][state]
        boolean[][] skippedStates=new boolean[wordLength][alphabet.length];
        int testedWordCount=1;
        for (int i = 0; i < wordLength; i++) {
            int testedStates=alphabet.length; //a priori, all states should be tested
            for (int j = 0; j < alphabet.length; j++) {
                if (ppSet[i][j]<stateTreshold) {
                    skippedStates[i][j]=true;
                    testedStates--; //descrese the number of states that will be tested at this site
                }
            }
            testedWordCount=testedWordCount*testedStates;
        }
        //set the array to the number of words that will be tested
        //initial capacity of the arraylist is arbitrary, but tests showed that
        //only ~ 10000 words survived the 2 tresholds if both set at 1e-6
        //table realloc will occur if we get behind this size,
        //I was initally using testedWordCount as the maximum size, but this is
        //a memory killer for nothing...
        int initalCapacity=10000;
        ArrayList<byte[]> words=new ArrayList(initalCapacity);
        if (generateProbas) {
            wordProbas2=new ArrayList<>(initalCapacity);
        }
        
        
//            for (int j = 0; j < skippedStates.length; j++) {
//                for (int i = 0; i <  skippedStates[0].length; i++) {
//                    System.out.print(skippedStates[j][i]+ " "); 
//                }
//                System.out.println("");
//            }
        
        
        
        //launch algo
        // note that we generat the words from right to left,
        //i.e. for k=4 0000,0001,0002,0003,0010,0011,0012,0013,0020,0021 ...
        //for this example, skippedStates[2][1]=ture, the loop is shortened to
        //0000,0001,0002,0003,0020,0021,...,0100,0101;0102,0103,0120,0121,0122,0123,...
        //all xx1x are basically skipped
        
        //int loop=0;
        while(true) {
//            System.out.println("loop "+loop);
            
            
            byte[] word=new byte[wordLength];
            float proba=1.0f;
            //build the word and its prob product from the current index
            for (int i = 0; i<wordLength; ++i) { 
                word[i] = alphabet[index[i]];
                if (generateProbas) {
                    proba=proba*ppSet[i][index[i]];
                }
            }
            //register the word and associated proba if > to wordTreshold
            if (proba>=wordTreshold) {
                words.add(word);
                if (generateProbas) {
                    wordProbas2.add(proba);
                }
            }
//            System.out.print("WORD ");
//            for (int j = 0; j < word.length; j++) {
//                System.out.print(word[j]+" ");
//            }
//            System.out.println("");

            
//                System.out.print("INDEXA ");
//                for (int j = 0; j < index.length; j++) {
//                    System.out.print(index[j]+" ");
//                }
//                System.out.println("");
            
            for (int pos = wordLength-1; ; --pos) { 
//                System.out.println("pos- "+pos);
                //we have process all sites from to to first, so we return the list of words
                if (pos < 0) {
                    //attempt to explicitely flag some useless stuff for the garbage collector
                    index=null;
                    alphabet=null;
                    skippedStates=null;
                    return words;
                }
                //if the next position/state must be skipped we increase the index by 2
                //at this stage, index is 1 incrementation less than the value that may be skipped, so +1 when getting (index[pos]+1)
                if (index[pos]+1<alphabet.length && skippedStates[pos][index[pos]+1]) {
//                    System.out.println("SKIP "+pos+" "+index[pos]);
                    index[pos]=index[pos]+2;
                } else {
                    index[pos]++;
                }
                
                
//                System.out.print("INDEXB ");
//                for (int j = 0; j < index.length; j++) {
//                    System.out.print(index[j]+" ");
//                }
//                System.out.println("");
                if (index[pos] == alphabet.length) {
                    index[pos] = 0;
                } else {
                    break;
                }
            }
            //loop++;
        }
        
        
    }

    /**
     * generate words for this alphabet if the corresponding probabilities are meaningful\n
     * for a particular node, it may be useless to generate a whole batch of words if the probability of a particular\n
     * character, at a particular position is ~0
     * @param ppSet table of posterior probabilities  ppSet[nodeId][stateProba], ex: ppSet[23][3], pp of state c for node with id 23
     * @param stateProbaTreshold minimum pp value for a state being considered
     * @param wordProbaTreshold minimum pp* for returning a word.
     * @return the list of generated words, which survived the thresholds
     */
    public ArrayList<ProbabilisticWord> generateProbableWords3(int refPosition, float[][] ppSet, double stateProbaTreshold, double wordProbaTreshold) {
        
        //set alphabet directly from ppSet, for dna = {1,2,3,4}
        byte[] alphabet=new byte[ppSet[0].length];
        for (byte i = 0; i < alphabet.length; i++) {
            alphabet[i]=i;
        }
        //set wordLength from ppSet
        int wordLength=ppSet.length;
        //index used for the incremental generation of the words
        int[] index=new int[wordLength];
        Arrays.fill(index, 0);

        //table to settle once for all the states that should be skipped
        //true if the state at this position is below the given treshold
        //i.e. the state will be considered in word generation
        //boolean[site][state]
        boolean[][] skippedStates=new boolean[wordLength][alphabet.length];
        int testedWordCount=1;
        for (int i = 0; i < wordLength; i++) {
            int testedStates=alphabet.length; //a priori, all states should be tested
            for (int j = 0; j < alphabet.length; j++) {
                if (ppSet[i][j]<stateProbaTreshold) {
                    skippedStates[i][j]=true;
                    testedStates--; //descrese the number of states that will be tested at this site
                }
            }
            testedWordCount=testedWordCount*testedStates;
        }
        //set the array to the number of words that will be tested
        //initial capacity of the arraylist is arbitrary, but tests showed that
        //only ~ 10000 words survived the 2 tresholds if both set at 1e-6
        //table realloc will occur if we get behind this size,
        //I was initally using testedWordCount as the maximum size, but this is
        //a memory killer for nothing...
        int initalCapacity=10000;
        ArrayList<ProbabilisticWord> words=new ArrayList(initalCapacity);

        //launch algo
        // note that we generat the words from right to left,
        //i.e. for k=4 0000,0001,0002,0003,0010,0011,0012,0013,0020,0021 ...
        //for this example, skippedStates[2][1]=ture, the loop is shortened to
        //0000,0001,0002,0003,0020,0021,...,0100,0101;0102,0103,0120,0121,0122,0123,...
        //all xx1x are basically skipped
        
        while(true) {
            byte[] word=new byte[wordLength];
            float proba=1.0f;
            //build the word and its prob product from the current index
            for (int i = 0; i<wordLength; ++i) { 
                word[i] = alphabet[index[i]];
                proba=proba*ppSet[i][index[i]];
                if (proba<wordProbaTreshold) 
                    break;
            }
            //register the word and associated proba if > to wordTreshold
            if (proba>=wordProbaTreshold) {
                words.add(new ProbabilisticWord(word, proba,refPosition));
            }
            for (int pos = wordLength-1; ; --pos) { 
                //we have process all sites from to to first, so we return the list of words
                if (pos < 0) {
                    //attempt to explicitely flag some useless stuff for the garbage collector
                    index=null;
                    alphabet=null;
                    skippedStates=null;
                    return words;
                }
                //if the next position/state must be skipped we increase the index by 2
                //at this stage, index is 1 incrementation less than the value that may be skipped, so +1 when getting (index[pos]+1)
                if (index[pos]+1<alphabet.length && skippedStates[pos][index[pos]+1]) {
                    index[pos]=index[pos]+2;
                } else {
                    index[pos]++;
                }
                if (index[pos] == alphabet.length) {
                    index[pos] = 0;
                } else {
                    break;
                }
            }
            //loop++;
        }
        
        
    }


    /**
     * generate words for this alphabet if the corresponding probabilities are meaningful\n
     * for a particular node, it may be useless to generate a whole batch of words if the probability of a particular\n
     * character, at a particular position is ~0\n
     * in this version, the registry of the Colmer is directly filled to economize memory
     * @param cs the @ColmerSet in which to register the generated words
     * @param c the current colmer for which words are generated
     * @param n the current node for which probas are generated
     * @param ppSet table of posterior probabilities  ppSet[nodeId][stateProba], ex: ppSet[23][3], pp of state c for node with id 23
     * @param stateProbaTreshold minimum pp value for a state being considered
     * @param wordProbaTreshold minimum pp* for returning a word.
     * @return the list of generated words, which survived the tresholds
     */
    public int generateProbableWords4(ColmerSet cs, Colmer c, PhyloNode n, float[][] ppSet, double stateProbaTreshold, double wordProbaTreshold) {
        
        boolean test=false;
//        if ( (n.getId()==23 || n.getId()==36) && (c.getStartSite()>1947) && (c.getStartSite()<1950)  ) {
//            test=true;
//            System.out.println("CURRENT NODE: id="+n.getId()+" label="+n.getLabel());
//        }
        
        
        //set alphabet directly from ppSet, for dna = {1,2,3,4}
        byte[] alphabet=new byte[ppSet[0].length];
        for (byte i = 0; i < alphabet.length; i++) {
            alphabet[i]=i;
        }
        //set wordLength from ppSet
        int wordLength=ppSet.length;
        //index used for the incremental generation of the words
        int[] index=new int[wordLength];
        Arrays.fill(index, 0);
        
        if (test)
            for (int i = 0; i < ppSet.length; i++) {
                for (int j = 0; j < ppSet[i].length; j++) {
                    System.out.print(ppSet[i][j]+";");
                }
                System.out.println("");
            }

        //table to settle once for all the states that should be skipped
        //true if the state at this position is below the given treshold
        //i.e. the state will be considered in word generation
        //boolean[site][state]
        boolean[][] skippedStates=new boolean[wordLength][alphabet.length];
        int testedWordCount=1;
        for (int i = 0; i < wordLength; i++) {
            int testedStates=alphabet.length; //a priori, all states should be tested
            for (int j = 0; j < alphabet.length; j++) {
                if (ppSet[i][j]<stateProbaTreshold) {
                    skippedStates[i][j]=true;
                    if (test)
                        System.out.println("Will skip state (ppSet[nodeId][stateProba]): ["+i+"]["+j+"]="+ppSet[i][j]);
                    testedStates--; //descrese the number of states that will be tested at this site
                }
            }
            testedWordCount=testedWordCount*testedStates;
        }

        //launch algo
        // note that we generat the words from right to left,
        //i.e. for k=4 0000,0001,0002,0003,0010,0011,0012,0013,0020,0021 ...
        //for this example, skippedStates[2][1]=ture, the loop is shortened to
        //0000,0001,0002,0003,0020,0021,...,0100,0101;0102,0103,0120,0121,0122,0123,...
        //all xx1x are basically skipped
        int generatedWordsCount=0;
        while(true) {
            byte[] word=new byte[wordLength];
            double proba=1.0;
            //build the word and its prob product from the current index
            for (int i = 0; i<wordLength; ++i) { 
                word[i] = alphabet[index[i]];
                proba=proba*ppSet[i][index[i]];
            }
            
            if (test)
                System.out.println("Node_id"+n.getId()+" colmerStart: "+c.getStartSite()+" Generated: "+Arrays.toString(word)+" PP*="+proba);
            
            //register the word and associated proba if > to wordTreshold
            if (proba>=wordProbaTreshold) {
                generatedWordsCount++;
                if (test)
                    System.out.println("REGISTERED!");                
                cs.registerWord(new SimpleWord(word));
                cs.associateLocality(new SimpleWord(word),n,c,proba);
            }
            for (int pos = wordLength-1; ; --pos) { 
                //we have process all sites from to to first, so we return the list of words
                if (pos < 0) {
                    //attempt to explicitely flag some useless stuff for the garbage collector
                    index=null;
                    alphabet=null;
                    skippedStates=null;
                    return generatedWordsCount;
                }
                //if the next position/state must be skipped we increase the index by 2
                //at this stage, index is 1 incrementation less than the value that may be skipped, so +1 when getting (index[pos]+1)
                if (index[pos]+1<alphabet.length && skippedStates[pos][index[pos]+1]) {
                    index[pos]=index[pos]+2;
                } else {
                    index[pos]++;
                }
                if (index[pos] == alphabet.length) {
                    index[pos] = 0;
                } else {
                    break;
                }
            }
            //loop++;
        }
    }

   
    
    
}
