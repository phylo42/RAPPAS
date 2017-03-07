/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import core.older.Colmer;
import core.older.ColmerSet;
import core.DNAStates;
import core.ProbabilisticWord;
import core.older.PProbas;
import core.States;
import core.Word;
import core.algos.SequenceKnife;
import core.algos.WordGenerator;
import etc.Environement;
import etc.Infos;
import inputs.Fasta;
import inputs.FASTAPointer;
import inputs.InputManager;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.ExtendedTree;
import tree.Tree;

/**
 *
 * @author ben
 */
public class Session {
    
    int k=5;
    int minK=5;
    double stateThreshold=1e-5;
    double wordThreshold=1e-5;
            
    
    
    States states=null;
    Alignment align=null;
    Tree tree=null;
    PProbas parsedProbas=null;
    ExtendedTree relaxedTree=null; //optinnal, checked through a boolean flag
    
    ColmerSet cs=null;
    
    
    public Session(int k, int mink, double stateThreshold, double wordThreshold) {
        this.k=k;
        this.minK=mink;
        this.stateThreshold=stateThreshold;
        this.wordThreshold=wordThreshold;
    }
    
    public void associateStates(States s) {
        this.states=s;
    }
    
    public void associateInputs(InputManager im) {
        this.tree=im.getTree();
        if (im.getTree() instanceof ExtendedTree)
            this.relaxedTree=(ExtendedTree)im.getTree();
        else
            this.relaxedTree=null;
        this.align=im.getAlignment();
        this.parsedProbas=im.getPProbas();
    }
    
    
    public void generateColmerSetData(boolean enforceGC,int samplingMode) {
        Infos.println("wordLength: " + k);
        this.cs=new ColmerSet(align, tree, states,k,minK,samplingMode);
        Infos.println("Tree structure:");
        ArrayList<Integer> orderedNodes = tree.getNodeIdsByDFS();
        for (int i = 0; i < orderedNodes.size(); i++) {
            Integer get = orderedNodes.get(i);
            PhyloNode n=tree.getById(get);
            Infos.println("Node: "+n.toString());
        }
        Infos.println("# nodes for which Words will be generated: " + orderedNodes.size());
        Infos.println("Colmers sampled as "+Arrays.toString(Arrays.copyOfRange(cs.getMerOrder(),0,15))+" ...");
        Environement.printMemoryUsageDescription();
        Colmer c=null;
        long meanTime=0;
        while ((c=cs.getNextColmer())!=null) {
            Infos.println( "current colmer: "+c.getStartSite()+";"+c.getEndSite());
            long startTime = System.currentTimeMillis();
            int totalWords=0;
            for (int i = 0; i < orderedNodes.size(); i++) {
                PhyloNode n=tree.getById(orderedNodes.get(i));
                if (n.isLeaf()) {continue;}
                //Infos.println("  current node: "+n);
                if (!cs.isFull()) {
                    //all these words are registered automatically in the hash of the ColmerSet
                    WordGenerator wg=new WordGenerator();
                    totalWords+=wg.generateProbableWords4(  cs,
                            c,
                            n,
                            parsedProbas.getPPSet(n.getId(), c.getStartSite(), c.getEndSite()),
                            stateThreshold,
                            wordThreshold
                    );
                    wg=null;     
                }
            }
            //after having build the colmer, register if in the colmerset
            cs.addColmer(c);
            if (enforceGC) {
                System.gc();
            }
            long endTime = System.currentTimeMillis();
            long time=endTime - startTime;
            meanTime+=time;
            Infos.println("final avg node size (words) " +(totalWords/tree.getNodeCount())+"  (max="+Double.valueOf(Math.pow(4, k)).intValue()+")");
            Infos.println("Colmer generation " + time + " ms");
            Environement.printMemoryUsageDescription();
            
//                if (c.getStartSite()>300) {
//                    break;
//                }
        }
        Infos.println("Average Colmer generation time: " + (meanTime/cs.getColmerCount()) + " ms");
        Infos.println("Total words registered in the ColmerSet: " +cs.getSetSize());
    }
    
    public void generateColmerSetDataForInterval(boolean enforceGC,int samplingMode,int intervalStart,int intervalEnd) {
        Infos.println("wordLength: " + k);
        this.cs=new ColmerSet(align, tree, states,k,minK,samplingMode);
        Infos.println("Tree structure:");
        ArrayList<Integer> orderedNodes = tree.getNodeIdsByDFS();
        for (int i = 0; i < orderedNodes.size(); i++) {
            Integer get = orderedNodes.get(i);
            PhyloNode n=tree.getById(get);
            Infos.println("Node: "+n.toString());
        }
        Infos.println("# nodes for which Words will be generated: " + orderedNodes.size());
        Infos.println("Colmers sampled as "+Arrays.toString(Arrays.copyOfRange(cs.getMerOrder(),0,15))+" ...");
        Environement.printMemoryUsageDescription();
        Colmer c=null;
        long meanTime=0;
        while ((c=cs.getNextColmer())!=null) {
            Infos.println( "current colmer: "+c.getStartSite()+";"+c.getEndSite());
            
            if ( (c.getStartSite()<intervalStart) || (c.getStartSite()>intervalEnd) ) {
                Infos.println("Not in interval set by the user, skipped.");
                //add an empty colmer in which no words are registered
                cs.addColmer(c);
                continue;
            }
            
            long startTime = System.currentTimeMillis();
            int totalWords=0;
            for (int i = 0; i < orderedNodes.size(); i++) {
                PhyloNode n=tree.getById(orderedNodes.get(i));
                if (n.isLeaf()) {continue;}
                //Infos.println("  current node: "+n);
                if (!cs.isFull()) {
                    //all these words are registered automatically in the hash of the ColmerSet
                    WordGenerator wg=new WordGenerator();
                    totalWords+=wg.generateProbableWords4(  cs,
                            c,
                            n,
                            parsedProbas.getPPSet(n.getId(), c.getStartSite(), c.getEndSite()),
                            stateThreshold,
                            wordThreshold
                    );
                    wg=null;     
                }
            }
            //after having build the colmer, register if in the colmerset
            cs.addColmer(c);
            if (enforceGC) {
                System.gc();
            }
            long endTime = System.currentTimeMillis();
            long time=endTime - startTime;
            meanTime+=time;
            Infos.println("final avg node size (words) " +(totalWords/tree.getNodeCount())+"  (max="+Double.valueOf(Math.pow(4, k)).intValue()+")");
            Infos.println("Colmer generation " + time + " ms");
            Environement.printMemoryUsageDescription();
            
//                if (c.getStartSite()>300) {
//                    break;
//                }
        }
        Infos.println("Average Colmer generation time: " + (meanTime/cs.getColmerCount()) + " ms");
        Infos.println("Total words registered in the ColmerSet: " +cs.getSetSize());
    }

    
    public void testMethod() {
        
        File coronavirusReads=new File("/media/ben/STOCK/DATA/NCBI_VIROMES/blast_extractions/coronavirus/extract_ERR1360080_1e-25_95_80");

        FASTAPointer fpp=new FASTAPointer(coronavirusReads,true);
        
        int limit=10;
        int readCounter=0;
        Fasta fasta=null;
        while ((fasta=fpp.nextSequenceAsFastaObject())!=null) {
            
            if (readCounter>limit) {break;}
            
            System.out.println("Read: "+fasta.getHeader()+"\n"+fasta.getSequence());
            SequenceKnife knife=new SequenceKnife(fasta.getSequence(),k , minK, new DNAStates(), SequenceKnife.SAMPLING_LINEAR);
            System.out.println("Mer order"+Arrays.toString(knife.getMerOrder()));
            Word w=null;
            while (( w=knife.getNextWord())!=null) {
                System.out.println("  Query word: "+w);
                //List<QueryWord> mutatedWords = w.getMutatedWords(new DNAStates());
//                ArrayList<Integer> queryMatch = cs.getQueryMatch(w);
//                if (queryMatch!=null) {
//                    for (int i = 0; i < queryMatch.size(); i++) {
//                        Colmer co = cs.getColmerById(queryMatch.get(i));
//                        System.out.println("     Colmer:"+co.getId()+" ("+co.getStartSite()+"-"+co.getEndSite()+")");
////                        LinkedHashMap<Integer, Double> ancetralScores = co.getAncetralScores(w);
////                        for (Iterator<Integer> iterator = ancetralScores.keySet().iterator(); iterator.hasNext();) {
////                            PhyloNode next = tree.getById(iterator.next());
////                            System.out.println("       "+next.getLabel()+":"+ancetralScores.get(next));
////                        }
//                    }
//                }
                
            }
            readCounter++;
        }
        fpp.closePointer();
        System.exit(0);
        
    }
    
    public boolean store(File f) {
        try {
            long startTime = System.currentTimeMillis();
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
            oos.writeInt(k);
            oos.writeInt(minK);
            oos.writeDouble(stateThreshold);
            oos.writeDouble(wordThreshold);
            //flag for relaxedTree
            if (this.relaxedTree!=null)
                oos.writeBoolean(true);
            else
                oos.writeBoolean(false);
            Infos.println("Storing of States");
            oos.writeObject(states);
            Infos.println("Storing of Alignment");
            oos.writeObject(align);
            Infos.println("Storing of PhyloTree");
            oos.writeObject(tree);
            Infos.println("Storing of RelaxedTree");
            if (this.relaxedTree!=null)
                oos.writeObject(relaxedTree);
            Infos.println("Storing of PPStats");
            oos.writeObject(parsedProbas);
            Infos.println("Storing of ColmerSet");
            oos.writeObject(cs);
            oos.close();
            fos.close();
            long endTime = System.currentTimeMillis();
            Infos.println("Complete session storage " + (endTime - startTime) + " ms");
            Infos.println("Session stored in : "+f.getAbsolutePath());
        } catch (IOException ex) {
            ex.printStackTrace();
            return false;
        }
        return true;
    }
    
    public static Session load(File f) {
        try {
            long startTime = System.currentTimeMillis();
            FileInputStream fis = new FileInputStream(f);
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,4096));
            int k=ois.readInt();
            int minK=ois.readInt();
            double stateThreshold=ois.readDouble();
            double wordThreshold=ois.readDouble();
            boolean relaxedTree=ois.readBoolean();
            Session s=new Session(k, minK, stateThreshold, wordThreshold);
            Infos.println("Loading States");
            s.states = (States)ois.readObject();
            Infos.println("Loading Alignment");
            s.align = (Alignment)ois.readObject();
            Infos.println("Loading PhyloTree");
            s.tree = (PhyloTree)ois.readObject();
            Infos.println("Loading RelaxedTree");
            if (relaxedTree)
                s.tree = (ExtendedTree)ois.readObject();
            Infos.println("Loading of PPStats");
            s.parsedProbas = (PProbas)ois.readObject();
            Infos.println("Loading of ColmerSet");
            s.cs = (ColmerSet)ois.readObject();
            ois.close();
            fis.close();
            long endTime = System.currentTimeMillis();
            Infos.println("Complete session loading " + (endTime - startTime) + " ms");
            return s;
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
            return null;
        }
    }



    
    
    @Deprecated
    private List<ProbabilisticWord> generateWords2(PhyloNode n,PProbas stats,int siteStart, int siteEnd) {

            //long startTime = System.currentTimeMillis();
            //Infos.println("Word generation for node: "+n.getLabel()+" ("+n.getId()+")");
            
            int wordLength=siteEnd-siteStart+1;
            
            WordGenerator wg=new WordGenerator();
            
            final ArrayList<byte[]> probableWords2 = wg.generateProbableWords2(stats.getPPSet(n.getId(), siteStart, siteStart+wordLength-1),1e-6,true,1e-6);
            //System.out.println("Length probableWords2: "+probableWords2.size());
            final ArrayList<Float> generatedProbas = wg.getAssociatedProbas2();
            
            final ArrayList<ProbabilisticWord> wordSet=new ArrayList<>(generatedProbas.size());
            AtomicInteger ii=new AtomicInteger();
            generatedProbas.stream().forEachOrdered((d) ->  {
                                                            wordSet.add(new ProbabilisticWord(probableWords2.get(ii.get()), d));
                                                            ii.incrementAndGet();
                                                            }
                                                        );
            Collections.sort(wordSet);
            
            wg=null;
            
            return wordSet.stream().filter((w) -> w.getPpStarValue()> 1e-6 ).collect(Collectors.toCollection(ArrayList::new));
            //long endTime = System.currentTimeMillis();
            //Infos.println("Generation took: " + (endTime - startTime) + " ms ("+selectedWords.size()+" words)");  
                        
    }
    

    
    
    private void printFirst(List<ProbabilisticWord> l,int n) {
        if (n>l.size()){n=l.size();}
        for (int i = 0; i < n; i++) {
            System.out.println(l.get(i).toStringNice(states));
        }
    }
 
    
    @Deprecated
    private List<ProbabilisticWord> generateWords(PhyloNode n,PProbas stats,int siteStart, int siteEnd) {

            //long startTime = System.currentTimeMillis();
            //Infos.println("Word generation for node: "+n.getLabel()+" ("+n.getId()+")");
            
            int wordLength=siteEnd-siteStart+1;
            
            WordGenerator wg=new WordGenerator();
                        
            
            final byte[][] probableWords2 = wg.generateProbableWords(stats.getPPSet(n.getId(), siteStart, siteStart+wordLength-1),1e-250,true,1e-250);
            System.out.println("Length probableWords2: "+probableWords2.length);
            final float[] generatedProbas = wg.getAssociatedProbas();
            
            final ArrayList<ProbabilisticWord> wordSet=new ArrayList<>(generatedProbas.length);
            AtomicInteger ii=new AtomicInteger();
            
            
//            Arrays.stream(generatedProbas).forEachOrdered((d) ->  {
//                                                                    wordSet.add(new ProbabilisticWord(probableWords2[ii.get()], d));
//                                                                    ii.incrementAndGet();
//                                                                    }
//                                                        );
            Collections.sort(wordSet);
            
            return wordSet.stream().filter((w) -> w.getPpStarValue()> 1e-6 ).collect(Collectors.toCollection(ArrayList::new));
            //long endTime = System.currentTimeMillis();
            //Infos.println("Generation took: " + (endTime - startTime) + " ms ("+selectedWords.size()+" words)");  
                        
    }
    
    @Deprecated
    private void generateWordsTests() {
// byte[][] probableWords=null;
//            int wordLength=12;
//            
////            startTime = System.currentTimeMillis();
////            probableWords = WordGenerator.getAllWords(wordLength, 4);
////            endTime = System.currentTimeMillis();
////            Infos.println("Generating all words took " + (endTime - startTime) + " ms");
//            
//            
//            startTime = System.currentTimeMillis();
//            double[][] ppSet = parsedProbas.getPPSet(tree.getByName("N2").getId(), 0, wordLength);
//            WordGenerator wg=new WordGenerator();
//            probableWords = wg.generateProbableWords(ppSet,0.00001,false,1e-50);
//            endTime = System.currentTimeMillis();
//            Infos.println("Generating probable words " + (endTime - startTime) + " ms");
////            for (int i = 0; i < 10; i++) {
////                System.out.println(Arrays.toString(probableWords[i]));
////            }
//            
//            startTime = System.currentTimeMillis();
//            ppSet = parsedProbas.getPPSet(tree.getByName("N2").getId(), 0, wordLength);
//            wg=new WordGenerator();
//            probableWords = wg.generateProbableWords(ppSet,0.00001,true,1e-50);
//            double[] generatedProbas = wg.getAssociatedProbas();
//            endTime = System.currentTimeMillis();
//            Infos.println("Generating probable words + probas (simultaneously) " + (endTime - startTime) + " ms");
////            for (int i = 0; i < 10; i++) {
////                System.out.println(Arrays.toString(probableWords[i])+" "+generatedProbas[i]);
////            }
//            
//            
//            startTime = System.currentTimeMillis();
//            ppSet = parsedProbas.getPPSet(tree.getByName("N2").getId(), 0, wordLength);
//            probableWords = wg.generateProbableWords(ppSet,0.00001,false,1e-50);
//            generatedProbas = ProbaGenerator.generateProbas(probableWords, ppSet);
//            endTime = System.currentTimeMillis();
//            Infos.println("Generating probable words + probas (a posteriori) " + (endTime - startTime) + " ms");
////            for (int i = 0; i < 10; i++) {
////                System.out.println(Arrays.toString(probableWords[i])+" "+generatedProbas[i]);
////            }
//            
//            startTime = System.currentTimeMillis();
//            ppSet = parsedProbas.getPPSet(tree.getByName("N2").getId(), 0, wordLength);
//            final byte[][] probableWords2 = wg.generateProbableWords(ppSet,0.00001,true,1e-5);
//            generatedProbas = wg.getAssociatedProbas();
//            endTime = System.currentTimeMillis();
//            Infos.println("Generating probable words + probas (simultaneously) + global treshold " + (endTime - startTime) + " ms");
////            for (int i = 0; i < 10; i++) {
////                System.out.println(Arrays.toString(probableWords[i])+" "+generatedProbas[i]);
////            }
//            
//
//            Infos.println("Word list size:"+probableWords2.length);  
    }
        
    @Deprecated
    private void wrappingAndSortingTests() {
//            startTime = System.currentTimeMillis();
//            List dList = new ArrayList<Double>();
//            for(double dValue : generatedProbas) {
//                dList.add(dValue);
//            }
//            endTime = System.currentTimeMillis();
//            Infos.println("List<Double> loop wrapper " + (endTime - startTime) + " ms");
//            print5First(dList);
//            
//            startTime = System.currentTimeMillis();
//            ArrayList<Double> list = Arrays.stream(generatedProbas).boxed().collect(Collectors.toCollection(ArrayList::new));  
//            endTime = System.currentTimeMillis();
//            Infos.println("ArrayList<Double> stream wrapper" + (endTime - startTime) + " ms");    
//            print5First(list);
//            
//            startTime = System.currentTimeMillis();
//            ArrayList<Word> wordSet=new ArrayList<>(generatedProbas.length);
//            for (int i = 0; i < probableWords2.length; i++) {
//                 wordSet.add(new Word(probableWords2[i],generatedProbas[i]));
//            }
//            endTime = System.currentTimeMillis();
//            Infos.println("ArrayList<Word> loop wrapper took " + (endTime - startTime) + " ms");   
//            print5First(wordSet);
//            
//            startTime = System.currentTimeMillis();
//            final ArrayList<Word> wordSet1=new ArrayList<>(generatedProbas.length);
//            AtomicInteger ii=new AtomicInteger();
//            Arrays.stream(generatedProbas).forEachOrdered((d) ->    {
//                                                                    Word w=new Word(probableWords2[ii.get()], d);
//                                                                    wordSet1.add(w);
//                                                                    ii.incrementAndGet();
//                                                                    }
//                                                         
//                                                        );
//            
//            endTime = System.currentTimeMillis();
//            Infos.println("ArrayList<Word> stream wrapper took " + (endTime - startTime) + " ms");
//            print5First(wordSet1);
//            
//            
//            System.out.println("\n-----------------------------------");
//            System.out.println("-----------------------------------\n");
//            
//            
//
//            
////            startTime = System.currentTimeMillis();  //SLOW
////            //QuickSort.quicksort2(list);  
////            endTime = System.currentTimeMillis();
////            Infos.println("Quicksort2 " + (endTime - startTime) + " ms");    
////            
////            startTime = System.currentTimeMillis();  //SLOW
////            //QuickSort.quicksort4(list); 
////            endTime = System.currentTimeMillis();
////            Infos.println("Quicksort3 " + (endTime - startTime) + " ms");    
////            startTime = System.currentTimeMillis();  //SLOW
////            //QuickSort.quicksort4(wordSet1);
////            endTime = System.currentTimeMillis();
////            Infos.println("QuickSort4(List) " + (endTime - startTime) + " ms");
//            
//            startTime = System.currentTimeMillis();
//            Word[] array=new Word[wordSet1.size()];
//            SmoothSort.sort(wordSet1.toArray(array), 0, wordSet1.size()-1);
//            ArrayList<Word> selectedWords=Arrays.stream(array).filter((w) -> w.getValue()> 1e-5 ).collect(Collectors.toCollection(ArrayList::new));
//            endTime = System.currentTimeMillis();
//            Infos.println("Smoothsort(Word[]) " + (endTime - startTime) + " ms");
//            print5First(selectedWords);
//
//            wordSet=new ArrayList<>(generatedProbas.length);
//            for (int i = 0; i < probableWords2.length; i++) {
//                 wordSet.add(new Word(probableWords2[i],generatedProbas[i]));
//            }
//            startTime = System.currentTimeMillis();
//            HeapSort.heapSort(wordSet);
//            selectedWords=wordSet.stream().filter((w) -> w.getValue()> 1e-5 ).collect(Collectors.toCollection(ArrayList::new));
//            endTime = System.currentTimeMillis();
//            Infos.println("Heapsort(List) " + (endTime - startTime) + " ms");
//            print5First(selectedWords);
//
//            //THIS ONE DON'T WORK
////            wordSet=new ArrayList<>(generatedProbas.length);
////            for (int i = 0; i < probableWords2.length; i++) {
////                 wordSet.add(new Word(probableWords2[i],generatedProbas[i]));
////            }
////            startTime = System.currentTimeMillis();
////            QuickSort.Quicksort1(wordSet, (a,b) -> Double.valueOf(a.getValue()-b.getValue()).intValue() );
////            selectedWords=wordSet.stream().filter((w) -> w.getValue()> 1e-5 ).collect(Collectors.toCollection(ArrayList::new));
////            endTime = System.currentTimeMillis();
////            Infos.println("QuickSort1(List) with lambda expression " + (endTime - startTime) + " ms");  
////            print5First(selectedWords);
//
//            wordSet=new ArrayList<>(generatedProbas.length);
//            for (int i = 0; i < probableWords2.length; i++) {
//                 wordSet.add(new Word(probableWords2[i],generatedProbas[i]));
//            }
//            startTime = System.currentTimeMillis();
//            Collections.sort(wordSet);
//            selectedWords=wordSet.stream().filter((w) -> w.getValue()> 1e-5 ).collect(Collectors.toCollection(ArrayList::new));
//            endTime = System.currentTimeMillis();
//            Infos.println("Java Collections sort " + (endTime - startTime) + " ms");  
//            print5First(selectedWords);
//
//            
//            wordSet=new ArrayList<>(generatedProbas.length);
//            for (int i = 0; i < probableWords2.length; i++) {
//                 wordSet.add(new Word(probableWords2[i],generatedProbas[i]));
//            }
//            startTime = System.currentTimeMillis();
//            selectedWords=wordSet.stream().sorted().filter((w) -> w.getValue()> 1e-5 ).collect(Collectors.toCollection(ArrayList::new));
//            endTime = System.currentTimeMillis();
//            Infos.println("Java Stream sort " + (endTime - startTime) + " ms");  
//            print5First(selectedWords);
//
//            wordSet=new ArrayList<>(generatedProbas.length);
//            for (int i = 0; i < probableWords2.length; i++) {
//                 wordSet.add(new Word(probableWords2[i],generatedProbas[i]));
//            }
//            startTime = System.currentTimeMillis();
//            selectedWords=wordSet.parallelStream().sorted().filter((w) -> w.getValue()> 1e-5 ).collect(Collectors.toCollection(ArrayList::new));
//            endTime = System.currentTimeMillis();
//            Infos.println("Java ParallelStream sort " + (endTime - startTime) + " ms");  
//            print5First(selectedWords);
//            
//            System.exit(0);
//
//            
//            //SmoothSort.sort(args, wordLength, wordLength);
//            
//            //ArrayList<Double> list=new ArrayList<Double>(Arrays.asList(generatedProbas));
//            //list.
//            
//            //QuickSort.quickSort3(Arrays.);
//            
////            startTime = System.currentTimeMillis();
////            ConcurrentSkipListMap<Word,Boolean> wordSet=null;
////            wordSet=new ConcurrentSkipListMap(new WordComparator());
////            for (int i = 0; i < generatedProbas.length; i++) {
////                wordSet.put(new Word(probableWords[i], generatedProbas[i]), Boolean.FALSE);
////            }
////            endTime = System.currentTimeMillis();
////            Infos.println("Generating skiplist took " + (endTime - startTime) + " ms");
//
//
//            
////            for (Iterator<Word> iterator = wordSet.descendingKeySet().iterator(); iterator.hasNext();) {
////                Word next = iterator.next();
////                for (int i = 0; i < next.getKey().length; i++) {
////                    System.out.print(next.getKey()[i]+"; ");
////                }
////                System.out.println(next.getValue());
////                
////            }
//            
//            
//            //display probas/words
////            for (int i = 0; i < ppSet.length; i++) {
////                double[] states = ppSet[i];
////                for (int j = 0; j < states.length; j++) {
////                    double b = states[j];
////                    System.out.print(b+"  ");
////                }
////                System.out.println("");
////            }
////            NumberFormat nf = new DecimalFormat("0.#####E0");
////            nf.setMinimumFractionDigits(5);
////            for (int i = 0; i < probableWords.length; i++) {
////                byte[] word = probableWords[i];
////                for (int j = 0; j < word.length; j++) {
////                    byte b = word[j];
////                    System.out.print(b+"; ");
////                }
////                System.out.println(associatedProbas[i]+" "+generateProbas[i]);
////            }
    }
    
}
