/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.DNAStates;
import core.PProbasSorted;
import core.SiteProba;
import core.States;
import core.older.PProbas;
import etc.Infos;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import main.SessionNext;
import tree.NewickReader;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class PHYMLWrapper implements DataWrapper {
    
    Alignment align=null;
    PhyloTree tree=null;
    States states=null;
    
    /**
     * 
     * @param ARInputAlign input alignment used at ancestral reconstruction
     * @param states
     */

    public PHYMLWrapper(Alignment ARInputAlign, States states) {
        this.align=ARInputAlign;
        this.states=states;
    }
    
    
    
    

    /**
     * parse the tree from the basic.phylip_phyml_tree.txt file
     * @param input
     * @return
     * @throws IOException 
     */
    public PhyloTree parseTree(InputStream input) throws IOException {
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        String originalTreeString=null;
        while ((line=br.readLine())!=null) {
            if (line.trim().length()<1) {continue;}
            originalTreeString=line;
        }
        this.tree=NewickReader.parseNewickTree2(originalTreeString, false);
        return this.tree;
    }
    

    /**
     * parse the posterior probas themselves
     * @param input
     * @param sitePPThreshold
     * @param asLog10
     * @param debugNodeLimit limit PP parsing to few nodes
     * @return
     * @throws IOException 
     */
    @Override
    public PProbasSorted parseSortedProbas(InputStream input, float sitePPThreshold, boolean asLog10, int debugNodeLimit){
        
        BufferedReader br=null;
        int lineNumber=0;
        PProbasSorted matrix=null;
        try {
            matrix=new PProbasSorted(tree.getNodeCount(), align.getLength(), states.getNonAmbiguousStatesCount());
            br = new BufferedReader(new InputStreamReader(input,"UTF-8"));
            String line=null;
            boolean start=false;
            int currentPP=0;
            String currentNode=null;
            int nodeId=-1;
            int previousId=-1;
            int nodeCount=0;
            char[] stateOrder=null;
            ArrayList<SiteProba> probasPerSite=new ArrayList<>(states.getNonAmbiguousStatesCount());
            for (int i = 0; i < states.getNonAmbiguousStatesCount(); i++) {
                probasPerSite.add(new SiteProba());
            }   Infos.println("Starting to parse PHYML posterior probas...");
            while ((line=br.readLine())!=null) {
                lineNumber++;
                if (lineNumber%500000==0) {
                    Infos.println("Line: >"+lineNumber);
                }
                
                if (line.startsWith("Site\tNode")) { //start of section
                    String[] stateStrings=line.split("\t");
                    stateOrder=new char[stateStrings.length-2]; //2 first columns (site/node)
                    for (int i = 0; i < stateOrder.length; i++) {
                        stateOrder[i]=stateStrings[i+2].charAt(0);
                    }
                    System.out.println("States found:"+Arrays.toString(stateOrder));
                    for (byte i = 0; i < states.getStateCount(); i++) {
                        System.out.println(i+":"+states.byteToState(i));
                    }                    
                    start=true;
                    continue;
                }
                if (line.trim().equals("")) { //useless lines
                    continue;
                }

                if (start) {

                    //groups:
                    //Site\tNode\tA\tC\tG\tT                    
                    String[] data=line.split("\t");

                    //node parsing
                    currentNode=data[1].trim();
                    nodeId=tree.getByName("x"+currentNode).getId();
                    if (nodeId!=previousId) {
                        nodeCount++;
                        //DEBUG
                        if (nodeCount>debugNodeLimit) {
                            break;
                        //DEBUG
                        }
                    }
                    previousId=nodeId;
                    
                    //site and probas parsing
                    int site=Integer.parseInt(data[0].trim());
                    for (int i = 0; i < stateOrder.length; i++) {
                        char s = stateOrder[i];
                        SiteProba sp=new SiteProba();
                        sp.state=states.stateToByte(s);
                        sp.proba=Float.parseFloat(data[i+2]);
                        if (sp.proba<sitePPThreshold)
                            sp.proba=sitePPThreshold;
                        if (asLog10)
                            sp.proba=(float)Math.log10(sp.proba);
                        probasPerSite.set(i,sp);
                    }
                    //sort 
                    Collections.sort(probasPerSite);
                    System.out.println(probasPerSite);
                    //register
                    matrix.setStates(nodeId, site-1, probasPerSite);
                    
                        
                    currentPP++;
                }
            }
            Infos.println( "Number of (site x nodes) for which pp were parsed: "+currentPP);
            Infos.println( "Number of (sites) for which pp were parsed: "+(0.0+currentPP/(tree.getNodeCount()-tree.getLeavesCount())));
            
        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (java.lang.NumberFormatException ex) {
            Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states.getStateCount()));                    
            ex.printStackTrace();
            System.exit(1);
        } catch (IOException ex) {
            Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return matrix;
    }
    
    
    
    public static void main(String[] args) {
        
        try {
            States s=new DNAStates();
            File a=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/JC69_based_comparison/phyml/basic.fasta");
            File t=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/JC69_based_comparison/phyml/basic.phylip_phyml_tree.txt");
            File p=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/JC69_based_comparison/phyml/basic.phylip_phyml_ancestral_seq");
            
            //////////////////////
            //LOAD ORIGINAL ALIGNMENT
            FASTAPointer fp=new FASTAPointer(a, false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(fastas);
            Infos.println(align.describeAlignment(false));
            fp.closePointer();
            
            
            //tree
            PHYMLWrapper pw=new PHYMLWrapper(align,s);
            PhyloTree tree=pw.parseTree(new FileInputStream(t));
            //tree.displayTree();
            //probas
            PProbasSorted probas = pw.parseSortedProbas(new FileInputStream(p),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            System.out.println(probas.getStateCount());
            System.out.println(probas.getSiteCount());

            int nodeId=tree.getByName("x5").getId();
            System.out.println("nodeId:"+nodeId );
            
            for (int i = 0; i < probas.getSiteCount(); i++) {
                for (int j = 0; j < probas.getStateCount(); j++) {
                    
                    System.out.print(s.byteToState(probas.getState(nodeId, i, j))+"="+probas.getPP(nodeId, i, j));
                    System.out.print("\t");
                    
                }
                System.out.println("");
            }
            System.out.println("----------------------------------------");

            System.out.println("");
            for (int i = 0; i < probas.getSiteCount(); i++) {
                System.out.print(i+":");
                for (int j = 0; j < probas.getStateCount(); j++) {
                    System.out.print(s.byteToState(probas.getState(nodeId, i, j))+"="+Math.pow(10,probas.getPP(tree.getByName("x5").getId(), i, j)));
                    System.out.print("\t");
                    
                }
                System.out.println("");
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(PHYMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(PHYMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }

    
}
