/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.older.PProbas;
import core.PProbasSorted;
import core.States;
import etc.Infos;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.NewickReader;
import tree.PhyloTree;
import tree.Tree;

/**
 * class allowing to switch between software inputs that are a priori used\n
 * for the ancestral reconstruction
 * @author ben
 */
public class InputManagerNext {
    
    public static final int SOURCE_PAML=1;
    public static final int SOURCE_FASTML=2;
    public static final int SOURCE_PHYML=3;
    
    private int currentActiveSource=-1;
    
    //the parsed data
    private States s=null;
    private Alignment align=null;
    private Tree tree=null;
    private PProbasSorted probas=null;
    
    
    
    /**
     * non-associated manager, needs to be associated to sources through @associate()
     */
    public InputManagerNext() {}
    
    /**
     * determines the source registered to this manager at instantiation
     * @param source the source, as defined by the class SOURCE_* variables
     * @param alignment
     * @param Tree
     * @param probas
     */
    public InputManagerNext(int source,File alignmentFile, File treeFile, File probasFile, States s) {
        this.s=s;
        this.currentActiveSource=source;
        try {
            parse(alignmentFile, treeFile, probasFile);
        } catch (IOException ex) {
            Logger.getLogger(InputManagerNext.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
     * associate sources if not done at instantiation
     * @param source
     * @param alignment
     * @param Tree
     * @param probas
     */
    public void associate(int source,File alignmentFile, File treeFile, File probasFile, States s) {
        this.s=s;
        this.currentActiveSource=source;
        try {
            parse(alignmentFile, treeFile, probasFile);
        } catch (IOException ex) {
            Logger.getLogger(InputManagerNext.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    
    public Alignment getAlignment() {
        return align;
    }
    
    public PProbasSorted getPProbas() {
        return probas;
    }
    
    public Tree getTree() {
        return tree;
    }
    
    private void parse(File alignmentFile, File treeFile, File probasFile) throws IOException {

        if (this.currentActiveSource==SOURCE_PAML) {
            //alignment
            long startTime = System.currentTimeMillis();
            FASTAPointer fpp=new FASTAPointer(alignmentFile,false);
            Fasta f=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((f=fpp.nextSequenceAsFastaObject())!=null) {
                fastas.add(f);
            }
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of alignment used " + (endTime - startTime) + " ms");
            this.align=new Alignment(fastas);
            //tree
            startTime = System.currentTimeMillis();
            PAMLWrapper pw=new PAMLWrapper(align,s);
            if (treeFile==null)
                this.tree=pw.parseTree(new FileInputStream(probasFile));
            else
                this.tree=pw.parseTree(new FileInputStream(treeFile));
            endTime = System.currentTimeMillis();
            Infos.println("Loading of tree used " + (endTime - startTime) + " ms");
            //probas
            startTime = System.currentTimeMillis();
            this.probas = pw.parseSortedProbas(new FileInputStream(probasFile),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            endTime = System.currentTimeMillis();
            Infos.println("Loading of probas used " + (endTime - startTime) + " ms");
        } else if (this.currentActiveSource==SOURCE_PHYML) {
            //alignment
            long startTime = System.currentTimeMillis();
            FASTAPointer fpp=new FASTAPointer(alignmentFile,false);
            Fasta f=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((f=fpp.nextSequenceAsFastaObject())!=null) {
                fastas.add(f);
            }
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of alignment used " + (endTime - startTime) + " ms");
            this.align=new Alignment(fastas);
            //tree
            startTime = System.currentTimeMillis();
            PHYMLWrapper pw=new PHYMLWrapper(align,s);
            this.tree=pw.parseTree(new FileInputStream(treeFile));
            endTime = System.currentTimeMillis();
            Infos.println("Loading of tree used " + (endTime - startTime) + " ms");
            //probas
            startTime = System.currentTimeMillis();
            this.probas = pw.parseSortedProbas(new FileInputStream(probasFile),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            endTime = System.currentTimeMillis();
            Infos.println("Loading of probas used " + (endTime - startTime) + " ms");
        } else if (this.currentActiveSource==SOURCE_FASTML) {
            //alignment
            long startTime = System.currentTimeMillis();
            FASTAPointer fpp=new FASTAPointer(alignmentFile,false);
            Fasta f=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((f=fpp.nextSequenceAsFastaObject())!=null) {
                fastas.add(f);
            }
            this.align=new Alignment(fastas);
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of alignment used " + (endTime - startTime) + " ms");
            //tree
            startTime = System.currentTimeMillis();
            BufferedReader br=new BufferedReader(new FileReader(treeFile));
            String treeString=null;
            while((treeString=br.readLine())!=null) {}
            this.tree = NewickReader.parseNewickTree2(treeString);
            endTime = System.currentTimeMillis();
            Infos.println("Loading of tree used " + (endTime - startTime) + " ms");
            //probas
            startTime = System.currentTimeMillis();
            FASTMLWrapper fw=new FASTMLWrapper(align, tree, s);
            
            //!!!!!!!!!!!!!!!
            //this.probas = fw..parseProbas(new FileInputStream(probasFile),1e-6f,false);
            Infos.println("return of  PProbasSorted not implemented yet for FASTML !");
            System.exit(1);
            //!!!!!!!!!!!!!!!!
            
            endTime = System.currentTimeMillis();
            Infos.println("Loading of probas used " + (endTime - startTime) + " ms");
        } 

        
    }
    
    
    
}
