/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.PProbasSorted;
import core.States;
import etc.Infos;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.ExtendedTree;
import tree.PhyloTree;

/**
 * class allowing to switch between software inputs that are used for AR ; the
 * AR softwares are changing the labels of internal nodes or even tips, here 
 * the class makes the correspondence between extended tree and modified 
 * extended trees.
 * Note that it's the ExtendedTree object itself that already contains
 * the correspondence between original PhyloTree and ExtendedTree.
 * @author ben
 */
public class ARResults {
    
    //the parsed data
    private States s=null;
    //original extended alignment
    private Alignment extendedAlign=null;
    //original tree (before Extension)
    private PhyloTree originalTree=null;
    //original extended tree (before AR)
    private ExtendedTree extendedTree=null;
    //extended tree, after AR. Labels modified by AR software
    private PhyloTree ARTree=null;
    //posterior probas associated to AR
    private PProbasSorted probas=null;
    //mapping of the node id before and after AR. map(ExtendedTree NodeID)=ARTree NodeID
    private HashMap<Integer, Integer> ARTreeToExtendedTreeNodeMapping = null;
    //Launcher used to perfomr the AR
    private ARProcessLauncher arpl=null;
    
    /**
     * non-associated manager, needs to be associated to sources through @associate()
     */
    public ARResults() {}
    
    /**
     * determines the source registered to this manager at instantiation
     * @param source the source, as defined in the class ARProcessLauncher ; one of ARProcessLauncher.AR_???
     * @param alignment
     * @param Tree
     * @param probas
     */
    public ARResults(ARProcessLauncher arpl, Alignment extendedAlign, PhyloTree originalTree, ExtendedTree extendedTree, States s) {
        this.arpl=arpl;
        this.s=s;
        this.extendedAlign=extendedAlign;
        this.originalTree=originalTree;
        this.extendedTree=extendedTree;
        //step 1: parse AR results
        try {
            parseResults();
        } catch (IOException ex) {
            Logger.getLogger(ARResults.class.getName()).log(Level.SEVERE, null, ex);
        }
        //step 2: Matches ARTree node labels and
        //takes into account if the AR software unrooted/rooted things
        ARTreeToExtendedTreeNodeMapping = ARTree.mapNodes(extendedTree);
    }
    
    
    public Alignment getExtendedAlignment() {
        return extendedAlign;
    }
    
    public PProbasSorted getPProbas() {
        return probas;
    }
    
    public PhyloTree getOriginalTree() {
        return originalTree;
    }
    
    public ExtendedTree getExtendedTree() {
        return extendedTree;
    }
    
    public PhyloTree getARTree() {
        return ARTree;
    }

    /**
     * map(ARTreeNodeID)=ExtendedTreeNodeID
     * @return 
     */
    public HashMap<Integer,Integer> getTreeMapping() {
        return ARTreeToExtendedTreeNodeMapping;
    }
    
    
    public int getOriginalExtendedTreeNodeId(int ARExtencedTreeNodeId) {
        return ARTreeToExtendedTreeNodeMapping.get(ARExtencedTreeNodeId);
    }
    
    
    private void parseResults() throws IOException {

        if (this.arpl.currentProg==ARProcessLauncher.AR_PAML) {
            //with PAML both tree and probas are in the rst file
            File rst = new File(arpl.ARPath.getAbsolutePath()+File.separator+"rst");
            long startTime = System.currentTimeMillis();
            PAMLWrapper pw=new PAMLWrapper(extendedAlign,s);
            this.ARTree=pw.parseTree(new FileInputStream(rst));
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of PAML modified tree used " + (endTime - startTime) + " ms");
            //probas
            startTime = System.currentTimeMillis();
            this.probas = pw.parseSortedProbas(new FileInputStream(rst),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            endTime = System.currentTimeMillis();
            Infos.println("Loading of PAML Posterior Probas used " + (endTime - startTime) + " ms");
        } else if (this.arpl.currentProg==ARProcessLauncher.AR_PHYML) {
            //with PHYML both tree and probas files are alignment name + extension
            File tree = new File(arpl.ARPath.getAbsolutePath()+File.separator+arpl.alignPath.getName()+"_phyml_ancestral_tree.txt");
            long startTime = System.currentTimeMillis();
            PHYMLWrapper pw=new PHYMLWrapper(extendedAlign,s);
            this.ARTree=pw.parseTree(new FileInputStream(tree),false);
            Infos.println("AR tree rooted? :"+ARTree.isRooted());
            //phyml AR unroot input tree even if they are rooted
            //we need to reroot the tree as is was in the inital input
            //What doeos phyml ?
            //input:        ((C1,C2)node,S3)root;
            //phyml output: (C3,C1,C2)newick_root;
            //then we go back to the original input configuration
            this.ARTree=pw.parseTree(new FileInputStream(tree),true);
            
            
            Infos.println("PhyML unroots trees, rerooting this tree. AR tree rooted? :"+ARTree.isRooted());
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of PHYML modified tree used " + (endTime - startTime) + " ms");
            //probas
            File align = new File(arpl.ARPath.getAbsolutePath()+File.separator+arpl.alignPath.getName()+"_phyml_ancestral_seq.txt");
            startTime = System.currentTimeMillis();
            this.probas = pw.parseSortedProbas(new FileInputStream(align),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            endTime = System.currentTimeMillis();
            Infos.println("Loading of PHYML Postrerior Probas used " + (endTime - startTime) + " ms");
            
        } else if (this.arpl.currentProg==ARProcessLauncher.AR_FASTML) {
            //probas
            long startTime = System.currentTimeMillis();
            //!!!!!!!!!!!!!!!
            //this.probas = fw..parseProbas(new FileInputStream(probasFile),1e-6f,false);
            Infos.println("return of  PProbasSorted not implemented yet for FASTML !");
            System.exit(1);
            //!!!!!!!!!!!!!!!!
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of probas used " + (endTime - startTime) + " ms");
        } 

        
    }
    
    
    
}
