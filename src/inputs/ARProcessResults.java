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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import main_v2.ARProcessLauncher;
import tree.ExtendedTree;
import tree.NewickReader;
import tree.PhyloTree;
import tree.Tree;

/**
 * class allowing to switch between software inputs that are used for AR ; the
 * AR softwares are changing the labels of internal nodes or even tips, here 
 * the class makes the correspondence between extended tree and modified 
 * extended trees.
 * Note that it's the ExtendedTree object itself that already contains
 * the correspondence between original PhyloTree and ExtendedTree.
 * @author ben
 */
public class ARProcessResults {
    
    //the parsed data
    private States s=null;
    private Alignment extendedAlign=null;
    private ExtendedTree extendedTree=null;
    private PhyloTree ARExtendedTree=null;
    private PProbasSorted probas=null;
    File ARWorkDir=null;
    
    ARProcessLauncher arpl=null;
    
    /**
     * non-associated manager, needs to be associated to sources through @associate()
     */
    public ARProcessResults() {}
    
    /**
     * determines the source registered to this manager at instantiation
     * @param source the source, as defined in the class ARProcessLauncher ; one of ARProcessLauncher.AR_???
     * @param alignment
     * @param Tree
     * @param probas
     */
    public ARProcessResults(ARProcessLauncher arpl, Alignment extendedAlign, ExtendedTree extendedTree, States s, File ARWorkDir) {
        this.arpl=arpl;
        this.s=s;
        this.extendedAlign=extendedAlign;
        this.extendedTree=extendedTree;
        this.ARWorkDir=ARWorkDir;
        try {
            parseResults();
        } catch (IOException ex) {
            Logger.getLogger(ARProcessResults.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    
    public Alignment getExtendedAlignment() {
        return extendedAlign;
    }
    
    public PProbasSorted getPProbas() {
        return probas;
    }
    
    public PhyloTree getARExtendedTree() {
        return ARExtendedTree;
    }
    
    private void parseResults() throws IOException {

        if (this.arpl.currentProg==ARProcessLauncher.AR_PAML) {
            //with PAML both tree and probas are in the rst file
            File rst = new File(ARWorkDir.getAbsolutePath()+File.separator+"rst");
            long startTime = System.currentTimeMillis();
            PAMLWrapper pw=new PAMLWrapper(extendedAlign,s);
            this.ARExtendedTree=pw.parseTree(new FileInputStream(rst));
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of PAML modified tree used " + (endTime - startTime) + " ms");
            //probas
            startTime = System.currentTimeMillis();
            this.probas = pw.parseSortedProbas(new FileInputStream(rst),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            endTime = System.currentTimeMillis();
            Infos.println("Loading of PAML Posterior Probas used " + (endTime - startTime) + " ms");
        } else if (this.arpl.currentProg==ARProcessLauncher.AR_PHYML) {
            //with PHYML both tree and probas files are alignment name + extension
            File tree = new File(ARWorkDir.getAbsolutePath()+File.separator+arpl.alignPath.getName()+"_phyml_tree.txt");
            long startTime = System.currentTimeMillis();
            PHYMLWrapper pw=new PHYMLWrapper(extendedAlign,s);
            this.ARExtendedTree=pw.parseTree(new FileInputStream(tree));
            long endTime = System.currentTimeMillis();
            Infos.println("Loading of PHYML modified tree used " + (endTime - startTime) + " ms");
            //probas
            File align = new File(ARWorkDir.getAbsolutePath()+File.separator+arpl.alignPath.getName()+"_phyml_stats.txt");
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
