/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.PProbasSorted;
import core.States;
import java.io.InputStream;
import tree.PhyloTree;

/**
 * 
 * @author ben
 */
public class FASTMLWrapper implements ARWrapper {

    Alignment align=null;
    PhyloTree tree=null;
    int states=4;
    
    /**
     * Used FASTML input/output files + FASTML results to fill the alignment/tree objects
     * @param inputAlign input alignment used in FASTML 
     * @param inputTree tree.newick.text output from FASTML
     * @param states
     */
    public FASTMLWrapper(Alignment inputAlign,PhyloTree inputTree,States states) {
        this.align=inputAlign;
        this.tree=inputTree;
        this.states=states.getStateCount();
    }

    @Override
    public PProbasSorted parseSortedProbas(InputStream input, float sitePPThreshold, boolean asLog10, int debugNodeLimit) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }


    
}
