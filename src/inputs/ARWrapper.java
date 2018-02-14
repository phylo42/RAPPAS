/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import core.PProbasSorted;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

/**
 *
 * @author ben
 */
public interface ARWrapper {

    /**
     * a wrapper load the posterior probas from the output of any ancestral state reconstruction software
     * @param input
     * @param sitePPThreshold
     * @param asLog10
     * @param debugNodeLimit
     * @return a matrix @PPStats of posterior probas associated to nodes and sites 
     */    
    public PProbasSorted parseSortedProbas(InputStream input, float sitePPThreshold, boolean asLog10, int debugNodeLimit);

    
    
}
