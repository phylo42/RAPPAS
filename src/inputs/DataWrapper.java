/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import core.PProbas;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import tree.PhyloNode;

/**
 *
 * @author ben
 */
public interface DataWrapper {

    /**
     * a wrapper load the posterior probas from the output of any ancestral state reconstruction software
     * @param input
     * @return a matrix @PPStats of posterior probas associated to nodes and sites
     * @throws FileNotFoundException
     * @throws IOException 
     */    
    public PProbas parseProbas(InputStream input, double sitePPThreshold) throws IOException;
}
