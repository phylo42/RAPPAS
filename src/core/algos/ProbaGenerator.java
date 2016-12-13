/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

/**
 *
 * @author ben
 */
public class ProbaGenerator {

    double[] wordProbas=null;

    
    public static double[] generateProbas(byte[][] words, double[][] ppSet) {
        double[] probas=new double[words.length];
        double val=1.0;
        for (int i = 0; i < words.length; i++) {
            val=1.0;
            for (int site = 0; site < words[i].length; site++) {
                val*=ppSet[site][words[i][site]];
            }
            probas[i]=val;
        }

        return probas;
    }
    
    
    
}
