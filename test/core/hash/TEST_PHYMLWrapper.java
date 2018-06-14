/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import alignement.Alignment;
import core.DNAStates;
import core.PProbasSorted;
import core.States;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import inputs.PHYMLWrapper;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class TEST_PHYMLWrapper {

    
    public static void main(String[] args) {
        
        try {
            States s=new DNAStates();
            File a=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/JC69_based_comparison/phyml/basic.fasta");
            File t=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/JC69_based_comparison/phyml/basic.phylip_phyml_tree.txt");
            File p=new File("/media/ben/STOCK/DATA/ancestral_reconstruct_tests/JC69_based_comparison/phyml/basic.phylip_phyml_ancestral_seq");
//              File a=new File("/home/yann/Bureau/Data/test_220518/extended_align.fasta");
//              File t=new File("/home/yann/Bureau/Data/test_220518/extended_align.phylip_phyml_ancestral_tree.txt");
//              File p=new File("/home/yann/Bureau/Data/test_220518/extended_align.phylip_phyml_ancestral_seq.txt");
            
            //////////////////////
            //LOAD ORIGINAL ALIGNMENT
            FASTAPointer fp=new FASTAPointer(a, false);
            Fasta fasta=null;
            ArrayList<Fasta> fastas=new ArrayList<>();
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fastas.add(fasta);
            }
            Alignment align=new Alignment(s,fastas);
            Infos.println(align.describeAlignment(false));
            fp.closePointer();
            
            
            //tree
            PHYMLWrapper pw=new PHYMLWrapper(align,s);
            PhyloTree tree=pw.parseTree(new FileInputStream(t),false);
            //tree.displayTree();
            //probas
            PProbasSorted probas = pw.parseSortedProbas(new FileInputStream(p),Float.MIN_VALUE,true,Integer.MAX_VALUE);
            System.out.println(probas.getStateCount());
            //System.out.println(probas.getState(0, 0, 3)+" ; "+probas.getPP(0, 0, 3));
            System.out.println(probas.getSiteCount());

            int nodeId=tree.getByName("x5").getId();
            //int nodeId=tree.getByName("38").getId();
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
                    //System.out.print(s.byteToState(probas.getState(nodeId, i, j))+"="+Math.pow(10,probas.getPP(tree.getByName("38").getId(), i, j)));
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
