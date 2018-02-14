package misc;


import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class FastaAlignedToFasta {
    
    public static void main(String[] args) {
        
        
        args=new String[2];
        args[0]="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/p4z1r36_query_only.fasta";
        args[1]="/media/ben/STOCK/DATA/ancestral_reconstruct_tests/paml/pplacer_refpkg/vaginal_16s_ORIGINAL/mod_p4z1r36_query_only.fasta";
        
        
        FileWriter fw=null;
        FASTAPointer fp=null;
        try {
            System.out.println("PARAMS: file_to_convert new_file_name");
            fp=new FASTAPointer(new File(args[0]), true);
            fw = new FileWriter(new File(args[1]));
            Fasta fasta=null;
            while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
                fw.write(">"+fasta.getHeader()+"\n");
                fw.write(fasta.getSequence(false)+"\n");
            }
            
            
        } catch (IOException ex) {
            Logger.getLogger(FastaAlignedToFasta.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fp.closePointer();
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(FastaAlignedToFasta.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }
    
    
}
