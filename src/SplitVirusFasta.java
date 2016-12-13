
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import inputs.FASTAPointer;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ben
 */
public class SplitVirusFasta {
    
    public static void main(String[] args) {
        
        args=new String[1];
        args[0]= "/media/ben/STOCK/DATA/NCBI_VIROMES/Poxviridae_sources/all_translations.fasta";

        try {
            System.out.println("Split a VIPR fasta by species. (First id before SPACE)");
            System.out.println("A unique id is added before the seq description. ");

            FASTAPointer fpp=new FASTAPointer(new File(args[0]),true);
            
            //will check the size to confirm that the last id attributed to the fasta is matching.
            int size=fpp.getContentSize();
            fpp.resetPointer();
            
            String previousGenomeIdentifier="";
            String currentGenomeIdentifier="";
            
            StringBuffer fastaEntry=null;
            int counterTotal=0;
            int counterOrganism=0;
            FileWriter fw=null;
            
            while ((fastaEntry=fpp.nextSequenceAsFasta())!=null) {
                                
                currentGenomeIdentifier=fastaEntry.substring(0, fastaEntry.indexOf("\n")).split(" ")[0].substring(1);
                counterTotal++;
                
                if (previousGenomeIdentifier.equals("")) { //1st file
                    fw=new FileWriter("/media/ben/STOCK/DATA/NCBI_VIROMES/Poxviridae_sources/test/"+currentGenomeIdentifier+".fasta");
                } else if (!previousGenomeIdentifier.equals(currentGenomeIdentifier) ) {
                    //close previous organism file
                    fw.close();
                    System.out.println(counterOrganism+" seqs added to "+previousGenomeIdentifier+".fasta");
                    counterOrganism=0;

                    //build new organism file
                    fw=new FileWriter("/media/ben/STOCK/DATA/NCBI_VIROMES/Poxviridae_sources/test/"+currentGenomeIdentifier+".fasta");
                }
                
                fw.append(fastaEntry.replace(0, currentGenomeIdentifier.length()+1, ">"+currentGenomeIdentifier+"|"+(counterOrganism+1)).toString().replaceAll(" \\(modified\\)", ""));
                previousGenomeIdentifier=currentGenomeIdentifier;
                counterOrganism++;
            }
            
            fw.close();
            
            
            System.out.println("Total fasta parsed: "+counterTotal);
            
            
            fpp.closePointer();
            
            
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SplitVirusFasta.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SplitVirusFasta.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
                
        
        
    }
}
