package misc;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
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
public class DoOrganismXML {
    
    public static void main(String[] args) {
        
        args=new String[1];
        args[0]="/home/ben/Desktop/Genomes_defined_for_orthoinspector.csv";
        
        try {
            FileWriter fw = new FileWriter(new File("/media/ben/STOCK/DATA/NCBI_VIROMES/Poxviridae_sources/data_for_orthoinspector_selection_of_best_loci/organism.xml"));
            BufferedWriter bw =new BufferedWriter(fw);
            
            bw.append( "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
            bw.append( "<orthoinspector>\n");
            
            FileReader fr=new FileReader(args[0]);
            BufferedReader br=new BufferedReader(fr);
            
            String line=null;
            while ((line=br.readLine())!=null) {                
                String[] infos = line.split("\t"); //identifier taxid genusspecies taxonomy
                
                bw.append( "<organism>\n" );
                String genus = infos[2].split(" ")[0];
                String species = infos[2].substring(infos[2].split(" ")[0].length()+1);
                bw.append( "<genus>"+genus+"</genus>\n");
                bw.append( "<species>"+species+"</species>\n");
                bw.append( "<identifier>"+infos[0]+"</identifier>\n");
                bw.append( "<phylum>"+infos[3]+"</phylum>\n");
                bw.append( "<source>\n" +
                            "<bank>\n" +
                            "<bk_date>2016-007-01</bk_date>\n" +
                            "<bk_identifier>Genebank</bk_identifier>\n" +
                            "<bk_description>Genebank</bk_description>\n" +
                            "<bk_taxinomy>\n" +
                            "<taxid>"+infos[1]+"</taxid>\n" +
                            "</bk_taxinomy>\n" +
                            "</bank>\n" +
                            "</source>\n");
                bw.append(  "<fastafile>"+infos[0]+".fasta</fastafile>\n");
                bw.append( "</organism>\n");    
                        
            }
            bw.append( "</orthoinspector>\n");
            
            bw.close();
            fw.close();
            
        } catch (IOException ex) {
            Logger.getLogger(DoOrganismXML.class.getName()).log(Level.SEVERE, null, ex);
        }
        
   
    }
    
    
}
