/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ben
 */
public class FASTQPointer implements SequencePointer {

private StringBuffer sb=new StringBuffer();
    boolean markPositionned=false;
    String lastHeader="";
    File myFile;
    LineNumberReader lnr=null;
    boolean first=true;
    boolean last=false;
    boolean sequence=false; //by opposition to quality section

    int currentLine=-1;

    /**
     * replaced characters during the parsing
     */
    public static final String[] eliminatedCharacters= {};

    //nombre d'entete fasta lues
    int fastqChecked=0;
    //limiteHaute du nombre de fastas lus
    int limitHigh=0;

    //nombre de fastas dans le fichier
    int size=0;

    public FASTQPointer (File f) {
        this.myFile=f;
        //ouverture du reader si pas déjà fait
        try {
            size=checkSize();
            lnr = new LineNumberReader(new FileReader(myFile));
        } catch (IOException ex) {
            Logger.getLogger(FASTAPointer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * return the next fasta
     * @return
     */
    @Override
    public StringBuffer nextSequenceAsFasta() {
        try {

            sb = new StringBuffer();
            //ajout de l'header repéré à la précédante itération
            sb.append(lastHeader);
            
            String line;
            //positionnement à la mark, si elle a déjà été dèf
            if (markPositionned) {
                lnr.setLineNumber(currentLine);
            }
            // lit le fichier normalement
            String lineRead = "";
            while ((lineRead = lnr.readLine()) != null) {
                if (lineRead.isEmpty()) {
                    continue;
                }
                if (lineRead.equals("\n") | lineRead.startsWith("#")) {
                    continue;
                }
                //elimination des character spéciaux
                /*for (String s:eliminatedCharacters) {
                    lineRead=lineRead.replaceAll(s, ":");
                }*/
                //si c une ligne avec un '>'
                if (lineRead.charAt(0)=='@' & first == false) {
                    fastqChecked++;
                    lastHeader = ">"+lineRead.substring(1) + "\n";
                    sequence=true;
                    return sb;
                }
                // Mark the position in the file.
                currentLine=lnr.getLineNumber();
                //si c le premier '>'
                if (lineRead.startsWith("@") & first == true) {
                    fastqChecked++;
                    first = false;
                    sb.append(">"+lineRead.substring(1) + "\n");
                    sequence=true;
                    continue;
                }
                if (lineRead.startsWith("+")) {
                    sequence=false;
                }
                
                //ajoute la ligne au stringbuffer
                if (sequence) {
                    sb.append(lineRead + "\n");
                }
                lineRead=null;
            }
            if (lineRead==null & last==false) {
                last=true;
                return sb;
            }

            lnr.close();
            return null;
        } catch (IOException ex) {
            Logger.getLogger(FASTAPointer.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }

     }
    
    /**
     *
     * @return
     */
    public Fasta nextProteinAsFastaObject() {
        StringBuffer seq=nextSequenceAsFasta();
        if (seq!=null) {
            int linereturnIndex=seq.indexOf("\n");
            return new Fasta(seq.substring(1, linereturnIndex).toString(), seq.substring(linereturnIndex).trim());
        } else {
            return null;
        }
    }

    

    public void setlimitLow(int limitLow) {
        try {
            int n = 0;
            String lineRead = "";
            while ((lineRead = lnr.readLine()) != null) {
                if (lineRead.startsWith(">")) {
                    n++;
                    if (n == limitLow) {
                        lnr.mark((int) myFile.length());
                        markPositionned=true;
                    }
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(FASTAPointer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void setLimitHigh(int limitHigh) {
        this.limitHigh=limitHigh;
    }

    public void resetPointer() {
        try {
            sb=new StringBuffer();
            markPositionned=false;
            lastHeader="";
            lnr.close();
            lnr=null;
            size = 0;
            first=true;
            last=false;
            lnr = new LineNumberReader(new FileReader(myFile));
            fastqChecked=0;
        } catch (IOException ex) {
            Logger.getLogger(FASTAPointer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    /**
     * return the number of fasta contained in the associated file
     * @return
     */
    @Override
    public int getContentSize() {
        return this.size;
    }

    /**
     * call it at the end of the parsing
     */
    @Override
    public void closePointer() {
        try {
            this.lnr.close();
            this.sb=null;
        } catch (IOException ex) {
            Logger.getLogger(FASTAPointer.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private int checkSize() throws FileNotFoundException, IOException {
        int n=0;
        String lineRead = "";
        BufferedReader buf=new BufferedReader(new FileReader(myFile));
        while ((lineRead = buf.readLine()) != null) {
            if (lineRead.isEmpty()) {
                continue;
            }
            if (lineRead.charAt(0)=='>') {
                n++;
            }
        }
        buf.close();
        buf=null;
        return n;
    }

    public void setPointerPosition(int fastaNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getContentMean() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Fasta nextSequenceAsFastaObject() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
