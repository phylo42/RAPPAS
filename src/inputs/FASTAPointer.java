/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package inputs;

import inputs.Fasta;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * allow to parse a file with multiple blast inside
 * @author linard
 */
public class FASTAPointer implements SequencePointer {
    private StringBuffer sb=new StringBuffer();
    boolean markPositionned=false;
    String lastHeader="";
    File myFile;
    LineNumberReader lnr=null;
    boolean first=true;
    boolean last=false;
    
    boolean gapsRemoved=false;

    int currentLine=-1;

    /**
     * replaced characters during the parsing
     */
    public static final String[] eliminatedCharacters= {"\\(", "\\)", "\\,", "\'","\""};

    //nombre d'entete fasta lues
    int fastaChecked=0;
    //limiteHaute du nombre de fastas lus
    int limitHigh=0;

    //nombre de fastas dans le fichier
    int size=0;

    public FASTAPointer (File f,boolean gapsRemoved) {
        this.gapsRemoved=gapsRemoved;
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
                if (lineRead.charAt(0)=='>' & first == false) {
                    fastaChecked++;
                    lastHeader = lineRead + "\n";
                    return sb;
                }
                // Mark the position in the file.
                currentLine=lnr.getLineNumber();
                //si c le premier '>'
                if (lineRead.startsWith(">") & first == true) {
                    fastaChecked++;
                    first = false;
                    if (gapsRemoved) {
                        sb.append(lineRead.replaceAll("-", "")+"\n");
                    } else {
                        sb.append(lineRead + "\n");
                    }
                    continue;
                }
                //ajoute la ligne au stringbuffer
                if (gapsRemoved) {
                    sb.append(lineRead.replaceAll("-", "")+"\n");
                } else {
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
    public Fasta nextSequenceAsFastaObject() {
        StringBuffer seq=nextSequenceAsFasta();
        if (seq!=null) {
            int linereturnIndex=seq.indexOf("\n");
            if (gapsRemoved) {
                return new Fasta(seq.substring(1, linereturnIndex).toString(), seq.substring(linereturnIndex).replaceAll("\n", "").replaceAll("-", "").trim());
            } else {
                return new Fasta(seq.substring(1, linereturnIndex).toString(), seq.substring(linereturnIndex).replaceAll("\n", "").trim());
            }
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
            fastaChecked=0;
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
}
