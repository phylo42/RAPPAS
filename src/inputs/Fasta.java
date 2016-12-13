/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

/**
 * Fasta representation, not that the '' is not in the header
 * @author ben
 */
public class Fasta {
    
    String sequence=null;
    String header=null;

    public Fasta(String header, String sequence) {
        this.sequence=sequence;
        this.header=header;
    }

    public String getHeader() {
        return header;
    }

    public String getSequence() {
        return sequence;
    }

    public void setHeader(String header) {
        this.header = header;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    
    public String getFormatedFasta() {
        return ">"+header+"\n"+sequence;
    }
    
    
    
}
