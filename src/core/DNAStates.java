/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 * simple object describing the conversion between bytes and residues
 * @author ben
 */
public class DNAStates extends AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6001L;
    private int ambigousStates=2;

    public DNAStates() {
        
        states =new LinkedHashMap<>();
        bytes =new LinkedHashMap<>();
    
        states.put((byte)0, 'A');
        states.put((byte)1, 'T');
        states.put((byte)2, 'C');
        states.put((byte)3, 'G');
        states.put((byte)4, 'N');
        states.put((byte)5, '-');
        bytes.put('A',(byte)0);
        bytes.put('T',(byte)1);
        bytes.put('C',(byte)2);
        bytes.put('G',(byte)3);
        //ambigous states which are allowed
        ambigousStates=2;
        bytes.put('N',(byte)4);
        bytes.put('-',(byte)5);
    }
    
    @Override
    public char byteToState(byte b) {
        return states.get(b);
    }
    
    @Override
    public byte stateToByte(char c) {
        try {
            return bytes.get(c);
        } catch (NullPointerException ex) {
            //ex.printStackTrace();
            System.out.println("Unexpected state in the sequence, replaced with N. (char='"+String.valueOf(c)+"')");
        }
        return 'N';
    }
    

    @Override
    public String getSequence(byte[] bytes) {
        char[] c=new char[bytes.length];
        for (int i = 0; i < c.length; i++) {
            c[i]=states.get(bytes[i]);
            
        }
        return new String(c);
    }

    @Override
    public int getStateCount() {
        return states.size();
    }
    
    @Override
    public int getNonAmbiguousStatesCount() {
        return states.size()-ambigousStates;
    }

    
}
