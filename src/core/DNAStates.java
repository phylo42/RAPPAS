/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.Infos;
import etc.exceptions.NonIUPACStateException;
import java.io.Serializable;
import java.util.HashMap;

/**
 * simple object describing the conversion between bytes and residues
 * @author ben
 */
public class DNAStates extends AbstractStates implements Serializable {
    
    private static final long serialVersionUID = 6001L;
    
    char[]states = {'A','T','C','G','N','-'};
    HashMap<Character,Boolean> convertedToN=new HashMap<>();
    byte[] bytes = {(byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03,(byte)0x04,(byte)0x05}; 

    public DNAStates() {
        //ambigous states which are allowed
        ambigousStates=2;
        //fill hashmap that correspond to IUPAC code
        convertedToN.put('R', Boolean.TRUE);convertedToN.put('r', Boolean.TRUE);
        convertedToN.put('Y', Boolean.TRUE);convertedToN.put('y', Boolean.TRUE);
        convertedToN.put('S', Boolean.TRUE);convertedToN.put('s', Boolean.TRUE);
        convertedToN.put('W', Boolean.TRUE);convertedToN.put('w', Boolean.TRUE);
        convertedToN.put('K', Boolean.TRUE);convertedToN.put('k', Boolean.TRUE);
        convertedToN.put('M', Boolean.TRUE);convertedToN.put('m', Boolean.TRUE);
        convertedToN.put('B', Boolean.TRUE);convertedToN.put('b', Boolean.TRUE);
        convertedToN.put('D', Boolean.TRUE);convertedToN.put('d', Boolean.TRUE);
        convertedToN.put('H', Boolean.TRUE);convertedToN.put('h', Boolean.TRUE);
        convertedToN.put('V', Boolean.TRUE);convertedToN.put('v', Boolean.TRUE);
    }

    @Override
    protected byte charToByte(char c) throws NonIUPACStateException {
        byte b=-1;
        switch (c) {
            case 'a':
                b=0x00; break;
            case 'A':
                b=0x00; break;
            case 't':
                b=0x01; break;
            case 'T':
                b=0x01; break;
            case 'c':
                b=0x02; break;
            case 'C':
                b=0x02; break;
            case 'g':
                b=0x03; break;
            case 'G':
                b=0x03; break;
            case 'N':
                b=0x04; break;
            case '-':
                b=0x05; break;
            default:
                if (convertedToN.containsKey(c)) {
                    Infos.println("Unexpected IUPAC state in the sequence (not ATUCGN-.), replaced with N. (char='"+String.valueOf(c)+"')");
                    b=0x04;
                } else {
                    throw new NonIUPACStateException(this, c);
                }
                break; //put N if other IUPAC base
        }
        return b;
    }
    
    
    
    @Override
    public char byteToState(byte b) {
        return states[b];
    }
    
    @Override
    public byte stateToByte(char c) throws NonIUPACStateException {
        if (c=='U') {
            c='T';
        }
        try {
            return bytes[charToByte(c)];
        } catch (NullPointerException ex) {
            //ex.printStackTrace();
            Infos.println("Unexpected (not ATUCGN) state in the sequence, replaced with N. (char='"+String.valueOf(c)+"')");
        }
        return 'N';
    }
    

    @Override
    public String getSequence(byte[] bytes) {
        char[] c=new char[bytes.length];
        for (int i = 0; i < c.length; i++) {
            c[i]=states[bytes[i]];
        }
        return new String(c);
    }

    @Override
    public int getStateCount() {
        return states.length;
    }
       
    @Override
    public int getNonAmbiguousStatesCount() {
        return states.length-ambigousStates;
    }

    @Override
    public int stateToInt(char c) throws NonIUPACStateException {
        if (c=='U') {
            c='T';
        }
        return bytes[charToByte(c)];
    }

    @Override
    public byte[] compressMer(byte[] bytes) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public char[] expandMer(byte[] mer, int k) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }


}
