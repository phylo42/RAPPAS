/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.Infos;
import etc.exceptions.NonSupportedStateException;
import java.io.Serializable;
import java.util.HashMap;

/**
 *
 * @author ben
 */
public class DNAStatesShifted extends AbstractStates implements Serializable {
 
    private static final long serialVersionUID = 6003L;
    
    //compression of 4 DNA residues in 1 byte, 2 bits per residue.
    //if 4 bases, in 1 byte:  00|00|00|00
    //          4 positions:   3  2  1  0
    //00000011 mask for position 0
    //00001100 mask for position 1
    //00110000 mask for position 2
    //11000000 mask for position 3
    byte[] maskArray={(byte)0x03,(byte)0x0C,(byte)0x30,(byte)0xC0};
    
    char[] states = {'A','T','C','G','N','-','.'};
    byte[] bytes = {(byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03,(byte)0x04,(byte)0x05,(byte)0x06};     
    //char[] states = {'A','T','C','G'};
    //byte[] bytes = {(byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03};     
    HashMap<Character,Boolean> ambiguousState=new HashMap<>(30);

    public DNAStatesShifted() {
        //ambigous states which are allowed
        ambigousStatesCount=24;
        //fill hashmap that correspond to IUPAC code
        ambiguousState.put('R', Boolean.TRUE);ambiguousState.put('r', Boolean.TRUE);
        ambiguousState.put('Y', Boolean.TRUE);ambiguousState.put('y', Boolean.TRUE);
        ambiguousState.put('S', Boolean.TRUE);ambiguousState.put('s', Boolean.TRUE);
        ambiguousState.put('W', Boolean.TRUE);ambiguousState.put('w', Boolean.TRUE);
        ambiguousState.put('K', Boolean.TRUE);ambiguousState.put('k', Boolean.TRUE);
        ambiguousState.put('M', Boolean.TRUE);ambiguousState.put('m', Boolean.TRUE);
        ambiguousState.put('B', Boolean.TRUE);ambiguousState.put('b', Boolean.TRUE);
        ambiguousState.put('D', Boolean.TRUE);ambiguousState.put('d', Boolean.TRUE);
        ambiguousState.put('H', Boolean.TRUE);ambiguousState.put('h', Boolean.TRUE);
        ambiguousState.put('V', Boolean.TRUE);ambiguousState.put('v', Boolean.TRUE);
        ambiguousState.put('N', Boolean.TRUE);ambiguousState.put('n', Boolean.TRUE);
        ambiguousState.put('.', Boolean.TRUE);
        ambiguousState.put('-', Boolean.TRUE);
    }
    
    
    
    
    /**
     * compress DNA mer into a mer.size/4 byte array, compressing 4 residues per byte
     * @param bytes
     * @param bytesCompressed
     * @return 
     */
    @Override
    public byte[] compressMer(byte[] bytes, byte[] bytesCompressed) {
        //chose number of bytes necessary to represent this DNA
        //we can code 4 bases in 1 byte
        //int byteCount=bytesCompressed.length;
        //int byteCount=new Double(Math.ceil((0.0+bytes.length)/4)).intValue();
        //the kmer compressed in byteCount bytes
        //byte[] bytesCompressed=new byte[byteCount];
        
        byte fourBasesByte=0x00;
        int insertions=0;
        for (int i = 0; i < bytes.length; i++) {
            if ( (i>0) & (i%4==0) ) {      
                //instance copy this byte in byte array
                bytesCompressed[(i/4)-1]=new Byte(fourBasesByte).byteValue();
                //reset pivot byte
                fourBasesByte=0x00;
            }
            //shift the 2 bits coding the i-th base, in the i%4 byte
            fourBasesByte= (byte) ( fourBasesByte | (bytes[i] << (2*(i%4))) ) ;            
        }
        //last mer fragment was last than 4 bases
        //use a complete byte, and at expansion, k will be used to know how
        //many bit pairs will need to be read.
        if (insertions<bytesCompressed.length) {
            bytesCompressed[bytesCompressed.length-1]=fourBasesByte;
        }

        return bytesCompressed;
        
    }
    
    /**
     * expand compressed DNA mer into char array
     * @param mer
     * @param k
     * @return 
     */
    @Override
    public char[] expandMer(byte[] mer, int k) {
        
        char[] charMer=new char[k];
        
        for (int i = 0; i < mer.length; i++) {
            byte b = mer[i];
            for (int j = 0; j < 4; j++) {
                //to skip last unused bits
                if ( (i*4)+j > k-1) {
                    break;
                }
                int residue=b;
                residue= (byte) (b & maskArray[j]) ;
                //masking by 0x00FF necessary to avoid that the head bits of integer
                //shifted too (when shift > 5), as in java they are 1 (positive/negative interval of int) by default
                residue=residue & 0x00FF; 
                residue= (byte) (residue >> (j*2));
                charMer[(i*4)+j]=byteToState((byte)residue);
            }
        }
        return charMer;
    }

    /**
     * 
     * @param c
     * @return
     * @throws NonSupportedStateException 
     */
    @Override
    protected byte charToByte(char c) throws NonSupportedStateException {
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
            case 'u':
                b=0x01; break;
            case 'U':
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
            case 'n':
                b=0x04; break;
            case '-':
                b=0x05; break;
            case '.':
                b=0x06; break;  
            default:
                b=0x05; 
                
        }
        return b;
    }

    @Override
    public boolean isAmbiguous(char c) {
        return ambiguousState.containsKey(c);
    }

    @Override
    public char byteToState(byte b) {
        return states[b];
    }
    
    @Override
    public byte stateToByte(char c) throws NonSupportedStateException{
        
        try {
            return bytes[charToByte(c)];
        } catch (NullPointerException ex) {
            throw new NonSupportedStateException(this, c);
        }
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
        return 4;
    }

    @Override
    public int stateToInt(char c) throws NonSupportedStateException {
        return bytes[charToByte(c)];
    }
    
    
    
    
}
