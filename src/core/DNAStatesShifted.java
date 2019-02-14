/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.Infos;
import etc.exceptions.NonSupportedStateException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ben
 */
public final class DNAStatesShifted extends AbstractStates implements Serializable {
 
    private static final long serialVersionUID = 6003L;
    
    //compression of 4 DNA residues in 1 byte, 2 bits per residue.
    //if 4 bases, in 1 byte:  00|00|00|00
    //          4 positions:   3  2  1  0
    //00000011 mask for position 0
    //00001100 mask for position 1
    //00110000 mask for position 2
    //11000000 mask for position 3
    byte[] maskArray={(byte)0x03,(byte)0x0C,(byte)0x30,(byte)0xC0};
    
    char[] states = {'A','T','C','G'};
    byte[] bytes = {(byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03};     
    //char[] states = {'A','T','C','G'};
    //byte[] bytes = {(byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03};     

    public DNAStatesShifted() {
        //ambigous states which are allowed
        ambigousStatesCount=24;
        //fill hashmap that correspond to IUPAC code
        ambiguousState.put('R', new byte[2]);ambiguousState.put('r', new byte[2]);
        ambiguousState.put('Y', new byte[2]);ambiguousState.put('y', new byte[2]);
        ambiguousState.put('S', new byte[2]);ambiguousState.put('s', new byte[2]);
        ambiguousState.put('W', new byte[2]);ambiguousState.put('w', new byte[2]);
        ambiguousState.put('K', new byte[2]);ambiguousState.put('k', new byte[2]);
        ambiguousState.put('M', new byte[2]);ambiguousState.put('m', new byte[2]);
        ambiguousState.put('B', new byte[3]);ambiguousState.put('b', new byte[3]);
        ambiguousState.put('D', new byte[3]);ambiguousState.put('d', new byte[3]);
        ambiguousState.put('H', new byte[3]);ambiguousState.put('h', new byte[3]);
        ambiguousState.put('V', new byte[3]);ambiguousState.put('v', new byte[3]);
        ambiguousState.put('N', new byte[4]);ambiguousState.put('n', new byte[4]);
        //ambiguousState.put('Z', new ArrayList<>());ambiguousState.put('z', new ArrayList<>()); //the rarely used Zero nucleotide
        ambiguousState.put('.', new byte[4]);
        ambiguousState.put('-', new byte[4]);
        
        try {
            //purine/pyrimide
            ambiguousState.get('R')[0]=stateToByte('A');
            ambiguousState.get('R')[1]=stateToByte('G');
            ambiguousState.get('Y')[0]=stateToByte('C');
            ambiguousState.get('Y')[1]=stateToByte('T');
            //strong/weak
            ambiguousState.get('S')[0]=stateToByte('C');
            ambiguousState.get('S')[1]=stateToByte('G');
            ambiguousState.get('W')[0]=stateToByte('A');
            ambiguousState.get('W')[1]=stateToByte('T');
            //keto/amino
            ambiguousState.get('K')[0]=stateToByte('G');
            ambiguousState.get('K')[1]=stateToByte('T');
            ambiguousState.get('M')[0]=stateToByte('A');
            ambiguousState.get('M')[1]=stateToByte('C');
            //not A
            ambiguousState.get('B')[0]=stateToByte('C');
            ambiguousState.get('B')[1]=stateToByte('G');     
            ambiguousState.get('B')[2]=stateToByte('T');
            //not C
            ambiguousState.get('D')[0]=stateToByte('A');
            ambiguousState.get('D')[1]=stateToByte('G');        
            ambiguousState.get('D')[2]=stateToByte('T');
            //not G
            ambiguousState.get('H')[0]=stateToByte('A');
            ambiguousState.get('H')[1]=stateToByte('C');   
            ambiguousState.get('H')[2]=stateToByte('T');
            //not T
            ambiguousState.get('V')[0]=stateToByte('A');
            ambiguousState.get('V')[1]=stateToByte('C');
            ambiguousState.get('V')[2]=stateToByte('G');
            //any
            ambiguousState.get('N')[0]=stateToByte('A');
            ambiguousState.get('N')[1]=stateToByte('C');
            ambiguousState.get('N')[2]=stateToByte('G');
            ambiguousState.get('N')[3]=stateToByte('T');
            ambiguousState.get('.')[0]=stateToByte('A');
            ambiguousState.get('.')[1]=stateToByte('C');         
            ambiguousState.get('.')[2]=stateToByte('G');
            ambiguousState.get('.')[3]=stateToByte('T');
            ambiguousState.get('-')[0]=stateToByte('A');
            ambiguousState.get('-')[1]=stateToByte('C');     
            ambiguousState.get('-')[2]=stateToByte('G');
            ambiguousState.get('-')[3]=stateToByte('T');
            
        } catch (NonSupportedStateException ex) {
            Logger.getLogger(DNAStatesShifted.class.getName()).log(Level.SEVERE, null, ex);
        }

        
        
    }
    
    
    
    
    /**
     * compress DNA mer into a mer.size/4 byte array, compressing 4 residues per byte
     * @param bytes
     * @return 
     */
    @Override
    public byte[] compressMer(byte[] bytes) {
        //chose number of bytes necessary to represent this DNA
        //we can code 4 bases in 1 byte
        int byteCount=new Double(Math.ceil((0.0+bytes.length)/4)).intValue();
        //the kmer compressed in byteCount bytes
        byte[] kmer=new byte[byteCount];
        
        byte fourBasesByte=0x00;
        int insertions=0;
        for (int i = 0; i < bytes.length; i++) {
            if ( (i>0) & (i%4==0) ) {      
                //instance copy this byte in byte array
                kmer[(i/4)-1]=new Byte(fourBasesByte).byteValue();
                //reset pivot byte
                fourBasesByte=0x00;
            }
            //shift the 2 bits coding the i-th base, in the i%4 byte
            fourBasesByte= (byte) ( fourBasesByte | (bytes[i] << (2*(i%4))) ) ;            
        }
        //last mer fragment was last than 4 bases
        //use a complete byte, and at expansion, k will be used to know how
        //many bit pairs will need to be read.
        if (insertions<kmer.length) {
            kmer[kmer.length-1]=fourBasesByte;
        }

        return kmer;
        
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
            default:
                throw new NonSupportedStateException(this, c);
        }
        return b;
    }

    @Override
    public boolean isAmbiguous(char c){
       return ambiguousState.containsKey(c);
    }
    
    /**
     * get list of states equivalent to this ambiguity
     * @param c
     * @return 
     * @throws etc.exceptions.NonSupportedStateException 
     */
    @Override
    public byte[] ambiguityEquivalence(char c) throws NonSupportedStateException {
        if (!ambiguousState.containsKey(c)) {
            throw new NonSupportedStateException(this, c);
        }
        return ambiguousState.get(c);
    }

    @Override
    public char byteToState(byte b) {
        return states[b];
    }
    
    @Override
    public byte stateToByte(char c) throws NonSupportedStateException{
        try {
            return bytes[charToByte(c)];
        } catch (Exception ex) {
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
