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
import java.util.LinkedHashMap;

/**
 *
 * @author ben
 */
public class AAStates extends AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6002L;
    
    
    char[] states = {   
                        'R','H','K','D','E',
                        'S','T','N','Q','C',
                        'G','P','A','I','L',
                        'M','F','W','Y','V',
                        '-'
                    };
    byte[] bytes =  {    
                        (byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03,(byte)0x04,
                        (byte)0x05,(byte)0x06,(byte)0x07,(byte)0x08,(byte)0x09,
                        (byte)0x0A,(byte)0x0B,(byte)0x0C,(byte)0x0D,(byte)0x0E,
                        (byte)0x0F,(byte)0x10,(byte)0x11,(byte)0x12,(byte)0x13,
                        (byte)0x14
                    };
    
    LinkedHashMap<Byte,Character> s =new LinkedHashMap<>();
    LinkedHashMap<Character, Byte> b =new LinkedHashMap<>();
    
    HashMap<Character,Boolean> ambiguousState=new HashMap<>(30);
    
    /**
     *
     * @param convertUOX the value of convertUOX
     */
    public AAStates(boolean convertUOX) {
       
        //positives
        s.put((byte)0, 'R');s.put((byte)0, 'r');
        s.put((byte)1, 'H');s.put((byte)1, 'h');
        s.put((byte)2, 'K');s.put((byte)2, 'k');
        //Negatives
        s.put((byte)3, 'D');s.put((byte)3, 'd');
        s.put((byte)4, 'E');s.put((byte)4, 'e');
        //Polar uncharged
        s.put((byte)5, 'S');s.put((byte)5, 's');
        s.put((byte)6, 'T');s.put((byte)6, 't');
        s.put((byte)7, 'N');s.put((byte)7, 'n');
        s.put((byte)8, 'Q');s.put((byte)8, 'q');
        //Other
        s.put((byte)9, 'C');s.put((byte)9, 'c');
        s.put((byte)10, 'G');s.put((byte)10, 'g');
        s.put((byte)11, 'P');s.put((byte)11, 'p');
        //Hydrophobic
        s.put((byte)12, 'A');s.put((byte)12, 'a');
        s.put((byte)13, 'I');s.put((byte)13, 'i');
        s.put((byte)14, 'L');s.put((byte)14, 'l');
        s.put((byte)15, 'M');s.put((byte)15, 'm');
        s.put((byte)16, 'F');s.put((byte)16, 'f');
        s.put((byte)17, 'W');s.put((byte)17, 'w');
        s.put((byte)18, 'Y');s.put((byte)18, 'y');
        s.put((byte)19, 'V');s.put((byte)19, 'v');
        s.put((byte)20, '-');
        
        b.put('R',(byte)0);b.put('r',(byte)0);
        b.put('H',(byte)1);b.put('h',(byte)1);
        b.put('K',(byte)2);b.put('k',(byte)2);
        b.put('D',(byte)3);b.put('d',(byte)3);
        b.put('E',(byte)4);b.put('e',(byte)4);
        b.put('S',(byte)5);b.put('s',(byte)5);
        b.put('T',(byte)6);b.put('t',(byte)6);
        b.put('N',(byte)7);b.put('n',(byte)7);
        b.put('Q',(byte)8);b.put('q',(byte)8);
        b.put('C',(byte)9);b.put('c',(byte)9);
        b.put('G',(byte)10);b.put('g',(byte)10);
        b.put('P',(byte)11);b.put('p',(byte)11);
        b.put('A',(byte)12);b.put('a',(byte)12);
        b.put('I',(byte)13);b.put('i',(byte)13);
        b.put('L',(byte)14);b.put('l',(byte)14);
        b.put('M',(byte)15);b.put('m',(byte)15);
        b.put('F',(byte)16);b.put('f',(byte)16);
        b.put('W',(byte)17);b.put('w',(byte)17);
        b.put('Y',(byte)18);b.put('y',(byte)18);
        b.put('V',(byte)19);b.put('v',(byte)19);
        b.put('-',(byte)20);
        
        //ambigous states which are allowed
        ambigousStatesCount=2;
        //ambiguousState.put('-', true);
        ambiguousState.put('*', true); //stop codons can be present in the middle of protein translations
        ambiguousState.put('!', true); //codon containing a frameshift, used in MACSE aligner
        
        //special AA and unknown
        if (convertUOX) {
            s.put((byte)9, 'U');s.put((byte)9, 'u'); // U to C
            s.put((byte)14, 'O');s.put((byte)14, 'o'); // O to L
            s.put((byte)20, 'X');s.put((byte)20, 'x'); // X to -
            b.put('U',(byte)9);b.put('u',(byte)9);
            b.put('O',(byte)14);b.put('o',(byte)14);
            b.put('X',(byte)20);b.put('x',(byte)20);
        }
        
    }

    @Override
    public boolean isAmbiguous(char c) {
        return ambiguousState.containsKey(c);
    }
    
    
    
    

    @Override
    protected byte charToByte(char c) {
        return b.get(c);
    }

    @Override
    public char byteToState(byte b) {
        return states[b];
    }
    
    @Override
    public byte stateToByte(char c) throws NonSupportedStateException {
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
        return 20;
    }

    @Override
    public int stateToInt(char c) {
        return bytes[charToByte(c)];
    }

    @Override
    public byte[] compressMer(byte[] bytes) {
        return bytes;
    }

    @Override
    public char[] expandMer(byte[] mer, int k) {
        char[] c=new char[bytes.length];
        for (int i = 0; i < c.length; i++) {
            c[i]=states[bytes[i]];
        }
        return c;
    }
    
}
