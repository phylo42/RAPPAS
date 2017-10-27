/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.Infos;
import java.io.Serializable;
import java.util.LinkedHashMap;

/**
 *
 * @author ben
 */
public class AAStates extends AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6002L;
    
    
    char[] states = {   'R','H','K','D','E',
                        'S','T','N','Q','C',
                        'G','P','A','I','L',
                        'M','F','W','Y','V',
                        '-'
                    };
    byte[] bytes = {    
                        (byte)0x00,(byte)0x01,(byte)0x02,(byte)0x03,(byte)0x04,
                        (byte)0x05,(byte)0x06,(byte)0x07,(byte)0x08,(byte)0x09,
                        (byte)0x0A,(byte)0x0B,(byte)0x0C,(byte)0x0D,(byte)0x0E,
                        (byte)0x0F,(byte)0x10,(byte)0x11,(byte)0x12,(byte)0x13,
                        (byte)0x14
                    };
    
    LinkedHashMap<Byte,Character> s =new LinkedHashMap<>();
    LinkedHashMap<Character, Byte> b =new LinkedHashMap<>();
    
    public AAStates() {
       
        //positives
        s.put((byte)0, 'R');
        s.put((byte)1, 'H');
        s.put((byte)2, 'K');
        //Negatives
        s.put((byte)3, 'D');
        s.put((byte)4, 'E');
        //Polar uncharged
        s.put((byte)5, 'S');
        s.put((byte)6, 'T');
        s.put((byte)7, 'N');
        s.put((byte)8, 'Q');
        //Other
        s.put((byte)9, 'C');
        s.put((byte)10, 'G');
        s.put((byte)11, 'P');
        //Hydrophobic
        s.put((byte)12, 'A');
        s.put((byte)13, 'I');
        s.put((byte)14, 'L');
        s.put((byte)15, 'M');
        s.put((byte)16, 'F');
        s.put((byte)17, 'W');
        s.put((byte)18, 'Y');
        s.put((byte)19, 'V');
        
        b.put('R',(byte)0);
        b.put('H',(byte)1);
        b.put('K',(byte)2);
        b.put('D',(byte)3);
        b.put('E',(byte)4);
        b.put('S',(byte)5);
        b.put('T',(byte)6);
        b.put('N',(byte)7);
        b.put('Q',(byte)8);
        b.put('C',(byte)9);
        b.put('G',(byte)10);
        b.put('P',(byte)11);
        b.put('A',(byte)12);
        b.put('I',(byte)13);
        b.put('L',(byte)14);
        b.put('M',(byte)15);
        b.put('F',(byte)16);
        b.put('W',(byte)17);
        b.put('Y',(byte)18);
        b.put('V',(byte)19);
        //ambigous states which are allowed
        ambigousStates=1;
        s.put((byte)20, '-');
        b.put('-',(byte)20);
        
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
    public byte stateToByte(char c) {
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
