/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.Infos;
import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 *
 * @author ben
 */
public class AAStates extends AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6002L;
    

    public AAStates() {
        
        states =new LinkedHashMap<>();
        bytes =new LinkedHashMap<>();
        
        //positives
        states.put((byte)0, 'R');
        states.put((byte)1, 'H');
        states.put((byte)2, 'K');
        //Negatives
        states.put((byte)3, 'D');
        states.put((byte)4, 'E');
        //Polar uncharged
        states.put((byte)5, 'S');
        states.put((byte)6, 'T');
        states.put((byte)7, 'N');
        states.put((byte)8, 'Q');
        //Other
        states.put((byte)9, 'C');
        states.put((byte)10, 'G');
        states.put((byte)11, 'P');
        //Hydrophobic
        states.put((byte)12, 'A');
        states.put((byte)13, 'I');
        states.put((byte)14, 'L');
        states.put((byte)15, 'M');
        states.put((byte)16, 'F');
        states.put((byte)17, 'W');
        states.put((byte)18, 'Y');
        states.put((byte)19, 'V');
        
        bytes.put('R',(byte)0);
        bytes.put('H',(byte)1);
        bytes.put('K',(byte)2);
        bytes.put('D',(byte)3);
        bytes.put('E',(byte)4);
        bytes.put('S',(byte)5);
        bytes.put('T',(byte)6);
        bytes.put('N',(byte)7);
        bytes.put('Q',(byte)8);
        bytes.put('C',(byte)9);
        bytes.put('G',(byte)10);
        bytes.put('P',(byte)11);
        bytes.put('A',(byte)12);
        bytes.put('I',(byte)13);
        bytes.put('L',(byte)14);
        bytes.put('M',(byte)15);
        bytes.put('F',(byte)16);
        bytes.put('W',(byte)17);
        bytes.put('Y',(byte)18);
        bytes.put('V',(byte)19);
        //ambigous states which are allowed
        ambigousStates=2;
        states.put((byte)20, 'N');
        states.put((byte)21, '-');
        bytes.put('N',(byte)20);
        bytes.put('-',(byte)21);
        
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
            Infos.println("Unexpected state in the sequence, replaced with N. (char='"+String.valueOf(c)+"')");
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

    @Override
    public int stateToInt(char c) {
        return bytes.get(c).intValue();
    }
    
}
