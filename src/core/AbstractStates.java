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
public abstract class AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6000L;
    
    protected LinkedHashMap<Byte, Character> states =null;
    protected LinkedHashMap<Character, Byte> bytes =null;
    protected int ambigousStates=2;

    
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


    
}
