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
 *
 * @author ben
 */
public abstract class AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6000L;
    
    LinkedHashMap<Byte, Character> states =null;
    LinkedHashMap<Character, Byte> bytes =null;
    
    @Override
    public char byteToState(byte b) {
        return states.get(b);
    }
    
    @Override
    public byte stateToByte(char c) {
        return bytes.get(c);
    }

    @Override
    public int stateToInt(char c) {
        return Byte.toUnsignedInt(bytes.get(c));
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
    public byte[] getBytes(char[] sequence) {
        byte[] b=new byte[sequence.length];
        for (int i = 0; i < b.length; i++) {
            b[i]=bytes.get(sequence[i]);
            
        }
        return b;
    }

    @Override
    public int getStateCount() {
        return states.size();
    }
    
    @Override
    public int getNonAmbiguousStatesCount() {
        return states.size()-1;
    }
    
}
