/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.Infos;
import etc.exceptions.NonSupportedStateException;
import java.io.Serializable;

/**
 *
 * @author ben
 */
public abstract class AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6000L;
    
    protected int ambigousStatesCount=2;

    
    /**
     * important method to implement to effectively match char and byte
     * @param c
     * @return 
     * @throws etc.exceptions.NonSupportedStateException 
     */
    protected byte charToByte(char c) throws NonSupportedStateException {
        return -1;
    }

    @Override
    public byte[] compressMer(byte[] chars, byte[] bytesCompressed) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public char[] expandMer(byte[] mer, int k) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    
}
