/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

/**
 *
 * @author ben
 */
public interface States {
    
    /**
     * from a byte, rturn the corresponding character
     * @param b
     * @return 
     */
    public char byteToState(byte b);

    /**
     * from a character, return the corresponding byte
     * @param c
     * @return 
     */
    public byte stateToByte(char c);
    /**
     * from a character, returns the correposnding byte, but as an int.
     * @param c
     * @return 
     */
    public int stateToInt(char c); //simply convert the byte to unsigned int
    public String getSequence(byte[] bytes);
    public int getStateCount();
    public int getNonAmbiguousStatesCount();
    
    /**
     * allows the possibility of mer compression
     * @param bytes
     * @return 
     */
    public byte[] compressMer(byte[] bytes);

    /**
     * allows the possibility of mer expansion as chars
     * @param mer
     * @param k
     * @return 
     */
    public char[] expandMer(byte[] mer, int k);
    
}
