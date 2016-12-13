/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

/**
 * //simple word representation through a byte[]
 * @author ben
 */
public interface Word {

    public byte[] getWord();
    
    /**
     * Here the concept is that a word can be compared through its byte[] value\n
     * Two @Word holding the same sequence of byte[] must return\n
     * when @equals() is called\n
     * Note overriding @hashcode() obliges to override equals()\n
     * If note, many functions of the JDK will behave inconsistently
     * @param obj
     * @return 
     */
    @Override
    public boolean equals(Object obj);
    /**
     * Here the concept is that a word can be compared through its byte[] value\n
     * Two @Word holding the same sequence of byte[] must return\n
     * when @equals() is called\n
     * Note: overriding @hashcode() obliges to override equals()\n
     * If note, many functions of the JDK will behave inconsistently
     * @return 
     */
    @Override
    public int hashCode();
}
