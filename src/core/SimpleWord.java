/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.io.Serializable;
import java.util.Arrays;

/**
 * //just used for the hashes key system
 * @author ben
 */
public class SimpleWord extends AbstractWord implements Word,Serializable {

    private static final long serialVersionUID = 4023L;
    
    public SimpleWord(byte[] w) {
        this.word=w;
    }
    
    public SimpleWord(Word w) {
        this.word=Arrays.copyOf(w.getWord(), w.getWord().length);
    }

    @Override
    public String toString() {
        return Arrays.toString(this.word);
    }
    
    
    
}
