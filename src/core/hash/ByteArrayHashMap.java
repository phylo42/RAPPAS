/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author ben
 */
public class ByteArrayHashMap extends HashMap<byte[], HashPointer> {

    public ByteArrayHashMap(int initalCapacity, float loadFactor) {
        super(initalCapacity, loadFactor);
    }
    
    static final int hash(Object key) {
        int h;
        return (key == null) ? 0 : (h = Arrays.hashCode((byte[])key));
    }
    
    
}
