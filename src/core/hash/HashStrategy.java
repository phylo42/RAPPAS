/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import it.unimi.dsi.fastutil.Hash;
import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public class HashStrategy implements Hash.Strategy<byte[]>, Serializable{
    
    private static final long serialVersionUID = 7001L;

    public HashStrategy() {}
    
    @Override
    public int hashCode(byte[] object) {
        return Arrays.hashCode(object);
    }
    @Override
    public boolean equals(byte[] o1, byte[] o2) {
        return Arrays.equals(o1, o2);
    }
                            
}
