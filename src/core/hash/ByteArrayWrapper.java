/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public final class ByteArrayWrapper implements Serializable {
    
    private static final long serialVersionUID = 1000001L;
    
    private byte[] data;

    public ByteArrayWrapper(byte[] data)
    {
        if (data == null)
        {
            throw new NullPointerException();
        }
        this.data = data;
    }
    
    public byte[] getArray() {
        return data;
    }
    
    public ByteArrayWrapper setArray(byte[] data) {
        this.data=data;
        return this;
    }

    @Override
    public boolean equals(Object other)
    {
        if (!(other instanceof ByteArrayWrapper))
        {
            return false;
        }
        return Arrays.equals(data, ((ByteArrayWrapper)other).data);
    }

    @Override
    public int hashCode()
    {
        return Arrays.hashCode(data);
    }
}
