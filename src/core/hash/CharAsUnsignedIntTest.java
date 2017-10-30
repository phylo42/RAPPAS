/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.hash;

/**
 *
 * @author ben
 */
public class CharAsUnsignedIntTest {
 
    public static void main(String[] args) {
        //int i=65535;
        int i=65537;
        System.out.println("i:"+i);
        char c=Character.forDigit(i, 10);
        System.out.println("c:"+c);
        char c2='\uFFFF';
        System.out.println("c2:"+c2);
        char c3=(char)i;
        int i2=(int)c3;
        int i3=(int)c2;
        System.out.println("c3:"+c3);
        System.out.println("i2:"+i2);
        System.out.println("i3:"+i3);
        
    }
    
    
}
