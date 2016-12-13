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
public class AAStates extends AbstractStates implements States,Serializable {
    
    private static final long serialVersionUID = 6002L;
    

    public AAStates() {
        
        states =new LinkedHashMap<>();
        bytes =new LinkedHashMap<>();
        
        //positives
        states.put((byte)0, 'R');
        states.put((byte)1, 'H');
        states.put((byte)2, 'K');
        //Negatives
        states.put((byte)3, 'D');
        states.put((byte)4, 'E');
        //Polar uncharged
        states.put((byte)5, 'S');
        states.put((byte)6, 'T');
        states.put((byte)7, 'N');
        states.put((byte)8, 'Q');
        //Other
        states.put((byte)9, 'C');
        states.put((byte)10, 'G');
        states.put((byte)11, 'P');
        //Hydrophobic
        states.put((byte)12, 'A');
        states.put((byte)13, 'I');
        states.put((byte)14, 'L');
        states.put((byte)15, 'M');
        states.put((byte)16, 'F');
        states.put((byte)17, 'W');
        states.put((byte)18, 'Y');
        states.put((byte)19, 'V');
        states.put((byte)20, '?');
        
        bytes.put('R',(byte)0);
        bytes.put('H',(byte)1);
        bytes.put('K',(byte)2);
        bytes.put('D',(byte)3);
        bytes.put('E',(byte)4);
        bytes.put('S',(byte)5);
        bytes.put('T',(byte)6);
        bytes.put('N',(byte)7);
        bytes.put('Q',(byte)8);
        bytes.put('C',(byte)9);
        bytes.put('G',(byte)10);
        bytes.put('P',(byte)11);
        bytes.put('A',(byte)12);
        bytes.put('I',(byte)13);
        bytes.put('L',(byte)14);
        bytes.put('M',(byte)15);
        bytes.put('F',(byte)16);
        bytes.put('W',(byte)17);
        bytes.put('Y',(byte)18);
        bytes.put('V',(byte)19);
        bytes.put('?',(byte)20);
    }

}
