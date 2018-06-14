/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import etc.exceptions.NonIUPACStateException;

/**
 *
 * @author ben
 */
public class TEST_DNAStatesShifted {

    
    public static void main(String[] args) throws NonIUPACStateException {
        
        //tests affichage
        char[] stat = {'A','T','C','G','N','-'};
        for (int i = 0; i < stat.length; i++) {
            char state = stat[i];
            System.out.println((int)state);
        }
        DNAStates states=new DNAStates();
        //print test
        byte[] t= { (byte)0, (byte)1, (byte)2, (byte)3};
        System.out.println(printBytes(t));
        System.out.println(printByte((byte)0xFF));
        System.out.println( printByte((byte)0x0F));
        
        //if 4 bases, in 1 byte:  00|00|00|00
        //          4 positions:   3  2  1  0

        //00000011 mask for position 0
        //00001100 mask for position 1
        //00110000 mask for position 2
        //11000000 mask for position 3
        byte[] maskArray={(byte)0x03,(byte)0x0C,(byte)0x30,(byte)0xC0};
        // maskArray: [ 00000011 00001100 00110000 11000000 ]
        
        System.out.println("------------------------------");
        String s3="ATCGATCGGCTAGCTAGCATCA";
        String tested=s3;
        System.out.println("QUERY: "+s3);
        System.out.println("QUERY LENGTH: "+s3.length());
        System.out.println("------------------------------");
        System.out.println("-- TEST COMPRESSION");
        
        
        //writing to DNA to bytes
        
        //chose number of bytes necessary to represent this DNA
        //we can code 4 bases in 1 byte
        int byteCount=new Double(Math.ceil(((0.0+tested.length())/4))).intValue();
        System.out.println("byteCount:"+byteCount);
        byte[] kmer=new byte[byteCount];
        
        byte fourBasesByte=0x00;
        int insertions=0;
        for (int i = 0; i < tested.length(); i++) {
            if ( (i>0) & (i%4==0) ) {
                System.out.println("i:"+i);
                System.out.println("i/4:"+(i/4));
                System.out.println("res at "+i+":"+printByte((byte)fourBasesByte));        
                //instance copy this byte in byte array
                kmer[(i/4)-1]=new Byte(fourBasesByte).byteValue();
                //reset pivot byte
                fourBasesByte=0x00;
            }
            
            char c=tested.charAt(i);
            System.out.println("c:"+c);
            byte bc=states.stateToByte(c);
            fourBasesByte= (byte) ( fourBasesByte | (bc << (2*(i%4))) ) ;            
        }
        //last mer fragment was last than 4 bases
        //just add the 
        if (insertions<kmer.length) {
            System.out.println("i:AFTER");
            System.out.println("res at AFTER:"+printByte((byte)fourBasesByte)); 
            kmer[kmer.length-1]=fourBasesByte;
        }
        
        System.out.println("res-total:"+printBytes(kmer));        
        
        System.out.println("---");
        
        //reading bytes to DNA
        byte res=(byte)0b11100100; //=ATCG
        System.out.println("p1:"+printByte((byte) ((res & 0x03) >> 0) ));
        System.out.println("p2:"+printByte((byte) ((res & 0x0C) >> 2) ));
        System.out.println("p3:"+printByte((byte) ((res & 0x30) >> 4) ));
        System.out.println("p4:"+printByte((byte) ((res & 0xC0) >> 6) ));
        
        System.out.println("-- TEST EXPANSION");
        int k=s3.length();
        
        StringBuilder sb=new StringBuilder(k);
        
        for (int i = 0; i < kmer.length; i++) {
            byte b = kmer[i];
            System.out.println("expanding next byte:"+printByte((byte)b)); 
            for (int j = 0; j < 4; j++) {
                //to skip last unused bits
                if ( (i*4)+j > k-1) {
                    break;
                }
                int residue=b;
                residue= (byte) (b & maskArray[j]) ;
                //masking by 0x00FF necessary to avoid that the head bits of integer
                //shifted too (when shift > 5), as in java they are 1 (positive/negative interval of int) by default
                residue=residue & 0x00FF; 
                residue= (byte) (residue >> (j*2));
                System.out.println("p      "+((i*4)+j)+": "+printByte( (byte)residue ));
                sb.append( states.byteToState((byte)residue)   );
            }
        }
        System.out.println("------------------------------");
        System.out.println("RESULT: "+sb.toString());
        System.out.println("SIZE  : "+sb.toString().length());
        System.out.println("QUERY EQUALS RESULT ? : "+tested.equals(sb.toString()));
        
    }
    
    public static String printByte(byte b) {
        StringBuilder sb = new StringBuilder(12);
        sb.append("[ ");
        sb.append( String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0'));
        sb.append(" ");
        sb.append("]");
        return sb.toString();
    }
    
    public static String printBytes(byte[] bytes) {
        StringBuilder sb = new StringBuilder(12);
        sb.append("[ ");
        for (byte b : bytes) {
            sb.append( String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0'));
            sb.append(" ");
        }
        sb.append("]");
        return sb.toString();
    }
    
}
