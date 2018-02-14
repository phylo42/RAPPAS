/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignment;

import alignement.Alignment;
import etc.Infos;
import inputs.FASTAPointer;
import inputs.Fasta;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ben
 */
public class TEST_Alignment {

    
    public static void main(String[] args) {
        
        Fasta f1=new Fasta("1", "-ATC-TG--GT---");
        Fasta f2=new Fasta("2", "-A-C-T---GT---");
        Fasta f3=new Fasta("3", "-AT--T-C-GT---");
        
        List<Fasta> align=new ArrayList<>(); 
        align.add(f1);
        align.add(f2);
        align.add(f3);
        
        System.out.println("HERE !");
        
        Alignment a=new Alignment(align);
        
        System.out.println(a.describeAlignment(true));
        System.out.println("ORIGINAL:");
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        System.out.println(Arrays.toString(a.getGapIntervals()));
        
        a.addSequence("4", new String("-A---TCC-GT--A").toCharArray());
        System.out.println("ADD 1");
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        System.out.println(Arrays.toString(a.getGapIntervals()));
        
        System.out.println("ADD 2");
        String[] labels={"5","6"};
        ArrayList<char[]> seqs= new ArrayList<>();
        seqs.add(new String("-A---TCC-GT--A").toCharArray());
        seqs.add(new String("AAT---CA-GT-AA").toCharArray());
        a.addAllSequences(labels, seqs);
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        System.out.println(Arrays.toString(a.getGapIntervals()));

        System.out.println("REMOVE 1");
        a.removeSequence("3");
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        System.out.println(a.describeAlignment(true));
        System.out.println(Arrays.toString(a.getGapIntervals()));

        System.out.println("REDUCTION");
        a.reduceAlignment(0.80);
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        System.out.println(a.describeAlignment(true));
        System.out.println(Arrays.toString(a.getGapIntervals()));
        
        
        System.out.println("END");
        
        Infos.println("Loading Alignment...");
        FASTAPointer fp=new FASTAPointer(new File("/home/ben/Dropbox/viromeplacer/test_datasets/accuracy_tests/R_analysis/BOLD_matk_v0.95/A14_older.mfa"), false);
        Fasta fasta=null;
        ArrayList<Fasta> fastas=new ArrayList<>();
        while ((fasta=fp.nextSequenceAsFastaObject())!=null) {
            fastas.add(fasta);
        }
        a=new Alignment(fastas);
        System.out.println(a.describeAlignment(true));
        double[] gapProportions1 = a.getGapProportions();
//        for (int i = 0; i < gapProportions1.length; i++) {
//            double d = gapProportions1[i];
//            System.out.println("i:"+i+ " prop="+d);
//        }
        a.reduceAlignment(1.0);
        gapProportions1 = a.getGapProportions();
//        for (int i = 0; i < gapProportions1.length; i++) {
//            double d = gapProportions1[i];
//            System.out.println("i:"+i+ " prop="+d);
//        }
        System.out.println(a.describeAlignment(true));
        fp.closePointer();
        
        try {
            a.writeAlignmentAsFasta(new File("/home/ben/test.fasta"));
        } catch (IOException ex) {
            Logger.getLogger(Alignment.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("------------");
        Alignment a2=a.copy();
        gapProportions1 = a2.getGapProportions();
//        for (int i = 0; i < gapProportions1.length; i++) {
//            double d = gapProportions1[i];
//            System.out.println("i:"+i+ " prop="+d);
//        }
        a2.removeSequence("SDH3502-15.matK_Rosales_1-790/1-790");
        System.out.println(a2.describeAlignment(true));
        gapProportions1 = a2.getGapProportions();
//        for (int i = 0; i < gapProportions1.length; i++) {
//            double d = gapProportions1[i];
//            System.out.println("i:"+i+ " prop="+d);
//        }
        a2.removeSequence("DBMPP221-14.matK_Rosales_1-784/1-784");
        System.out.println(a2.describeAlignment(true));
        gapProportions1 = a2.getGapProportions();
//        for (int i = 0; i < gapProportions1.length; i++) {
//            double d = gapProportions1[i];
//            System.out.println("i:"+i+ " prop="+d);
//        }
        System.out.println("Reduce at 1.0");
        a2.reduceAlignment(1.0);
        System.out.println(a2.describeAlignment(true));
        
        try {
            a2.writeAlignmentAsFasta(new File("/home/ben/test2.fasta"));
        } catch (IOException ex) {
            Logger.getLogger(Alignment.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
}
