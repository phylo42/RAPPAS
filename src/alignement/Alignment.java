/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignement;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import inputs.Fasta;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author ben
 */
public class Alignment implements Serializable {
    
    private static final long serialVersionUID = 1000L;

    //base def of the alignment
    private char[][] charMatrix=null; //char[row][column]
    private int[] colPartitionIds=null;
    private String[] rowLabels=null;

    //metadata
    private ArrayList<Partition> partitions=null;
    private boolean reduced=false;
    private int reducedColumnCount=0;
    private double[] gapProportions=null;
    private double reductionThreshold=0.995;
    

    /**
     * simplest copy constructor
     * @param charMatrix
     * @param colLabels
     * @param rowLabels 
     */
    private Alignment(char[][] charMatrix, String[] rowLabels, int [] colPartitionIds, boolean reduced,double[] gapProportions, double reductionThreshold) {
        this.charMatrix=charMatrix;
        this.rowLabels=rowLabels;
        this.colPartitionIds=colPartitionIds;
        this.reduced=reduced;
        this.reductionThreshold=reductionThreshold;
    }
    
    /**
     * constructor copy 
     * @return 
     */
    public Alignment copy(){
        return new Alignment(charMatrix, rowLabels,colPartitionIds,reduced,gapProportions,reductionThreshold);
    }

    /**
     * build a new alignment from a list of fasta objects,\n
     * order is kept 
     * @param fastas 
     */
    public Alignment(List<Fasta> fastas) {
        fillAlignment(fastas);
    }

    /**
     * build a new alignment from a list of fasta objects,\n
     * order is kept, partitions define the column labels
     * @param fastas 
     */
    public Alignment(List<Fasta> fastas, ArrayList<Partition> partitions) {
        fillAlignment(fastas);
        this.partitions=partitions;
        //fill column names with partitions
        for (int i = 0; i < partitions.size(); i++) {
            Partition p=partitions.get(i);
            for (int j = p.getStart(); j < p.getEnd()+1; j++) {
                colPartitionIds[j]=i;
            }         
        }
    }

    private void fillAlignment(List<Fasta> fastas) {
        //to calculate gaps proportions in all columns
        gapProportions=new double[fastas.get(0).getSequence(false).length()];
        
        for (int i = 0; i < fastas.size(); i++) {
            Fasta f = fastas.get(i);
            //init matrix 
            if (i==0) {
                charMatrix=new char[fastas.size()][f.getSequence(false).length()];
                colPartitionIds=new int[f.getSequence(false).length()];
                rowLabels=new String[fastas.size()];
            }
            //fill matrix
            for (int j = 0; j < f.getSequence(false).length(); j++) {
                char c=f.getSequence(false).charAt(j);
                charMatrix[i][j]=c;
                if (c=='-' || c=='.') {
                    gapProportions[j]++;
                }
            }
            rowLabels[i]=f.getHeader();
        }
        //gap proportions 
        for (int j = 0; j < fastas.get(0).getSequence(false).length(); j++) {
            gapProportions[j]/=0.0+fastas.size();
        }
    }
    
    /**
     * reduces alignment by deleting all columns containg a proportion of gaps 
     * (dash or dot) >= to the given ratio (note that this operation copy the whole array). 
     * @param ratio 
     */
    public void reduceAlignment(double ratio) {
        this.reductionThreshold=ratio;
        boolean[] toRemove=new boolean[gapProportions.length];
        for (int j = 0; j < gapProportions.length; j++) {
            if (gapProportions[j]>=reductionThreshold) {
                toRemove[j]=true;
                reducedColumnCount++;
            }
        }
        
        char[][] newCharMatrix=new char[rowLabels.length][charMatrix[0].length-reducedColumnCount];
        //fill new matrix
        for (int i = 0; i < charMatrix.length; i++) {
            int shift=0;
            for (int j = 0; j < charMatrix[0].length; j++) {
                if (toRemove[j]) {
                    shift++;
                    continue;
                }
                newCharMatrix[i][j-shift]=charMatrix[i][j];
            }
        }
        char[][] old=charMatrix;
        charMatrix=newCharMatrix;
        old=null;
        System.gc();//to free memory in case these matrices are large
        
        reduced=true;
    }
    
    public double[] getGapProportions() {
        return this.gapProportions;
    }
    
    /**
     * add a single sequence to an existing alignment
     * @param label
     * @param seq 
     */
    public void addSequence(String label, char[] seq) {
        //reinstantiate table with a new line
        char[][] newCharMatrix=new char[charMatrix.length+1][charMatrix[0].length];
        for(int i=0; i<charMatrix.length; i++)
            System.arraycopy(charMatrix[i], 0, newCharMatrix[i], 0, charMatrix[i].length);
        System.arraycopy(seq, 0, newCharMatrix[charMatrix.length], 0, charMatrix[0].length);
        charMatrix=newCharMatrix;
        //reinstantiate labels with a new element
        String[] newRowLabels=new String[rowLabels.length+1];
        System.arraycopy(rowLabels, 0, newRowLabels, 0, rowLabels.length);
        newRowLabels[rowLabels.length]=label;
        rowLabels=newRowLabels;
        
        //update gap proportions: divide by previous #seqs and multiply by new #seqs
        for (int j = 0; j < charMatrix[0].length; j++) {
            char c=charMatrix[charMatrix.length-1][j];
            int increment=0;
            if (c=='-' || c=='.') {increment=1;}
            gapProportions[j]= ((gapProportions[j]*(charMatrix.length-1))+increment)/charMatrix.length;
        }

        
        
    }
    
    /**
     * add several sequences to an existing alignment
     * @param labels
     * @param seqs 
     */
    public void addAllSequences(String[] labels, ArrayList<char[]> seqs) {
        //reinstantiate table with a new line
        char[][] newCharMatrix=new char[charMatrix.length+labels.length][charMatrix[0].length];
        for(int i=0; i<charMatrix.length; i++)
            newCharMatrix[i]=charMatrix[i];
        for (int i = 0; i < seqs.size(); i++) {
            newCharMatrix[charMatrix.length+i]=seqs.get(i);
        }
        
        charMatrix=newCharMatrix;
        //register new gaps
        int[] gapsCount = new int[charMatrix[0].length]; //#gaps by site
        for(int i=0; i<seqs.size(); i++) {
            for (int j = 0; j < seqs.get(0).length; j++) {
                char c=seqs.get(i)[j];
                if (c=='-' || c=='.') {gapsCount[j]+=1;}
            }
        }
        //reinstantiate table with a new line
        String[] newRowLabels=new String[rowLabels.length+labels.length];
        for(int i=0; i<rowLabels.length; i++)
            newRowLabels[i]=rowLabels[i];
        for(int i=0; i<labels.length; i++)
            newRowLabels[rowLabels.length+i]=labels[i];
        rowLabels=newRowLabels;
        
        //update gap proportions: divide by previous #seqs and multiply by new #seqs
        for (int j = 0; j < charMatrix[0].length; j++) {
            gapProportions[j]= ((gapProportions[j]*(charMatrix.length-seqs.size()))+gapsCount[j])/charMatrix.length;
        }
    }
    
    
    
    /**
     * Removes a specific sequence: use with parsimony, this reinstanciates all arrays.
     * @param label 
     */
    public void removeSequence(String label) {
        boolean found=false;
        //copy table, remove only the concerned sequence
        char[][] newMatrix=new char[charMatrix.length-1][charMatrix[0].length];
        String[] newLabels=new String[rowLabels.length-1];
        int shift=0;//line index shift
                
        for (int i = 0; i < rowLabels.length; i++) {
            String rowLabel = rowLabels[i];
            if (rowLabel.equals(label)) {
                found=true;
                shift++;
                //correct gap proportion before removal of this sequence
                for (int j = 0; j < charMatrix[0].length; j++) {
                    char c=charMatrix[i][j];
                    int decrement=0;
                    if (c=='-' || c=='.') {
                        decrement=1;
                    }
                    gapProportions[j]= ((gapProportions[j]*charMatrix.length)-decrement)/newMatrix.length;
                }
                continue;
            }
            //copy data in table where sequence removed
            newLabels[i-shift]=rowLabels[i];
            for (int k = 0; k < charMatrix[i].length; k++) {
                newMatrix[i-shift][k]=charMatrix[i][k];
            }
        }
        //update this object
        rowLabels=null;
        charMatrix=null;
        rowLabels=newLabels;
        charMatrix=newMatrix;

        assert found==true;
    }
    
    /**
     * return a particular sequence as a fasta (keeping gaps)
     * @param label
     * @return 
     */
    public Fasta getFasta(String label, boolean withGaps) {
        Fasta f=null;
        for (int i = 0; i < rowLabels.length; i++) {
            if (rowLabels[i].equals(label)) {
                String seq=new String(charMatrix[i]);
                if (withGaps) {seq.replaceAll("-", "");}
                f=new Fasta(label, seq);
                break;
            }
            
        }
        return f;
    }
    
    /**
     * get this alignment as a list of Fasta objects
     * @param withGaps
     * @return 
     */
    public List<Fasta> getAllFasta(boolean withGaps) {
        List<Fasta> l=new ArrayList<>(rowLabels.length);
        for (int i = 0; i < rowLabels.length; i++) {
            String seq=new String(charMatrix[i]);
            if (withGaps) {seq.replaceAll("-", "");}
            l.add(new Fasta(rowLabels[i], seq));  
        }
        return l;
    }
    
    
    public char[][] getCharMatrix() {
        return charMatrix;
    }

    public String[] getRowLabels() {
        return rowLabels;
    }

    public void setRowLabels(String[] rowLabels) {
        this.rowLabels = rowLabels;
    }
    
    public int getLength() {
        return charMatrix[0].length;
    }
    

    /**
     * output alignment as a Fasta file
     * @param f
     * @throws IOException 
     */
    public void writeAlignmentAsFasta(File f) throws IOException {
        BufferedWriter br=new BufferedWriter(new FileWriter(f),4096);
        for (int i = 0; i < charMatrix.length; i++) {
            br.append(">");
            br.append(rowLabels[i]);
            br.append("\n");
            br.append(String.copyValueOf(charMatrix[i]));
            br.append("\n");
        }
        br.close();
    }
    
    /**
     * output alignment as a Phylip file
     * @param f
     * @throws IOException 
     */
    public void writeAlignmentAsPhylip(File f) throws IOException {
        
        
        //here choose 50 because PAML don't allow more in original sources !
        int allowedLabelSize=250; //in fact(48 +2 spaces)
        //this can be changed in paml source, rising #define LSPNAME value
        //changed in baseml.c then recompiled
        
        //paml also wants 2 spaces between label and sequence
        int numColumns=2;
        
        BufferedWriter br=new BufferedWriter(new FileWriter(f),4096);
        br.append(String.valueOf(charMatrix.length)+" "+String.valueOf(charMatrix[0].length)+"\n");     
        for (int i = 0; i < charMatrix.length; i++) {
            //add label
            String label=rowLabels[i];
            if (label.length()>allowedLabelSize) {
                label=label.substring(0, allowedLabelSize-2);
            }
            br.append(label);

            for (int j = 0; j < (allowedLabelSize-label.length()); j++)  {
                br.write(' ');
            }
            //add sequence
            for (int j=0;j<charMatrix[0].length;j++) {
                if (j!=0 && (j%allowedLabelSize)==0)
                    if (allowedLabelSize*numColumns==0)
                        br.write('\n');
                    else
                        br.write(' ');
                br.write(charMatrix[i][j]);
            }
            br.append("\n");
        }
        br.close();
    }
    
    /**
     * simple textual description of this alignment
     * @param extended
     * @return 
     */
    public String describeAlignment(boolean extended) {
        StringBuilder sb=new StringBuilder();
        sb.append("Dimension: "+charMatrix[0].length+"x"+charMatrix.length+" (colxline)");
        if (extended) {
            if (partitions!=null) {
                sb.append("\nPartitions: "+partitions.size()+"\n");
                for (int i = 0; i < partitions.size(); i++) {
                    sb.append("P"+i+" :"+partitions.get(i)+"\n");
                }
                sb.append("P assigned: ");
                if (colPartitionIds!=null) {
                    sb.append("YES\n");
                } else {
                    sb.append("NO\n");
                } 
            }else {
                sb.append(" Partitions: none");
            }
            if (reduced) {
                sb.append(" reduced:"+reducedColumnCount);
            } else {
                sb.append(" reduced:false");
            }
        } 
        return sb.toString();
    }
    
    protected void printAlignment() {
        for (int i = 0; i < charMatrix.length; i++) {
            for (int j = 0; j < charMatrix[0].length; j++) {
                System.out.print(charMatrix[i][j]);
            }
            System.out.println("");
        }
    }
    
    
    public static void main(String[] args) {
        
        Fasta f1=new Fasta("1", "ATC-TG--GT---");
        Fasta f2=new Fasta("2", "A-C-T---GT---");
        Fasta f3=new Fasta("3", "AT--T-C-GT---");
        
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
        
        a.addSequence("4", new String("A---TCC-GT--A").toCharArray());
        System.out.println("ADD 1");
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        
        System.out.println("ADD 2");
        String[] labels={"5","6"};
        ArrayList<char[]> seqs= new ArrayList<>();
        seqs.add(new String("A---TCC-GT--A").toCharArray());
        seqs.add(new String("AT---CA-GT-AA").toCharArray());
        a.addAllSequences(labels, seqs);
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        
        System.out.println("REMOVE 1");
        a.removeSequence("3");
        a.printAlignment();
        System.out.println(Arrays.toString(a.getGapProportions()));
        System.out.println(a.describeAlignment(true));
        
        System.out.println("REDUCTION");
        a.reduceAlignment(0.995);
        System.out.println(a.describeAlignment(true));
        a.printAlignment();
       
        
        System.out.println("END");
        
    }
    
    
    
}
