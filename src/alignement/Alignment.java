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
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.LinkedList;
import core.States;
import etc.Infos;
import etc.exceptions.NonSupportedStateException;
import java.util.HashMap;

/**
 *
 * @author ben
 */
public class Alignment implements Serializable {
    
    private static final long serialVersionUID = 1000L;

    //states associated to this object
    States s=null;
    
    //base def of the alignment
    private char[][] charMatrix=null; //char[row][column]
    private int[] colPartitionIds=null;
    private String[] rowLabels=null;

    //metadata
    private ArrayList<Partition> partitions=null;
    private boolean reduced=false;
    private int reducedColumnCount=0;
    //% of gap in each site
    private double[] gapProportions=null;
    private double reductionThreshold=0.99;
    //gap intervals
    private ArrayList<Integer>[] gapIntervals=null;
    

    /**
     * simplest copy constructor
     * @param charMatrix
     * @param colLabels
     * @param rowLabels 
     */
    private Alignment(States s, char[][] charMatrix, String[] rowLabels, int [] colPartitionIds, boolean reduced,int reducedColumnCount, double[] gapProportions, double reductionThreshold,ArrayList<Integer>[] gapIntervals) {
        this.s=s;
        this.charMatrix=charMatrix;
        this.rowLabels=rowLabels;
        this.colPartitionIds=colPartitionIds;
        this.reduced=reduced;
        this.reducedColumnCount=reducedColumnCount;
        this.reductionThreshold=reductionThreshold;
        this.gapProportions=gapProportions;
        this.gapIntervals=gapIntervals;
    }
    
    /**
     * constructor copy 
     * @return 
     */
    public Alignment copy(){
        return new Alignment(s,charMatrix, rowLabels,colPartitionIds,reduced,reducedColumnCount,gapProportions,reductionThreshold,gapIntervals);
    }

    /**
     * build a new alignment from a list of fasta objects,\n
     * order is kept 
     * @param fastas 
     */
    public Alignment(States s, List<Fasta> fastas) {
        this.s=s;
        fillAlignment(fastas);
    }

    /**
     * build a new alignment from a list of fasta objects,\n
     * order is kept, partitions define the column labels
     * @param fastas 
     */
    @Deprecated
    public Alignment(States s, List<Fasta> fastas, ArrayList<Partition> partitions) {
        this.s=s;
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

    /**
     * fills alignment at instanciation
     * @param fastas 
     */
    private void fillAlignment(List<Fasta> fastas) {
        //to calculate gaps proportions in all columns
        gapProportions=new double[fastas.get(0).getSequence(false).length()];
        
        //list of gap intervals
        gapIntervals=new ArrayList[fastas.get(0).getSequence(false).length()];
        
        //register how many ambiguous characters were found
        HashMap<Character,Integer> ambiguityInventory=new HashMap<>();        
        for (int i = 0; i < fastas.size(); i++) {
            Fasta f = fastas.get(i);
            //init matrix 
            if (i==0) {
                charMatrix=new char[fastas.size()][f.getSequence(false).length()];
                colPartitionIds=new int[f.getSequence(false).length()];
                rowLabels=new String[fastas.size()];
            }
            //test seq length, should be the same
            checkLength(f.getHeader(),f.getSequence(false).toCharArray());
            //prepare interval def
            int firstGapIndex=-1;
            char previousChar='n';
            //fill matrix
            for (int j = 0; j < f.getSequence(false).length(); j++) {
                char c=f.getSequence(false).charAt(j);
                //test char per char at read
                if (c!='-') {
                    if (!ambiguityInventory.containsKey(c)) {
                        ambiguityInventory.put(c, 1);
                    } else {
                        ambiguityInventory.put(c,ambiguityInventory.get(c)+1);
                    }
                }
                //either this is a known ambiguity or it raises an exception
                try {
                    //either this is a known ambiguity definition
                    if (!s.isAmbiguous(c)) {
                    //or this is a valid state or it raises an exception
                        s.stateToByte(c);                        
                    }
                } catch (NonSupportedStateException ex) {
                    ex.printStackTrace(System.err);
                    System.out.println("Reference alignment contains a non supported ambiguous state.");
                    System.exit(1); //do not exit here, AR will take care of transforming them to gaps
                }
                charMatrix[i][j]=c;
                if (c=='-') {
                    gapProportions[j]++;
                    if (previousChar!='-') {
                        //activate gap counter
                        if (firstGapIndex==-1) { firstGapIndex=j;}
                    }
                } else {
                    if (firstGapIndex!=-1) {
                        if (gapIntervals[firstGapIndex]==null)
                            gapIntervals[firstGapIndex]=new ArrayList<>(5);
                        int length=j-firstGapIndex;
                        if (!gapIntervals[firstGapIndex].contains(length)) {
                            gapIntervals[firstGapIndex].add(length);
                        }
                        firstGapIndex=-1;
                    }
                }
                previousChar=c;
            }
            //close last gap interval
            if (firstGapIndex!=-1) {
                gapIntervals[firstGapIndex]=new ArrayList<>(5);
                int length=f.getSequence(false).length()-firstGapIndex;
                if (!gapIntervals[firstGapIndex].contains(length)) {
                    gapIntervals[firstGapIndex].add(length);
                }
                firstGapIndex=-1;
            }
            
            //sequence labels
            rowLabels[i]=f.getHeader();
        }
        //warn about ambiguities
        if (ambiguityInventory.keySet().size()>0) {
            System.out.println("Some ambiguous states were found in the alignment.");
            for (Iterator<Character> iterator = ambiguityInventory.keySet().iterator(); iterator.hasNext();) {
                Character key = iterator.next();
                Infos.println("Ambiguous state: char='"+key+"' occurences="+ambiguityInventory.get(key));
            }
        }

        
        
        //gap proportions 
        for (int j = 0; j < fastas.get(0).getSequence(false).length(); j++) {
            gapProportions[j]/=0.0+fastas.size();
        }
        
        
        
    }
    
    /**
     * test if new seq is same length as alignment seqs
     * @param label
     * @param seq 
     */
    private void checkLength(String label, char[] seq){
        //if has already at least 1 seq
        if (rowLabels[0]!=null) {
            // Check if sequences contains same number of sites
            if(seq.length!=charMatrix[0].length){
                System.err.println("Error: Sequences in the input alignment don't have same number of sites!");
                System.err.println("This sequence has a length different from 1st sequence: "+label+" (1st is "+rowLabels[0]+")");
                System.exit(1);
            }
        }
    }
    
    
    /**
     * fill the gapInterval table
     */
    private void updateGapIntervals() {
        gapIntervals=new ArrayList[charMatrix[0].length];
        for (int i = 0; i < charMatrix.length; i++) {
            //prepare interval def
            int firstGapIndex=-1;
            char previousChar='n';
            //fill matrix
            for (int j = 0; j < charMatrix[i].length; j++) {
                char c=charMatrix[i][j];
                if (c=='-') {
                    if (previousChar!='-') {
                        //activate gap counter
                        if (firstGapIndex==-1) { firstGapIndex=j;}
                    }
                } else {
                    if (firstGapIndex!=-1) {
                        if (gapIntervals[firstGapIndex]==null)
                            gapIntervals[firstGapIndex]=new ArrayList<>(5);
                        int length=j-firstGapIndex;
                        if (!gapIntervals[firstGapIndex].contains(length)) {
                            gapIntervals[firstGapIndex].add(length);
                        }
                        firstGapIndex=-1;
                    }
                }
                previousChar=c;
            }
        }
        
    }
    
    
    
    /**
     * reduces alignment by deleting all columns containg a proportion of gaps 
     * (dash or dot) >= to the given ratio (note that this operation copy the whole array). 
     * @param ratio 
     */
    public void reduceAlignment(double ratio) {
        this.reductionThreshold=ratio;
        
        reducedColumnCount=0;
        boolean[] toRemove=new boolean[gapProportions.length];
        for (int j = 0; j < gapProportions.length; j++) {
            if (gapProportions[j]>=reductionThreshold) {
                toRemove[j]=true;
                reducedColumnCount++;
            }
        }

        
        //fill new "reduced" matrix        
        char[][] newCharMatrix=new char[rowLabels.length][charMatrix[0].length-reducedColumnCount];
        //update gap proportions 
        gapProportions=new double[newCharMatrix[0].length];

        for (int i = 0; i < charMatrix.length; i++) {
            int shift=0;
            for (int j = 0; j < charMatrix[0].length; j++) {
                if (toRemove[j]) {
                    shift++;
                    continue;
                }
                newCharMatrix[i][j-shift]=charMatrix[i][j];
                if (newCharMatrix[i][j-shift]=='-') {
                    gapProportions[j-shift]++;
                }
            }
        }
       
        for (int j = 0; j < newCharMatrix[0].length; j++) {
            gapProportions[j]/=newCharMatrix.length;
        }
        
        charMatrix=newCharMatrix;
        
        System.gc();//to free memory in case these matrices are large
       
        
        //update gap intervals
        updateGapIntervals();
        
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
        //test if new seq length is correct
        checkLength(label, seq);
        
        //reinstantiate table with a new line
        char[][] newCharMatrix=new char[charMatrix.length+1][charMatrix[0].length];
        for(int i=0; i<charMatrix.length; i++) {
            newCharMatrix[i]=charMatrix[i];
        }
        newCharMatrix[charMatrix.length]=seq;
        charMatrix=newCharMatrix;
        //reinstantiate labels with a new element
        String[] newRowLabels=new String[rowLabels.length+1];
        System.arraycopy(rowLabels, 0, newRowLabels, 0, rowLabels.length);
        newRowLabels[rowLabels.length]=label;
        rowLabels=newRowLabels;
        
        //reset gap proportions
        gapProportions=new double[charMatrix[0].length];
        
        //update gap proportions: divide by previous #seqs and multiply by new #seqs
        //update gap interval, just add eventual new intervals represented by this new seqeunce
        int firstGapIndex=-1;
        char previousChar='n';
        for (int j = 0; j < charMatrix[0].length; j++) {
            char c=charMatrix[charMatrix.length-1][j];
            //update proportions
            int increment=0;
            if (c=='-') {
                increment=1;
            }
            gapProportions[j]= ((gapProportions[j]*(charMatrix.length-1))+increment)/charMatrix.length;
            
            //update intervals
            if (c=='-') {
                if (previousChar!='-') {
                    //activate gap counter
                    if (firstGapIndex==-1) { firstGapIndex=j;}
                }
            } else {
                if (firstGapIndex!=-1) {
                    if (gapIntervals[firstGapIndex]==null)
                        gapIntervals[firstGapIndex]=new ArrayList<>(5);
                    int length=j-firstGapIndex;
                    if (!gapIntervals[firstGapIndex].contains(length)) {
                        gapIntervals[firstGapIndex].add(length);
                    }
                    firstGapIndex=-1;
                }
            }
            previousChar=c;
                     
            
        }  
    }
    
    /**
     * add several sequences to an existing alignment
     * @param labels
     * @param seqs 
     */
    public void addAllSequences(String[] labels, ArrayList<char[]> seqs) {
        
        //test if new seq length is correct
        for (int i = 0; i < labels.length; i++) {
            checkLength(labels[i], seqs.get(i));
        }
        
        //reinstantiate sequence table with a new line
        char[][] newCharMatrix=new char[charMatrix.length+labels.length][charMatrix[0].length];
        for(int i=0; i<charMatrix.length; i++)
            newCharMatrix[i]=charMatrix[i];
        for (int i = 0; i < seqs.size(); i++) {
            newCharMatrix[charMatrix.length+i]=seqs.get(i);
        }
        charMatrix=newCharMatrix;
        //reset gap proportions
        gapProportions=new double[charMatrix[0].length];
        //register new gaps
        for(int i=0; i<newCharMatrix.length; i++) {
            for (int j = 0; j < newCharMatrix[i].length; j++) {
                if (newCharMatrix[i][j]=='-') {gapProportions[j]+=1;}
            }
        }
        //reinstantiate labels table with a new line
        String[] newRowLabels=new String[rowLabels.length+labels.length];
        for(int i=0; i<rowLabels.length; i++)
            newRowLabels[i]=rowLabels[i];
        for(int i=0; i<labels.length; i++)
            newRowLabels[rowLabels.length+i]=labels[i];
        rowLabels=newRowLabels;
        
        //update gap proportions: divide by previous #seqs and multiply by new #seqs
        for (int j = 0; j < charMatrix[0].length; j++) {
            gapProportions[j]/=charMatrix.length;
        }
        
        //update gap intervals
        updateGapIntervals();
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
        //reset gap proportions
        gapProportions=new double[newMatrix[0].length];
        
        for (int i = 0; i < rowLabels.length; i++) {
            String rowLabel = rowLabels[i];
            if (rowLabel.equals(label)) {
                found=true;
                shift++;
                continue;
            }
            //copy data in table where sequence removed
            newLabels[i-shift]=rowLabels[i];
            for (int j = 0; j < charMatrix[i].length; j++) {
                newMatrix[i-shift][j]=charMatrix[i][j];
                if (charMatrix[i][j]=='-') {
                    gapProportions[j]+=1.0;
                }
            }
        }
        for (int j = 0; j < charMatrix[0].length; j++) {
            gapProportions[j]/=newMatrix.length;
        }
        
        
        
        //update this object
        rowLabels=null;
        charMatrix=null;
        rowLabels=newLabels;
        charMatrix=newMatrix;

        assert found==true;
        
        //update gap intervals
        updateGapIntervals();
        
    }
    
    /**
     * return a particular sequence as a fasta (keeping gaps)
     * @param label
     * @param withGaps
     * @return 
     */
    public Fasta getFasta(String label, boolean withGaps) {
        Fasta f=null;
        for (int i = 0; i < rowLabels.length; i++) {
            if (rowLabels[i].equals(label)) {
                String seq=new String(charMatrix[i]);
                if (!withGaps) {seq=seq.replaceAll("-", "");}
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
    
    /**
     * gap intervals, defined by a fixed array of lists of variable length:
     * 
     *     - - A T C G - T
     *     A - - T G G - C
     * 
     *    [0|1|2| | | | |align_length]
     *     | | |       |
     *     v v v       v
     * [0] 2 1 1       1
     * [1]   2
     * [2]
     * 
     * @return 
     */
    public ArrayList<Integer>[] getGapIntervals() {
        return this.gapIntervals;
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
    
    @Deprecated
    public Object2IntOpenHashMap<char[]> getKmerRedundancy(int k) {
        
        Object2IntOpenHashMap<char[]> map=new Object2IntOpenHashMap <>();
        map.defaultReturnValue(-1);
        LinkedList<Character> list=new LinkedList<>();
        char[] mer=new char[k];
        for (int i = 0; i < charMatrix.length; i++) {
            char[] cs = charMatrix[i];
            for (int j = 0; j < cs.length; j++) {
                //mer = [j,j+k[
                if (j<k) {
                    list.addLast(cs[j]);
                } else {
                    list.removeFirst();
                    list.addLast(cs[j]);
                    int pos=0;
                    for (Iterator<Character> iterator = list.iterator(); iterator.hasNext();) {
                        mer[pos]=list.get(pos);
                        pos++;
                    }
                    int count=map.getInt(mer);
                    if (count==-1) {
                        map.put(mer, 0);
                    } else {
                        map.put(mer, count+1);
                    }
                }  
            }
        }
        return map;
        
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
                sb.append(" reduced: #cols="+reducedColumnCount);
            } else {
                sb.append(" reduced:false");
            }
        } 
        return sb.toString();
    }
    
    public void printAlignment() {
        for (int i = 0; i < charMatrix.length; i++) {
            for (int j = 0; j < charMatrix[0].length; j++) {
                System.out.print(charMatrix[i][j]);
            }
            System.out.println("");
        }
    }
    
    
        
}
