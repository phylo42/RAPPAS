/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import alignement.Alignment;
import core.PProbasSorted;
import core.States;
import core.hash.CustomHash;
import core.hash.CustomHash_Triplet;
import core.hash.CustomHash_v4_FastUtil81;
//import core.hash.CustomHash_Triplet;
import etc.Infos;
import inputs.ARResults;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import tree.ExtendedTree;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class SessionNext {
    
    private final static int bufferSize=2097152; // buffer of 2mo
    
    
    public int k=-1;
    public int minK=-1;
    public int branchPerEdge=-1;
    public float stateThreshold=Float.MIN_VALUE;
    public float PPStarThreshold=Float.MIN_VALUE;
    public float PPStarThresholdAsLog10=Float.NEGATIVE_INFINITY;
    public float alpha=1.0f;
    
    public States states=null;
    public Alignment align=null;
    public PhyloTree originalTree=null;
    public ExtendedTree extendedTree=null;
    public PhyloTree ARTree=null;
    /**
     * map(ARTree NodeID)= extended tree NodeID
     */
    public HashMap<Integer,Integer> nodeMapping=null;
    public PProbasSorted parsedProbas=null;    
    public CustomHash hash=null;
    //public CustomHash_Triplet hash=null;
    public boolean onlyFakes=false;
    public int hashType=CustomHash.NODES_UNION;
    public Float calibrationNormScore=null;
    
    /**
     *
     */
    public SessionNext() {}

    public int getK() {
        return this.k;
    }

    public int getMinK() {
        return this.minK;
    }

    public float getAlpha() {
        return this.alpha;
    }

    public int getBranchPerEdge() {
        return this.branchPerEdge;
    }

    public float getStateThreshold() {
        return this.stateThreshold;
    }

    public float getPPStarThreshold() {
        return this.PPStarThreshold;
    }

    public float getPPStarThresholdAsLog10() {
        return this.PPStarThresholdAsLog10;
    }
    
    
    /**
     *
     * @param k the value of k
     * @param mink the value of mink
     * @param alpha the value of alpha
     * @param branchPerEdge the number of branches injected on each edge
     * @param stateThreshold the value of stateThreshold
     * @param PPStarThreshold the value of PPStarThreshold
     * @param PPStarThresholdAsLog10 the value of PPStarThresholdAsLog10
     */
    public void associateParameters(int k, int mink, float alpha, int branchPerEdge, float stateThreshold, float PPStarThreshold, float PPStarThresholdAsLog10) {
        this.k=k;
        this.minK=mink;
        this.alpha=alpha;
        this.branchPerEdge=branchPerEdge;
        this.stateThreshold=stateThreshold;
        this.PPStarThreshold=PPStarThreshold;
        this.PPStarThresholdAsLog10=PPStarThresholdAsLog10;
    }
    
    public void associateStates(States s) {
        this.states=s;
    }
    
    public void associateInputs(ARResults arpl) {
        this.originalTree=arpl.getOriginalTree();
        this.extendedTree=arpl.getExtendedTree();
        this.ARTree=arpl.getARTree();
        this.nodeMapping=arpl.getTreeMapping();
        this.align=arpl.getExtendedAlignment();
        this.parsedProbas=arpl.getPProbas();
    }
    
    public void associateHash(CustomHash hash, boolean onlyFakes, int hashType) {
        this.hash=hash;
        this.onlyFakes=onlyFakes;
        this.hashType=hashType;
    }
    
    public void associateCalibrationScore(float score) {
        this.calibrationNormScore=score;
    }

    public boolean storeHash(File f) {
        try {
            long startTime = System.currentTimeMillis();
            
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,bufferSize));
            oos.writeInt(hashType);
            oos.writeInt(k);
            oos.writeInt(minK);
            oos.writeFloat(alpha);
            oos.writeInt(branchPerEdge);
            oos.writeFloat(stateThreshold);
            oos.writeFloat(PPStarThreshold);
            oos.writeFloat(PPStarThresholdAsLog10);
            Infos.println("Storing of States");
            oos.writeObject(states);
            Infos.println("Storing of Alignment");
            oos.writeObject(align);
            Infos.println("Storing of Original Tree");
            oos.writeObject(originalTree);
            Infos.println("Storing of Extended Tree");
            oos.writeObject(extendedTree);
            Infos.println("Storing of AR Tree");
            oos.writeObject(ARTree);
            Infos.println("Storing of AR node mappings");
            oos.writeObject(nodeMapping);
//            Infos.println("Storing of PPStats");
//            oos.writeObject(parsedProbas);
            Infos.println("Storing of Calibration");
            oos.writeFloat(calibrationNormScore); 
            Infos.println("Storing of Hash");
            oos.writeBoolean(onlyFakes);
            oos.writeObject(hash);           
            oos.close();
            fos.close();
            
            long endTime = System.currentTimeMillis();
            Infos.println("Complete session storage " + (endTime - startTime) + " ms");
            Infos.println("Session stored in : "+f.getAbsolutePath());
        } catch (IOException ex) {
            ex.printStackTrace();
            return false;
        }
        return true;
    }
    
    /**
     *
     * @param f
     * @param loadHash
     */
    public void load(File f,boolean loadHash) {
        try {
            long startTime = System.currentTimeMillis();
            FileInputStream fis = new FileInputStream(f);
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,bufferSize));
            this.hashType=ois.readInt();
            this.k=ois.readInt();
            this.minK=ois.readInt();
            this.alpha=ois.readFloat();
            this.branchPerEdge=ois.readInt();
            this.stateThreshold=ois.readFloat();
            this.PPStarThreshold=ois.readFloat();
            this.PPStarThresholdAsLog10=ois.readFloat();
            associateParameters(k, minK, alpha, branchPerEdge, stateThreshold, PPStarThreshold, PPStarThresholdAsLog10);
            Infos.println("Loading States");
            states = (States)ois.readObject();
            Infos.println("Loading Alignment");
            align = (Alignment)ois.readObject();
            Infos.println("Loading Original Tree");
            originalTree = (PhyloTree)ois.readObject();
            Infos.println("Loading Extended Tree");
            extendedTree = (ExtendedTree)ois.readObject();
            Infos.println("Loading AR Tree");
            ARTree = (PhyloTree)ois.readObject();
            Infos.println("Loading of AR node mappings");
            nodeMapping = (HashMap<Integer,Integer>)ois.readObject();
//            Infos.println("Loading of PPStats");
//            s.parsedProbas = (PProbasSorted)ois.readObject();
            Infos.println("Loading of calibration");
            calibrationNormScore=ois.readFloat();
            if (loadHash) {
                Infos.println("Loading Hash");
                onlyFakes=ois.readBoolean();
                if (onlyFakes) {
                    //System.out.println("Loaded DB only contain ancestral kmers associated to fake nodes.");
                } else {
                    System.out.println("Loaded DB only also contain ancestral kmers associated to original nodes.");
                }
                hash = (CustomHash)ois.readObject();
                if (hash instanceof CustomHash_v4_FastUtil81) {
                    Infos.println("HashType: NODES_UNION");
                } else if (hash instanceof CustomHash_Triplet) {
                    Infos.println("HashType: NODES_TRIPLET");
                }
            }

            ois.close();
            fis.close();
            long endTime = System.currentTimeMillis();
            Infos.println("Complete session loading " + (endTime - startTime) + " ms");
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    
}
