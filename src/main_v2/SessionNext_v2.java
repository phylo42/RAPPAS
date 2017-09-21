/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import alignement.Alignment;
import core.PProbasSorted;
import core.States;
import core.hash.CustomHash_v2;
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
public class SessionNext_v2 {
    
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
    public CustomHash_v2 hash=null;
    public Float calibrationNormScore=null;
    
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
    public SessionNext_v2(int k, int mink, float alpha, int branchPerEdge, float stateThreshold, float PPStarThreshold, float PPStarThresholdAsLog10) {
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
    
    public void associateHash(CustomHash_v2 hash) {
        this.hash=hash;
    }
    
    public void associateCalibrationScore(float score) {
        this.calibrationNormScore=score;
    }

    public boolean storeHash(File f) {
        try {
            long startTime = System.currentTimeMillis();
            
            Infos.println("Storing of hash");
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,bufferSize));
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
    
    
    
    public static SessionNext_v2 load(File f,boolean loadHash) {
        try {
            long startTime = System.currentTimeMillis();
            FileInputStream fis = new FileInputStream(f);
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,bufferSize));
            int k=ois.readInt();
            int minK=ois.readInt();
            float alpha=ois.readFloat();
            int branchPerEdge=ois.readInt();
            float stateThreshold=ois.readFloat();
            float PPStarThreshold=ois.readFloat();
            float PPStarThresholdAsLog10=ois.readFloat();
            SessionNext_v2 s=new SessionNext_v2(k, minK, alpha, branchPerEdge, stateThreshold, PPStarThreshold, PPStarThresholdAsLog10);
            Infos.println("Loading States");
            s.states = (States)ois.readObject();
            Infos.println("Loading Alignment");
            s.align = (Alignment)ois.readObject();
            Infos.println("Loading Original Tree");
            s.originalTree = (PhyloTree)ois.readObject();
            Infos.println("Loading Extended Tree");
            s.extendedTree = (ExtendedTree)ois.readObject();
            Infos.println("Loading AR Tree");
            s.ARTree = (PhyloTree)ois.readObject();
            Infos.println("Loading of AR node mappings");
            s.nodeMapping = (HashMap<Integer,Integer>)ois.readObject();
//            Infos.println("Loading of PPStats");
//            s.parsedProbas = (PProbasSorted)ois.readObject();
            Infos.println("Loading of calibration");
            s.calibrationNormScore=ois.readFloat();
            if (loadHash) {
                Infos.println("Loading Hash");
                s.hash = (CustomHash_v2)ois.readObject();
            }

            ois.close();
            fis.close();
            long endTime = System.currentTimeMillis();
            Infos.println("Complete session loading " + (endTime - startTime) + " ms");
            return s;
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
            return null;
        }
    }
    
    
}
