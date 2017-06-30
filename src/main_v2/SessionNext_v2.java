/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

import alignement.Alignment;
import core.PProbasSorted;
import core.States;
import core.hash.SimpleHash_v2;
import etc.Infos;
import inputs.ARProcessResults;
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
    
    int k=-1;
    int minK=-1;
    float stateThreshold=Float.MIN_VALUE;
    float wordThreshold=Float.MIN_VALUE;
    float factor=1.0f;
    
    
    States states=null;
    Alignment align=null;
    PhyloTree originalTree=null;
    ExtendedTree extendedTree=null;
    PhyloTree ARTree=null;
    /**
     * map(ARTree NodeID)= extended tree NodeID
     */
    HashMap<Integer,Integer> nodeMapping=null;
    PProbasSorted parsedProbas=null;    
    SimpleHash_v2 hash=null;
    
    
    public SessionNext_v2(int k, int mink, float factor, float stateThreshold, float wordThreshold) {
        this.k=k;
        this.minK=mink;
        this.factor=factor;
        this.stateThreshold=stateThreshold;
        this.wordThreshold=wordThreshold;
    }
    
    public void associateStates(States s) {
        this.states=s;
    }
    
    public void associateInputs(ARProcessResults arpl) {
        this.originalTree=arpl.getOriginalTree();
        this.extendedTree=arpl.getExtendedTree();
        this.ARTree=arpl.getARTree();
        this.nodeMapping=arpl.getTreeMapping();
        this.align=arpl.getExtendedAlignment();
        this.parsedProbas=arpl.getPProbas();
    }
    
    public void associateHash(SimpleHash_v2 hash) {
        this.hash=hash;
    }
    
    
    public boolean storeFullHash(File f) {
        try {
            long startTime = System.currentTimeMillis();
            
            Infos.println("Storing of \"FULL\" hash");
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
            oos.writeInt(k);
            oos.writeInt(minK);
            oos.writeFloat(factor);
            oos.writeFloat(stateThreshold);
            oos.writeFloat(wordThreshold);
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
            Infos.println("Storing of PPStats");
            oos.writeObject(parsedProbas);
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
    
    /**
     * should be called AFTER storeFullHash
     * @param f
     * @return 
     */
    public boolean storeMediumHash(File f) {
        try {
            long startTime = System.currentTimeMillis();
            
            Infos.println("Storing of \"MEDIUM\" hash");
            hash.reduceToMediumHash();
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
            oos.writeInt(k);
            oos.writeInt(minK);
            oos.writeFloat(factor);
            oos.writeFloat(stateThreshold);
            oos.writeFloat(wordThreshold);
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
            Infos.println("Storing of PPStats");
            oos.writeObject(parsedProbas);
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
    
    /**
     * should be called AFTER storeSmallHash
     * @param f
     * @return 
     */
    public boolean storeSmallHash(File f, int X) {
        try {
            long startTime = System.currentTimeMillis();
            
            Infos.println("Storing of \"SMALL\" hash");
            hash.reducetoSmallHash(X);
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
            oos.writeInt(k);
            oos.writeInt(minK);
            oos.writeFloat(factor);
            oos.writeFloat(stateThreshold);
            oos.writeFloat(wordThreshold);
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
            Infos.println("Storing of PPStats");
            oos.writeObject(parsedProbas);
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
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,4096));
            int k=ois.readInt();
            int minK=ois.readInt();
            float factor=ois.readFloat();
            float stateThreshold=ois.readFloat();
            float wordThreshold=ois.readFloat();
            SessionNext_v2 s=new SessionNext_v2(k, minK, factor, stateThreshold, wordThreshold);
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
            Infos.println("Loading of PPStats");
            s.parsedProbas = (PProbasSorted)ois.readObject();
            if (loadHash) {
                Infos.println("Loading Hash");
                s.hash = (SimpleHash_v2)ois.readObject();
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
