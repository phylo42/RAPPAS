/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import core.PProbasSorted;
import core.hash.CustomHash;
import core.States;
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
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class SessionNext {
    
    int k=-1;
    int minK=-1;
    float stateThreshold=Float.MIN_VALUE;
    float wordThreshold=Float.MIN_VALUE;
    float factor=1.0f;
    
    
    States states=null;
    Alignment align=null;
    PhyloTree tree=null;
    PProbasSorted parsedProbas=null;    
    CustomHash hash=null;
    
    
    public SessionNext(int k, int mink, float factor, float stateThreshold, float wordThreshold) {
        this.k=k;
        this.minK=mink;
        this.factor=factor;
        this.stateThreshold=stateThreshold;
        this.wordThreshold=wordThreshold;
    }
    
    public void associateStates(States s) {
        this.states=s;
    }
    
    public void associateInputs(ARResults im) {
        this.tree=im.getARTree();
        this.align=im.getExtendedAlignment();
        this.parsedProbas=im.getPProbas();
    }
    
    public void associateHash(CustomHash hash) {
        this.hash=hash;
    }
    
    
    public boolean store(File f) {
        try {
            long startTime = System.currentTimeMillis();
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
            Infos.println("Storing of Tree");
            oos.writeObject(tree);
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
    
    public static SessionNext load(File f,boolean loadHash) {
        try {
            long startTime = System.currentTimeMillis();
            FileInputStream fis = new FileInputStream(f);
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,4096));
            int k=ois.readInt();
            int minK=ois.readInt();
            float factor=ois.readFloat();
            float stateThreshold=ois.readFloat();
            float wordThreshold=ois.readFloat();
            SessionNext s=new SessionNext(k, minK, factor, stateThreshold, wordThreshold);
            Infos.println("Loading States");
            s.states = (States)ois.readObject();
            Infos.println("Loading Alignment");
            s.align = (Alignment)ois.readObject();
            Infos.println("Loading Tree");
            s.tree = (PhyloTree)ois.readObject();
            Infos.println("Loading PPStats");
            s.parsedProbas = (PProbasSorted)ois.readObject();
            if (loadHash) {
                Infos.println("Loading Hash");
                s.hash = (CustomHash)ois.readObject();
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
