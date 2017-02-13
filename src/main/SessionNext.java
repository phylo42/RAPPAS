/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import alignement.Alignment;
import core.Colmer;
import core.ColmerSet;
import core.DNAStates;
import core.ProbabilisticWord;
import core.PProbas;
import core.PProbasSorted;
import core.SimpleHash;
import core.States;
import core.Word;
import core.algos.SequenceKnife;
import core.algos.WordGenerator;
import etc.Environement;
import etc.Infos;
import inputs.Fasta;
import inputs.FASTAPointer;
import inputs.InputManager;
import inputs.InputManagerNext;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import tree.PhyloNode;
import tree.PhyloTree;
import tree.RelaxedTree;

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
    RelaxedTree relaxedTree=null; //optinnal, checked through a boolean flag
    
    SimpleHash hash=null;
    
    
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
    
    public void associateInputs(InputManagerNext im) {
        this.tree=im.getTree();
        if (im.getTree() instanceof RelaxedTree)
            this.relaxedTree=(RelaxedTree)im.getTree();
        else
            this.relaxedTree=null;
        this.align=im.getAlignment();
        this.parsedProbas=im.getPProbas();
    }
    
    
    public boolean store(File f) {
        try {
            long startTime = System.currentTimeMillis();
            FileOutputStream fos = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(fos,4096));
            oos.writeInt(k);
            oos.writeInt(minK);
            oos.writeDouble(factor);
            oos.writeDouble(stateThreshold);
            oos.writeDouble(wordThreshold);
            //flag for relaxedTree
            if (this.relaxedTree!=null)
                oos.writeBoolean(true);
            else
                oos.writeBoolean(false);
            Infos.println("Storing of States");
            oos.writeObject(states);
            Infos.println("Storing of Alignment");
            oos.writeObject(align);
            Infos.println("Storing of PhyloTree");
            oos.writeObject(tree);
            Infos.println("Storing of RelaxedTree");
            if (this.relaxedTree!=null)
                oos.writeObject(relaxedTree);
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
    
    public static SessionNext load(File f) {
        try {
            long startTime = System.currentTimeMillis();
            FileInputStream fis = new FileInputStream(f);
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(fis,4096));
            int k=ois.readInt();
            int minK=ois.readInt();
            float factor=ois.readFloat();
            float stateThreshold=ois.readFloat();
            float wordThreshold=ois.readFloat();
            boolean relaxedTree=ois.readBoolean();
            SessionNext s=new SessionNext(k, minK, factor, stateThreshold, wordThreshold);
            Infos.println("Loading States");
            s.states = (States)ois.readObject();
            Infos.println("Loading Alignment");
            s.align = (Alignment)ois.readObject();
            Infos.println("Loading PhyloTree");
            s.tree = (PhyloTree)ois.readObject();
            Infos.println("Loading RelaxedTree");
            if (relaxedTree)
                s.tree = (RelaxedTree)ois.readObject();
            Infos.println("Loading of PPStats");
            s.parsedProbas = (PProbasSorted)ois.readObject();
            Infos.println("Loading of Hash");
            s.hash = (SimpleHash)ois.readObject();
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
