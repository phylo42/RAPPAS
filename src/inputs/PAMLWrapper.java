/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.PProbas;
import core.PProbasSorted;
import core.SiteProba;
import core.States;
import etc.Infos;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import tree.NewickReader;
import tree.PhyloTree;

/**
 * PAML has a strange way to give number to the leaves and internal nodes.
 * It requires to 1 st parse the original tree (#1st tree in rst file),
 * then to parse a second tree (3rd tree in rst file ) in the same way
 * (we will use DFS pre-order) to make the correspondance between
 * original tree nodes and the node ids of PAML.
 * @author ben
 */
public class PAMLWrapper implements DataWrapper {

    Alignment align=null;
    PhyloTree tree=null;
    States states=null;
    
    //to accelerate parsin, will use this regexp
    //careful, in paml outputscientic number with negative, but also positive(!) powers.
    //ex: T(1.000000E+00) C(1.050222E-12) A(2.261192E-13) G(1.568578E-13)
    String pattern="^\\s+([0-9]+)\\s+[0-9]+\\s+[^:]+:\\s+([A-Z])\\(([0-9E\\.+-]+)\\)\\s+([A-Z])\\(([0-9E\\.+-]+)\\)\\s+([A-Z])\\(([0-9E\\.+-]+)\\)\\s+([A-Z])\\(([0-9E\\.+-]+)\\)\\s+";
    Pattern p=Pattern.compile(pattern);
    
    /**
     * Used Paml 'rst' file (all infos are in there)
     * @param inputAlign input alignment used in FASTML 
     * @param inputTree tree.newick.text output from FASTML
     * @param states
     */

    public PAMLWrapper(Alignment inputAlign, States states) {
        this.align=inputAlign;
        this.states=states;
    }
    

    
    /**
     * pare the tree from the rst file (it contains several version of the tree,
     * none with complete information, and the post probas)
     * @param input
     * @return
     * @throws IOException 
     */
    public PhyloTree parseTree(InputStream input) throws IOException {
     
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        String originalTreeString=null;
        String thirdTreeString=null;
        boolean startOriginal=false;
        boolean start3rdTree=false;
        while ((line=br.readLine())!=null) {
            if (line.trim().length()<1) {continue;}
            if (line.startsWith("Ancestral reconstruction by BASEML")) {
                startOriginal=true;
                continue;
            }
            if (line.startsWith("tree with node labels for Rod Page's TreeView")) {
                start3rdTree=true;
                continue;
            }
            if (startOriginal) {
                originalTreeString=line.replaceAll(" ", "");
                Infos.println("Original tree found: "+originalTreeString);
                startOriginal=false;
            }
            if (start3rdTree) {
                thirdTreeString=line.replaceAll(" ", "");
                Infos.println("PAML node numbering found: "+thirdTreeString);
                start3rdTree=false;
                break; //to not parse the whole file
            }        
        }
        br.close();
        //System.out.println(originalTreeString);
        //System.out.println(thirdTreeString);
        
        //the tree parsing will be done in 2 phases, the 1st tree in the rst file
        //(has branch length), then the 3rd tree in the rst file (has node ids
        //from PAML, necessary to associate posterior probas).
        
        //the 1st tree is the original input tree, with branch length
        NewickReader np=new NewickReader();
        if (originalTreeString!=null) {
            this.tree=np.parseNewickTree(originalTreeString);; //keep a reference, which will be used by parseProbas()
            //tree.displayTree();
        } 
        //now we parse the second tree
        PhyloTree thirdTree=null;
        np=new NewickReader();
        if (originalTreeString!=null) {
            thirdTree=np.parseNewickTree(thirdTreeString);; //keep a reference, which will be used by parseProbas()
            //tree.displayTree();
        } 
        //do do the same DFS pre-order in both trees to associate PAML
        //3rd tree node ids to original tree.
        Infos.println("Nodes, original Tree: "+this.tree.getLabelsByDFS());
        Infos.println("Nodes, PAML numbering: "+thirdTree.getLabelsByDFS());
        
        Pattern p= Pattern.compile("^([0-9]+)(_.*$)?");
        for (int i=0;i<tree.getNodeIdsByDFS().size();i++) {
            Integer nodeId = tree.getNodeIdsByDFS().get(i);
            String PAMLString=thirdTree.getLabelsByDFS().get(i);
            Matcher m=p.matcher(PAMLString);
            if (m.matches()) {
                this.tree.getById(nodeId).setExternalId(Integer.parseInt(m.group(1)));
                if (m.group(2)==null)
                    this.tree.getById(nodeId).setLabel(m.group(1));
            }
            //m.group(2), not used but would be '_label'
        }
        //we changed the tree ids, rebuild the indexes !
        this.tree.initIndexes();
                
        return this.tree;
    }
    
    
    public PProbas parseProbas(InputStream input, float sitePPThreshold, boolean asLog10) throws IOException {
        
        PProbas matrix=new PProbas(tree.getNodeCount(), align.getLength(), states.getNonAmbiguousStatesCount());
        
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        boolean start=false;
        int currentPP=0;
        int lineNumber=0;
        String currentNode=null;
        
        while ((line=br.readLine())!=null) {
            lineNumber++;
            if (line.startsWith("Prob distribs at nodes,")) { //start of section
                start=true;
                continue;
            }
            if (line.trim().equals("") || line.trim().equals("site  Freq   Data")) { //useless lines
                continue;
            }
            if (line.equals("Prob of best state at each node, listed by site")) { //end of section
                break;
            }
            if (start) {
                
                if (line.startsWith("Prob distribution at node ")) {
                    currentNode=line.split(" ")[4].replaceAll(",", "");
                    //System.out.println("Current node: "+currentNode);
                    continue;
                }
                
                String[] infos=line.split("[ ]+");
                
                int nodeId=tree.getByName(currentNode).getId();
                int site=Integer.parseInt(infos[1])-1;
                
                try {
                    for (int i=4;i<infos.length;i++) {
                        int stateIndex=states.stateToInt(infos[i].split("\\(")[0].charAt(0));
                        //substring(2,length-1) to remove A( and last )
                        float val= Float.parseFloat(infos[i].substring(2, infos[i].length()-1)); //substring(2,length-1) to remove last )
                        //System.out.println("Infos parsed: "+nodeId+","+site+","+stateIndex+","+val);
                        if (val<sitePPThreshold)
                            val=sitePPThreshold;
                        if (asLog10)
                            matrix.setState(nodeId, site, stateIndex, (float)Math.log10(val));
                        else
                            matrix.setState(nodeId, site, stateIndex, val);
                    }
                } catch (java.lang.NumberFormatException ex) {
                    Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states.getNonAmbiguousStatesCount()));
                    //for now, simple hack, set probabilities to 1/states (DNA)
                    for (int i=0;i<states.getNonAmbiguousStatesCount();i++) {
                        if (asLog10)
                            matrix.setState(nodeId, site, i, (float)Math.log10(1.0/states.getNonAmbiguousStatesCount()));
                        else
                            matrix.setState(nodeId, site, i, (1.0f/states.getNonAmbiguousStatesCount()));
                    }
                    //ex.printStackTrace();
                    //br.close();
                    //System.exit(1);
                }
                
                currentPP++;
            }
        }
        Infos.println( "Number of (site x nodes) for which pp were parsed: "+currentPP);
        Infos.println( "Number of (sites) for which pp were parsed: "+(0.0+currentPP/(tree.getNodeCount()-tree.getLeavesCount())));
        br.close();
        return matrix;
    }
    
    
    
    
    public PProbasSorted parseSortedProbas(InputStream input, float sitePPThreshold, boolean asLog10, int debugNodeLimit) throws IOException {
        
        PProbasSorted matrix=new PProbasSorted(tree.getNodeCount(), align.getLength(), states.getNonAmbiguousStatesCount());
        
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        boolean start=false;
        int currentPP=0;
        int lineNumber=0;
        String currentNode=null;
        
        int nodeId=-1;
        int nodeCount=0;
        ArrayList<SiteProba> probasPerSite=new ArrayList<>(states.getNonAmbiguousStatesCount());
        for (int i = 0; i < states.getNonAmbiguousStatesCount(); i++) {
            probasPerSite.add(new SiteProba());
        }

        
        while ((line=br.readLine())!=null) {
            lineNumber++;
            if (lineNumber%500000==0) {
                Infos.println("Line: >"+lineNumber);
            }
            
            if (line.startsWith("Prob distribs at nodes,")) { //start of section
                start=true;
                continue;
            }
            if (line.trim().equals("") || line.trim().equals("site  Freq   Data")) { //useless lines
                continue;
            }
            if (line.equals("Prob of best state at each node, listed by site")) { //end of section
                break;
            }
            if (start) {
                
                               
                if (line.startsWith("Prob distribution at node ")) {
                    currentNode=line.split(" ")[4].replaceAll(",", "");
                    //System.out.println("Current node: "+currentNode);
                    nodeId=tree.getByName(currentNode).getId();
                    //DEBUG
                    nodeCount++;
                    if (nodeCount>debugNodeLimit)
                        break;
                    //DEBUG
                    continue;
                }
                
                try {

                    //groups:
                    //1: 1 (int of site)
                    //2: T
                    //3: 2.295358E-01
                    //4: C
                    //5: 3.009531E-01
                    //6: A
                    //7: 2.168515E-01
                    //8: G
                    //9: 2.526595E-01
                    Matcher matcher = p.matcher(line);
                    int site=-1;

                    if(matcher.matches()) {
                        site=Integer.parseInt(matcher.group(1));

                        SiteProba sp1=new SiteProba();
                        sp1.state=states.stateToByte(matcher.group(2).charAt(0));
                        sp1.proba=Float.parseFloat(matcher.group(3));

                        SiteProba sp2=new SiteProba();
                        sp2.state=states.stateToByte(matcher.group(4).charAt(0));
                        sp2.proba=Float.parseFloat(matcher.group(5));

                        SiteProba sp3=new SiteProba();
                        sp3.state=states.stateToByte(matcher.group(6).charAt(0));
                        sp3.proba=Float.parseFloat(matcher.group(7));

                        SiteProba sp4=new SiteProba();
                        sp4.state=states.stateToByte(matcher.group(8).charAt(0));
                        sp4.proba=Float.parseFloat(matcher.group(9));

                        if (sp1.proba<sitePPThreshold)
                            sp1.proba=sitePPThreshold;
                        if (sp2.proba<sitePPThreshold)
                            sp2.proba=sitePPThreshold;
                        if (sp3.proba<sitePPThreshold)
                            sp3.proba=sitePPThreshold;
                        if (sp4.proba<sitePPThreshold)
                            sp4.proba=sitePPThreshold;

                        if (asLog10) {
                            sp1.proba=(float)Math.log10(sp1.proba);
                            sp2.proba=(float)Math.log10(sp2.proba);
                            sp3.proba=(float)Math.log10(sp3.proba);
                            sp4.proba=(float)Math.log10(sp4.proba);
                        }

                        probasPerSite.set(0,sp1);
                        probasPerSite.set(1,sp2);
                        probasPerSite.set(2,sp3);
                        probasPerSite.set(3,sp4);
                        
                        Collections.sort(probasPerSite);
                        matrix.setStates(nodeId, site-1, probasPerSite);                    
                    }

                } catch (java.lang.NumberFormatException ex) {
                    Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states.getStateCount()));
                    ex.printStackTrace();
                    br.close();
                    System.exit(1);
                }
                
                currentPP++;
            }
        }
        Infos.println( "Number of (site x nodes) for which pp were parsed: "+currentPP);
        Infos.println( "Number of (sites) for which pp were parsed: "+(0.0+currentPP/(tree.getNodeCount()-tree.getLeavesCount())));
        br.close();
        return matrix;
    }
    
    
}
