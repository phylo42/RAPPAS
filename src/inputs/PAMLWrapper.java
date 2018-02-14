/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.PProbasSorted;
import core.SiteProba;
import core.States;
import etc.Infos;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import tree.NewickReader;
import tree.PhyloTree;

/**
 * PAML has a strange way to give number to the leaves and internal nodes.
 * It requires to 1 st parse the original tree (#1st tree in rst file),
 * then to parse a second tree (3rd tree in rst file ) in the same way
 * (we will use DFS pre-order) to make the match between
 * original tree nodes and the node ids modified by PAML.
 * @author ben
 */
public class PAMLWrapper implements ARWrapper {

    Alignment align=null;
    PhyloTree tree=null;
    States states=null;
    
    //intermediat tree loaded during parsing (PAML output modified trees)
    PhyloTree firstTree=null;
    PhyloTree thirdTree=null;
    Pattern pNodeMatch= Pattern.compile("^([0-9]+)(_.*$)?"); //used to match PAML 1rst and 3rd tree
    
    
    //to accelerate posterior probas parsing (rst file), we use this regexp
    //careful, in paml scientic number have negative, but also positive powers.('+' symbol present!)
    //ex: '  1  2  dfd-g: T(1.000000E+00) C(1.050222E-12) A(2.261192E-13) G(1.568578E-13)'
    //groups: 2=site 3=freq 4=state(char) 5=PP
    //note: 4 and 5 null after 1st state
    String pattern="(^\\s+([0-9]+)\\s+([0-9]+).*:\\s+)?([A-Z])\\(([0-9E\\.+-]+)\\)";
    Pattern p=Pattern.compile(pattern);
    
    /**
     * Used Paml 'rst' file (all infos are in there)
     * @param inputAlign input alignment used in FASTML 
     * @param states
     */

    public PAMLWrapper(Alignment inputAlign, States states) {
        this.align=inputAlign;
        this.states=states;
    }
    

    
    /**
     * parse the tree from the rst file (it contains several version of the tree,
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
            if (line.startsWith("Ancestral reconstruction by")) {
                startOriginal=true;
                continue;
            }
            if (line.startsWith("tree with node labels for Rod Page's TreeView")) {
                start3rdTree=true;
                continue;
            }
            if (startOriginal) {
                originalTreeString=line.replaceAll(" ", "");
                //Infos.println("1st PAML tree found: "+originalTreeString);
                Infos.println("1st PAML tree found");
                startOriginal=false;
            }
            if (start3rdTree) {
                thirdTreeString=line.replaceAll(" ", "");
                //Infos.println("3rd PAML tree found: "+thirdTreeString);
                Infos.println("3rd PAML tree found");
                start3rdTree=false;
                break; //to not parse the whole file
            }        
        }
        br.close();
        
        //the tree parsing will be done in 2 phases, the 1st tree in the rst file
        //(has branch length), then the 3rd tree in the rst file (has node ids
        //from PAML, necessary to associate posterior probas).
        
        //the 1st tree is the original input tree, with branch length but 
        //no internal nodes and renaming of the tips
        if (originalTreeString!=null) {
            this.firstTree=NewickReader.parseNewickTree2(originalTreeString, false, false); //keep a reference, which will be used by parseProbas()
            //tree.displayTree();
        } 
        //now we parse the second tree
        if (thirdTreeString!=null) {
            this.thirdTree=NewickReader.parseNewickTree2(thirdTreeString, false, false); //keep a reference, which will be used by parseProbas()
            //tree.displayTree();
        } 
        //do do the same DFS pre-order in both trees to associate PAML
        //3rd tree node ids to original tree.
        //Infos.println("Nodes, first tree: "+firstTree.getLabelsByDFS());
        //Infos.println("Nodes, third tree: "+thirdTree.getLabelsByDFS());
        

        for (int i=0;i<firstTree.getNodeIdsByDFS().size();i++) {
            Integer nodeId = firstTree.getNodeIdsByDFS().get(i);
            String PAMLString=thirdTree.getLabelsByDFS().get(i);
            Matcher m=pNodeMatch.matcher(PAMLString);
            if (m.matches()) {
                //this.firstTree.getById(nodeId).setExternalId(Integer.parseInt(m.group(1)));
                if (m.group(2)==null)
                    this.firstTree.getById(nodeId).setLabel(m.group(1));
            }
            //m.group(2), not used but would be '_label'
        }
        //we changed the tree ids, rebuild the indexes !
        this.firstTree.initIndexes();    
        
        //we set the modified 1st tree as the tree supported by this wrapper
        this.tree=firstTree;
        return this.tree;
    }
    
    /**
     * parse the posterior probas themselves
     * @param input
     * @param sitePPThreshold
     * @param asLog10
     * @param debugNodeLimit limit PP parsing to few nodes
     * @return
     * @throws IOException 
     */
    public PProbasSorted parseSortedProbas(InputStream input, float sitePPThreshold, boolean asLog10, int debugNodeLimit){
        
        BufferedReader br=null;
        int lineNumber=0;
        PProbasSorted matrix=null;
        try {
            matrix=new PProbasSorted(tree.getNodeCount(), align.getLength(), states.getNonAmbiguousStatesCount());
            br = new BufferedReader(new InputStreamReader(input,"UTF-8"));
            String line=null;
            boolean start=false;
            int currentPP=0;
            String currentNode=null;
            int nodeId=-1;
            int nodeCount=0;
            ArrayList<SiteProba> probasPerSite=new ArrayList<>(states.getNonAmbiguousStatesCount());
            for (int i = 0; i < states.getNonAmbiguousStatesCount(); i++) {
                probasPerSite.add(new SiteProba());
            }   Infos.println("Starting to parse PAML posterior probas...");
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
                    

                    //groups: see Pattern declaration
                    Matcher matcher = p.matcher(line);
                    int site=-1;
                    int stateIdx=0;
                    while (matcher.find()) {
                        if (matcher.group(2)!=null) {
                            site=Integer.parseInt(matcher.group(2));
                        }
                        SiteProba sp=new SiteProba();
                        sp.state=states.stateToByte(matcher.group(4).charAt(0));
                        sp.proba=Float.parseFloat(matcher.group(5));  
                        if (sp.proba<sitePPThreshold)
                            sp.proba=sitePPThreshold;
                        
                        if (asLog10) {
                            sp.proba=(float)Math.log10(sp.proba);
                        }

                        probasPerSite.set(stateIdx,sp);
                        stateIdx++;
                        
                    }
                    Collections.sort(probasPerSite);
                    matrix.setStates(nodeId, site-1, probasPerSite);

                    /*OLD VERSION FOR DNA ONLY
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
                    }*/
                        
                    currentPP++;
                }
            }
            Infos.println( "Number of (site x nodes) for which pp were parsed: "+currentPP);
            Infos.println( "Number of (sites) for which pp were parsed: "+(0.0+currentPP/(tree.getNodeCount()-tree.getLeavesCount())));
            
        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (java.lang.NumberFormatException ex) {
            Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states.getStateCount()));                    
            ex.printStackTrace();
            System.exit(1);
        } catch (IOException ex) {
            Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return matrix;
    }
    
    
}
