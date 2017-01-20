/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import alignement.Alignment;
import core.PProbas;
import core.States;
import etc.Infos;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
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
    
    @Override
    public PProbas parseProbas(InputStream input, double sitePPThreshold) throws IOException {
        
        PProbas matrix=new PProbas(tree.getNodeCount(), align.getLength(), states.getStateCount());
        
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
                        double val=Double.parseDouble(infos[i].substring(2, infos[i].length()-1)); //substring(2,length-1) to remove last )
                        //System.out.println("Infos parsed: "+nodeId+","+site+","+stateIndex+","+val);
                        if (val<sitePPThreshold)
                            val=sitePPThreshold;
                        matrix.setState(nodeId, site, stateIndex, val);
                    }
                } catch (java.lang.NumberFormatException ex) {
                    Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states.getStateCount()));
                    //for now, simple hack, set probabilities to 1/states (DNA)
                    for (int i=0;i<states.getStateCount();i++) {
                        matrix.setState(nodeId, site, i, (1.0/states.getStateCount()));
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
    
    
    
    
}
