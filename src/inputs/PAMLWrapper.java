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
import tree.NewickParser;
import tree.PhyloTree;

/**
 *
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
    
    public PhyloTree parseTree(InputStream input) throws IOException {
     
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        String treeString=null;
        boolean start=false;
        int lineNumber=0;
        while ((line=br.readLine())!=null) {
            lineNumber++;
            if (line.startsWith("tree with node labels for Rod Page's TreeView")) {
                start=true;
                continue;
            }
            if (start) {
                treeString=line.replaceAll(" ", "");
                Infos.println("Tree found: "+treeString);
                break;
            }
        }
        br.close();
        NewickParser np=new NewickParser();
        if (treeString!=null) {
            this.tree=np.parseNewickTree(treeString);; //keep a reference, which will be used by parseProbas()
            //tree.displayTree();
            return tree;
        } else {
            return null;
        }
    }
    
@Override
    public PProbas parseProbas(InputStream input) throws IOException {
        
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
