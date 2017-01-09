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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import tree.PhyloTree;

/**
 * 
 * @author ben
 */
public class FASTMLWrapper implements DataWrapper {

    Alignment align=null;
    PhyloTree tree=null;
    int states=4;
    
    /**
     * Used FASTML input/output files + FASTML results to fill the alignment/tree objects
     * @param inputAlign input alignment used in FASTML 
     * @param inputTree tree.newick.text output from FASTML
     * @param treeAncestor tree.ancestor.txt output from FASTML
     * @param probaMarginal prob.marginal.txt from FASTML
     */
    public FASTMLWrapper(Alignment inputAlign,PhyloTree inputTree,States states) {
        this.align=inputAlign;
        this.tree=inputTree;
        this.states=states.getStateCount();
    }
    
    
    public PProbas parseProbas(InputStream input) throws IOException {
        
        PProbas matrix=new PProbas(tree.getNodeCount(), align.getLength(), states);
        
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        boolean start=false;
        int currentPP=0;
        int lineNumber=0;
        
        while ((line=br.readLine())!=null) {
            lineNumber++;
            if (line.startsWith("node,site,")) {
                start=true;
                continue;
            }
            if (line.trim().equals("")) {
                continue;
            }
            if (line.equals("++++++++++++++++++++++++ marginal log likelihood +++++++++++++++++++++++++++++++")) {
                break;
            }
            if (start) {
                String[] infos=line.split(",");
                String nodeName=infos[0];
                //System.out.println("NodeName:"+nodeName);
                int nodeId=tree.getByName(nodeName).getId();
                int site=Integer.parseInt(infos[1])-1;
                try {
                    for (int i=2;i<infos.length;i++) {
                        //thestate order from DNStates is the same as in the FastMl output, no need to change it
                        matrix.setState(nodeId, site, i-2, Double.parseDouble(infos[i]));
                    }
                } catch (java.lang.NumberFormatException ex) {
                    Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states));
                    //for now, simple hack, set probabilities to 1/states (DNA)
                    for (int i=0;i<states;i++) {
                        matrix.setState(nodeId, site, i, (1.0/states));
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
