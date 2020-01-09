package inputs;

import alignement.Alignment;
import core.PProbasSorted;
import core.SiteProba;
import core.States;
import etc.Infos;
import etc.exceptions.NonSupportedStateException;
import tree.NewickReader;
import tree.PhyloTree;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

public class RAXMLNGWrapper implements ARWrapper {

    Alignment align=null;
    PhyloTree tree=null;
    States states=null;

    /**
     *
     * @param ARInputAlign input alignment used at ancestral reconstruction
     * @param states
     */

    public RAXMLNGWrapper(Alignment ARInputAlign, States states) {
        this.align=ARInputAlign;
        this.states=states;
    }


    /**
     * parse the tree from the raxml.ancestralTree file
     * @param input
     * @param rerootRAxMLngTree
     * @return
     * @throws IOException
     */
    public PhyloTree parseTree(InputStream input,boolean rerootRAxMLngTree) throws IOException {
        BufferedReader br=new BufferedReader(new InputStreamReader(input,"UTF-8"));
        String line=null;
        String originalTreeString=null;
        while ((line=br.readLine())!=null) {
            if (line.trim().length()<1) {continue;}
            originalTreeString=line;
        }
        //if necessary reroot phyml tree
        if (rerootRAxMLngTree) {
            Infos.println("Reversing phyML unrooting...");

            //first change newick string to get root sons
            //order from (C3,C1,C2); to (C1,C2,C3);

            //1st, define root
            int cladeClosingIndex=-1;
            for (int i = originalTreeString.length()-1; i >= 0; i--) {
                if (originalTreeString.charAt(i)==')') {
                    cladeClosingIndex=i;
                    break;
                }
            }
            //System.out.println("Closing index: "+cladeClosingIndex);

            //2nd, extract C1 to C3
            String[] clades=new String[4]; //the 4th contains the root node
            int depth=0;
            int cladeStart=1;
            int cladeCounter=0;
            for (int i = 0; i < originalTreeString.length(); i++) {
                char c=originalTreeString.charAt(i);
                if (c=='(') { depth++; }
                if (c==')') { depth--; }
                //we arrive or return to clade attached to roto
                if ( (depth==1 && c==',') || (depth==0 && i==cladeClosingIndex) ) {
                    //store previous clade string
                    if (i>0) {
                        clades[cladeCounter]=originalTreeString.substring(cladeStart,i);
                        //System.out.println(clades[cladeCounter]);
                        cladeCounter++;
                    }
                    cladeStart=i+1;

                    continue;
                }

            }
            //last one = the root
            clades[cladeCounter]=originalTreeString.substring(cladeStart,originalTreeString.length());
            //System.out.println(clades[cladeCounter]);

            //reorder
            originalTreeString="("+clades[1]+","+clades[2]+","+clades[0]+")"+clades[3];
            //System.out.println("---\n"+originalTreeString+"\n---\n");


            //then force rooting, which goes back to the original extended tree rooting
            this.tree=NewickReader.parseNewickTree2(originalTreeString, true, false);

        } else {
            this.tree=NewickReader.parseNewickTree2(originalTreeString, false, false);
        }
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
    @Override
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
            int previousId=-1;
            int nodeCount=0;
            List<Character> stateOrder=new ArrayList<>(states.getNonAmbiguousStatesCount());
            ArrayList<SiteProba> probasPerSite=new ArrayList<>(states.getNonAmbiguousStatesCount());
            for (int i = 0; i < states.getNonAmbiguousStatesCount(); i++) {
                probasPerSite.add(new SiteProba());
            }
            Infos.println("Starting to parse RAXMLNG posterior probas...");
            while ((line=br.readLine())!=null) {
                lineNumber++;
                if (lineNumber%500000==0) {
                    Infos.println("Line: >"+lineNumber);
                }
                if (line.startsWith("Node\tSite\tState")) { //start of section
                    //header is Node\tSite\tState\tp_A\tp_C\tp_G\tp_T\n
                    String[] stateStrings=line.split("\t");
                    for (int i = 3; i < stateStrings.length; i++) { //skip 1st to 3rd columns
                        if (!stateStrings[i].trim().equals("MPEE")) {
                            stateOrder.add(stateStrings[i].trim().charAt(2));
                        }
                    }
                    Infos.println("States found:"+stateOrder);
                    start=true;
                    continue;
                }
                if (line.trim().isEmpty()) { //useless lines
                    continue;
                }

                if (start) {

                    //groups:
                    //Node\tSite\tState\tp_A\tp_C\tp_G\tp_T\n
                    String[] data=line.split("\t");

                    //node parsing
                    currentNode=data[0].trim();
                    nodeId=tree.getByName(currentNode).getId();
                    if (nodeId!=previousId) {
                        nodeCount++;
                        //DEBUG
                        if (nodeCount>debugNodeLimit) {
                            break;
                            //DEBUG
                        }
                    }
                    previousId=nodeId;

                    //site and probas parsing
                    int site=Integer.parseInt(data[1].trim());
                    if (site>align.getLength()) {
                        System.out.println("It seems the RAxML-ng AR output contains more sites than the input reference alignment.");
                        System.out.println("Was the AR performed on the same alignment ?");
                        System.exit(1);
                    }

                    //System.out.println("NODE/SITE: "+nodeId+" "+site);
                    for (int i = 0; i < stateOrder.size(); i++) {
                        SiteProba sp=new SiteProba();
                        try {
                            sp.state=states.stateToByte(stateOrder.get(i));
                        } catch ( NonSupportedStateException ex) {
                            ex.printStackTrace(System.err);
                            System.out.println("RAxML-ng wrapper encountered a non supported state. (state="+stateOrder.get(i)+")");
                            System.exit(1);
                        }
                        sp.proba=Float.parseFloat(data[i+3]);
                        //System.out.print(" "+sp.state+":"+sp.proba);
                        if (sp.proba<sitePPThreshold)
                            sp.proba=sitePPThreshold;
                        if (asLog10)
                            sp.proba=(float)Math.log10(sp.proba);
                        probasPerSite.set(i,sp);
                    }
                    //System.out.println("");
                    //sort
                    Collections.sort(probasPerSite);
                    //System.out.println(probasPerSite);
                    //register
                    matrix.setStates(nodeId, site-1, probasPerSite);

                    currentPP++;
                }
            }
            Infos.println( "Number of (site x nodes) for which pp were parsed: "+currentPP);
            Infos.println( "Number of (sites) for which pp were parsed: "+(0.0+currentPP/(tree.getNodeCount()-tree.getLeavesCount())));

        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (java.lang.NumberFormatException ex) {
            Infos.println("Parsing error (line "+lineNumber+"): all states will have pp="+(1.0/states.getNonAmbiguousStatesCount()));
            ex.printStackTrace();
            System.exit(1);
        } catch (IOException ex) {
            Logger.getLogger(PAMLWrapper.class.getName()).log(Level.SEVERE, null, ex);
        }
        return matrix;
    }


}
