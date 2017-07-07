/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.Enumeration;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;

/**
 *
 * @author ben
 */
public class NewickWriter {
    
    private int currentNodeIndex=0;
    private int level=0;
    private int treeCount=0;
    
    private boolean branchLength=true;
    private boolean internalNodeNames=true;
    // if true, a {X} labelling edge is added to each node, this represents the
    // edge linking the node to its parent. Same ID (integer) is used for the
    //node and this parent edge.
    private boolean jplaceBranchLabels=false; 
    
    private NumberFormat format=null;
    
    private Writer w=null;

    public NewickWriter() {
        setupWriter(w);
    }

    public NewickWriter(Writer w) {
        setupWriter(w);
    }

    public NewickWriter(File f) throws IOException {
        setupWriter(new FileWriter(f));
    }
    
    /**
     * init basic stuff in the writer, such as number formats for branch length
     * @param w 
     */
    private void setupWriter(Writer w) {
        this.w=w;
        //makes all branch length as a decimal number, no scientific number
        format=NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(12);
        format.setMinimumFractionDigits(12);
        format.setParseIntegerOnly(false);
    }
    
    
 /**
     * write a phylotree as a newick tree file, can be called several times to write
     * different trees in the same newick output. For jplaceBranchLabels, the id
     * of the node is reused transfered to give an id on the parent branch,
     * following the JSON specification of jplace format (i.e. node:0.123{branch_id} ).
     * @param tree 
     * @param withInternalNodeNames 
     * @param withJplaceBranchLabels 
     * @param withBranchLength 
     * @throws java.io.IOException 
     */
    public void writeNewickTree(PhyloTree tree,boolean withBranchLength, boolean withInternalNodeNames, boolean withJplaceBranchLabels) throws IOException {
        this.branchLength=withBranchLength;
        this.internalNodeNames=withInternalNodeNames;
        this.jplaceBranchLabels=withJplaceBranchLabels;
        if (treeCount>1) {w.append('\n');}
        currentNodeIndex=0;
        level=0;
        PhyloNode root=(PhyloNode)tree.getModel().getRoot();
        StringBuilder sb = new StringBuilder();
        sb = writerDFS(root,sb);
        w.append(sb);
        treeCount++;
    }
    
 /**
     * return a phylotree as a newick tree string. One tree at each call.
     * For jplaceBranchLabels, the id of the node is reused transfered to give
     * an id on the parent branch, following the JSON specification of jplace
     * format (i.e. node:0.123{branch_id} ).
     * @param tree 
     * @param withInternalNodeNames 
     * @param withJplaceBranchLabels 
     * @param withBranchLength 
     * @throws java.io.IOException 
     */
    public String getNewickTree(PhyloTree tree,boolean withBranchLength, boolean withInternalNodeNames, boolean withJplaceBranchLabels) throws IOException {
        this.branchLength=withBranchLength;
        this.internalNodeNames=withInternalNodeNames;
        this.jplaceBranchLabels=withJplaceBranchLabels;
        currentNodeIndex=0;
        level=0;
        PhyloNode root=(PhyloNode)tree.getModel().getRoot();
        StringBuilder sb = new StringBuilder();
        sb = writerDFS(root,sb);
        return sb.toString();
    }
    

    /**
     * depth first search used by the newick writer
     * @param node 
     * @param sb
     * @param withBranchLength
     */
    private StringBuilder writerDFS(PhyloNode node, StringBuilder sb) {
        
        //start this level
        sb.append("(");
        int childrenLeft=node.getChildCount();
        //System.out.println("   current node: "+node);
        for (Enumeration e=node.children();e.hasMoreElements();) {
            //System.out.println("   childrenLeft"+childrenLeft);
            childrenLeft-=1; 
            PhyloNode currentNode=(PhyloNode)e.nextElement();
            if (currentNode.isLeaf()) {
                sb.append(currentNode.getLabel());
                if(branchLength) {
                    sb.append(':');
                    sb.append(format.format(currentNode.getBranchLengthToAncestor()));
                }
                if (jplaceBranchLabels) {
                    sb.append('{');
                    sb.append(currentNode.getId());
                    sb.append('}');
                }
                
            } else {
                //System.out.println("LEVEL "+level+" TO "+(++level));
                writerDFS(currentNode,sb);
                //System.out.println("RETURN TO "+(--level) +"  (node:"+node+")");
            }
            //return from recusion or simple leaf,
            //add ',' if more children
            //close level with current node info if no children left
            if (childrenLeft>0) {
                sb.append(',');
            } else {
                sb.append(')');
                if(internalNodeNames) {
                    sb.append(node.getLabel());
                }
                if (branchLength) {
                    sb.append(':');
                    sb.append(format.format(node.getBranchLengthToAncestor()));
                }
                if (jplaceBranchLabels) {
                    sb.append('{');
                    sb.append(currentNode.getId());
                    sb.append('}');
                }            }
        }
        //close the string with a ';' after all root children were passed
        if (node.isRoot() && childrenLeft<1) {
            sb.append(';'); 
        }
        
        //System.out.println(sb+"\n");
        return sb;
    }
    
    /**
     * close writer 
     * @throws IOException 
     */
    public void close() throws IOException {
        w.flush();
        w.close();
    }
    
    

    
    /**
     * test main
     * @param args 
     */
    public static void main(String[] args) {
        
        String treeFASTML=  "((JN856008:0.035546,AF124992:0.018812)N2:0.011794,(AY874541:0.011802,((((AB907632:0.012757,AB907633:0.002090)N7:0.010082,(((EU769559:0" +
                            ".006293,(EU769560:0.002622,AB907631:0.012287)N11:0.002049)N10:0.008871,AF124986:0.005430)N9:0.008257,((AB781796:0.020279,((AB907634:0." +
                            "022758,((((KJ473819:2.180129,((KJ473805:0.000001,KJ473804:0.020913)N21:0.142716,KF843851:0.119877)N20:0.116703)N19:0.099975,(KF843858:" +
                            "0.484028,KF843855:0.147765)N22:0.086147)N18:0.416505,AB781795:0.002254)N17:0.039562,AB781793:0.013682)N16:0.005580)N15:0.005506,KF5301" +
                            "23:0.029973)N14:0.011016)N13:0.081245,AF516906:0.022711)N12:0.006863)N8:0.002245)N6:0.001145,(DQ431016:0.005662,DQ431014:0.005897)N23:" +
                            "0.011966)N5:0.008440,((AB907625:0.002438,((AB907630:0.006924,(AB907628:0.002390,AB781792:0.002455)N28:0.002897)N27:0.005056,AB781791:0" +
                            ".005370)N26:0.004328)N25:0.000001,AY878324:0.002363)N24:0.005962)N4:0.004490)N3:0.001528,KC175339:0.014500)N1:0.0;";

        String tree="(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;";

        String treePAML=    "(5_KJ473819,((13_KF843851,(10_KJ473805,4_KJ473804)34)33,((12_KF843855,14_KF843858)36,(26_AB781795,"
                + "(23_AB781793,(17_AB907634,(2_KF530123,(25_AB781796,(11_AF516906,((6_AF124986,(8_EU769559,(9_EU769560,21_AB907631)46)45)44,"
                + "((20_AB907632,22_AB907633)48,((30_DQ431016,29_DQ431014)50,((27_AY874541,(1_KC175339,(3_JN856008,7_AF124992)54)53)52,"
                + "(28_AY878324,(16_AB907625,(15_AB781791,(19_AB907630,(18_AB907628,24_AB781792)59)"
                + "58)57)56)55)51)49)47)43)42)41)40)39)38)37)35)32)31;";
        
        

        
        //test parsing
        long startTime = System.currentTimeMillis();
        PhyloTree t=new NewickReader().parseNewickTree2(treeFASTML, false);
        long endTime = System.currentTimeMillis();
        System.out.println("Parsing took " + (endTime - startTime) + " milliseconds");

        //test writing
        startTime = System.currentTimeMillis();
        try {
            NewickWriter nw=new NewickWriter(new File("test_fake.tree"));
            nw.writeNewickTree(t,true,true,true);
            System.out.println(nw.getNewickTree(t, true, true, true));
            nw.close();
            
        } catch (IOException ex) {
            Logger.getLogger(NewickReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        endTime = System.currentTimeMillis();
        System.out.println("Writing took " + (endTime - startTime) + " milliseconds");
        
        
        
        JFrame f=new JFrame();
        t.setSize(700,700);
        t.setEnabled(true);
        t.setVisible(true);
        f.setSize(800, 800);
        f.add(t);
        f.setVisible(true);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
    
}
