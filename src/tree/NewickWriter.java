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
import java.util.Locale;

/**
 *
 * @author ben
 */
public class NewickWriter {
    
    private int currentNodeIndex=0;
    private int level=-1;
    private int treeCount=0;
    
    private boolean branchLength=true;
    private boolean internalNodeNames=true;
    // if true, a {X} labelling edge is added to each node, this represents the
    // edge linking the node to its parent. Same ID (integer) is used for the
    //node and this parent edge.
    private boolean jplaceBranchLabels=false; 
    //if true, the programmatic nodeid is added as a prefix to node labels
    //ex: (8BV-4:0.023443493992,3BV-52:0.010550861247):0.011659557000)
    //    (__12__8BV-4:0.023443493992,__13__3BV-52:0.010550861247)__11__:0.011659557000)
    private boolean nodeIdPrefix=false;
    
    private NumberFormat format=null;
    
    private Writer w=null;

    public NewickWriter() {
        setupWriter(null);
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
        //and with . as fraction separator
        format=NumberFormat.getNumberInstance(Locale.UK);
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
     * @param withBranchLength 
     * @param withInternalNodeNames 
     * @param withJplaceBranchLabels 
     * @param withNodeIdPrefix 
     * @throws java.io.IOException 
     */
    public void writeNewickTree(PhyloTree tree, boolean withBranchLength, boolean withInternalNodeNames, boolean withJplaceBranchLabels, boolean withNodeIdPrefix) throws IOException {
        this.branchLength=withBranchLength;
        this.internalNodeNames=withInternalNodeNames;
        this.jplaceBranchLabels=withJplaceBranchLabels;
        this.nodeIdPrefix=withNodeIdPrefix;
        if (treeCount>1) {w.append('\n');}
        currentNodeIndex=0;
        //the level variable is used to not set any branch length or
        //jplace edge ids at the virtual root of an unrooted newick
        //i.e.  (A,B,C):0.0{0};  => not valid
        //      (A,B,C);         => valid
        if (!tree.isRooted()) {
            level=-1;
        } else {
            level=0;
        }
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
     * @param withBranchLength 
     * @param withInternalNodeNames 
     * @param withJplaceBranchLabels
     * @param withNodeIdPrefix the value of withNodeIdPrefix 
     * @throws java.io.IOException
     * @return the java.lang.String 
     */
    public String getNewickTree(PhyloTree tree, boolean withBranchLength, boolean withInternalNodeNames, boolean withJplaceBranchLabels, boolean withNodeIdPrefix) throws IOException {
        this.branchLength=withBranchLength;
        this.internalNodeNames=withInternalNodeNames;
        this.jplaceBranchLabels=withJplaceBranchLabels;
        this.nodeIdPrefix=withNodeIdPrefix;
        currentNodeIndex=0;
        //the level variable is used to not set any branch length or
        //jplace edge ids at the virtual root of an unrooted newick
        //i.e.  (A,B,C):0.0{0};  => not valid
        //      (A,B,C);         => valid
        if (!tree.isRooted()) {
            level=-1;
        } else {
            level=0;
        }
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
        //System.out.println("sb: "+sb);
        //start this level
        sb.append("(");
        int childrenLeft=node.getChildCount();
        Enumeration e=node.children();
        while (e.hasMoreElements()) {
            //System.out.println("   childrenLeft"+childrenLeft);
            childrenLeft-=1; 
            PhyloNode currentNode=(PhyloNode)e.nextElement();
            //System.out.println("currentNode "+currentNode);
            if (currentNode.isLeaf()) {
                if (nodeIdPrefix) {
                    sb.append("__");
                    sb.append(currentNode.getId());
                    sb.append("__");
                }
                sb.append(currentNode.getLabel());
                if(branchLength) {
                    sb.append(':');
                    sb.append(format.format(currentNode.getBranchLengthToAncestor()));
                }
                if (jplaceBranchLabels) {
                    sb.append('{');
                    sb.append(currentNode.getJplaceEdgeId());
                    sb.append('}');
                }
            } else {
                //System.out.println("LEVEL "+level+" TO "+(++level));
                level++;
                writerDFS(currentNode,sb);
                level--;
                //System.out.println("RETURN TO "+(--level) +"  (node:"+node+")");
            }
            //return from recursion or simple leaf,
            //add ',' if more children
            //close level with current node info if no children left
            if (childrenLeft>0) {
                sb.append(',');
            } else {
                sb.append(')');
                if (nodeIdPrefix) {
                    sb.append("__");
                    sb.append(currentNode.getId());
                    sb.append("__");
                }
                if(internalNodeNames) {
                    sb.append(node.getLabel());
                }
                if (branchLength && (level>-1) ) {
                    sb.append(':');
                    sb.append(format.format(node.getBranchLengthToAncestor()));
                }
                if (jplaceBranchLabels && (level>-1)) {
                    sb.append('{');
                    sb.append(node.getJplaceEdgeId());
                    sb.append('}');
                }            
            }

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
        if (w!=null) {
            w.flush();
            w.close();
        }
    }
  
}
