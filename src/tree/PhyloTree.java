/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import etc.Infos;
import java.awt.Dimension;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;

/**
 * A tree of PhyloNodes, based on the original java structure.
 * @author ben
 */
public class PhyloTree extends JTree implements Serializable {
    
    private static final long serialVersionUID = 2000L;
    
    private boolean isRooted=false;
    private int nodeCount=0;
    private int leavesCount=0;
    //simple map of pointers to directly access the nodes by name/id
    //whithout having to run along the tree again
    private HashMap<String,PhyloNode> indexByName=null;
    private HashMap<Integer,PhyloNode> indexById=null;
    
    
    //filled by a single DF search after calling init(),
    //avoids to do it again in other program modules
    private ArrayList<Integer> orderedLeavesIds=null;
    private ArrayList<Integer> orderedNodesIds=null;
    private ArrayList<Integer> orderedInternalNodesIds=null;
    private ArrayList<String> orderedNodesLabels=null;

    //necessary to use the Extended tree specialization
    public PhyloTree() {}
    

    protected PhyloTree(TreeModel newModel, boolean rooted) {
        super(newModel);
        this.isRooted=rooted;
        
    }
    
    public PhyloNode getByName(String nodeName) {
        return indexByName.get(nodeName);
    }
    
    public PhyloNode getById(int id) {
        return indexById.get(id);
    }
    
    /**
     * Carefull! don't forget to relaunch initIndexes() after making several 
     * ids updates
     * @param oldId
     * @param newId 
     */
    public void updateId(int oldId, int newId) {
        if (!indexById.containsKey(oldId)) {
            Infos.println("Id don't exists and cannot be updated.");
            return;
        }
        PhyloNode n=getById(oldId);
        
    }
    
    
    /**
     * count including internal nodes and leaves
     * @return 
     */
    public int getNodeCount() {
        return nodeCount;
    }
    
    /**
     * count including leaves only
     * @return 
     */
    public int getLeavesCount() {
        return leavesCount;
    }
    
    public PhyloNode getRoot() {
        return (PhyloNode)this.getModel().getRoot();
    }
    
    
    
    /**
     * as seen during newick parsing (2 or 3 sons at highest newick level)
     * @return 
     */
    public boolean isRooted() {
        return this.isRooted;
    }

    /**
     * for graphical representation only
     * @param tree
     * @param expand 
     */
    public void expandTree(JTree tree, boolean expand) {
        TreeNode root = (TreeNode) tree.getModel().getRoot();
        expandAll(tree, new TreePath(root), expand);
    }

    /**
     * expand all nodes of the JTree
     * @param tree
     * @param path
     * @param expand 
     */
    protected void expandAll(JTree tree, TreePath path, boolean expand) {
        TreeNode node = (TreeNode) path.getLastPathComponent();

        if (node.getChildCount() >= 0) {
            Enumeration enumeration = node.children();
            while (enumeration.hasMoreElements()) {
                TreeNode n = (TreeNode) enumeration.nextElement();
                TreePath p = path.pathByAddingChild(n);

                expandAll(tree, p, expand);
            }
        }

        if (expand) {
            tree.expandPath(path);
        } else {
            tree.collapsePath(path);
        }
    }
    
    /**
     * internal nodes ordered through depth-first search from the root
     * @return a list of nodeIds
     */
    public ArrayList<Integer> getInternalNodesByDFS() {
        return orderedInternalNodesIds;
    }
    
    /**
     * leaves ordered through depth-first search from the root
     * @return 
     */
    public ArrayList<Integer> getLeavesByDFS() {
        return orderedLeavesIds;
    }
    
    /**
     * all nodes, ordered through depth-first search from the root
     * @param node
     * @return 
     */
    public ArrayList<Integer> getNodeIdsByDFS() {
        return orderedNodesIds;
    }
    
    /**
     * all nodes, ordered through depth-first search from the root
     * @param node
     * @return 
     */
    public ArrayList<String> getLabelsByDFS() {
        return orderedNodesLabels;
    }
    
    /**
     * init the indexes (node by name, id...)
     */
    public void initIndexes() {
        Infos.println("Building tree indexation");
        this.indexByName=new HashMap<>();
        this.indexById=new HashMap<>();
        this.orderedLeavesIds = new ArrayList<>();
        this.orderedNodesIds = new ArrayList<>();
        this.orderedNodesLabels = new ArrayList<>();
        this.orderedInternalNodesIds = new ArrayList<>();
        this.nodeCount=0;
        this.leavesCount=0;
        dfs((PhyloNode)this.getModel().getRoot());
    }

    /**
     * depth first search from v
     * @param node 
     */
    private void dfs(PhyloNode node) {
        nodeCount++;
        indexById.put(node.getId(), node);
        indexByName.put(node.getLabel(), node);
        //report leaves encountered
        if (node.isLeaf()) {
            leavesCount++;
            orderedLeavesIds.add(node.getId());
        } else {
            orderedInternalNodesIds.add(node.getId());
        }
        orderedNodesIds.add(node.getId());
        orderedNodesLabels.add(node.getLabel());
        //go down recursively
        for (int i=0;i<node.getChildCount();i++) {
            PhyloNode child=(PhyloNode) node.getChildAt(i);
            dfs(child);
        }
    }
    
    
    /**
     * very basic tree visualizer, using the native swing libraries
     * @param t
     */
    public void displayTree() {
        JFrame f=new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        JScrollPane sp=new JScrollPane(this, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        sp.setPreferredSize(new Dimension(800, 800));
        setSize(700,700);
        setEnabled(true);
        setVisible(true);
        setRootVisible(true);
        expandTree(this, true);
        f.setSize(800, 800);
        f.add(sp);
        f.setVisible(true);
        
    }
    
    
    /**
     * return a map making the link between the nodeIds of this tree and another
     * tree provided by the user: the mapping is based on the leaves labels,
     * and a Depth First Traversals search. Nodes are mapped when going up
     * @param treeOriginal
     * @param <error>
     * @return 
     */
    
    
    //--------------------------------------------------------------------------
    //This block of methods is used to map the internal lables of a tree
    //of similar topology than this tree (rooted or unrooted) 
    HashMap<Integer,Integer> nodeMapping=new HashMap<>();
    
    /**
     * map(this.tree nodeId)=otherTree nodeID
     * @param otherTree
     * @return 
     */
    public HashMap<Integer,Integer> mapNodes(PhyloTree otherTree) {
        DFSforMapping(this.getRoot(),otherTree);
        //finally, if both tree rooted on same edge, associate the root ids
        //both trees rooted ?
        if (this.isRooted() && otherTree.isRooted) {
            System.out.println("both trees are rooted");
            //rooted on same edge ? if yes roots can be mapped between
            //both trees, if not no mapping is done for the root
            PhyloNode root=this.getRoot();
            Enumeration e=root.children();
            boolean sameChildren=true;
            while (e.hasMoreElements()) {
                PhyloNode child=(PhyloNode)e.nextElement();
                int idInOtherTree=nodeMapping.get(child.getId());
//                System.out.println("test "+child.getId()+" "+idInOtherTree);
//                System.out.println("A:"+otherTree.getById(idInOtherTree).getParent());
//                System.out.println("B:"+otherTree.getRoot());
                if (otherTree.getById(idInOtherTree).getParent() != otherTree.getRoot()) {
                    sameChildren=false;
                }
            }
            if (sameChildren)
                nodeMapping.put(this.getRoot().getId(), otherTree.getRoot().getId());
        }
        
        return nodeMapping;
    }
    
    /**
     * do a depth first transversal to go to all nodes, starts from a bait node 
     * @param node the bait node
     */
    private void DFSforMapping(PhyloNode node,PhyloTree otherTree) {
        
        //go down recursively is some children
        Enumeration e=node.children();        
        while (e.hasMoreElements()) {
            PhyloNode n=(PhyloNode)e.nextElement();
            DFSforMapping(n,otherTree);
        }
        //if root, no parent to associate
        if (node.isRoot()) {
            System.out.println("SKIP ROOT");
            return;
        }
        //System.out.println("IN: "+node);
        
        //returns up, time to associate associate nodes, 3 possible case:
        //1. node child of root, 
        //we will associate parent only if both trees are rooted on same branch
        //this is verified after completing the mapping, in th mapnodes() method
        if (((PhyloNode)node.getParent())==this.getRoot()) {
            if (node.isLeaf()) {
                PhyloNode otherLeaf=otherTree.getByName(node.getLabel());
                nodeMapping.put(node.getId(), otherLeaf.getId());           
            }
        //2. a leaf, not problem to map parent
        } else if (node.isLeaf()) {
            //associate leaf
            PhyloNode otherLeaf=otherTree.getByName(node.getLabel());
            nodeMapping.put(node.getId(), otherLeaf.getId());
            //associate parent of leaf
            PhyloNode parent=(PhyloNode)node.getParent();
            PhyloNode otherParent=(PhyloNode)otherLeaf.getParent();
            nodeMapping.put(parent.getId(), otherParent.getId());
        //3. internal node, reuse previous mapping to associate parent
        } else {
            int nodeId=node.getId();
            int otherNodeId=nodeMapping.get(nodeId);
            PhyloNode parent=(PhyloNode)node.getParent();
            PhyloNode otherParent=(PhyloNode)otherTree.getById(otherNodeId).getParent();
            nodeMapping.put(parent.getId(), otherParent.getId());
        }
        
        //System.out.println(nodeMapping);
    }
    //--------------------------------------------------------------------------

    
}
