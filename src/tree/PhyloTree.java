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
    
    protected int nodeCount=0;
    protected int leavesCount=0;
    protected HashMap<String,PhyloNode> indexByName=null;
    protected HashMap<Integer,PhyloNode> indexById=null;
    
    
    //filled by DFS searches
    protected ArrayList<Integer> orderedLeavesIds=null;
    protected ArrayList<Integer> orderedNodesIds=null;
    protected ArrayList<Integer> orderedInternalNodesIds=null;
    protected ArrayList<String> orderedNodesLabels=null;

    public PhyloTree() {
    }

    public PhyloTree(TreeModel newModel) {
        super(newModel);
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
    
    
    
    
    
}
