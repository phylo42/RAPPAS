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
import java.util.LinkedHashMap;
import java.util.List;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;

/**
 * A tree of PhyloNodes, based on the original java structure.
 * @author ben
 */
public class PhyloTree extends JTree implements Serializable {
    
    private static final long serialVersionUID = 2000L;
    
    protected boolean isRooted=false;
    protected boolean isJplaceType=false;
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
    
    //only if tree is jplace-related tree
    //this map will be not null and contain the mapping between the jplace
    //edge ids and the current PhyloTree node ids (edge data holded by son node)
    //map(jplace_edge_id)=PhyloTree_node_id
    HashMap<Integer, Integer> jPlaceEdgeMappingJPToNodeID=null;
    //map(PhyloTree_node_id)=jplace_edge_id
    HashMap<Integer, Integer> jPlaceEdgeMappingNodeIDToJP=null;
    int jplaceIdsCounter=0;

    /**
     * necessary to use the ExtendedTree specialization
     */
    public PhyloTree() {}
    
    /**
     * constructor copy, used only to copy tree structure before pruning experiments
     *  or at rerooting/unrooting of trees
     * @param newModel the value of newModel
     * @param isRooted the value of isRooted
     * @param isFromJPlace the value of isFromJPlace
     */
    public PhyloTree(TreeModel newModel, boolean isRooted, boolean isFromJPlace) {
        super(newModel);
        this.isRooted=isRooted;
        this.isJplaceType=isFromJPlace;
    }
    
    /**
     * return node using name mapping
     * @param nodeName
     * @return 
     */
    public PhyloNode getByName(String nodeName) {
        return indexByName.get(nodeName);
    }
    /**
     * return node using id mapping
     * @param id
     * @return 
     */
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
     * return root of the TreeModel
     * @return 
     */
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
     * true if this PhyloTree was build from a jplace file (presence of {x} 
     * edge labels
     * @return 
     */
    public boolean isFromJplace() {
        return isJplaceType;
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
        this.jPlaceEdgeMappingJPToNodeID=new HashMap<>();
        this.jPlaceEdgeMappingNodeIDToJP=new HashMap<>();
        dfs((PhyloNode)this.getModel().getRoot());     
    }

    /**
     * depth first search, starting from node
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
        if (isJplaceType)
            jPlaceEdgeMappingJPToNodeID.put(node.getJplaceEdgeId(), node.getId());
        //go down recursively
        for (int i=0;i<node.getChildCount();i++) {
            PhyloNode child=(PhyloNode) node.getChildAt(i);
            dfs(child);
        }
    }
    
    
    /**
     * very basic tree visualiser, using the native swing libraries
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
     * CAN BE LAUNCHED ONLY AFTER initindexes() is called,
     * update the jplace edge ids to be similar to EPA/PPlacer orders
     * 
     */
    public void resetJplaceEdgeIds() {
        jplaceIdsCounter=0;
        jPlaceEdgeMappingJPToNodeID=new HashMap<>();
        jPlaceEdgeMappingNodeIDToJP=new HashMap<>();
        jplaceEdgeIdsUpdaterDFS(this.getRoot());
    }
    
     /**
     * depth first search to update the jplace ids
     * @param node 
     */
    private void jplaceEdgeIdsUpdaterDFS(PhyloNode node) {
        //start this level
        int childrenLeft=node.getChildCount();
        Enumeration e=node.children();
        while (e.hasMoreElements()) {
            childrenLeft-=1; 
            PhyloNode currentNode=(PhyloNode)e.nextElement();
            if (currentNode.isLeaf()) {
                jPlaceEdgeMappingJPToNodeID.put(jplaceIdsCounter, currentNode.getId());
                jPlaceEdgeMappingNodeIDToJP.put(currentNode.getId(),jplaceIdsCounter);
                currentNode.setJPlaceEdgeId(jplaceIdsCounter++);
            } else {
                jplaceEdgeIdsUpdaterDFS(currentNode);
            }
            if (childrenLeft>0) {
            } else {
                jPlaceEdgeMappingJPToNodeID.put(jplaceIdsCounter, currentNode.getId());
                jPlaceEdgeMappingNodeIDToJP.put(currentNode.getId(),jplaceIdsCounter);
                node.setJPlaceEdgeId(jplaceIdsCounter++);           
            }
        }
    }
    
    /**
     * used only if parsed tree is from jplace, i.e. edges are assigned to ids
     * with the {x} annotation on the right of branch length. map(jplaceEdgeId)=nodeId
     * @param jplaceEdgeId 
     * @return the nodeId holding the equivalent edge
     */
    public int getJplaceMappingJPToNodeID(int jplaceEdgeId) {
        return jPlaceEdgeMappingJPToNodeID.get(jplaceEdgeId);
    }
    
    /**
     * used only if parsed tree is from jplace, i.e. edges are assigned to ids
     * with the {x} annotation on the right of branch length. map(nodeId)=jplaceEdgeId
     * @param jplaceEdgeId 
     * @return the nodeId holding the equivalent edge
     */
    public int getJplaceMappingNodeIdToJP(int nodeId) {
        return jPlaceEdgeMappingNodeIDToJP.get(nodeId);
    }
    
    /**
     * return all jplace mapping in a map (see @getJplaceMapping)
     * @return 
     */
    public HashMap<Integer,Integer> getAllJPlaceMappingsJPToNodeID() {
        return jPlaceEdgeMappingJPToNodeID;
    }

    @Override
    public String toString() {
        StringBuilder sb=new StringBuilder();
        sb.append("#nodes="+this.getNodeCount());
        sb.append(" #leaves="+this.getLeavesCount());
        sb.append(" rooted="+this.isRooted());
        return sb.toString();
    }
    
    
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //This block of methods is used to map the internal lables of a tree
    //of similar topology than this tree (rooted or unrooted) 
    LinkedHashMap<Integer,Integer> nodeMapping=new LinkedHashMap<>();
    
    /**
     * map(this.tree nodeId)=otherTree nodeID;
     * not that compared PhyloTree need to be both rooted or both unrooted
     * @param otherTree
     * @return 
     */
    public LinkedHashMap<Integer,Integer> mapNodes(PhyloTree otherTree) {
        mapNodesByDFS(this.getRoot(),otherTree);
        //both rotted or both unrooted
        if ( (this.isRooted && otherTree.isRooted) || (!this.isRooted && !otherTree.isRooted)) {
//            System.out.println("nodeMapping: both trees are rooted");
//            //test if rooted on same edges ? 
//            PhyloNode root=this.getRoot();
//            Enumeration e=root.children();
//            boolean sameChildren=true;
//            while (e.hasMoreElements()) {
//                PhyloNode child=(PhyloNode)e.nextElement();
//                int idInOtherTree=nodeMapping.get(child.getId());
////                System.out.println("test "+child.getId()+" "+idInOtherTree);
////                System.out.println("A:"+otherTree.getById(idInOtherTree).getParent());
////                System.out.println("B:"+otherTree.getRoot());
//                if (otherTree.getById(idInOtherTree).getParent() != otherTree.getRoot()) {
//                    sameChildren=false;
//                }
//            }
//            //if yes roots can be mapped between both trees
//            if (sameChildren)
                nodeMapping.put(this.getRoot().getId(), otherTree.getRoot().getId());
            // if not, no mapping is done for the root, this is return as -1
//            else
//                nodeMapping.put(this.getRoot().getId(), -1);
        } else {
            System.out.println("Mapping of rooted to unrooted has no good solution... this case should not occur !");
            System.exit(1);
        }
        
        return nodeMapping;
    }
    
    /**
     * do a depth first transversal to go to all nodes, starts from a bait node 
     * @param node the bait node
     */
    private void mapNodesByDFS(PhyloNode node,PhyloTree otherTree) {
        
        //go down recursively in children
        Enumeration e=node.children();        
        while (e.hasMoreElements()) {
            PhyloNode n=(PhyloNode)e.nextElement();
            mapNodesByDFS(n,otherTree);
        }
        //returning up, from the sons
        //IF ROOT, no parent to associate
        if (node.isRoot()) {
            //System.out.println("SKIP ROOT");
            return;
        }
        //System.out.println("IN: "+node);
        //System.out.println("    parent to associate="+node.getParent());
        
//        //returns up, time to associate nodes, 3 possible case:
        //1. node is child of root, and just a leaf
        //we will associate only the leaf, not its parent
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
        
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    

    /**
     * return shortest path between nodes A and B (non null nodes)
     * @param root
     * @param a
     * @param b
     * @return 
     */
    public Path shortestPath(PhyloNode root, PhyloNode a, PhyloNode b) {

        assert null!=root;
        assert null!=a;
        assert null!=b;
        
        
        List<PhyloNode> l =new ArrayList<>();
                
        TreeNode[] pathAToRoot = a.getPath(); //[[10]added_root:0.000, [0]W:0.100, [2]B:0.200]
        TreeNode[] pathBToRoot = b.getPath(); //[[10]added_root:0.000, [0]W:0.100, [1]A:0.100]
        //LCA will be last common node on path
        PhyloNode LCA=null;
        int LCAIndex=-1;
        int shortestPathLength=pathAToRoot.length;
        if (pathBToRoot.length<shortestPathLength)
            shortestPathLength=pathBToRoot.length;

        for (int i = 0; i < shortestPathLength; i++) {
            if (pathAToRoot[i]!=pathBToRoot[i]) {  //[[10]added_root:0.000, [7]Z:0.300, [9]F:0.200]
                LCA=(PhyloNode)pathAToRoot[i-1];   //[[10]added_root:0.000, [0]W:0.100, [1]A:0.100]
                LCAIndex=i-1;
                //System.out.println("Diff at: "+i);
                break;
            }
        }
        //if we reach the end of the shortest path 
        if (LCAIndex==-1) {
            LCA=(PhyloNode)pathAToRoot[shortestPathLength-1];  //[[10]added_root:0.000, [7]Z:0.300, [9]F:0.200]
            LCAIndex=shortestPathLength-1;                     //[[10]added_root:0.000, [7]Z:0.300]
        }
        
        
//        System.out.println("A_to_root:"+Arrays.toString(pathAToRoot));        
//        System.out.println("B_to_root:"+Arrays.toString(pathBToRoot));
//        System.out.println("LCA      :"+LCA);
//        System.out.println("LCA_index:"+LCAIndex);
        
        //now merge both list, from LCAIndex to end of list, with inversion of A list
        float branchDist=0.0f;
        int nodeDist=0;
        boolean added_root=false;
        for (int i = pathAToRoot.length-1; i>LCAIndex; i--) {
            //System.out.println("i:"+i);
            //System.out.println("add rl "+pathAToRoot[i]);
            l.add((PhyloNode)pathAToRoot[i]);
            if (((PhyloNode)pathAToRoot[i]).getLabel().equals("added_root"))
                added_root=true;
            branchDist+=((PhyloNode)pathAToRoot[i]).getBranchLengthToAncestor();
            //System.out.println("branchDist:"+branchDist);
            if (i>LCAIndex && (i!=pathAToRoot.length-1)) {
                nodeDist++;
                //System.out.println("nodeDist ++  "+i);
            }
        }
        for (int i = LCAIndex; i<pathBToRoot.length; i++) {
            //System.out.println("add lr "+pathBToRoot[i]);
            l.add((PhyloNode)pathBToRoot[i]);
            if (((PhyloNode)pathBToRoot[i]).getLabel().equals("added_root"))
                added_root=true;
            if (i>LCAIndex) {
                branchDist+=((PhyloNode)pathBToRoot[i]).getBranchLengthToAncestor();
                //System.out.println("branchDist:"+branchDist);

            }
            if (i>LCAIndex && i!=pathBToRoot.length-1) { 
                nodeDist++;
                //System.out.println("nodeDist ++  "+i);
            }
        } 
        //if this path >1 (at least neighboors nodes)
        //then add +1 to node count, because LCA node was not counted above
        if (l.size()>1)
            nodeDist++;
        
//        System.out.println("path:"+l);
//        System.out.println("branchDist:"+branchDist);
//        System.out.println("nodeDist:"+nodeDist);
        
        Path p=new Path(l,nodeDist,branchDist,added_root);
        
        return p;
        
        
    }
    
    public class Path {
        public float branchDistance=-1.0f;
        public int nodeDistance=-1;
        public List<PhyloNode> path=null;
        public boolean withAddedRoot=false;

        public Path(List<PhyloNode> p, int nd, float bd, boolean withAddedRoot) {
            this.branchDistance=bd;
            this.nodeDistance=nd;
            this.path=p;
            this.withAddedRoot=withAddedRoot;
        }

        public boolean isWithAddedRoot() {
            return withAddedRoot;
        }
        
        @Override
        public String toString() {
            return "nd="+nodeDistance+" bd="+branchDistance+" "+path;
        }
 
    }
    
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    

}
