/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import alignement.Alignment;
import etc.Infos;
import inputs.Fasta;
import java.awt.Dimension;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;
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
    

    public PhyloTree(TreeModel newModel, boolean rooted) {
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
        
        //for RMQ-LCA algo, we init the corresponding tables
        // nodeCount is the highest value of node in our tree
        euler = new int[2 * nodeCount - 1]; // for euler tour sequence
        level = new int[2 * nodeCount - 1]; // level of nodes in tour sequence
        f_occur= new int[2 * nodeCount - 1]; // to store 1st occurance of nodes
        sc = new St_class();        
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
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
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
        mapNodesByDFS(this.getRoot(),otherTree);
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
    private void mapNodesByDFS(PhyloNode node,PhyloTree otherTree) {
        
        //go down recursively is some children
        Enumeration e=node.children();        
        while (e.hasMoreElements()) {
            PhyloNode n=(PhyloNode)e.nextElement();
            mapNodesByDFS(n,otherTree);
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
    

    /**
     * return a copy of this tree as a completely independant deep copy;
     * the tree is copied node per node through a preorder traversal,
     * not keeping any reference to the current tree.
     * This is useful for doing tree backups before pruning
     * experiments.
     * @param tree
     * @return 
     */
    public PhyloTree copyPhyloTree(PhyloTree tree) {
        
        PhyloNode root = tree.getRoot();
        
        
        return null;
    }
    
    private PhyloNode copyNode(PhyloNode originalNode) {
        if (originalNode==null) {
            return null;
        }
        //copy node, using constructor copy 
        PhyloNode copiedNode=new PhyloNode();
        Enumeration<PhyloNode> children = originalNode.children();
        while (children.hasMoreElements()) {
            PhyloNode nextElement = children.nextElement();
            copiedNode.add(copyNode(nextElement));
        }
        return copiedNode;
        
    }
    
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //This block of methods is used to furnish effecient calculation of
    //shortest path between 2 nodes (better than the bottom-up apporach 
    //of the standard DefaultMutableNode.getSharedAncestor(node)
    //this uses RMQ+LCA: as described in geeksforgeeks.org/find-lca-in-binary-tree-using-rmq/
    // it is O(n) for preprocessing, then O(log n) for finding LCA of 2 nodes
    
    private class St_class {
        int st;
        int stt[] = new int[10000];
    }
    
    
    //these tables are init to a specific size by initIndexes()
    int euler[] = null; // for euler tour sequence
    int level[] = null; // level of nodes in tour sequence
    int f_occur[] = null; // to store 1st occurance of nodes
    int fill; // variable to fill euler and level arrays
    St_class sc = new St_class();
  
    // log base 2 of x
    private int Log2(int x) {
        int ans = 0;
        int y = x >>= 1;
        while (y-- != 0)
            ans++;
        return ans;
    }
  
    private int swap(int a, int b) {
        return a;
    }
  
    /*  A recursive function to get the minimum value in a given range
     of array indexes. The following are parameters for this function.
   
     st    --> Pointer to segment tree
     index --> Index of current node in the segment tree. Initially
     0 is passed as root is always at index 0
     ss & se  --> Starting and ending indexes of the segment represented
     by current node, i.e., st[index]
     qs & qe  --> Starting and ending indexes of query range */
    private int RMQUtil(int index, int ss, int se, int qs, int qe, St_class st) {
        // If segment of this node is a part of given range, then return
        //  the min of the segment
        if (qs <= ss && qe >= se)
            return st.stt[index];
  
        // If segment of this node is outside the given range
        else if (se < qs || ss > qe)
            return -1;
  
        // If a part of this segment overlaps with the given range
        int mid = (ss + se) / 2;
  
        int q1 = RMQUtil(2 * index + 1, ss, mid, qs, qe, st);
        int q2 = RMQUtil(2 * index + 2, mid + 1, se, qs, qe, st);
  
        if (q1 == -1)
            return q2;
        else if (q2 == -1)
            return q1;
  
        return (level[q1] < level[q2]) ? q1 : q2;
    }
  
    // Return minimum of elements in range from index qs (quey start) to
    // qe (query end).  It mainly uses RMQUtil()
    private int RMQ(St_class st, int n, int qs, int qe) {
        // Check for erroneous input values
        if (qs < 0 || qe > n - 1 || qs > qe) 
        {
            System.out.println("Invalid input");
            return -1;
        }
  
        return RMQUtil(0, 0, n - 1, qs, qe, st);
    }
  
    // A recursive function that constructs Segment Tree for array[ss..se].
    // si is index of current node in segment tree st
    private void constructSTUtil(int si, int ss, int se, int arr[], St_class st) {
        // If there is one element in array, store it in current node of
        // segment tree and return
        if (ss == se)
            st.stt[si] = ss;
        else
        {
            // If there are more than one elements, then recur for left and
            // right subtrees and store the minimum of two values in this node
            int mid = (ss + se) / 2;
            constructSTUtil(si * 2 + 1, ss, mid, arr, st);
            constructSTUtil(si * 2 + 2, mid + 1, se, arr, st);
  
            if (arr[st.stt[2 * si + 1]] < arr[st.stt[2 * si + 2]])
                st.stt[si] = st.stt[2 * si + 1];
            else
                st.stt[si] = st.stt[2 * si + 2];
        }
    }
  
    /* Function to construct segment tree from given array. This function
     allocates memory for segment tree and calls constructSTUtil() to
     fill the allocated memory */
    private int constructST(int arr[], int n) {
        // Allocate memory for segment tree
        // Height of segment tree
        int x = Log2(n) + 1;
          
        // Maximum size of segment tree
        int max_size = 2 * (1 << x) - 1;  //  2*pow(2,x) -1
  
        sc.stt = new int[max_size];
  
        // Fill the allocated memory st
        constructSTUtil(0, 0, n - 1, arr, sc);
          
        // Return the constructed segment tree
        return sc.st;
    }
  
    // Recursive version of the Euler tour of T
    private void eulerTour(PhyloNode node, int l) 
    {
        euler[fill] = node.getId()+1; // insert in euler array  //ben:node id shift by 1 to be in [1;nodeCount]
        level[fill] = l;         // insert l in level array
        fill++;                  // increment index

        /* if unvisited, mark first occurrence */
        if (f_occur[node.getId()+1] == -1)            //ben:node id shift by 1 to be in [1;nodeCount]
            f_occur[node.getId()+1] = fill - 1;       //ben:node if shift by 1 to be in [1;nodeCount]
        //ben: below left/right recursion modified to handle more that 2 sons
        //i.e. for unrooted tree, where the Jtree "root" has 3 sons
//            /* tour left subtree if exists, and remark euler
//               and level arrays for parent on return */
//            if (node.getChildAt(0) != null) 
//            {
//                eulerTour(node.left, l + 1);
//                euler[fill] = node.getId()+1;  //ben:node id shift by 1 to be in [1;nodeCount]
//                level[fill] = l;
//                fill++;
//            }
//            /* tour right subtree if exists, and remark euler
//               and level arrays for parent on return */
//            if (node.getChildAt(1) != null) 
//            {
//                eulerTour(node.right, l + 1);
//                euler[fill] = node.getId()+1;  //ben:node id shift by 1 to be in [1;nodeCount]
//                level[fill] = l;
//                fill++;
//            }
        Enumeration<PhyloNode> children = node.children();
        while (children.hasMoreElements()) {
            PhyloNode nextChild = children.nextElement();
            eulerTour(nextChild, l + 1);
            euler[fill] = node.getId()+1;  //ben:node id shift by 1 to be in [1;nodeCount]
            level[fill] = l;
            fill++;
        }

            
        
    }
  
    /**
     * returns LCA of node u and v assuming they are present in tree; note that
     * this will be launch only after the tree get more than 2 leaves
     * (3 nodes) and that initIndexes() was called.
     * @param node
     * @param a node a
     * @param b node b
     * @return 
     */
    public PhyloNode findLCA(PhyloNode a, PhyloNode b) 
    {
        //note that this function can be executed only if 
        //the tree was indexed with initIndexes()
        assert nodeCount>=3;
        
        //ben: note that this algo considers node ids between 1 and nodeCount
        //as our program uses ids between 0 and nodeCount-1, we shift these ids
        //we unshift the LCA node id when returned
        int u=a.getId()+1;
        int v=b.getId()+1;
        
        assert u>0 && u<=nodeCount;
        assert v>0 && v<=nodeCount;
        
        /* Mark all nodes unvisited.  Note that the size of
           firstOccurrence is 1 as node values which vary from
           1 to 9 are used as indexes */
        Arrays.fill(f_occur, -1);
  
        /* To start filling euler and level arrays from index 0 */
        fill = 0;
  
        /* Start Euler tour with root node on level 0 */
        eulerTour(this.getRoot(), 0);
        
        System.out.println("Euler tour : "+Arrays.toString(euler));
        System.out.println("Euler level: "+Arrays.toString(level));
        System.out.println("First occur: "+Arrays.toString(f_occur));
         
        /* construct segment tree on level array */
        sc.st = constructST(level, 2 * v - 1);
        
        System.out.println(sc.st);
        System.out.println(Arrays.toString(sc.stt));
          
        /* If v before u in Euler tour.  For RMQ to work, first
         parameter 'u' must be smaller than second 'v' */
        System.out.println("u="+u+" v="+v);
        if (f_occur[u] > f_occur[v]) {
            //we swap
            int t=u;
            u=v;
            v=t;
            System.out.println("swapped! u="+u+" v="+v);
        }
  
        // Starting and ending indexes of query range
        int qs = f_occur[u];
        int qe = f_occur[v];
        
        System.out.println("qs="+qs+" qe="+qe);
  
        // query for index of LCA in tour
        int index = RMQ(sc, 2 * nodeCount - 1, qs, qe);
  
        /* return LCA node */
        return this.getById(euler[index]);
  
    }
    
    //--------------------------------------------------------------------------


    /**
     * return shortest path between non null nodes
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
        
        
        //System.out.println("A_to_root:"+Arrays.toString(pathAToRoot));        
        //System.out.println("B_to_root:"+Arrays.toString(pathBToRoot));
        //System.out.println("LCA      :"+LCA);
        //System.out.println("LCA_index:"+LCAIndex);
        
        //now merge both list, from LCAIndex to end of list, with inversion of A list
        float branchDist=0.0f;
        int nodeDist=0;
        for (int i = pathAToRoot.length-1; i>LCAIndex; i--) {
            //System.out.println("i:"+i);
            //System.out.println("add rl "+pathAToRoot[i]);
            l.add((PhyloNode)pathAToRoot[i]);
            branchDist+=((PhyloNode)pathAToRoot[i]).getBranchLengthToAncestor();
            if (i>LCAIndex && (i!=pathAToRoot.length-1)) {
                nodeDist++;
                //System.out.println("nodeDist ++  "+i);
            }
        }
        for (int i = LCAIndex; i<pathBToRoot.length; i++) {
            //System.out.println("add lr "+pathBToRoot[i]);
            l.add((PhyloNode)pathBToRoot[i]);
            if (i>LCAIndex)
                branchDist+=((PhyloNode)pathBToRoot[i]).getBranchLengthToAncestor();
            if (i>LCAIndex && i!=pathBToRoot.length-1) { 
                nodeDist++;
                //System.out.println("nodeDist ++  "+i);
            }
        } 
        //if this path >2 (not neighboors are same nodes)
        //then add +1 to node count, because LCA node was not counted above
        if (l.size()>2)
            nodeDist++;
        
        //System.out.println("path:"+l);
        //System.out.println("branchDist:"+branchDist);
        //System.out.println("nodeDist:"+nodeDist);
        
        Path p=new Path(l,nodeDist,branchDist);
        
        return p;
        
        
    }
    
    public class Path {
        public float branchDistance=-1.0f;
        public int nodeDistance=-1;
        public List<PhyloNode> path=null;

        public Path(List<PhyloNode> p, int nd, float bd) {
            this.branchDistance=bd;
            this.nodeDistance=nd;
            this.path=p;
        }  
    }
    
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    
    
    public static void main(String[] args) {
        //test that RMQ-LCA approach was correctly adapted
        
        //load a test tree
        String treeString="(A:0.1,B:0.2,((C:0.1,D:0.2)Y:0.1,(E:0.1,F:0.2)Z:0.2)X:0.2)W:0.0;";
        PhyloTree tree = NewickReader.parseNewickTree2(treeString, false);
        tree.initIndexes();
        
        for (int i = 0; i < tree.getNodeIdsByDFS().size(); i++) {
            System.out.println(tree.getById(tree.getNodeIdsByDFS().get(i)));
        }
        
        
        System.out.println("A^B:"+tree.findLCA(tree.getByName("A"), tree.getByName("B")));
        

        
        System.out.println("B^Y:"+tree.findLCA(tree.getByName("B"), tree.getByName("Y")));
        
    }
    

}
