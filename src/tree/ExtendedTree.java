/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import etc.Infos;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.LinkedHashMap;
import javax.swing.tree.TreeModel;

/**
 * builds a tree from a PhyloTree with added fake nodes
 * @author ben
 */
public class ExtendedTree extends PhyloTree {
    
    private static final long serialVersionUID = 2100L;
    
    //operational variables
    public static final int BRANCHING_ON_NODE=1;
    public static final int BRANCHING_ON_BRANCH=2;
    public static final int BL_FROM_SUBTREE_MINMAX=100; //deprecated
    public static final int BL_FROM_SUBTREE_MEAN=101;
    
    //default values
    public static final float DEFAULT_BRANCHBREAK_LENGTH=-1.0f;
    public static final int DEFAULT_N=1;
    private int branchingMode=BRANCHING_ON_NODE;
    private int newBranchLengthMode=BL_FROM_SUBTREE_MEAN;
    private float branchbreakThreshold=DEFAULT_BRANCHBREAK_LENGTH;
    private int N=DEFAULT_N;
    
    //variables used only in the search of longest/shortest branch
    //or add new nodes
    private int fakeNodeCounter=0;
    private int level=0;
    private float l_DFS_cumul=0.0f;
    private float longest=0.0f;
    private float shortest=Float.MAX_VALUE; //init it to the highest possible float
    private int sum_B_leaves=0;
    private float l_sum_B_subtree=0.0f;
    
    //list of new terminal nodes (will be used to add gap-only fake sequences
    //to the extended alignment)
    private int lastOriginalId=-1;
    private ArrayList<PhyloNode> newLeaves=null;
    private ArrayList<PhyloNode> newInternalNodes=null;
    
    
    //associative maps describing the edges before and after tree extension
    /////////////////////////////////////////////////
    //map linking the new fake nodes to the original branches they come from.
    //for now, all fake nodes injected in a branch will be mapped to the 
    //son node of the original branch
    //note: in jplace output format, the integer used as branch ids (node{id})
    //are the nodeId of the PhyloTree object.
    //so here, we associated nodeId of (fake nodes injected in the branch)
    //to nodeId of (son of the original branch before fake nodes injection)
    //map(fakeNode Id)=nodeId of original node, son of the branch
    private HashMap<Integer,Integer> extendedNodesToOriginalNodes=null;
    // index=edgeId ; key=son,val=parent
    public LinkedHashMap<PhyloNode,PhyloNode> originalEdges=new LinkedHashMap<>();
    // key=edgeId ; key=son,vals=parent
    public LinkedHashMap<PhyloNode,PhyloNode> extendedEdges=new LinkedHashMap<>();
    
    
    
    
    /**
     * Build a extended tree with new N nodes on branches (BRANCHING_ON_BRANCH mode)
     * or one new node on nodes (not root)
     * @param tree
     * @param branchingMode one of BRANCHING_ON_NODE, BRANCHING_ON_BRANCH
     */
    public ExtendedTree(PhyloTree tree, int branchingMode ) {
        super(tree.getModel(),tree.isRooted());
        this.branchingMode=branchingMode;
        initRelaxedTree(tree,DEFAULT_BRANCHBREAK_LENGTH,DEFAULT_N);
    }

    /**
     * Build a extended tree with N new nodes on branches (BRANCHING_ON_BRANCH mode)
     * @param tree
     * @param branchbreackThreshold branch length below which no nodes are added on the branch
     * @param N number of fake nodes to add
     */
    public ExtendedTree(PhyloTree tree,float branchbreackThreshold,int N) {
        setModel(tree.getModel());
        this.branchingMode=BRANCHING_ON_BRANCH;
        this.N=N;
        this.isRooted=tree.isRooted;
        initRelaxedTree(tree, branchbreackThreshold,N);
        
    }
    
    /**
     * return all the leaves that were created in this extended tree.
     * @return 
     */
    public ArrayList<PhyloNode> getFakeLeaves() {
        return this.newLeaves;
    }
    
    /**
     * return all the internal nodes that were created in this extended tree.
     * @return 
     */
    public ArrayList<PhyloNode> getFakeInternalNodes() {
        return this.newInternalNodes;
    }    
    /**
     * for the given fake node nodeId, get the nodeId of the original node 
     * (son of the branch to which was originally injected the fake node),
     *  map(fakeNode Id)=nodeId of original node(son of the modified branch)
     * or the same id if this is an original node (same id iin original and 
     * extended tree).
     * @param nodeId
     * @return 
     */
    public Integer getFakeToOriginalId(int nodeId) {
        return extendedNodesToOriginalNodes.get(nodeId);
    }
    
    /**
     * for the given fake node nodeId, get the nodeId of the original node 
     * (son of the branch to which was originally injected the fake node),
     *  map(fakeNode Id)=nodeId of original node(son of the modified branch)
     * @return map of the mappings
     */
    public HashMap<Integer,Integer> getFakeNodeMapping() {
        return extendedNodesToOriginalNodes;
    }    
    /**
     * init the extended tree, used in constructors
     * @param tree
     * @param branchbreakThreshold 
     * @param 
     */
    private void initRelaxedTree(PhyloTree tree,float branchbreakThreshold,int N) {
        this.fakeNodeCounter=tree.getNodeCount();
        this.lastOriginalId=fakeNodeCounter;
        this.branchbreakThreshold=branchbreakThreshold;
        this.newLeaves=new ArrayList<>();
        this.newInternalNodes=new ArrayList<>();
        this.extendedNodesToOriginalNodes=new HashMap<>();
        Infos.println("# nodes in tree before extension: "+fakeNodeCounter);
        switch (branchingMode) {
            case BRANCHING_ON_NODE:
                populateWithFakeBranchToNodes_DFS((PhyloNode)tree.getModel().getRoot());
                break;
            case BRANCHING_ON_BRANCH:
                populateWithFakeBranchToEdges_DFS((PhyloNode)tree.getModel().getRoot());
                break;
            default:
                Infos.println("Cannot instanciate extended tree: unknown mode.");
                break;
        }
        Infos.println("# nodes in tree after extension: "+fakeNodeCounter);
    }

    /**
     * add a node X and its 2 sons FakeX1 and FakeX2 to each internal node 
     * of the original tree
     * (to the exception of the root)
     */
    @Deprecated
    private void populateWithFakeBranchToNodes_DFS(PhyloNode node) {
        // "<lastOriginalId" verify that node is nt one of the new fake nodes
        if (!node.isLeaf() && node.getId()<lastOriginalId) {
            fakeNodeCounter+=3;
            PhyloNode X0 =new PhyloNode(fakeNodeCounter-2, (fakeNodeCounter-2)+"_FakeX0", 0.0f);
            PhyloNode fakeX1 = new PhyloNode(fakeNodeCounter-1, (fakeNodeCounter-1)+"_FakeX1", 0.01f);
            PhyloNode fakeX2 = new PhyloNode(fakeNodeCounter, (fakeNodeCounter)+"_FakeX2", 0.01f);
            //add the 2 fakes to X
            X0.add(fakeX1);
            X0.add(fakeX2);
            newLeaves.add(fakeX1);
            newLeaves.add(fakeX2);
           
            float lengthNew=0.0f;
            if (newBranchLengthMode==BL_FROM_SUBTREE_MINMAX) {
                //calculate the shortest/longest branch length of X children subtrees
                //first init the shortest/longest with branch length of X to children
                l_DFS_cumul=0.0f;
                longest=0.0f;
                shortest=Float.MIN_VALUE;
                getBLFromMinMax_DFS(node,0);
                //System.out.println("    shortest after DFS search: "+shortest);
                //System.out.println("    longest after DFS search: "+longest);
                lengthNew=(longest-shortest)/2;
            } else if (newBranchLengthMode==BL_FROM_SUBTREE_MEAN){
                //mean build from cumulated path length divided by number
                //of encountered leaves
                l_DFS_cumul=0.0f;
                l_sum_B_subtree=0;
                sum_B_leaves=0;
                //System.out.println("### LAUNCH FROM "+node);
                getBLFromMean_DFS(node,0);
                lengthNew=(l_sum_B_subtree)/sum_B_leaves;
                //System.out.println("    --> lengthNew="+lengthNew);

            }
            X0.setBranchLengthToAncestor(lengthNew);
            
            //finally, add X to current node
            if (!node.isRoot())
                node.add(X0);
        }
        
        //go down recursively
        for (int i=0;i<node.getChildCount();i++) {
            PhyloNode child=(PhyloNode) node.getChildAt(i);
            if (!child.getLabel().contains("FAKE"))         //TO DO, change this to a node attribute
                populateWithFakeBranchToNodes_DFS(child);
        }
        
    }
        

    /**
     * add N nodes XN1 and their respective children XN2, FakeXN1 and FakeXN2
     * to each internal branch of the original tree
     * (to the exception of the root)
     */
    private void populateWithFakeBranchToEdges_DFS(PhyloNode node) {
        
//        System.out.println("### LAUNCH FROM "+node);
//        System.out.println("lastOriginalId:"+lastOriginalId);
//        System.out.println("node.getId()<lastOriginalId:"+(node.getId()<lastOriginalId));
//        System.out.println("!(node.getBranchLengthToAncestor()<=branchbreakThreshold):"+(!(node.getBranchLengthToAncestor()<=branchbreakThreshold)));
        

        //memorize current parent and children
        PhyloNode A=(PhyloNode)node.getParent();
        PhyloNode B=node;
//            System.out.println("   A:"+A);
//            System.out.println("   B:"+B);
        //register them in the node mapping

        //go down recursively, after memorizing which were the inital children
        Enumeration enumChildren = B.children();
        ArrayList<PhyloNode> initialChildren=new ArrayList<>();
        while (enumChildren.hasMoreElements()) {
            PhyloNode Bi= (PhyloNode) enumChildren.nextElement();
            initialChildren.add(Bi);  
        }
//        System.out.println("### LIST of B children: "+initialChildren);
        for (int i = 0; i < initialChildren.size(); i++) {
            PhyloNode Bi = initialChildren.get(i);
//            System.out.println("### DFS LAUNCHED on "+Bi);
            populateWithFakeBranchToEdges_DFS(Bi);
//            System.out.println("### DFS RETURNS from "+Bi);  
        }


        //before conditions below (branchbreakThreshold), register edges in original map
        if (!B.isRoot())
            originalEdges.put(B, A);

        if ( (!B.isRoot()) &&
             (((PhyloNode)B.getParent()).getId()<lastOriginalId) && 
             (B.getId()<lastOriginalId) &&
             (!(B.getBranchLengthToAncestor()<branchbreakThreshold)) ) {

            //define l_init and l_b for the current edge
            float l_init=B.getBranchLengthToAncestor();
            float l_b=(0.0f+l_init)/(N+1);
//                System.out.println("  l_init:"+l_init+" l_b:"+l_b);

            //cut parent A from children B
            A.remove(B);
//            System.out.println("B:"+B+" cut from A:"+A);
            //add tehm in node mappings
            extendedNodesToOriginalNodes.put(A.getId(), A.getId());
            extendedNodesToOriginalNodes.put(B.getId(), B.getId());

            //build and attach N fake nodes to parent and subsequent Xi
            PhyloNode currentParent=A;
            for (int j = 0; j < N; j++) {
//                System.out.println("currentParent:"+currentParent);
                //built subtree of the jth X0
                fakeNodeCounter+=4;
                PhyloNode X0 =new PhyloNode(fakeNodeCounter-3, (fakeNodeCounter-3)+"_X0", 0.01f);
                PhyloNode X1 =new PhyloNode(fakeNodeCounter-2, (fakeNodeCounter-2)+"_X1", 0.01f);
                PhyloNode X2 = new PhyloNode(fakeNodeCounter-1, (fakeNodeCounter-1)+"_X2", 0.01f);
                PhyloNode X3 = new PhyloNode(fakeNodeCounter, fakeNodeCounter+"_X3", 0.01f);
                X1.add(X2);
                X1.add(X3);
                X0.add(X1);
                newLeaves.add(X2);
                newLeaves.add(X3);
                newInternalNodes.add(X0);
                newInternalNodes.add(X1);
                //register the internal nodes assoxiation in the map
                extendedNodesToOriginalNodes.put(X0.getId(), B.getId());
                extendedNodesToOriginalNodes.put(X1.getId(), B.getId());
                extendedEdges.put(X2, X1);
                extendedEdges.put(X3, X1);
                extendedEdges.put(X1, X0);

                //define length left from this X0 to B
                float l_XO_B=l_init-l_b*(j+1);
//                System.out.println("l_XO_B:"+l_XO_B);
                //define X0-X1 branch length
                float l_new=0.0f;
                if (newBranchLengthMode==BL_FROM_SUBTREE_MINMAX) {
                    //calculate the shortest/longest branch length of X children subtrees
                    //first init the shortest/longest with branch length of X to children
                    l_DFS_cumul=0.0f;
                    longest=0.0f;
                    shortest=Float.MIN_VALUE;
                    getBLFromMinMax_DFS(B,0);
                    l_new=((longest+N*l_b)-(shortest+N*l_b))/2;
                } else if (newBranchLengthMode==BL_FROM_SUBTREE_MEAN){
                    //mean build from cumulated path length divided by number
                    //of encountered leaves
                    l_DFS_cumul=0.0f;
                    l_sum_B_subtree=0;
                    sum_B_leaves=0;
                    //if this branch is from internal node to leaf
                    //we will not define the XO-X1 length as the mean branch
                    //length of the subtree, but as l_b.
                    if (!B.isLeaf()) {
                        getBLFromMean_DFS(B,0);
                        l_new=(sum_B_leaves*l_XO_B+l_sum_B_subtree)/sum_B_leaves;
                    } else {
                        sum_B_leaves=1;
                        l_new=l_b;
                    }

                }
//                System.out.println("l_new:"+l_new);
                //set branch lengths of X0 and X1
                X1.setBranchLengthToAncestor(l_new);
                X0.setBranchLengthToAncestor(l_b);
                
                //set branch length to nodes of original tree 
                X1.setBranchLengthToOriginalAncestor((j+1)*l_b+l_new);
                X1.setBranchLengthToOriginalSon(l_XO_B+l_new);
                X0.setBranchLengthToOriginalAncestor((j+1)*l_b);
                X0.setBranchLengthToOriginalSon(l_XO_B);
//                System.out.println("X1 blta:"+((j+1)*l_b+l_new));
//                System.out.println("X1 blts:"+(l_XO_B+l_new));
//                System.out.println("X0 blta:"+((j+1)*l_b));
//                System.out.println("X0 blts:"+l_XO_B);
                
                //attach X0 to parent
                currentParent.add(X0);
                extendedEdges.put(X0, currentParent);
                //and will move to next X0
                currentParent=X0;
            }

            //when all X0 are added, attach the last X0 to B
            currentParent.add(B);
            B.setBranchLengthToAncestor(l_init-l_b*N);  

        }

        //System.out.println("### AFTER NEW BRANCHING "+node);

    }

    
    /**
     * calculate a mean branch length from the subtree of a given node
     * @param node 
     */
    private void getBLFromMean_DFS(PhyloNode node, int level) {
        if(node.isLeaf()) {
            //System.out.println("          LEAF:"+node);
            l_sum_B_subtree+=l_DFS_cumul+node.getBranchLengthToAncestor();
            sum_B_leaves+=1;
            //System.out.println("          meanSum: "+meanSum);
            //System.out.println("          meanLeaves: "+meanLeaves);
        } else {
            if (level>0) {
                l_DFS_cumul+=node.getBranchLengthToAncestor();
                //System.out.println("       cumul+:"+cumul);
            } else {
                l_DFS_cumul=0;
                l_sum_B_subtree=0;
                //System.out.println("       cumul to 0");
            }
            //go down recursively
            for (int i=0;i<node.getChildCount();i++) {
                PhyloNode child=(PhyloNode) node.getChildAt(i);
                getBLFromMean_DFS(child,(level+1));
            }
            if (level>0) {
                l_DFS_cumul-=node.getBranchLengthToAncestor();
                //System.out.println("       cumul-:"+cumul);
            }
        }
    }
    
    
    
    
    /**
     * retrieve longest and shortest path to leaves from a given node
     * @param node
     * @param level 
     */
    @Deprecated
    private void getBLFromMinMax_DFS(PhyloNode node, int level) {
        //System.out.println("       DFS_shortlong: "+node);
        if(node.isLeaf()) {
            //System.out.println("          LEAF:"+node);
            float pathLength=l_DFS_cumul+node.getBranchLengthToAncestor();
            //System.out.println("          pathLength: "+pathLength);
            if ( pathLength > longest ) {
                longest=pathLength;
                //System.out.println("          UPDATE longest:"+pathLength);
            }
            if (pathLength < shortest ) {
                shortest=pathLength;
                //System.out.println("          UPDATE shortest:"+pathLength);
            }
        } else {
            if (level>0) {
                l_DFS_cumul+=node.getBranchLengthToAncestor();
                //System.out.println("       cumul+:"+cumul);
            }
            //go down recursively
            for (int i=0;i<node.getChildCount();i++) {
                PhyloNode child=(PhyloNode) node.getChildAt(i);
                getBLFromMinMax_DFS(child,(level+1));
            }
            if (level>0) {
                l_DFS_cumul-=node.getBranchLengthToAncestor();
                //System.out.println("       cumul-:"+cumul);
            }
        }
    }
        
    
    
}
