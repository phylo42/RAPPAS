/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import etc.Infos;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Stack;
import javax.swing.JFrame;

/**
 * 
 * @author ben
 */
public class NewickReader {
    
    /**
     * build a phylotree from a newick string, and init its indexes
     * @param s
     * @param forceRooting if the newick describes an unrooted tree (3 sons at
     * top level), then force the return of a rooted tree; not that the root
     * will always be placed as follows: 
     * (son1,son2,son3)newick_root; -->  ((son1,son2)newick_root,son3)added_root;
     * @param considerJplaceEdgeIds consider the {x} edge ids when this
     * newick tree is a jplace style tree. If not activated the presence of this 
     * labelling will explicitely throw an error. Example:
     * NumberFormatException: For input string: "0.00000100000050002909{0}" ;
     * if activated, before parsing branch length, a test is done and a map 
     * of jplace edgeIds to PhyloTree nodeIds is build (edge assigned to son node) ;
     * this mapping can be called with PhyloTree.getJPlaceMapping() 
     * @return 
     */
    
    public static PhyloTree parseNewickTree2(String s, boolean forceRooting, boolean considerJplaceEdgeIds) {
        
        if (s==null) {
            System.out.println("Cannot read tree, string is null");
            return null;
        }
        
        
        
        //counter to build internal nodeIds
        //(different from labels read in the newick)
        int currentNodeIndex=-1;
        //tree depth
        int depth=0;
        //pile of parent nodes
        Stack<PhyloNode> stackedParents=new Stack<>();
        HashMap<Integer,ArrayList<PhyloNode>> nodesPerDepth=new HashMap<>();
        nodesPerDepth.put(depth, new ArrayList<>()); //init highest depth (=0)
        
        StringBuilder sb=new StringBuilder();
        boolean descending=false;
        boolean ascending=false;
        PhyloNode bufferedNode=null;
        //the newick is red from left to rigth, parenthesis are used to 
        //jump dig in and out at different depths
        for (int i = 0; i < s.length(); i++) {
//            System.out.println("i:"+i);
//            System.out.println("--------------------------------------------");
//            System.out.println("-BEFORE-------------------------------------");
//            System.out.println("char:"+s.charAt(i));
//            System.out.println("depth:"+depth);
//            System.out.println("stackedParents:"+stackedParents);
//            System.out.println("nodesPerDepth:"+nodesPerDepth);
//            System.out.println("sb:"+sb);
//            System.out.println("currentNodeIndex:"+currentNodeIndex);
//            System.out.println("bufferedNode:"+bufferedNode);
//            System.out.println("descending:"+descending);
//            System.out.println("ascending:"+ascending);
            //init new node when (
            if (s.charAt(i)=='(') {
                stackedParents.push(new PhyloNode(++currentNodeIndex));
                depth++;
                nodesPerDepth.put(depth, new ArrayList<>());
                sb.delete(0, sb.length());
                descending=true;
                ascending=false;
            //fill init but empty node when )
            } else if (s.charAt(i)==')') {
                //make node with previous node (if exists), at depth before )
                if (sb.length()>0) {
                    //take potential content of sb
                    //to fill node
                    String[] data=sb.toString().split(":");
                    String label=data[0];
                    Float bl=0.0f;
                    int jplaceEdgeId=-1;
                    if (data.length>1) {
                        //if jplace file, consider {x} as jplace nodeids
                        if (considerJplaceEdgeIds) {
                            int openBracketIndex=data[1].indexOf('{');
                            jplaceEdgeId=Integer.valueOf(data[1].substring(openBracketIndex+1, data[1].indexOf('}')));
                            bl=Float.parseFloat(data[1].substring(0, openBracketIndex));
                        } else {
                            bl=Float.parseFloat(data[1]);
                        }
                    }
                    if (!ascending)
                        nodesPerDepth.get(depth).add(new PhyloNode(++currentNodeIndex, label, bl, jplaceEdgeId, false));
                    else {
                        bufferedNode.setLabel(label);
                        bufferedNode.setBranchLengthToAncestor(bl);
                        bufferedNode.setJPlaceEdgeId(jplaceEdgeId);
                    }
                }
                sb.delete(0, sb.length());
                
                //get parent which is on top of stack
                bufferedNode=stackedParents.pop();
                //associate nodes to this parent
                for (PhyloNode son:nodesPerDepth.get(depth)) {
                    bufferedNode.add(son);
                    //System.out.println("DO: parent="+bufferedNode+"  son="+son);
                }
                nodesPerDepth.get(depth).clear();
                
                //now moving to depth after )
                //not forgetting to register the paretn to this new depth
                depth--;
                descending=false;
                ascending=true;
                nodesPerDepth.get(depth).add(bufferedNode);

            //add to current depth when ,    
            } else if (s.charAt(i)==',') {
                //take potential content of sb
                //to fill node
                String[] data=sb.toString().split(":");
                String label=data[0];
                Float bl=0.0f;
                int jplaceEdgeId=-1;
                if (data.length>1) {
                    //if jplace file, consider {x} as jplace nodeids
                    if (considerJplaceEdgeIds) {
                        int openBracketIndex=data[1].indexOf('{');
                        jplaceEdgeId=Integer.valueOf(data[1].substring(openBracketIndex+1, data[1].indexOf('}')));
                        bl=Float.parseFloat(data[1].substring(0, openBracketIndex));
                    } else {
                        bl=Float.parseFloat(data[1]);
                    }
                }
                //if comes from upper depth, start list of nodes for this depth
                if (descending) {
                    nodesPerDepth.get(depth).add(new PhyloNode(++currentNodeIndex, label, bl, jplaceEdgeId, false));
                //if ascending, come up from subtree, we use the unstacked node
                } else if (ascending) {
                    bufferedNode.setLabel(label);
                    bufferedNode.setBranchLengthToAncestor(bl);
                    bufferedNode.setJPlaceEdgeId(jplaceEdgeId);
                    //levelList.get(depth).add(bufferedNode);
                //neither ascending or descending (case happening when unrooted
                //tree, 1st level appears with 3 nodes, i.e. ((),(),());
                } else {
                    nodesPerDepth.get(depth).add(new PhyloNode(++currentNodeIndex, label, bl, jplaceEdgeId, false));
                }
                sb.delete(0, sb.length());
                 
                descending=false;
                ascending=false;
                 
            //end of newick, take last node as root
            } else if (s.charAt(i)==';') {
                //take potential content of sb
                //to fill node
                String[] data=sb.toString().split(":");
                String label=data[0];
                Float bl=-1.0f; 
                int jplaceEdgeId=-1;
                if (data.length>1) {
                    //if jplace file, consider {x} as jplace nodeids
                    if (considerJplaceEdgeIds) {
                        int openBracketIndex=data[1].indexOf('{');
                        jplaceEdgeId=Integer.valueOf(data[1].substring(openBracketIndex+1, data[1].indexOf('}')));
                        bl=Float.parseFloat(data[1].substring(0, openBracketIndex));
                    } else {
                        bl=Float.parseFloat(data[1]);
                    }
                }
                
                bufferedNode.setJPlaceEdgeId(jplaceEdgeId);
                bufferedNode.setLabel(label);
                bufferedNode.setBranchLengthToAncestor(bl);
            //else, simply fill the stringbuilder, 
            //which content will be used when 
            //new nodes needs to be filled with
            //node name and branch length
            } else {
                sb.append(s.charAt(i));
            }
            
            
//            System.out.println("-AFTER----");
//            System.out.println("char:"+s.charAt(i));
//            System.out.println("depth:"+depth);
//            System.out.println("stack:"+stackedParents);
//            System.out.println("nodesPerDepth:"+nodesPerDepth);
//            System.out.println("sb:"+sb);
//            System.out.println("currentNodeIndex:"+currentNodeIndex);
//            System.out.println("bufferedNode:"+bufferedNode);
//            System.out.println("descending:"+descending);
//            System.out.println("ascending:"+ascending);
        }
        
        //if the root contain 2 sons, it's rooted, if three it's unrooted
        boolean rooted=false;
        int rootSons=0;
        for (Enumeration e=bufferedNode.children();e.hasMoreElements();) {
            PhyloNode pn= (PhyloNode)e.nextElement();
            if (e!=null)
                rootSons++;
        }
        //rooted
        if (rootSons<3) {
            rooted=true;
        }
        
        
        //root this unrooted tree if asked by the user
        PhyloTree tree=null;
        if (!rooted && forceRooting) {
            Infos.println("Rooting of input unrooted Tree !");
            //rooting will be done on the edge linking the newick root 
            //and the 3 son: 
            //(son1,son2,son3)newick_root; -->  ((son1,son2)newick_root,son3)added_root;
            //
            //   newick_root                      added_root
            //     / | \bl=1.5  ==>          bl=0/   \ bl=1.5
            //    /  |  \               newick_root   \
            // son1 son2 son3                /  \      son3
            //                            son1  son2
            //
            PhyloNode newick_root=bufferedNode;
            //PhyloNode son1=bufferedNode.getChildAt(0);
            //PhyloNode son2=bufferedNode.getChildAt(1);
            PhyloNode son3=bufferedNode.getChildAt(2);
            PhyloNode added_root=new PhyloNode(++currentNodeIndex, "added_root", 0.0f, -1, false);
            //unlink sons3
            float son3_bl=son3.getBranchLengthToAncestor();
            son3.removeFromParent();
            //set new branch lengths
            son3.setBranchLengthToAncestor(son3_bl);
            newick_root.setBranchLengthToAncestor(0.0f);
            //link son3 and newick_root to added_root
            added_root.add(newick_root);
            added_root.add(son3);
            //build tree
            tree=new PhyloTree(new PhyloTreeModel(added_root),true, considerJplaceEdgeIds);
        } else {
            //last bufferedNode is the root, i.e. ([sons])bufferedNode; in the newick
            tree=new PhyloTree(new PhyloTreeModel(bufferedNode),rooted, considerJplaceEdgeIds);
        }
        
        
        
        
        //init indexes related to internal/leaves stats
        tree.initIndexes();
        
        
        return tree;
        
    }
    
    

    
    
    
    
    /**
     * test main
     * @param args 
     */
    public static void main(String[] args) throws InterruptedException {
        
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
        
        String t_basic_unrooted="((L1:1,L2:2)X:4,(L3:2,L4:2)Y:1,L5:4)Z:4;"; //unrooted
        String t_basic_rZX="((L1:1,L2:2)X:4,((L3:2,L4:2)Y:1,L5:4)Z:4)root:0;";//rooted on Z-X
        String t_basic_rZX_inv="(((L3:2,L4:2)Y:1,L5:4)Z:4,(L1:1,L2:2)X:4)root:0;";//rooted on Z-X, inversed children of root
        String t_basic_rZL5="(L5:4,((L3:2,L4:2)Y:1,(L1:1,L2:2)X:4)Z:4)root:0;"; //rooted on Z-L5


        
        String t_basic2="((L1:2,L2:2)I:2,L3:4)root:0;";
        String t_unrooted="(L3:4,L4:3,(L1:2,L2:2)I:2);";
        String t_unrooted2="(L3:4,(L1:2,L2:2)I:2,L4:3);";
        String t_unrooted3="((L1:2,L2:2)I:2,L3:4,L4:3);";
        
        System.out.println("START");
        System.out.println(t_basic_rZL5);
        System.out.println(t_basic_rZX);
        
        PhyloTree tree1 = NewickReader.parseNewickTree2(t_basic_rZX_inv, false, false);
        System.out.println("t_basic parsed!");
        tree1.displayTree();
        System.out.println("isRooted:"+tree1.isRooted());
        System.out.println("Struct root:"+tree1.getRoot());
        
        System.out.println("START");
        PhyloTree tree2 = NewickReader.parseNewickTree2(t_basic_rZX, false, false);
        System.out.println("t_basic_unrooted parsed!");
        tree2.displayTree();
        System.out.println("isRooted:"+tree2.isRooted());
        System.out.println("Struct root:"+tree2.getRoot());
        
        
        //node mapping test
        System.out.println(tree1.mapNodes(tree2));
        
        
        
        
        Thread.sleep(60000);
        
        System.exit(1);
        

        
        //test parsing
        long startTime = System.currentTimeMillis();
        PhyloTree t=new NewickReader().parseNewickTree2(treeFASTML, false, false);
        long endTime = System.currentTimeMillis();
        System.out.println("Parsing took " + (endTime - startTime) + " milliseconds");
        
        
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
