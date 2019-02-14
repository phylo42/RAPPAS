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

/**
 * 
 * @author ben
 */
public class NewickReader {
    
    /**
     * build a phylotree from a newick string, and init its indexes
     * @param s
     * @param forceRooting if the newick describes an unrooted tree (3 sons at
     * top level), then forces the return of a rooted tree; note that the root
     * will always be inserted as follows: 
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
        
        try {

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
            //this is where the jPlaceEdgeMappings of PhyloTree are updated
            tree.initIndexes();


            return tree;
            
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("Something went wrong during tree parsing.");
            System.out.println("Please check that you input tree is in newick format.");
            System.out.println("Please avoid (),:; characters in your labels.");
            System.out.println("If tree contains {x}, i.e. jplace edge labels, configure parser accordingly.");
            System.exit(1);
            return null;
        }
        
    }
    
}
