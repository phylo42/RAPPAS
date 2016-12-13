/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import etc.Infos;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.tree.DefaultTreeModel;

/**
 * 
 * @author ben
 */
public class NewickParser {
    
    private int currentElementIndex=1;
    //get root info name:length;
    private Pattern rootPattern= Pattern.compile("\\)([^\\)]*);$"); //*? to make it non greedy    
    //internal nodes infos are on the right of a subtree (subtree)name:length
    private Pattern internalNodePattern= Pattern.compile("\\)([^\\)]*?):([0-9\\.Ee-]*?)$|\\)([^\\)]*?)$"); //*? to make it non greedy    
    //simple name:length
    private Pattern leafPattern =Pattern.compile("(.*?):([0-9\\.Ee-]*?)$|(.*)$");
    
    private int level=0;
    
    /**
     * simple newick parser, assumes that the root is the first element of the newick sting representation
     * @param newickTree
     * @return 
     */
    public PhyloTree parseNewickTree (String newickTree) {
        currentElementIndex=1;
        level=0;
//        long startTime = System.currentTimeMillis();
        //we will initialize the root immediately during the newick parsing
        PhyloNode root=null;
        
        //get root info
        Matcher matcher = rootPattern.matcher(newickTree);

        int rootStringLen=0;
        String name="";
        String len="";
        String[] infos=null;
        if (matcher.find()) {
            infos=matcher.group(1).split(":");
            name=infos[0];
            if (infos.length>1) {len=infos[1];}
        }

        //the last element, with ';' is the root infos (label, branch length)
        if (name.equals("") && len.equals("")) {
            root=new PhyloNode(1, "root",0.0);
            rootStringLen=1; //only ;
        } else if (len.equals("")) {
            root=new PhyloNode(1, name,0.0);
            rootStringLen=name.length()+1; //name;
        } else if (name.equals("")) {
            root=new PhyloNode(1, "root", Double.parseDouble(len));
            rootStringLen=1+len.length()+1; //:len;
        } else {
            root=new PhyloNode(1, "root", 0.0);
            rootStringLen=name.length()+1+len.length()+1; //name:len;
        }
        
        if (root==null) {System.out.println("The newickTree input may be missing the terminal ';'.");return null;};
        
        //this removes the info of the root node and send the subtree to the recursivity
        //-1 just to remove the last parenthesis
        String subtree=newickTree.substring(0, newickTree.length()-rootStringLen-1);
        //System.out.println("from root subtree:"+subtree);
        addSubTrees(root,subtree);
        
//        long endTime = System.currentTimeMillis();
//        System.out.println("Newick parsing used " + (endTime - startTime) + " ms");

        PhyloTree tree=new PhyloTree(new PhyloTreeModel(root));
        tree.initIndexes();
        return tree;
    }
    
    /**
     * assumes that the newick string representation is the 1st line of the file
     * @param f
     * @return
     * @throws FileNotFoundException 
     */
    public PhyloTree parseNewickTree (File f) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line=null;
        PhyloTree tree=null;
        while ((line=br.readLine())!=null) {
            tree=parseNewickTree(line.trim());
            break;
        }
        
        br.close();
        return tree;
    }
    
    /**
     * write a phylotree as a newick tree
     * @param tree
     * @param fileName 
     * @param withBranchLength 
     * @throws java.io.IOException 
     */
    public void writeNewickTree(PhyloTree tree, File fileName,boolean withBranchLength, boolean writeInternalNodeNames) throws IOException {
        currentElementIndex=1;
        level=0;
        Infos.println("Writing newick tree...");
        try (FileWriter fw = new FileWriter(fileName)) {
            PhyloNode root=(PhyloNode)tree.getModel().getRoot();
            StringBuilder sb = new StringBuilder();
            sb = writerDFS(root,sb,withBranchLength,writeInternalNodeNames);
            fw.append(sb);
        }
        Infos.println("Tree saved: "+fileName.getAbsolutePath());
        
    }
    
    /**
     * depth first search used by the newick writer
     * @param node 
     * @param sb
     * @param withBranchLength
     */
    private StringBuilder writerDFS(PhyloNode node, StringBuilder sb, boolean withBranchLength, boolean writeInternalNodeNames) {
        
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
                if(withBranchLength) {
                    sb.append(":");
                    sb.append(currentNode.getBranchLengthToAncestor());
                }
            } else {
                //System.out.println("LEVEL "+level+" TO "+(++level));
                writerDFS(currentNode,sb,withBranchLength,writeInternalNodeNames);
                //System.out.println("RETURN TO "+(--level) +"  (node:"+node+")");
            }
            //return from recusion or simple leaf,
            //add ',' if more children
            //close level with current node info if no children left
            if (childrenLeft>0) {
                sb.append(",");
            } else {
                sb.append(")");
                if(writeInternalNodeNames) {
                    sb.append(node.getLabel());
                }
                if (withBranchLength) {
                    sb.append(":");
                    sb.append(node.getBranchLengthToAncestor());
                }
            }
        }
        //close the string with a ';' after all root children were passed
        if (node.isRoot() && childrenLeft<1) {
            sb.append(";"); 
        }
        
        //System.out.println(sb+"\n");
        return sb;
    }
    
    
    
    private void addSubTrees(PhyloNode parent, String parentSubtree) {
        
        //System.out.println("CALL "+parent+"   "+parentSubtree);
        
        //search all ',' of nestidness 1, parsing from left to right
        ArrayList<Integer> nestedCommas=new ArrayList();
        int nestedness=0;
        for (int i=0;i<parentSubtree.length();i++) {
            if (parentSubtree.charAt(i)=='(') {
                nestedness++;
                continue;
            } else if (parentSubtree.charAt(i)==')') {
                nestedness--;
                continue;
            }
            if (parentSubtree.charAt(i)==',' && nestedness==1) {
                nestedCommas.add(i);
            }
        }
        //System.out.println(nestedCommas);
        //System.out.println("parent subtrees before split: "+parentSubtree);
        //split the current level by ',' of nesteness 1
        int indexStart=1; //=1 to zap the '(' of the current level
        ArrayList<String> subtrees=new ArrayList<>();
        for (int i = 0; i < nestedCommas.size(); i++) {
            Integer index = nestedCommas.get(i);
            subtrees.add(parentSubtree.substring(indexStart, index));
            indexStart=index+1;
        }
        //System.out.println("subTree size after split: "+subtrees.size());
        //note there is 2 scenarios here: 
        //1. the last element is a subtree which root has a name, (xxx)NAME:LEN
        //2. the last element is a leaf NAME:LEN
        //in second case, we remove the last ')' so that the internalNodePattern matches
        subtrees.add(parentSubtree.substring(indexStart,parentSubtree.length()));
        
        //System.out.println("next subtrees: "+subtrees);
        //for each split element, recursive call, or create a leaf if not a tree
        for (int i = 0; i < subtrees.size(); i++) {
            String subTree=subtrees.get(i);
            Matcher m=null;
            if (subTree.startsWith("(")) {//subtree
                //if internal node has ')label:length' 
                //Remove the last ) for the pattern to work
                if (subTree.charAt(subTree.length()-1)==')') {
                    //System.out.println("Remove last ) !");
                    subTree=subTree.substring(0, subTree.length()-1);
                }
                //System.out.println("before internalnode pattern"+subTree);
                m=internalNodePattern.matcher(subTree);
                if (m.find()) {
                    String name=m.group(1);
                    String len=m.group(2);
                    String nameAlt=m.group(3);//alternative for when the leaf name is a number and no :
                    //System.out.println(name+"::"+len+"::"+nameAlt);
                    if (name==null && len==null) {
                        name=nameAlt;
                        len="";
                    }
                    PhyloNode internalNode=null;
                    if (!len.equals("")) {
                        internalNode=new PhyloNode(++currentElementIndex,name,Double.parseDouble(len));
                    } else {
                        internalNode=new PhyloNode(++currentElementIndex,name,0.1);
                    }
                    parent.add(internalNode);
                    //remove the node info, just send the subtree to next recursive call
                    //System.out.println("   INTERNnodeInfo: name='"+name+"':len='"+len+"'");
                    int nodeInfoLength=-1;
                    if (len.equals("")) {
                        nodeInfoLength=name.length();
                    } else {
                        nodeInfoLength=name.length()+1+len.length();
                    }                    
                    //substring(1,end-1) to send the next 
                    String nextSubTree=subTree.substring(0,subTree.length()-nodeInfoLength);
                    addSubTrees(internalNode, nextSubTree);
                } else {
                    System.out.println("This subtree seems incorrect:\n"+subTree);
                }
            } else {//leaf
                //if internal node was 'leaf)'
                //Remove the last ) for the pattern to work
                if (subTree.charAt(subTree.length()-1)==')') {
                    //System.out.println("Remove last ) !");
                    subTree=subTree.substring(0, subTree.length()-1);
                }                
                //System.out.println("before leaf pattern: "+subTree);
                m=leafPattern.matcher(subTree);
                if (m.find()) {
                    String name=m.group(1);
                    String len=m.group(2);
                    String nameAlt=m.group(3);//alternative for when the leaf name is a number and no :
                    if (name==null && len==null) {
                        name=nameAlt;
                        len="";
                    }
                    PhyloNode leaf=null;
                    if (!len.equals("")) {
                        leaf=new PhyloNode(++currentElementIndex,name,Double.parseDouble(len));
                    } else {
                        leaf=new PhyloNode(++currentElementIndex,name,0.1);
                    }
                    //System.out.println("   LEAFnodeInfo: name='"+name+"':len='"+len+"'");
                    parent.add(leaf);
                } else {
                    //System.out.println("This leaf seems incorrect:\n"+subTree);
                }
            }
            
        }
    
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
        PhyloTree t=new NewickParser().parseNewickTree(treeFASTML);
        long endTime = System.currentTimeMillis();
        System.out.println("Parsing took " + (endTime - startTime) + " milliseconds");

        //test writing
        startTime = System.currentTimeMillis();
        try {
            new NewickParser().writeNewickTree(t, new File("test_fake.tree"),true,true);
        } catch (IOException ex) {
            Logger.getLogger(NewickParser.class.getName()).log(Level.SEVERE, null, ex);
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
