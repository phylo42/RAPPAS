/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import java.util.ArrayList;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.TreeCellRenderer;
import javax.swing.tree.TreeModel;

/**
 * A tree of PhyloNodes, partially based on the original JTree class
 * @author ben
 */
public interface Tree {
    
    public TreeModel getModel();
 
    public PhyloNode getByName(String nodeName);
    
    public PhyloNode getById(int id);
    
    public int getNodeCount();
    
    public int getLeavesCount();

    public ArrayList<Integer> getInternalNodesByDFS();
    
    public ArrayList<Integer> getLeavesByDFS();
    
    public ArrayList<Integer> getNodeIdsByDFS();
   
    public ArrayList<String> getLabelsByDFS();
    
    public TreeCellRenderer getCellRenderer();
    
    public void addTreeSelectionListener(TreeSelectionListener tsl);
    
    public void displayTree();
}
