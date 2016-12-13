/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tree;

import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;

/**
 *
 * @author ben
 */
public class PhyloTreeModel extends DefaultTreeModel {
    
    public PhyloTreeModel(TreeNode root) {
        super(root);
    }
    
    public PhyloTreeModel(TreeNode root, boolean asksAllowsChildren) {
        super(root, asksAllowsChildren);
    }
    
}
