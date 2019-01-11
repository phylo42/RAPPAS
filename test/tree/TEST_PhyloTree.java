/*
 * Copyright (C) 2019 ben
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package tree;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ben
 */
public class TEST_PhyloTree {
    
    public static void main(String[] args) {
        
        System.setProperty("debug.verbose", "1");
        
        String t_basic_unrooted="((L1:1,L2:2)X:3,(L3:2,L4:2)Y:1,L5:4)Z:4;"; //unrooted
        String t_basic_rZX="((L1:1,L2:2)X:4,((L3:2,L4:2)Y:1,L5:4)Z:4)root:0;";//rooted on Z-X
        String t_basic_rZX_inv="(((L3:2,L4:2)Y:1,L5:4)Z:4,(L1:1,L2:2)X:4)root:0;";//rooted on Z-X, inversed children of root
        String t_basic_rZL5="(L5:4,((L3:2,L4:2)Y:1,(L1:1,L2:2)X:4)Z:4)root:0;"; //rooted on Z-L5
        
        String treeBase="(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927)homo:0.08386):0.06124):0.15057)primates:0.54939,Mouse:1.21460):0.10;";

        String treePAML=    "(5_KJ473819,((13_KF843851,(10_KJ473805,4_KJ473804)34)33,((12_KF843855,14_KF843858)36,(26_AB781795,"
                + "(23_AB781793,(17_AB907634,(2_KF530123,(25_AB781796,(11_AF516906,((6_AF124986,(8_EU769559,(9_EU769560,21_AB907631)46)45)44,"
                + "((20_AB907632,22_AB907633)48,((30_DQ431016,29_DQ431014)50,((27_AY874541,(1_KC175339,(3_JN856008,7_AF124992)54)53)52,"
                + "(28_AY878324,(16_AB907625,(15_AB781791,(19_AB907630,(18_AB907628,24_AB781792)59)"
                + "58)57)56)55)51)49)47)43)42)41)40)39)38)37)35)32)31;";
        
        String treePAMLUnrooted= "(5_KJ473819,(13_KF843851,(10_KJ473805,4_KJ473804)34)33,((12_KF843855,14_KF843858)36,(26_AB781795,"
                + "(23_AB781793,(17_AB907634,(2_KF530123,(25_AB781796,(11_AF516906,((6_AF124986,(8_EU769559,(9_EU769560,21_AB907631)46)45)44,"
                + "((20_AB907632,22_AB907633)48,((30_DQ431016,29_DQ431014)50,((27_AY874541,(1_KC175339,(3_JN856008,7_AF124992)54)53)52,"
                + "(28_AY878324,(16_AB907625,(15_AB781791,(19_AB907630,(18_AB907628,24_AB781792)59)"
                + "58)57)56)55)51)49)47)43)42)41)40)39)38)37)35)32;";
        
        PhyloTree tree = NewickReader.parseNewickTree2(t_basic_rZX, false, false);
        tree.initIndexes();
        
        PhyloNode Y = tree.getByName("Y");
        //PhyloNode L5 = tree.getByName("L5");        
        tree.rerootTree(Y,true);

//        PhyloNode Y = tree.getByName("homo");
//        tree.rerootTree(Y, true);

//        PhyloNode Y = tree.getByName("58");
//        //PhyloNode L5 = tree.getByName("L5");        
//        tree.rerootTree(Y);


        try {
            System.out.println(new NewickWriter().getNewickTree(tree, true, true, false, false));
        } catch (IOException ex) {
            Logger.getLogger(TEST_PhyloTree.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        tree.displayTree();
        
    }
    
}
