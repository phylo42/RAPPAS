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

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;

/**
 *
 * @author ben
 */
public class TEST_NewickWriter {


    
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
        PhyloTree t=new NewickReader().parseNewickTree2(tree, false, false);
        long endTime = System.currentTimeMillis();
        System.out.println("Parsing took " + (endTime - startTime) + " milliseconds");

        System.out.println(t.getById(0));
        System.out.println(t.getById(10));
        System.out.println(t.getById(8));
        
        //test writing
        startTime = System.currentTimeMillis();
        try {
            NewickWriter nw=new NewickWriter(new File("test_fake.tree"));
            //nw.writeNewickTree(t,true,true,true);
            System.out.println(nw.getNewickTree(t, true, true, true, false));
            nw.close();
            
        } catch (IOException ex) {
            Logger.getLogger(NewickReader.class.getName()).log(Level.SEVERE, null, ex);
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
