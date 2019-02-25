/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos.sorting;

import core.DNAStatesShifted;
import core.QueryWord;
import core.algos.SequenceKnife;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author ben
 */
public class TEST_SequenceKnife {
    public static void main(String[] args) {
        SequenceKnife knife=new SequenceKnife("ATCGCTGATCGATCGA", 7, 4, new DNAStatesShifted(), SequenceKnife.SAMPLING_SEQUENTIAL);
        System.out.println(Arrays.toString(knife.getMerOrder()));
        knife=new SequenceKnife("ATCGCTGATCGATCGA", 7, 4, new DNAStatesShifted(), SequenceKnife.SAMPLING_STOCHASTIC);
        knife.forceSeed(12345);
        System.out.println(Arrays.toString(knife.getMerOrder()));
        knife=new SequenceKnife("ATCGCTGATCGATCGA", 7, 4, new DNAStatesShifted(), SequenceKnife.SAMPLING_LINEAR);
        System.out.println(Arrays.toString(knife.getMerOrder()));
        
        QueryWord w=null;
        while ((w=knife.getNextWord())!=null) {
            System.out.println(w);
            List<QueryWord> mutatedWords = w.getMutatedWords(new DNAStatesShifted());
            for (int i = 0; i < mutatedWords.size(); i++) {
                System.out.println("   "+ mutatedWords.get(i));
                
            }
        }
    }
}
