/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author ben
 */
public class QueryWord extends AbstractWord implements Word, Serializable {
    
    private static final long serialVersionUID = 4022L;
    
    private int originalPosition=-1;

    public QueryWord(byte[] word, int originalPosition) {
        this.word=word;
        this.originalPosition=originalPosition;
    }
    
    /**
     * Currently allows only 1 mutation in the original word \n
     * Very basic approach...
     * @param s
     * @return 
     */
    public List<QueryWord> getMutatedWords(States s) {
        ArrayList<QueryWord> list=new ArrayList<>();
        for (int i = 0; i < word.length; i++) {
            for (byte j = 0; j < s.getStateCount(); j++) {
                if (word[i]==j) {
                    continue;
                }
                byte[] newWord=Arrays.copyOf(word, word.length);
                newWord[i]=j;
                list.add(new QueryWord(newWord, originalPosition));
            }
        }
        return list;
    }

    @Override
    public String toString() {
        return Arrays.toString(word)+":"+originalPosition;
    }
    
   
    public int getOriginalPosition() {
        return this.originalPosition;
    }
    
    
}
