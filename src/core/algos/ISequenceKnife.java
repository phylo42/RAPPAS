package core.algos;

import core.QueryWord;
import inputs.Fasta;

public interface ISequenceKnife {
    public QueryWord getNextWord();
    
    public void init(String seq);
    public void init(Fasta f);

    public byte[] getNextByteWord();
    
    public int getMerCount();
}
