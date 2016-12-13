/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignement;

/**
 * an alignment's partition, generally defining a gene\n
 *
 * @author ben
 */
public class Partition {

    private int length=-1;
    private int start=-1;
    private int end=-1;
    private String name;

    public Partition(String name, int start, int end) {
        this.name = name;
        this.start=start;
        this.end=end;
        this.length=end-start+1;
    }

    public int getEnd() {
        return end;
    }

    public String getName() {
        return name;
    }

    public int getStart() {
        return start;
    }

    public void setName(String name) {
        this.name = name;
    }

    @Override
    public String toString() {
        return "name:"+name+" start:"+start+" end:"+end+" "+super.toString();
    }

    
    
    
    
}
