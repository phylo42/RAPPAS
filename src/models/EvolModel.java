/*
 * Copyright (C) 2018 ben
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
package models;

import core.AAStates;
import core.States;
import java.util.HashMap;

/**
 * manages modelId used during ancestral sequence reconstruction
 * @author ben
 */
public class EvolModel {
    
    //parameters related to AR models
    //note only models supported by both phyml and paml are supported in RAPPAS
    //values 0 and 20 are let free for potential custom modelId support in the future.
    //nucl models: 7 commons
    public static final int MODEL_JC69=1; //phyml=JC69 ; baseml=JC69
    public static final int MODEL_K80=2; //phyml=K80  ; baseml=K80
    public static final int MODEL_F81=3; //phyml=F81  ; baseml=F81
    public static final int MODEL_F84=4; //phyml=F84  ; baseml=F84
    public static final int MODEL_HKY85=5; //phyml=HKY85  ; baseml=HKY85
    public static final int MODEL_TN93=6; //phyml=TN93  ; baseml=TN93
    public static final int MODEL_GTR=7;  //phyml=GTR  ; baseml=REV
    //prot models 9 commons
    public static final int MODEL_LG=21; //phyml=LG ; baseml=lg.dat
    public static final int MODEL_WAG=22; //phyml=WAG  ; baseml=wag.dat
    public static final int MODEL_JTT=23; //phyml=JTT  ; baseml=jones.dat
    public static final int MODEL_DAYHOFF=24; //phyml=Dayhoff  ; baseml=dayhoff.dat
    public static final int MODEL_DCMUT=25; //phyml=DCMut  ; baseml=dayhoff_dimut.dat
    public static final int MODEL_CPREV=26; //phyml=CpREV  ; baseml=cpREV10.dat
    public static final int MODEL_MTMAM=27;  //phyml=MtMam  ; baseml=mtmam.dat
    public static final int MODEL_MTREV=28;  //phyml=mMtREV  ; baseml=mtREV24.dat
    public static final int MODEL_MTART=29;  //phyml=MtArt  ; baseml=mtart.dat
    
    //correspondance to paml modelId files
    HashMap<Integer, String> pamlModelFiles=new HashMap<>();
        
    //default parameters
    public int modelId=MODEL_GTR;
    public String modelString="GTR"; //string representation of models directly callable in phyml command-line
    public float alpha=1.0f;
    public int categories=4;

    
    
    

    /**
     * default constructor, return a basic GTR or LG modelId
     * @param s 
     */
    private EvolModel(States s) {
        pamlProtModelsFiles();
        if (s instanceof AAStates) {
            this.modelId=MODEL_LG;
            this.modelString="LG";
        } else {
            this.modelId=MODEL_GTR;
            this.modelString="GTR";
        }
    }
    /**
     * 
     * @param model
     * @param alpha
     * @param categories 
     */
    public EvolModel(String model, float alpha, int categories) {
        pamlProtModelsFiles();
        this.alpha=alpha;
        this.categories=categories;
        chooseFromString(model);            
    }
    
    /**
     * 
     * @param model one of Model.MODEL_JC69 , Model.MODEL_JC69 etc
     * @param alpha
     * @param categories 
     */
    public EvolModel(int model, float alpha, int categories) {
        pamlProtModelsFiles();
        this.alpha=alpha;
        this.categories=categories;
        this.modelId=model;
    }
    
    public boolean isProteinModel() {
        return (this.modelId>19);
    }
    
    
    /**
     * using a string, choose the correct modelId definition
     * @return 
     */
    private void chooseFromString(String val) {
        switch (val) {
            //nucleotidic models: 
            case "JC69":
                this.modelId=MODEL_JC69;
                this.modelString="JC69";
                break;
            case "HKY85":
                this.modelId=MODEL_HKY85;
                this.modelString="HKY85";
                break;
            case "K80":
                this.modelId=MODEL_K80;
                this.modelString="K80";
                break;
            case "F81":
                this.modelId=MODEL_F81;
                this.modelString="F81";
                break;
            case "TN93":
                this.modelId=MODEL_TN93;
                this.modelString="TN93";
                break;
            case "GTR":
                this.modelId=MODEL_GTR;
                this.modelString="GTR";
                break;
            //proteic models
            case "LG":
                this.modelId=MODEL_LG;
                this.modelString="LG";
                break;
            case "WAG":
                this.modelId=MODEL_WAG;
                this.modelString="WAG";
                break;
            case "JTT":
                this.modelId=MODEL_JTT;
                this.modelString="JTT";
                break;
            case "Dayhoff":
                this.modelId=MODEL_DAYHOFF;
                this.modelString="Dayhoff";
                break;
            case "DCMut":
                this.modelId=MODEL_DCMUT;
                this.modelString="DCMut";
                break;
            case "CpREV":
                this.modelId=MODEL_CPREV;
                this.modelString="CpREV";
                break;
            case "MtMam":
                this.modelId=MODEL_MTMAM;
                this.modelString="MtMam";
                break;
            case "MtREV":
                this.modelId=MODEL_MTREV;
                this.modelString="mMtREV";
                break;
            case "MtArt":
                this.modelId=MODEL_MTART;
                this.modelString="MtArt";
                break;   
            default:
                System.out.println("Model name not recognized. Please use one of:");
                System.out.println("*nucl : JC69, HKY85, K80, F81, TN93, GTR");
                System.out.println("*amino: LG, WAG, JTT, Dayhoff, DCMut, CpREV, mMtREV, MtMam, MtArt");  
                System.exit(1);
        }

    }
    
    /**
     * defines equivalents for paml control (*.ctl) files
     */
    private void pamlProtModelsFiles() {
        //nucl
        //0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
        pamlModelFiles.put(MODEL_JC69,"0");
        pamlModelFiles.put(MODEL_HKY85,"4");
        pamlModelFiles.put(MODEL_K80,"1");
        pamlModelFiles.put(MODEL_F81,"2");
        pamlModelFiles.put(MODEL_TN93,"6");
        pamlModelFiles.put(MODEL_GTR,"7");
        //prot
        pamlModelFiles.put(MODEL_LG,"lg.dat");
        pamlModelFiles.put(MODEL_WAG,"wag.dat");
        pamlModelFiles.put(MODEL_JTT,"jones.dat");
        pamlModelFiles.put(MODEL_DAYHOFF,"dayhoff.dat");
        pamlModelFiles.put(MODEL_DCMUT,"dayhoff_dimut.dat");
        pamlModelFiles.put(MODEL_CPREV,"cpREV10.dat");
        pamlModelFiles.put(MODEL_MTMAM,"mtmam.dat");
        pamlModelFiles.put(MODEL_MTREV,"mtREV24.dat");
        pamlModelFiles.put(MODEL_MTART,"mtart.dat");
    }
    
    
    
    public String getPAMLEquivalent() {
        return pamlModelFiles.get(modelId);
    }
    
    
    /**
     * build default model during arguments parsing.
     * @param s one of States.STATES_DNA or States.STATES_PROTEIN
     * @return 
     */
    public static EvolModel getDefault(States s) {
        return new EvolModel(s);
    }

    @Override
    public String toString() {
        return "m="+modelString+";a="+alpha+";c="+categories+";";
    }
    
    
    
    
}
