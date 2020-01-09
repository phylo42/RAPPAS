/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main_v2;

import etc.Environement;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * parse the arguments when command_line version is used
 * @author linard
 */
public class ArgumentsParser_v2 {

    public static final int DBBUILD_PHASE=0;
    public static final int PLACEMENT_PHASE=1;
    
    public static final int DB_FULL=0;
    public static final int DB_MEDIUM=0;
    public static final int DB_SMALL=0;
    
    public static final int STATES_DNA=1;
    public static final int STATES_PROTEIN=2;
    
    private HashMap<Integer,String> argsMap=null;
    
    String version="debug";


    //RAPPAS general parameters
    public int phase=DBBUILD_PHASE;
    public File workingDir=null;//current directory by default, see below
    public int verbose=0; //default, no verbosity
    public int states=STATES_DNA;
    boolean convertUO=false; //must be true so that U/O/X are converted in C/L/-
    
    //RAPPAS parameters for alignment reduction
    public boolean reduction=true;
    public File reducedAlignFile=null;
    public double reductionRatio=0.99;
    
    //RAPPAS parameters for DB build
    public int k=8; //default=8
    public float omega=1.5f; //default=1., quantity allowing to modulate the threshold
    public int ghostsAmount=1;  //default =1
    public File alignmentFile=null;
    public File treeFile=null;
    public boolean forceRooting=false;
    public String modelString=null; //if not set by options, default will be used (see Main)
    public float alpha=1.0f;
    public int categories=4;
    public String arparameters=null;
    public String dbfilename=null;
    public int threads=4; // #threads for AR, currently only used by raxml-ng

    //passed for debugging in DB_build
    public File ARBinary=new File("phyml"); //default = phyml command line
    public File ARDirToUse=null;  // used to load directly an already computed AR
    public File exTreeDir=null; //not functionnal anymore
    public boolean onlyAR=false;
    public boolean onlyARInput=false;
    public boolean builddbfull=false; //default=false, as dbfull is quite useless with the current algo
    public boolean noCalibration=true; //skip calibration step
    public boolean dbInRAM=false; //do not write DB in a file and immediately places reads passed with -q
    public boolean unionHash=true; //if false, use old positionnal hash
    public boolean onlyFakeNodes=true; //if false, uses ancestral kmers of original nodes
    public boolean doGapJumps=false; //take gap jumps into account when building kmers
    public boolean limitTo1Jump=true; //only allow a 1st jump, not jump combinations
    public float gapJumpThreshold=0.3f; //gap jumps are activated if >30% gaps in the ref alignment
    public boolean onlyX1Nodes=false;
    public boolean jsondb=false;
    public boolean acceptUnrootedRefTree=false;
    
    //RAPPAS parameters for placement
    public int minOverlap=100; //used in entropy computation
    public List<File> queriesFiles=null;
    public File databaseFile=null;
    public Float nsBound=null; //do not consider threshold calculated by formula but this particular threshold, null if not set
    public int keepAtMost=7; //as in pplacer
    public float keepFactor=0.01f; //as in pplacer
    public boolean guppyCompatible=false;
    public boolean treatAmbiguities=true; // treat ambiguities 
    public boolean treatAmbiguitiesWithMax=false; //default: treated with means

    
    //call string
    public String callString=null;
    
    public ArgumentsParser_v2(String[] args,String version) {
        this.version=version;
        argsMap=new HashMap<Integer,String>();
        StringBuilder sb=new StringBuilder();
        //args with priority
        for (int i=0;i<args.length;i++) {
            if (args[i].equals("--help") || args[i].equals("-h"))
                showHelpAndExit();
            sb.append(" "+args[i]);
            argsMap.put(i, args[i]);
        }
        callString=sb.toString();
        sb=null;
        
        //if no phase, bad
        if ((!argsMap.containsValue("-p")) && (!argsMap.containsValue("--phase")) ) {
            System.out.println("Cannot find 'phase' (option -p). See help (-h) for details.");
            System.exit(1);
        }
        
        try {
            //if arguments values are correct
            loadParameters();
        } catch (Exception ex) {
            Logger.getLogger(ArgumentsParser_v2.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("Unexpected parameters !");
            System.out.println("Something went wrong with the given command line parameters... ");
            System.exit(1);
        }

    }


    /**
     * check if arguments assoiations and values are correct
     */
    private void loadParameters() throws Exception {

        //general loop, must pass a 1st time on all arguments to set basal
        //options which will influence other options (ex: protein vs nucleotides analysis)
        boolean wGiven=false;
        boolean statesGiven=false;
        boolean modelGiven=false;
        for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
            int index=it.next();
            
            //check if --phase associated to correct minimum values
            if (argsMap.get(index).equals("-p") || argsMap.get(index).equals("--phase")) {
                if (argsMap.get(index+1).equalsIgnoreCase("b")) {this.phase=DBBUILD_PHASE;}
                if (argsMap.get(index+1).equalsIgnoreCase("p")) {this.phase=PLACEMENT_PHASE;}
            }
            
            //check verbosity
            if (argsMap.get(index).equals("-v") || argsMap.get(index).equals("--verbosity")) {
                try {
                    this.verbose=Integer.parseInt(argsMap.get(index+1));
                    if (verbose>0) { //for now, only one verbosity level
                        System.setProperty("debug.verbose", "1");
                    }
                } catch (NumberFormatException ex) {
                    System.out.println("Cannot parse '-v (--verbosity)' as an integer value.");
                    System.exit(1);
                }
            }
            //check working directory
            if (argsMap.get(index).equals("-w") || argsMap.get(index).equals("--workdir")) {
                    System.out.println("wordDir: "+(new File(argsMap.get(index+1)).getAbsolutePath()));
                    File wd=new File(argsMap.get(index+1));
                    if (wd.isDirectory() && wd.canRead() && wd.canWrite()) {
                        this.workingDir=wd;
                        wGiven=true;
                    } else {
                        System.out.println("Cowardly refusing to use a working directory that does not exist !");
                        System.out.println("Will not read/write from to this directory: "+wd.getAbsolutePath());
                        System.exit(1);
                    }
            }        
            
            //test --convertUO parameter
            if (argsMap.get(index).equals("--convertUO")) {
                this.convertUO=true;
            }
            
            //test -m modelId: 
            if (argsMap.get(index).equals("--model") || argsMap.get(index).equals("-m")) {
                String fVal=argsMap.get(index+1);
                this.modelString=fVal;
                modelGiven=true;
            }
            

        }
        if ( !(phase==DBBUILD_PHASE) && !(phase==PLACEMENT_PHASE) ) {
            System.out.println("Unexpected -p (--phase) value, must be 'b' (build_db) or 'p' (place) !");
            System.exit(1);
        }
        
        if (!wGiven) {
            this.workingDir=Environement.getCurrentDirectory().toFile();
        }
        
        ////from here, we know in which phase the program will be launched///////
        ////////////////////////////////////////////////////////////////////////
        //check if -a ,-t, -k given when phase=b
        switch (phase) {
            
            case DBBUILD_PHASE:
                
                boolean kGiven=false;
                boolean omegaGiven=false;
                boolean reductionGiven=false;
                boolean fGiven=false;
                boolean alignGiven=false;
                boolean treeGiven=false;
                
                for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
                    int index=it.next();
                    
                    //test -s parameter, this must be tested before many other options as 
                    //analyses differ when doing protein or nucleotides analysis
                    if (argsMap.get(index).equals("--states") || argsMap.get(index).equals("-s")) {
                        if (argsMap.get(index+1).equalsIgnoreCase("nucl")) {
                            statesGiven=true;
                            this.states=STATES_DNA;
                        } else if (argsMap.get(index+1).equalsIgnoreCase("amino")) {
                            statesGiven=true;
                            this.states=STATES_PROTEIN;
                        } else {
                            System.out.println("Unexpected -s (--states) value, must be one of [nucl|amino].");
                            System.exit(1);
                        }
                    }
                    
                    //test -r parameter
                    if (argsMap.get(index).equals("--refalign") || argsMap.get(index).equals("-r")) {
                        File alignment=new File(argsMap.get(index+1));
                        if (alignment.isFile() && alignment.canRead()) {
                            this.alignmentFile=alignment;
                            alignGiven=true;
                        } else {
                            System.out.println(alignment.getAbsolutePath());
                            System.out.println("Cannot open alignment: Not a file or no read permission.");
                            System.exit(1);
                        }
                    }
                    //test -t parameter
                    if (argsMap.get(index).equals("--reftree") || argsMap.get(index).equals("-t")) {
                        File tree=new File(argsMap.get(index+1));
                        if (tree.isFile() && tree.canRead()) {
                            this.treeFile=tree;
                            treeGiven=true;
                        } else {
                            System.out.println(tree.getAbsolutePath());
                            System.out.println("Cannot open tree: Not a file or no read permission.");
                            System.exit(1);
                        }
                    }
                    //test -k parameter
                    if (argsMap.get(index).equals("--k") || argsMap.get(index).equals("-k")) {
                        String kVal=argsMap.get(index+1);
                        try {
                            this.k=Integer.parseInt(kVal);
                            kGiven=true;
                            if (this.k<3) {
                                this.k=3;
                                System.out.println("--k set to 3 (minimum authorised values) .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-k (--k)' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    //test -a parameter
                    if (argsMap.get(index).equals("--omega")) {
                        String omegaVal=argsMap.get(index+1);
                        try {
                            this.omega=Float.parseFloat(omegaVal);
                            omegaGiven=true;
                            if (this.omega<0) {
                                this.omega=1.0f;
                                System.out.println("--omega set to 1.0 .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--omega' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    
                    
                    //test -f parameter
                    if (argsMap.get(index).equals("--ghosts") || argsMap.get(index).equals("-g")) {
                        String fVal=argsMap.get(index+1);
                        try {
                            this.ghostsAmount=Integer.parseInt(fVal);
                            fGiven=true;
                            if (this.ghostsAmount<1) {
                                this.ghostsAmount=1;
                                System.out.println("--ghosts set to 1.");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-g (--ghosts)' as an integer value.");
                            System.exit(1);
                        }
                       
                    }
                    
                    //test -b parameter
                    if (argsMap.get(index).equals("--arbinary") || argsMap.get(index).equals("-b")) {
                        File ARBinaryFile=new File(argsMap.get(index+1));
                        System.out.println("Using AR binary provided by user: "+ARBinaryFile.getAbsolutePath());
                        if (ARBinaryFile.isFile() && ARBinaryFile.canExecute()) {
                            this.ARBinary=ARBinaryFile;
                        } else {
                            System.out.println("Cannot execute binary loaded through option --arbinary: Not a file or no execution permission.");
                            System.exit(1);
                        }
                    }
                    
                    //test -a parameter
                    if (argsMap.get(index).equals("--alpha") || argsMap.get(index).equals("-a")) {
                        String fVal=argsMap.get(index+1);
                        try {
                            this.alpha=Float.parseFloat(fVal);
                            if (this.alpha<0.0) {
                                this.alpha=1.0f;
                                System.out.println("--alpha set to 1.0 .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-a (--alpha)' as a float value.");
                            System.exit(1);
                        }
                       
                    }
                    
                    //test -c parameter
                    if (argsMap.get(index).equals("--categories") || argsMap.get(index).equals("-c")) {
                        String fVal=argsMap.get(index+1);
                        try {
                            this.categories=Integer.parseInt(fVal);
                            if (this.categories<0) {
                                this.categories=4;
                                System.out.println("--categories set to 4 .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-c (--categories)' as an integer value.");
                            System.exit(1);
                        }
                       
                    }
                    
                    //test --dbfilename
                    if (argsMap.get(index).equals("--dbfilename")) {
                        this.dbfilename=argsMap.get(index+1);
                        System.out.println("Output DB will be renamed as: "+this.dbfilename);
                    }
                  
                    
                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS
                    
                    //test --no-reduction
                    if (argsMap.get(index).equals("--no-reduction")) {
                        System.out.println("Original alignment will be used (no gapped columns removed).");
                        this.reduction=false;
                    }
                    
                    //test --write-reduction
                    if (argsMap.get(index).equals("--write-reduction")) {
                        File f=new File(argsMap.get(index+1));
                        if (f.isFile() && f.canWrite()) {
                            this.reducedAlignFile=f;
                        } else {
                            System.out.println("Cannot write file given via option --write-reduction: Not a file or no write permission.");
                            System.exit(1);
                        }
                    }
                    
                    //test --ratio-reduction
                    if (argsMap.get(index).equals("--ratio-reduction")) {
                        String alphaVal=argsMap.get(index+1);
                        try {
                            this.reductionRatio=Double.parseDouble(alphaVal);
                            reductionGiven=true;
                            if (this.reductionRatio>1) {
                                this.reductionRatio=1;
                                System.out.println("--ratio-reduction set to 1.0 .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--ratio-reduction' as a double value.");
                            System.exit(1);
                        }
                        
                    }
                    
                    //test --ardir parameter
                    if (argsMap.get(index).equals("--ardir")) {
                        File ARDir=new File(argsMap.get(index+1));
                        System.out.println("AR provided by user: "+ARDir.getAbsolutePath());
                        if (ARDir.isDirectory() && ARDir.canRead()) {
                            this.ARDirToUse=ARDir;
                        } else {
                            System.out.println("Cannot open directory given through option --ardir: Not a directory or no read permission.");
                            System.exit(1);
                        }
                    }
                   
                    
                    //test --extree parameter
                    if (argsMap.get(index).equals("--extree")) {
                        File exTreeDir=new File(argsMap.get(index+1));
                        System.out.println("Using extended trees provided by user: "+exTreeDir.getAbsolutePath());
                        System.out.println("exTreeDir:"+exTreeDir.isDirectory());
                        System.out.println("exTreeDir:"+exTreeDir.canRead());
                        if (exTreeDir.isDirectory()&& exTreeDir.canRead()) {
                            this.exTreeDir=exTreeDir;
                        } else {
                            System.out.println("Cannot open directory given through option --extree: Not a directory or no read permission.");
                            System.exit(1);
                        }
                    }
                    
                    //test --dbfull parameter
                    if (argsMap.get(index).equals("--dbfull")) {
                        this.builddbfull=true;
                    }
                    
                    //test --force-root parameter
                    if (argsMap.get(index).equals("--force-root")) {
                        this.forceRooting=true;
                    }
                    
                    //test --dbinram parameter
                    if (argsMap.get(index).equals("--dbinram")) {
                        this.dbInRAM=true;
                    }
                    
                    //test -q parameter (to use with --dbinram)
                    if (argsMap.get(index).equals("--queries") || argsMap.get(index).equals("-q")) {
                        this.queriesFiles=new ArrayList<>();
                        //split eventually if these are filenames separated by ','
                        String[] elts=argsMap.get(index+1).split(",");
                        for (int i = 0; i < elts.length; i++) {
                            String elt = elts[i];
                            File query=new File(elt);
                            if (query.isFile() && query.canRead()) {
                                queriesFiles.add(new File(elt));
                            } else {
                                System.out.println(query.getAbsolutePath());
                                System.out.println("Cannot open query file: Not a file or no read permission.");
                                System.out.println("File: "+elt);
                                System.exit(1);
                            }
                        }
                    }                  
                    
                    //test --nsbound parameter (to use with --dbinram)
                    if (argsMap.get(index).equals("--nsbound")) {
                        String nsBoundVal=argsMap.get(index+1);
                        try {
                            this.nsBound=Float.parseFloat(nsBoundVal);
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--nsbound' as a float value.");
                            System.exit(1);
                        }
                    }
                    
                    //test --nocalib parameter
                    if (argsMap.get(index).equals("--calibration")) {
                        this.noCalibration=false;
                    }
                    
                    //test --poshash parameter
                    if (argsMap.get(index).equals("--poshash")) {
                        this.unionHash=false;
                    }
                    
                    //test --original-nodes
                    if (argsMap.get(index).equals("--original-nodes")) {
                        System.out.println("DB will also contain ancestral kmers associated to original nodes (--original-nodes).");
                        this.onlyFakeNodes=false;
                    }
                    
                    //test --keep-at-most
                    if (argsMap.get(index).equals("--keep-at-most")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.keepAtMost=Integer.parseInt(val);
                            if (this.keepAtMost<1) {
                                this.keepAtMost=1;
                                System.out.println("--keep-at-most set to 1 .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--keep-at-most' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    //test --keep-factor
                    if (argsMap.get(index).equals("--keep-factor")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.keepFactor=Float.parseFloat(val);
                            if (this.keepFactor>1) {
                                this.keepFactor=1;
                                System.out.println("--keep-factor set to 1.0 .");
                            }                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--keep-factor' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    //test --no-gap-jumps
                    if (argsMap.get(index).equals("--force-gap-jump")) {
                        this.doGapJumps=true;
                        System.out.println("Forced gap jump activation.");
                    }
                    //test --do-n-jumps
                    if (argsMap.get(index).equals("--do-n-jumps")) {
                        this.limitTo1Jump=false;
                        System.out.println("N gaps interval will be considered.");
                    }
                    //test --gap-jump-thresh
                    if (argsMap.get(index).equals("--gap-jumps-thresh")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.gapJumpThreshold=Float.parseFloat(val);
                            if (this.gapJumpThreshold>1.0f) {
                                this.gapJumpThreshold=1.0f;
                                System.out.println("--gap-jump-thresh set to 1.0");
                            }
                            if (this.gapJumpThreshold<0.0f) {
                                this.gapJumpThreshold=0.0f;
                                System.out.println("--gap-jump-thresh set to 0.0");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--gap-jump-thresh' as a float value.");
                            System.exit(1);
                        }
                    }             
                    
                    //test --arparameters
                    if (argsMap.get(index).equals("--arparameters")) {
                        String val=argsMap.get(index+1);
                        this.arparameters=val;
                        System.out.println("AR paramaters fixed by user : '"+arparameters+"'");
                        if (this.arparameters.length()<1) {
                            System.out.println("AR parameters string is less than 1 character, please check.");
                            System.exit(1);
                        }
                        //forbidden params, because set by RAPPAS
                        if (    this.arparameters.contains("--ancestral ") ||
                                this.arparameters.contains("-i ") ||
                                this.arparameters.contains("--input ") ||
                                this.arparameters.contains("-u ") ||
                                this.arparameters.contains("--inputtree ")
                                ) {
                            System.out.println("AR parameters string set with '--arparameters' contains a forbidden option (--ancestral,-i,--input,-u,--inputtree).\nThose are always managed by RAPPAS, please remove them.");
                            System.exit(1);
                        }
                    }
                    
                    //test --onlyX1
                    if (argsMap.get(index).equals("--onlyX1")) {
                        this.onlyX1Nodes=true;
                        System.out.println("Only X1 nodes will be considered.");
                    }
                    
                    //test --jsondb
                    if (argsMap.get(index).equals("--jsondb")) {
                        this.jsondb=true;
                        System.out.println("DB will be written as json.");
                    }
                    
                    //test --use_unrooted 
                    if (argsMap.get(index).equals("--use_unrooted")) {
                        this.acceptUnrootedRefTree=true;
                        System.out.println("User confirmed to use an unrooted tree.");
                    }
                    
                    //test --aronly 
                    if (argsMap.get(index).equals("--aronly")) {
                        this.onlyAR=true;
                        System.out.println("Only extended_tree and AR will be built.");
                    }
                    
                    //test --arinputonly 
                    if (argsMap.get(index).equals("--arinputonly")) {
                        this.onlyARInput=true;
                        System.out.println("Only extended_tree and AR inputs will be built.");
                    }

                    //test --threads parameter
                    if (argsMap.get(index).equals("--threads")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.threads=Integer.parseInt(val);
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--threads' as an integer value.");
                            System.exit(1);
                        }
                    }
                    
                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS END HERE
                    
                    
                }
                
                ///////////////////////////////////////////////
                //do not continue if the type of analysis is not set
                ///////////////////////////////////////////////
                if (!statesGiven) {System.out.println("Analysis states not found. Use option -s (--states) and one of 'nucl' or 'amino'.");System.exit(1);}
                ///////////////////////////////////////////////
                //do not continue if align/tree not set (-r / -t)
                ///////////////////////////////////////////////
                if (!alignGiven) {System.out.println("Reference alignment not set correctly (use -r) ."); System.exit(1);}
                if (!treeGiven) {System.out.println("Reference tree not set correctly (use -t) ."); System.exit(1);}
                ///////////////////////////////////////////////
                //use defaults if -k,-a,-f not set
                ///////////////////////////////////////////////
                if (!kGiven) {System.out.println("Default k="+this.k+" will be used.");}
                if (!omegaGiven) {System.out.println("Default omega="+this.omega+" will be used.");}
                if (!reductionGiven) {System.out.println("Default ratio-reduction="+this.reductionRatio+" will be used.");}
                if (!fGiven) {
                    System.out.println("Default injectionPerBranch="+this.ghostsAmount+" will be used.");
                }
                
                break;
                
            //check if -q given when phase=p    
            case PLACEMENT_PHASE:
                
                boolean dGiven=false;
                boolean qGiven=false;
                
                for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
                    int index=it.next();
                    
                    //test -q parameter
                    if (argsMap.get(index).equals("--queries") || argsMap.get(index).equals("-q")) {
                        this.queriesFiles=new ArrayList<>();
                        //split eventually if these are filenames separated by ','
                        String[] elts=argsMap.get(index+1).split(",");
                        for (int i = 0; i < elts.length; i++) {
                            String elt = elts[i];
                            File query=new File(elt);
                            if (query.isFile() && query.canRead()) {
                                queriesFiles.add(new File(elt));
                            } else {
                                System.out.println(query.getAbsolutePath());
                                System.out.println("Cannot open query file: Not a file or no read permission.");
                                System.out.println("File: "+elt);
                                System.exit(1);
                            }
                        }
                        if (queriesFiles.size()>0) {
                            qGiven=true;
                        }
                    }
                    
                    //test -d parameter
                    if (argsMap.get(index).equals("--database") || argsMap.get(index).equals("-d")) {
                        File database=new File(argsMap.get(index+1));
                        if (database.isFile() && database.canRead()) {
                            this.databaseFile=database;
                            dGiven=true;
                        } else {
                            System.out.println(database.getAbsolutePath());
                            System.out.println("Cannot open database: Not a file or no read permission.");
                            System.exit(1);
                        }
                    }
                    
                    

                    
                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS
                    
                    //test --nsbound parameter
                    if (argsMap.get(index).equals("--nsbound")) {
                        String nsBoundVal=argsMap.get(index+1);
                        try {
                            this.nsBound=Float.parseFloat(nsBoundVal);
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--nsbound' as a float value.");
                            System.exit(1);
                        }
                    }
                    
                    //test --keep-at-most
                    if (argsMap.get(index).equals("--keep-at-most")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.keepAtMost=Integer.parseInt(val);
                            if (this.keepAtMost<1) {
                                this.keepAtMost=1;
                                System.out.println("--keep-at-most set to 1 .");
                            }
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--keep-at-most' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    //test --keep-factor
                    if (argsMap.get(index).equals("--keep-factor")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.keepFactor=Float.parseFloat(val);
                            if (this.keepFactor>1) {
                                this.keepFactor=1;
                                System.out.println("--keep-factor set to 1.0 .");
                            }                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--keep-factor' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    //test --guppy-compat
                    if (argsMap.get(index).equals("--guppy-compat")) {
                        this.guppyCompatible=true;
                        System.out.println("Jplace format changed to be guppy-compatible.");
                    }
                    //test --noamb
                    if (argsMap.get(index).equals("--noamb")) {
                        this.treatAmbiguities=false;
                        System.out.println("Ambiguities will not be treated.");
                    }
                    //test --ambwithmax
                    if (argsMap.get(index).equals("--ambwithmax")) {
                        this.treatAmbiguitiesWithMax=true;
                        System.out.println("Ambiguities will be treated with max(w'), not mean(w').");
                    }
                    
                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS END HERE
                    
                }
                
                //check if both -d and -q were set
                if (!dGiven) {System.out.println("User did not provided a reference database (-d).");System.exit(1);}
                if (!qGiven) {System.out.println("User did not provided a query (-q).");System.exit(1);}
                break;
                
            //mode was not set, we must exit
            default:
                System.out.println("You should never arrive here. If yes, welcome to a world of ... don't know, probably another dimension.");
                break;
        }
        

    }

    /**
     * print help and exit execution
     */
    private void showHelpAndExit() {
        System.out.println("################################################");
        System.out.println("## RAPPAS v"+version);
        System.out.println("## ---------------------------------------------");
        System.out.println("## Rapid Alignment-free Phylogenetic Placement ");
        System.out.println("## via Ancestral Sequences");
        System.out.println("## Linard B, Swenson KM, Pardi F.");
        System.out.println("## LIRMM, Univ. of Montpellier, CNRS, France");
        System.out.println("## Citation:");
        System.out.println("## https://doi.org/10.1093/bioinformatics/btz068");
        System.out.println("## benjamin/dot/linard/at/lirmm/dot/fr");
        System.out.println("################################################");
        System.out.print(
        "\n--------------------------------------------------------------------\n" +
        " Minimum usage:\n\n"+
        " 1. For building the phylo-kmers database:\n"+
        "    java -jar RAPPAS.jar -p b -s [nucl|amino] -b ARbinary \n" +
        "    -w workdir -r alignment.fasta -t tree.newick\n" +
        "   \n" + 
        " 2. For placing sequences, using the database (DB) built in step 1:\n"+
        "    java -jar RAPPAS.jar -p p -d DB.union -q queries.fa \n"+ 
        "    \n"+ 
        " Note: For large references or high values of k, allocate more RAM :\n" +
        "       ex: java -Xms1024m -Xmx16g -jar RAPPAS.jar [options] \n"+
        "       -Xms -> memory allocated at startup. (m=MegaByte, g=GigaByte)\n"+
        "       -Xmx -> maximum memory allocated to the process.  \n"+
        "---------------------------------------------------------------------\n"+
        "\n" +
        "Main options:     Default values are in [].\n" +
        "---------------------------------------------------------------------\n"+
        "-b (--arbinary)   [file] Binary for marginal AR, currently 'phyml' and \n" +
        "                  'baseml' (from PAML) are supported. (b phase)\n" +
        "-d (--database)   [file] The database of ancestral kmers. (b|p phase) \n"+
        "-p (--phase)      [b|p] One of 'b' for \"Build\" or 'p' for \"Place\"\n" +
        "                  b: Build DB of phylo-kmers (done 1 time). \n" +
        "                  p: Phylogenetic placement itself (done n times)\n"+
        "                     requires the DB generated during 'build' phase.\n" +
        "-r (--refalign)   [file] Reference alignment in fasta format.\n" +
        "                  It must be the multiple alignment from which was \n" +
        "                  built the reference tree loaded with -t. (b phase) \n"+        
        "-s (--states)     ['nucl'|'amino'] States used in analysis. (b|p phase) \n" +    
        "-t (--reftree)    [file] Reference tree, in newick format.\n"+
        "-q (--queries)    [file[,file,...]] Fasta queries to place on the tree.\n" +
        "                  Can be a list of files separated by ','. (b|p phase)\n"+
        "-v (--verbosity)  [0] Verbosity level: -1=none ; 0=default ; 1=high\n" +  
        "-w (--workdir)    [path] Working directory for temp files. (b|p phase)\n" +  
        "\n" +
        "Outputs options:  Jplace, log files...  \n" +
        "---------------------------------------------------------------------\n"+
        "--keep-at-most    [7] Max number of placement per query kept in \n" +
        "                  the jplace output. (p phase)\n" +
        "--keep-factor     [0.01] Report placement with likelihood_ratio higher\n" +
        "                  than (factor x best_likelihood_ratio). (p phase)\n" +      
        "--write-reduction [file] Write reduced alignment to file. (b phase)\n" +
        "--guppy-compat    [] Ensures output is Guppy compatible. (p phase)\n" +
        "--dbfilename      [string] Set DB filename. (b phase)\n" +
        "\n" +
        "Algo options:     Use only if you know what you are doing...    \n" +
        "---------------------------------------------------------------------\n"+
        "-a (--alpha)      [1.0] Gammma shape parameter used in AR . (b phase)\n" +     
        "-c (--categories) [4] # categories used in AR . (b phase)\n" +   
        "-g (--ghosts)     [1] # ghost nodes injected per branches. (b phase)\n"+
        "-k (--k)          [8] k-mer length used at DB build. (b mode)\n" +   
        "-m (--model)      [GTR|LG] Model used in AR, one of the following:\n" +   
        "                  nucl  : JC69, HKY85, K80, F81, TN93, GTR \n" +  
        "                  amino : LG, WAG, JTT, Dayhoff, DCMut, CpREV,\n" +
        "                          mMtREV, MtMam, MtArt \n" +  
        "--arparameters    [string] Parameters passed to the software used for\n" +
        "                  anc. seq. reconstuct. Overrides -a,-c,-m options.\n" +
        "                  Value must be quoted by ' or \". Do not set options\n" +
        "                  -i,-u,--ancestral (managed by RAPPAS). (b phase)\n" +
        "                  PhyML example: \"-m HIVw -c 10 -f m -v 0.0 --r_seed 1\"\n" +     
        "--convertUO       [] U,O amino acids are converted to C,L. (b|p phase)\n"+        
        "--force-root      [] Root input tree (if unrooted) by adding a root\n"+
        "                  node on righmost branch of the trifurcation.(b phase)\n" +
        "--gap-jump-thresh [0.3] Gap ratio above which gap jumps are activated.\n" +
        "--no-reduction    [] Do not operate alignment reduction. This will \n" +
        "                  keep all sites of input reference alignment and \n" +
        "                  may produce erroneous ancestral k-mers. (b phase)\n" +
        "--ratio-reduction [0.99] Ratio for alignment reduction, e.g. sites \n" +
        "                  holding >99% gaps are ignored. (b phase)\n" +
        "--omega           [1.0] Modifier levelling the threshold used during\n"+
        "                  phylo-kmer filtering, T=(omega/#states)^k .(b phase)\n" +
        "--use_unrooted    [] Confirms you accept to use an unrooted reference\n"+
        "                  tree (option -t). The trifurcation described by the\n"+
        "                  newick file will be considered as root. Be aware that\n" +
        "                  meaningless roots may impact accuracy. (b phase)\n" +
        "\n" +
        "Debug options:    Use only if you know what you are doing...    \n" +
        "---------------------------------------------------------------------\n"+
        "--ambwithmax      [] Treat ambiguities with max, not mean. (p phase)\n" +
        "--ardir           [dir] Skip ancestral sequence reconstruction, and \n"+
        "                  uses outputs from the specified directory. (b phase)\n" +
        "--arinputonly     [] Generate only AR inputs. (b phase)\n" +
        "--aronly          [] Launch AR, but not DB build. (b phase)\n" +
        "--dbinram         [] Build DB, do not save it to a file, but directly\n" +
        "                     place queries given via -q instead.\n" +
        "--do-n-jumps      [] Shifts from 1 to n jumps. (b phase) \n" +
        "--force-gap-jump  [] Forces gap jump even if %gap<thresh. (b phase) \n" +
        "--jsondb          [] DB written as json. (careful, huge file outputs!)\n" +
        "--noamb           [] Do not treat ambiguous states. (p phase)\n" +
        "--threads         [4] #threads used in AR (if raxml-ng). (b phase)\n" +
        "\n" +
        "Final note:\n" +
        "---------------------------------------------------------------------\n"+
        "When you use this software, please cite RAPPAS and the binary used for \n" +
        "ancestral reconstruction, e.g. one of :\n" +
        " * phyml: Oliva et al, 2019. doi: 10.1093/bioinformatics/btz249\n" +
        " * paml: Yang Z, 2007. doi: 10.1093/molbev/msm088\n" +
        " * raxml-ng: Kzlov et al, 2019. doi: 10.1093/bioinformatics/btz305\n" +

        );
       System.exit(0);
    }  

    @Override
    public String toString() {
        return argsMap.toString();
    }
    
    
    
    
}
