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
import models.EvolModel;

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
    
    String version="";


    //RAPPAS general parameters
    public int phase=DBBUILD_PHASE;
    public File workingDir=null;//current directory by default, see below
    public int verbose=0; //default, no verbosity
    public int states=STATES_DNA;
    boolean convertUOX=false; //must be true so that U/O/X are converted in C/L/-
    
    //RAPPAS parameters for alignment reduction
    public boolean reduction=true;
    public File reducedAlignFile=null;
    public double reductionRatio=0.999;
    
    //RAPPAS parameters for DB build
    public int k=8; //default=8
    public float omega=1.0f; //default=1.0, quantity allowing to modulate the threshold
    public int ghostsAmount=1;  //default =1
    public File alignmentFile=null;
    public File treeFile=null;
    public boolean forceRooting=false;
    public String modelString=null; //if not set by options, default will be used (see Main)
    public float alpha=1.0f;
    public int categories=4;
    public String arparameters=null;

    //passed for debugging in DB_build
    public File ARBinary=new File("phyml"); //default = phyml command line
    public File ARDirToUse=null;  // used to load directly an already computed AR
    public File exTreeDir=null; //not functionnal anymore
    public boolean builddbfull=false; //default=false, as dbfull is quite useless with the current algo
    public boolean noCalibration=true; //skip calibration step
    public boolean dbInRAM=false; //do not write DB in a file and immediately places reads passed with -q
    public boolean unionHash=true; //if false, use old positionnal hash
    public boolean onlyFakeNodes=true; //if false, uses ancestral kmers of original nodes
    public boolean doGapJumps=false; //take gap jumps into account when building kmers
    public boolean limitTo1Jump=true; //only allow a 1st jump, not jump combinations
    public float gapJumpThreshold=0.3f; //gap jumps are activated if >30% gaps in the ref alignment
    
    //RAPPAS parameters for placement
    public int minOverlap=100; //used in entropy computation
    public List<File> queriesFiles=null;
    public File databaseFile=null;
    public Float nsBound=Float.NEGATIVE_INFINITY; //do not consider threshold calculated by formula but this particular threshold
    public int keepAtMost=7; //as in pplacer
    public float keepFactor=0.01f; //as in pplacer
    public boolean guppyCompatible=false;
    
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
            //test --convertUOX parameter
            if (argsMap.get(index).equals("--convertUOX")) {
                this.convertUOX=true;
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
            //this.workingDir=Environement.getExecutablePathWithoutFilename(this.getClass());
            this.workingDir=Environement.getCurrentDirectory().toFile();
            System.out.println("Default working directory (current directory) will be used.");
            System.out.println("workDir="+this.workingDir.getAbsolutePath());
        }
        
        if (!statesGiven) {
            System.out.println("Analysis states not found. Use option -s (--states) and one of 'nucl' or 'amino'.");
            System.exit(1);
        }

        
        
        ////from here, we know in which phase the program will be launched///////
        ////////////////////////////////////////////////////////////////////////
        //check if -a ,-t, -k given when phase=b
        switch (phase) {
            
            case DBBUILD_PHASE:
                
                boolean kGiven=false;
                boolean omegaGiven=false;
                boolean fGiven=false;
                boolean alignGiven=false;
                boolean treeGiven=false;
                
                for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
                    int index=it.next();
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
                            omegaGiven=true;
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
                        System.out.println("Using AR provided by user: "+ARDir.getAbsolutePath());
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

                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS END HERE
                    
                    
                }
                

                ///////////////////////////////////////////////
                //do not continue if align/tree not set (-r / -t)
                ///////////////////////////////////////////////
                if (!alignGiven) {System.out.println("Reference alignment not set correctly (use -r) ."); System.exit(1);}
                if (!treeGiven) {System.out.println("Reference tree not set correctly (use -t) ."); System.exit(1);}
                ///////////////////////////////////////////////
                //use defaults if -k,-a,-f not set
                ///////////////////////////////////////////////
                if (!kGiven) {System.out.println("Default k="+this.k+" will be used.");}
                if (!omegaGiven) {System.out.println("Default alpha="+this.omega+" will be used.");}
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
        System.out.println("## Linard B, Swenson KM, Pardi F");
        System.out.println("## LIRMM, Univ. of Montpellier, CNRS, France");
        System.out.println("## https://doi.org/10.1101/328740");
        System.out.println("## benjamin/dot/linard/at/lirmm/dot/fr");
        System.out.println("################################################");
        System.out.print(
        "\n--------------------------------------------------------------------\n" +
        " Minimum usage:\n\n"+
        " 1. For building the phylo-kmers database:\n"+
        "    java -jar viromplacer.jar -p b -s [nucl|prot] -b ARbinary \n" +
        "    -w workdir -s nucl -r alignment.fasta -t tree.newick\n" +
        "   \n" + 
        " 2. For placing the query reads, using the database built in 1. :\n"+
        "    java -jar viromplacer.jar -p p -d database.union -q queries.fasta \n"+ 
        "    \n"+ 
        " Note: For larger alignments or values of k, allocate more memory:\n" +
        "       ex: java -jar -Xms1024m -Xmx16g viromplacer.jar [...] \n"+
        "       -Xms -> memory allocated at startup. (m=MegaByte, g=GigaByte)\n"+
        "       -Xmx -> maximum allocation allowed.  \n"+
        "---------------------------------------------------------------------\n"+
        "\n" +
        "Main options:     Default values are in [].\n" +
        "---------------------------------------------------------------------\n"+
        "-b (--arbinary)   [file] Binary for marginal AR, currently 'phyml' and \n" +
        "                  'baseml' (from PAML) are supported. (b phase)\n" +
        "-d (--database)   [file] The database of ancestral kmers. (b|p phase) \n"+
        "-p (--phase)       ['b'|'p'] One of 'b' for \"Build\" or 'p' for \"Place\"\n" +
        "                   * b: Build DB of phylo-kmers (done 1 time). \n" +
        "                   * p: Phylogenetic placement itself (done n times)\n"+
        "                        requires the DB generated during 'build' phase.\n" +
        "-r (--refalign)   [file] Reference alignment in fasta format.\n" +
        "                  It must be the multiple alignment from which was \n" +
        "                  inferred the reference tree (option -t). (b phase) \n"+        
        "-s (--states)     ['nucl'|'amino'] States used in analysis. (b|p phase) \n" +    
        "-t (--reftree)    [file] Reference tree, in newick format.\n"+
        "                  reconstruction and DB build (b phase).\n" +
        "-q (--queries)    [file[,file,...]] Fasta queries to place on the tree.\n" +
        "                  Can be a list of files separated by ','. (b|p phase)\n"+
        "                  be placed if filenames are separated by ','.\n" +
        "-v (--verbosity)  [0] Verbosity level: -1=null ; 0=low ; 1=high\n" +  
        "-w (--workdir)    [.] Path to the working directory (b|p phase).\n" +  
        "\n" +
        "Outputs options:  Jplace, log files...  \n" +
        "---------------------------------------------------------------------\n"+
        "--keep-at-most    [7] Max number of placement per query kept in \n" +
        "                  the jplace output. (p phase)\n" +
        "--keep-factor     [0.01] Report placement with likelihood_ratio higher\n" +
        "                  than (factor x best_likelihood_ratio). (p phase)\n" +      
        "--write-reduction [file] Write reduced alignment to file. (b phase)\n" +
        "--guppy-compat    [] Ensures output is Guppy compatible. (p phase)\n" +
        "\n" +
        "Algo options:     Use only if you know what you are doing...    \n" +
        "---------------------------------------------------------------------\n"+
        "-a (--alpha)      [1.0] Shape parameter used in AR . (b phase)\n" +     
        "-c (--categories) [4] # categories used in AR . (b phase)\n" +   
        "-g (--ghosts)     [1] # ghost nodes injected per branches. (b phase)\n"+
        "-k (--k)          [8] k-mer length used at DB build. (b mode)\n" +   
        "-m (--model)      [GTR|LG] Model used in AR, one of the following:\n" +   
        "                  *nucl  : JC69, HKY85, K80, F81, TN93, GTR \n" +  
        "                  *amino : LG, WAG, JTT, Dayhoff, DCMut, CpREV,\n" +
        "                           mMtREV, MtMam, MtArt \n" +  
        "--convertUOX      [] U,O,X amino acids become C,L,- (b|p phase).\n"+        
        "--force-root      [] Root input tree if non rooted. (b phase)\n" +
        "--ratio-reduction [0.999] Ratio for alignment reduction, e.g. sites \n" +
        "                  holding >99.9% gaps are ignored. (b phase)\n" +
        "--arparameters    [string] Parameters passed to the software used for\n" +
        "                  anc. seq. reconstuct. Overrides -a,-c,-m options.\n" +
        "                  Value must be quoted by ' or \". Do not set options\n" +
        "                  -i,-u,--ancestral (managed by RAPPAS). (b phase)\n" +
        "                  PhyML example: \"-m HIVw -c 10 -f m -v 0.0 --r_seed 1\"\n" +
        "--no-reduction    [] Do not operate alignment reduction. This will \n" +
        "                  keep all sites of input reference alignment and \n" +
        "                  may produce erroneous ancestral k-mers. (b phase)\n" +
        "--gap-jump-thresh [0.3] Gap ratio above which gap jumps are activated.\n" +
        "--omega           [1.0] Modifier levelling the threshold used during\n"+
        "                  phylo-kmers filtering, T=(omega/#states)^k .(b phase)\n" +
        "\n" +
        "Debug options:    Use only if you know what you are doing...    \n" +
        "---------------------------------------------------------------------\n"+
        "--ardir           [dir] Skip ancestral sequence reconstruction, and \n"+
        "                  uses outputs from the specified directory. (b phase)\n" +
        "--extree          [dir] Skip phantom nodes injection, and use already\n"+
        "                  injected trees from the specified directory.(b phase)\n" +
        "--nsbound         [float] Force normalized score bound. (p phase)\n" +
        "--dbinram         [] Operate B mode, but whitout saving DB to files and\n" +
        "                  directly place queries given via -q .\n" +
        "--calibration     [] Prototype calib. on random anc. kmers. (b phase).\n" +
        "--do-n-jumps      [] Shifts from 1 to n jumps. (b phase) \n" +
        "--force-gap-jump  [] Forces gap jump even if %gap<thresh. (b phase) \n" +
        "\n"
        );
       System.exit(0);
    }  

    @Override
    public String toString() {
        return argsMap.toString();
    }
    
    
    
    
}
