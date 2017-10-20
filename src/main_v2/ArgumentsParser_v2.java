/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main_v2;

import main.*;
import etc.Environement;
import java.io.File;
import java.nio.file.Files;
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

    public static final int DBBUILD_MODE=0;
    public static final int PLACEMENT_MODE=1;
    
    public static final int DB_FULL=0;
    public static final int DB_MEDIUM=0;
    public static final int DB_SMALL=0;
    
    public static final int TYPE_DNA=1;
    public static final int TYPE_PROTEIN=2;

    
    private HashMap<Integer,String> argsMap=null;
    
    //general parameters
    public int mode=DBBUILD_MODE;
    public File workingDir=null;//current directory by default, see below
    int verbose=0; //default, no verbosity
    int analysisType=TYPE_DNA;
    
    //parameters for alignment reduction
    public boolean reduction=true;
    public File reducedAlignFile=null;
    double reductionRatio=0.999;
    
    //parameters for DB build
    public int k=8; //default=8
    public float alpha=1.4f; //default=1.5
    public int fakeBranchAmount=1;  //default =1
    public File alignmentFile=null;
    public File treeFile=null;
    public boolean builddbfull=false; //default=false, as dbfull is quite useless with the current algo
    public boolean forceRooting=false;
    public boolean noCalibration=false;
    public boolean onlyFakeNodes=false;
    //eventual directories passed for debugging
    public File ARBinary=new File("phyml"); //default = phyml command line
    public File ARDirToUse=null;
    public File exTreeDir=null;
    public boolean dbInRAM=false;
    public boolean unionHash=true;
    
    //parameters for placement
    public int minOverlap=100; //default =100
    public List<File> queriesFiles=null;
    public File databaseFile=null;
    public int dbsize=DB_MEDIUM;
    public Float nsBound=null;
    public int keepAtMost=7;
    public float keepFactor=0.1f;
    
    //call string
    public String callString=null;
    
    public ArgumentsParser_v2(String[] args) {
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
        
        //if no mode, bad
        if ((!argsMap.containsValue("-m")) && (!argsMap.containsValue("--mode")) ) {
            System.out.println("Cannot find 'mode' (option -m). See help (-h) for details.");
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

        //general loop, must pass a 1st time on all arguments to set the mode for sure
        boolean wGiven=false;
        for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
            int index=it.next();
            
            //check if --mode associated to correct minimum values
            if (argsMap.get(index).equals("-m") || argsMap.get(index).equals("--mode")) {
                if (argsMap.get(index+1).equalsIgnoreCase("b")) {this.mode=DBBUILD_MODE;}
                if (argsMap.get(index+1).equalsIgnoreCase("p")) {this.mode=PLACEMENT_MODE;}
            }
            
            //check verbosity
            if (argsMap.get(index).equals("-v") || argsMap.get(index).equals("--verbose")) {
                try {
                    this.verbose=Integer.parseInt(argsMap.get(index+1));
                    if (verbose>0) { //for now, only one verbosity level
                        System.setProperty("debug.verbose", "1");
                    }
                } catch (NumberFormatException ex) {
                    System.out.println("Cannot parse '-v (--verbose)' as an integer value.");
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
            //test -s parameter
            if (argsMap.get(index).equals("--states") || argsMap.get(index).equals("-s")) {
                if (argsMap.get(index+1).equalsIgnoreCase("nucl")) {
                    this.analysisType=TYPE_DNA;
                } else if (argsMap.get(index+1).equalsIgnoreCase("prot")) {
                    this.analysisType=TYPE_PROTEIN;
                } else {
                    System.out.println("Unexpected -s (--states) value, must be one of [nucl|prot].");
                    System.exit(1);
                }
            }
            
            
            
        }
        if ( !(mode==DBBUILD_MODE) && !(mode==PLACEMENT_MODE) ) {
            System.out.println("Unexpected -m (--mode) value, must be 'b' (build_db) or 'p' (place) !");
            System.exit(1);
        }
        
        if (!wGiven) {
            //this.workingDir=Environement.getExecutablePathWithoutFilename(this.getClass());
            this.workingDir=Environement.getCurrentDirectory().toFile();
            System.out.println("Default working directory (current directory) will be used.");
            System.out.println("workDir="+this.workingDir.getAbsolutePath());
        }
        
        
        
        ////from here, we know in which mode the program will be launched///////
        ////////////////////////////////////////////////////////////////////////
        //check if -a ,-t, -k given when mode=b
        switch (mode) {
            
            case DBBUILD_MODE:
                
                boolean kGiven=false;
                boolean alphaGiven=false;
                boolean fGiven=false;
                
                for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
                    int index=it.next();
                    //test -a parameter
                    if (argsMap.get(index).equals("--input") || argsMap.get(index).equals("-i")) {
                        File alignment=new File(argsMap.get(index+1));
                        if (alignment.isFile() && alignment.canRead()) {
                            this.alignmentFile=alignment;
                        } else {
                            System.out.println(alignment.getAbsolutePath());
                            System.out.println("Cannot open alignment: Not a file or no read permission.");
                            System.exit(1);
                        }
                    }
                    //test -t parameter
                    if (argsMap.get(index).equals("--tree") || argsMap.get(index).equals("-t")) {
                        File tree=new File(argsMap.get(index+1));
                        if (tree.isFile() && tree.canRead()) {
                            this.treeFile=tree;
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
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-k (--k)' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    //test -a parameter
                    if (argsMap.get(index).equals("--alpha") || argsMap.get(index).equals("-a")) {
                        String alphaVal=argsMap.get(index+1);
                        try {
                            this.alpha=Float.parseFloat(alphaVal);
                            alphaGiven=true;
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-a (--alpha)' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    
                    
                    //test -f parameter
                    if (argsMap.get(index).equals("--fakebranch") || argsMap.get(index).equals("-f")) {
                        String fVal=argsMap.get(index+1);
                        try {
                            this.fakeBranchAmount=Integer.parseInt(fVal);
                            fGiven=true;
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '-f (--fakebranch)' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    
                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS
                    
                    //test --skipredu
                    if (argsMap.get(index).equals("--skip-redu")) {
                        System.out.println("Original alignment will be used (no gapped columns removed).");
                        this.reduction=false;
                    }
                    
                    //test --writeredu
                    if (argsMap.get(index).equals("--write-redu")) {
                        File f=new File(argsMap.get(index+1));
                        if (f.isFile() && f.canWrite()) {
                            this.reducedAlignFile=f;
                        } else {
                            System.out.println("Cannot open directory given through option --write-reduced-alignment: Not a file or no write permission.");
                            System.exit(1);
                        }
                    }
                    
                    //test --ratioredu
                    if (argsMap.get(index).equals("--ratioredu")) {
                        String alphaVal=argsMap.get(index+1);
                        try {
                            this.reductionRatio=Double.parseDouble(alphaVal);
                            alphaGiven=true;
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--ratioredu' as a double value.");
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
                    
                    //test --arbinary parameter
                    if (argsMap.get(index).equals("--arbinary")) {
                        File ARBinaryFile=new File(argsMap.get(index+1));
                        System.out.println("Using AR binary provided by user: "+ARBinaryFile.getAbsolutePath());
                        if (ARBinaryFile.isFile() && ARBinaryFile.canExecute()) {
                            this.ARBinary=ARBinaryFile;
                        } else {
                            System.out.println("Cannot execute binary loaded through option --arbinary: Not a file or no execution permission.");
                            System.exit(1);
                        }
                    }
                    
                    //test --extendedtree parameter
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
                    
                    //test --builddbfull parameter
                    if (argsMap.get(index).equals("--dbfull")) {
                        this.builddbfull=true;
                    }
                    
                    //test --froot parameter
                    if (argsMap.get(index).equals("--froot")) {
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
                    if (argsMap.get(index).equals("--nocalib")) {
                        this.noCalibration=true;
                    }
                    
                    //test --unihash parameter
                    if (argsMap.get(index).equals("--poshash")) {
                        this.unionHash=false;
                    }
                    
                    //test --onlyfake
                    if (argsMap.get(index).equals("--onlyfakes")) {
                        System.out.println("DB will only contain ancestral kmers associated to fake nodes.");
                        this.onlyFakeNodes=true;
                    }
                    
                    //test --keep-at-most
                    if (argsMap.get(index).equals("--keep-at-most")) {
                        String val=argsMap.get(index+1);
                        try {
                            this.keepAtMost=Integer.parseInt(val);
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
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--keep-factor' as an integer value.");
                            System.exit(1);
                        }
                        
                    }
                    
                    
                    //////////////////////////////////////
                    //////////////////////////////////////
                    //DEBUG OPTIONS END HERE
                    
                    
                }
                
                //use defaults if -k,-a,-f not set
                if (!kGiven) {System.out.println("Default k="+this.k+" will be used.");}
                if (!alphaGiven) {System.out.println("Default alpha="+this.alpha+" will be used.");}
                if (!fGiven) {
                    System.out.println("Default injectionPerBranch="+this.fakeBranchAmount+" will be used.");
                }
                
                break;
                
            //check if -q given when mode=p    
            case PLACEMENT_MODE:
                
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
                        } catch (NumberFormatException ex ) {
                            System.out.println("Cannot parse '--keep-factor' as an integer value.");
                            System.exit(1);
                        }
                        
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
        System.out.print(
        "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" +
        " Minimum usage,\n"+
        " 1.For building the ancestral kmers database:\n"+
        "   java -jar viromplacer.jar -m B -b ARbinary -w workdir \\ \n"+
        "   -s dna -i alignment.fasta -t tree.newick\n" +
        " 2.For placing the query reads, using the database built in 1. :\n"+
        "   java -jar viromplacer.jar -m P -q queries.fasta \n"+
                
        " Note: Do not hesistate to allocate lots of memory at the first step,\n" +
        "       as increasing k rapidly brings to large requirements.\n"+
        "       ex: java -jar -Xms1024m -Xmx16g viromplacer.jar [...] \n"+
        "       -Xms -> memory allocated at startup. (m=MegaByte, g=GigaByte)\n"+
        "       -Xmx -> maximum allocation allowed.  \n"+
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"+
        "-a (--alpha)      [float] Alpha modifier levelling the proba threshold \n"+
        "                  used in ancestral words filtering. (B mode only)\n" +
        "-b (--arbinary)   [file] Binary for marginal AR, currently phyml and \n" +
        "                  baseml (from PAML) are supported. (B mode only)\n" +
        "-d (--database)   [file] The database of ancestral kmers. (B/P mode) \n"+
        "-f (--fakebranch) [int] # fake nodes to add along ref tree branches. \n"+
        "-i (--input)      Input sequences, in fasta format. When building the \n"+
        "                  database (-d mode), it is the multiple alignment from\n"+
        "                  which was inferred the phylogenetic tree. \n"+
        "-k (--k)          [int] Word length used for the DB build. (B mode only)\n" +
        "-m (--mode)       One of 'B' for \"Build\" or 'P' for \"Place\"\n" +
        "                   * B: Build DB of ancestral words and associated \n"+
        "                        probabilites. \n" +
        "                   * P: Placement of query sequences, using a DB of\n"+
        "                        ancestral words prevously built with mode B.\n" +
        "-s (--states)     [nucl|prot] States used in analysis (B mode). \n" +    
        "-t (--tree)       [file] Reference tree, in newick format.\n"+
        "                  reconstruction and DB build (B mode only).\n" +
        "-q (--queries)    [file[,file,...]] Fasta queries to place on the tree." +
        "                  Can be a list of files separated by ','. (B/P mode)\n"+
        "                  be placed if filenames are separated by ','.\n" +
        "-v (--verbose)    Verbosity level: -1=null ; 0=basic ; 1=full\n" +
        "-w (--workdir)    [dir] Path to the working directory (B/P mode).\n\n" +
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"+
        "Debug options: Use only if you know what you are doing...    \n" +
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"+
        "--ardir           [dir] Skip AR reconstruction, searching its results\n"+
        "                  in the specified directory (B mode only)\n" +
        "--extree          [dir] Skip fake nodes injection, and use files present\n"+
        "                  in the specified directory instead (B mode only)\n" +
        "--dbfull          [] Save full DB (unused in algo). (B mode only)\n" +      
        "--froot           [] If input tree is unrooted, root it. (B mode only).\n" +
        "--nsbound         [float] Force normalized score bound. (P mode only)\n" +
        "--dbinram         [] Operate B mode, but whitout saving DB to files and\n" +
        "                  directly place queries given via -q .\n" +
        "--nocalib         [] No calibration, use threshold formula. (B mode only).\n" +
        "--poshash         [] Places using older deprecated hash. (B mode only)\n" +
        "--skip-redu       [] Do not operate alignment reduction. By default,\n" +
        "                  reference alignment columns with more than 'ratio'\n" +
        "                  gaps are ignored (see --ratio-redu). (B mode only)\n" +
        "--write-redu      [file] Write reduced alignment to file. (B mode only)\n" +
        "--ratio-redu      [float] Ratio for alignment reduction, default=0.999 \n" +
        "                  i.e. sites with >99.9% gaps are ignored. (B mode only)\n" +
        "--onlyfakes       [] Ancestral k-mers only for fake nodes. (B mode only)\n" +
        "--keep-at-most    [int] Max number of placement per query we keep in\n" +
        "                  the jplace output (default=7). (B/P mode)\n" +
        "--keep-factor     [float] Report placement with likelihood_ratio higher\n" +
        "                  than (factor x best_likelihood_ratio) (default=0.1).\n" +
        "                  (B/P mode)" +
        "\n\n"
        );
       System.exit(1);
    }  

    @Override
    public String toString() {
        return argsMap.toString();
    }
    
    
    
    
}
