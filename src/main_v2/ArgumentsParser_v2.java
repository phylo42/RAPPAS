/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package main_v2;

import main.*;
import etc.Environement;
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
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
    
    private HashMap<Integer,String> argsMap=null;
    
    //general parameters
    public int mode=DBBUILD_MODE;
    public File workingDir=null;//current directory by default, see below
    int verbose=0; //default, no verbosity

    
    //parameters for DB build
    public int k=8; //default=8
    public float alpha=1.5f; //default=1.5
    public int fakeBranchAmount=1;  //default =1
    public File alignmentFile=null;
    public File treeFile=null;
    //eventual directories passed for debugging
    public File ARBinary=new File("phyml"); //default = phyml command line
    public File ARDirToUse=null;
    public File exTreeDir=null;
    
    //parameters for placement
    public int minOverlap=100; //default =100
    public File queriesFile=null;
    public File databaseFile=null;
    public int dbsize=DB_FULL;
    
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
                        System.out.println("Cannot read from / write to this directory: "+wd.getAbsolutePath());
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
                        if (exTreeDir.isDirectory()&& exTreeDir.canRead()) {
                            this.exTreeDir=exTreeDir;
                        } else {
                            System.out.println("Cannot open directory given through option --extree: Not a directory or no read permission.");
                            System.exit(1);
                        }
                    }
                }
                
                //use defaults if -k,-a,-f not set
                if (!kGiven) {System.out.println("Default k="+this.k+" will be used.");}
                if (!alphaGiven) {System.out.println("Default alpha="+this.alpha+" will be used.");}
                if (!fGiven) {System.out.println("Default fakeNodePerBranch="+this.fakeBranchAmount+" will be used.");}
                
                break;
                
            //check if -q given when mode=p    
            case PLACEMENT_MODE:
                
                boolean dGiven=false;
                boolean qGiven=false;
                
                for (Iterator<Integer> it=argsMap.keySet().iterator();it.hasNext();) {
                    int index=it.next();
                    //test -q parameter
                    if (argsMap.get(index).equals("--queries") || argsMap.get(index).equals("-q")) {
                        File queries=new File(argsMap.get(index+1));
                        if (queries.isFile() && queries.canRead()) {
                            this.queriesFile=queries;
                            qGiven=true;
                        } else {
                            System.out.println(queries.getAbsolutePath());
                            System.out.println("Cannot open queries: Not a file or no read permission.");
                            System.exit(1);
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
                    //test -s parameter
                    if (argsMap.get(index).equals("--dbsize") || argsMap.get(index).equals("-s")) {
                        if (argsMap.get(index+1).equalsIgnoreCase("full")) {
                            this.dbsize=DB_FULL;
                        } else if (argsMap.get(index+1).equalsIgnoreCase("medium")) {
                            this.dbsize=DB_MEDIUM;
                        } else if (argsMap.get(index+1).equalsIgnoreCase("small")) {
                            this.dbsize=DB_SMALL;
                        }
                    }
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
        " 1.for building the initial database:\n"+
        "   java -jar viromplacer.jar -m B -i alignment.fasta -t tree.newick\n"+
        " 2.for placing the reads:\n"+
        "   java -jar viromplacer.jar -m P -q queries.fasta \n"+
        " Note: Do not hesistate to allocate more memory to the process.\n"+
        "       ex: java -jar -Xms8000m -Xmx16000m viromplacer.jar [...] \n"+
        "       Xms -> memory allocated at start ; Xmx maximum memory allocated  \n"+
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"+
        "-a (--alpha)      Alpha modificator levelling the score threshold. \n"+
        "                  ancestral reconstruction and DB build (B mode only).\n" +
        "-d (--database)   The database of ancestral words (P mode only). \n"+
        "-f (--fakebranch) # of fake branches to add on each original branch. \n"+
        "-i (--input)      Input sequences, in fasta format. When building the \n"+
        "                  database (-d mode), it is the multiple alignment from\n"+
        "                  which was inferred the phylogenetic tree. \n"+
        "-k (--k)          Word length used for the DB build (B mode only).\n" +
        "-l (--minoverlap) Option not implemented yet.\n" +
        "-m (--mode)       One of 'B' for \"Build\" or 'P' for \"Place\"\n" +
        "                   * B: Build DB of ancestral words and associated \n"+
        "                        probabilites. \n" +
        "                   * P: Placement of query sequences, using a DB of\n"+
        "                        ancestral words prevously built with mode B.\n" +
        "-s (--dbsize)     DB size to load for the placement (P mode only).\n" +    
        "                  One of [full|medium|small] \n" +    
        "-t (--tree)       Reference tree, in newick format. Used for ancestral \n"+
        "                  reconstruction and DB build (B mode only).\n" +
        "-q (--queries)    Sequences to place, in fasta format (P mode only)\n" +
        "-v (--verbose)    Verbosity level: 0=none ; 1=basic ; 2=full\n" +
        "-w (--workdir)    Path of the working directory (default= current dir).\n\n" +
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"+
        "Debug options: Use them only if you know what you are doing...    \n" +
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"+
        "--ardir           Skip AR reconstruction, searching its results in the\n"+
        "                  specified directory (B mode only).\n" +
        "--arbinary        Binary used for AR, ex: phyml. (B mode only).\n" +
        "--extree          Skip fake nodes injection, and use files present in\n"+
        "                  the specified directory instead (B mode only).\n" +
        "\n\n"
        );
       System.exit(1);
    }  

    @Override
    public String toString() {
        return argsMap.toString();
    }
    
    
    
    
}
