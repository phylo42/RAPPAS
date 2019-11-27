/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package inputs;

import core.AAStates;
import core.States;
import etc.Infos;
import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import models.EvolModel;

/**
 * class used to setup the maginal ancerstral reconstruction (AR),
 using an external ARBinary. Currently PAML and PHYML are supported.
 * @author ben
 */
public class ARProcessLauncher {
    
    public static final int AR_PAML=1;
    public static final int AR_PHYML=2;
    public static final int AR_RAXMLNG=3;
    
    public int currentProg=AR_PHYML;
    public File ARBinary=null;
    public File ARPath=null;
    public File alignPath=null;
    public File treePath=null;
    
    private boolean verboseAR=true;
    private States s=null;
    private EvolModel model=null;
    private int threads;
    private String ARParameters;
    private String phymlVersion=null;
    private boolean phymlAcceptsDuplicates=true;
    
    /**
     * prepare the marginal AR
     * @param ARBinary
     * @param verboseAR
     * @param s
     * @param model
     * @param ARParameters
     * @param threads
     */
    public ARProcessLauncher(File ARBinary, boolean verboseAR, States s, EvolModel model, String ARParameters,int threads) {
        this.ARBinary=ARBinary;
        this.verboseAR=verboseAR;
        this.s=s;
        this.ARParameters = ARParameters;
        this.model=model;
        this.threads = threads;
        //test if ARBinary is something supported
        if (ARBinary.getName().contains("phyml")) {
            System.out.println("I guess, from binary name, we are using PhyML (phyml).");
            this.currentProg=AR_PHYML;
        } else if(ARBinary.getName().contains("baseml")) {
            System.out.println("I guess, from binary name, we are using PAML (baseml).");
            //check tha baseml or codeml is called according to states
            if (s instanceof AAStates) {
                System.out.println("You set 'baseml' as AR binary, but states are amino acids. Please use 'codeml' binary instead.");
                System.exit(1);
            }
            this.currentProg=AR_PAML;
        } else if (ARBinary.getName().contains("codeml")) {
            System.out.println("I guess, from binary name, we are using PAML (codeml).");
            //check tha baseml or codeml is called according to states
            if (!(s instanceof AAStates)) {
                System.out.println("You set 'codeml' as AR binary, but states are nucleotides. Please use 'baseml' binary instead.");
                System.exit(1);
            }
            this.currentProg=AR_PAML;
        } else if (ARBinary.getName().contains("raxml-ng")) {
            System.out.println("I guess, from binary name, we are using RAXML-NG.");
            this.currentProg=AR_RAXMLNG;
        } else {
                System.out.println("AR binary is unknown, currently RAPPAS support only phyml, raxml-ng & baseml+codeml (from paml package).");
            System.exit(1);
        }

    }
    
    
    
    
    
    /**
     * execute immediately marginal AR
     * @param ARPath path in which the external AR program will work
     * @param alignPath alignment in Phylip format, used for the marginal AR
     * @param treePath tree in Newick format, used for the marginal AR
     */
    public void launchAR(File ARPath, File alignPath, File treePath) {
        this.ARPath=ARPath;
        this.alignPath=alignPath;
        this.treePath=treePath;
        
        if (ARBinary==null) {
            System.out.println("Path to executable used for ancestral reconstruction is not set correctly.");
            System.exit(1);
        }
        if ( (!ARBinary.isFile()) || (!ARBinary.canExecute()) ) {
            System.out.println("The set AR binary is not a file or do not have execution rights.");
            System.out.println("AR binary: "+ARBinary.getAbsolutePath());
            System.exit(1);
        }
        
        switch (currentProg) {
            case AR_PAML:
                Infos.println("PAML AR was selected.");
                launchPAML();
                break;
            case AR_PHYML:
                Infos.println("PHYML AR was selected.");
                {
                    try {
                        testPHYMLVersion();
                    } catch (IOException ex) {
                        Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                launchPHYML();
                break;
            case AR_RAXMLNG:
                Infos.println("RAXML-NG AR was selected.");
                launchRAXMLNG();
                break;
            default:
                break;
        }
    }
    /**
     * load already existing marginal AR, will verify if files found in 
     * ARPath are compatible with the AR_PROG set at instantiation 
     * @param ARPath path in which the external AR program will work
     * @param alignPath alignment in Phylip format, used for the marginal AR
     * @param treePath tree in Newick format, used for the marginal AR
     */
    public void loadExistingAR(File ARPath, File alignPath, File treePath) {
        this.ARPath=ARPath;
        this.alignPath=alignPath;
        this.treePath=treePath;
        switch (currentProg) {
            case AR_PAML:
                //we expect the rst file to be present in ARPath
                File rst=new File(ARPath.getAbsolutePath()+File.separator+"rst");
                if (!rst.exists() || !rst.canRead()) {
                    System.out.println("rst file do not exists or cannot be read in "+rst.getAbsolutePath());
                    System.exit(1);
                }
                break;
            case AR_PHYML:
                //we expect 2 files to be present in ARPath
                //alignName_phyml_stats.txt
                //alignName_phyml_tree.txt
                File stats=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName()+"_phyml_stats.txt");
                File tree=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName()+"_phyml_ancestral_tree.txt");
                File seq=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName()+"_phyml_ancestral_seq.txt");
                if (!stats.exists() || !stats.canRead()) {
                    System.out.println(stats.getAbsolutePath()+" do not exists or cannot be read.");
                    System.exit(1);
                }
                if (!tree.exists() || !tree.canRead()) {
                    System.out.println(tree.getAbsolutePath()+" do not exists or cannot be read.");
                    System.exit(1);
                }
                if (!seq.exists() || !seq.canRead()) {
                    System.out.println(seq.getAbsolutePath()+" do not exists or cannot be read.");
                    System.exit(1);
                }
                break;
            case AR_RAXMLNG:
                //we expect 2 files to be present in ARPath
                //alignName.raxml.ancestralProbs
                //alignName.raxml.ancestralTree
                File statsraxml=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName()+".raxml.log");
                File treeraxml=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName()+".raxml.ancestralTree");
                File probaraxml=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName()+".raxml.ancestralProbs");
                if (!statsraxml.exists() || !statsraxml.canRead()) {
                    System.out.println(statsraxml.getAbsolutePath()+" do not exists or cannot be read.");
                    System.exit(1);
                }
                if (!treeraxml.exists() || !treeraxml.canRead()) {
                    System.out.println(treeraxml.getAbsolutePath()+" do not exists or cannot be read.");
                    System.exit(1);
                }
                if (!probaraxml.exists() || !probaraxml.canRead()) {
                    System.out.println(probaraxml.getAbsolutePath()+" do not exists or cannot be read.");
                    System.exit(1);
                }
                break;
        }
    }
    
    /**
     * do not execute marginal AR, but prepare its config files (if required,
     * like in PAML) and directories and returns he command to execute to launch
     * the AR from a console.
     * @param ARPath path in which the external AR program will work
     * @param alignPath alignment in Phylip format, used for the marginal AR
     * @param treePath tree in Newick format, used for the marginal AR
     */
    public String prepareAR(File ARPath, File alignPath, File treePath) {
        this.ARPath=ARPath;
        this.alignPath=alignPath;
        this.treePath=treePath;
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }
        StringBuilder sb=new StringBuilder();
        List<String> com=null;
        switch (this.currentProg) {
            case AR_PAML:
                Infos.println("PAML AR was selected.");
                //build PHYML program command-line for the AR, without execution but built its .ctl file.
                FileWriter fw=null;
                StringBuilder sb2=null;
                try {
                    //paml launch command itself
                    //launch paml externally to build the posterior probas on the extended tree
                    com=buildPAMLCommand(fw);
                    //execution
                } catch (IOException ex) {
                    Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        fw.close();
                        return sb2.toString();
                    } catch (IOException ex) {
                        Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                break;
            case AR_PHYML:
                Infos.println("PHYML AR was selected.");
                {
                    try {
                        testPHYMLVersion();
                    } catch (IOException ex) {
                        Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                com=buildPhyMLCommand();
                break;
            case AR_RAXMLNG:
                Infos.println("RAXML-NG AR was selected.");
                com=buildRAXMLNGCommand();
                break;
            default:
                break;
        }
        for (int i = 0; i < com.size(); i++) {
            sb.append(" ");
            sb.append(com.get(i));
        }
        Infos.println("Ancestral reconstruct command: "+sb.toString());
        return sb.toString();
    }

    /**
     * execute RAXML-NG program for the AR.
     * important: input tree should NOT have node labels !
     */
    private void launchPHYML() {
        if ((!ARPath.isDirectory()) || (!ARPath.canWrite())) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }

        try {
            List<String> com = buildPhyMLCommand();
            Infos.println("Ancestral reconstruct command: " + com);
            //execution
            executeProcess(com);
            //phyml is written all data files near the input alignment file...
            //we move them to the AR directory
            //files are:
            // 1. alignName_phyml_ancestral_seq.txt         (used)
            // 2. alignName_phyml_stats.txt                 (unused)
            // 3. alignName_phyml_ancestral_tree.txt        (used)
            // 4. alignName_phyml_tree.txt                  (unused)
            File stats = new File(alignPath.getAbsolutePath() + "_phyml_stats.txt");
            File tree = new File(alignPath.getAbsolutePath() + "_phyml_ancestral_tree.txt");
            File seq = new File(alignPath.getAbsolutePath() + "_phyml_ancestral_seq.txt");
            File oriTree = new File(alignPath.getAbsolutePath() + "_phyml_tree.txt");
            //check that they were correctly created
            if (!stats.exists() || !tree.exists() || !seq.exists() || !oriTree.exists()) {
                System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                System.out.println("!!! Phyml outputs are missing, the process may have failed... !!!");
                System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                System.out.println("Some clues may be found in the AR/AR_sdterr.txt  or AR/AR_sdtour.txt files.\n");
                System.out.println("A common reason: phyml limits are set during compilation (defaults are <4000 taxa, label<1000 characters...)");
                System.out.println("You can increase these limits in the src/utilities.h source file and recompile phyml.");
                System.out.println("Changing these lines to these new values should be enough in most cases: ");
                System.out.println("  #define  T_MAX_FILE          2000");
                System.out.println("  #define  T_MAX_NAME          5000");
                System.out.println("  #define  N_MAX_OTU         262144");
                System.exit(1);
            }
            File statsNew = new File(alignPath.getParent().replace("/extended_trees", "/AR") + File.separator + alignPath.getName() + "_phyml_stats.txt");
            File treeNew = new File(alignPath.getParent().replace("/extended_trees", "/AR") + File.separator + alignPath.getName() + "_phyml_ancestral_tree.txt");
            File seqNew = new File(alignPath.getParent().replace("/extended_trees", "/AR") + File.separator + alignPath.getName() + "_phyml_ancestral_seq.txt");
            File oriTreeNew = new File(alignPath.getParent().replace("/extended_trees", "/AR") + File.separator + alignPath.getName() + "_phyml_tree.txt");

            boolean move = stats.renameTo(statsNew);
            boolean move2 = tree.renameTo(treeNew);
            boolean move3 = seq.renameTo(seqNew);
            boolean move4 = oriTree.renameTo(oriTreeNew);

            if (!move || !move2 || !move3 || !move4) {
                System.out.println("Cannot move phyml results from /extended_tree to /AR directory");
                System.exit(1);
            }

        } catch (IOException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    
    /**
     * execute RAXMLNG program for the AR.
     */
    private void launchRAXMLNG() {
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }
        
        try {
            List<String> com=buildRAXMLNGCommand();
            Infos.println("Ancestral reconstruct command: "+com);
            //execution
            executeProcess(com);
            //raxml-ng writes all data files in same directory than the alignment file...
            //we move results to the AR directory
            //files are:
            // 1. alignName.raxml.ancestralProbs             (used)
            // 2. alignName.raxml.log                        (unused)
            // 3. alignName.raxml.ancestralTree              (used)
            // 4. alignName.raxml.startTree                  (unused)
            // 5. alignName.raxml.ancestralStates            (unused)
            // 6. alignName.raxml.rba                        (unused)

            File stats=new File(alignPath.getAbsolutePath()+".raxml.log");
            File tree=new File(alignPath.getAbsolutePath()+".raxml.ancestralTree");
            File seq=new File(alignPath.getAbsolutePath()+".raxml.ancestralProbs");
            File oriTree=new File(alignPath.getAbsolutePath()+".raxml.startTree");
            File as=new File(alignPath.getAbsolutePath()+".raxml.ancestralStates");
            File rba=new File(alignPath.getAbsolutePath()+".raxml.rba");

            //check that they were correctly created
            if (!stats.exists() || !tree.exists() || !seq.exists() || !oriTree.exists()) {
                System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                System.out.println("!!! Raxml-ng outputs are missing, the process may have failed... !!!");
                System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                System.out.println("Some clues may be found in the AR/AR_sdterr.txt  or AR/AR_sdtour.txt files.\n");
                System.exit(1);
            }
            File statsNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+".raxml.log");
            File treeNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+".raxml.ancestralTree");
            File seqNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+".raxml.ancestralProbs");
            File oriTreeNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+".raxml.startTree");
            File asNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+".raxml.ancestralStates");
            File rbaNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+".raxml.rba");
            
            boolean move=stats.renameTo(statsNew);
            boolean move2=tree.renameTo(treeNew);
            boolean move3=seq.renameTo(seqNew);
            boolean move4=oriTree.renameTo(oriTreeNew);
            boolean move5=as.renameTo(asNew);
            boolean move6=rba.renameTo(rbaNew);


            if (!move || ! move2 || !move3 || !move4 || !move5 || !move6) {
                System.out.println("Cannot move raxml-ng results from /extended_tree to /AR directory");
                System.exit(1);
            }

        } catch (IOException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }

    
    private void launchPAML() {
        
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }
        FileWriter fw=null;
        try {
            //paml launch command itself
            //launch paml externally to build the posterior probas on the extended tree
            List<String> com=buildPAMLCommand(fw);
            //launch paml externally to build the posterior probas on the extended tree
            Infos.println("Ancestral reconstruct command: "+com);
            executeProcess(com);
            //here should add code to check that expected output files were created correctly
            //TODO

        } catch (IOException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }

    /**
     *
     * @return
     */
    private List<String> buildPhyMLCommand() {
        List<String> com=new ArrayList<>();
        com.add(ARBinary.getAbsolutePath());
        com.add("--ancestral"); //marginal reconstruct
        com.add("--no_memory_check"); //no interactive questions for mem usage
        com.add("-i"); //align
        com.add(alignPath.getAbsolutePath());
        com.add("-u"); //tree
        com.add(treePath.getAbsolutePath());
        if (ARParameters==null) {
            com.add("-m"); //model
            com.add(model.modelString);
            if (model.isProteinModel()) {
                com.add("-d"); //analysis type
                com.add("aa");
            }
            com.add("-c"); //number of relative substitution rate categories
            com.add(String.valueOf(model.categories));
            com.add("-b"); //neither approximate likelihood ratio test nor bootstrap values are computed
            com.add("0");
            com.add("-v"); //proportion of invariable sites
            com.add("0.0");
            com.add("-o"); //rate parameters are optimised
            com.add("r");
            com.add("-a"); //gamma shape param
            com.add(String.valueOf(model.alpha));
            com.add("-f"); //base frequencies based on aligned
            com.add("e"); 
            //com.add("--quiet"); //no interactive questions
            //phymlversion 3.3.20180621 do not accept anymore duplicate seqs
            if (phymlAcceptsDuplicates) {
                com.add("--leave_duplicates");
            }
        } else {
            //if parameters given by user via --arparameters, forget previous 
            //command and use this one instead.
            com.addAll(Arrays.asList(ARParameters.split(" ")));
        }
        
        return com;
    }

    /**
     *
     * @return
     */
    private List<String> buildRAXMLNGCommand() {
        List<String> com=new ArrayList<>();
        com.add(ARBinary.getAbsolutePath());
        com.add("--ancestral"); //marginal reconstruct
        com.add("--msa"); //align
        com.add(alignPath.getAbsolutePath());
        com.add("--tree"); //tree
        com.add(treePath.getAbsolutePath());
        com.add("--threads");
        com.add(Integer.toString(threads));
        com.add("--redo"); //no error when launched previously
        com.add("--precision");
        com.add("9");
        com.add("--seed");
        com.add("1");
        com.add("--force");
        com.add("msa"); //do not check for duplicate sequences and empty columns, controlled on rappas side
        if (ARParameters==null) {
            com.add("--data-type"); //analysis type
            if (model.isProteinModel()) {
                com.add("AA");
            } else {
                com.add("DNA");
            }
            com.add("--model"); //model
            StringBuilder sb=new StringBuilder();
            sb.append(model.modelString);
            sb.append("+G");
            sb.append(model.categories);
            sb.append("{");
            sb.append(model.alpha);
            sb.append("}");
            sb.append("+IU{0}");
            sb.append("+FC");
            com.add(sb.toString());
            com.add("--blopt");
            com.add("nr_safe");
            com.add("--opt-model");
            com.add("on");
            com.add("--opt-branches");
            com.add("on");
        } else {
            //if parameters given by user via --arparameters, forget previous
            //command and use this one instead.
            com.addAll(Arrays.asList(ARParameters.split(" ")));
        }
        return com;
    }

    
    /**
     * execute paml AR, by 1st creating config files, then execution.
     * Important, input tree should NOT have node labels !
     * @param fwCTLFile paml ctl file
     * @return
     * @throws IOException 
     */
    private List<String> buildPAMLCommand(FileWriter fwCTLFile) throws IOException {
        List<String> com=new ArrayList<>();
        com.add(ARBinary.getAbsolutePath());
        //first of all, write model file
        String modelString=model.getPAMLEquivalent();
        
        if (ARParameters==null) {
            //content of paml control file (*.ctl)
            StringBuilder sb=new StringBuilder();
            if(!(s instanceof AAStates)) {
                //convert model to PAML model identifier
                sb.append("seqfile = "+alignPath.getAbsolutePath()+"\n");
                sb.append("treefile = "+treePath.getAbsolutePath()+"\n");
                sb.append("outfile = "+ARPath.getAbsolutePath()+File.separator+"paml_output"+"\n");
                sb.append("noisy = 3   * 0,1,2,3: how much rubbish on the screen\n");
                sb.append("verbose = 2   * set to 2 to output posterior proba distribution\n");
                sb.append("runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI\n");
                sb.append("model = "+modelString+"   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n");
                sb.append("Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n");
                sb.append("* ndata = 100\n");
                sb.append("clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n");
                sb.append("fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below\n");
                sb.append("kappa = 5  * initial or fixed kappa\n");
                sb.append("fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below\n");
                sb.append("alpha = "+String.valueOf(model.alpha)+"   * initial or fixed alpha, 0:infinity (constant rate)\n");
                sb.append("Malpha = 0   * 1: different alpha's for genes, 0: one alpha\n");
                sb.append("ncatG = "+String.valueOf(model.categories)+"   * # of categories in the dG, AdG, or nparK models of rates\n");
                sb.append("nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK\n");
                sb.append("nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n");
                sb.append("getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n");
                sb.append("RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n");
                sb.append("Small_Diff = 7e-6\n");
                sb.append("cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n");
                sb.append("* icode = 0  * (with RateAncestor=1. try \"GC\" in data,model=4,Mgene=4)\n");
                sb.append("fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                sb.append("method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n");
                fwCTLFile = new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"baseml.ctl"));
                Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"baseml.ctl");
                fwCTLFile.append(sb);
                fwCTLFile.close();
                com.add(ARPath.getAbsolutePath()+File.separator+"baseml.ctl");

            } else {   
                //if protein analysis, empiric models defined in files
                //need to copy them
                File modelFile = new File(ARPath.getAbsolutePath()+File.separator+modelString);
                //then write qsub array shell script in the working directory
                Infos.println("Loading PAML model file: "+modelString);
                InputStream resourceAsStream = this.getClass().getClassLoader().getResourceAsStream("models/"+modelString);
                Files.copy(resourceAsStream, modelFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                resourceAsStream.close();
              
                sb.append("      seqfile = "+alignPath.getAbsolutePath()+"\n");
                sb.append("     treefile = "+treePath.getAbsolutePath()+"\n");
                sb.append("      outfile = "+ARPath.getAbsolutePath()+File.separator+"paml_output"+"\n");
                sb.append("        noisy = 3   * 0: concise; 1: detailed, 2: too much\n");
                sb.append("      verbose = 2   * set to 2 to output posterior proba distribution\n");
                sb.append("      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI  -2: pairwise\n");
                sb.append("      seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs\n");
                sb.append("    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n");
                sb.append("      * ndata = 100\n");
                sb.append("        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n");
                sb.append("       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n");
                sb.append("   aaRatefile = "+modelFile.getAbsolutePath()+"  * only used for aa seqs with model=empirical(_F)\n");
                sb.append("                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\n");
                sb.append("        model = 2\n");
                sb.append("                   * models for codons:\n");
                sb.append("                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n");
                sb.append("                   * models for AAs or codon-translated AAs:\n");
                sb.append("                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F\n");
                sb.append("                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n");
                sb.append("      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n");
                sb.append("                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n");
                sb.append("                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n");
                sb.append("                   * 13:3normal>0\n");
                sb.append("        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n");
                sb.append("        Mgene = 0\n");
                sb.append("                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n");
                sb.append("                   * AA: 0:rates, 1:separate\n");
                sb.append("    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n");
                sb.append("        kappa = 2  * initial or fixed kappa\n");
                sb.append("    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n");
                sb.append("        omega = .4 * initial or fixed omega, for codons or codon-based AAs\n");
                sb.append("    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n");
                sb.append("        alpha = "+String.valueOf(model.alpha)+" * initial or fixed alpha, 0:infinity (constant rate)\n");
                sb.append("       Malpha = 0  * different alphas for genes\n");
                sb.append("        ncatG = "+String.valueOf(model.categories)+"  * # of categories in dG of NSsites models\n");
                sb.append("        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates\n");
                sb.append(" RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n");
                sb.append("   Small_Diff = .5e-6\n");
                sb.append("    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n");
                sb.append("  fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                sb.append("       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n");
                sb.append("        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n");
                sb.append(" RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n");
                sb.append("   Small_Diff = 7e-6\n");
                sb.append("    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n");
                sb.append("      * icode = 0  * (with RateAncestor=1. try \"GC\" in data,model=4,Mgene=4)\n");
                sb.append("  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                sb.append("       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n");
                //write paml ctl file 
                fwCTLFile = new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"codeml.ctl"));
                Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"codeml.ctl");
                fwCTLFile.append(sb);
                fwCTLFile.close();
                com.add(ARPath.getAbsolutePath()+File.separator+"codeml.ctl");

            }  
        } else {
            //if parameters given by user via --arparameters, forget previous 
            //command and use this one instead.
            com.addAll(Arrays.asList(ARParameters.split(" ")));
        }
        //Note: PAML is doing something weird with the stdout buffer at the very beginning of its C code. 
        //to avoid this and restore default behaviour of printf, we use the --stdout-no-buf hidden option
        com.add("--stdout-no-buf"); 
        
        return com;
    }
    
    
    /**
     * execution itself
     * @param com 
     */
    private void executeProcess(List<String> com) throws IOException {
        
        ProcessBuilder pb = new ProcessBuilder(com);
        //pb.environment().entrySet().stream().forEach((e) ->{ System.out.println(e.getKey()+"="+e.getValue()); });
        //env.put("VAR1", "myValue"); env.remove("OTHERVAR");
        pb.directory(ARPath);
        Infos.println("Current directory:"+pb.directory().getAbsolutePath());
        pb.redirectErrorStream(false);
        pb.redirectOutput(ProcessBuilder.Redirect.PIPE);
        pb.redirectInput(ProcessBuilder.Redirect.PIPE);
        Infos.println("External process operating reconstruction is logged in: "+new File(ARPath.getAbsolutePath()+File.separator+"AR_sdtout.txt").getAbsolutePath());
        Infos.println("Launching ancestral reconstruction (go and take a coffee, it might take hours if > 5000 leaves!) ...");
        System.out.println("Output from external software:");
        
        Process p = pb.start();
        assert pb.redirectInput() == ProcessBuilder.Redirect.PIPE;
        assert p.getInputStream().read() == -1;
        //redirect sdtout/stdin to files
        FileOutputStream STDOUTOutputStream=new FileOutputStream(new File(ARPath.getAbsolutePath()+File.separator+"AR_sdtout.txt"));
        FileOutputStream STDERROutputStream=new FileOutputStream(new File(ARPath.getAbsolutePath()+File.separator+"AR_sdterr.txt"));
        if (verboseAR)
            inputStreamToOutputStream(p.getInputStream(), System.out);
        inputStreamToOutputStream(p.getInputStream(), STDOUTOutputStream);
        inputStreamToOutputStream(p.getErrorStream(), STDERROutputStream);

        try {
            p.waitFor();
            Thread.sleep(100);
        } catch (InterruptedException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.flush();
        STDOUTOutputStream.flush();
        STDERROutputStream.flush();
        STDOUTOutputStream.close();
        STDERROutputStream.close();
        System.out.println(""); //this line ensures line return after the external process output
        System.out.println("Ancestral reconstruction finished. Return to RAPPAS process.");
    }
    
    
    /**
     * utility function to redirect a stream
     * @param inputStream
     * @param out 
     */
    public static void inputStreamToOutputStream(final InputStream inputStream, final OutputStream out) {
        Thread t = new Thread(new Runnable() {
            @Override
            public void run() {
                try {
                    int d;
                    BufferedInputStream bis=new BufferedInputStream(inputStream);
                    while ((d = bis.read()) != -1) {
                        out.write(d);
                    }
                } catch (IOException ex) {
                    Infos.println(ex.getCause());
                }
            }
        });
        t.setDaemon(true);
        t.start();
    }
    
    /**
     * If phyml, needs to test which version as --ancestral option is recent
     * ; no other way than launching it and getting it from the help header
     */
    private void testPHYMLVersion() throws UnsupportedEncodingException, IOException {
        ProcessBuilder pb = null;
        System.out.print("Testing phyml version :");
        ArrayList<String> phymlHelpCom=new ArrayList<>();
        phymlHelpCom.add(ARBinary.getAbsolutePath());
        phymlHelpCom.add("-h");
        pb = new ProcessBuilder(phymlHelpCom);
        pb.directory(ARPath);
        ByteArrayOutputStream output=new ByteArrayOutputStream();
        Process p = pb.start();
        inputStreamToOutputStream(p.getInputStream(), output);

        try {
            p.waitFor();
            Thread.sleep(100);
        } catch (InterruptedException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        String res=output.toString("UTF-8");
        Pattern pattern = Pattern.compile("- PhyML [0-9a-z+:\\.-]+ -");
        Matcher matcher = pattern.matcher(res);
        while(matcher.find()) {
            this.phymlVersion=matcher.group(0);
            System.out.println(" found "+this.phymlVersion);
        }          
        if (    (!phymlVersion.equals("- PhyML 3.3.20180214 -")) 
             && (!phymlVersion.equals("- PhyML 3.3.20180621 -"))
             && (!phymlVersion.equals("- PhyML 3.3.20190321 -"))
             && (!phymlVersion.equals("- PhyML 3.3.20190909 -"))   
            ){
            System.out.println("Please, use one of these PhyML releases: 3.3.20180214, 3.3.20180621, 3.3.20190321, 3.3.20190909 .");
            System.out.println("They can be downloaded from:  https://github.com/stephaneguindon/phyml/releases");
            System.exit(1);
        }
        //testing --leave_duplicates behaviour because  
        //PhyML 3.3.20180621 shows version 3.3.20180214 in its help !!!
        System.out.print("Testing phyml 'duplicates' option :");
        phymlHelpCom=new ArrayList<>();
        phymlHelpCom.add(ARBinary.getAbsolutePath());
        phymlHelpCom.add("--leave_duplicates");
        pb = new ProcessBuilder(phymlHelpCom);
        pb.directory(ARPath);
        output=new ByteArrayOutputStream();
        p = pb.start();
        //inputStreamToOutputStream(p.getInputStream(), output);
        inputStreamToOutputStream(p.getErrorStream(), output);
        try {
            p.waitFor();
            Thread.sleep(100);
        } catch (InterruptedException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        res=output.toString("UTF-8");
        if (res.contains("unrecognized option '--leave_duplicates'")) {
            System.out.println(" UNSUPPORTED");
            this.phymlAcceptsDuplicates=false;
        } else {
            System.out.println(" SUPPORTED");
        }
        
    }
    
}
