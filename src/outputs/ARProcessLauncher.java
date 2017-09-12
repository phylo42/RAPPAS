/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package outputs;

import etc.Infos;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import main_v2.ArgumentsParser_v2;

/**
 * class used to setup the maginal ancerstral reconstruction (AR),
 using an external ARBinary. Currently PAML and PHYML are supported.
 * @author ben
 */
public class ARProcessLauncher {
    
    public static final int AR_PAML=1;
    public static final int AR_PHYML=2;
    public static final int AR_FASTML=3;
    
    public int currentProg=AR_PHYML;
    public File ARBinary=null;
    public File ARPath=null;
    public File alignPath=null;
    public File treePath=null;
    boolean verboseAR=true;
    int analysisType=ArgumentsParser_v2.TYPE_DNA;
    
    /**
     * prepare the marginal AR launcher
     * @param AR_PROG one of ARProcessLauncher.AR_PAML,
     * ARProcessLauncher.AR_PHYML ...etc...
     * @param ARBinary
     * @param verboseAR
     * @param analysisType One of ArgumentsParser_v2.TYPE_DNA, ArgumentsParser_v2.TYPE_PROTEIN
     */
    public ARProcessLauncher(File ARBinary, boolean verboseAR, int analysisType) {
        this.ARBinary=ARBinary;
        this.verboseAR=verboseAR;
        this.analysisType=analysisType;
        //test if ARBinary is something supported
        if (ARBinary.getName().contains("phyml")) {
            System.out.println("I guess, from binary name, we are using PhyML (phyml).");
            this.currentProg=AR_PHYML;
        } else if(ARBinary.getName().contains("baseml")) {
            System.out.println("I guess, from binary name, we are using PAML (baseml).");
            this.currentProg=AR_PAML;
        } else if (ARBinary.getName().contains("codeml")) {
            System.out.println("I guess, from binary name, we are using PAML (codeml).");
            this.currentProg=AR_PAML;
        } else {
            System.out.println("AR binary could not be associated to know AR configuration...");
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
        
        switch (this.currentProg) {
            case AR_PAML:
                Infos.println("PAML AR was selected.");
                launchPAML();
                break;
            case AR_PHYML:
                Infos.println("PHYML AR was selected.");
                launchPHYML();
                break;
            case AR_FASTML:
                Infos.println("FASTML AR was selected.");
                launchFASTML();
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
        switch (this.currentProg) {
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
                File stats=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName().split("\\.")[0]+"_phyml_stats.txt");
                File tree=new File(ARPath.getAbsolutePath()+File.separator+alignPath.getName().split("\\.")[0]+"_phyml_tree.txt");
                if (!stats.exists() || !stats.canRead()) {
                    System.out.println("alignName_phyml_stats.txt file do not exists or cannot be read in "+stats.getAbsolutePath());
                    System.exit(1);
                }
                if (!tree.exists() || !tree.canRead()) {
                    System.out.println("alignName_phyml_tree.txt file do not exists or cannot be read in "+tree.getAbsolutePath());
                    System.exit(1);
                }                
                break;
            case AR_FASTML:
                System.out.println("prepareAR() for FastML not supported yet...");
                System.exit(1);
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
        String command=null;
        switch (this.currentProg) {
            case AR_PAML:
                Infos.println("PAML AR was selected.");
                command=preparePAML();
                break;
            case AR_PHYML:
                Infos.println("PHYML AR was selected.");
                command=preparePHYML();
                break;
            case AR_FASTML:
                Infos.println("FASTML AR was selected.");
                command=prepareFASTML();
                break;
            default:
                break;
        }
        return command;
    }
    
    /**
     * build PHYML program command-line for the AR, without execution.
     * @return 
     */
    private String preparePHYML() {
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }

        //launch paml externally to build the posterior probas on the extended tree
        List<String> com=new ArrayList<>();
        com.add(ARBinary.getAbsolutePath());
        com.add("--ancestral"); //marginal reconstruct
        com.add("-i"); //align
        //System.out.println(alignPath.getAbsolutePath());
        com.add(alignPath.getAbsolutePath());
        com.add("-u"); //tree
        com.add(treePath.getAbsolutePath());
        com.add("-m"); //model
        com.add("GTR");
        com.add("-c"); //number of relative substitution rate categories
        com.add("4");
        com.add("-b"); //neither approximate likelihood ratio test nor bootstrap values are computed
        com.add("0");
        com.add("-v"); //proportion of invariable sites
        com.add("0.0");
        com.add("-o"); //no parameter is optimised
        com.add("n");
        com.add("-a"); //gamma shape param
        com.add("0.5");
        com.add("--quiet"); //no interactive questions
        com.add("--no_memory_check"); //no interactive questions for mem usage

        Infos.println("Ancestral reconstruct command: "+com);

        StringBuilder sb=new StringBuilder(com.get(0));
        for (int i = 1; i < com.size(); i++) {
            sb.append(" ");
            sb.append(com.get(i));
        }
        return sb.toString();
    }
    
    /**
     * build FASTML program command-line for the AR, without execution.
     * @return 
     */
    private String prepareFASTML() {
        System.out.println("FASTML AR not supported yet.");
        System.exit(1);
        return null;
    }    

    /**
     * build PHYML program command-line for the AR, without execution but built its .ctl file. 
     */
    private String preparePAML() {
        
        FileWriter fw=null;
        StringBuilder sb2=null;
        try {
            if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
                System.out.println("AR path is not a directory or do not have read rights.");
                System.exit(1);
            }
            //content of paml control file (*.ctl)
            StringBuilder sb=new StringBuilder();
            //paml launch command itself
            //launch paml externally to build the posterior probas on the extended tree
            List<String> com=new ArrayList<>();
            com.add(ARBinary.getAbsolutePath());
            
            switch(analysisType) {
                
                case ArgumentsParser_v2.TYPE_DNA:
                    
                    sb.append("seqfile = "+alignPath.getAbsolutePath()+"\n");
                    sb.append("treefile = "+treePath.getAbsolutePath()+"\n");
                    sb.append("outfile = "+ARPath.getAbsolutePath()+File.separator+"paml_output"+"\n");
                    sb.append("noisy = 3   * 0,1,2,3: how much rubbish on the screen\n");
                    sb.append("verbose = 2   * set to 2 to output posterior proba distribution\n");
                    sb.append("runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI\n");
                    sb.append("model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n");
                    sb.append("Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n");
                    sb.append("* ndata = 100\n");
                    sb.append("clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n");
                    sb.append("fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below\n");
                    sb.append("kappa = 5  * initial or fixed kappa\n");
                    sb.append("fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below\n");
                    sb.append("alpha = 0.433838   * initial or fixed alpha, 0:infinity (constant rate)\n");
                    sb.append("Malpha = 0   * 1: different alpha's for genes, 0: one alpha\n");
                    sb.append("ncatG = 4   * # of categories in the dG, AdG, or nparK models of rates\n");
                    sb.append("nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK\n");
                    sb.append("nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n");
                    sb.append("getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n");
                    sb.append("RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n");
                    sb.append("Small_Diff = 7e-6\n");
                    sb.append("cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n");
                    sb.append("* icode = 0  * (with RateAncestor=1. try \"GC\" in data,model=4,Mgene=4)\n");
                    sb.append("fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                    sb.append("method = 1  * Optimization method 0: simultaneous; 1: one branch a time\n");
                    fw = new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"baseml.ctl"));
                    Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"baseml.ctl");
                    fw.append(sb);
                    fw.close();
                    com.add(ARPath.getAbsolutePath()+File.separator+"baseml.ctl");
                    
                case ArgumentsParser_v2.TYPE_PROTEIN:
                    
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
                    sb.append("   aaRatefile = dat/wag.dat  * only used for aa seqs with model=empirical(_F)\n");
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
                    sb.append("        alpha = 0.5 * initial or fixed alpha, 0:infinity (constant rate)\n");
                    sb.append("       Malpha = 0  * different alphas for genes\n");
                    sb.append("        ncatG = 4  * # of categories in dG of NSsites models\n");
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
                    sb.append("  fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                    sb.append("       method = 1  * Optimization method 0: simultaneous; 1: one branch a time\n");
                    
                    fw = new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"codeml.ctl"));
                    Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"codeml.ctl");
                    fw.append(sb);
                    fw.close();
                    com.add(ARPath.getAbsolutePath()+File.separator+"codeml.ctl");
                    
            }
                    
            //Note: PAML is doing something weird with the stdout buffer at the very beginning of its C code. 
            //to avoid this and restore default behaviour of printf, we use the --stdout-no-buf hidden option
            com.add("--stdout-no-buf"); 
            Infos.println("Ancestral reconstruct command: "+com);     
            //execution
            sb2=new StringBuilder(com.get(0));
            for (int i = 1; i < com.size(); i++) {
                sb2.append(" ");
                sb2.append(com.get(i));
            }
            
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
        return null;
    }
    
    

    private void launchPHYML() {
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }
        
        try {
           
            //launch paml externally to build the posterior probas on the extended tree
            List<String> com=new ArrayList<>();
            com.add(ARBinary.getAbsolutePath());
            com.add("--ancestral"); //marginal reconstruct
            com.add("-i"); //align
            //System.out.println(alignPath.getAbsolutePath());
            com.add(alignPath.getAbsolutePath());
            com.add("-u"); //tree
            com.add(treePath.getAbsolutePath());
            com.add("-m"); //model
            com.add("GTR");
            com.add("-c"); //number of relative substitution rate categories
            com.add("4");
            com.add("-b"); //neither approximate likelihood ratio test nor bootstrap values are computed
            com.add("0");
            com.add("-v"); //proportion of invariable sites
            com.add("0.0");
            com.add("-o"); //no parameter is optimised
            com.add("n");
            com.add("-a"); //gamma shape param
            com.add("0.5");
            com.add("--quiet"); //no interactive questions
            
            Infos.println("Ancestral reconstruct command: "+com);
            
            //execution
            executeProcess(com);
            
            //phyml is written all data files near the input aignment file...
            //we move them to the AR directory
            //files are:
            // 1. alignName_phyml_ancestral_seq   (unused)
            // 2. alignName_phyml_stats.txt       (used)
            // 3. alignName_phyml_tree.txt        (used)
            File stats=new File(alignPath.getAbsolutePath()+"_phyml_stats.txt");
            File tree=new File(alignPath.getAbsolutePath()+"_phyml_tree.txt");
            File seq=new File(alignPath.getAbsolutePath()+"_phyml_ancestral_seq");
            //check that they were correctly created
            if (!stats.exists() || !stats.exists()) {
                System.out.println("Phyml outputs are missing, the process may have failed...");
                System.exit(1);
            }
            File statsNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+"_phyml_stats.txt");
            File treeNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+"_phyml_tree.txt");
            File seqNew=new File(alignPath.getParent().replace("/extended_trees", "/AR")+File.separator+alignPath.getName()+"_phyml_ancestral_seq");
            
            boolean move=stats.renameTo(statsNew);
            boolean move2=tree.renameTo(treeNew);
            boolean move3=seq.renameTo(seqNew);
            
            if (!move || ! move2 || !move3) {
                System.out.println("Could not move phyml results to /AR directory");
                System.exit(1);
            }
            

        } catch (IOException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    
    private void launchFASTML() {
        System.out.println("FASTML AR not supported yet.");
        System.exit(1);
    }    

    /**
     * execute PAML program for the AR. First build the .ctl control file
     * then execute PAML in the same directory.
     * @param ARPath
     * @param alignPath
     * @param treePath
     */
    private void launchPAML() {
        
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }
        FileWriter fw=null;
        try {
            //content of paml control file (*.ctl)
            StringBuilder sb=new StringBuilder();
            //paml launch command itself
            //launch paml externally to build the posterior probas on the extended tree
            List<String> com=new ArrayList<>();
            com.add(ARBinary.getAbsolutePath());

            switch(analysisType) {
                
                case ArgumentsParser_v2.TYPE_DNA:
                    
                    sb.append("seqfile = "+alignPath.getAbsolutePath()+"\n");
                    sb.append("treefile = "+treePath.getAbsolutePath()+"\n");
                    sb.append("outfile = "+ARPath.getAbsolutePath()+File.separator+"paml_output"+"\n");
                    sb.append("noisy = 3   * 0,1,2,3: how much rubbish on the screen\n");
                    sb.append("verbose = 2   * set to 2 to output posterior proba distribution\n");
                    sb.append("runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI\n");
                    sb.append("model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n");
                    sb.append("Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n");
                    sb.append("* ndata = 100\n");
                    sb.append("clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n");
                    sb.append("fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below\n");
                    sb.append("kappa = 5  * initial or fixed kappa\n");
                    sb.append("fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below\n");
                    sb.append("alpha = 0.433838   * initial or fixed alpha, 0:infinity (constant rate)\n");
                    sb.append("Malpha = 0   * 1: different alpha's for genes, 0: one alpha\n");
                    sb.append("ncatG = 4   * # of categories in the dG, AdG, or nparK models of rates\n");
                    sb.append("nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK\n");
                    sb.append("nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n");
                    sb.append("getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n");
                    sb.append("RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n");
                    sb.append("Small_Diff = 7e-6\n");
                    sb.append("cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?\n");
                    sb.append("* icode = 0  * (with RateAncestor=1. try \"GC\" in data,model=4,Mgene=4)\n");
                    sb.append("fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                    sb.append("method = 1  * Optimization method 0: simultaneous; 1: one branch a time\n");
                    fw=new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"baseml.ctl"));
                    Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"baseml.ctl");
                    fw.append(sb);
                    fw.close();
                    
                case ArgumentsParser_v2.TYPE_PROTEIN:
                    
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
                    sb.append("   aaRatefile = dat/wag.dat  * only used for aa seqs with model=empirical(_F)\n");
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
                    sb.append("        alpha = 0.5 * initial or fixed alpha, 0:infinity (constant rate)\n");
                    sb.append("       Malpha = 0  * different alphas for genes\n");
                    sb.append("        ncatG = 4  * # of categories in dG of NSsites models\n");
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
                    sb.append("  fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed\n");
                    sb.append("       method = 1  * Optimization method 0: simultaneous; 1: one branch a time\n");
                    
                    fw = new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"codeml.ctl"));
                    Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"codeml.ctl");
                    fw.append(sb);
                    fw.close();
                    
            }
            
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
        Infos.println("External process operating reconstruction is logged in: "+new File(ARPath.getAbsolutePath()+File.separator+"AR_sdtout.txt").getAbsolutePath());
        Infos.println("Launching ancestral reconstruction (go and take a coffee, it might take hours!) ...");
        try {
            p.waitFor();
            Thread.sleep(1000);
        } catch (InterruptedException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.flush();
        STDOUTOutputStream.flush();
        STDERROutputStream.flush();
        STDOUTOutputStream.close();
        STDERROutputStream.close();
        System.out.println(""); //this line ensures line return after the external process output
        System.out.println("Ancestral reconstruction finished.");
    }
    
    
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
    
}
