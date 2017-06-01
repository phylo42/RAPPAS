/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main_v2;

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

/**
 * class used to setup the maginal ancerstral reconstruction (AR),
 * using an external executable. Currently PAML and PHYML are supported.
 * @author ben
 */
public class ARProcessLauncher {
    
    public static final int AR_PAML=1;
    public static final int AR_PHYML=2;
    public static final int AR_FASTML=3;
    
    public int currentProg=AR_PHYML;
    public File executable=null;
    File ARPath=null;
    File alignPath=null;
    File treePath=null;
    boolean verboseAR=true;
    
    /**
     * prepare the marginal AR launcher
     * @param AR_PROG one of ARProcessLauncher.AR_PAML,
     * ARProcessLauncher.AR_PHYML ...etc...
     */
    public ARProcessLauncher(int AR_PROG, File pathToExecutable,boolean verboseAR) {
        this.currentProg=AR_PROG;
        this.executable=pathToExecutable;
        this.verboseAR=verboseAR;
        
        if (executable==null) {
            System.out.println("Path to executable used for ancestral reconstruction is not set correctly.");
            System.exit(1);
        }
        if ( (!pathToExecutable.isFile()) || (!pathToExecutable.canExecute()) ) {
            System.out.println("AR executable is not a file or do not have execution rights.");
            System.exit(1);
        }
        
        
    }
    
    /**
     * launch marginal AR
     * @param ARPath path in which the external AR program will work
     * @param alignPath alignment in Phylip format, used for the marginal AR
     * @param treePath tree in Newick format, used for the marginal AR
     */
    public void launchAR(File ARPath, File alignPath, File treePath) {
        this.ARPath=ARPath;
        this.alignPath=alignPath;
        this.treePath=treePath;
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
    

    private void launchPHYML() {
        if ( (!ARPath.isDirectory()) || (!ARPath.canWrite()) ) {
            System.out.println("AR path is not a directory or do not have read rights.");
            System.exit(1);
        }
        
        try {
           
            //launch paml externally to build the posterior probas on the extended tree
            List<String> com=new ArrayList<>();
            com.add(executable.getAbsolutePath());
            com.add("--ancestral"); //marginal reconstruct
            com.add("-i"); //align
            System.out.println(alignPath.getAbsolutePath());
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
        
        try {
            StringBuilder sb=new StringBuilder();
            
            sb.append("seqfile = "+alignPath.getAbsolutePath()+"\n");
            sb.append("treefile = "+treePath.getAbsolutePath()+"\n");

            sb.append("outfile = "+ARPath+File.separator+"paml_output"+"\n");
            sb.append("noisy = 2   * 0,1,2,3: how much rubbish on the screen\n");
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
            
            FileWriter fw=new FileWriter(new File(ARPath.getAbsolutePath()+File.separator+"baseml.ctl"));
            Infos.println("Ancestral reconstruciton parameters written in: "+ARPath.getAbsolutePath()+File.separator+"baseml.ctl");
            fw.append(sb);
            fw.close();
            
            //launch paml externally to build the posterior probas on the extended tree
            List<String> com=new ArrayList<>();
            com.add(executable.getAbsolutePath());
            com.add(ARPath.getAbsolutePath()+File.separator+"baseml.ctl");
            Infos.println("Ancestral reconstruct command: "+com);
            
            //execution
            executeProcess(com);
            
            
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
        Infos.println("Launching ancestral reconstruction (go and take a coffee, it mights take hours!) ...");
        try {
            p.waitFor();
        } catch (InterruptedException ex) {
            Logger.getLogger(ARProcessLauncher.class.getName()).log(Level.SEVERE, null, ex);
        }

        STDOUTOutputStream.close();
        STDERROutputStream.close();
        Infos.println("Ancestral reconstruction finished.");
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
