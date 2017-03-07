package etc;

import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLDecoder;
import java.text.DecimalFormat;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * contains several tools for system environement variables
 * @author linard
 */
public class Environement {

    /**
     * print all system environement variables to stdout
     */
    public static void printEnvironementVariables() {
        java.util.Enumeration liste = System.getProperties().propertyNames();
        String cle;
        while( liste.hasMoreElements() ) {
                cle = (String)liste.nextElement();
                System.out.println( cle + " = " + System.getProperty(cle) );
        }
    }

    /**
     * get path of the jar or the class that call this method (name of file included)
     * @param c
     * @return
     * @throws UnsupportedEncodingException 
     */
    public static String getExecutablePath(Class c) throws UnsupportedEncodingException {
        String path = "/" + c.getName().replace('.', '/') + ".class";
        URL url = c.getClass().getResource(path);
        path = URLDecoder.decode(url.toString(), "UTF-8");
        // suppression de  la classe ou du jar du path de l'url
        int index = path.lastIndexOf("/");
        path = path.substring(0, index);
        if (path.startsWith("jar:file:")) {
            // suppression de jar:file: de l'url d'un jar
            // ainsi que du path de la classe dans le jar
            index = path.indexOf("!");
            path = path.substring(9, index);
            return path;
        } else {
            // suppresion du file: de l'url si c'est une classe en dehors d'un jar
            // et suppression du path du package si il est présent.
            path = path.substring(5, path.length());
            Package pack = c.getClass().getPackage();
            if (null != pack) {
                String packPath = pack.toString().replace('.', '/');
                if (path.endsWith(packPath)) {
                    path = path.substring(0, (path.length() - packPath.length()));
                }
            }
            return path;
        }
    }
    /**
     * get path of the jar or the class that call this method (only directory)
     * @param c
     * @return
     * @throws UnsupportedEncodingException
     */
    public static File getExecutablePathWithoutFilename(Class c) {
        try {
            String path=Environement.getExecutablePath(c);
            return new File(path).getParentFile();
        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(Environement.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    /**
     * print memory usage in MB format( max memory, free memory and used memory)
     */
    public static void printMemoryUsageDescription() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(2);
        df.setMinimumFractionDigits(2);
        double allocated=new Double(Runtime.getRuntime().totalMemory()/1048576);
        double allowed=new Double(Runtime.getRuntime().maxMemory()/1048576);
        double free=new Double(Runtime.getRuntime().freeMemory()/1048576);
        double actuallyUsed=allocated-free;
        System.out.println("MEMORY: allocated: " + df.format(allocated) + "MB of " + df.format(allowed)+ "MB (" + df.format(free) +"MB free) used:"+df.format(actuallyUsed)+" MB");
    }

    /**
     * print memory usage in MB format( max memory, free memory and used memory)
     */
    public static String getMemoryUsage() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(2);
        df.setMinimumFractionDigits(2);
        double allocated=new Double(Runtime.getRuntime().totalMemory()/1048576);
        double allowed=new Double(Runtime.getRuntime().maxMemory()/1048576);
        double free=new Double(Runtime.getRuntime().freeMemory()/1048576);
        double actuallyUsed=allocated-free;
        return "MEMORY: allocated: " + df.format(allocated) + "MB of " + df.format(allowed)+ "MB (" + df.format(free) +"MB free) used:"+df.format(actuallyUsed)+" MB";
    }
    
    /**
     * get actual memory usage, which corresponds to (memory_allocated-memory_free)
     * @return 
     */
    public static Double getMemoryUsageAsMB() {
        return new Double((Runtime.getRuntime().totalMemory()/1048576)-(Runtime.getRuntime().freeMemory()/1048576));
    }

    /**
     * return file size as MB, or null if file not exists
     * @param f
     * @return 
     */
    public static Double getFileSize(File f) {
        if(f.exists()){

            double bytes = f.length();
            double kilobytes = (bytes / 1024);
            double megabytes = (kilobytes / 1024);
            return megabytes;
//            double gigabytes = (megabytes / 1024);
//            double terabytes = (gigabytes / 1024);
//            double petabytes = (terabytes / 1024);
//            double exabytes = (petabytes / 1024);
//            double zettabytes = (exabytes / 1024);
//            double yottabytes = (zettabytes / 1024);
        }else{
            Infos.println("File does not exists!");
            return null;
        }    
    }
    
    
}
