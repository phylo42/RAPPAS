/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package etc;

import sun.reflect.Reflection;

/**
 * provides customized outputs if debug mode enabled
 * @author linard
 */
public class Infos {

    /**
     * customized output with new line
     * @param c
     * @param s
     */
    public static void println(Object s) {
        if (System.getProperty("debug.verbose")!=null)
            if (System.getProperty("debug.verbose").equals("1")) {
                System.out.println(getCallerClassName()+"-->  "+s+"");
            }
    }

    /**
     * customized output without new line
     * @param c
     * @param s
     */
    public static void print(Object s) {
        if (System.getProperty("debug.verbose")!=null)
            if (System.getProperty("debug.verbose").equals("1")) {
                System.out.print(getCallerClassName()+"-->  "+s+"");
            }
    }
    
    
    public static String getCallerClassName() { 
        StackTraceElement[] stElements = Thread.currentThread().getStackTrace();
        for (int i=1; i<stElements.length; i++) {
            StackTraceElement ste = stElements[i];
            if (!ste.getClassName().equals(Infos.class.getName()) && ste.getClassName().indexOf("java.lang.Thread")!=0) {
                return ste.getClassName();
            }
        }
        return null;
     }
    
    public static String getCallerCallerClassName() { 
        StackTraceElement[] stElements = Thread.currentThread().getStackTrace();
        String callerClassName = null;
        for (int i=1; i<stElements.length; i++) {
            StackTraceElement ste = stElements[i];
            if (!ste.getClassName().equals(Infos.class.getName())&& ste.getClassName().indexOf("java.lang.Thread")!=0) {
                if (callerClassName==null) {
                    callerClassName = ste.getClassName();
                } else if (!callerClassName.equals(ste.getClassName())) {
                    return ste.getClassName();
                }
            }
        }
        return null;
    }
}
