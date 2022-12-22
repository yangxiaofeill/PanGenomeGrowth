/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenomegrowth;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author fei
 */
public class CPCGrowth {
    
    public  static String versionString = "v1.0";
    
    private static HashMap<String, ArrayList<Region>> getExcludeRegion(String regionExcludeFile){
        System.out.println("[---in getExcludeRegion function---]");
        HashMap<String, ArrayList<Region>> regionExclude = new HashMap<>();   
        try {
            BufferedReader brer = new BufferedReader(new FileReader(regionExcludeFile));
            String line = "";
            while((line = brer.readLine()) != null){
                String[] theEles = line.split("\t");
                String chr = theEles[1];
                int start = new Integer(theEles[2]);
                int end = new Integer(theEles[3]);
                Region oneRegion = new Region();
                oneRegion.chr = chr;
                oneRegion.start = start;
                oneRegion.end = end;
                if(regionExclude.containsKey(chr)){
                    regionExclude.get(chr).add(oneRegion);
                }else{
                    ArrayList<Region> oneChrRegions = new ArrayList<>();
                    oneChrRegions.add(oneRegion);
                    regionExclude.put(chr, oneChrRegions);
                }
            }
            brer.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        }
        return regionExclude;
    }
    
    private static boolean isCen(String chr, int start, HashMap<String, ArrayList<Region>> regionExclude){
        boolean isInCen = false;
        ArrayList<Region> excludeRegion = regionExclude.get(chr);
        for(Region oneRegin : excludeRegion){
            if(start >= oneRegin.start && start <= oneRegin.end){
                // in centromere region
                isInCen = true;
                break;
            }
        }
        return  isInCen;
    }
    
    private static int getLen(String ref, String alt){
        int sameStartLen = 0;
        for(int snpi = 0; ; ++snpi){
            if(snpi >= ref.length() || snpi >= alt.length()){
                break;
            }
            if(ref.charAt(snpi) == alt.charAt(snpi)){
                sameStartLen++;
            }
        }
        int len = 0;
        if(alt.length() > ref.length()){
            len = alt.length() - sameStartLen;
        }else if(alt.length() < ref.length()){
            len = -1*(ref.length() - sameStartLen);
        }else{
            len = ref.length() - sameStartLen;
        }
        return len;
    }
    
    private static int getNonZeroMin (ArrayList<Integer> theRawLen){
        int minV = Integer.MAX_VALUE;
        for(int oneLen : theRawLen){
            if(oneLen > 0){
                if(oneLen < minV){
                    minV = oneLen;
                }
            }
        }
        return minV;
    }
    
    private static ArrayList<ArrayList<Integer>> processTheLen(ArrayList<Integer> theRawLen){
        ArrayList<ArrayList<Integer>> theProLen = new ArrayList<>();
        ArrayList<Integer> currLen = theRawLen;
        while(true){
            int minLen = getNonZeroMin(currLen);
            if(minLen == Integer.MAX_VALUE){
                // all value is zero
                break;
            }
            ArrayList<Integer> minEqualArray = new ArrayList<>();
            ArrayList<Integer> newOne = new ArrayList<>();
            currLen.forEach(oneLen -> {
                if(oneLen > 0){
                    minEqualArray.add(minLen);
                    newOne.add(oneLen - minLen);
                }else{
                    newOne.add(0);
                    minEqualArray.add(0);
                }
            });
            currLen = newOne;
            theProLen.add(minEqualArray);
        }
        return theProLen;
    }
    
    private static ArrayList<ArrayList<Integer>> processALT(String ref, String alt, String[] GTs){
        String[] alts = null;
        if(alt.contains(",")){
            alts = alt.split(",");
        }else{
            alts = new String[1];
            alts[0] = alt;
        }
        
        ArrayList<Integer> theRawLen = new ArrayList<>();
        for(String oneGT : GTs){
            String[] GTEles = oneGT.split("\\|");
            for(String oneGTEle : GTEles){
                if(oneGTEle.equals(".") || oneGTEle.equals("0")){
                    theRawLen.add(0);
                }else{
                    String newAlt = alts[new Integer(oneGTEle) - 1];
                    int onelen = getLen(ref, newAlt);
                    if(onelen > 3000000){
                        onelen = 0;
                    }
                    if(onelen > 0)
                        theRawLen.add(onelen);
                    else
                        theRawLen.add(0);  // remove deletion
                }
            }
        }
        return processTheLen(theRawLen);
    }
    
    private static void processALTAndGT(String line, int linenum, int sampleNum, 
            int sampleStartPosInVcf, BufferedWriter bw, float corePer, float commonPer){
//        System.out.println("[---processing ALT and GTs ----]");
        String[] theEles = line.split("\t");
        String ref = theEles[3];
        String alt = theEles[4];

        String[] GTs = new String[sampleNum];
        for(int i = sampleStartPosInVcf; i < theEles.length; ++i){
            GTs[i - sampleStartPosInVcf] = theEles[i];
        }
        ArrayList<ArrayList<Integer>> theLengths = processALT(ref, alt, GTs);
        
        
        try{
            for(ArrayList<Integer> oneLengthArray : theLengths){
                boolean isLargerThan1Mb = false;
                boolean isLargerThan200Kb = false;
                String lineClass = "single";
                int nonZeroNum = 0;
                for(int oneLen : oneLengthArray){
                    if(oneLen > 0){
                        ++nonZeroNum;
                        if(oneLen > 1000000){  
                            isLargerThan1Mb = true;
                        }
                        if(oneLen > 200000){ 
                            isLargerThan200Kb = true;
                        }
                    }
                    bw.write(oneLen + "\t");
                }
                bw.write(linenum + "\t");
                if(nonZeroNum >= 2 && nonZeroNum < (sampleNum * 2 * commonPer)){
                    lineClass = "poly";
                }
                if(nonZeroNum >= (sampleNum * 2 * commonPer) && nonZeroNum < (sampleNum * 2 * corePer)){
                    lineClass = "common";
                }
                if(nonZeroNum >= (sampleNum * 2 * corePer)){
                    lineClass = "core";
                }
                if(isLargerThan200Kb){
                    lineClass = "poly";
                }
                if(isLargerThan1Mb){
                    lineClass = "single";
                }
                
                bw.write(lineClass + "\n");
            }
        }catch (FileNotFoundException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public static ArrayList<String> getChrList(String chrList){
         ArrayList<String> theChrs = new ArrayList<>();
        try {
            BufferedReader brChr = new BufferedReader(new FileReader(chrList));
            String line = "";
            while((line = brChr.readLine()) != null){
                String[] theEles = line.split("\t");
                theChrs.addAll(Arrays.asList(theEles));
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        }
        return theChrs;
    }
    
    public static void getHapDetph(String infile, int sampleNum, 
            int idx, String chrListFile, 
            float corePer, float commonPer, 
            int superLT, String regionExcludeFile, 
            boolean isRemoveN){

        // infile is the gziped vcf file
        // sampleNum is the number of samples
        // sampleStartPosInVcf is the start column in (from 0) of the first sample in vcf file
        try {
            
            String line = "";
            BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile))));
            BufferedWriter bw = new BufferedWriter(new FileWriter(infile + ".hapDepth_PGG_" + versionString));
            int linenum = 0;
            
            ArrayList<String> theChrsList = null;
            boolean isFilterChr = false;
            if(chrListFile != null){
                theChrsList = getChrList(chrListFile);
//                for(String oneChr : theChrsList){
//                    System.out.println(oneChr);
//                }
                isFilterChr = true;
            }
            
            while((line = br.readLine()) != null){
                if(line.startsWith("##"))
                    continue;
                if(line.startsWith("#")){
                    // get the sample ID
                    String[] theIDs = line.split("\t");
                    for(int i = idx; i < theIDs.length; ++i){
                        bw.write(theIDs[i] + "_H1" + "\t");
                        bw.write(theIDs[i] + "_H2" + "\t");
                    }
                    bw.write("lineIndex\tlineClass\n");
                    continue;
                }
                ++linenum;
                
                String[] theEles = line.split("\t");
                
                if(isFilterChr){
                    if(!theChrsList.contains(theEles[0])){  // remain the chr1-Y
                        System.out.println("line " + linenum + ", the Chr: " + theEles[0] + " should be filter");
                        continue;
                    }
                }
                

                // remain the N region in reference
                // these code use to remove N region anc cen region
                if(isRemoveN){  // remain N or not
                    if(theEles[3].contains("N")){  
                        continue;  
                    }
                }
                
                String chr = theEles[0];
                int start = new Integer(theEles[1]);
                
                if(regionExcludeFile != null){ // remain cen or not
                    HashMap<String, ArrayList<Region>> regionExclude = getExcludeRegion(regionExcludeFile);
                    boolean isInCen = isCen(chr, start, regionExclude);
                    if(isInCen){
                        System.out.println(linenum + "\t" + chr + "\t" + start + "\t in cen");
                        continue;
                    }
                }
                
                processALTAndGT(line, linenum, sampleNum, idx, bw, corePer, commonPer);
                if(linenum %10000 == 0){
                    bw.flush();
                    System.out.println("process " + linenum + " lines");
                }
            }
            br.close();
            bw.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(CPCGrowth.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    public static void getHelpMessage(){
        System.out.println("The version: " + versionString);
        System.out.println("java -jar PanGenomeGrowth.jar "
                + "-f vcfFile.gz -n sampleNum -idx index -chr chrlist "
                + "-core 0.95 -common 0.05 -superLT 3000000 -r region.bed "
                + "-N T -h");
        
        System.out.println("*********The parameters*******************");
        System.out.println("-f\t\t[required] The VCF File in gzipped format");
        System.out.println("-n\t\t[required] The number of samples");
        System.out.println("-idx\t\t[optional] The column index (start from 0) of the first sample in vcf, default: 9");
        System.out.println("-chr\t\t[optional] The file of list of chrs to be kept in results. "
                           + "All chrs in one line and seperating by tab, default: all chromosomes");
        System.out.println("-core\t\t[optional] The percentage of all haplotyes represents the core, default: 0.95");
        System.out.println("-common\t\t[optional] The percentage of all haplotypes represents the common, default: 0.05");
        System.out.println("-superLT\t\t[optional] The length threshold of super larger segment that need to remove, default: 3000000");
        System.out.println("-r\t\t[optional] The regions stored in bed file format that need to remove, default: no");
        System.out.println("-N\t\t[optional] True(T) or False(F) remove regions with N or not, default: False"); 
        System.out.println("-h\t\t print this help message"); 
        System.out.println("*******************************************");
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        for(int i = 0; i < args.length; ++i){
            if(args[i].equals("-h") || args[i].equals("--help")){
                getHelpMessage();
                System.exit(0);
            }
        }
        getHelpMessage();
        
        String theVCFFile = null;    // required
        int sampleNum = -1;          // required
        int idx = 9;                 // optional
        String chrListFile = null;   // optional 
        float coreP = 0.95f;         // optional 
        float commonP = 0.05f;       // optional
        int superLT = 3000000;       // optional
        String regionExcludeFile = null; // optional
        boolean isRemoveN = false;   // optional
        
        if(args.length == 0){
            getHelpMessage();
        }
        
        int argsi = 0;
        while(argsi < args.length){
            if(args[argsi].equals("-f")){
                theVCFFile = args[argsi + 1];
            }
            if(args[argsi].equals("-n")){
                sampleNum = new Integer(args[argsi+1]);
            }
            if(args[argsi].equals("-idx")){
                idx = new Integer(args[argsi + 1]);
            }
            
            if(args[argsi].equals("-chr")){
                chrListFile = args[argsi + 1];
            }
            
            if(args[argsi].equals("-core")){
                coreP = new Float(args[argsi + 1]);
            }
            
            if(args[argsi].equals("-common")){
                commonP = new Float(args[argsi + 1]);
            }
            
            if(args[argsi].equals("-superLT")){
                commonP = new Float(args[argsi + 1]);
            }
            
            if(args[argsi].equals("-N")){
                if(args[argsi + 1].toUpperCase().equals("T") || args[argsi + 1].toUpperCase().equals("TRUE")){
                    isRemoveN = true;
                }
            }
            if(args[argsi].equals("-r")){
                regionExcludeFile = args[argsi + 1];
            }
            argsi += 2;
        }
        
        if(theVCFFile == null || sampleNum == -1){
            System.err.println("ERROR: -f vcfFile.gz -n sampleNum are two required parameters, "
                    + "please provided both of them");
            System.exit(1);
        }
        
        System.out.println("*******Everything is ok, begin to run**********");
        System.out.println("The output file is " + theVCFFile + ".hapDepth_PGG_" + versionString);
        
        getHapDetph(theVCFFile, sampleNum, idx, chrListFile, 
                coreP, commonP, superLT, regionExcludeFile, isRemoveN);
        System.out.println("*******The END**********");
    }
    
}
