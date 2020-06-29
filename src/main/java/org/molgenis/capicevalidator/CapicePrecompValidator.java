package org.molgenis.capicevalidator;

import net.sf.samtools.util.BlockCompressedInputStream;
import org.molgenis.capice.Helper;
import org.molgenis.genotype.Allele;
import org.molgenis.vcf.VcfReader;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMeta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.sql.SQLOutput;
import java.util.*;

/**
 * CapicePrecompValidator tool.
 * Validations for precomputed CAPICE SNV file.
 */
public class CapicePrecompValidator {

    /*
     * Class variables
     */
    private File input;

    /*
     * Static variables
     */
    private static final String version = "v0.0.1";

    /*
     * Constructor
     */
    public CapicePrecompValidator(File input) {
        this.input = input;
    }

    /**
     * Run the CapicePrecompValidator after constructing.
     */
    void run() throws Exception
    {

        /*
         * Initialize the VCF reader
         */
        BlockCompressedInputStream is = new BlockCompressedInputStream(input);

        /*
         * Keep track of which line we are at
         */
        int lineNr = 0;

        /*
         * Remember current and previous position, must find 3 lines for 1 pos,
         * then next one must be +1 unless next chrom
         */
        int previousPos = -1;

        /*
         * Must find 3 lines exactly for each pos
         */
        int nrOfLinesForPos = 0;

        /*
         * Keep track of all chromosomes that have been seen and that we do
         * not expect to see again (e.g. chrom 1 and 2 when we are at 3)
         */
        Set<String> staleChroms = new HashSet<>();

        /*
         * Keep track of the active chromosome.
         * When this changes, make it stale.
         * Also store min and max pos per chrom.
         */
        String previousChrom = null;
        HashMap<String, Integer> chromMinPos = new HashMap<>();
        HashMap<String, Integer> chromMaxPos = new HashMap<>();

        /*
         * Keep track of unique ref and alts found at specific position
         */
        Set<String> refs = new HashSet<>();
        Set<String> alts = new HashSet<>();

        /*
         * Did we encounter first data line yet, starting with chr 1 ?
         */
        boolean startOfDataFound = false;

        while(true) {
            String line = is.readLine();

            /*
             * Stop when there is no more content
             */
            if (line == null) {
                break;
            }

            /*
             * Loop over any lines that don't start with '1' until we find it
             */
            if (!startOfDataFound && line.startsWith("1")) {
                startOfDataFound = true;
            }
            if (!startOfDataFound) {
                continue;
            }

            /*
             * Check if split length is what we expect
             */
            String[] split = line.split("\t", -1);
            if (split.length != 5) {
                throw new Exception(
                        "Expected length == 5 but found " + split.length +
                                " for line: " + line);
            }

            String currentChrom = split[0];
            int currentPos = Integer.parseInt(split[1]);
            String ref = split[2];
            String alt = split[3];
            Double score = Double.parseDouble(split[4]);

            /*
             * Some sanity checks:
             * CAPICE score is in expected range
             * Position is positive
             * Ref and alt are 1 char
             */
            if(score < 0.0 || score > 1.0)
            {
                throw new Exception("CAPICE score outside 0-1 range: " + score + " at line: " + line);
            }
            if(currentPos < 0){
                throw new Exception("Position negative: " + currentPos + " at line: " + line);
            }
            if(ref.length() != 1)
            {
                throw new Exception("Ref not 1 char at line: " + line);
            }
            if(alt.length() != 1){
                throw new Exception("Alt not 1 char at line: " + line);
            }
            if(!ref.equals("A") && !ref.equals("T") && !ref.equals("G") && !ref.equals("C"))
            {
                throw new Exception("Ref does not equal A, T, G or C at line:" +
                        " " + line);
            }
            if(!alt.equals("A") && !alt.equals("T") && !alt.equals("G") && !alt.equals("C"))
            {
                throw new Exception("Alt does not equal A, T, G or C at line:" +
                        " " + line);
            }

            /*
             * Detect chromosome change. Check if new one was stale.
             * Reset the position for further checks.
             * Add current position as min for current chromosome and
             * previous position as max for previous chromosome
             */
            if(previousChrom != null && !currentChrom.equals(previousChrom))
            {
                if(staleChroms.contains(currentChrom))
                {
                    throw new Exception("Current chrom seen before, is your " +
                            "ordering correct? at line: " + line);
                }
                previousPos = -1;
                staleChroms.add(previousChrom);
                chromMinPos.put(currentChrom, currentPos);
                chromMaxPos.put(previousChrom, previousPos);
            }

            /*
             * Quick check if positions are at least not counting down
             */
            if(previousPos != -1 && currentPos < previousPos){
                throw new Exception("Current pos precedes previous pos: " + currentPos + " < " + previousPos + " at line: " + line);
            }

            /*
             * Detect position change. Do all kinds of checks.
             */
            if(previousPos != -1 && currentPos != previousPos)
            {

                /*
                 * Check if increment is 1
                 */
                if(currentPos - previousPos != 1)
                {
                    throw new Exception("Position increment not 1 around " +
                            "line:" + line);
                }

                /*
                 * Check if previously, there were 3 lines for this position
                 */
                if(nrOfLinesForPos != 3)
                {
                    throw new Exception("Expecting 3 lines per unique " +
                            "position but found " + nrOfLinesForPos + " " +
                            "just before line: " + line);
                }

                /*
                 * Check if ref was unique
                 */
                if(refs.size() != 1)
                {
                    throw new Exception("Non-unique ref around line " + line);
                }

                /*
                 * Check if there were 3 unique alts
                 */
                if(alts.size() != 3)
                {
                    throw new Exception("Expected 3 alts but found " + alts.size() + " around line " + line);
                }

                /*
                 * Reset variables that keep track of things
                 */
                refs = new HashSet<>();
                alts = new HashSet<>();
                nrOfLinesForPos = 0;
            }

            refs.add(ref);
            alts.add(alt);
            previousPos = currentPos;
            previousChrom = currentChrom;
            nrOfLinesForPos++;
            lineNr++;

            if(lineNr % 1000000 == 0)
            {
                System.out.println("Processed " + lineNr + " lines...");
            }
        }

        System.out.println("Done checking " + lineNr + " lines (excl. header)");

        for(String chrom : chromMinPos.keySet())
        {
            System.out.println(chrom + " -> " + chromMinPos.get(chrom) + " " +
                    "to" + chromMaxPos.get(chrom));
        }






    }
}
