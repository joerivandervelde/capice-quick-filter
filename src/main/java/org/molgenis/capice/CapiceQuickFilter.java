package org.molgenis.capice;

import net.sf.samtools.util.BlockCompressedInputStream;
import org.molgenis.genotype.Allele;
import org.molgenis.vcf.VcfReader;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMeta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.*;

/**
 * CapiceQuickFilter tool. With only CAPICE and GnomAD annotations, perform
 * basic but effective filtering and reporting of potentially clinically
 * interesting variants.
 */
public class CapiceQuickFilter {

    /*
     * Class variables
     */
    private File input;
    private File output;
    private double capiceThreshold;
    private double gnomadThreshold;
    private String caseSampleID;
    private List<String> controlSampleIDs;

    /*
     * Static variables
     */
    private static final String version = "v0.0.1";
    private static final String DE_NOVO = "Potential de novo/uncontrolled hetzygote: ";
    private static final String HOM_ALT = "Potential homozygous                    : ";
    private static final String NON_AUT = "Potential non-autosomal                 : ";
    private static final String COMPHET = "Potential compound heterozygote         : ";

    /*
     * Constructor
     */
    public CapiceQuickFilter(File input, File output, double capiceThreshold, double gnomadThreshold, String caseSampleID, List<String> controlSampleIDs) {
        this.input = input;
        this.output = output;
        this.capiceThreshold = capiceThreshold;
        this.gnomadThreshold = gnomadThreshold;
        this.caseSampleID = caseSampleID;
        this.controlSampleIDs = controlSampleIDs;
    }

    /**
     * Run the CapiceQuickFilter after constructing.
     */
    void run() throws Exception
    {
        /*
         * Counters for reporting
         */
        int totalVariantCount = 0;
        int droppedByGnomAD = 0;
        int droppedByCAPICE = 0;
        int droppedByNullOrRefCaseGeno = 0;
        int droppedByHomZygAltControlGeno = 0;
        int droppedByHetZygAltNoHetComp = 0;
        int variantWithoutGnomAD = 0;
        int variantWithoutCAPICE = 0;

        /*
         * Initialize the VCF reader
         */
        BlockCompressedInputStream is = new BlockCompressedInputStream(input);
        VcfReader r = new VcfReader(is);
        VcfMeta vm = r.getVcfMeta();

        /*
         * Prepare objects to store output
         */
        HashMap<String, List<String>> geneToHetZyg = new HashMap<>();
        HashMap<String, List<String>> reportedVariants = new HashMap<>();
        reportedVariants.put(DE_NOVO, new ArrayList<>());
        reportedVariants.put(HOM_ALT, new ArrayList<>());
        reportedVariants.put(NON_AUT, new ArrayList<>());
        reportedVariants.put(COMPHET, new ArrayList<>());

        /*
         * Create output file writer
         */
        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);

        /*
         * Get the sample names from the VCF meta-data
         * TODO: verify that order is guaranteed
         */
        List<String> sampleNames = new ArrayList<>();
        for(String sample: vm.getSampleNames()){
            sampleNames.add(sample);
        }

        /*
         * Sanity checks: are the sample and control IDs present in the VCF?
         * Also, store indices of case and control samples for later use and
         * a list of all indices for convenience
         */
        if(!sampleNames.contains(caseSampleID))
        {
            throw new Exception("index sample id not found: " + caseSampleID);
        }
        int caseSampleIndex = sampleNames.indexOf(caseSampleID);
        List<Integer> controlSampleIndices = new ArrayList<Integer>();
        for(String control : controlSampleIDs)
        {
            if(!sampleNames.contains(control))
            {
                throw new Exception("control sample id not found: " + control);
            }
            controlSampleIndices.add(sampleNames.indexOf(control));
        }
        List<Integer> allIndices = controlSampleIndices;
        allIndices.add(caseSampleIndex);

        /*
         * Start iterating over the input VCF file
         */
        Iterator<VcfRecord> vcfIt = r.iterator();
        vcfIterator:
        while(vcfIt.hasNext())
        {
            VcfRecord vr = vcfIt.next();
            totalVariantCount++;

            /*
             * Retrieve GnomAD and CAPICE values if present
             */
            Double lowestGnomAD =
                    Helper.getLowestGnomAD(vr.getInformation().iterator());
            Double highestCapice = Helper.getHighestCapice(vr.getInformation().iterator());

            /*
             * Keep track of missing GnomAD and CAPICE values
             */
            if(highestCapice == null)
            {
                variantWithoutCAPICE++;
            }
            if(lowestGnomAD == null)
            {
                variantWithoutGnomAD++;
            }

            /*
             * If not missing, we have reasons to drop variants
             */
            if(highestCapice != null && highestCapice < capiceThreshold)
            {
                droppedByCAPICE++;
                continue;
            }
            if(lowestGnomAD != null && lowestGnomAD > gnomadThreshold)
            {
                droppedByGnomAD++;
                continue;
            }

            /*
             * Now we need to investigate genotypes
             * Store case and control genotype in objects
             */
            List<Allele> caseGeno = null;
            HashMap<String, List<Allele>> controlGeno = new HashMap<>();

            /*
             * Iterate over samples and extract the case and control genotypes
             */
            Iterator<VcfSample> vi = vr.getSamples().iterator();
            VcfSample nextSample;
            int sampleIndex = 0;
            while(vi.hasNext())
            {
                nextSample = vi.next();
                if(sampleIndex == caseSampleIndex)
                {
                    caseGeno = nextSample.getAlleles();
                }
                else if(controlSampleIndices.contains(sampleIndex))
                {
                    controlGeno.put(sampleNames.get(sampleIndex),
                            nextSample.getAlleles());
                }
                sampleIndex++;
            }

            /*
             * Drop variant if case genotype consists of only reference
             * alleles and/or missing alleles
             */
            int caseAltCount = 0;
            for(Allele a: caseGeno)
            {
                if(!a.equals(Allele.ZERO) && !a.equals(vr.getReferenceAllele())) {
                    caseAltCount++;
                }
            }
            if(caseAltCount == 0)
            {
                droppedByNullOrRefCaseGeno++;
                continue;
            }

            /*
             * Analyse further. If one control sample is homozygous, drop the
             * variant. If not, track if 1+ control(s) are heterozygous.
             */
            boolean atLeastOneCtrlWithOneAlt = false;
            for(String key : controlGeno.keySet())
            {
                List<Allele> controlAlleles = controlGeno.get(key);
                int controlAltCount = 0;
                for(Allele a: controlAlleles)
                {
                    if(!a.equals(Allele.ZERO) && !a.equals(vr.getReferenceAllele())) {
                        controlAltCount++;
                    }
                }
                if(controlAltCount == 2)
                {
                    droppedByHomZygAltControlGeno++;
                    continue vcfIterator;
                }
                else if(controlAltCount == 1)
                {
                    atLeastOneCtrlWithOneAlt = true;
                }
            }


            /*
             * There are no homozygous controls. So if the case is homozygous
             * alternative, report it and continue.
             */
            if(caseAltCount == 2)
            {
                reportedVariants.get((Helper.isAutosomal(vr) ? HOM_ALT : NON_AUT)).add(Helper.retainIndices(vr.toString(), allIndices));
                continue;
            }

            /*
             * If case has 1 alt, and there are no controls with alt alleles,
             * it is either de novo (controls present), or an 'uncontrolled'
             * heterozygote (controls not present). Report and continue.
             */
            if(caseAltCount == 1 && !atLeastOneCtrlWithOneAlt)
            {
                reportedVariants.get((Helper.isAutosomal(vr) ? DE_NOVO : NON_AUT)).add(Helper.retainIndices(vr.toString(), allIndices));
                continue;
            }

            /*
             * If case has 1 alt but there are also controls with 1 alt, it
             * could still be compound heterozygous. But we can only tell
             * after we have seen all variants from this gene. Save for later.
             * Exception is variants on allosomes, always report these.
             */
            Set<String> genes = Helper.getGenes(vr.getInformation().iterator());
            if(caseAltCount == 1)
            {
                if(!Helper.isAutosomal(vr))
                {
                    reportedVariants.get(NON_AUT).add(Helper.retainIndices(vr.toString(), allIndices));
                    continue;
                }
                else
                {
                    for(String gene : genes)
                    {
                        if(!geneToHetZyg.containsKey(gene))
                        {
                            geneToHetZyg.put(gene, new ArrayList<>());
                        }
                        geneToHetZyg.get(gene).add(vr.toString());
                    }
                      continue;
                }
            }

            /*
             * We should have covered all states when looping over all
             * variants in the input VCF. If not, crash the program.
             */
            throw new Exception("Bad state: all possibilities should be covered by now. Offending variant: " + vr.toString());
        }


        /*
         * Iterate over the heterozygous variants that may become compound
         * heterozygous variants if there are two in one gene. Keep track
         * which are reported to prevent duplicates.
         */
        Set<String> hasBeenReported = new HashSet<>();
        Set<String> hasBeenDropped = new HashSet<>();
        for(String gene : geneToHetZyg.keySet()) {
            if (geneToHetZyg.get(gene).size() > 1) {
                for (String rec : geneToHetZyg.get(gene)) {
                    if (!hasBeenReported.contains(rec)) {
                        reportedVariants.get(COMPHET).add(Helper.retainIndices(rec, allIndices));
                        hasBeenReported.add(rec);
                    }
                }
            }
        }

        /*
         * Also iterate over the leftovers to make sure all numbers add up.
         * For this, we must also considered those already reported for a
         * different gene.
         */
        for(String gene : geneToHetZyg.keySet()) {
            if (geneToHetZyg.get(gene).size() == 1)
            {
                String rec = geneToHetZyg.get(gene).get(0);
                if(!hasBeenDropped.contains(rec) && !hasBeenReported.contains(rec)) {
                    droppedByHetZygAltNoHetComp++;
                    hasBeenDropped.add(rec);
                }
            }
        }

        /*
         * Count total reported and total dropped
         */
        int totalRep =
                reportedVariants.get(HOM_ALT).size() + reportedVariants.get(DE_NOVO).size() + reportedVariants.get(COMPHET).size() + reportedVariants.get(NON_AUT).size();
        int totalDrop =
                droppedByGnomAD + droppedByCAPICE + droppedByNullOrRefCaseGeno + droppedByHomZygAltControlGeno + droppedByHetZygAltNoHetComp;

        /*
         * Print the header with information in theo utput VCF file.
         * TODO: retain original header for proper meta-data
         * TODO: assign VCF version
         * TODO: perhaps sort output by chrom/pos instead of category
         */
        bw.write("## Output of CapiceQuickFilter " + version + "\n");
        bw.write("## Settings:\n");
        bw.write("## - Input file: " + input.getAbsolutePath() + "\n");
        bw.write("## - Output file: " + output.getAbsolutePath() + "\n");
        bw.write("## - CAPICE threshold: " + capiceThreshold + "\n");
        bw.write("## - GnomAD threshold: " + gnomadThreshold + "\n");
        bw.write("## - Case sample ID: " + caseSampleID + "\n");
        bw.write("## - Control sample IDs: " + controlSampleIDs + "\n");
        bw.write("## Total number of variants processed: " + totalVariantCount + "\n");
        bw.write("## Total number of potential candidates found: " + totalRep + "\n");
        bw.write("## Breakdown of potential candidates by type:" + "\n");
        bw.write("## - " + HOM_ALT + reportedVariants.get(HOM_ALT).size() + "\n");
        bw.write("## - " + DE_NOVO + reportedVariants.get(DE_NOVO).size() + "\n");
        bw.write("## - " + COMPHET + reportedVariants.get(COMPHET).size() + "\n");
        bw.write("## - " + NON_AUT + reportedVariants.get(NON_AUT).size() + "\n");
        bw.write("## Total number of variants dropped: " + totalDrop + "\n");
        bw.write("## Breakdown of dropped variants by reason:" + "\n");
        bw.write("## - CAPICE score below threshold = " + droppedByCAPICE + "\n");
        bw.write("## - GnomAD allele frequency over threshold = " + droppedByGnomAD + "\n");
        bw.write("## - Case genotype null or reference = " + droppedByNullOrRefCaseGeno + "\n");
        bw.write("## - Homozygous control was present = " + droppedByHomZygAltControlGeno + "\n");
        bw.write("## - Flagged for compound but no second hit: " + droppedByHetZygAltNoHetComp + "\n");
        bw.write("## Additional information:" + "\n");
        bw.write("## - Variants without GnomAD annotation: " + variantWithoutGnomAD + "\n");
        bw.write("## - Variants without CAPICE annotation: " + variantWithoutCAPICE + "\n");
        bw.write("## Potential candidates categorized by type (full info below, can be copy-pasted side by side):" + "\n");
        for(String key : reportedVariants.keySet())
        {
            for(String variant : reportedVariants.get(key))
            {
                bw.write("## " + key + (variant.length() > 50 ?
                        variant.substring(0, 50).replace("\t", " ") :
                        variant.replace("\t", " ")) + "\n");
            }
        }

        /*
         * Print the VCF columns with sample names and then all variant data.
         * By sorting the indices first, we maintain the same order used in
         * retainIndices() to print the genotypes.
         */
        Collections.sort(allIndices);
        StringBuilder sb = new StringBuilder();
        for(int index : allIndices)
        {
            sb.append("\t");
            sb.append(sampleNames.get(index));
        }
        bw.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" + sb.toString() + "\n");
        for(String key : reportedVariants.keySet())
        {
            for(String variant : reportedVariants.get(key))
            {
                bw.write(variant + "\n");
            }
        }

        /*
         * Flush and close file writer.
         */
        bw.flush();
        bw.close();
    }
}
