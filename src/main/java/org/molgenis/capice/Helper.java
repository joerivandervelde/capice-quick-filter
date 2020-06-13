package org.molgenis.capice;

import org.molgenis.vcf.VcfInfo;
import org.molgenis.vcf.VcfRecord;

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Helper class with some self-contained functions used in CapiceQuickFilter.
 */
public class Helper {

    /**
     * Retrieve genes symbols are annoted by VEP.
     */
    static Set<String> getGenes(Iterator<VcfInfo> infoIter)
    {
        Set<String> genes = new HashSet<String>();
        while(infoIter.hasNext())
        {
            VcfInfo info = infoIter.next();
            String key = info.getKey();
            if(key.equals("CSQ"))
            {
                String val = info.getValRaw();
                String[] valSplit = val.split(",", -1);
                for (String forEachAltAndTranscript : valSplit) {
                    String[] csqField =
                            forEachAltAndTranscript.split("\\|", -1);
                    String gene = csqField[3];
                    if (!gene.isEmpty()) {
                        genes.add(gene);
                    }
                }
            }
        }
        return genes;
    }

    /**
     * Get highest CAPICE score, or NULL if CAPICE is not present.
     * Note that we are not matching exact alt allele here. If one variant
     * has a high enough score, it passes for further interpretation.
     */
    static Double getHighestCapice(Iterator<VcfInfo> infoIter)
    {
        Double highestCapice = null;
        while(infoIter.hasNext())
        {
            VcfInfo info = infoIter.next();
            String key = info.getKey();
            if(key.equals("CAPICE"))
            {
                String val = info.getValRaw();
                String[] valSplit = val.split(",", -1);
                for(String CS : valSplit){
                    double CSdouble = Double.parseDouble(CS);
                    if(highestCapice == null || CSdouble > highestCapice)
                    {
                        highestCapice = CSdouble;
                    }
                }
            }
        }
        return highestCapice;
    }

    /**
     * Get lowest GnomAD allele frequency, or NULL if GnomAD is not present.
     * Note that we are not matching exact alt allele here. If one variant is
     * rare enough, it passes for further interpretation.
     */
    static Double getLowestGnomAD(Iterator<VcfInfo> infoIter)
    {
        Double lowestGnomAD = null;
        while(infoIter.hasNext())
        {
            VcfInfo info = infoIter.next();
            String key = info.getKey();
            if(key.equals("CSQ"))
            {
                String val = info.getValRaw();
                String[] valSplit = val.split(",", -1);
                for (String forEachAltAndTranscript : valSplit) {
                    String[] csqField =
                            forEachAltAndTranscript.split("\\|", -1);
                    String gnomadAFStr = csqField[26];
                    if (!gnomadAFStr.isEmpty()) {
                        double gnomadAF = Double.parseDouble(gnomadAFStr);
                        if (lowestGnomAD == null || gnomadAF < lowestGnomAD) {
                            lowestGnomAD = gnomadAF;
                        }
                    }
                }
            }
        }
        return lowestGnomAD;
    }

    /**
     * Return a VCF line but only genotypes for selected indices.
     */
    static String retainIndices(String variant,
                                          List<Integer> indices)
    {
        StringBuilder out = new StringBuilder();
        String[] split = variant.split("\t", -1);
        for(int i = 0; i < split.length; i++)
        {
            if(i < 9 || indices.contains(i-9))
            {
                out.append(split[i]);
                out.append("\t");
            }

        }
        out.deleteCharAt(out.length()-1);
        return out.toString();
    }

    /**
     * Check if a chromosome is autosomal or not.
     */
    static boolean isAutosomal(VcfRecord vr)
    {
        return vr.getChromosome().matches("\\d+(\\.\\d+)?");
    }
}
