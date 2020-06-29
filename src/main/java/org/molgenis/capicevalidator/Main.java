package org.molgenis.capicevalidator;

import org.molgenis.capice.CapiceQuickFilter;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Main class for running CapiceQuickFilter from command-line.
 */
public class Main
{
    public static void main(String args[]) throws Exception
    {
        /*
         * Print help if incorrect number of arguments are given
         */
        if(args.length != 1)
        {
            System.out.println("Please supply 1 arguments:");
            System.out.println("- File location of your CAPICE precomputed scores file.");
            System.exit(0);
        }

        /*
         * Input .VCF.GZ file
         */
        File input = new File(args[0]);
        if(!input.getName().endsWith(".gz"))
        {
            System.out.println("Input CAPICE precomputed scores file name '" + input.getName() + "' does not end in '.gz'. Are you sure this is a valid input?");
            System.exit(0);
        }
        if(!input.exists())
        {
            System.out.println("Input CAPICE precomputed scores file not found at " + input.getAbsolutePath()+".");
            System.exit(0);
        }

        /*
         * Argument parsing done, start CapicePrecompValidator
         */
        System.out.println("Arguments OK. Starting...");
        long start = System.nanoTime();
        CapicePrecompValidator cpv = new CapicePrecompValidator(input);
        cpv.run();
        System.out.println("...completed in " + ((System.nanoTime()-start)/1000000)+"ms.");
    }
}
