using System;

namespace DNA
{
    class MainClass
    {
        // Method for counting occurence of nucleotide in sequence.
        static int CountNucl(char nucl, string seq)
        {
            int count = 0;
            foreach (char x in seq)
            {
                if (x == nucl)
                {
                    count++;
                }
            }
            return count;
        }

        public static void Main(string[] args)
        {
            // Load in DNA sequence.
            string seq = System.IO.File.ReadAllText(
                @"/Users/david/GitHub/Rosalind/DNA/DNA/" +
                "rosalind_dna_1_dataset.txt");

            // Write out nucleotide counts to console.
            Console.WriteLine($"" +
                $"{CountNucl('A', seq)} " +
                $"{CountNucl('C', seq)} " +
                $"{CountNucl('T', seq)} " +
                $"{CountNucl('G', seq)} ");
        }
    }
}
