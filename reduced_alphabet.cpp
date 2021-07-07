#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <set>

using namespace seqan;
using namespace std;

struct ModifyStringOptions 
{
   CharString queryFileName;
   CharString referenceFileName;
   CharString outputFileName;
   int kmer_size;
   bool donot_reduce_alphabet;
   CharString type;
};

seqan::ArgumentParser::ParseResult parseCommandLine(
	ModifyStringOptions & options, 
	int argc, 
	char const ** argv) {

   ArgumentParser parser("reduced_alphabet");
   addOption(parser, ArgParseOption("q", "query-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE,
                                    "IN"));
   setRequired(parser, "query-file");

   addOption(parser, ArgParseOption("r", "reference-file", 
                                    "Path to the input filter file", 
                                    ArgParseArgument::INPUT_FILE, 
                                    "IN"));
   setRequired(parser, "reference-file");

   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE,
                                    "OUT"));
   setRequired(parser, "output-file");

   addOption(parser, ArgParseOption("k", "kmer-size", "k-mer size",
             ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "kmer-size", "3");

   addOption(parser, ArgParseOption("nr", "no-reduce",
             "Do not reduce the alphabet. Calculate all kmers, useful \
              for comparing the reduced alphabet with a standard\
               calculation."));

   addOption(parser, ArgParseOption("t", "distance-type",
                                    "The method of calculating the distance \
                                    between two sequences. For descriptions of \
                                    distance please refer to the wiki. ",
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "distance-type",
                  "euclid manhattan bc ngd canberra normalised_canberra cosine");
   setDefaultValue(parser, "distance-type", "cosine");

   setShortDescription(parser, "reduced_alphabet");
   setVersion(parser, "0.0.1");
   setDate(parser, "July 2021");
   addUsageLine(parser, "-q query.fa -r reference.fa -o output.txt \
                         [\\fIOPTIONS\\fP] ");
   addDescription(parser, "Alignment free but using a reduced alphabet.");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   if (res != ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.queryFileName, parser, "query-file");
   getOptionValue(options.referenceFileName, parser, "reference-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.kmer_size, parser, "kmer-size");
   options.donot_reduce_alphabet = isSet(parser, "no-reduce");
   getOptionValue(options.type, parser, "distance-type");

   return ArgumentParser::PARSE_OK;
}

std::set<CharString> getAllKmers(map<CharString, int> &qry_counts, 
                                 map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers;
   for(auto i : qry_counts)
      kmers.insert(i.first);
   for(auto i : ref_counts)
      kmers.insert(i.first);

   return kmers;
}

double cosine(map<CharString, int> &qry_counts, map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers = getAllKmers(qry_counts, ref_counts);

   long long int sumqCrC = 0;
   long long int sumqC2 = 0;
   long long int sumrC2 = 0;

   for(auto i : kmers)
   {
      sumqCrC += (long long int)qry_counts[i] * (long long int)ref_counts[i];
      sumqC2 += (long long int)qry_counts[i] * (long long int)qry_counts[i];
      sumrC2 += (long long int)ref_counts[i] * (long long int)ref_counts[i];
   }

   double score = sumqCrC / (sqrt(sumqC2) * sqrt(sumrC2));
   return (1-score);
}

double euler(map<CharString, int> &qry_counts, map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers = getAllKmers(qry_counts, ref_counts);

   double score = 0.0;
   long long int rN = 0;
   long long int qN = 0;

   for(auto i : kmers)
   {
      rN = rN + qry_counts[i];
      qN = qN + ref_counts[i];
   }

   for(auto i : kmers)
   {
      double rF = qry_counts[i] / (double)rN;
      double qF = ref_counts[i] / (double)qN;
      score = score + (pow((rF - qF), 2));
   }

   return pow(score, 0.5);
}

double bray_curtis_distance(map<CharString, int> &qry_counts, map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers = getAllKmers(qry_counts, ref_counts);

   double sumMinus = 0.0;
   double sumPlus = 0.0;

   for(auto i : kmers)
   {
      sumMinus = sumMinus + abs((long long int)ref_counts[i] - (long long int)qry_counts[i]);
      sumPlus = sumPlus + abs((long long int)ref_counts[i] + (long long int)qry_counts[i]);
   }

   double result = (double)sumMinus/(double)sumPlus;
   return result;
}

double normalised_google_distance(map<CharString, int> &qry_counts, map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers = getAllKmers(qry_counts, ref_counts);

   double score = 0.0;
   double sumqC = 0.0;
   double sumrC = 0.0;
   double sum_min_qr = 0.0;

   for(auto i : kmers)
   {
      sumrC += qry_counts[i];
      sumqC += ref_counts[i];

      if(qry_counts[i] < ref_counts[i])
         sum_min_qr += qry_counts[i];
      else
         sum_min_qr += ref_counts[i];
   }

   double sum_max, sum_min;

   if(sumqC > sumrC)
   {
      sum_max = sumqC;
      sum_min = sumrC;
   }
   else
   {
      sum_max = sumrC;
      sum_min = sumqC;
   }

   double sum_all = sumqC + sumrC;

   return (sum_max - sum_min_qr) / (sum_all - sum_min);
}

double canberra(map<CharString, int> &qry_counts, map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers = getAllKmers(qry_counts, ref_counts);

   double score = 0.0;
   for(auto i : kmers)
   {
      double p1 = abs((long long int)qry_counts[i] - (long long int)ref_counts[i]);
      double p2 = abs((long long int)qry_counts[i]) + abs((long long int)ref_counts[i]);
      if(p1 == 0 && p2 == 0)
         score += 0.0;
      else if(p2 == 0)
         score += 0.0;
      else
         score += (double)p1/(double)p2;
   }
   return score;
}

double normalised_canberra(map<CharString, int> &qry_counts, map<CharString, int> &ref_counts)
{
   std::set<CharString> kmers = getAllKmers(qry_counts, ref_counts);

   double score = 0.0;
   for(auto i : kmers)
   {
      double p1 = abs((long long int)qry_counts[i] - (long long int)ref_counts[i]);
      double p2 = abs((long long int)qry_counts[i]) + abs((long long int)ref_counts[i]);
      if(p1 == 0 && p2 == 0)
         score += 0.0;
      else if(p2 == 0)
         score += 0.0;
      else
         score += (double)p1/(double)p2;
   }
   return score / length(kmers);
}

CharString sortSeq(CharString seq)
{
   string copySeq = toCString(seq);
   std::sort(copySeq.begin(), copySeq.end());
   CharString converted = copySeq;
   return copySeq;
}

int count(CharString seq, int kmer_size, map<CharString, int> &counts,
          ModifyStringOptions options)
{
   for(int i = 0; i < length(seq)-kmer_size; i++)
   {
      CharString inf = infix(seq, i, i+kmer_size);
      if(options.donot_reduce_alphabet == true)
         counts[inf]++;
      else
         counts[sortSeq(inf)]++;
   }
   return 0;
}

int precalculateRefCounts(vector<map<CharString, int>> &referenceCounts, 
                          StringSet<CharString> &ids,
                          StringSet<CharString> &seqs,
                          int kmer_size,
                          ModifyStringOptions options)
{
   for(auto seq : seqs)
   {
      map<CharString, int> counts;
      count(seq, kmer_size, counts, options);
      referenceCounts.push_back(counts);
   }
   return 0;
}

// A basic template to get up and running quickly
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options,argc, argv);
   if (res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // this will begin super simple
/*
    > Count only the kmers that exist
       > sort the kmer, this will reduce the kmer size
    > do a distance calculation


*/

   // read the query file into memory
   SeqFileIn queryFileIn;
   if(!open(queryFileIn, (toCString(options.queryFileName))))
   {
      cerr << "Error: could not open query file ";
      cerr << toCString(options.queryFileName) << endl;
      return 1;
   }

   // read the reference into memory
   SeqFileIn referenceFileIn;
   if(!open(referenceFileIn, (toCString(options.referenceFileName))))
   {
      cerr << "Error: could not open reference file ";
      cerr << toCString(options.referenceFileName) << endl;
      return 1;
   }   

   StringSet<CharString> ids;
   StringSet<CharString> seqs;
   readRecords(ids, seqs, referenceFileIn);

   vector<map<CharString, int>> referenceCounts;
   precalculateRefCounts(referenceCounts, ids, seqs, options.kmer_size, options);

   while(!atEnd(queryFileIn))
   {
      CharString id;
      CharString seq;
      readRecord(id, seq, queryFileIn);

      map<CharString, int> counts;
      count(seq, options.kmer_size, counts, options);
      int j = 0;
      for(auto i : referenceCounts)
      {
         cout << id << "\t" << ids[j] << "\t";

         if(options.type == "cosine")
            cout << cosine(counts, i);
         else if(options.type == "canberra")
            cout << canberra(counts, i);
         else if(options.type == "normalised_canberra")
            cout << normalised_canberra(counts, i);
         else if(options.type == "ngd")
            cout << normalised_google_distance(counts, i);
         else if(options.type == "bc")
            cout << bray_curtis_distance(counts, i);
         else if(options.type == "euclid")
            cout << euler(counts, i);

         cout << endl;
         j++;
      }
   }

   close(referenceFileIn);
   close(queryFileIn);

   return 0;
}
