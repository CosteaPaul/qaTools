/*                                                                                                                                                                                                                                          
    qaTools - Just more qa tools.
    Copyright (C) 2011  P. Costea(paul.igor.costea@scilifelab.se)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.                
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string>
#include "sam.h"

using namespace std;

#define EXIT_IF_NULL(P) \
  if (P == NULL) \
    return 1;

/**                                                                                   
 * Check if read is properly mapped                                                      
 * @return true if read mapped, false otherwise                                           
 */
static bool is_mapped(const bam1_core_t *core, int minQual)
{

  if ((core->flag&BAM_FUNMAP) || (core->qual < minQual)) {
    return false;
  }

  return true;
}

static int print_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: removeUnmapped \n");
    fprintf(stderr, "Contact: Paul Costea <paul.igor.costea@scilifelab.se>\n\n");
    fprintf(stderr, "Usage:   removeUnmapped [options] <in.bam/sam> <out.bam/sam>\n\n");
    fprintf(stderr, "Options: -q INT        quality threshold (strictly smaller than) [30]\n");
    fprintf(stderr, "         -i INT        minimum insert size [0]. Does not apply to single end libs\n");
    fprintf(stderr, "         -s            keep good quality unpaired reads\n");
    fprintf(stderr, "         -l            library is single end\n");
    fprintf(stderr, "         -k STR        write removed to file. If STR ends with .fastq, output will be in fastq format. Otherwise, SAM format\n");
    fprintf(stderr, "         -g STR        write good pairs to file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: Paired-end sam/bam must contain paires as subsequent entries. Position sorted files will not be processed!\n\n");
    return 1;
}

/** 
 * Open a .sam/.bam file.       
 * @returns NULL is open failed.
 */
samfile_t * open_alignment_file(std::string path, void* aux = NULL)
{
  samfile_t * fp = NULL;
  std::string flag = (aux==NULL) ? "r" : "w";
  if (path.substr(path.size()-3).compare("bam") == 0) {
    //BAM file!
    flag += "b";
  } else {
    flag += "h"; //Force write header!
  }
  if ((fp = samopen(path.c_str(), flag.c_str() , aux)) == 0) {
    fprintf(stderr, "Failed to open file %s\n", path.c_str());
  }
  return fp;
}


void print_bam_to_fastq(bam1_t *b, FILE* fastq, int cutoff=10000)
{
  char seq[200];
  char qual[200];
  //Get bam core.                                                                                                                                                                                                                             
  const bam1_core_t *core = &b->core;
  //Get the sequence                                                                                                                                                                                                                          
  int i=0;
  for (i; (i<core->l_qseq) && (i<cutoff); ++i) {
    seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)];
    //ADD 64 to bring it back to Illumina standard! silly hack, i'm sure it will cause problems later..
    qual[i] = bam1_qual(b)[i]+64;
  }
  seq[i] = '\0';
  qual[i] = '\0';
  //Print                                                                                                                                                                                                    
  fprintf(fastq,"@%s\n%s\n+\n%s\n",bam1_qname(b),seq,qual);
}

/**
 * Main of app
 */
int main(int argc, char *argv[])
{
    samfile_t *fp;
    samfile_t *out;
    samfile_t *out_rem = NULL;
    FILE * out_rem_f = NULL;
    FILE * out_rem_f1 = NULL;
    FILE * link = NULL;

    string removed_file = "";
    //Default to SAM
    bool rem_to_fastq = false;
    char* good_file = NULL;
    int minQual = 30;
    int minInsert = 0;
    bool keepSingle = false;
    bool isSingleEnd = false;
    int arg;
    //Get args
    while ((arg = getopt(argc, argv, "q:i:slk:g:")) >= 0) {
      switch (arg) {
      case 'q': minQual = atoi(optarg); break;
      case 'i': minInsert = atoi(optarg); break;
      case 'k':
	removed_file = optarg;
	if (removed_file.substr(removed_file.size()-5).compare("fastq") == 0) {
	  rem_to_fastq = true;
	}
	break;
      case 'g':
	good_file = new char[ strlen(optarg) ];
	strcpy(good_file, optarg);
	break;
      case 's':
	keepSingle = true;
	break;
      case 'l':
	isSingleEnd = true;
	break;
      default:
	fprintf(stderr,"Read wrong arguments! \n");
	break;
      }
    }
  
    //Check quality threshold
    if (minQual > 60) {
        //There is no qual higher than 60
      fprintf(stderr,"Minimum quality threshold too high\n");
        return 1;
    }

    if (argc-optind != 2) {
      print_usage();
      //Give up
      return 1;
    }

    fp = open_alignment_file(argv[optind]);
    EXIT_IF_NULL(fp);

    out = open_alignment_file(argv[optind+1],fp->header);
    EXIT_IF_NULL(out);

    if (!removed_file.empty()) {
    //Keep the removed ones somewhere
      if (rem_to_fastq) {
	if (isSingleEnd){//Only one output file!
	  out_rem_f = fopen(removed_file.c_str(),"w");
	  EXIT_IF_NULL(out_rem_f);
	} else {
	  out_rem_f = fopen((removed_file+"_1").c_str(),"w");
	  EXIT_IF_NULL(out_rem_f);
	  out_rem_f1 = fopen((removed_file+"_2").c_str(),"w");
	  EXIT_IF_NULL(out_rem_f1);
	}
      }else {
	out_rem = open_alignment_file(removed_file,fp->header);
      }
    removed_file = "";
  }

  if (good_file != NULL) {
    //Write good pair
    if ((link = fopen(good_file,"w")) == 0) {
      fprintf(stderr, "removeUnmapped: Filed to create output file %s\n", good_file);
      return 1;
    }
    delete[] good_file;
    good_file = NULL;
  }

  bam1_t *b = bam_init1();
  bam1_t *c = bam_init1();

  if (isSingleEnd) {
    while (samread(fp,b) >= 0) {
      if (is_mapped(&b->core, minQual)) {//Good, mapped!
	samwrite(out,b);
      } else if (out_rem_f || out_rem) {//Write removed?
	//Write removed here!
	if (rem_to_fastq) {
	  print_bam_to_fastq(b,out_rem_f);
	} else {
	  samwrite(out_rem,b);
	}
      }
    }

  } else {
    while (samread(fp, b) >= 0) {
      samread(fp,c);

      //Are these mapings paired
      if (b->core.isize == 0) {
	//Compute insert size
	b->core.isize = abs(b->core.pos - c->core.pos);
      }

      if (is_mapped(&b->core, minQual) && is_mapped(&c->core, minQual) && ( (abs(b->core.isize) >= minInsert) || (b->core.tid != c->core.tid) )) {
	string b_name = bam1_qname(b);
	string c_name = bam1_qname(c);
	if (b_name.find('/') != -1) {
	  //This is a bowtie generated file. removed "pair" information
	  b_name.erase(b_name.find('/'),2);
	  c_name.erase(c_name.find('/'),2);
	}
	samwrite(out,b);
	samwrite(out,c);
	
	if (b_name.compare(c_name) != 0)
	  fprintf(stderr,"Bam file is not properly ordered: %s - %s,\n or might your lib be single end? Use -l then.\n",b_name.c_str(),c_name.c_str());
	
	if (link != NULL) {
	  //Write to link file.
	  fprintf(link,"%s\t%s\t%d\t%d\n", b_name.c_str(), fp->header->target_name[b->core.tid], b->core.pos, b->core.pos + b->core.l_qseq);
	  fprintf(link,"%s\t%s\t%d\t%d\n", c_name.c_str(), fp->header->target_name[c->core.tid], c->core.pos, c->core.pos + c->core.l_qseq);
	}
      } else if (keepSingle && (b->core.flag&BAM_FUNMAP || c->core.flag&BAM_FUNMAP)){
	//      printf("Found one\n");
	//One of them is unmapped
	if (b->core.flag&BAM_FUNMAP && (b->core.qual >= (int)minQual/2)) {
	  //Write only c
	  samwrite(out,c);
	  fprintf(stderr,"Writing c: %s\t%s\t%d\t%d\t%d\n", bam1_qname(c),fp->header->target_name[c->core.tid], b->core.qual, b->core.flag&BAM_FUNMAP, c->core.flag&BAM_FUNMAP);
	} else if (c->core.flag&BAM_FUNMAP && (c->core.qual >= (int)minQual/2)) {
	  //Write only b
	  fprintf(stderr,"Writing b: %s\t%s\t%d\t%d\t%d\n", bam1_qname(b), fp->header->target_name[b->core.tid] ,c->core.qual, b->core.flag&BAM_FUNMAP, c->core.flag&BAM_FUNMAP);
	  samwrite(out,b);
	}
      } else if (out_rem || out_rem_f) {
	if (rem_to_fastq) {
	  print_bam_to_fastq(b,out_rem_f);
	  print_bam_to_fastq(c,out_rem_f1);
	}else {
	  samwrite(out_rem, b);
	  samwrite(out_rem, c);
	}
      }
    }
  }

  bam_destroy1(b);
  bam_destroy1(c);
  
  samclose(fp);
  samclose(out);
  if (out_rem)
    samclose(out_rem);
  else if (out_rem_f) {
    if (out_rem_f != NULL)
      fclose(out_rem_f);
    if (out_rem_f1 != NULL)
      fclose(out_rem_f1);
  }
  if (link != NULL)
    fclose(link);
  return 0;
}
	
      
