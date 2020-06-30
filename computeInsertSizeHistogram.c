/* 
    qaTools - Just more qa tools.
    Copyright (C) 2011  P. Costea(paul.igor.costea@embl.de)

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
#include <getopt.h>
#include "sam.h"

using namespace std;

#define MIN(x,y) \
  (x < y) ? (x) : (y)

#define MAX(x,y) \
  (x > y) ? (x) : (y)

/**                                                                                   
 * Check if read is properly mapped                                                      
 * @return true if read mapped, false otherwise                                           
 */
static bool is_mapped(const bam1_core_t *core, int minQual)
{

  if ((core->flag&BAM_FUNMAP) || (int(core->qual) < minQual)) {
    return false;
  }

  return true;
}

static int print_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: computeInsertSizeHistogram \n");
    fprintf(stderr, "Version: 1.2\n");
    fprintf(stderr, "Contact: Paul Costea <paul.igor.costea@embl.de>\n\n");
    fprintf(stderr, "Usage:   computeInsertSizeHistogram [options] <in.bam/sam> <out.hist>\n\n");
    fprintf(stderr, "Options: -q INT        minimum quality mapping to consider in counting distribution [60]\n");
    fprintf(stderr, "         -l INT        maximum insert size (size of distribution) [1000]\n");
    fprintf(stderr, "         -v            split histogram by FR,RF,T\n");
    fprintf(stderr, "         -s            assume file is sorted");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: Input file must contain paires as subsequent entries, unless -s was specified \n\n");
    return 1;
}

/**
 * Main of app
 */
int main(int argc, char *argv[])
{
    samfile_t *fp;
    FILE * out = NULL;

    //Minimum quality to consider mapping
    int minQual = 60;
    int max = 1000;
    bool sorted = false;
    bool verbose = false;
    int arg;
    //Get args
    while ((arg = getopt(argc, argv, "q:l:vs")) >= 0) {
      switch (arg) {
      case 'q': minQual = atoi(optarg); break;
      case 'l': max = atoi(optarg); break;
      case 's': sorted = true; break;
      case 'v': verbose = true; break;
      default:
    	  fprintf(stderr,"Read wrong arguments! \n");
    	  break;
      }
    }

    if (argc-optind != 2) {
      print_usage();
      //Give up
      return -1;
    }

    int32_t **hist = NULL;
    if (verbose) {
      //Need a bit more memory...
      hist = new int32_t*[3];
      hist[0] = new int32_t[max];
      hist[1] = new int32_t[max];
      hist[2] = new int32_t[max];
    } else {
      hist = new int32_t*[1];
      hist[0] = new int32_t[max];
      for (int i=0; i<max;++i) hist[0][i]=0;
    }

    string alignFile(argv[optind]);
    //Check if this is sam or bam file
    string flag = "r";
    if (alignFile.substr(alignFile.size()-3).compare("bam") == 0) {
      //BAM file!
      flag += "b";
    }

    if ((fp = samopen(alignFile.c_str(), flag.c_str() , 0)) == 0) {
      fprintf(stderr, "Fail to open file %s\n", alignFile.c_str());
      return 1;
    }
    
    if ((out = fopen(argv[optind+1], "w")) == 0) {
      fprintf(stderr, "Filed to create output file %s\n", argv[optind+1]);
      return 1;
    }


  bam1_t *b = bam_init1();
  bam1_t *c = bam_init1();

  //Test for proper "sorting" (we can afford to chuck the first two reads
  samread(fp, b);
  samread(fp, c);
  if (string(bam1_qname(b)).compare(string(bam1_qname(c))) == 0) {
	  if (sorted) {
		  fprintf(stderr,"File doesn't seem to be sorted...Don't use -s \n");
	  }
  } else {
	  if (!sorted) {
		  fprintf(stderr,"File seems to be sorted...Use -s \n");
	  }
  }

  while (samread(fp,b) >= 0) {
      if (!sorted) {
	samread(fp,c);
      }

      if (is_mapped(&b->core, minQual)) { //First is mapped...
	if (sorted) {//This is a sorted file! Forget the second one, it's not here
	  if (!(b->core.flag&BAM_FMUNMAP) //Mate is also mapped!
	      && (b->core.pos < b->core.mpos)) {//Count pair only once! Thus, do this only for leftmost in pair	    
	    int i_size = 0;
	    if (b->core.tid != b->core.mtid)
	      i_size = 1;
	    else
	      i_size = abs(b->core.pos - b->core.mpos)+b->core.l_qseq;
	    if (i_size < max) {
	      if (verbose) {
		if (b->core.tid != b->core.mtid) {//Forget this one...
		  ++hist[0][i_size];
		  continue;
		}
		//Check FR, RF, T
		//Get flag of the leftmost read!
		int flag = b->core.flag;
		if ( (flag&(BAM_FMREVERSE)) && !(flag&BAM_FREVERSE) ) {//ForwardReverse
		  ++hist[0][i_size];
		} else if ( (flag&BAM_FREVERSE) && !(flag&BAM_FMREVERSE) ){//ReverseForward!
		  ++hist[2][i_size];
		} else if (((flag&(BAM_FMREVERSE|BAM_FREVERSE))==(BAM_FMREVERSE|BAM_FREVERSE))
			   || ((flag&(BAM_FMREVERSE|BAM_FREVERSE))==0)) {//T                                                                                                                                                              
		  ++hist[1][i_size];
		} else {
		  fprintf(stderr,"Cannot clasify %s %d\n",bam1_qname(b),b->core.flag);
		}
	      } else {
		++hist[0][i_size];
	      /*if (abs(b->core.isize)==1 && b->core.tid == b->core.mtid) {
		fprintf(stderr,"%s\n",bam1_qname(b));
		}*/
	      }
	    }
	  } 
	} else if (is_mapped(&c->core,minQual)) {//This is a consecutive entry list!
	  //Compute the "actual" insert size. SAM is being funny here.
	  int32_t i_size = 0;
	  if (b->core.tid != c->core.tid) //Diff chromosomes!
	    i_size = 1;
	  else {
	    int32_t left = MIN(b->core.pos,c->core.pos);
	    int32_t right = MAX(((c->core.pos)+(c->core.l_qseq)),((b->core.pos)+(b->core.l_qseq)));
	    i_size = right - left;
	    /*if (i_size < 300) {
	      fprintf(stderr,"%s\t%d\t%d\n",bam1_qname(b),i_size,b->core.isize);
	      }*/
	    /* Adding the length of the b sequence is not necessarily the best way to go,
	       because the other read may actually have a different length. However, this
	       should not prove to be a problem overall.
	     */
	  }
	  if (abs(i_size) < max) {
	    if (verbose) {
	      if (b->core.tid != c->core.tid) {//Forget this one...
		++hist[0][i_size];
		continue;
	      }
	      //Check FR, RF, T
	      //Get flag of the leftmost read!
	      int flag = (b->core.pos < c->core.pos) ? (b->core.flag) : (c->core.flag);
	      if ( (flag&(BAM_FMREVERSE)) && !(flag&BAM_FREVERSE) ) {//ForwardReverse
		++hist[0][i_size];
	      } else if ( (flag&BAM_FREVERSE) && !(flag&BAM_FMREVERSE) ){//ReverseForward!
		//fprintf(stderr,"%s %d\n",bam1_qname(b),b->core.flag);
		++hist[2][i_size];
	      } else if (((b->core.flag&(BAM_FMREVERSE|BAM_FREVERSE))==(BAM_FMREVERSE|BAM_FREVERSE)) 
		      || ((b->core.flag&(BAM_FMREVERSE|BAM_FREVERSE))==0)) {//T
		++hist[1][i_size];
	      } else {
		fprintf(stderr,"Cannot clasify %s %d\n",bam1_qname(b),b->core.flag);
	      }
	    } else {
	      ++hist[0][i_size];
	    }
	  }
	}
      }

  }
 
  bam_destroy1(b);
  bam_destroy1(c);

  if (verbose) {
    //Print header!
    fprintf(out,"FR\tT\tRF\n");
    for (int i=0; i<max; ++i) {
      fprintf(out,"%d\t%d\t%d\n",hist[0][i],hist[1][i],hist[2][i]);
    } 
  } else {
    for (int i=0; i<max; ++i)
      fprintf(out,"%d\n", hist[0][i]);
  }

  samclose(fp);
  fclose(out);

  return 0;
}
	
      
