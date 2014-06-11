#ifndef _PLOT_H
#define _PLOT_H

class plot {
public: 
  static int pts;
  static int count;
  static int neighbor_length;
  static string format;
  static string term;
  static string convert;
  static string folder;
  static bool save;
  
  static int minq;
  static int min_baseQ;
  
  static int bam_tid;
  static int pe_insert;
  static int pe_insert_sd;
  
  ~plot(){};
};

struct plot_st {
  string region;
  string title;
  string psfile;
  string datafile;
  string gpfile;
  int tid;
  int beg;
  int end;
  int pbeg;
  int pend;
  plot_st(): region(""),
	     title(""),
	     psfile(""),
	     datafile(""),
	     gpfile(""),
	     tid(-1),
	     beg(-1),
	     end(-1),
	     pbeg(-1),
	     pend(-1){};
} ;

struct pair_st {
  int tid;
  int F1,F2,R1,R2;
  pair_st(): F1(-1),
	     F2(-1),
	     R1(-1),
	     R2(-1) {};
} ;


float gnuplot_version();
void plot_cnv_from_bam(vector<plot_st> &plots, samfile_t *fp_in, bam_index_t *bamidx);

#endif
