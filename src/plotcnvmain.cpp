/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
//#include <sys/types.h>
//#include <cstdlib>
using namespace std;

/* samtools header */
#include "samfunctions.h"

/**** user headers ****/
#include "wufunctions.h"
#include "plotcnv.h"


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int plot_usage() 
{
  cerr << "Usage:\n"
       << "   bamdepth -b BAM -r REGION -t TITLE -o OUTPUT\n"
       << "   bamdepth -b BAM -rf FILE -p FOLDER\n"
       << "\nOptions:\n"
       << "   REGION   in samtools format\n"
       << "   FILE     tab delimited BED format for each line, first 3 \n"
       << "            fields are required\n"
       << "            CHR BEGIN END COMMENT\n"
       << "   TITLE    optional, plot title, default to $REGION\n"
       << "   FOLDER   optional, folder for ps files, default to ./plots\n"
       << "   OUTPUT   optional, default to $FOLDER/$REGION.ps\n"
       << "   -q  INT  minimum map quality, INT=0\n"
       << "   -Q  INT  minimum base quality, INT=0\n"
       << "   -d  INT  plot d bases before and after region, INT=2*REGION\n"
       << "   -S       save intermediate dat and gnuplot files, default NOT\n"
       << "   -PE INT INT provide paired end insert and sd, otherwise\n"
       << "               calculate from bam file\n" 
       << "\nNote:\n"
       << "This program requires gnuplot to plot the figures. The figures are\n"
       << "saved in ./plots directory in postscript format. If ImageMagic is\n"
       << "available, the ps files will be converted to the png format.\n"
       << "GC contents are not adjusted.\n"
       << "\nExample:\n"
       << "   bamdepth -b $BAM -r  chr1:1000-2000 \n"
       << "   bamdepth -b $BAM -rf $CNV -p $PLOTFOLDER\n"
       << endl;
  
  return 0;
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int main(int argc, char* argv[])
{
  string mycommand="",RNAME;
  
  if ( gnuplot_version() < 0 ) {
    cerr << "gnuplot not found" << endl;
    return 0;
  }
  if ( argc < 2 ) return plot_usage();
  
  string bamfile="", cnvfile="", outfile="";
  string folder="rdplots", region="", title="", output="";   
  
  vector<string> ARGV(0);
  for(int i=0;i<argc;++i) ARGV.push_back(string(argv[i]));
  for(int i=1;i<(int)ARGV.size();++i) {
#define _next2 ARGV[i]=""; ARGV[i+1]=""; continue;
#define _next1 ARGV[i]=""; continue;
    if ( ARGV[i]=="-b" ) { bamfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-p" ) { folder=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-rf" ) { cnvfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-r" ) { region=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-t" ) { title=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-d" ) { plot::neighbor_length=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-o" ) { outfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-q" ) { plot::minq=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-Q" ) { plot::min_baseQ=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-S" ) { plot::save=true; _next1; }
    if ( ARGV[i]=="-PE" ) { 
      plot::pe_insert=atoi(ARGV[i+1].c_str());
      plot::pe_insert_sd=atoi(ARGV[i+2].c_str());
      ARGV[i+2]="";
      _next2; 
    }
  }
  
  plot_st iplot;
  vector<plot_st> plots(0);
  
  samfile_t *fp_in=samopen(bamfile.c_str(), "rb", 0);
  if ( ! fp_in ) { cerr << "Cannot open " << bamfile << endl; exit(0); }
  bam_index_t *bamidx=bam_index_load(bamfile.c_str()); 
  if ( !bamidx ) { cerr << "Cannot open idx file" << endl; exit(0); }
  
  
  if ( region!="" ) {
    iplot.region=region;
    iplot.title= title!="" ? title : region;
    size_t semicolon=region.find(':');
    if ( semicolon!=string::npos ) region[semicolon]='_';
    iplot.psfile= outfile!="" ? 
      outfile : 
      folder+"/"+region+".ps";
    plots.push_back(iplot);
  }
  
  if ( cnvfile != "" ) {
    ifstream INP(cnvfile.c_str());
    if ( !INP ) { cerr << "Cannot open " << cnvfile << endl; exit(0); }
    while ( !INP.eof() ) {
      string tmps, tmps1;
      getline(INP,tmps);
      istringstream iss(tmps);
      if (tmps.length() < 5 || tmps[0] == '#') continue;
      string chr,sbeg,send,comment="";
      if ( ! ( iss >> chr >> sbeg >> send ) ) {
	cerr << "expecting CHR START END in each line" << endl;
	continue;
      }
      if ( !(iss >> comment) ) comment="";
      if ( comment.length()>20 ) comment=comment.substr(0,20);
      replace(comment.begin(), comment.end(), ' ', '_');
      replace(comment.begin(), comment.end(), '\t', '_');
      iplot.region=chr+":"+sbeg+"-"+send;
      iplot.title=chr+":"+sbeg+"-"+send+" "+comment;
      iplot.psfile=folder+"/"+chr+"_"+sbeg+"-"+send+"_"+comment+".ps";
      plots.push_back(iplot);
    }
    INP.close();
  }    
  
  if ( plots.size()==0 ) {
    cerr << "no cnv found" << endl;
    plot_usage();
    goto CLOSEBAM;
  }
  
  for(size_t i=0; i<plots.size(); ++i) {
    if ( plots[i].region.find(":") == string::npos ) plots[i].region+=":";
    int ref=-1, beg=-1, end=0x7fffffff;
    int is_solved=bam_parse_region(fp_in->header, plots[i].region.c_str(), 
				   &ref, &beg, &end); 
    if ( is_solved<0 || ref<0 || ref>=fp_in->header->n_targets ||
	 beg<0 || end<0 ) {
      cerr << "Cannt resolve region " << plots[i].region << "\n"
	   << ref << "\t" << beg << "\t" << end << endl;
      goto CLOSEBAM;
    }
    if ( beg > end ) swap(beg, end);
    plots[i].tid=ref;
    plots[i].beg=beg;
    plots[i].end=end;
    plots[i].title+=" " + to_string(end-beg);
  }
  
  /*
    for(size_t i=0; i<plots.size(); ++i) {
    cerr << plots[i].region << "\n" << plots[i].title << "\n"
    << plots[i].tid << "\t" << plots[i].beg << "\t" << plots[i].end << "\n"
    << plots[i].psfile << "\n"
    << endl;
    }
  */
  
  plot_cnv_from_bam(plots, fp_in, bamidx);
  
  
 CLOSEBAM:
  if ( fp_in) samclose(fp_in);
  if ( bamidx) bam_index_destroy(bamidx);
  
  return 0;
} 
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
