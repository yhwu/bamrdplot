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

/**** user headers ****/
#include "pstream.h"  // http://pstreams.sourceforge.net/
#include "wu2.h"
#include "wufunctions.h"
#include "alglibinterface.h"
/**** user headers ****/

#include "sam.h"
#include "samfunctions.h"

#include "plotcnv.h"

int plot::pts=30000;
int plot::count=0;
int plot::neighbor_length=0;
string plot::format="ps";
string plot::term="png xffffff x222222";
string plot::convert="convert -limit thread 1 -limit area 256MB -limit disk 512MB -density 72 -rotate 90 -background white -render -antialias -flatten ";
int plot::minq=0;
int plot::min_baseQ=0;
int plot::bam_tid=-1;
int plot::pe_insert=-1;
int plot::pe_insert_sd=-1;
bool plot::save=false;

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
std::string exec(const char* cmd) 
{
  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[1024];
  std::string result = "";
  while( !feof(pipe) ) {
    if(fgets(buffer, 128, pipe) != NULL) result+=buffer;
  }
  pclose(pipe);
  return result;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
float gnuplot_version(){
  string tmps="gnuplot -V 2>/dev/null | cut -d' ' -f2 ";
  float v=-1.0;
  string tmp1=exec( tmps.c_str() );
  istringstream iss(tmp1);
  
  if ( iss >> v ) return v;
  else {
    cerr << "gnuplot not found" << endl;
    return(-1.0);
  }
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_void(plot_st icnv, Array<int>& RD)
{
  if ( plot::format=="ps" ) plot::term="postscript color enhanced solid";
  if ( plot::format=="eps" ) plot::term="postscript eps enhanced solid";
  if ( plot::format=="png" ) plot::term="png xffffff x222222";
  
  ofstream FDAT(icnv.datafile.c_str());
  if ( !FDAT ) {
    cerr << "file " << icnv.datafile << " can't be created" << "\n"
	 << "plot_void()" << endl;
    return;
  }
  
  FDAT << "#" << icnv.title << endl;
  
  // dat file for gnuplot
  // data block 0
  double step=(double)(icnv.pend-icnv.pbeg+1)/(double)plot::pts;
  if ( step < 1.0 ) step=1.0;
  for (double ir=(double)icnv.pbeg+0.01; 
       ir<(double)icnv.pend+0.011; 
       ir+=step) {
    int i=(int) ir;
    int ic=i-1;
    if ( ic<0 || ic>=RD.size() ) FDAT << i << "\t0\n";
    else FDAT << i << "\t" << RD[ic] << "\n";
  }
  FDAT.close();

  icnv.title+=" Cross border";
  
  // mark xtics
  string xtics="("+to_string(icnv.pbeg)+" 0";
  xtics +=  ", " + to_string( (icnv.pbeg+icnv.beg+1)/2 ) + " 0";
  xtics +=  ", " + to_string( icnv.beg+1 ) + " 0";
  if ( icnv.end!=icnv.beg+1 ) {
    if ( (icnv.end-icnv.beg) < (icnv.pend-icnv.pbeg)/10 ) 
      xtics +=  ", \"\" " + to_string( icnv.end ) + " 0";
    else xtics +=  ", " + to_string( icnv.end ) + " 0";
  }
  xtics +=  ", " + to_string( (icnv.pend+icnv.end)/2 ) + " 0";
  xtics +=  ", " + to_string( icnv.pend ) + " 0";
  xtics +=  ")";
  
  ofstream FGP( icnv.gpfile.c_str() );
  string offset="offset";
  
  FGP << "#f=\"" << icnv.datafile << "\"" << endl
      << "#info=\"" <<  icnv.title << " \"" << endl
      << "set terminal " << plot::term << endl
      << "set output \"" << icnv.psfile << "\"" << endl
      << "set title \"" << icnv.title << "\"" << offset << " 0,-0.5" << endl
      << "set xtics font \"Helvetica bold,18\"" << endl
      << "set ytics font \"Helvetica bold,18\"" << endl
      << "set xtics " << xtics << "\n"
      << "set ytics nomirror" << endl
      << "set format x \"%.0f\" " << endl
      << "plot \""
      << icnv.datafile << "\" in 0 u 1:2 w p pt 7 ps 0.5 lt 3 not\n";
  FGP << "set output" << endl
      << "quit"
      << endl;
  FGP.close();
  
  string plotcommand="gnuplot < " + icnv.gpfile;
  redi::ipstream proc(plotcommand);
  string tmps;
  while (getline(proc, tmps))  std::cerr << tmps << endl;
  proc.close();
  //  system(plotcommand.c_str());
  
  if ( ! plot::save ) {
    remove(icnv.datafile.c_str());
    remove(icnv.gpfile.c_str());
  }
  
  return;  
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_icnv(plot_st icnv, 
	       double q0, 
	       Array<int>& RD, 
	       vector<pair_st>& plotpairs)
{
  int pts=plot::pts;
  
  if ( icnv.beg>icnv.end ) swap(icnv.beg, icnv.end) ;
  
  if ( icnv.beg > RD.size() || icnv.end >= RD.size() ) {
    cerr << "CNV exceeds reference length" << endl;
    plot_void(icnv, RD);
    return;
  }
  if ( icnv.pbeg<0 ) icnv.pbeg=0;
  if ( icnv.pend>RD.size()-1 ) icnv.pend=RD.size()-1;
  
  if ( plot::format=="ps" ) plot::term="postscript color enhanced solid";
  if ( plot::format=="eps" ) plot::term="postscript eps enhanced solid";
  if ( plot::format=="png" ) plot::term="png xffffff x222222";
  
  // sample points if there are too many
  int i1=icnv.pbeg;
  int i2=icnv.pend;
  if ( i1<1 ) i1=1;
  if ( i1>RD.size()-1 ) i1=RD.size()-1;
  if ( i2>RD.size()-1 ) i2=RD.size()-1;
  
  int i1c=i1-1;
  int i2c=i2-1;
  if ( i1c<0 ) i1c=0;
  if ( i1c>=RD.size() ) i1c=RD.size()-1;
  if ( i2c<0 ) i2c=0;
  if ( i2c>=RD.size() ) i2c=RD.size()-1;
  int c1=icnv.beg-1;
  int c2=icnv.end-1;
  if ( c1<0 ) c1=0;
  if ( c1>=RD.size() ) c1=RD.size()-1;
  if ( c2>=RD.size() ) c2=RD.size()-1;
  
  if ( i2c-i1c != i2-i1 ) {
    cerr << "coordinate converting error\n"
	 << i1 << "\t" << i1c << "\n"
	 << i2 << "\t" << i2c << endl;
  }
  if ( c2-c1 != icnv.end-icnv.beg ) {
    cerr << "coordinate converting error\n"
	 << icnv.beg << "\t" << c1 << "\n"
	 << icnv.end << "\t" << c2 << endl;
  }

  Array<int> ref(abs(icnv.beg-i1+i2-icnv.end));
  int i,k;
  double sum=0;
  for(i=i1, k=0; i<icnv.beg; ++i,++k) {
    int ic=i-1;
    if ( ic>=RD.size() ) break;
    if ( ic<0 || ic>=RD.size() ) ref[k]=0;
    else ref[k]=RD[ic];
    sum+=ref[k];
  }
  for(i=icnv.end+1; i<=i2; ++i,++k) {
    int ic=i-1;
    if ( ic>=RD.size() ) break;
    if ( ic<0 || ic>=RD.size() ) ref[k]=0;
    else ref[k]=RD[ic];
    sum+=ref[k];
  }
  if ( k!=ref.size() ) { 
    cerr << "size error\n" << "plot_icnv()" << endl;
    return;
  }
  
  int RDmax=0;
  for(int i=icnv.pbeg; i<icnv.pend; ++i) if ( RD[i]>RDmax ) RDmax=RD[i];
  
  // double refmed=sum/double(k);
  
  sum=0.0;
  for (i=c1; i<=c2 && i<RD.size(); ++i) sum+=RD[i];
  // double cnvmed=sum/double(c2-c1+1);
  
  // get basic statistics
  double icnv_ref_med,icnv_ref_lqt,icnv_ref_uqt,icnv_cnv_med,icnv_cnv_lqt,icnv_cnv_uqt;
  int n=c2-c1+1;
  if (n<0) {
    cerr << c1 << "\t" << c2 << "\n"
	 << "Cannot plot it" << endl
      	 << "plot_icnv()" << endl;
    return;
  }
  icnv_cnv_med=_median(&RD[c1],n);
  icnv_cnv_lqt=_lowerquartile(&RD[c1],n);
  icnv_cnv_uqt=_upperquartile(&RD[c1],n);
  
  Array<float> refrunmean(ref.size());
  Array<int> refrunmed(ref.size());
  int d=abs(icnv.end-icnv.beg)/2+1;
  int band=d+((d+1)%2)*1;
  if ( band < 50 ) band=50;
  if ( band > ref.size() ) band=ref.size()/2+((ref.size()/2+1)%2)*1;
  while ( band > ref.size() ) band-=2;
  if ( band<0 ) {
    cerr << "negative band size " << band << "\t" << ref.size() << "\n"
	 << "plot_icnv()" << endl;
    return;
  }
  runmean(ref, refrunmean, ref.size(), band, 1);
  icnv_ref_med=_median(&(ref[0]),ref.size());
  icnv_ref_lqt=_lowerquartile(&(ref[0]),ref.size());
  icnv_ref_uqt=_upperquartile(&(ref[0]),ref.size());
  
  // check range of plot
  int ymax=icnv_ref_med*2.25;
  // double y2max=(double)ymax / refmed / 2.0;
  if ( icnv_cnv_med < icnv_ref_med ) ymax=icnv_ref_med*2;
  else { ymax= icnv_cnv_uqt + 1.5*(icnv_ref_uqt-icnv_ref_lqt); }
  if ( ymax>RDmax ) ymax=RDmax+RDmax/10;
  ymax=RDmax+RDmax/20;
  if ( ymax<2 ) ymax=2;
  //ymax=plot_uqt+iqr*2;
  //if ( ymax>RDmax ) ymax=RDmax+RDmax/10;
  
  double r_upper=0;
  for(int i=(icnv.end+icnv.pend)/2; i<icnv.pend; ++i)  
    r_upper += ( RD[i]*2 > ymax );
  r_upper /= (icnv.pend - (icnv.end+icnv.pend)/2 + 0.05 );
  
  double l_upper=0;
  for(int i=icnv.pbeg; i<(icnv.beg+icnv.pbeg)/2; ++i)  
    l_upper += ( RD[i]*2 > ymax );
  l_upper /= ((icnv.beg+icnv.pbeg)/2 - icnv.pbeg + 0.05 );
  
  // mark xtics
  string xtics="("+to_string(icnv.pbeg)+" 0";
  xtics +=  ", " + to_string( (icnv.pbeg+icnv.beg+1)/2 ) + " 0";
  xtics +=  ", " + to_string( icnv.beg+1 ) + " 0";
  if ( icnv.end!=icnv.beg+1 ) {
    if ( (icnv.end-icnv.beg) < (icnv.pend-icnv.pbeg)/10 ) 
      xtics +=  ", \"\" " + to_string( icnv.end ) + " 0";
    else xtics +=  ", " + to_string( icnv.end ) + " 0";
  }
  xtics +=  ", " + to_string( (icnv.pend+icnv.end)/2 ) + " 0";
  xtics +=  ", " + to_string( icnv.pend ) + " 0";
  xtics +=  ")";
  
  ofstream FDAT(icnv.datafile.c_str());
  if ( !FDAT ) {
    cerr << "file " << icnv.datafile << " can't be created\n"
	 << "plot_cnv()" << endl;
    exit(0);
  }
  
  FDAT << "#" << icnv.title << '\n';
  
  // dat file for gnuplot
  // data block 0
  // before CNV  : pos RD refrunmean //
  double step=(double)(i2-i1+1)/(double)pts;
  if ( step < 1.0 ) step=1.0;
  // before cnv
  for (double ir=(double)icnv.pbeg+0.00001; ir<(double)icnv.beg+0.000011; ir+=step) {
    int i=(int) ir;
    if ( i<0 || i>=RD.size() ) FDAT << i+1 << "\t0\tNaN" << "\n";
    else FDAT << i+1 << "\t" << RD[i] << "\t" << refrunmean[i-i1] << "\n";
  }
  FDAT << "\n\n";
  
  // data block 1
  // cnv  : pos RD NaN //
  for (double ir=(double)icnv.beg+0.00001; ir<(double)icnv.end+0.000011; ir+=step) {
    int i=(int) ir;
    if ( i<0 || i>=RD.size() ) FDAT << i+1 << "\t0\tNaN" << "\n";
    else FDAT << i+1 << "\t" << RD[i] << "\t" << "NaN" << "\n";
  }
  FDAT << "\n\n";
  
  // data block 2
  // after cnv  : pos RD refrunmean //
  d=icnv.end-icnv.beg+1;
  for (double ir=(double)icnv.end+1.00001; ir<=(double)i2+0.000011; ir+=step) {
    int i=(int) ir;
    if ( i<0 || i>=RD.size() ) FDAT << i << "\t0\tNaN" << "\n";
    else FDAT << i+1 << "\t" << RD[i] << "\t" << refrunmean[i-i1-d] << "\n";
  }
  FDAT << "\n\n";
  
  // data block 3
  // dotted line to connect runmean //
  int imid1=icnv.beg-1-i1;
  int imid2=icnv.end+1-i1-d;
  if ( imid1<0 ) imid1=0;
  if ( imid2>=refrunmean.size() ) imid2=refrunmean.size()-1;
  FDAT << icnv.beg+1 << "\t" << refrunmean[imid1] << "\n"
       << icnv.end << "\t" << refrunmean[imid2] << "\n";
  FDAT << "\n\n";
  
  // data block 4
  // line for median of CNV //
  FDAT << icnv.beg << "\t" << icnv_cnv_med << "\n"
       << icnv.end << "\t" << icnv_cnv_med << "\n";
  FDAT << "\n\n";
  
  // data block 5
  // lower quartile of CNV //
  FDAT << icnv.beg << "\t" << icnv_cnv_lqt << "\n"
       << icnv.end << "\t" << icnv_cnv_lqt << "\n";
  FDAT << "\n\n";
  
  // data block 6
  // upper quartile of CNV //
  FDAT << icnv.beg << "\t" << icnv_cnv_uqt << "\n"
       << icnv.end << "\t" << icnv_cnv_uqt << "\n";
  FDAT << "\n\n";
  
  // data block 7
  // line for median of neighbor //
  FDAT << i1 << "\t" << icnv_ref_med << "\n"
       << i2 << "\t" << icnv_ref_med << "\n";
  FDAT << "\n\n";
  
  // data block 8
  // line for lower quartile of neighbor //
  FDAT << i1 << "\t" << icnv_ref_lqt << "\n"
       << i2 << "\t" << icnv_ref_lqt << "\n";
  FDAT << "\n\n";
  
  // data block 9
  // line for upper of neighbor //
  FDAT << i1 << "\t" << icnv_ref_uqt << "\n"
       << i2 << "\t" << icnv_ref_uqt << "\n";
  FDAT << "\n\n";
  
  // data block 10
  // nothing //
  FDAT << "NaN\tNaN\n"
       << "NaN\tNaN\n";
  FDAT << "\n\n";
  
  // data block 11
  // pairs //
  double base,top, dy=ymax/50;
  if ( icnv_cnv_med < icnv_ref_med ) {
    base=icnv_cnv_uqt+(RDmax-icnv_cnv_uqt)/10;
    if ( base < ymax/10 ) base=ymax/10;
    top=RDmax;
  }
  else {
    base=icnv_cnv_lqt-RDmax/20;
    top=RDmax/10;
    if ( base > ymax*2/3 ) base=ymax*2/3;
    if ( base < ymax/3 ) {
      base=icnv_cnv_uqt;
      top=RDmax;
    }
  }
  dy=(top-base)/(plotpairs.size()+1);
  if ( dy<ymax/100 ) {
    base=RDmax;
    top=RDmax/10;
  }
  dy=(top-base)/(plotpairs.size()+1);
  if ( abs(dy)*25 > ymax ) dy=( dy>0 ? 1:-1 )*ymax/25;
  
  for(size_t i=0; i<plotpairs.size(); ++i) {
    float y=base+(float)(plotpairs.size()-1-i)*dy ;
    FDAT << plotpairs[i].F1 << "\t" << y << "\t"
	 << plotpairs[i].F2-plotpairs[i].F1 << "\t" << 0 << "\t"
      // F1 |--| F2
	 << plotpairs[i].F2 << "\t" << y << "\t"
	 << plotpairs[i].R1-plotpairs[i].F2 << "\t" << 0 << "\t"
      // F2 --| R1
	 << plotpairs[i].R1 << "\t" << y << "\t"
	 << plotpairs[i].R2-plotpairs[i].R1 << "\t" << 0
      // R1 --> R2
	 << endl;
  }
  FDAT << "NaN\tNaN\n"
       << "NaN\tNaN\n";
  FDAT << "\n\n";
  FDAT.close(); // save data
  
  ofstream FGP(icnv.gpfile.c_str());
  float version=gnuplot_version();
  //version=4.0;
  string offset="offset";
  if ( version < 4.19 ) offset="";
  
  string comment="mapq=0: " + to_string( int(q0*100.0) ) + "%";
  comment+="\\n";
  comment+="minq="+to_string(plot::minq);
  comment+="\\n";
  comment+="minQ="+to_string(plot::min_baseQ);
  if ( plotpairs.size() > 0 ) {
    comment+="\\n";
    comment+="Pair="+to_string(plotpairs.size());
  }
  comment = "\"" + comment + "\"";
  
  int margin=(icnv.pend-icnv.pbeg)/40;
  
  // keys and labels
  if ( r_upper>0.5 ) FGP << "set key right bottom" << endl;
  FGP << "set label " << comment << " font \"Courier bold, 16\"" ;
  if ( l_upper>0.5 ) FGP << " at graph  0.05, graph  0.18 \n";
  else FGP << " at graph  0.05, graph  0.95 \n";
  
  // formats
  FGP << "f=\"" << icnv.datafile << "\" \n"
      << "set datafile missing \'NaN\' \n"
      << "info=\"" <<  icnv.title << " \" \n"
      << "set terminal " << plot::term << "\n"
      << "set output \"" << icnv.psfile << "\" \n"
      << "#set nokey \n"
      << "##set label 1 info at graph  0.25, graph  0.9 \n"
      << "set title info offset 0,-0.5 \n"
      << "set title info font \"Helvetica bold,18\" \n"
      << "set xtics font \"Helvetica bold,16\" \n"
      << "set ytics font \"Helvetica bold,16\" \n"
      << "set format x \"%.0f\" \n"
      << "set ytics nomirror \n"
      << "set xrange [" << i1-margin << ":" << i2+margin << "] \n"
      << "set xtics " << xtics << "\n"
      << "set yrange [0:" << ymax << "] \n";
  
  // plot RD data
  FGP << "plot \\" << endl
      << "f i 1 u 1:2 w p pt 7 ps 0.5 lt rgb \"red\" t \"REGION\" ,\\\n"
      << "f i 4 u 1:2 w l lt 1 lw 8 lc rgb \"cyan\" t \"REGION med\" ,\\\n"
      << "f i 5 u 1:2 w l lt 0 lw 8 lc rgb \"cyan\" t \"REGION 25%, 75%\" ,\\\n"
      << "f i 6 u 1:2 w l lt 0 lw 8 lc rgb \"cyan\" not ,\\\n"
      << "f i 0 u 1:2 w p pt 7 ps 0.5 lt rgb \"blue\" t \"Neighbor\" ,\\\n"
      << "f i 2 u 1:2 w p pt 7 ps 0.5 lt rgb \"blue\" not ,\\\n"
    //<< "f i 0:2:2 u 1:2 smooth bezier lw 2 lt rgb \"pink\" not, \\\n"
      << "f i 0 u 1:3 w l lt 1 lw 8 lc rgb \"green\" t \"Runmean\" ,\\\n"
      << "f i 2 u 1:3 w l lt 1 lw 8 lc rgb \"green\" not ,\\\n"
      << "f i 3 u 1:2 w l lt 0 lw 8 lc rgb \"green\" not ,\\\n"
      << icnv_ref_lqt << " w l lt 0 lw 8 lc rgb \"green\" t \"Neighbor 25%, 75%\" ,\\\n"
      << icnv_ref_uqt << " w l lt 0 lw 8 lc rgb \"green\" not ,\\\n";
  
  // plot pairs
  if ( plotpairs.size() > 0 ) {
    // F read |---|
    FGP << "f i 11 u 1:2:3:4  with vectors heads size screen 0.005,90 lw 2 lc rgb \"brown\" t \"Abnormal Pairs\" ,\\\n"
      // read gap --->
	<< "f i 11 u 5:6:7:8  with vectors head size screen 0.008,15 lw 2 lc rgb \"gray40\" not ,\\\n"
      // RC read |---|
	<<  "f i 11 u 9:10:11:12  with vectors heads size screen 0.005,90  lw 2 lc rgb \"brown\" not ,\\\n";
  }
  FGP << "NaN not \n";
  
  FGP << "set output\nquit" << endl;
  FGP.close();
  
  string plotcommand="gnuplot < " + icnv.gpfile ;
  
  //  redi::ipstream proc(plotcommand, redi::pstreams::pstderr);
  redi::ipstream proc(plotcommand);
  string tmps;
  while (getline(proc, tmps))  std::cerr << tmps << endl;
  proc.close();
  //  system(plotcommand.c_str());
  
  if ( ! plot::save ) {
    remove(icnv.datafile.c_str());
    remove(icnv.gpfile.c_str());
  }  

  return;  
}

// load a region
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void load_data_from_bam(samfile_t *fp_in, bam_index_t *bamidx, 
			plot_st icnv, 
			Array<int>& RD, 
			vector<pair_st>& pairs,
			double& q0) 
  
{
  q0=-1.0;
  pairs.clear();
  
  bam1_t *b=NULL;   
  bam_iter_t iter=0;
  b = bam_init1();
  
  if ( icnv.beg>icnv.end ) swap(icnv.beg, icnv.end);
  if ( icnv.pbeg>icnv.pend ) swap(icnv.pbeg, icnv.pend);
  
  if ( RD.size() != (int)fp_in->header->target_len[icnv.tid] ) {
    RD.resize(fp_in->header->target_len[icnv.tid]);
    RD.assign(0);
  }
  
  if ( icnv.pbeg<0 ) icnv.pbeg=0;
  if ( icnv.pend>RD.size() ) icnv.pend=RD.size();
  
  for(int i=icnv.pbeg; i<icnv.pend; ++i ) RD[i]=0;
  
  double count=0, q0_count=0, in_count=0;
  iter = bam_iter_query(bamidx, icnv.tid, icnv.pbeg, icnv.pend);
  while( bam_iter_read(fp_in->x.bam, iter, b)>0 ) {
    count+=1;
    if ( (size_t)count % 1000000 ==0 ) 
      cerr << "#processed " << (size_t)count/1000000 << "M reads" << endl;
    if ( (int)b->core.pos == 0 ) continue;
    if ( (int)b->core.tid < 0 ) continue;
    if ( (int)b->core.qual < plot::minq ) continue;
    if ( b->core.flag & BAM_DEF_MASK ) continue;
    
    POSCIGAR_st bm;
    resolve_cigar_pos(b, bm);  
    
    // get read depth; RD is 0 based
    for(int k=0;k<(int)bm.op.size();++k) {
      if ( bm.op[k] != BAM_CMATCH && bm.op[k] != BAM_CEQUAL ) continue;
      int p1=bm.cop[k]-1;
      int q1=bm.qop[k];
      for(size_t i=0; i<bm.nop[k] && p1<RD.size(); ++i, ++p1, ++q1)      
	if ( bam1_qual(b)[q1]>=plot::min_baseQ ) ++RD[p1];
    }
    
    // get q0
    if ( (b->core.pos+1 >=icnv.beg && b->core.pos < icnv.end) ||
	 (b->core.pos+b->core.l_qseq >=icnv.beg && b->core.pos+b->core.l_qseq < icnv.end) ) {
      in_count+=1;
      if ( b->core.qual == 0 ) q0_count+=1;
    }
    
    // get abnormal pairs
    if ( b->core.mtid != b->core.tid ) continue;
    if ( bool(b->core.flag&BAM_FREVERSE) ==
	 bool(b->core.flag&BAM_FMREVERSE) ) continue;
    if ( b->core.flag & BAM_FREVERSE ) continue; // only use F reads
    if ( b->core.isize > plot::pe_insert/3 && 
	 b->core.isize < plot::pe_insert+4*plot::pe_insert_sd ) continue;
    
    // 1 based position
    pair_st ipair;
    ipair.tid=icnv.tid; 
    ipair.F1=b->core.pos+1; 
    ipair.F2=bam_calend(&b->core, bam1_cigar(b))+1;
    ipair.R1=b->core.mpos+1;
    ipair.R2=b->core.pos+b->core.isize;
    pairs.push_back(ipair);
  }
  bam_destroy1(b);
  bam_iter_destroy(iter);
  
  q0=q0_count/(in_count+0.01);
  
  // check read depth
  double rdmean=0.0;
  for(int k=icnv.pbeg; k<icnv.pend; ++k) rdmean+=RD[k];
  rdmean/=(double)(icnv.pend-icnv.pbeg+0.000000001);
  
  // check read depth
  vector<bool> is_overlap(pairs.size(), true);
  for(size_t i=0; i<pairs.size(); ++i) {
    int F1=pairs[i].F1;
    int R2=pairs[i].R2;
    if ( F1>R2 ) swap(F1, R2);
    if ( max(F1, icnv.beg) > min(R2, icnv.end) ) is_overlap[i]=false;
    //if ( F1 < icnv.pbeg || R2>icnv.pend )  is_overlap[i]=false;
  }
  int k=0;
  for(size_t i=0; i<pairs.size(); ++i) {
    if ( is_overlap[i] ) pairs[k]=pairs[i];	
    k+=is_overlap[i];
  }
  pairs.erase(pairs.begin()+k, pairs.end());
  
  cerr << "Loaded length\t" << icnv.pend-icnv.pbeg << "\t" << rdmean << endl;
  cerr << "q0=" << q0 << endl;
  cerr << "Abnormal pairs:\t" << pairs.size() << endl;
  // cerr << "Note: GC content not adjusted" << endl;
  return;
}



/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_cnv_from_bam(vector<plot_st> &plots, 
		       samfile_t *fp_in, 
		       bam_index_t *bamidx)
{
  bool ImageMagick = true;
  string convertcom="convert -version | grep Image";
  plot::convert="convert -limit thread 1 -limit area 256MB -limit disk 512MB -density 72 -rotate 90 -background white -flatten ";
  if (exec(convertcom.c_str()).find("ImageMagick")==string::npos) {
    ImageMagick=false;
  }
  
  
  Array<int> RD(1);
  
  for(int i=0; i<(int)plots.size(); ++i) {
    plot_st icnv=plots[i];
    
    if ( plot::pe_insert<0 || plot::bam_tid !=icnv.tid ) {
      get_pairend_info(fp_in, bamidx, 
		       icnv.tid, plot::pe_insert, plot::pe_insert_sd);
      plot::bam_tid=icnv.tid;
    }
    
    int d=icnv.end-icnv.beg;
    if ( d < 10 ) d=20;
    else if ( d < 20 ) d=50;
    else if ( d < 100 ) d=100;
    if ( plot::neighbor_length>0 ) d=plot::neighbor_length;
    
    icnv.pbeg=max(icnv.beg-2*d,1);
    icnv.pend=plots[i].end+2*d;
    
    vector<pair_st> plotpairs(0);
    double q0=-1;
    
    
    cerr << icnv.region << "\n" << icnv.title << "\n" << icnv.psfile << endl;
    
    load_data_from_bam(fp_in, bamidx, icnv,
		       RD, plotpairs, q0 );
    
    
    icnv.datafile=plots[i].psfile+".dat";
    icnv.gpfile=plots[i].psfile+".gp";
    
    size_t p=plots[i].psfile.find('/');
    if ( p!=string::npos ) {
      redi::pstream proc("mkdir -p " + plots[i].psfile.substr(0,p) );
      proc.close();
    }
    
    plot_icnv(icnv, q0, RD, plotpairs);
    
    if ( ! ImageMagick ) continue;
    
    string png=plots[i].psfile.substr(0,icnv.psfile.size()-3)+".png";
    cerr << "converting: " << icnv.psfile << " to " << png << "\n" << endl;
    convertcom=plot::convert + icnv.psfile + " " + png ; 
    redi::pstream proc(convertcom);
    string tmps;
    while (getline(proc, tmps))  std::cerr << tmps << endl;
    proc.close();
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
