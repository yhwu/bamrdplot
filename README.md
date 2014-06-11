![The pair end algorithm](https://raw.githubusercontent.com/yhwu/bamrdplot/master/src/sampledup.png)
	
BAMRDPLOT plots the read depth and abnormal pairs in red for a given region. Regions twice the length of the given region before and after are plotted as well in blue. If the bam file is paired end, abnormal pairs are defined as 4 s.d. away from the average insert size. The pairs are plotted as an arrow with 3 segments

|---|-------->|---|

where the first segment is the Forward orientation read colored red, the second segment is colored gray for the unread part of the template, and the third segment colored red is for the Reverse Compliment read. So, something looks like 

|<--|---------|---|

indicates duplication. 

To install, simply download and make.

## Usage
```
[yhwu@node63 bamdepth]$ ./bamrdplot 
Usage:
   bamdepth -b BAM -r REGION -t TITLE -o OUTPUT
   bamdepth -b BAM -rf FILE -p FOLDER

Options:
   REGION   in samtools format
   FILE     tab delimited BED format for each line, first 3 
            fields are required
            CHR BEGIN END COMMENT
   TITLE    optional, plot title, default to $REGION
   FOLDER   optional, folder for ps files, default to ./plots
   OUTPUT   optional, default to $FOLDER/$REGION.ps
   -q  INT  minimum map quality, INT=0
   -Q  INT  minimum base quality, INT=0
   -d  INT  plot d bases before and after region, INT=2*REGION
   -S       save intermediate dat and gnuplot files, default NOT
   -PE INT INT provide paired end insert and sd, otherwise
               calculate from bam file

Note:
This program requires gnuplot to plot the figures. The figures are
saved in ./plots directory in postscript format. If ImageMagic is
available, the ps files will be converted to the png format.
GC contents are not adjusted.

Example:
   bamdepth -b $BAM -r  chr1:1000-2000 
   bamdepth -b $BAM -rf $CNV -p $PLOTFOLDER
```

