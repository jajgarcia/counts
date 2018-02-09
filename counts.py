#!/usr/local/bin/python
#
# counts.py
#
# Print the counts in different bands of a Chandra HETG spectra
# (broad band, O, Fe, Ne, Mg, Si, S)
#
# Requires: xspec
#
import sys 
import pyfits
from xspec import *
from optparse import OptionParser
import os,os.path
import glob
#
# ------------------------------------------------------------------------------
#
# MAIN PROGRAM
#
#
#
version='0.2a'
date='- Wed Mar 12 13:18:24 EDT 2014 -'
author='Javier Garcia'
#
ul=[]
ul.append("usage: %prog [options] PREFIX")
ul.append("")
ul.append("Get total counts in different bands for a given observation")
ul.append("PREFIX can be a single PHA file or a group (e.g. *.pha)")
usage=""
for u in ul: usage+=u+'\n'

parser=OptionParser(usage=usage)
parser.add_option("-v","--version",action="store_true",dest="version",default=False,help="show version number")
parser.add_option("-i","--input",dest="inpfile",default="",help="specify spectrum file (default: meg_-1.pha.gz)")
parser.add_option("-o","--output",dest="outfile",default="",help="specify alternative output file (default: counts.out)")

(options,args)=parser.parse_args()

if options.version:
  print 'phacorr.py version:',version,date
  print 'Author:',author
  sys.exit()

if len(args) == 0:
  parser.print_help()
  sys.exit(0)

# If outfile is not define use the default
outfile=options.outfile
if outfile == "":
   outfile='counts.out'
f = open(outfile, 'w') 
output=[]
output.append('Observation          Broad-Band  O-Band  Fe-Band  Ne-Band  Mg-Band  Si-band  S-Band')

# If input is not define use the default
input=options.inpfile
if input == "":
   input='meg_-1.pha.gz'

# No chatter
Xset.chatter = 0


#-----
files=[]
allcounts=[]
for obspath in args:

  # Change dir
  print obspath
  os.chdir(obspath)

  # Get list of spectrum files
  files=glob.glob(input)

  # Get object name
  counts=[]
  hdulist = pyfits.open(files[0])
  object = hdulist[1].header['OBJECT']
  counts.append(object)

  cbb=0.
  cob=0.
  cfeb=0.
  cneb=0.
  cmgb=0.
  csib=0.
  csb=0.
  for specfile in files:
    # Check if specfile exist
    if not os.path.isfile(specfile):
      print 'Warning: spectrum file',specfile,'does not exist!'
      print 'Skiping...'

    else:

      # Load data
      s1 = Spectrum(specfile);

      # Exposure time
      et = s1.exposure
    
      # Ignore/notice data
      s1.ignore("**")
      s1.notice("0.49-6.19") # Just an example!!!

      # Count rate
      cr = s1.rate

      # Total counts
      cbb+=cr[0]*et

      # O-Band
      s1.ignore("**")
      s1.notice("0.51-0.59")   # Just an example !!!
      cr = s1.rate
      cob+=cr[0]*et

      # Fe-Band
      s1.ignore("**")
      s1.notice("0.67-0.72")   # Just an example !!!
      cr = s1.rate
      cfeb+=cr[0]*et

      # Ne-Band
      s1.ignore("**")
      s1.notice("0.82-0.95")   # Just an example !!!
      cr = s1.rate
      cneb+=cr[0]*et

      # Mg-Band
      s1.ignore("**")
      s1.notice("1.12-1.54")   # Just an example !!!
      cr = s1.rate
      cmgb+=cr[0]*et

      # Si-Band
      s1.ignore("**")
      s1.notice("1.77-2.06")   # Just an example !!!
      cr = s1.rate
      csib+=cr[0]*et

      # S-Band
      s1.ignore("**")
      s1.notice("2.4-2.6")   # Just an example !!!
      cr = s1.rate
      csb+=cr[0]*et


      # Unload data
      AllData -= s1

  # Output
  output.append(object+'   '+str(cbb)+'  '+str(cob)+'  '+str(cfeb)+' '+str(cneb)+'  '+str(cmgb)+'  '+str(csib)+'  '+str(csb))

  counts.append(cbb)
  counts.append(cob)
  counts.append(cfeb)
  counts.append(cneb)
  counts.append(cmgb)
  counts.append(csib)
  counts.append(csb)
  allcounts.append(counts)

  # Return to working dir
  os.chdir("../../")
#
# Write output
for o in output:
    f.write(o+'\n')

f.close
sys.exit()
# ------------------------------------------------------------------------------
