'''

***Eventually, this file should probably be merged into the io/pilatus.py file.***

This file is designed to read the log file produced by the PILATUS instrument, and then generate 
a single directory containing all of the information needed by the extract_iofqt() function in 
io/pilatus.py.  As input, you provide the following arguments, in order

sys.argv[1]	the time inverval between each image (enter '0' if you are going to customize later)
sys.argv[2]	the name of the log file (e.g., "0deg")
sys.argv[3]	the prefix common to all image files (e.g., "0deg_PIL1"
sys.argv[4:]	the list of SCAN NUMBERS YOU WANT TO COLLECT (e.g., 10 11 12 13 14 15)

Given these inputs, this script will parse the log file, copy the images in proper order into 
a new directory (called "output"), and then write a file containing the flux associated with 
each image.  No attempt to scale the values in the images is made in this script; that is 
handled in the function extract_iofqt() in io/pilatus.py.

This script also makes a number of assumptions that should be noted:
* The scan numbers must be listed in the order they appear in the file (i.e., numerical order)
* Images are renamed upon placement into the new directory, in increasing sequential order
* Consistent with the above, it is assumed that every image corresponds to a new ion beam time
* Again, NO ATTEMPT IS MADE TO AVERAGE ANY SCANS, and I ignored initial scans that needed to be averaged

Obviously, this file was written before we started collecting mutliple scans per data point, and 
to deal with that some extra code will need to be written.  In addition, I see that this script 
could probably be combined with the function extract_iofqt() in io/pilatus.py.

'''



import sys
import os
import shutil
from collections import deque

# time interval between each image; enter '0' if this is non-trivial and you need to customize later
simple_time_interval = float(sys.argv[1])
simple_time = 0.0

# open and read PILATUS log file
plog = sys.argv[2]
f = open(plog, "r")
lines = f.readlines()
f.close()

# prefix common to all image files
tiff_prefix = sys.argv[3]

# this is a list of all of the scan numbers we want to incorporate
target_scans = deque(sys.argv[4:])

# some placeholder information
lnum = -1
current_scan = None
global_reading = 0
fluxdata = []

# report
target_scan = int(target_scans.popleft())
print "looking for scan ", target_scan
sys.stdout.flush()




# Main Text Processing Segment
# -----------------------------------
'''
This is a moderately-complicated text processing routine that goes through the 
PILATUS log file
'''




done_with_all_scans = False
while (not done_with_all_scans):

  # increment the line number
  lnum += 1

  # if we are out of lines to read, then set the exit flag and continue
  if (lnum >= len(lines)):
    print "finished reading file ", plog, ". exiting"
    done_with_all_scans = True
    continue

  # otherwise, read the line and split it into words
  words = lines[lnum].split()

  # if the line is empty, just skip it
  if len(words) == 0:
    continue ;

  # if we are entering a new scan number, store it for future reference and continue
  if (words[0] == '#S'):
    current_scan = int(words[1])
    continue ;

  # keep skipping lines until we match a time series for a scan of interest
  if (words[0] == '#L') and (words[1] == 'Time') and (current_scan == target_scan):

    # okay, if we are in this loop it means that we found a match, 
    # and the *next line* will contain flux information.  But first we
    # need to do some setup work.

    # First, report that we found a match (for diagnostic purposes)
    print "found timeseries for scan ", current_scan
    sys.stdout.flush()

    # next, start a new increment for the fluxes read in this scan series
    local_reading = 0

    # Now begin an inner loop that gets all of the readings
    done_with_this_scan = False
    while (not done_with_this_scan):

      # move to the next line and split into words
      lnum += 1

      # if we have passed the end of the entire file, then we are *definitely* done with this scan, so set flag and exit
      if lnum >= len(lines):
        done_with_this_scan = True
        continue

      # otherwise, read the line and split into words
      vals = lines[lnum].split()

      # if the line has no entries or moves on to a new type of information, we are done.  Set flag and continue
      if len(vals) == 0 or vals[0] == "#C":
        done_with_this_scan = True
        continue

      # otherwise, then get the next flux value and add it to the list
      flux = int(float(vals[3]))	#  This is the value 'I2'
      inst_time = float(vals[0])	#  This is the value 'time'
      srcfile = "%s_%03d_%04d.tiff" % (tiff_prefix, target_scan, local_reading)
      fluxdata.append([simple_time, target_scan, local_reading, srcfile, flux, inst_time])

      # finally, increment the local and global flux reading counts
      local_reading += 1
      global_reading += 1
      simple_time += simple_time_interval


    # at this point, we have finished with one scan.  Report that
    print "finished with scan ", target_scan

    # If the deque of target scans is empty, we are done (successfully)
    if len(target_scans) == 0:
      done_with_all_scans = True
      continue

    # Otherwise, we should move on to the next target scan
    else:
      target_scan = int(target_scans.popleft())
      print "looking for scan ", target_scan
      sys.stdout.flush()



# Now we are totally done, and all images are transferred.  
# The only thing left to do is to write the flux data file
thefile = open('pilatus-config.csv', 'w')
for item in fluxdata:
  thefile.write("%f,%d,%d,%s,%d,%f\n" % (item[0], item[1], item[2], item[3], item[4], item[5]))
thefile.close()



      


  

