#!/anaconda/bin/python

import sys,os,glob

gui = 'tk'
from numpy import *

try:
   from IPython.terminal.embed import InteractiveShellEmbed
except:
   try:
      from IPython.frontend.terminal.embed import InteractiveShellEmbed
   except:
      print "You don't have an up-to-date version of ipython."
      print "SNooPy requires version 0.12 or greater"
      sys.exit(1)

ipshell=InteractiveShellEmbed(banner1='Welcome to InteractiveShellEmbed')
ipshell.enable_pylab(gui=gui)
ipshell()
