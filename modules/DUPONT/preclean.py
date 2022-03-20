PmPyeasy.iraf.hedit.unlearn()
PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'OBSERVAT', '(@"SITENAME")', add='yes', addonly='yes', verify='no')

#PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'ccdsec,datasec', delete='yes', verify='no', update='yes')
#PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'EPOCH', "(EQUINOX)", add='yes', addonly='yes', verify='no')

PmPyeasy.iraf.ccdproc.unlearn()
PmPyeasy.iraf.hedit(PmPyeasy.WorkDir+'*_SQ???.fits', 'ccdsec,datasec', delete='yes', verify='no', update='yes')
