


# IRAF

IRAF, the Image Reduction and Analysis Facility, is one of the more historic way that astronomers reduce and analyze imaging and spectroscopic data.


[The Harvey Mudd College Intro to IRAF](1)
[Documentation from IRAF website](2) Use [FTP](3)


[1]: http://www.physics.hmc.edu/Astronomy/Imain.html
[2]: http://iraf.noao.edu/docs/
[3]: http://iraf.noao.edu/iraf/ftp/docs/






#PyRAF

Command language for running IRAF tasks that is based on the Python scripting language. As of October 1, 2019, STScI no longer supports PyRAF. Users will still be able to access an installation of PyRAF through [Astroconda](1). However, STScI will no longer answer help calls related to PyRAF installation, bugs, or other usage questions.

Known installation problems:
(1) AttributeError: Undefined IRAF task `chkupdate’ 
Solution: create new login.cl by mkiraf and search for 'chkupdate' and comment out it

Please read the [documentation](2) carefully before use this.

import iraf task

> from pyraf import iraf

Some tasks can import directly from iraf, for example:

 > from iraf import imcombine
 > from iraf import apphot

Some task can't be imported before its parent package is imported. For example, if you want to load task *colbias*, plase see below for how it works.

 > from iraf import imred
 > from iraf import bias
 > from iraf import colbias

[1]: https://astroconda.readthedocs.io/en/latest/
[2]: http://stsdas.stsci.edu/pyraf/doc.old/pyraf_guide/





#Useful Tasks

- noao.obsutil.psfmeasure
- 
