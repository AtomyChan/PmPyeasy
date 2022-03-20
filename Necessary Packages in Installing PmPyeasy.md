# Necessary Packages in Installing PmPyeasy

## Prerequisites

0. install anaconda


   refer to: https://docs.anaconda.com/anaconda/install/linux/


1. Install iraf27 environment

   ```shell
   conda config --add channels http://ssb.stsci.edu/astroconda
   conda create -n iraf27 python=2.7 iraf-all pyraf-all stsci && source activate iraf27
   ```

2. Possible necessary packages

   ```shell
   pip install photutils ccdproc==1.3.0 wget astroquery==0.3.9 astroobs ephem
   ```
**Notes**``:
For MacOS, make sure Xcode work.

## A Quick Start

The most fundamental class is `photcore.photometry`.

```python
from photcore import photometry
```

**Notes**:
0. Issue encounterd on ubuntu 18.04 system with iraf/pyraf:

   It is the same error as reported here: https://iraf.net/forum/viewtopic.php?showtopic=1469803
   The instruction below didn't solve the poster's problem but did sovlve mine!!
   $sudo apt-get update
   $sudo apt-get install libc6:i386 libz1:i386 libncurses5:i386 libbz2-1.0:i386 libuuid1:i386 libxcb1:i386 libxmu6:i386


1. If the `photcore` module has trouble importing `iraf.stsdas` with errors: undefined variable 'uparm' in string 'uparm$clecl.par'. You can first try by explicitly imoprt stsdas as follow:

   ```python
   from iraf import stsdas
   ```

   If no issue from above, then 'mkiraf' in the PmPyeasy folder. The reference is here: https://maravelias.info/2011/05/login-cl-iraf-pyraf/

2. If there is problem in importing pyds9 during loading photcore, try to import pyds9 explicitly first.
    ```python
    import pyds9
    import photcore
    ```

The `photometry` class is initialized by automatically reading `.fits` files from your data directory, and generates ... at a result directory. With the parameter `telescope='XXX'` it custmizes the nomenclature of headers for telescopes listed in `./conf/telescopes.info`.

```python
sn = photometry(sn=SOURCE_NAME, telescope='XXX', data_repository_dir=YOUR_DATA_DIR, result_dir=YOUR_RESULT_DIR, pipe=False)
```

Information of these sources are stored in a dictionary `sn.photometry_info`, but it is easier to check `sn.photometry_record_table`, a more human-readable table.

>name                 realimg                 flt       obstime       camera exptime bitpix   NX    NY  ... instmagerr  relmag relmagerr  magzpt magzpterr  calmag calmagerr  drop
> str10                  str100                str10      float64       str30  float64 int32  int32 int32 ...  float64   float64  float64  float64  float64  float64  float64  int32
>-------- ------------------------------------ ----- ------------------ ------ ------- ------ ----- ----- ... ---------- ------- --------- ------- --------- ------- --------- -----

### External Packages

- `DIAPL`

  1. Run `./install.csh`

**Notes** if /bin/tcsh not available, install with "sudo apt install tcsh" on ubuntu

  2. enter fwhm_ping2 folder and
     ```shell
     make
     mv fwhm ../DIAPL_BIN/fwhm_ping2
     ```
     Notes: fwhm function in fwhm_ping2 enable file-based input parameter

  3. Add `DIAPL-BIN` path to your `.bashrc/.zshrc`

     ```shell
     export PATH="XXX/PmPyeasy/extern/DIAPL/DIAPL_BIN:$PATH"
     ```

- `fitsh-0.9.1`

  1. Run `./configure`

  2. Run `make`

  3. Add `src` path to your `.bashrc/.zshrc`

     ```shell
     export PATH="XXX/PmPyeasy/extern/fitsh-0.9.1/src:$PATH"
     ```

     Notes: if fitsh-0.9.1 fails, try 0.9.4 https://fitsh.net/wiki/Releases. (0.9.1 failed on our new work station with ubuntu but 0.9.4 complied smoothly)

- `CFITSIO`
  1. download from https://heasarc.gsfc.nasa.gov/fitsio
  2. ./configure --prefix=/usr/local; make; make install

  Notes: --prefix=/usr/local require write permission to /usr/local
         if zlib is missing, sudo apt install zlib1g-dev


  Notes: one error encountered when using hotpants
	"...ERROR: Mismatch in the CFITSIO_SONAME value in the fitsio.h include file that was used to build the CFITSIO library, and the value in the include fiel that was used when compiling the application program: version used to build CFITSIO library = 5, version included by the application program =9"
	Issue identified: apt list | grep cfitsio --> libcfitsio5 was installed and was used to compile hotpants
	Solution: apt remove libcfitsio5 and recompile hotpants


- `CCDproc`
  1. Run `make`

  Notes: modify the path to the HEASARC CFITSIO libraries in the `MAKEFILE` at `PmPyeasy/extern/CCDProc` (`CFITSIO=YOUR_CFITSIO_LIB`).
  if g++ not availale, install with "sudo apt install g++" on ubuntu


- `dophot`
  1. make dophot (see README file)

  Notes: gfortran is required; gfortran-7 works; gfortran-9 failed...
         Makefile cfitsio library "XXX/libcfitsio.so" works, "XXX/libcfitsio.a" failed
  Known issue: error while loading shared libraries libcfitsio.so.9: cannot open shared object file: No such file or directory
  Solution: Your library is a dynamic library. You need to tell the operating system where it can locate it at runtime.

    To do so, we will need to do those easy steps:

      (1) Find where the library is placed if you don't know it.

        $ sudo find / -name the_name_of_the_file.so

      (2) Check for the existence of the dynamic library path environment variable(LD_LIBRARY_PATH)

        $ echo $LD_LIBRARY_PATH

      if there is nothing to be displayed, add a default path value (or not if you wish to)

        $ LD_LIBRARY_PATH=/usr/local/lib

      (3) We add the desire path, export it and try the application.

        Note that the path should be the directory where the path.so.something is. So if path.so.something is in /my_library/path.so.something it should be :

          $ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/my_library/
          $ export LD_LIBRARY_PATH
          $ ./my_app

- `hotpants`
  1. make



