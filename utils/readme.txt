Python scripts used to perform photometry on images obtained from LCOGT telescopes


augement philosophy: result file and  log files during each task/subroutine are stored. When new images are added to the directory containing all original images, the pipeline will try to find the new added images and work on those new members. The control is realized by control message "new" and "add"


1. image_preprocess_pip.py
usage: python image_preprocess_pip.py target_name mode vebose

** Look up to filename_obstime_dump.dat for newly added images
** Get the image names 'data_dir' and get the obstime from fits files 
** preprocess the images: fix the saturated pixels which is set by a critical value; convert the data format from float32 to int16 which is required by the functions in diapl package

output file:

filename_obstime_dump.dat 
the names of all image files and corresponding observe time (JD) are recored in ascending sequence against JD.

filename_obstime_registration.txt
ascii version of the above file

fix_bad_pixels.log
log file for fixing saturated pixels task

flt2int.log
log file for flt2int task

2. registration_pip.py 
usage: python registration_pip.py target_name mode verbose
Note: if you want to display poor quality image and store them for later check, let verbose be 2

** Look up to image_info_dump.dat to find new images
** Get image information on filter, exposure time, background, fwhm
** The sfind task is used to find the sources on each image and the results are stored in data_dir/xys

output file: 
image_info_dump.dat
format: image_name  obstime  exptime     flts      bkg     fwhm    nstar  starnum 

image_info.txt 
ASCII version of image_info_dump.dat

good_image_info_dump.dat & poor_quality_image_info_dump.dat
The good or poor quality is defined by whether the backgound and fwhm information is retrived from the image

get_fwhm_bkg.log
get_star_number_rough.log

3. template_pip.py 
usage: python template_pip.py target_name filters verbose
Note: if you want to display template image and store them for later check, let verbose be 2

** Background value, fwhm and source number obtained in task 2 above are evaluated and the best image is used as the template
** The criteria on determining the 'best' is defined by parameter 'factors' which give bkg,fwhm and starnum weights

output file:
template.log 
The file name of  original image which is used as template for each filter is given 

3.5 mannual work 
add the coordinates(image coordinates or physical coordinate in equatorial FK5) of the object in interest to template.log 


4. photometry_aperture_pip.py
usage: python photometry_aperture_pip.py config_file verbose

4.5 mannual work (current using online app) 
calibration of the template images 

5. get_stdref_stars_pip.py
usage: usage: python get_stdref_stars_pip.py target_name ra dec image_size

6. stdcal_mag_pip.py
usage: usage: python stdcal_mag_pip.py target_name aperture_size filters vebose






