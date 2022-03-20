Here I list the steps on how to perform image subtraction with the photometry pipeline:



(0) prepare the template image

(1) source detection: self.__get_star_number_rough_single_apphot_daofind

(2) match images: self.match_sources_on_two_images

(3) interpolate and resample image: self.interpolate_resample 














Some tips/notes or tricks:
(1) for Iowa images, the template images and science images for old targets have differecnt sizes (change of camera?)
The template images have size of 2048*2048 while the science images have size of 1536 x 1023

what to do? cut the images with smaller size to that of big size

it doesn't match well... what to do? do the matching and resampling again on the result product from previous procedure 
