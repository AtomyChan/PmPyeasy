procedure bias4km (images)
#
# A script for overscan subtraction for the MDM 4k
# 
#   Paul Martini
#   2010 April 15
#
# Fixed typo in line 69 to allow for rectangular images -RS 19 May 2010
#

string images{ prompt="List of images to process " }
string suffix{"b", prompt="Suffix for bias-subtracted output "} 
struct *imglist

begin

  int i
  int overx, overy, naxis1, naxis2, binx, biny, quadx, quady, datax, datay
  string img, imgfile, d1, d2, d3, d4, d5, d6

# Expand the image template into a text file list

  imgfile = mktemp ("tmp$ctr")
  sections ( images, option="fullname", > imgfile )
  
  imglist = imgfile

  if (! defpac ("noao")) { noao }
  if (! defpac ("noao.imred")) { imred }
  if (! defpac ("noao.imred.bias")) { bias }

# Bias subtract each frame 
  
  while (fscan (imglist, img) != EOF ) {

    i = strlen (img)
    if ( substr (img, i-4, i) == ".fits")
      img = substr (img, 1, i-5)

    imgets(img, "i_naxis1", >&"dev$null") 
    naxis1 = int(imgets.value) 
    imgets(img, "i_naxis2", >&"dev$null") 
    naxis2 = int(imgets.value) 
    imgets(img, "overscnx", >&"dev$null") 
    overx = int(imgets.value) 
    imgets(img, "overscny", >&"dev$null") 
    overy = int(imgets.value) 
    imgets(img, "ccdxbin", >&"dev$null") 
    binx = int(imgets.value) 
    imgets(img, "ccdybin", >&"dev$null") 
    biny = int(imgets.value) 
    datax = naxis1-2*overx
    datay = naxis2-2*overy
    quadx = datax/2
    quady = datay/2
    overx = overx/binx
    overy = overy/biny

    print ("bias4k: processing "//img//".fits to output file "//img//suffix//".fits") 
    if(access(img//suffix//".fits")) {imdel(img//suffix)} 
    if(access(img//suffix//"1.fits")) {imdel(img//suffix//"1")} 
    if(access(img//suffix//"2.fits")) {imdel(img//suffix//"2")} 
    if(access(img//suffix//"3.fits")) {imdel(img//suffix//"3")} 
    if(access(img//suffix//"4.fits")) {imdel(img//suffix//"4")} 
 
    colbias(img, img//suffix//"1", bias="[1:"//overx//",1:"//quady//"]",trim="["//overx+1//":"//overx+quadx//",1:"//quady//"]",median-,interac-,function="legendre",order=4)
    colbias(img, img//suffix//"2", bias="[1:"//overx//","//quady+1//":"//2*quady//"]",trim="["//overx+1//":"//overx+quadx//","//quady+1//":"//2*quady//"]",median-,interac-,function="legendre",order=4)
    colbias(img, img//suffix//"3", bias="["//naxis1-overx+1//":"//naxis1//",1:"//quady//"]",trim="["//quadx+overx+1//":"//naxis1-overx//",1:"//quady//"]",median-,interac-,function="legendre",order=4)
    colbias(img, img//suffix//"4",bias="["//naxis1-overx+1//":"//naxis1//","//quady+1//":"//2*quady//"]",trim="["//quadx+overx+1//":"//naxis1-overx//","//quady+1//":"//2*quady//"]",median-,interac-,function="legendre",order=4)
    imarith(img//"["//1+overx//":"//naxis1-overx//",1:"//2*quady//"]", "*", 1.0, img//suffix, ver-) 
    imcopy(img//suffix//"1[1:"//quadx//",1:"//quady//"]", img//suffix//"[1:"//quadx//",1:"//quady//"]", ver-) 
    imcopy(img//suffix//"2[1:"//quadx//",1:"//quady//"]", img//suffix//"[1:"//quadx//","//quady+1//":"//2*quady//"]", ver-) 
    imcopy(img//suffix//"3[1:"//quadx//",1:"//quady//"]", img//suffix//"["//quadx+1//":"//2*quadx//",1:"//quady//"]", ver-) 
    imcopy(img//suffix//"4[1:"//quadx//",1:"//quady//"]", img//suffix//"["//quadx+1//":"//2*quadx//","//quady+1//":"//2*quady//"]", ver-)
    date | scan(d1, d2, d3, d4, d5, d6) 
    hedit (img//suffix, fields="BIAS4K", value=d1//" "//d2//" "//d3//" "//d4//" "//d5//" "//d6, add+, del-, ver-, show-, upd+) 
    hedit (img//suffix, fields="OVERSCAN", value="[1:32,*] and [4097:4128,*]", add+, del-, ver-, show-, upd+) 
    hedit (img//suffix, fields="CCDSEC", value="[33:4096,1:4064]", add+, del-, ver-, show-, upd+) 
  }

# Clean up

  delete (imgfile, ver-) 
  if(access(img//suffix//"1.fits")) {imdel(img//suffix//"1")} 
  if(access(img//suffix//"2.fits")) {imdel(img//suffix//"2")} 
  if(access(img//suffix//"3.fits")) {imdel(img//suffix//"3")} 
  if(access(img//suffix//"4.fits")) {imdel(img//suffix//"4")} 

end
