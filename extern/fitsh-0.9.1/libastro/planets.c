/*****************************************************************************/
/* ed_plan.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Calculates the heliocentric ecliptical Cartesian coordiants of planets.   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Written by Kospal, A. (kospal@szofi.elte.hu)	and Pal, A. 		     */
/* (apal@szofi.elte.hu), maintained by Pal, A.				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>

#include <astro/astro.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841968
#endif

#define		EQNX1900
/* #define		EQNX2000 */

#ifdef		EQNX1900

/************************ EQUINOX 1900.0, T_0=1900.0 *************************/
/* Meeus page 95-97 */

static double	merc[]={				/* L,a,e,i,o,O */

    178.179078 , 0.3870986, 0.20561421 , 7.002881   ,28.753753   ,47.145944,
 149474.07078  , 0.0      , 0.00002046 , 0.0018608  , 0.3702806  , 1.1852083,
      0.0003011, 0.0      ,-0.000000030,-0.0000183  , 0.0001208  , 0.0001739,
      0.0      , 0.0      , 0.0        , 0.0        , 0.0        , 0.0
};

static double	venu[]={
   342.767053 , 0.7233316, 0.00682069 , 3.393631  ,54.384186   ,75.779647,
 58519.21191  , 0.0      ,-0.00004774 , 0.0010058 , 0.5081861  , 0.8998500,
     0.0003097, 0.0      , 0.000000091,-0.0000010 ,-0.0013864  , 0.0004100,
     0.0      , 0.0      , 0.0        , 0.0       , 0.0        , 0.0
};

static double	mars[]={
   293.737334 , 1.5236883, 0.09331290 , 1.850333   ,285.431761   , 48.786442,
 19141.69551  , 0.0      , 0.000092064,-0.0006750  ,  1.0697667  ,  0.7709917,
     0.0003107, 0.0      ,-0.000000077, 0.0000126  ,  0.0001313  , -0.0000014,
     0.0      , 0.0      , 0.0        , 0.0        ,  0.00000414 , -0.00000533
};

static double	jupi[]={
  238.049257  , 5.202561, 0.04833475  , 1.308736 , 273.277558  , 99.443414,
 3036.301986  , 0.0     , 0.000164180 ,-0.0056961,   0.5994317 ,  1.0105300,
    0.0003347 , 0.0     ,-0.0000004676, 0.0000039,   0.00070405,  0.00035222,
   -0.00000165, 0.0     ,-0.0000000017, 0.0      ,   0.00000508, -0.00000851
};

static double	satu[]={
  266.564377 ,9.554747, 0.05589232   , 2.492519  ,338.307800  ,112.790414,
 1223.509884 ,0.0     ,-0.00034550   ,-0.0039189 ,  1.0852207 ,  0.8731951,
    0.0003245,0.0     ,-0.000000728  ,-0.00001549,  0.00097854, -0.00015218,
   -0.0000058,0.0     , 0.00000000074, 0.00000004,  0.00000992, -0.00000531
};


#else
#ifdef EQNX2000

/************************ EQUINOX 2000.0, T_0=1900.0 *************************/
/* Meeus page 95-97 + 100-101 */

static double	merc[]={				/* L,a,e,i,o,O */

    178.179078 , 0.3870986, 0.20561421 , 7.010678   ,28.839814   ,48.456876,
 149474.07078  , 0.0      , 0.00002046 ,-0.0059556  , 0.2842765  ,-0.1254715,
      0.0003011, 0.0      ,-0.000000030, 0.00000069 , 0.00007445 ,-0.00008844,
      0.0      , 0.0      , 0.0        ,-0.000000035, 0.000000043,-0.000000068
};

static double	venu[]={
   342.767053 , 0.7233316, 0.00682069 , 3.395459  ,54.602827   ,76.957740,
 58519.21191  , 0.0      ,-0.00004774 ,-0.0007913 , 0.2892764  ,-0.2776656,
     0.0003097, 0.0      , 0.000000091,-0.0000325 ,-0.00114464 ,-0.0001401,
     0.0      , 0.0      , 0.0        , 0.00000018,-0.000000794, 0.000000769
};

static double	mars[]={
   293.737334 , 1.5236883, 0.09331290 , 1.857866   ,285.762379   , 49.852347,
 19141.69551  , 0.0      , 0.000092064,-0.0081565  ,  0.7387251  , -0.2941821,
     0.0003107, 0.0      ,-0.000000077,-0.00002304 ,  0.00046556 , -0.00064344,
     0.0      , 0.0      , 0.0        ,-0.000000044,  0.000006939, -0.000008159
};

static double	jupi[]={
  238.049257  , 5.202561, 0.04833475  , 1.305288   , 273.829584   ,100.287838,
 3036.301986  , 0.0     , 0.000164180 ,-0.0022374  ,   0.0478404  ,  0.1659357,
    0.0003347 , 0.0     ,-0.0000004676, 0.00002942 ,  -0.00021857 ,  0.00096672,
   -0.00000165, 0.0     ,-0.0000000017, 0.000000127,   0.000008999, -0.00001246
};

static double	satu[]={
  266.564377 ,9.554747, 0.05589232   , 2.486204   ,338.571353   ,113.923406,
 1223.509884 ,0.0     ,-0.00034550   , 0.0024449  ,  0.8220515  , -0.2599254,
    0.0003245,0.0     ,-0.000000728  ,-0.00005017 ,  0.00070747 , -0.00018997,
   -0.0000058,0.0     , 0.00000000074, 0.000000002,  0.000006177, -0.000001589
};


#else

#error	"Invalid equinox"

#endif
#endif

/************ EQUINOX 2000.0, Uranus, Neptune, Pluto, T_0=2000.0 *************/

static double	uran[]={			/* L, a, e, i, ~o (!), O */
  313.23218, 19.19126, 0.0471677 , 0.76986, 170.96424, 74.22988,
  428.48549, 0.00152, -0.00019155,-0.00058, 0.3646, -0.4670565,
 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0
};

static double	nept[]={
  304.88003 , 30.068963, 0.00858587, 1.76917    , 44.97135, 131.72169  ,
 218.4581139, -0.001252, 0.00000251, -0.00101111, -0.23456, -0.04201389,
 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0
};


static double	plut[]={
238.92881, 39.48168677, 0.24880766, 17.14175, 224.06676, 110.30347,
145.20775,0,0,0,0,0,
0,0,0,0,0,0,
0,0,0,0,0,0 
};


/* nept: */
/* 30.06896348, 0.00858587, 1.76917, 131.72169,  44.97135, 304.88003, */
/* -0.00125196,  0.0000251,   -3.64,   -151.25,   -844.43, 786449.21  */

/* plut: */
/* 39.48168677, 0.24880766,17.14175, 110.30347, 224.06676, 238.92881, */
/* -0.00076912, 0.00006465,   11.07,    -37.33,   -132.25, 522747.90  */


static double pol_jupi[]={
 0.331364, -0.010281,-0.0004692,
 0.003228, -0.064436, 0.002075,
-0.003083, -0.000275, 0.000489,
 0.002472,  0,        0,
 0.013619,  0,        0,
 0.018472,  0,        0,
 0.006717,  0,        0,
 0.002775,  0,        0,
 0.007275, -0.001253, 0,
 0.006417,  0,        0,
 0.002439,  0,        0,
-0.033839, -0.001125, 0,
-0.003767,  0,        0,
-0.035681, -0.001208, 0,
-0.004261,  0,        0,
 0.002178,  0,        0,
-0.006333,  0.001161, 0,
-0.006675,  0,        0,
-0.002664,  0,        0,
-0.002572,  0,        0,
-0.003567,  0,        0,
 0.002094,  0,        0,
 0.003342,  0,        0,  /* perturbations in the mean longitude (Meeus page 106) */

     3606,  130, -43,
     1289, -580,   0,
    -6764,    0,   0,
    -1110,    0,   0,
     -224,    0,   0,
     -204,    0,   0,
     1284,  116,   0,
      188,    0,   0,
     1460,  130,   0,
      224,    0,   0,
     -817,    0,   0,
     6074,    0,   0,
      992,    0,   0,
      508,    0,   0,
      230,    0,   0,
      108,    0,   0,
     -956,  -73,   0,
      448,    0,   0,
      137,    0,   0,
     -997,  108,   0,
      480,    0,   0,
      148,    0,   0,
     -956,   99,   0,
      490,    0,   0,
      158,    0,   0,
      179,    0,   0,
     1024,   75,   0,
     -437,    0,   0,
     -132,    0,   0,  /* perturbatin in the eccentricity (Meeus page 106-107) */

 0.007192, -0.003147, 0,
-0.020428, -0.000675, 0.000197,
 0.007269,  0.000672, 0,
-0.004344,  0,        0,
 0.034036,  0,        0,
 0.005614,  0,        0,
 0.002964,  0,        0, 
 0.037761,  0,        0, 
 0.006158,  0,        0,
-0.006603,  0,        0,
-0.005356,  0,        0,
 0.002722,  0,        0,
 0.004483,  0,        0,
-0.002642,  0,        0,
 0.004403,  0,        0,
-0.002536,  0,        0,
 0.005547,  0,        0,
-0.002689,  0,        0, /* perturbations in the perihelion (Meeus page 107) */

 -263, 0, 0,
  205, 0, 0,
  693, 0, 0,
  312, 0, 0,
  147, 0, 0,
  299, 0, 0, 
  181, 0, 0, 
  204, 0, 0, 
  111, 0, 0,
 -337, 0, 0,
 -111, 0, 0 /* perturbations is the semimajor axis (Meeus page 108) */
};

static double pol_satu[]={

 -0.814181, 0.018150	,-0.016714,
 -0.010497, 0.160906	,-0.004100,
  0.007581, 0		, 0,
 -0.007986, 0		, 0,
 -0.148811, 0		, 0,
 -0.040786, 0		, 0,
 -0.015208, 0		, 0,
 -0.006339, 0		, 0,
 -0.006244, 0		, 0,
  0.008931, 0.002728	, 0,
 -0.016500, 0		, 0,
 -0.005775, 0		, 0,
  0.081344, 0.003206	, 0,
  0.015019, 0		, 0,
  0.085581, 0.002494	, 0,
  0.025328,-0.003117	, 0,
  0.014394, 0		, 0,
  0.006319, 0		, 0,
  0.006369, 0		, 0,
  0.009156, 0		, 0,
  0.007525, 0		, 0,
 -0.005236, 0		, 0,
 -0.007736, 0		, 0,
 -0.007528, 0		, 0,  /* perturbations in the mean longitude (Meeus page 109) */

 -7927, 2548, 91,
 13381, 1226,-253,
   248, -121, 0,
  -305,  -91, 0, 
   412,    0, 0,
 12415,    0, 0,
   390, -617, 0,
   165, -204, 0,
 26599,    0, 0,
 -4687,    0, 0,
 -1870,    0, 0,
  -821,    0, 0,
  -377,    0, 0,
   497,    0, 0,
   163, -611, 0,
-12696,    0, 0,
 -4200,    0, 0,
 -1503,    0, 0,
  -619,    0, 0,
  -268,    0, 0,
  -282,-1306, 0,
   -86,  230, 0,
   461,    0, 0,
  -350,    0, 0,
  2211, -286, 0,
 -2208,    0, 0,
  -568,    0, 0,
  -346,    0, 0,
 -2780, -222, 0,
  2022,  263, 0,
   248,    0, 0,
   242,    0, 0,
   467,    0, 0,
  -490,    0, 0,
 -2842, -279, 0,
   128,  226, 0,
   224,    0, 0,
 -1594,  282, 0,
  2162, -207, 0,
   561,    0, 0,
   343,    0, 0,
   469,    0, 0,
  -242,    0, 0,
  -205,    0, 0,
   262,    0, 0, 
   208,    0, 0, 
  -271,    0, 0, 
  -382,    0, 0, 
  -376,    0, 0,  /* petrurbations in the eccentricity (Meeus page 109-110) */

  0.077108, 0.007186	, -0.001533,
  0.045803,-0.014766	, -0.000536,
 -0.007075, 0		, 0,
 -0.075825, 0		, 0,
 -0.024839, 0		, 0,
 -0.008631, 0		, 0,
 -0.072586, 0		, 0,
 -0.150383, 0		, 0,
  0.026897, 0		, 0,
  0.010053, 0		, 0,
 -0.013597,-0.001719	, 0,
 -0.007742, 0.001517	, 0,
  0.013586,-0.001375	, 0,
 -0.013667, 0.001239	, 0,
  0.011981, 0		, 0,
  0.014861, 0.001136	, 0,
 -0.013064,-0.001628	, 0,  /* perturbations in the perihelion (Meeus page 111) */

   572,   0, 0,
  2933,   0, 0,
 33629,   0, 0,
 -3081,   0, 0,
 -1423,   0, 0,
  -671,   0, 0,
  -320,   0, 0,
  1098,   0, 0,
 -2812,   0, 0,
   688,   0, 0,
  -393,   0, 0,
  -228,   0, 0,
  2138,   0, 0,
  -999,   0, 0,
  -642,   0, 0,
  -325,   0, 0,
  -890,   0, 0,
  2206,   0, 0,
 -1590,   0, 0,
  -647,   0, 0,
  -344,   0, 0,
  2885,   0, 0,
  2172, 102, 0,
   296,   0, 0,
  -267,   0, 0,
  -778,   0, 0,
   495,   0, 0,
   250,   0, 0,
  -856,   0, 0,
   441,   0, 0,
   296,   0, 0,
   211,   0, 0,
  -427,   0, 0,
   398,   0, 0,
   344,   0, 0,
  -427,   0, 0 /* perturbations in the semimajor axis (Meeus page 111) */
};


/* 1900 */
/* static double	uran[]={
 244.197470  ,19.21814, 0.0463444  ,0.772464 ,98.071581  ,73.477111,
 429.863546  , 0.0    ,-0.00002658 ,0.0006253, 0.9857650 , 0.4986678,
   0.0003160 , 0.0    , 0.000000077,0.0000395,-0.0010745 , 0.0013117,
  -0.00000060, 0.0    , 0.0        ,0.0      ,-0.00000061, 0.0
};

static double	nept[]={
  84.457994  ,30.10957, 0.00890704 , 1.779242 ,276.045975   ,130.681389,
 219.885914  , 0.0    , 0.000006330,-0.0095436,  0.3256394  ,  1.0989350,
   0.0003205 , 0.0    ,-0.000000002,-0.0000091,  0.00014095 ,  0.00024987,
  -0.00000060, 0.0    , 0.0        , 0.0      ,  0.000004113, -0.000004718
}; */

/* 2000 */
/* static double	uran[]={
 244.197470  ,19.21814, 0.0463444  , 0.774950   ,99.021587   ,73.923501,
 429.863546  , 0.0    ,-0.00002658 ,-0.0017660  , 0.0337219  , 0.0545828,
   0.0003160 , 0.0    , 0.000000077,-0.00000027 ,-0.00049812 , 0.00042674,
  -0.00000060, 0.0    , 0.0        , 0.000000123, 0.000013904,-0.000014536
};

static double	nept[]={
  84.457994  ,30.10957, 0.00890704 , 1.769715   ,276.335328   ,131.788486,
 219.885914  , 0.0    , 0.000006330,-0.0000144  ,  0.0368127  , -0.0084187,
   0.0003205 , 0.0    ,-0.000000002,-0.00000227 ,  0.00003849 ,  0.00004428,
  -0.00000060, 0.0    , 0.0        , 0.000000018,  0.000002226, -0.000002858
};
*/
                              


static	void sum_jupiter_perturbations(	double *ax,double *ex,double *el,
					double *mm,double *om,double jc1)
{
 double	nu,nu2,pe,ku,es,ve,zeta,a,b,ex1,ax1;

	nu  =jc1/5+0.1;
	nu2 =nu*nu;

	pe  =237.47555+3034.9061*jc1;
	ku  =265.91650+1222.1139*jc1;
	es  =243.51721+ 428.4677*jc1;

	ve  =5*ku-2*pe;
	zeta=  ku-  pe;

 a=
	(pol_jupi[ 0]+pol_jupi[ 1]*nu+pol_jupi[ 2]*nu2)*sina(ve)+
 	(pol_jupi[ 3]+pol_jupi[ 4]*nu+pol_jupi[ 5]*nu2)*cosa(ve)+
	(pol_jupi[ 6]+pol_jupi[ 7]*nu+pol_jupi[ 8]*nu2)*sina(2*ve)+
	(pol_jupi[ 9])*sina(2*pe-6*ku+3*es)+
	(pol_jupi[12])*sina(zeta)+
	(pol_jupi[15])*sina(2*zeta)+
	(pol_jupi[18])*sina(3*zeta);
	a+=
	(pol_jupi[21])*sina(4*zeta)+
	(pol_jupi[24]+pol_jupi[25]*nu)*sina(zeta)*sina(ku)+
	(pol_jupi[27])*sina(2*zeta)*sina(ku)+
	(pol_jupi[30])*sina(3*zeta)*sina(ku)+
	(pol_jupi[33]+pol_jupi[34]*nu)*cosa(zeta)*sina(ku)+
	(pol_jupi[36])*cosa(2*zeta)*sina(ku);
	a+=
	(pol_jupi[39]+pol_jupi[40]*nu)*sina(zeta)*cosa(ku)+
	(pol_jupi[42])*sina(2*zeta)*cosa(ku)+
	(pol_jupi[45])*cosa(ku)+
	(pol_jupi[48]+pol_jupi[49]*nu)*cosa(zeta)*cosa(ku)+
	(pol_jupi[51])*cosa(2*zeta)*cosa(ku)+
	(pol_jupi[54])*cosa(3*zeta)*cosa(ku);
	a+=
	(pol_jupi[57])*sina(zeta)*sina(2*ku)+
	(pol_jupi[60])*sina(2*zeta)*sina(2*ku)+
	(pol_jupi[63])*cosa(zeta)*cosa(2*ku)+
	(pol_jupi[66])*cosa(2*zeta)*cosa(2*ku);

	ex1=
	(pol_jupi[ 69]+pol_jupi[ 70]*nu+pol_jupi[ 71]*nu2)*sina(ve)+
	(pol_jupi[ 72]+pol_jupi[ 73]*nu)*cosa(ve)+
	(pol_jupi[ 75])*sina(zeta)*sina(ku)+
	(pol_jupi[ 78])*sina(2*zeta)*sina(ku)+
	(pol_jupi[ 81])*sina(3*zeta)*sina(ku)+
	(pol_jupi[ 84])*sina(ku)+
	(pol_jupi[ 87]+pol_jupi[ 88]*nu)*cosa(zeta)*sina(ku)+
	(pol_jupi[ 90])*cosa(2*zeta)*sina(ku);
	ex1+=
	(pol_jupi[ 93]+pol_jupi[ 94]*nu)*sina(zeta)*cosa(ku)+
	(pol_jupi[ 96])*sina(2*zeta)*cosa(ku)+
	(pol_jupi[ 99])*cosa(ku)+
	(pol_jupi[102])*cosa(zeta)*cosa(ku)+
	(pol_jupi[105])*cosa(2*zeta)*cosa(ku)+
	(pol_jupi[108])*cosa(3*zeta)*cosa(ku);
	ex1+=
	(pol_jupi[111])*cosa(4*zeta)*cosa(ku)+
	(pol_jupi[114])*cosa(5*zeta)*cosa(ku)+
	(pol_jupi[117]+pol_jupi[118]*nu)*sina(zeta)*sina(2*ku)+
	(pol_jupi[120])*sina(2*zeta)*sina(2*ku)+
	(pol_jupi[123])*sina(3*zeta)*sina(2*ku)+
	(pol_jupi[126]+pol_jupi[127]*nu)*cosa(zeta)*sina(2*ku)+
	(pol_jupi[129])*cosa(2*zeta)*sina(2*ku)+
	(pol_jupi[132])*cosa(3*zeta)*sina(2*ku);
	ex1+=
	(pol_jupi[135]+pol_jupi[136]*nu)*sina(zeta)*cosa(2*ku)+
	(pol_jupi[138])*sina(2*zeta)*cosa(2*ku)+
	(pol_jupi[141])*sina(3*zeta)*cosa(2*ku)+
	(pol_jupi[144])*cosa(2*ku)+
	(pol_jupi[147]+pol_jupi[148]*nu)*cosa(zeta)*cosa(2*ku)+
	(pol_jupi[150])*cosa(2*zeta)*cosa(2*ku)+
	(pol_jupi[153])*cosa(3*zeta)*cosa(2*ku);
	ex1*=1e-7;

	b=
	(pol_jupi[156]+pol_jupi[157]*nu)*sina(ve)+
	(pol_jupi[159]+pol_jupi[160]*nu+pol_jupi[161]*nu2)*cosa(ve)+
	(pol_jupi[162]+pol_jupi[163]*nu)*sina(zeta)*sina(ku)+
	(pol_jupi[165])*sina(ku)+
	(pol_jupi[168])*cosa(zeta)*sina(ku)+
	(pol_jupi[171])*cosa(2*zeta)*sina(ku)+
	(pol_jupi[174])*cosa(3*zeta)*sina(ku)+
	(pol_jupi[177])*sina(zeta)*cosa(ku)+
	(pol_jupi[180])*sina(2*zeta)*cosa(ku);
	b+=
	(pol_jupi[183])*cosa(zeta)*cosa(ku)+
	(pol_jupi[186])*sina(zeta)*sina(2*ku)+
	(pol_jupi[189])*sina(2*zeta)*sina(2*ku)+
	(pol_jupi[192])*cosa(zeta)*sina(2*ku)+
	(pol_jupi[195])*cosa(2*zeta)*sina(2*ku)+
	(pol_jupi[198])*sina(zeta)*cosa(2*ku)+
	(pol_jupi[201])*sina(2*zeta)*cosa(2*ku)+
	(pol_jupi[204])*cosa(zeta)*cosa(2*ku)+
	(pol_jupi[207])*cosa(2*zeta)*cosa(2*ku);
	
	ax1=(
	(pol_jupi[210])*cosa(ve)+
	(pol_jupi[213])*cosa(zeta)+
	(pol_jupi[216])*cosa(2*zeta)+
	(pol_jupi[219])*cosa(3*zeta)+
	(pol_jupi[222])*cosa(4*zeta)+
	(pol_jupi[225])*sina(zeta)*sina(ku)+
	(pol_jupi[228])*cosa(2*zeta)*sina(ku)+
	(pol_jupi[231])*sina(2*zeta)*cosa(ku)+
	(pol_jupi[234])*sina(3*zeta)*cosa(ku)+
	(pol_jupi[237])*cosa(zeta)*cosa(ku)+
	(pol_jupi[240])*cosa(2*zeta)*cosa(ku)
	)*1e-6;

 (*ax)+=ax1;
 (*el)+=a;
 (*om)+=b;
 (*mm)+=a-b/(*ex);
 (*ex)+=ex1;
}
static	void sum_saturn_perturbations(	double *ax,double *ex,double *el,
					double *mm,double *om,double jc1)
{
 double	nu,nu2,pe,ku,es,ve,zeta,psi,a,b,ax1,ex1;

	nu  =jc1/5+0.1;
	nu2 =nu*nu;

	pe  =237.47555+3034.9061*jc1;
	ku  =265.91650+1222.1139*jc1;
	es  =243.51721+ 428.4677*jc1;

	ve  =5*ku-2*pe;
	zeta=  ku-  pe;
	psi =  es-  ku;

	a=
	(pol_satu[ 0]+pol_satu[ 1]*nu+pol_satu[ 2]*nu2)*sina(ve)+
	(pol_satu[ 3]+pol_satu[ 4]*nu+pol_satu[ 5]*nu2)*cosa(ve)+
	(pol_satu[ 6])*sina(2*ve)+
	(pol_satu[ 9])*sina(2*pe-6*ku+3*es)+
	(pol_satu[12])*sina(zeta)+
	(pol_satu[15])*sina(2*zeta)+
	(pol_satu[18])*sina(3*zeta)+
	(pol_satu[21])*sina(4*zeta)+
	(pol_satu[24])*sina(ku)+
	(pol_satu[27]+pol_satu[28]*nu)*sina(zeta)*sina(ku)+
	(pol_satu[30])*sina(2*zeta)*sina(ku)+
	(pol_satu[33])*sina(3*zeta)*sina(ku)+
	(pol_satu[36]+pol_satu[37]*nu)*cosa(zeta)*sina(ku);
	a+=
	(pol_satu[39])*cosa(2*zeta)*sina(ku)+
	(pol_satu[42]+pol_satu[43]*nu)*sina(zeta)*cosa(ku)+
	(pol_satu[45]+pol_satu[46]*nu)*cosa(zeta)*cosa(ku)+
	(pol_satu[48])*cosa(2*zeta)*cosa(ku)+
	(pol_satu[51])*cosa(3*zeta)*cosa(ku)+
	(pol_satu[54])*sina(zeta)*sina(2*ku)+
	(pol_satu[57])*sina(2*zeta)*sina(2*ku)+
	(pol_satu[60])*sina(3*psi)*sina(2*ku)+
	(pol_satu[63])*cosa(zeta)*cosa(2*ku)+
	(pol_satu[66])*cosa(2*zeta)*cosa(2*ku)+
       	(pol_satu[69])*cosa(3*psi)*cosa(2*ku);
	
	ex1=
	(pol_satu[ 72]+pol_satu[ 73]*nu+pol_satu[ 74]*nu2)*sina(ve)+
	(pol_satu[ 75]+pol_satu[ 76]*nu+pol_satu[ 77]*nu2)*cosa(ve)+
	(pol_satu[ 78]+pol_satu[ 79]*nu)*sina(2*ve)+
	(pol_satu[ 81]+pol_satu[ 82]*nu)*cosa(2*ve)+
	(pol_satu[ 84])*sina(2*zeta)+
	(pol_satu[ 87])*sina(ku)+
	(pol_satu[ 90]+pol_satu[ 91]*nu)*sina(zeta)*sina(ku)+
	(pol_satu[ 93]+pol_satu[ 94]*nu)*sina(2*zeta)*sina(ku)+
	(pol_satu[ 96])*cosa(zeta)*sina(ku);
	ex1+=
	(pol_satu[ 99])*cosa(2*zeta)*sina(ku)+
	(pol_satu[102])*cosa(3*zeta)*sina(ku)+
	(pol_satu[105])*cosa(4*zeta)*sina(ku)+
	(pol_satu[108])*cosa(5*zeta)*sina(ku)+
	(pol_satu[111])*cosa(2*psi)*sina(ku)+
	(pol_satu[114]+pol_satu[115]*nu)*cosa(ku);
	ex1+=
	(pol_satu[117])*sina(zeta)*cosa(ku)+
	(pol_satu[120])*sina(2*zeta)*cosa(ku)+
	(pol_satu[123])*sina(3*zeta)*cosa(ku)+
	(pol_satu[126])*sina(4*zeta)*cosa(ku)+
	(pol_satu[129])*sina(5*zeta)*cosa(ku)+
	(pol_satu[132]+pol_satu[133]*nu)*cosa(zeta)*cosa(ku)+
	(pol_satu[135]+pol_satu[136]*nu)*cosa(2*zeta)*cosa(ku)+
	(pol_satu[138])*sina(2*psi)*cosa(ku)+
	(pol_satu[141])*sina(2*ku);
	ex1+=
	(pol_satu[144]+pol_satu[145]*nu)*sina(zeta)*sina(2*ku)+
	(pol_satu[147])*sina(2*zeta)*sina(2*ku)+
	(pol_satu[150])*sina(3*zeta)*sina(2*ku)+
	(pol_satu[153])*sina(4*zeta)*sina(2*ku)+
	(pol_satu[156]+pol_satu[157]*nu)*cosa(zeta)*sina(2*ku)+
	(pol_satu[159]+pol_satu[160]*nu)*cosa(2*zeta)*sina(2*ku)+
	(pol_satu[162])*cosa(3*zeta)*sina(2*ku);
	ex1+=
	(pol_satu[165])*sina(3*psi)*sina(2*ku)+
	(pol_satu[168])*cosa(3*psi)*sina(2*ku)+
	(pol_satu[171])*cosa(2*ku)+
	(pol_satu[174]+pol_satu[175]*nu)*sina(zeta)*cosa(2*ku)+
	(pol_satu[177]+pol_satu[178]*nu)*sina(2*zeta)*cosa(2*ku)+
	(pol_satu[180])*sina(3*zeta)*cosa(2*ku);
	ex1+=
	(pol_satu[183]+pol_satu[184]*nu)*cosa(zeta)*cosa(2*ku)+
	(pol_satu[186]+pol_satu[187]*nu)*cosa(2*zeta)*cosa(2*ku)+
	(pol_satu[189])*cosa(3*zeta)*cosa(2*ku)+
	(pol_satu[192])*cosa(4*zeta)*cosa(2*ku)+
	(pol_satu[195])*sina(3*psi)*cosa(2*ku)+
	(pol_satu[198])*cosa(3*psi)*cosa(2*ku);
	ex1+=
	(pol_satu[201])*sina(zeta)*sina(3*ku)+
	(pol_satu[204])*sina(3*zeta)*sina(3*ku)+
	(pol_satu[207])*cosa(zeta)*cosa(3*ku)+
	(pol_satu[210])*cosa(3*zeta)*cosa(3*ku)+
	(pol_satu[213])*cosa(3*zeta)*sina(4*ku)+
	(pol_satu[216])*sina(3*zeta)*cosa(4*ku);
	ex1*=1e-7;

	b=
	(pol_satu[219]+pol_satu[220]*nu+pol_satu[221]*nu2)*sina(ve)+
	(pol_satu[222]+pol_satu[223]*nu+pol_satu[224]*nu2)*cosa(ve)+
	(pol_satu[225])*sina(zeta)+
	(pol_satu[228])*sina(zeta)*sina(ku)+
	(pol_satu[231])*sina(2*zeta)*sina(ku)+
	(pol_satu[234])*sina(3*zeta)*sina(ku);
	b+=
	(pol_satu[237])*cosa(ku)+
	(pol_satu[240])*cosa(zeta)*cosa(ku)+
	(pol_satu[243])*cosa(2*zeta)*cosa(ku)+
	(pol_satu[246])*cosa(3*zeta)*cosa(ku);
	b+=
	(pol_satu[249]+pol_satu[250]*nu)*sina(zeta)*sina(2*ku)+
	(pol_satu[252]+pol_satu[253]*nu)*cosa(zeta)*sina(2*ku)+
	(pol_satu[255]+pol_satu[256]*nu)*cosa(2*zeta)*sina(2*ku)+
	(pol_satu[258]+pol_satu[259]*nu)*sina(zeta)*cosa(2*ku)+
	(pol_satu[261])*sina(2*zeta)*cosa(2*ku)+
	(pol_satu[264]+pol_satu[265]*nu)*cosa(zeta)*cosa(2*ku)+
	(pol_satu[267]+pol_satu[268]*nu)*cosa(2*zeta)*cosa(2*ku);

	ax1=
	(pol_satu[270])*sina(ve)+
	(pol_satu[273])*cosa(ve)+
	(pol_satu[276])*cosa(zeta)+
	(pol_satu[279])*cosa(2*zeta)+
	(pol_satu[282])*cosa(3*zeta)+
	(pol_satu[285])*cosa(4*zeta)+
	(pol_satu[288])*cosa(5*zeta)+
	(pol_satu[291])*sina(ku)+
	(pol_satu[294])*sin(zeta)*sina(ku)+
	(pol_satu[297])*sina(2*zeta)*sina(ku)+
	(pol_satu[300])*sina(3*zeta)*sina(ku)+
	(pol_satu[303])*sina(4*zeta)*sina(ku);
	ax1+=
	(pol_satu[306])*cosa(zeta)*sina(ku)+
	(pol_satu[309])*cosa(2*zeta)*sina(ku)+
	(pol_satu[312])*cosa(3*zeta)*sina(ku)+
	(pol_satu[315])*cosa(4*zeta)*sina(ku)+
	(pol_satu[318])*cosa(ku)+
	(pol_satu[321])*sina(zeta)*cosa(ku)+
	(pol_satu[324])*sina(2*zeta)*cosa(ku);
	ax1+=
	(pol_satu[327])*sina(3*zeta)*cosa(ku)+
	(pol_satu[330])*sina(4*zeta)*cosa(ku)+
	(pol_satu[333])*cosa(zeta)*cosa(ku)+
	(pol_satu[336]+pol_satu[337]*nu)*cosa(2*zeta)*cosa(ku)+
	(pol_satu[339])*cosa(3*zeta)*cosa(ku)+
	(pol_satu[342])*sina(2*zeta)*sina(2*ku);
	ax1+=
	(pol_satu[345])*cosa(zeta)*sina(2*ku)+
	(pol_satu[348])*cosa(2*zeta)*sina(2*ku)+
	(pol_satu[351])*cosa(3*zeta)*sina(2*ku)+
	(pol_satu[354])*sina(zeta)*cosa(2*ku)+
	(pol_satu[357])*sina(2*zeta)*cosa(2*ku)+
	(pol_satu[360])*cosa(2*zeta)*cosa(2*ku);
	ax1+=
	(pol_satu[363])*cosa(3*zeta)*cosa(2*ku)+
	(pol_satu[366])*sina(zeta)*sina(3*ku)+
	(pol_satu[369])*cosa(3*zeta)*sina(3*ku)+
	(pol_satu[372])*cosa(zeta)*cosa(3*ku)+
	(pol_satu[375])*cosa(3*zeta)*cosa(3*ku);
	ax1*=1e-6;

 *ax+=ax1;
 *el+=a;
 *om+=b;
 *mm+=a-b/(*ex);
 *ex+=ex1;
}
int get_planet_coords(double jdate, int n, double *rx, double *ry, double *rz)
{
 double	jc1,jc2,jc3,jckk,x,y,z;
 double	el,ax,ex,in,om,oo;
 double	m0,m1,m2,m4,m5,m6,mm,vv;
 double	e,er,la,be;
 double	*pd;

 jc1 =(jdate-2415020.0)/36525.0;	  /* JCent, 1900.01.01 12:00 UT */
 jc2 =jc1*jc1;
 jc3 =jc2*jc1;
 jckk=(jdate-2451545.0)/36525.0;	  /* JCent, 2000.01.01 12:00 UT */

 switch(n)
  {	case 1: pd=(double*)merc;break;
	case 2: pd=(double*)venu;break;
	case 4: pd=(double*)mars;break;
	case 5: pd=(double*)jupi;break;
	case 6: pd=(double*)satu;break;
	case 7: pd=(double*)uran;break;
	case 8: pd=(double*)nept;break;
	case 9: pd=(double*)plut;break;
	default: pd=NULL;break;
  };

 if ( pd==NULL )	return(-1);

 if ( n<=6 )
  {	el=pd[ 0]+pd[ 6]*jc1+pd[12]*jc2+pd[18]*jc3;  /* L, mean longitude */
	ax=pd[ 1]+pd[ 7]*jc1+pd[13]*jc2+pd[19]*jc3;  /* a, semimajor axis */
	ex=pd[ 2]+pd[ 8]*jc1+pd[14]*jc2+pd[20]*jc3;  /* e, eccentricity */
	in=pd[ 3]+pd[ 9]*jc1+pd[15]*jc2+pd[21]*jc3;  /* i, inclination */
	om=pd[ 4]+pd[10]*jc1+pd[16]*jc2+pd[22]*jc3;  /* \omega, argument of perihelion */
	oo=pd[ 5]+pd[11]*jc1+pd[17]*jc2+pd[23]*jc3;  /* \Omega, longitude of ascending node */
  }
 else
  {	el=pd[ 0]+pd[ 6]*jckk;  /* L, mean longitude */
	ax=pd[ 1]+pd[ 7]*jckk;  /* a, semimajor axis */
	ex=pd[ 2]+pd[ 8]*jckk;  /* e, eccentricity */
	in=pd[ 3]+pd[ 9]*jckk;  /* i, inclination */
	om=pd[ 4]+pd[10]*jckk;  /* \varpi, longitude of perihelion */
	oo=pd[ 5]+pd[11]*jckk;  /* \Omega, longitude of ascending node */
	om=om-oo;
  }
    
/* Sun's mean anomaly: */
    
 m0=358.47583+ 35999.04975*jc1-0.000150*jc2-0.0000033*jc3;
    
/* Planets' mean anomaly: */
    
 m1=102.27938+149472.51529*jc1+0.000007*jc2;	/* merc */
 m2=212.60322+ 58517.80387*jc1+0.001286*jc2;	/* venu */
 m4=319.51913+ 19139.85475*jc1+0.000181*jc2;	/* mars */
 m5=225.32833+  3034.69202*jc1-0.000722*jc2;	/* jupi */
 m6=175.46622+  1221.55147*jc1-0.000502*jc2;	/* satu */
    
/* Periodic perturbations in the elements of the planet's orbit: */
    
 switch(n)
  {  case 1:
	mm=m1;
	break;
     case 2: 
	mm=m2;
	el+=	+0.00077*sina(237.24+150.27*jc1);
	break;
     case 4: 
	mm=m4;
	el+=	-0.01133*sina(3*m5-8*m4+4*m0)
		-0.00933*cosa(3*m5-8*m4+4*m0);
	break;
     case 5:
	mm=m5;

	sum_jupiter_perturbations(&ax,&ex,&el,&mm,&om,jc1);

	break;

     case 6: 
	mm=m6;

	sum_saturn_perturbations(&ax,&ex,&el,&mm,&om,jc1);


	break;
     case 7: case 8: case 9:
	mm=el-om-oo;	/* Default mean anomaly: M=L-\Omega-\omega=L-\varpi */
	break;
    };

 e =solve_kepler_equ(mm*M_PI/180.0,ex);
 er=ax*(1-ex*cos(e));
 vv=getpcoords(cos(e)-ex, sqrt(1-ex*ex)*sin(e));

 la=el+vv-mm-oo;

 x=er*cosa(la);
 y=er*sina(la);
 z=0;
 rotate(&y,&z,in);
 rotate(&x,&y,oo);
 be=asina(z/er);
 la=getpcoords(x,y);
 
 switch(n)
  {  case 1: 
	la+=	+0.00204*cosa(5*m2-2*m1+ 12.220)
		+0.00103*cosa(2*m2-  m1-160.692)
		+0.00091*cosa(2*m5-  m1- 37.003)
		+0.00078*cosa(5*m2-3*m1+ 10.137);

	er+=	+0.000007525*cosa(2*m5-  m1+ 53.013)
		+0.000006802*cosa(5*m2-3*m1-259.918)
		+0.000005457*cosa(2*m2-2*m1- 71.188)
		+0.000003569*cosa(5*m2-  m1- 77.75);
		
	break;
     case 2:
	la+=	+0.00313*cosa(2*m0-2*m2-148.225)
		+0.00198*cosa(3*m0-3*m2+  2.565)
		+0.00136*cosa(  m0-  m2-119.107)
		+0.00096*cosa(3*m0-2*m2-135.912)
		+0.00082*cosa(  m5-  m2-208.087);

	er+=	+0.000022501*cosa(2*m0-2*m2- 58.208)
		+0.000019045*cosa(3*m0-3*m2+ 92.577)
		+0.000006887*cosa(  m5-  m2-118.090)
		+0.000005172*cosa(  m0-  m2- 29.110)
		+0.000003620*cosa(5*m0-4*m2-104.208)
		+0.000003283*cosa(4*m0-4*m2+ 63.513)
		+0.000003074*cosa(2*m5-2*m2- 55.167);

	break;
     case 4:
	la+=	+0.00705*cosa(  m5-  m4- 48.958)
		+0.00607*cosa(2*m5-  m4-188.380)
		+0.00445*cosa(2*m5-2*m4-191.897)
		+0.00388*cosa(  m0-2*m4+ 20.495)
		+0.00238*cosa(  m0-  m4+ 35.097)
		+0.00204*cosa(2*m0-3*m4+158.638)
		+0.00177*cosa(3*m4-  m2- 57.602)
		+0.00136*cosa(2*m0-4*m4+154.093)
		+0.00104*cosa(  m5     + 17.618);

	er+=	+0.000053227*cosa(  m5-  m4+ 41.1306)
		+0.000050989*cosa(2*m5-2*m4-101.9847)
		+0.000038278*cosa(2*m5-  m4- 98.3292)
		+0.000015996*cosa(  m0-  m4- 55.5550)
		+0.000014764*cosa(2*m0-3*m4+ 68.6220)
		+0.000008966*cosa(  m5-2*m4+ 43.6150)
		+0.000007914*cosa(3*m5-2*m4-139.7370)
		+0.000007004*cosa(2*m5-3*m4-102.8880)
		+0.000006620*cosa(  m0-2*m4+113.2020)
		+0.000004930*cosa(3*m5-3*m4- 76.2430)
		+0.000004693*cosa(3*m0-5*m4+190.6030)
		+0.000004571*cosa(2*m0-4*m4+244.7020)
		+0.000004409*cosa(3*m5-  m4-115.8280);

	break;
     case 6:
	x=237.47555+3034.9061*jc1;	/* PE */
	y=265.91650+1222.1139*jc1;	/* KU */
	x=y-x;				/* ZETA = KU-PE */

	be+=	+0.000747*cosa(  x)*sina(  y)
		+0.001069*cosa(  x)*cosa(  y)
		+0.002108*sina(2*x)*sina(2*y)
		+0.001261*cosa(2*x)*sina(2*y)
		+0.001236*sina(2*x)*cosa(2*y)
		-0.002075*cosa(2*x)*cosa(2*y);

	break;

 };

 *rx=er*cosa(be)*cosa(la);
 *ry=er*cosa(be)*sina(la);
 *rz=er*sina(be);

 return(0);
}

/*****************************************************************************/
                                     
