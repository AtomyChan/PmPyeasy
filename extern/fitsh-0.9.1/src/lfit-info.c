/*****************************************************************************/
/* lfit-info.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Symbolic fitting & arithmetic evaluating utility. 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 1996, 2002, 2004-2005, 2006, 2007-2008, 2009; 			     */
/* Pal, A. (apal@szofi.net)						     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "longhelp.h"

#include <lfit/lfit.h>

#include "lfit-info.h"

/*****************************************************************************/

int fprint_lfit_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tlfit [-o <output>] [-h|--help|--version|--examples|--function-list]\n"
"\t[-V|--verbose] [--quiet]\n");
 fprintf(fw,
"Common parameters and arguments for fitting and/or regression analysis:\n"
"\t-v <fitvariable>[[:]=<initial>[:<error>][[<min>]:[<max>]]] [,...]\n"
"\t-g <derivedvariable>=<expression>[,...]\n"
"\t[-F|--format <fitvariable>=<numeric-printf-like-format>[,...]]\n"
"\t[-C|--correlation-format <numeric-printf-like-format>\n"
"\t[-q|--difference <fitvariable>=<difference>[,...]]\n"
"\t[-t|--constraint <linear constraint>=<value>[,...]]\n"
"\t[-x|--define|--macro <function>(v1,v2,...)=<definition> [-x ...]]\n");
 fprintf(fw,
"Parameters for simple fitting (single data block):\n"
"\t[<file>|-] -c <column names> -f <function>[=<dependent>]\n"
"\t[-y <dependent>] [-e <error>|-w <weight>]\n");
 fprintf(fw,
"Parameters for fitting more than one data blocks:\n"
"\t-i<key1> <file>|- -c<key1> <colum names 1> -f<key1> <funct1>[=<dep1>]\n"
"\t[-y<key1> <dependent1>] [-e<key1> <error1>|-w<key1> <weight1>]\n"
"\t-i<key2> <file>|- -c<key2> <colum names 2> -f<key2> <funct2>[=<dep2>]\n"
"\t[-y<key2> <dependent2>] [-e<key2> <error2>|-w<key2> <weight2>]\n"
"\t[...]\n");
 fprintf(fw,
"Fit methods:\tderivatives & linearity\n"
"\t[-L|--clls]\tyes\tyes\t"	"# Classic linear least squares\n"
"\t[-N|--nllm]\tyes\tno\t"	"# Nonlinear Levenberg-Marquardt\n"
"\t[-E|--emce]\topt.\topt.\t"	"# Refitting to synthetic data sets\n"
"\t[-M|--mcmc]\tno\tno\t"	"# Markov Chain Monte-Carlo\n"
"\t[-K|--mchi]\tno\tno\t"	"# Grid mapping of Chi^2\n"
"\t[-X|--xmmc]\tyes\tno\t"	"# Extended Markov Chain Monte-Carlo\n"
"\t[-U|--lmnd]\tno\tno\t"	"# Levenberg-Marquardt + num. derivatives\n"
"\t[-D|--dhsx]\topt.\tno\t"	"# Downhill simplex\n"
"\t[-A|--fima]\tyes\tno\t"	"# Fisher Information Matrix Analysis\n");
 fprintf(fw,
"Fine-tune parameters [-P|--params]:\n"
"\t-N -P\t[default],[lambda=<l>],[multiply=<m>],[iterations=<i>]\n"
"\t-U -P\t[default],[lambda=<l>],[multiply=<m>],[iterations=<i>] -q <...>\n"
"\t-D -P\t[default],[fisher]\n"
"\t-E -P\t[default],[skip],{linear|nonlinear|lmnd|dhsx[fisher]|mc}\n"
"\t-M -P\t[default],[[non]accepted],[gibbs]\n"
"\t-X -P\t[default],[[non]accepted],[skip],\n"
"\t\t[adaptive],[iterations=<N>],[window=<W>]\n"
"\t-A -P\t[default],[[no]orig],[[no]errors],[[no]corr],[mc|montecarlo]\n"
"\t[-U|-E [-q <...>]]\n");
 fprintf(fw,
"\t[-i|--{emce|mcmc|xmmc|fima}-iterations <n>] [-s|--seed <seed>|-1]\n"
"\t[--errors|--error-line|--error-columns] [--residual]\n"
"\t[--perturbations <key1>=<noise1>[,<key2>=<noise2>[,...]]]\n"
"\t[-r <rejection-level> -n <num-of-rejections> [--weighted-sigma]]\n"
"\t[-k|--separate [@],[!]<fitvariable>,...]\n");
#ifdef	LFIT_ENABLE_DYNAMIC_EXTENSIONS
 fprintf(fw,
"Dynamically loaded external libraries and functions:\n"
"\t[-d|--dynamic <library.so>:<array>[,<array2>,...] [-d <...>]]\n");
#endif
 fprintf(fw,
"General output specification:\n"
"\t[-o|--output <output-file-name]  [-z|--columns-output <>]\n");
 fprintf(fw,
"More outputs:\n"
"\t[-u <dump-used-to>] [-j <dump-rejected-to>]\n"
"\t[-a <dump-all-to> [--delta|--delta-comment]]\n"
"\t[-p <dump-expression-to>] [-l <dump-fitted-variables-to>]\n");
 return(0);
}

longhelp_entry lfit_long_help[]=
{
 LONGHELP_OPTIONS,

 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help",
	"Gives a detailed list of command line options." },
 { "--version",
 	"Gives some version information about the program." },
 { "--functions, --list-functions, --function-list", 
	"Lists the available arithmetic operations and built-in functions "
	"supported by the program." },
 { "--examples",
	"Prints some very basic examples for the program invocation." },

 { "Common options for regression analysis:", NULL },

 { "-v, --variable, --variables <list-of-variables>",
	"Comma-separated list of regression variables. In case of non-linear "
	"regression analysis, all of these fit variables are expected to have "
	"some initial values (specified as <name>=<value>), otherwise the "
	"initial values are set to be zero. Note that in the case of some of "
 	"the regression/analysis methods, additional parameters should be "
	"assigned to these fit/regression variables. See the section "
	"``Regression analysis methods'' for additional details." },

 { "-c, --column, --columns <independent>[:<column index>],...",
	__extension__ 
	"Comma-separated list of independet variable names as read from the "
	"subsequent columns of the primary input data file. If the independent "
	"variables are not in sequential order in the input file, the optional "
	"column indices should be defined for each variable, by separating the "
	"column index with a colon after the name of the variable. "
	"In the case of multiple input files and data blocks, the user "
	"should assign the individual independent variables and the "
	"respective column names and definitions for each file "
	"(see later, Sec. ``Multiple data blocks'')." },
 { "-f, --function <model function>",
	__extension__ 
	"Model function of the analysis in a symbolic form. This expression "
	"for the model function "
	"should contain built-in arithmetic operators, built-in functions, "
	"user-defined macros (see -x, --define) or functions provided by the "
	"dynamically loaded external modules (see -d, --dynamic). The model "
	"function can depend on both the fit/regression variables "
	"(see -v, --variables) and the independent variables read from the "
	"input file (see -c, --columns). In the case of multiple input files " 
	"and data blocks, the user should assign the respective model functions "
	"for each data block (see later). Note that some of the analysis methods "
	"expects the model function to be either differentiable or linear in "
	"the fit/regression variables. See ``Regression analysis methods'' later "
	"on about more details." },
 { "-y, --dependent <dependent expression>",
	__extension__
	"The dependent variable of the regression analysis, in a form of an "
	"arithmetic expression. This expression for the dependent variable "
	"can depend only on the variables read from the input file "
	"(see -c, --columns). In the case of multiple input files " 
	"and data blocks, the user should assign the respective dependent "
	"expressions for each data block (see later)." },
 { "-o, --output <output file>",
	"Name of the output file into which the fit results (the values for the "
	"fit/regression variables) are written." },

 { "Common options for function evaluation:", NULL },

 { "-f, --function <function to evaluate>[...]",
	"List of functions to be evaluated. More expressions can be specified by "
	"either separating the subsequent expressions by a comma or by "
	"specifying more -f, --function options in the command line." },

 { "Note that the two basic modes of `lfit` are distinguished only by the "
   "presence or the absence of the -y, --dependent command line argument. "
   "In other words, there isn't any explicit command line argument which "
   "specify the mode of `lfit`. If the -y, --dependent command line argument "
   "is omitted, `lfit` runs in function evaluation mode, otherwise "
   "the program runs in regression analysis mode.", NULL }, { "",NULL },
 
 { "-o, --output <output file>",
	"Name of the output file in which the results of the function "
	"evaluation are written." },

 { "Regression analysis methods:", NULL },

 { "-L, --clls, --linear",
	__extension__
	"The default mode of `lfit`, the classical linear least squares (CLLS) method. "
	"The model functions specified after -f, --function are expected to be "
	"both differentiable and linear with respect to the fit/regression "
	"variables. Otherwise, `lfit` detects the non-differentiable and "
	"non-linear property of the model function(s) and refuses the analysis. "
	"In this case, other types of regression analysis methods can be applied "
	"depending our needs, for instance the Levenberg-Marquardtalgorithm "
	"(NLLM, see -N, --nllm) or the downhill simplex minimization "
	"(DHSX, see -D, --dhsx)." },

 { "-N, --nllm, --nonlinear",
	__extension__
	"This option implies a regression involving the nonlinear "
	"Levenberg-Marquardt (NLLM) minimization "
	"algorithm. The model function(s) specified after -f, --function are "
	"expected to be differentiable with respect to the fit/regression "
	"variables. Otherwise, `lfit` detects the non-differentiable property "
	"and refuses the analysis. There some fine-tune parameters of the " 
	"Levenberg-Marquardt algorithm, see also the secion "
	"``Fine-tuning of regression analysis methods'' for more details how "
	"these additional regression parameters can be set. Note that all of "
	"the fit/regression variables should have a proper initial value, "
	"defined in the command line argument -v, --variable (see also there)." } ,

 { "-U, --lmnd",
	__extension__
	"Levenberg-Marquardt minimization with numerical partial "
	"derivatives (LMND). Same as the NLLM method, with the exception of "
	"that the partial "
	"derivatives of the model function(s) are calculated numerically. " 
	"Therefore, the model function(s) may contain functions of which "
	"partial derivatives are not known in an analytic form. "
	"The differences used in the computations of the partial derivatives "
	"should be declared by the user, see also the command line option "
	"-q, --differences. " },

 { "-D, --dhsx, --downhill",
	__extension__ 
	"This option implies a regression involving the nonlinear "
	"downhill simplex (DHSX) minimization algorithm. "
	"The user should specify the proper inital values "
	"and their uncertainties as <name>=<initial>:<uncertainty>, unless "
	"the ``fisher'' option is passed to the -P, --parameters command line "
	"argument (see later in the section ``Fine-tuning of regression "
	"analysis methods''). In the first case, the initial size of the "
	"simplex is based on the uncertainties provided by the user while "
	"in the second case, the initial simplex is derived from the "
	"eigenvalues and eigenvectors of the Fisher covariance matrix. "
	"Note that the model functions must be differentiable in the "
	"latter case. " },


 { "-M, --mcmc",
	__extension__
	"This option implies the method of Markov Chain Monte-Carlo (MCMC). "
	"The model function(s) can be arbitrary in the point of "
	"differentiability. However, each of the fit/regression variables "
	"must have an initial assumption for their uncertainties which must "
	"be specified via the command line argument -v, --variable. "
	"The user should specify the proper inital values "
	"and uncertainties of these as <name>=<initial>:<uncertainty>. "
	"In the actual implementation of `lfit`, each variable has an "
	"uncorrelated Gaussian a priori distribution with the specified "
	"uncertainty. The MCMC algorithm has some fine-tune parameters, " 
	"see the section ``Fine-tuning of regression analysis methods'' "
	"for more details." },

 { "-K, --mchi, --chi2",
	__extension__
	"With this option one can perform a ``brute force'' Chi^2 minimization "
	"by evaluating the value of the merit function of Chi^2 on a grid of "
	"the fit/regression variables. In this case the grid size and resolution "
	"must be specified in a specific form "
	"after the -v, --variable command line argument. Namely each of the "
	"fit/regression variables intended to be varied on a grid must "
	"have a format of <name>=[<min>:<step>:<max>] while the other ones "
	"specified as <name>=<value> are kept fixed. The output of this "
	"analysis will be a series of lines with N+1 columns, where the "
	"values of fit/regression variables are followed by the value "
	"of the merit function. Note that all of the declared fit/regression "
	"variables are written to the output, including the ones which are fixed "
	"(therefore the output is somewhat redundant)." },

 { "-E, --emce",
	__extension__
	"This option implies the method of ``refitting to synthetic data sets'', "
	"or ``error Monte-Carlo estimation'' (EMCE). This method must have a "
	"primarily assigned minimization algorithm (that can be any of the "
	"CLLS, NLLM or DHSX methods). First, the program searches the best fit "
	"values for the fit/regression variables involving the assigned primary "
	"minimization algorithm and reports these best fit variables. "
	"Then, additional synthetic data sets are generated around this set "
	"of best fit variables and the minimization is repeated involving the same "
	"primary method. The synthetic data sets are generated independently "
	"for each input data block, taking into account the fit residuals. "
	"The noise added to the best fit data is generated from the power "
	"spectrum of the residuals." },

 { "-X, --xmmc",
	__extension__
	"This option implies an improved/extended version of the Markov Chain "
	"Monte-Carlo analysis (XMMC). The major differences between the classic MCMC "
	"and XMMC methods are the following. 1/ The transition distribution "
	"is derived from the Fisher covariance matrix. 2/ The program performs "
	"an initial minimization of the merit function involving the method "
	"of downhill simplex. 3/ Various sanity checks are performed in order "
	"to verify the convergence of the Markov chains (including the "
	"comparison of the actual and theoretical transition probabilities, "
	"the computation of the autocorrelation lengths of each "
	"fit/regression variable series and the comparison of the statistical "
	"and Fisher covariance). " },

 { "-A, --fima",
	__extension__
	"Fisher information matrix analysis (FIMA). With this analysis method " 
	"one can estimate the uncertainties and correlations of the "
	"fit/regression variables involving the method of Fisher matrix "
	"analysis. This method does not minimize the merit functions by "
	"adjusting the fit/regression variables, instead, the initial values "
	"(specified after the -v, --variables option) are expected to be the "
	"``best fit'' ones." },

 { "Fine-tuning of regression analysis methods:", NULL },

 { "-e, --error <error expression>",
	"Expression for the uncertainties. Note that zero or negative "
	"uncertainty is equivalent to zero weight, i.e. input lines with zero "
	"or negative errors are discarded from the fit." },
 
 { "-w, --weight <weight expression>",
	"Expression for the weights. The weight is simply the reciprocal of " 
	"the uncertainty. The default error/uncertainty (and therefore "
	"the weight) is unity. Note that most of the analysis/regression "
	"methods are rather sensitive to the uncertainties since the merit "
	"function also depends on these." },
 
 { "-P, --parameters <regression parameters>",
	"This option is followed by a set of optional fine-tune parameters, "
	"that is different for each primary regression analysis method:" },

 { "default, defaults",
	"Use the default fine-tune parameters for the given regression method." },
 { "clls, linear",
	"Use the classic linear least squares method as the primary minimization "
	"algorithm of the EMCE method. Like in the case of the CLLS regression "
	"analysis (see -L, --clls), the model function(s) must be both differentiable and linear "
	"with respect to the fit/regression variables. " },
 { "nllm, nonlinear",
	"Use the non-linear Levenberg-Marquardt minimization algorithm as "
	"the primary minimization algorithm of the EMCE method. Like in the case of the NLLM regression "
	"analysis (see -N, --nllm), the model function(s) must be differentiable "
	"with respect to the fit/regression variables. " },
 { "lmnd",
	"Use the non-linear Levenberg-Marquardt minimization algorithm as "
	"the primary minimization algorithm of the EMCE method. "
	"Like in the case of -U, --lmnd regression method, the parametric "
	"derivatives of the model function(s) are calculated by a "
	"numerical approximation (see also -U, --lmnd and -q, --differences for "
	"additional details)." },
 { "dhsx, downhill",
	"Use the downhill simplex (DHSX) minimization as the primary "
	"minimization algorithm of the EMCE method. Unless the additional "
	"'fisher' option is specified directly, like in the default case of "
	"the DHSX regression method, the user should specify the uncertainties "
	"of the fit/regression variables that are used as an initial size "
	"of the simplex." },

 { "mc, montecarlo",
	"Use a primitive Monte-Carlo diffusion minimization technique as the "
	"primary minimization algorithm of the EMCE method. The user should "
	"specify the uncertainties of the fit/regression variables which are "
	"then used to generate the Monte-Carlo transitions. This primary "
	"minimization technique is rather nasty (very slow), "
	"so its usage is not recommended. " },

 { "fisher",
	__extension__
	"In the case of the DHSX regression method or in the case of the EMCE "
	"method when the primary minimization is the downhill simplex algorithm, "
	"the initial size of the simplex is derived from the Fisher covariance "
	"approximation evaluated at the point represented by the initial "
	"values of the fit/regression variables. Since the derivation of the "
	"Fisher covariance requires the knowledge of the partial derivatives "
	"of the model function(s) with respect to the fit/regression variables, "
	"the(se) model function(s) must be differentiable. On the other hand, "
	"the user do not have to specify the initial uncertainties after the "
	"-v, --variables option since these uncertainties derived automatically "
	"from the Fisher covariance." },

 { "skip",
	"In the case of EMCE and XMMC method, the initial minimization "
	"is skipped. " },

 { "lambda=<value>",
	"Initial value for the ``lambda'' parameter of the Levenberg-Marquardt "
	"algorithm. " },
 { "multiply=<value>",
	"Value of the ``lambda multiplicator'' parameter of the "
	"Levenberg-Marquardt algorithm. " },
 { "iterations=<max.iterations>",
	"Number of iterations during the Levenberg-Marquardt algorithm. " },

 { "accepted",
	"Count the accepted transitions in the MCMC and XMMC methods (default)." },
 { "nonaccepted",
	"Count the total (accepted plus non-accepted) transitions in the MCMC "
	"and XMMC methods." },
 { "gibbs",
	"Use the Gibbs sampler in the MCMC method. " },
 { "adaptive",
	"Use the adaptive XMMC algorithm (i.e. the Fisher covariance is "
	"re-computed after each accepted transition). " },
 { "window=<window size>",
	"Window size for calculating the autocorrelation lengths for the "
	"Markov chains (these autocorrelation lengths are reported only "
	"in the case of XMMC method). The default value is 20, which is "
	"fine in the most cases since the typical autocorrelation lengths "
	"are between 1 and 2 for nice convergent chains." },

 { "-q, --difference <variablename>=<difference>[,...]",
	"The analysis method of LMND (Levenberg-Marquardt minimization using "
	"numerical derivatives, see -U, --lmnd) requires the differences that "
	"are used during the computations of the partial derivatives of the "
	"model function(s). With this option, one can specify these differences." },

 { "-k, --separate <variablename>[,...]",
	__extension__
	"In the case of non-linear regression methods (for instance, DHSX or XMMC) "
	"the fit/regression variables in which the model functions are linear "
	"can be separated from the nonlinear part and therefore make the "
	"minimization process more robust and reliable. Since the "
	"set of variables in which the model functions are linear is "
	"ambiguous, the user should explicitly specify this supposedly "
	"linear subset of regression variables. "
	"(For instance, the model function ``a*b*x+a*cos(x)+b*sin(x)+c*x^2'' "
	"is linear in both ``(a,c)'' and ``(b,c)'' parameter vectors but it "
	"is non-linear in ``(a,b,c)''.) " 
	"The program checks whether the specified subset of "
	"regression variables is a linear subset and reports "
	"a warning if not. " 
	"Note that the subset of separated linear variables (defined here) "
	"and the subset of the fit/regression variables affected by "
	"linear constraints (see also section ``Constraints'') "
	"must be disjoint." },

 { "--perturbations <noise level>, --perturbations <key>=<noise level>[,...]",
	"Additional white noise to be added to each EMCE synthetic data sets. " 
	"Each data block (referred here by the approprate data block keys, "
	"see also section ``Multiple data blocks'') may have different white "
	"noise levels. If there is only one data block, this command line "
	"argument is followed only by a single number specifying the "
	"white noise level." },

 { "Additional parameters for Monte-Carlo analysis:", NULL },
 { "-s, --seed <random seed>",
	__extension__
	"Seed for the random number generator. By default this seed is 0, thus "
	"all of the Monte-Carlo regression analyses (EMCE, MCMC, XMMC "
	"and the optional generator for the FIMA method) generate "
	"reproducible parameter distributions. A positive value after this "
	"option yields alternative random seeds while all negative values "
	"result in an automatic random seed (derived from various available "
	"sources, such as /dev/[u]random, system time, hardware MAC address "
	"and so), therefore distributions generated involving this kind of "
	"automatic random seed are not reproducible." },
 { "-i, --[mcmc,emce,xmmc,fima]-iterations <iterations>",
	"The actual number of Monte-Carlo iterations for the MCMC, EMCE, "
	"XMMC methods. Additionally, the FIMA method is capable to generate " 
	"a mock Gaussian distribution of the parameter with the same covariance "
	"as derived by the Fisher analysis. The number of points in this mock "
	"distribution is also specified by this command line option. "},

 { "Clipping outlier data points:", NULL },

 { "-r, --sigma, --rejection-level <level>",
	"Rejection level in the units of standard deviations." },
 { "-n, --iterations <number of iterations>",
	"Maximum number of iterations in the outlier clipping cycles. "
	"The actual number of outlier points can be traced by increasing the "
	"verbosity of the program (see -V, --verbose)." },
 { "--[no-]weighted-sigma",
	"During the derivation of the standard deviation, the contribution of "
	"the data points data points can be weighted by the respective "
	"weights/error bars (see also -w, --weight or -e, --error in the "
	"section ``Fine-tuning of regression analysis methods''). If no "
	"weights/error bars are associated to the data points (i.e. both "
	"-w, --weight or -e, --error options are omitted), this option will have "
	"no practical effect." },

 { "Note that in the actual version of `lfit`, only the CLLS, NLLM and LMND "
   "regression methods support the above discussed way of "
   "outlier clipping.", NULL }, { "", NULL },

 { "Multiple data blocks:", NULL },

 { "-i<key> <input file name>",
	"Input file name for the data block named as <key>." },
 { "-c<key> <independent>[:<column index>],...",
	"Column definitions (see also -c, --columns) for the given data " 
	"block named as <key>." },
 { "-f<key> <model function>",
	"Expression for the model function assigned to the data block named as <key>." },
 { "-y<key> <dependent expression>",
	"Expression of the dependent variable for the data block named as <key>." },
 { "-e<key> <errors>",
	"Expression of the uncertainties for the data block named as <key>." },
 { "-w<key> <weights>",
	"Expression of the weights for the data block named as <key>. Note that "
	"like in the case of -e, --errors and -w, --weights, only one of the "
	"-e<key>, -w<key> arguments should be specified." },
	
 { "Constraints:", NULL },

 { "-t, --constraint, --constraints <expression>{=<>}<expression>[,...]",
	__extension__
	"List of fit and domain constraints between the regression variables. "
	"Each fit constraint expression must be linear in the fit/regression variables. "
	"The program checks the linearity of the fit constraints and reports an "
	"error if any of the constraints are non-linear. "
	" A domain constraint can be any " 
	"expression involving arbitrary binary arithmetic relation (such as "
	"strict greater than: '>', strict less than: '<', "
	"greater or equal to: '>=' and less or requal to: '<='). "
	"Constraints can be "
	"specified either by a comma-separated list after a single command "
	"line argument of -t, --constraints or by multiple of these "
	"command line arguments. "
 },
 { "-v, --variable <name>:=<value>",
	"Another form of specifying constraints. The variable specifications "
	"after -v, --variable can also be used to define constraints by writing "
	"``:='' instead of ``='' between the variable name and initial value. "
	"Thus, -v <name>:=<value> is equivalent to -v <name>=<value> "
	"-t <name>=<value>." },

 { "User-defined functions:", NULL },

 { "-x, --define, --macro <name>(<parameters>)=<definition expression>",
	__extension__
	"With this option, the user can define additional functions "
	"(also called macros) on the top of the built-in functions and operators, "
	"dynamically loadaded functions and previously defined macros. "
	"Note that each such user-defined function must be stand-alone, i.e. " 
	"external variables (such as fit/regression variables and independent "
	"variables) cannot be part of the definition expression, only the "
	"parameters of these functions." },

 { "Dynamically loaded extensions and functions:", NULL },

 { "-d, --dynamic <library>:<array>[,...]",
	__extension__
	"Load the dynamically linked library (shared object) named <library> " 
	"and import the global `lfit`-compatible set of functions defined "
	"in the arrays specified after the name of the library. The arrays "
	"must have to be declared with the type of 'lfitfunction', as it is "
	"defined in the file ``lfit.h''. Each record in this array contains "
	"information about a certain imported function, namely "
	"the actual name of this function, flags specifying whether the "
	"function is differentiable and/or linear in its regression parameters, "
	"the number of regression variables and independent variables "
	"and the actual C subroutine that implements the evaulation of the "
	"function (and the optional computation of the partial derivatives). " 
	"The module 'linear.c' and 'linear.so' provides a simple example "
	"that implements the ``line(a,b,x)=a*x+b'' function. "
	"This example function has "
	"two regression variables (``a'' and ``b'') and one independent "
	"variable (``x'') and the function itself is linear in the regression "
	"variables." },

 { "More on outputs:", NULL },

 { "-z, --columns-output <column indices>",
	__extension__
	"Column indices where the results are written in evaluation mode. " 
	"If this option is omitted, the results of the function evaluation "
	"are written sequentally. Otherwise, the input file is written to "
	"the output and the appropriate columns (specified here) are replaced " 
	"by the respective results of the function evaluation. Thus, although "
	"the default column order is sequential, there is a significant "
	"difference between omitting this option and specifying ``-z 1,2,...,N''. "
	"In the first case, the output file contains only the results of the "
	"function evaluations, while in the latter case, the first N columns "
	"of the original file are replaced with the results. " },

 { "--errors, --error-line, --error-columns",
	"Print the uncertainties of the fit/regression variables." },

 { "-F, --format <variable name>=<format>[,...]",
	"Format of the output in printf-style for each fit/regression variable"
	"(see printf(3)). The default "
	"format is %12.6g (6 signifiant figures)." },

 { "-F, --format <format>[,...]",
	"Format of the output in evaluation mode. The default "
	"format is %12.6g (6 signifiant figures)." },

 { "-C, --correlation-format <format>",
	"Format of the correlation matrix elements. The default format "
	"is %6.3f (3 significant figures)." },

 { "-g, --derived-variable[s] <variable name>=<expression>[,...]", 
	"Some of the regression and analysis methods are capable to " 
	"compute the uncertainties and correlations for derived regression "
	"variables. These additional (and therefore not independent) "
	"variables can be defined with this command line option. "
	"In the definition expression one should use only the fit/regression "
	"variables (as defined by the -v, --variables command line argument). "
	"The output format of these variables can also be specified by the " 
	"-F, --format command line argument." },

 { "-u, --output-fitted <filename>",
	"Neme of an output file into which those lines of the input are "
	"written that were involved in the final regression. This option "
	"is useful in the case of outlier clipping in order to see what "
	"was the actual subset of input data that was used in the fit "
	"(see also the -n, --iterations and -r, --sigma options)." },

 { "-j, --output-rejected <filename>",
	"Neme of an output file into which those lines of the input are "
	"written that were rejected from the final regression. This option "
	"is useful in the case of outlier clipping in order to see what "
	"was the actual subset of input data where the dependent variable "
	"represented outlier points "
	"(see also the -n, --iterations and -r, --sigma options)." },

 { "-a, --output-all <filename>",
	"File containing the lines of the input file that were involved "
	"in the complete regression analysis. This file is simply the "	
	"original file, only the commented and empty lines are omitted. " },

 { "-p, --output-expression <filename>",
 	"In this file the model function is written in which the "
	"fit/regression variables are replaced by their best-fit values. " },

 { "-l, --output-variables <filename>",
	"List of the names and values of the fit/regression variables in the "
	"same format as used after the -v, --variables command line argument. "
	"The content of this file can therefore be passed to subsequent "
	"invocations of `lfit`." },

 { "--delta",
	"Write the individual differences between the independent "
	"variables and the evaluated best fit model function values for each "
	"line in the output files specified by the -u, --output-fitted, "
	"-j, --output-rejected and -a, --output-all command line options." },

 { "--delta-comment",
	"Same as --delta, but the differences are written as a comment " 
	"(i.e. separated by a '##' from the original input lines)." },

 { "--residual",
	"Write the final fit residual to the output file (after the list of " 
	"the best-fit values for the fit/regression variables)." },

 { NULL, NULL }
};

int fprint_lfit_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tlfit [method of analysis] [options] <input> [-o, --output <output>]\n"
"The program `lfit` is a standalone command line driven tool designed for\n"
"both interactive and batch processed data analysis and regression.\n");
 fprintf(fw,
"In principle, the program may run in two modes. First, `lfit` supports \n"
"numerous regression analysis methods that can be used to search for \n"
"``best fit'' parameters of model functions in order to model the input\n"
"data (which are read from one or more input files in tabulated form).\n");
 fprintf(fw,
"Second, `lfit` is capable to read input data and performs various \n"
"arithmetic operations as it is specified by the user. Basically this second\n"
"mode is used to evaluate the model functions with the parameters presumably\n"
"derived by the actual regression methods (and in order to complete this \n"
"evaluation, only slight changes are needed in the command line invocation \n"
"arguments).\n");
 fprintf(fw,"\n");

 longhelp_fprint(fw,lfit_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <apal@szofi.net>\n");

 return(0);
}

int fprint_lfit_examples(FILE *fw)
{
 fprintf(fw,
"Examples:\n"
"* linear regression in 1+1 dim:\n"
"\tlfit -v a,b -c x,y -f \"a*x+b\" -y y\n"
"* linear regression in 1+1 dim, taking into account errors (from 3d column):\n"
"\tlfit -v a,b -c x,y,yerr -f \"a*x+b\" -y y -e yerr\n"
"* fit a circle in 2 dim (find the center and radius of a circle which \n"
"  fits well to the points given in the first two columns of the input):\n"
"\tlfit -v u,v,uvr -c x,y -f \"2*x*u+2*y*v-uvr\" -y \"x*x+y*y\"\n"
"  The center of the circle is (u,v) and the raduis is r=sqrt(u*u+v*v-uvr).\n");
 fprintf(fw,
"Some notes: \n"
" - Avoid to put spaces or any nasty characters (such as wildchards: *, ?,\n"
"   backslash, brackets -- which are pre-interpreted by the shell) directly\n"
"   in the arguments of -v, -c, -f, -y and -e|-w. Simply put such arguments\n"
"   between quotation marks (or escape them), see the examples above.\n"
" - If the function 'function' is linear in the variables 'vars' are to be\n"
"   fitted, the standard linear regression algorithm will be used. Otherwise,\n");
 fprintf(fw,
"   the Levenberg-Marquardt method will be used (with the parameters optionally\n"
"   defined by the switches -l, -m and -i), if the switch -N is specified in the\n"
"   command line to force the non-linear method. In this case, the initial\n"
"   values of the variables 'vars' are important to be defined in the argument\n"
"   of -v. With the switch -V (--verbose) the evolution of the variables and\n"
"   the lambda parameter can be traced ('max_iter' lines are written to stderr).\n");
 return(0);
}

/*****************************************************************************/
                                                         
