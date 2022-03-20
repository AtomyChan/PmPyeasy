/*****************************************************************************/
/* psn.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program is free software; you can redistribute it and/or modify	     */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 2 of the License, or	     */
/* (at your option) any later version.					     */
/*                                                                           */
/* This library is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111, USA.    */
/*****************************************************************************/

#ifndef __PSN_H_INCLUDED
#define __PSN_H_INCLUDED

#define		PSN_MAX_PREC		32766

/*********** Definition of objects used by the psn_...() routines ************/

/* A term of a PSN sequence: */
typedef struct
 {	short	type;
	short	major;
	short	minor;
	short	cache;
 } psnterm;

/* Full information about a PSN sequence: */
typedef struct
{	psnterm	*terms;
	double	*cons;
	int	nterm;
	int	ncon;
	int	ntermalloc;
	int	nconalloc;
	int	nseq;
} psn;

/* Symbol: variable, operator or function: */
typedef struct
 {	int	type;
	int	major;
	char	*name;
	int	minor;
 } psnsym;

/* Properties of an operator or a function (major code, number of input	     */
/* parameters, percedency and associativity. For operators, TO_PREFIX,	     */
/* TO_INFIX or TO_SUFFIX should be used in the 'innum' field.		     */
typedef struct
 {	int	major;
	int	argnum;
	int	precedency;
	int	associativity;
 } psnprop;

/* Derivative rules: */
typedef struct
 {	int	major;
	short	*oplist;
 } psndiff;

/* Function descriptor, called by psn_double_calc(): */
typedef struct 
 {	int	major;
	int	(*funct)(double *top_of_the_stack_pointer);
 } psnfunct;

/* Rules for the description of the symbolic evaluation (thus, conversion    */
/* of the PSN sequences to any kind of character string), used by	     */
/* psn_convert_symbolic():						     */
typedef struct
 {	int	major;
	char	*string;
	int	strength;
	int	affixation;
 } psnsymeval;
	
/************************ Definition of global constants *********************/

#define	T_END		0x00	/* end-of-sequence			     */

#define	T_PAROPEN	0x10	/* parenthesis '('			     */
#define	T_SEP		0x11	/* separator 				     */
#define	T_PARCLOSE	0x12	/* parenthesis ')'			     */

#define	T_CONST		0x20	/* numerical constants			     */
#define	T_SCONST	0x21	/* special constants (=major)		     */
#define	T_VAR		0x22	/* variables				     */
#define	T_STACKVAR	0x23	/* variables on stack			     */

#define	T_OP		0x30	/* operators				     */
#define	T_FN		0x31	/* functions				     */
				/* in PSN sequences these are equivalent to  */
				/* T_OPs... the differece is between	     */
				/* associativity (functions are left-assoc)  */
				/* and precedence (maximal for functions)    */

#define	TO_PREFIX	(-1)	/* prefix operators (e.g. negation)	     */
#define	TO_INFIX	( 0)	/* infix  operators (e.g. addition,...)	     */
#define	TO_SUFFIX	( 1)	/* suffix operators (e.g. factorial, 'prime')*/

/********* Precedence and associativity definitions, psnprec->assoc **********/

#define	ASSOC_LEFT	( 0)	/* Left-associated operators (generic)	     */
#define	ASSOC_RIGHT	( 1)	/* Right-associated operators (e.g. power:^) */

/***** PSN stack indices, used in the description of derivative rules... *****/

#define	SS_REF		128	/* maximum number of stack elements, may be  */
				/* changed if functions with more than 9     */
				/* arguments are expected to be used...	     */
				/* If changed, the whole library should be   */
				/* re-compiled.				     */

#define	SS_D1	( -1)		/* refers to the top of the derivative stack */
#define	SS_D2	( -2)
#define	SS_D3	( -3)
#define	SS_D4	( -4)
#define	SS_D5	( -5)
#define	SS_D6	( -6)
#define	SS_D7	( -7)
#define	SS_D8	( -8)

#define	SS_N1	( -1-SS_REF)	/* refers to the top of the stack	     */
#define	SS_N2	(- 2-SS_REF)
#define	SS_N3	(- 3-SS_REF)
#define	SS_N4	(- 4-SS_REF)
#define	SS_N5	(- 5-SS_REF)
#define	SS_N6	(- 6-SS_REF)
#define	SS_N7	(- 7-SS_REF)
#define	SS_N8	(- 8-SS_REF)

#define	CON_0	(-0-2*SS_REF)	/* refers to 0, [T_SCONT, major=0]	     */
#define	CON_1	(-1-2*SS_REF)	/* refers to 1, [T_SCONT, major=1]	     */
#define	CON_2	(-2-2*SS_REF)	/* refers to 2, [T_SCONT, major=2]	     */
#define	CON_3	(-3-2*SS_REF)	/* refers to 3, [T_SCONT, major=3]	     */

/*************************** Simplification rules ****************************/

#define	S_EXP1	(-1)
#define	S_EXP2	(-2)
#define	S_EXP3	(-3)
#define	S_EXP4	(-4)
#define	S_EXP5	(-5)
#define	S_EXP6	(-6)
#define	S_EXP7	(-7)
#define	S_EXP8	(-8)

/***************************** psnerrno values *******************************/

extern	int	psnerrno;

#define	PEOK		0	/* Everything is ok...			     */
#define	PEALLOC		-1	/* Allocation error		  	     */
#define	PEINVALID	1	/* Unexpected character or term in a sequence*/
#define	PENOTFOUND	2	/* Required function/symbol not found        */
#define	PEPARSE		3	/* Parse error				     */
#define	PEANALYTICAL	4	/* failure due to analytical error	     */
#define	PENUMERICAL	5	/* failure due to numerical error	     */
#define	PENONDIFF	6	/* Non-differentable operator in a sequence  */
#define	PEOPTIMIZED	7	/* Sequence is optimized, but it must not be */
#define	PEMULTI		8	/* Sequence contains more sub-sequences      */

/* Note that the exact meanings of these pre-defined errors are depend on    */
/* the function which returned these values or set 'psnerrno'. See the 	     */
/* comments below and the full documentation in ./doc for more details.	     */

/*************************** Function prototypes *****************************/

/* psn_conv_string():
   Encodes the string 'in' using the symbols in 'symtable'. The symbol table
   'symtable' should contain a NULL-terminated list of pointers, these
   pointers point to {0,0,NULL}-terminated lists of symbols. If all operators 
   and alphabetic symbols can be resolved using the symbol table 'symtable', 
   the function returns a pointer to a newly allocated PSN structure, creates 
   the list 'sequence' and initializes the constant table 'cons'. If any of the
   operators or alphabetic symbols cannot be resolved, the encoding fails and 
   the function returns NULL. If the function returns NULL, the type of the
   error can be figured out by using the external variable 'psnerrno'.	     */
psn *	psn_conv_string(char *in,psnsym **symtable);

/* psn_init():
   Initializes the sequence 'seq'. If the initialization is successful, the
   function returns 0, otherwise the function returns a non-zero value and
   the type of the error can be found out from comparing it to the PE*
   pre-defined values.							     */
int	psn_init(psn *seq,psnprop *prop);

/* psn_cache_init(), psn_cache_clear():
   Initializes and clears the lookup-table (called cache) in the sequence 
   'seq'. This lookup table is used by psn_double_calc(). If the lookup table 
   was not initialized before the call of psn_double_calc() is not a problem,
   the first call of psn_double_calc() will do the initialization. Note that 
   in such (almost rare) cases, when the function list 'flist' changes between 
   two calls of psn_double_calc(), the cache must be cleared by 
   psn_cache_clear() or re-initialized by psn_cache_init() before the next 
   call of psn_double_calc(). The cache does not have to be cleared or 
   initialized before calling other funcions than psn_double_calc() which
   have any argument with the type of (psnfunct*), but some of the functions
   (e.g. psn_diff(), psn_simplify()) does clear the cache.		     */
int	psn_cache_init(psn *seq,psnfunct *flist);
int	psn_cache_clear(psn *seq);
int	psn_cache_search(int major,psnfunct *flist);

/* psn_conv():
   Re-orders the terms in the sequence 'in' to be a really suffix-notated
   list of operands and operators. The precendencies of the operators are
   stored in the {0,0,0}-terminated list 'plist'. If the conversion fails
   (e.g. due to syntax or parse error), the functions returns NULL, otherwise
   returns a pointer newly created PSN structure, so the input sequence 'in'
   is not destroyed. The new sequence inherits the constant table 'cons' of 
   the sequence 'in'. If the function returns NULL, the type of the
   error can be figured out by using the external variable 'psnerrno'.	     */
psn *	psn_conv(psn *in,psnprop *plist);

/* psn_append_term(), psn_append(), psn_cons_append():
   Append one or more psnterms or a constant ({type,major,minor,flag}, or the 
   'nterm'-lenght list 'terms') to the sequence 'seq'.			     */
void	psn_append_term(psn *seq,int type,int major,int minor,int flag);
void	psn_append(psn *seq,psnterm *terms,int nterm);
int	psn_cons_append(psn *seq,double con);

/* psn_fprint():
   Writes the sequence 'seq' to the stream 'f', using the symbol table 
   'symtable' to resolve the names of symbols (useful for debugging).	     */
void	psn_fprint(FILE *f,psn *seq,psnsym **symtable);

/* psn_get_length():
   Gets the number of sequence terms of the sequence 'seq'.		     */
int	psn_get_length(psn *seq);
/* psn_get_nseq():
   Gets how many individual subsequences are stored in the sequence 'seq'.   */
int	psn_get_nseq(psn *seq);

/* psn_get_nstack():
   Gets how many terms remain on the stack after the evaluation of the
   PSN sequence 'pseq'. If the sequence 'pseq' is incosistent, the function
   returns a negative value. If the optional argument 'ropt' is not NULL,
   this function stores the number of the optimized elements remaining on the 
   stack (note that this value is always less or equal to the value 
   returned).								     */
int	psn_get_nstack(psn *pseq,int *ropt);

/* psn_test():
   Tests the integrity of the PSN sequence 'seq': tries to evaluate the 
   expression coded in the sequence 'seq'. If the test fails, the function 
   returns a non-zero value (see PE* values), otherwise -- if the sequence is 
   consistent -- it returs 0. Note that before the call of psn_double_calc(), 
   all sequences should be checked by such a function like psn_test(). To make 
   the evaluation faster, psn_double_calc() does not do any checking and may 
   spectacularly fail (e.g. terminates with segmentation fault, SIGSEGV) if 
   the input sequence is incosistent.					     */
int	psn_test(psn *seq);

/* psn_double_calc():
   Evaluates the expression coded in the sequence 'seq' using the functions and
   operands listed in 'functlist'. In this library, this is the only function 
   related to the evaluation of PSN-based sequences and all operands are 
   treated as operands with the type of double. If sequences based on 
   non-real operands (or any other kind of operands which are not having a 
   type of double, e.g. complex variables, vectors, matrices) do want to be 
   used, the appropriate functions like this should be written.
   After the evaulation, the results are stored in 'result': this array  
   contains `psn_get_nseq(seq)` values. The array 'vars' should contain the 
   appropriate values of the variables (T_VARs). Note that the sequence
   'seq' should be checked (e.g. with psn_test()) before calling this 
   function, see also the comments near psn_test(). The functon returns
   0, if the caclulation was successful, otherwise it returns a non-zero
   value (see also the pre-defined error id's: PE*)			     */
int	psn_double_calc(psn *seq,psnfunct *functlist,double *result,
	double *vars);

/* psn_diff():
   Caluclates the partial derivative of 'seq' by the the variable 'var_index'
   using the derivative rules listed in the {0,0,NULL}-terminated list 
   'diffrules'. If the sequence 'seq' contains a non-differentable operator,
   or in any other cases, when the requested partial derivative cannot be
   calculated, the function returns NULL (and the type of the error is 
   indicated by 'psnerrno'), otherwise it returns a pointer which points to 
   a newly created sequence with the desired partial derivative.	     */
psn *	psn_diff(psn *seq,psnprop *plist,psndiff *diffrules,int var_index);

/* psn_simplify():
   Simplifies the sequence 'seq' using the rules listed in 'rules' and the 
   operator properties listed in 'plist'. If the simplification was successful,
   the function returns 0, otherwise it returns a non-zero value (see also
   the PE* values). Operators and functions with pure numerical (T_CONST, 
   T_SCONST) arguments can be evaluated (and replaced by the result) also, if 
   the list of the available functions is passed in the argument 'flist'. If 
   no such simplification is desired, pass NULL in the argument 'flist' to 
   psn_simplift().							     */
int	psn_simplify(psn *seq,short **rules,psnprop *plist,psnfunct *flist);

/* psn_diff_simplify():
   This function is subsequent combination of psn_diff() and psn_simplify().
   The call also checks the integrity of the resulted PSN sequence using
   psn_test(), i.e. if the return value is not NULL (a valid PSN sequence),
   there is no need to call psn_test() after this call. Note that if 'slist'
   is NULL, the function does not do any simplifications, in this case
   it is equivalent with psn_diff() and a psn_test() on its result.          */
psn *psn_diff_simplify(psn *seq,psnprop *plist,psndiff *dlist,
        short **slist,psnfunct *flist,int var_index);

/* psn_optimize():
   Tries to optimize the evaluation of the sequence 'seq'. If the sequence 
   'seq' has been optimized yet, or the optimization fails, the function 
   returns a non-zero value (see also the PE* values). If the optimization 
   was successful, the function returns 0.				     */
int	psn_optimize(psn *seq);

/* psn_replace():
   Replaces the operator/function in the sequence 'seq' with the major code 
   'major' with the sequence 'replacement'. The arguments of the given
   operator/functions are replaced by the variables (T_VARs) with the 
   codes specified in the array 'varlist'. If 'varlist' is NULL, then the 
   first argument is replaced by the variable with the code 0, the second
   is replaced by the variable with the code 1 and so on. The function
   returns a new (replaced) sequence or NULL if something is wrong.	     */
psn *	psn_replace(psn *seq,int major,psn *replacement,int *varlist);

/* psn_duplicate():
   Duplicates the sequence 'seq_to_be_duplicated'. If the duplication cannot
   be done, NULL is returned and the type of the error is indicated by 
   'psnerrno'.								     */
psn *	psn_duplicate(psn *seq_to_be_duplicated);

/* psn_concate():
   Concatenates the sequences 'seq2' to the sequence 'seq1'. After this call, 
   the sequence 'seq1' will contain `psn_get_nseq(seq1)+psn_get_nseq(seq2)` 
   sub-sequences. If the concatenation was successful, the function returns 0,
   otherwise it returns a non-zero value and the type of the error can be 
   figured out by using the PE* values.					     */
int	psn_concate(psn *seq1,psn *seq2);

/* psn_exctract():
   Extracts the nth (0th is the first) sub-sequence from the sequence 'seq' 
   even if it was previously optimized. A pointer to the newly created
   sequence is returned, if the extraction was successful, otherwise it 
   returns NULL and 'psnerrno' indicates the type of the error.		     */
psn *	psn_extract(psn *seq,int n);
/* psn_full_extract():
   Extinguishes all optimalizations from the sequence 'seq'. If the procedure
   was successful, the function returns a pointer to a newly created
   sequence, otherwise it returns NULL and 'psnerrno' indicates the type of 
   the error occured.							     */
psn *	psn_full_extract(psn *seq);

/* psn_convert_symbolic():
   Converts the PSN sequence 'seq' to a character string, using the conversion
   rules defined in the {0,NULL}-terminated array 'syt'. For the variables, the
   symbols found in 'sym' is used. The numerical values (constants) are 
   inserted into the sequence using sprintf() and the format string 'dformat'. 
   If 'dformat' is NULL, the format string "%g" will be used. The {0,0,0,0}- 
   terminated list of the property-array 'plist' also has to be passed to 
   this function to get the appropriate values of the precedencies of the 
   operators. The symbolic names of the functions and operators might be found 
   in 'sym' are not used in any case. The function returns a newly allocated
   character stream, which can be released by free(). If any step of the
   conversation fails (e.g. some of the required conversion rules are not found 
   in 'syt'), the function returns NULL.
   Note that the 'seq' should contain only one sub-sequence and should not
   be optimized. If it does not stand for the sequence 'seq', the function
   also returns NULL. Use psn_extract() and psn_full_extract() to extract
   any of the sub-sequences from a PSN sequence and/or extinguish the
   optimalizations. If the function returns NULL, the type of the error can
   be figured out by using 'psnerrno' and the PE* values.		     */
char *	psn_convert_symbolic(psn *seq,psnprop *plist,psnsymeval *syt,
	psnsym **sym,char *dformat);

/* psn_argument_chain():
   This function returns a dynamically allocated array which countains
   integers pointing to the operators or functions of which arguments are they.
   If something is wrong with 'seq', the function returns NULL and psnerrno
   is set appropriately. This function is useful to optimalize evaluation
   in some cases. E.g. when one wants to calculate A*B, and A is zero, 
   therefore there is no need to calculate B. The returned array can be
   released by free().							     */
int *	psn_argument_chain(psn *seq);

/* psn_free():
   Releases the PSN sequence 'seq'.					     */
void	psn_free(psn *seq);

/* psn_register_function():
   This function registers a new function (T_FN) with the given name 'name' to 
   the symbol list, property list, function list, differential rule list
   and symbolic evaluation list (named 'syms', 'props', 'functs', 'diffs' and
   'symevals', respectively). Note that these lists should be dynamically
   allocated lists (or NULLs). The new function has 'argnum' arguments,
   evaluated by the function 'funct', has a differential rule described
   by 'diffrule' and can be converted to human-readable string by 
   'symevalstr'. The latter two are optional, can be NULL, in this case
   'diffs' and 'symevals' also can be NULL (and vice versa, if these are NULL,
   'diffrule' and 'symevalstr' are ignored).				     */
int	psn_register_function(psnsym **syms,psnprop **props,
		psnfunct **functs,psndiff **diffs,psnsymeval **symevals,
		int major,char *name,int argnum,int (*funct)(double *),
		short *diffrule,char *symevalstr);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                                       
