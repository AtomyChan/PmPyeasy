/*****************************************************************************/
/* psn-general-ds.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A simple set of general, widely used functions and operators:	     */
/*	+, -, *, /, ^ (power)						     */
/*	sin, cos, tan, ctn (argument is in radians)			     */
/*	sina, cosa, tana, ctna (argument is in degrees)			     */
/*	exp, ln (log), (natural exp() and log())			     */
/*	sqrt, abs, sign, mod, pi (some other functions)			     */
/* These file contains simplification and differentation rules.		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2002-2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include	<stdio.h>

#include	<psn/psn.h>
#include	"psn-general.h"
#include	"psn-general-ds.h"

/* Derivative rules **********************************************************/

static	short	diffrule_o_add[]={ SS_D2,SS_D1,O_ADD,0 };
static	short	diffrule_o_sub[]={ SS_D2,SS_D1,O_SUB,0 };
static	short	diffrule_o_mul[]={ SS_D2,SS_N1,O_MUL,SS_D1,SS_N2,O_MUL,O_ADD,0 };
static	short	diffrule_o_div[]={ SS_D2,SS_N1,O_MUL,SS_D1,SS_N2,O_MUL,O_SUB,SS_N1,O_PSQ,O_DIV,0 };
static	short	diffrule_o_pow[]={ SS_N2,SS_N1,O_POW,SS_D1,O_MUL,SS_N2,F_LOG,O_MUL,SS_N2,SS_N1,CON_1,O_SUB,O_POW,SS_N1,O_MUL,SS_D2,O_MUL,O_ADD,0 };

static	short	diffrule_o_chs[]={ SS_D1,O_CHS,0 };
static	short	diffrule_o_psq[]={ SS_N1,CON_2,O_MUL,SS_D1,O_MUL,0 };
static	short	diffrule_o_rcp[]={ SS_N1,O_PSQ,O_RCP,O_CHS,SS_D1,O_MUL,0 };

static	short	diffrule_f_sqr[]={ SS_N1,F_SQR,CON_2,O_MUL,O_RCP,SS_D1,O_MUL,0 };
static	short	diffrule_f_abs[]={ SS_N1,F_SGN,SS_D1,O_MUL,0 };
static	short	diffrule_f_sgn[]={ CON_0,0 };
/* static	short	diffrule_f_fmod[]={ SS_D2,SS_N2,SS_N1,F_FMOD,SS_N1,O_DIV,SS_N2,SS_N1,O_PSQ,O_DIV,O_SUB,SS_D1,O_MUL,O_ADD,0 }; */
static	short	diffrule_f_fmod[]={ SS_D2,SS_N2,SS_N1,O_DIV,F_FINT,SS_D1,O_MUL,O_SUB,0 };
static	short	diffrule_f_fint[]={ CON_0,0 };
static	short	diffrule_f_vpi[]={ CON_0,0 };

static	short	diffrule_f_arg[]={ SS_N2,SS_D1,O_MUL,SS_N1,SS_D2,O_MUL,O_SUB,SS_N1,O_PSQ,SS_N2,O_PSQ,O_ADD,O_DIV,0 };

static	short 	diffrule_f_sin[]={ SS_N1,F_COS,SS_D1,O_MUL,0 };
static	short	diffrule_f_cos[]={ SS_N1,F_SIN,O_CHS,SS_D1,O_MUL,0 };
static	short	diffrule_f_tan[]={ SS_N1,F_COS,O_PSQ,O_RCP,SS_D1,O_MUL,0 };
static	short	diffrule_f_ctn[]={ SS_N1,F_SIN,O_PSQ,O_RCP,O_CHS,SS_D1,O_MUL,0 };

static	short 	diffrule_f_sina[]={ SS_N1,F_COS,SS_D1,O_MUL,0 };
static	short	diffrule_f_cosa[]={ SS_N1,F_SIN,O_CHS,SS_D1,O_MUL,0 };
static	short	diffrule_f_tana[]={ SS_N1,F_COS,O_PSQ,O_RCP,SS_D1,O_MUL,0 };
static	short	diffrule_f_ctna[]={ SS_N1,F_SIN,O_PSQ,O_RCP,O_CHS,SS_D1,O_MUL,0 };

static	short	diffrule_f_exp[]={ SS_N1,F_EXP,SS_D1,O_MUL,0 };
static	short	diffrule_f_log[]={ SS_N1,O_RCP,SS_D1,O_MUL,0 };

psndiff	psn_general_diff[]={	{O_ADD ,diffrule_o_add},
				{O_SUB ,diffrule_o_sub},
				{O_MUL ,diffrule_o_mul},
				{O_DIV ,diffrule_o_div},
				{O_POW ,diffrule_o_pow},

				{O_CHS ,diffrule_o_chs},
				{O_PSQ ,diffrule_o_psq},
				{O_RCP ,diffrule_o_rcp},

				{F_SQR ,diffrule_f_sqr},
				{F_ABS ,diffrule_f_abs},
				{F_SGN ,diffrule_f_sgn},
				{F_FMOD,diffrule_f_fmod},
				{F_FINT,diffrule_f_fint},
				{F_VPI ,diffrule_f_vpi},
				{F_ARG ,diffrule_f_arg},

				{F_SIN ,diffrule_f_sin},
				{F_COS ,diffrule_f_cos},
				{F_TAN ,diffrule_f_tan},
				{F_CTN ,diffrule_f_ctn},
				{F_SINA,diffrule_f_sina},
				{F_COSA,diffrule_f_cosa},
				{F_TANA,diffrule_f_tana},
				{F_CTNA,diffrule_f_ctna},

				{F_EXP ,diffrule_f_exp},
				{F_LOG ,diffrule_f_log},
				{0,NULL} };

/* Simplifying rules *********************************************************/
									/* Type */
		
static	short sim_o_add1[]={ S_EXP1,CON_0,O_ADD,0,S_EXP1,0 };		/* expr+0=expr		*/
static	short sim_o_add2[]={ CON_0,S_EXP1,O_ADD,0,S_EXP1,0 };		/* 0+expr=expr		*/

static	short sim_o_sub1[]={ S_EXP1,CON_0,O_SUB,0,S_EXP1,0 };		/* expr-0=expr		*/
static	short sim_o_sub2[]={ CON_0,S_EXP1,O_SUB,0,S_EXP1,O_CHS,0 };	/* 0-expr=-expr 	*/

static	short sim_o_mul1[]={ S_EXP1,CON_1 ,O_MUL,0,S_EXP1,0 };		/* expr*1=expr		*/
static	short sim_o_mul2[]={ CON_1 ,S_EXP1,O_MUL,0,S_EXP1,0 };		/* 1*expr=expr		*/
static	short sim_o_mul3[]={ S_EXP1,CON_0 ,O_MUL,0,CON_0,0  };		/* expr*0=0		*/
static	short sim_o_mul4[]={ CON_0 ,S_EXP1,O_MUL,0,CON_0,0  };		/* 0*expr=0		*/
static	short sim_o_mul5[]={ S_EXP1,S_EXP1,O_MUL,0,S_EXP1,O_PSQ,0 };	/* ex*ex=ex^(2)		*/

static	short sim_o_div1[]={ S_EXP1,CON_1 ,O_DIV,0,S_EXP1,0 };		/* expr/1=expr      	*/
static	short sim_o_div2[]={ CON_1 ,S_EXP1,O_DIV,0,S_EXP1,O_RCP,0 };	/* 1/expr=expr^(-1) 	*/
static	short sim_o_div3[]={ S_EXP1,CON_0 ,O_DIV,0,0 };			/* expr/0=UNDEF!     	*/
static	short sim_o_div4[]={ CON_0 ,S_EXP1,O_DIV,0,CON_0 ,0 };		/* 0/expr=0       	*/

static	short sim_o_pow1[]={ S_EXP1,CON_0,O_POW,0,CON_1 ,0 };		/* ex^0=1 		*/
static	short sim_o_pow2[]={ S_EXP1,CON_1,O_POW,0,S_EXP1,0 };		/* ex^1=ex		*/
static	short sim_o_pow3[]={ S_EXP1,CON_2,O_POW,0,S_EXP1,O_PSQ,0 };	/* ex^2=ex PSQ		*/

static	short sim_addmul_dis1[]={ S_EXP1,S_EXP2,O_MUL,S_EXP1,S_EXP3,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 };	/* a*b+a*c=a*(b+c) */
static	short sim_addmul_dis2[]={ S_EXP2,S_EXP1,O_MUL,S_EXP1,S_EXP3,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 }; /* a*b+c*a=a*(b+c) */
static	short sim_addmul_dis3[]={ S_EXP1,S_EXP2,O_MUL,S_EXP3,S_EXP1,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 }; /* b*a+a*c=a*(b+c) */
static	short sim_addmul_dis4[]={ S_EXP2,S_EXP1,O_MUL,S_EXP3,S_EXP1,O_MUL,O_ADD,0,S_EXP1,S_EXP2,S_EXP3,O_ADD,O_MUL,0 }; /* b*a+c*a=a*(b+c) */
static	short sim_submul_dis1[]={ S_EXP1,S_EXP2,O_MUL,S_EXP1,S_EXP3,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* a*b-a*c=a*(b-c) */
static	short sim_submul_dis2[]={ S_EXP2,S_EXP1,O_MUL,S_EXP1,S_EXP3,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* a*b-c*a=a*(b-c) */
static	short sim_submul_dis3[]={ S_EXP1,S_EXP2,O_MUL,S_EXP3,S_EXP1,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* b*a-a*c=a*(b-c) */
static	short sim_submul_dis4[]={ S_EXP2,S_EXP1,O_MUL,S_EXP3,S_EXP1,O_MUL,O_SUB,0,S_EXP1,S_EXP2,S_EXP3,O_SUB,O_MUL,0 }; /* b*a-c*a=a*(b-c) */
static	short sim_adddiv_dis []={ S_EXP2,S_EXP1,O_DIV,S_EXP3,S_EXP1,O_DIV,O_ADD,0,S_EXP2,S_EXP3,O_ADD,S_EXP1,O_DIV,0 }; /* b/a+c/a=(b+c)/a */
static	short sim_subdiv_dis []={ S_EXP2,S_EXP1,O_DIV,S_EXP3,S_EXP1,O_DIV,O_SUB,0,S_EXP2,S_EXP3,O_SUB,S_EXP1,O_DIV,0 }; /* b/a-c/a=(b-c)/a */
static	short sim_mulpow_dis []={ S_EXP2,S_EXP1,O_POW,S_EXP3,S_EXP1,O_POW,O_MUL,0,S_EXP2,S_EXP3,O_MUL,S_EXP1,O_POW,0 }; /* b^a*c^a=(b*c)^a */
static	short sim_divpow_dis []={ S_EXP2,S_EXP1,O_POW,S_EXP3,S_EXP1,O_POW,O_DIV,0,S_EXP2,S_EXP3,O_DIV,S_EXP1,O_POW,0 }; /* b^a/c^a=(b/c)^a */

static	short sim_tau1[]= { S_EXP1,S_EXP2,O_CHS,O_SUB,0,S_EXP1,S_EXP2,O_ADD,0 };	/* a-(-b)=a+b */

static  short sim_gon1[]= { S_EXP1,F_SIN,O_PSQ,S_EXP1,F_COS,O_PSQ,O_ADD,0,CON_1,0 };	/* sin^2(x)+cos^2(x)=1 */

short	*psn_general_simp[]=
 {	sim_o_add1,sim_o_add2,sim_o_sub1,sim_o_sub2,
	sim_o_mul1,sim_o_mul2,sim_o_mul3,sim_o_mul4,sim_o_mul5,
	sim_o_div1,sim_o_div2,sim_o_div3,sim_o_div4,
	sim_o_pow1,sim_o_pow2,sim_o_pow3,
	sim_addmul_dis1,sim_addmul_dis2,sim_addmul_dis3,sim_addmul_dis4,
	sim_submul_dis1,sim_submul_dis2,sim_submul_dis3,sim_submul_dis4,
	sim_adddiv_dis,sim_subdiv_dis,sim_mulpow_dis,sim_divpow_dis,
	sim_tau1,sim_gon1,
	NULL
 };

/* Symbolic rules ************************************************************/

psnsymeval	psn_general_symeval[] = {

	{ 0    , "(#0)",	  1,0	     },

	{ O_ADD, "#(1)+#(2)",	  0,TO_INFIX },
	{ O_SUB, "#(1)-#[2]",	  0,TO_INFIX },
	{ O_CHS, "(-#(1))",	  1,0        },
	{ O_MUL, "#(1)*#(2)",	  0,TO_INFIX },
	{ O_DIV, "#(1)/#[2]",	  0,TO_INFIX },
	{ O_POW, "#(1)^#(2)",	  1,0	     },
	{ O_RCP, "1.0/#[1]",	  0,TO_INFIX },
	{ O_PSQ, "#(1)^2",	  1,0	     },

	{ F_SIN, "sin(#1)",	  1,0	     },
	{ F_COS, "cos(#1)",	  1,0	     },
	{ F_TAN, "tg(#1)",	  1,0	     },
	{ F_CTN, "ctg(#1)",	  1,0	     },
	{ F_SINA, "sina(#1)",	  1,0	     },
	{ F_COSA, "cosa(#1)",	  1,0	     },
	{ F_TANA, "tga(#1)",	  1,0	     },
	{ F_CTNA, "ctga(#1)",	  1,0	     },

	{ F_ASIN, "asin(#1)",	  1,0	     },
	{ F_ACOS, "acos(#1)",	  1,0	     },
	{ F_ATAN, "atan(#1)",	  1,0	     },
	{ F_ACTN, "actn(#1)",	  1,0	     },
	{ F_ASINA, "asina(#1)",   1,0	     },
	{ F_ACOSA, "acosa(#1)",   1,0	     },
	{ F_ATANA, "atana(#1)",	  1,0	     },
	{ F_ACTNA, "actna(#1)",   1,0	     },

	{ F_EXP, "exp(#1)",   	  1,0	     },
	{ F_LOG, "log(#1)", 	  1,0	     },
	{ F_EXP, "exp10(#1)",     1,0	     },
	{ F_LOG, "log10(#1)", 	  1,0	     },
	{ F_SQR, "sqrt(#1)",	  1,0	     },
	{ F_ABS, "abs(#1)",	  1,0	     },
	{ F_SGN, "sign(#1)",	  1,0	     },
	{ F_FMOD, "mod(#1,#2)",	  2,0	     },
	{ F_FINT, "int(#1)",	  1,0	     },
	{ F_FDIV, "div(#1,#2)",	  2,0	     },

	{ F_VPI, "pi()",	  0,0	     },

	{ 0, NULL,0,0 }

   };

/*****************************************************************************/
                                                                   
