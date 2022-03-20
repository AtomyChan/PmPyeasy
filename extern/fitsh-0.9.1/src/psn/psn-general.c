/*****************************************************************************/
/* psn-general.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A simple set of general, widely used functions and operators:	     */
/*	+, -, *, /, ^ (power)						     */
/*	sin, cos, tan, ctn (argument is in radians)			     */
/*	sina, cosa, tana, ctna (argument is in degrees)			     */
/*	exp, ln (log), expt, lg (natural and 10-based exp() and log())	     */
/*	sqrt, abs, sign, mod, int, div, pi (some other functions)	     */
/*	==, !=, <, >, <=, >=, &&, || (relation and logical operators)	     */
/*	r, g (functions related to random number generation)		     */
/* These file contaions symbol definitions and evaluating function codes.    */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2002-2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>

#include	<psn/psn.h>
#include	"psn-general.h"

/* Definition of operators, basic symbols and functions, basic precedencies **/
/* and unary operators which have binary symbols also... *********************/

psnsym		psn_general_op[]= {	{ T_OP, O_ADD, "+", TO_INFIX  },
					{ T_OP, 0    , "+", TO_PREFIX },
					{ T_OP, O_SUB, "-", TO_INFIX  },
					{ T_OP, O_CHS, "-", TO_PREFIX },
					{ T_OP, O_MUL, "*", TO_INFIX  },
					{ T_OP, O_DIV, "/", TO_INFIX  },
					{ T_OP, O_POW, "^", TO_INFIX  },

					{ T_OP, O_RCP, "RCP", 0 },
					{ T_OP, O_PSQ, "PSQ", 0 },

					{0,0,NULL,0} };

psnsym		psn_general_xop[]= {	{ T_OP, O_EQU, "==", TO_INFIX },
					{ T_OP, O_NEQ, "!=", TO_INFIX },
					{ T_OP, O_LS , "<" , TO_INFIX },
					{ T_OP, O_GR , ">" , TO_INFIX },
					{ T_OP, O_LE , "<=", TO_INFIX },
					{ T_OP, O_GE , ">=", TO_INFIX },
					{ T_OP, O_AND, "&&", TO_INFIX },
					{ T_OP, O_OR , "||", TO_INFIX },
					{0,0,NULL,0} };

psnsym		psn_general_fn[]= {	
					{ T_FN, F_SIN  , "sin" ,1 },
					{ T_FN, F_COS  , "cos" ,1 },
					{ T_FN, F_TAN  , "tg"  ,1 },
					{ T_FN, F_TAN  , "tan" ,1 },
					{ T_FN, F_CTN  , "ctg" ,1 },
					{ T_FN, F_CTN  , "ctn" ,1 },

					{ T_FN, F_SINA , "sina",1 },
					{ T_FN, F_COSA , "cosa",1 },
					{ T_FN, F_TANA ,  "tga",1 },
					{ T_FN, F_TANA , "tana",1 },
					{ T_FN, F_CTNA , "ctga",1 },
					{ T_FN, F_CTNA , "ctna",1 },

					{ T_FN, F_ASIN, "arcsin",1 },
					{ T_FN, F_ASIN, "asin",1 },
					{ T_FN, F_ACOS, "arccos",1 },
					{ T_FN, F_ACOS, "acos",1 },
					{ T_FN, F_ATAN, "arctg" ,1 },
					{ T_FN, F_ATAN, "atan",1 },
					{ T_FN, F_ACTN, "arcctg",1 },
					{ T_FN, F_ACTN, "actan",1 },

					{ T_FN, F_ASINA, "arcsina",1 },
					{ T_FN, F_ASINA, "asina",1 },
					{ T_FN, F_ACOSA, "arccosa",1 },
					{ T_FN, F_ACOSA, "acosa",1 },
					{ T_FN, F_ATANA, "arctga" ,1 },
					{ T_FN, F_ATANA, "atana",1 },
					{ T_FN, F_ACTNA, "arcctga",1 },
					{ T_FN, F_ACTNA, "actana",1 },

					{ T_FN, F_EXP  , "exp" ,1 },
					{ T_FN, F_LOG  , "ln"  ,1 },
					{ T_FN, F_LOG  , "log" ,1 },
					{ T_FN, F_EXP10, "expt",1 },
					{ T_FN, F_LOG10, "lg"  ,1 },
					{ T_FN, F_EXP10, "exp10",1 },
					{ T_FN, F_LOG10, "log10",1 },

					{ T_FN, F_SQR  , "sqrt",1 },
					{ T_FN, F_ABS  , "abs" ,1 },
					{ T_FN, F_SGN  , "sign",1 },
					{ T_FN, F_FMOD , "mod" ,2 },
					{ T_FN, F_FINT , "int" ,1 },
					{ T_FN, F_FDIV , "div" ,2 },
					{ T_FN, F_VPI  , "pi"  ,0 },
					{ T_FN, F_ARG  , "arg" ,2 },

					{ 0,0,NULL,0 } };

psnsym		psn_general_fn_rnd[]= {	
					{ T_FN, F_RND  , "r",2 },
					{ T_FN, F_GAUSS, "g",2 },

					{ 0,0,NULL,0 } };


/* Precedence & associativity definitions ************************************/

psnprop		psn_general_prop[] = {	{ O_ADD, 2, 20, ASSOC_LEFT },
					{ O_SUB, 2, 20, ASSOC_LEFT },
					{ O_MUL, 2, 21, ASSOC_LEFT },
					{ O_DIV, 2, 21, ASSOC_LEFT },
					{ O_CHS, 1, 22 },
					{ O_POW, 2, 23, ASSOC_RIGHT },

					{ O_RCP, 1, 23 },
					{ O_PSQ, 1, 23 },

					{ O_EQU, 2, 10, ASSOC_LEFT },
					{ O_NEQ, 2, 10, ASSOC_LEFT },
					{ O_LS , 2, 10, ASSOC_LEFT },
					{ O_GR , 2, 10, ASSOC_LEFT },
					{ O_LE , 2, 10, ASSOC_LEFT },
					{ O_GE , 2, 10, ASSOC_LEFT },
					{ O_AND, 2,  5, ASSOC_LEFT },
					{ O_OR , 2,  4, ASSOC_LEFT },

					{ F_SIN  ,1, 0 },
					{ F_COS  ,1, 0 },
					{ F_TAN  ,1, 0 },
					{ F_TAN  ,1, 0 },
					{ F_CTN  ,1, 0 },
					{ F_CTN  ,1, 0 },

					{ F_SINA ,1, 0 },
					{ F_COSA ,1, 0 },
					{ F_TANA ,1, 0 },
					{ F_TANA ,1, 0 },
					{ F_CTNA ,1, 0 },
					{ F_CTNA ,1, 0 },

					{ F_EXP  ,1, 0 },
					{ F_LOG  ,1, 0 },
					{ F_EXP10,1, 0 },
					{ F_LOG10,1, 0 },

					{ F_SQR  ,1, 0 },
					{ F_ABS  ,1, 0 },
					{ F_SGN  ,1, 0 },
					{ F_FMOD ,2, 0 },
					{ F_FINT ,1, 0 },
					{ F_FDIV ,2, 0 },
					{ F_VPI  ,0, 0 },
					{ F_ARG  ,2, 0 },

					{ F_RND  ,2, 0 },
					{ F_GAUSS,2, 0 },

					{ 0, 0, 0, 0 } };

/* Definition of function codes **********************************************/

static int __fn_o_add(double *s) { *(s-2)=(*(s-2))+(*(s-1));return(0);	}
static int __fn_o_sub(double *s) { *(s-2)=(*(s-2))-(*(s-1));return(0);	}
static int __fn_o_mul(double *s) { *(s-2)=(*(s-2))*(*(s-1));return(0);	}
static int __fn_o_div(double *s) { if ( *(s-1)==0.0 ) return(1);
				   *(s-2)=(*(s-2))/(*(s-1));return(0);	}
static int __fn_o_pow(double *s) { /* if ( *(s-2)<=0.0 ) return(1); */
				   *(s-2)=pow(*(s-2),*(s-1));
				   return(0);				}

static int __fn_o_chs(double *s) { *(s-1)=-(*(s-1));return(0);		}

static int __fn_o_psq(double *s) { s--,(*s)*=*s;return(0);		}
static int __fn_o_rcp(double *s) { s--;if ( *s==0.0 ) return(1);
				   *s=1.0/(*s);return(0);		}

static int __fn_o_sqr(double *s) { s--;if ( *s<0.0 ) return(1);
				   *s=sqrt(*s);return(0);		}
static int __fn_o_abs(double *s) { s--;*s=fabs(*s);return(0);		}
static int __fn_o_sgn(double *s) { s--;	if ( *s>0 ) *s=1.0;
					else if ( *s<0 ) *s=-1.0;
					else *s=0.0;
				   return(0);				}
static int __fn_o_fmod(double *s) { double m,d;int k;
				   s-=2;m=*s,d=*(s+1);
				   if ( d<0.0 ) d=-d; else if ( d==0.0 ) return(1);
				   if ( m<0.0 ) k=(int)((-m)/d),m+=d*(double)(k+2);
				   k=(int)(m/d),m-=(double)k*d;
				   *s=m;return(0);			}
static int __fn_o_fint(double *s) { s--;*s=floor(*s);return(0);		}

static int __fn_o_fdiv(double *s) { double m,d;
				   s-=2;m=*s,d=*(s+1);
				   if ( d<0.0 ) d=-d; else if ( d==0.0 ) return(1);
				   m=floor(m/d);
				   *s=m;return(0);			}

#define			X_PI	3.1415926535897932

static int __fn_o_vpi(double *s) { *s=X_PI;return(0); }

static int __fn_f_arg(double *s) { s-=2; *s=atan2(*(s+1),*s);return(0); }


static int __fn_f_exp(double *s) { s--;*s=exp(*s);return(0);		}
static int __fn_f_log(double *s) { s--;if ( *s<=0.0 ) return(1);
				   *s=log(*s);return(0);		}
static int __fn_f_exp10(double *s) { s--;*s=exp(*s*(M_LN10));return(0);	}
static int __fn_f_log10(double *s) { s--;if ( *s<=0.0 ) return(1);
				   *s=log(*s)*M_LOG10E;return(0);	}

static int __fn_f_sin(double *s) { *(s-1)=sin(*(s-1));return(0);	}
static int __fn_f_cos(double *s) { *(s-1)=cos(*(s-1));return(0);	}
static int __fn_f_tan(double *s) { *(s-1)=tan(*(s-1));return(0);	}
static int __fn_f_ctn(double *s) { *(s-1)=1.0/tan(*(s-1));return(0);	}

static int __fn_f_sina(double *s) { *(s-1)=sin(*(s-1)*X_PI/180.0);return(0);	 }
static int __fn_f_cosa(double *s) { *(s-1)=cos(*(s-1)*X_PI/180.0);return(0);	 }
static int __fn_f_tana(double *s) { *(s-1)=tan(*(s-1)*X_PI/180.0);return(0);	 }
static int __fn_f_ctna(double *s) { *(s-1)=1.0/tan(*(s-1)*X_PI/180.0);return(0); }

static int __fn_f_asin(double *s) { *(s-1)=asin(*(s-1));return(0); }
static int __fn_f_acos(double *s) { *(s-1)=acos(*(s-1));return(0); }
static int __fn_f_atan(double *s) { *(s-1)=atan(*(s-1));return(0); }
static int __fn_f_actn(double *s) { *(s-1)=atan(1.0/(*(s-1)));return(0); }

static int __fn_f_asina(double *s) { *(s-1)=asin(*(s-1))*180.0/X_PI;return(0); }
static int __fn_f_acosa(double *s) { *(s-1)=acos(*(s-1))*180.0/X_PI;return(0); }
static int __fn_f_atana(double *s) { *(s-1)=atan(*(s-1))*180.0/X_PI;return(0); }
static int __fn_f_actna(double *s) { *(s-1)=atan(1.0/(*(s-1)))*180.0/X_PI;return(0); }

/* Definitions of additional operators: relations and boolean operators: *****/

static int __fn_o_equ(double *s) { if ( *(s-1)==*(s-2) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_neq(double *s) { if ( *(s-1)!=*(s-2) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_ls(double *s)  { if ( *(s-2)< *(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_gr(double *s)  { if ( *(s-2)> *(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_le(double *s)  { if ( *(s-2)<=*(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_ge(double *s)  { if ( *(s-2)>=*(s-1) ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_and(double *s) { if ( *(s-1)!=0 && *(s-2)!=0 ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }
static int __fn_o_or(double *s)  { if ( *(s-1)!=0 || *(s-2)!=0 ) *(s-2)=1.0; else *(s-2)=0.0; return(0); }

/* Definitions of functions related to generate random numbers: **************/

static int __fn_f_rnd  (double *s) { *(s-2)=*(s-2)+(*(s-1)-*(s-2))*drand48();return(0);	}
static int __fn_f_gauss(double *s) { *(s-2)=*(s-2)+(*(s-1))*sqrt(-2.0*log(drand48()))*cos(2.0*X_PI*drand48());return(0);	}

/*****************************************************************************/

psnfunct	psn_general_funct[] = {	{ O_ADD,__fn_o_add},
					{ O_SUB,__fn_o_sub},
					{ O_MUL,__fn_o_mul},
					{ O_DIV,__fn_o_div},
					{ O_POW,__fn_o_pow},

					{ O_CHS,__fn_o_chs},
					{ O_PSQ,__fn_o_psq},
					{ O_RCP,__fn_o_rcp},

					{ F_SQR,__fn_o_sqr},
					{ F_ABS,__fn_o_abs},
					{ F_SGN,__fn_o_sgn},
					{ F_FMOD,__fn_o_fmod},
					{ F_FINT,__fn_o_fint},
					{ F_FDIV,__fn_o_fdiv},
					{ F_VPI,__fn_o_vpi},

					{ F_ARG,__fn_f_arg},
		
					{ F_EXP,__fn_f_exp},
					{ F_LOG,__fn_f_log},
					{ F_EXP10,__fn_f_exp10},
					{ F_LOG10,__fn_f_log10},

					{ F_SIN,__fn_f_sin},
					{ F_COS,__fn_f_cos},
					{ F_TAN,__fn_f_tan},
					{ F_CTN,__fn_f_ctn},

					{ F_SINA,__fn_f_sina},
					{ F_COSA,__fn_f_cosa},
					{ F_TANA,__fn_f_tana},
					{ F_CTNA,__fn_f_ctna},

					{ F_ASIN ,__fn_f_asin },
					{ F_ACOS ,__fn_f_acos },
					{ F_ATAN ,__fn_f_atan },
					{ F_ACTN ,__fn_f_actn },

					{ F_ASINA,__fn_f_asina},
					{ F_ACOSA,__fn_f_acosa},
					{ F_ATANA,__fn_f_atana},
					{ F_ACTNA,__fn_f_actna},

					{ O_EQU,__fn_o_equ},
					{ O_NEQ,__fn_o_neq},
					{ O_LS ,__fn_o_ls},
					{ O_GR ,__fn_o_gr},
					{ O_LE ,__fn_o_le},
					{ O_GE ,__fn_o_ge},
					{ O_AND,__fn_o_and},
					{ O_OR ,__fn_o_or},

					{ F_RND  ,__fn_f_rnd  },
					{ F_GAUSS,__fn_f_gauss},

					{0,NULL} };

/*****************************************************************************/
                                    
