%{
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "gfparser.hh"
%}

%%
	/*** Rules ***/
#.*		/* ignore comments */
\(		return OPENP;
\)		return CLOSEP;
sum		return SUM;
inv		return INV;
def		return DEF;
eval		return EVAL;
\"[^"]+\"	yylval.string = strdup(yytext); return FUNCNAME;
[\+\-]		yylval.charac=yytext[0]; return PLUSMINUS;
[0-9]+		yylval.number=atoi(yytext); return DYNMAT;
[AaOoUuLl]	yylval.charac=yytext[0]; return PROJECTOR;
[A-Za-z_][A-Za-z_0-9]+|[B-KMNP-TV-Zb-kmnp-tv-z]	yylval.string = strdup(yytext); return SYMBOL;
[ \t\n]+	/*ignore whitespace*/
.		fprintf(stderr, "Parsing ERROR: Bad symbol: %s\n", yytext); exit(-1);

%%
