%{
#include <iostream>
#include <fstream>
#include <stack>
#include <map>
#include <string>
#include <vector>
#include "Expression.hh"
#include "SumExpression.hh"
#include "DynMatExpression.hh"
#include "InvertExpression.hh"
#include "ProjectorExpression.hh"
#include "PlusMinusExpression.hh"
#define DEBUG 0

extern "C"
{
	int yyparse(void);
	int yylex(void);
	int yywrap()
	{
		return 1;
	}
}


void yyerror(const char *str)
{
	std::cerr << "error: " << str << std::endl;
}

std::stack <Expression *> expstack;
std::map <std::string, Expression *> symbol_table;

#define READ_BUF_SIZE 2048

void writeProgramHeader(std::ostream& os, char *funcname) {
	os << "#include \"Matrix.hh\"\n"
	   << "#include \"DynMat.hh\"\n"
	   << "Matrix " << funcname
	   << "(DynMat &dynmat, double k[3]) {\n"
	   << "\tMatrix result(3u*dynmat.getNions(),3u*dynmat.getNions());\n"
	   << "\tresult = ";
}

void writeProgramFooter(std::ostream& os) {
	os << "\treturn result;\n" << "}" << std::endl;
}

int main()
{
	yyparse();
	return 0;
}


%}

%union
{
	int number;
	char *string;
	char charac;
}

%token OPENP CLOSEP DEF INV SUM EVAL

%token <string> SYMBOL
%token <string> FUNCNAME
%token <charac> PLUSMINUS
%token <charac> PROJECTOR
%token <number> DYNMAT

%type <number> commands

%%

program:
	| program define
	| program eval
	;

eval:
	OPENP EVAL FUNCNAME command CLOSEP
	{
		Expression *exp = expstack.top();
		expstack.pop();
#if defined(DEBUG) && DEBUG
		std::cerr << "Evaluation Statement: " << $3 << std::endl;
		std::cerr << "\tCommand: " << exp->getValue() << std::endl;
#endif
		int len = strlen($3);
		/* strip off the quotation marks */
		char *tmp = $3 + 1;
		tmp[len-2] = '\0'; 
		writeProgramHeader(std::cout, tmp);
		std::cout << exp->getValue() << ";\n";
		writeProgramFooter(std::cout);
		free($3);
	}
	;

define:
	OPENP DEF SYMBOL command CLOSEP
	{
		Expression *exp = expstack.top();
		expstack.pop();
#if defined(DEBUG) && DEBUG
		std::cerr << "Define statement: " << $3 << std::endl;
		std::cerr << "\tCommand: " << exp->getExpr() << std::endl;
#endif
		symbol_table[$3] = exp;
		free($3);
	}
	;

commands: command
	  {
	  	$$ = 1;
	  }
	| commands command
	  {
	  	$$ = $1 + 1;
	  }
	;

command:
	symbol
	|
	sum
	|
	invert
	|
	prod
	|
	dynmat
	|
	projector
	;

symbol:
	SYMBOL
	{
		expstack.push(symbol_table[$1]);
#if defined(DEBUG) && DEBUG
		std::cerr << "Symbol: " << $1 << "\tExpression: "
		          << expstack.top()->getExpr() << std::endl;
#endif
		free($1);
	}
	;

sum:
	OPENP SUM commands CLOSEP
	{
		std::vector <Expression *> v($3);
		for(int i=$3-1; i>=0; --i) {
			v[i] = expstack.top();
			expstack.pop();
		}
		expstack.push(new SumExpression(v));
#if defined(DEBUG) && DEBUG
		std::cerr << expstack.top()->getExpr() << std::endl;
#endif
	}
	;

invert:
	OPENP INV command CLOSEP
	{
		Expression *exp = new InvertExpression(expstack.top());
		expstack.pop();
		expstack.push(exp);
#if defined(DEBUG) && DEBUG
		std::cerr << expstack.top()->getExpr() << std::endl;
#endif
	}
	;



prod:
	OPENP PLUSMINUS commands CLOSEP
	{
		std::vector <Expression *> v($3);
		for(int i=$3-1; i>=0; --i) {
			v[i] = expstack.top();
			expstack.pop();
		}
		expstack.push(new PlusMinusExpression($2, v));
#if defined(DEBUG) && DEBUG
		std::cerr << expstack.top()->getExpr()
		          << std::endl;
#endif
	}
	;

dynmat:
	DYNMAT
	{
		expstack.push(new DynMatExpression($1));
#if defined(DEBUG) && DEBUG
		std::cerr << expstack.top()->getExpr() << std::endl;
#endif
	}
	;

projector:
	PROJECTOR
	{
		expstack.push(new ProjectorExpression($1));
#if defined(DEBUG) && DEBUG
		std::cerr << expstack.top()->getExpr() << std::endl;
#endif
	}
	;
