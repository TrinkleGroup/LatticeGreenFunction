include ../compilervars.mk

INCLUDES=-I/sw/include
LEX=flex
YACC=bison
YFLAGS=-d

all: gfparser-bin allmats.cc

.y.cc:
	$(YACC) $(YFLAGS) -o $@ $<

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

lex.yy.c: gfparser.l
	$(LEX) $(LFLAGS) $<

lex.yy.o: lex.yy.c gfparser.cc

gfparser.cc: gfparser.y

SumExpression.o: Expression.o
DynMatExpression.o: Expression.o
ProjectorExpression.o: Expression.o
InvertExpression.o: Expression.o
PlusMinusExpression.o: Expression.o

OBJS=lex.yy.o gfparser.o Expression.o SumExpression.o DynMatExpression.o \
     ProjectorExpression.o InvertExpression.o PlusMinusExpression.o

gfparser-bin: $(OBJS)
	$(CXX) $(CXXFLAGS) -o gfparser $(OBJS)

allmats.cc: allmats.gfc gfparser-bin
	./gfparser < allmats.gfc > ../matrixlib/allmats.cc 2> /dev/null

clean:
	-rm *.o lex.yy.c gfparser.cc gfparser.hh gfparser ../matrixlib/allmats.cc
