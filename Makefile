DIRS=expressiongen matrixlib buildlgf

all:
	-for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS)) ; done

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean) ; done
