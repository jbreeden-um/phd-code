SUPPORT = constraints.c utilities.c constants.c points.c verify.c

algorithm : $(SUPPORT) algorithm.c verify.h
	gcc algorithm.c $(SUPPORT) -o algorithm -lm -g -Wall

setup : $(SUPPORT) setup.c verify.h
	gcc setup.c $(SUPPORT) -o setup -lm -g -Wall

expert : $(SUPPORT) expert.c verify.h
	gcc expert.c $(SUPPORT) -o expert -lm -g -Wall