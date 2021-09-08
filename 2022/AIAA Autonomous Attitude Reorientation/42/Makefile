# This commit has a bug where it includes the wrong objects for the gui.
# This can be fixed by inserting the new GUIOBJ line below.

ifeq ($(42PLATFORM),__MSYS__)
   ifneq ($(strip $(GUIFLAG)),)
      GUIOBJ = $(OBJ)42gl.o $(OBJ)42glut.o $(OBJ)glkit.o
   endif
endif

# We also have to add OSQP to our linking options.

42 : $(42OBJ) $(GUIOBJ) $(SIMIPCOBJ) $(FFTBOBJ) $(SLOSHOBJ) $(KITOBJ) $(ACOBJ) $(GMSECOBJ) Makefile
	$(CC) $(LFLAGS) -L$(EXTERNDIR)/osqp/ -losqp $(GMSECBIN) -o $(EXENAME) $(42OBJ) $(GUIOBJ) $(FFTBOBJ) $(SLOSHOBJ) $(KITOBJ) $(ACOBJ) $(GMSECOBJ) $(SIMIPCOBJ) $(LIBS) $(GMSECLIB)

# We don't need to link OSQP to the 42fsw.c file, but we do need to include the h-file.

$(OBJ)42fsw.o       : $(SRC)42fsw.c $(INC)Ac.h $(INC)AcTypes.h
	$(CC) $(CFLAGS) -c $(SRC)42fsw.c -o $(OBJ)42fsw.o -I $(SRC)osqp/include/