# a few default targets for all makefiles to include

# missing PUT_VERSION because don't always want this
#depend :
#	makedepend -v -f ../tradesdep.mk -- $(CFLAGS1) -I$(INCPATH) $(VIBFLAG) -DX11 $(DEF) -- $(SRC)

$(TRAJLIBNAME) : $(TRAJLIBOBJS)
	-mkdir $(TRADESLIBS)
	ar r $(TRAJLIBNAME) $(TRAJLIBOBJS)
	ranlib $(TRAJLIBNAME) || true

$(FOLDTRAJLIBNAME) : $(FOLDTRAJLIBOBJS)
	-mkdir $(TRADESLIBS)
	ar r $(FOLDTRAJLIBNAME) $(FOLDTRAJLIBOBJS)
	ranlib $(FOLDTRAJLIBNAME) || true

$(MAKETRAJLIBNAME) : $(MAKETRAJLIBOBJS)
	-mkdir $(TRADESLIBS)
	ar r $(MAKETRAJLIBNAME) $(MAKETRAJLIBOBJS)
	ranlib $(MAKETRAJLIBNAME) || true

# DO NOT DELETE THIS LINE -- make  depend  depends  on it.	

$(SRCDIR)/vecmath.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/rotate.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/bbox.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/Rplot.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/rotlib.o : $(TRADESINC)/rotlib.h $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/potential.o : $(TRADESINC)/potential.h $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/crease.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/dft.o : $(TRADESINC)/dft.h
$(SRCDIR)/gil_gor.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/extaaseq.o : $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/randwalk.o : $(TRADESINC)/slriaccsurf.h $(TRADESINC)/geometry.h $(TRADESINC)/foldtrajlib.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/vismaketrj.o : $(TRADESINC)/vismaketrj.h $(TRADESINC)/foldtrajlib.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/maketg.o : $(TRADESINC)/foldtrajlib.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/GetTrjAngle.o : $(TRADESINC)/geometry.h $(TRADESINC)/foldtrajlib.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/trajtools.o : $(TRADESINC)/foldtrajlib.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/buildit.o : $(TRADESINC)/geometry.h $(TRADESINC)/mmdbtraj.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
$(SRCDIR)/newasn.o : $(TRADESINC)/foldtrajlib.h $(TRADESINC)/mmdbtraj_pub.h $(BZINC)/bzlib.h
