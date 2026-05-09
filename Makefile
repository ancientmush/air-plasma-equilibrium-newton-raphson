FC = ifx
LDFLAGS = -qmkl

SRCDIR = src
BUILDDIR = build
BINDIR = bin

BUILD ?= release

ifeq ($(BUILD), debug)
    FFLAGS = -O0 -g -check all -traceback -fpe0 -qmkl -module $(BUILDDIR) -I$(BUILDDIR)
else
    FFLAGS = -O3 -Warn all -xHOST -qmkl -module $(BUILDDIR) -I$(BUILDDIR)
endif


SRCS = $(SRCDIR)/mod_constants.f90 $(SRCDIR)/mod_types.f90 $(SRCDIR)/mod_types_functions.f90 $(SRCDIR)/mod_types_subroutines.f90 \
	 $(SRCDIR)/mod_newton_raphson.f90 $(SRCDIR)/main.f90

OBJS = $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.o, $(SRCS))

TARGET = $(BINDIR)/calc_pressure

all: $(TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDDIR) $(BINDIR):
	mkdir -p $@


$(BUILDDIR)/main.o: $(BUILDDIR)/mod_constants.o $(BUILDDIR)/mod_types.o $(BUILDDIR)/mod_newton_raphson.o 
$(BUILDDIR)/mod_types_functions.o: $(BUILDDIR)/mod_types.o
$(BUILDDIR)/mod_types_subroutines.o: $(BUILDDIR)/mod_types.o



.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(BINDIR)

