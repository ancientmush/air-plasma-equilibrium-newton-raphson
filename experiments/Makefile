FC = ifx
FFLAGS =  -O0 -g -check all -traceback -fpe0 -fpp
CPPFLAGS =

ifdef CUSTOM
    ifeq ($(A),)
        $(error "Error: no dummy argument for preprocessor are given")
    endif

    CPPFLAGS += -DUSE_CUSTOM_VALUE -DVAL_NO=$(A) -DVAL_NN=$(B)
endif

TARGET = a.out
SRC = filename_manipulation.f90

all: $(TARGET)
	./a.out

$(TARGET): $(SRC)
	$(FC) $(CPPFLAGS) $(FFLAGS) -o $@ $<

clean:
	rm -rf $(TARGET)

