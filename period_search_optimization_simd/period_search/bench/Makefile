# Banco aislado de CalcStrategyAsimd::bright(), para medir cambios en la ruta
# ASIMD sin tener que ejecutar un workunit completo.
#
#   make                          -> bench_bright + fdiv_probe (cross, aarch64)
#   make host                     -> bench_host (x86, ruta escalar; sirve para
#                                    validar el arnes, no para medir SIMD)
#   make EXTRA=-DPS_FAST_RECIP    -> variante alternativa de bright_asimd.cpp
#
# Las rutas siguen la convencion de Makefile_OpenWRT_aarch64 y se pueden
# sobreescribir desde la linea de ordenes:
#
#   make OPENWRT_DIR=/ruta/al/sdk
#
# Requiere libps_common.a construida antes:
#   cd ../../../period_search_common && make -f Makefile_OpenWRT_aarch64

OPENWRT_DIR ?= ../../../../openwrt
GCC_VERSION ?= 14.3.0
STAGING_DIR ?= $(OPENWRT_DIR)/staging_dir/toolchain-aarch64_cortex-a53_gcc-$(GCC_VERSION)_musl

PS_DIR        ?= ..
CM_DIR        ?= ../../../period_search_common
PS_COMMON_LIB ?= $(CM_DIR)/build/libps_common.a

CXX = $(STAGING_DIR)/bin/aarch64-openwrt-linux-musl-g++
CC  = $(STAGING_DIR)/bin/aarch64-openwrt-linux-musl-gcc

CXXFLAGS = -march=armv8-a -O3 -std=c++2a -DARM64 -I$(CM_DIR) -I$(PS_DIR)
EXTRA   ?=

# Fuentes de la app necesarias para levantar el estado que bright() consume.
# mrqcof_asimd y gauss_errc_asimd entran solo porque la vtable de
# CalcStrategyAsimd se emite en la TU de su primera funcion virtual (mrqcof).
SRCS = bright_asimd curv_asimd conv_asimd mrqcof_asimd gauss_errc_asimd \
       blmatrix matrix phasec areanorm sphfunc trifac
OBJS = $(SRCS:%=obj/%.o) obj/bench_bright.o

all: bench_bright fdiv_probe

obj/%.o: $(PS_DIR)/%.cpp
	@mkdir -p obj
	STAGING_DIR=$(STAGING_DIR) $(CXX) $(CXXFLAGS) $(EXTRA) -c $< -o $@

obj/bench_bright.o: bench_bright.cpp
	@mkdir -p obj
	STAGING_DIR=$(STAGING_DIR) $(CXX) $(CXXFLAGS) -c $< -o $@

bench_bright: $(OBJS)
	STAGING_DIR=$(STAGING_DIR) $(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(PS_COMMON_LIB) -lm

fdiv_probe: fdiv_probe.c
	STAGING_DIR=$(STAGING_DIR) $(CC) -march=armv8-a -O2 -o $@ $<

# --- build de validacion en el host (ruta escalar, sin NEON) ---
HOST_SRCS = bright curv conv mrqcof gauss_errc blmatrix matrix phasec \
            areanorm sphfunc trifac
HOST_OBJS = $(HOST_SRCS:%=objhost/%.o) objhost/dot_product.o \
            objhost/arrayHelpers.o objhost/bench_bright.o

host: bench_host

objhost/%.o: $(PS_DIR)/%.cpp
	@mkdir -p objhost
	$(HOSTCXX) -O2 -std=c++2a -I$(CM_DIR) -I$(PS_DIR) -c $< -o $@

objhost/%.o: $(CM_DIR)/%.cpp
	@mkdir -p objhost
	$(HOSTCXX) -O2 -std=c++2a -I$(CM_DIR) -I$(PS_DIR) -c $< -o $@

objhost/bench_bright.o: bench_bright.cpp
	@mkdir -p objhost
	$(HOSTCXX) -O2 -std=c++2a -I$(CM_DIR) -I$(PS_DIR) -c $< -o $@

bench_host: $(HOST_OBJS)
	$(HOSTCXX) -o $@ $(HOST_OBJS) -lm

HOSTCXX ?= g++

clean:
	rm -rf obj objhost bench_bright bench_host fdiv_probe

.PHONY: all host clean
