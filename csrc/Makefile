# Makefile for the C checker kernel.  Flags are part of the trusted base: no fast math, no FMA contraction
# (spec v0.3 section 1.2).  FLINT >= 3.0 with Arb (headers flint/arb.h) is required; set FLINT_PREFIX to a
# conda env or a from-source install if the headers are not on the default path.
CC      ?= cc
FLINT_PREFIX ?=
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -fno-fast-math -ffp-contract=off -fPIC
ifneq ($(FLINT_PREFIX),)
CFLAGS  += -I$(FLINT_PREFIX)/include
LDFLAGS += -L$(FLINT_PREFIX)/lib -Wl,-rpath,$(FLINT_PREFIX)/lib
endif
LDLIBS  := -lflint -lm
BUILD   := ../build

all: $(BUILD)/libkisscheck.so $(BUILD)/test_arith

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/libkisscheck.so: kc_arith.c kc_arith.h kc_model.c kc_model.h kc_replay.c kc_replay.h kc_leaf.c kc_leaf.h | $(BUILD)
	$(CC) $(CFLAGS) -shared -o $@ kc_arith.c kc_model.c kc_replay.c kc_leaf.c $(LDFLAGS) $(LDLIBS)

$(BUILD)/test_arith: test_arith.c kc_arith.c kc_arith.h | $(BUILD)
	$(CC) $(CFLAGS) -o $@ test_arith.c kc_arith.c $(LDFLAGS) $(LDLIBS)

test: $(BUILD)/test_arith
	$(BUILD)/test_arith 1000000

clean:
	rm -rf $(BUILD)
.PHONY: all test clean
