AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_CFLAGS = $(OGG_CFLAGS) $(CAIRO_CFLAGS)

EXTRA_DIST = \
	encoder_disabled.c \
	arm-intrinsics/armcpu.h \
	arm-intrinsics/armenc.h \
	arm-intrinsics/armint.h \
	c64x/c64xint.h \
	c64x/c64xdec.h \
	x86/mmxfrag.c \
	x86/mmxidct.c \
	x86/mmxloop.h \
	x86/mmxstate.c \
	x86/sse2idct.c \
	x86/x86cpu.c \
	x86/x86int.h \
	x86/x86state.c \
	x86/mmxencfrag.c \
	x86/mmxfdct.c \
	x86/sse2encfrag.c \
	x86/sse2fdct.c \
	x86/sse2trans.h \
	x86/x86enc.c \
	x86/x86enc.h \
	x86/x86enquant.c \
	x86/x86zigzag.h \
	x86_vc

lib_LTLIBRARIES = libtheoradec.la libtheoraenc.la libtheora.la

if THEORA_DISABLE_ENCODE
encoder_uniq_sources = \
	encoder_disabled.c

encoder_sources = \
	$(encoder_uniq_sources)
else
encoder_uniq_x86_sources = \
	x86/mmxencfrag.c \
	x86/mmxfdct.c \
	x86/sse2encfrag.c \
	x86/x86enquant.c \
	x86/x86enc.c

encoder_uniq_x86_64_sources = \
	x86/sse2fdct.c

encoder_shared_x86_sources = \
	x86/x86cpu.c \
	x86/mmxfrag.c \
	x86/mmxidct.c \
	x86/mmxstate.c \
	x86/sse2idct.c \
	x86/x86state.c

encoder_shared_x86_64_sources =

encoder_uniq_arm_sources = \
	arm-intrinsics/armencfrag.c \
	arm-intrinsics/armenc.c

encoder_shared_arm_sources = \
	arm-intrinsics/armcpu.c \
	arm-intrinsics/armfrag.c \
	arm-intrinsics/armidct.c \
	arm-intrinsics/armloop.c \
	arm-intrinsics/armstate.c

encoder_uniq_wasm_sources =

encoder_shared_wasm_sources = \
	wasm-simd/wav128frag.c

encoder_uniq_arch_sources =
encoder_shared_arch_sources =

if CPU_x86_64
encoder_uniq_arch_sources += \
 $(encoder_uniq_x86_sources) \
 $(encoder_uniq_x86_64_sources)
encoder_shared_arch_sources += \
 $(encoder_shared_x86_sources) \
 $(encoder_shared_x86_64_sources)
endif
if CPU_x86_32
encoder_uniq_arch_sources += $(encoder_uniq_x86_sources)
encoder_shared_arch_sources += $(encoder_shared_x86_sources)
endif
if CPU_arm
encoder_uniq_arch_sources += $(encoder_uniq_arm_sources)
encoder_shared_arch_sources += $(encoder_shared_arm_sources)
endif
if CPU_wasm
encoder_uniq_arch_sources += $(encoder_uniq_wasm_sources)
encoder_shared_arch_sources += $(encoder_shared_wasm_sources)
endif

encoder_uniq_sources = \
	analyze.c \
	fdct.c \
	encfrag.c \
	encapiwrapper.c \
	encinfo.c \
	encode.c \
	enquant.c \
	huffenc.c \
	mathops.c \
	mcenc.c \
	rate.c \
	tokenize.c \
	$(encoder_uniq_arch_sources)

encoder_sources = \
	apiwrapper.c \
	bitpack.c \
	dequant.c \
	fragment.c \
	idct.c \
	info.c \
	internal.c \
	state.c \
	quant.c \
	$(encoder_shared_arch_sources) \
	$(encoder_uniq_sources)

endif

decoder_x86_sources = \
	x86/x86cpu.c \
	x86/mmxidct.c \
	x86/mmxfrag.c \
	x86/mmxstate.c \
	x86/sse2idct.c \
	x86/x86state.c

decoder_arm_sources = \
	arm-intrinsics/armcpu.c \
	arm-intrinsics/armfrag.c \
	arm-intrinsics/armidct.c \
	arm-intrinsics/armloop.c \
	arm-intrinsics/armstate.c

decoder_c64x_sources = \
	c64x/c64xdec.c \
	c64x/c64xfrag.c \
	c64x/c64xidct.c \
	c64x/c64xstate.c

decoder_wasm_sources = \
	wasm-simd/wav128frag.c

decoder_arch_sources =

if CPU_x86_64
decoder_arch_sources += $(decoder_x86_sources)
endif
if CPU_x86_32
decoder_arch_sources += $(decoder_x86_sources)
endif
if CPU_arm
decoder_arch_sources += $(decoder_arm_sources)
endif
if CPU_c64x
decoder_arch_sources += $(decoder_c64x_sources)
endif
if CPU_wasm
decoder_arch_sources += $(decoder_wasm_sources)
endif

decoder_sources = \
	apiwrapper.c \
	bitpack.c \
	decapiwrapper.c \
	decinfo.c \
	decode.c \
	dequant.c \
	fragment.c \
	huffdec.c \
	idct.c \
	info.c \
	internal.c \
	quant.c \
	state.c \
	$(decoder_arch_sources)

noinst_HEADERS = \
	apiwrapper.h \
	bitpack.h \
	dct.h \
	decint.h \
	dequant.h \
	encint.h \
	enquant.h \
	huffdec.h \
	huffenc.h \
	huffman.h \
	internal.h \
	mathops.h \
	modedec.h \
	ocintrin.h \
	quant.h \
	state.h \
	arm-intrinsics/armcpu.h \
	c64x/c64xdec.h \
	c64x/c64xint.h \
	x86/mmxloop.h \
	x86/sse2trans.h \
	x86/x86cpu.h \
	x86/x86enc.h \
	x86/x86int.h \
	x86/x86zigzag.h

libtheoradec_la_SOURCES = \
	$(decoder_sources) \
	Version_script-dec theoradec.exp
libtheoradec_la_LDFLAGS = \
  -version-info @THDEC_LIB_CURRENT@:@THDEC_LIB_REVISION@:@THDEC_LIB_AGE@ \
  @THEORADEC_LDFLAGS@ @CAIRO_LIBS@ \
  -no-undefined

libtheoraenc_la_SOURCES = \
	$(encoder_sources) \
	Version_script-enc theoraenc.exp
libtheoraenc_la_LDFLAGS = \
  -version-info @THENC_LIB_CURRENT@:@THENC_LIB_REVISION@:@THENC_LIB_AGE@ \
  @THEORAENC_LDFLAGS@ $(OGG_LIBS) \
  -no-undefined

libtheora_la_SOURCES = \
	$(decoder_sources) \
	$(encoder_uniq_sources) \
	Version_script theora.exp
libtheora_la_LDFLAGS = \
  -version-info @TH_LIB_CURRENT@:@TH_LIB_REVISION@:@TH_LIB_AGE@ \
  @THEORA_LDFLAGS@ @CAIRO_LIBS@ $(OGG_LIBS) \
  -no-undefined

debug:
	$(MAKE) all CFLAGS="@DEBUG@" 

profile:
	$(MAKE) all CFLAGS="@PROFILE@"

# construct various symbol export list files
.def.exp : defexp.awk
	awk -f defexp.awk $< > $@
