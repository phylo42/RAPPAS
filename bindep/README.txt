##########################################
# RAPPAS ancestral reconstruction
# pre-compiled versions or PhyML and PAML 
##########################################

These are pre-compiled versions of PhyML 3.3 and PAML 4.9.
Those were compiled with static libraries, which should allow you to launch them on most unix architectures.
Their code was slightly modified to remove lots of outputs (disk access) which are not required by RAPPAS.

If you work on MAC OS or Windows, please download the last version of these software and do your own compilation.

PhyML is much faster than PAML but may require large amount of RAM for large trees.
PAML is slower but requires limited amounts of RAM.

!!! All versions are based on x64 libraries !!!

phyml           |   PhyML 3.3 without CPU instructions (slow but highest compatibility).
phyml_SSE       |   PhyML 3.3 using SSE CPU instructions (faster, available on most modern CPUs).
phyml_AVX       |   PhyML 3.3 using AVX CPU instructions (fastest, compatible with most recent multicore Intel CPUs).
baseml/codeml   |   binaries of PAML (slower than PhyML: codeml for DNA analysis, baseml for amino acid analysis)

