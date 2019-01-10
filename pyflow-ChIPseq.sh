#! /bin/bash


# invoke Snakemake in cluster mode with custom wrapper scripts

snakemake --rerun-incomplete -j 1000 --jobscript ./jobscript.sh \
		  --latency-wait 120  \
		  --cluster-config cluster.json --cluster './grun_cluster.py {dependencies}' \
		  "$@"
