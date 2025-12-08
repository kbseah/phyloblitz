SHELL := /bin/bash
DB=test_run/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta
READS=test_run/test1.sample.fastq.gz
OUTDIR=test_run/pbz_test
PREFIX=pbz

help:
	@echo "make test	Format code and run unit tests"
	@echo "make run	Test pipeline with read dataset"
	@echo "make clean	Delete test run output folder"
	@echo "make pclean	Delete test run output files except initial mapping, to test --resume mode"

test:
	black ./src
	black ./tests
	python3 -m unittest tests.utils
	python3 -m unittest tests.pipeline

run:
	phyloblitz --db ${DB} --outdir ${OUTDIR} --prefix ${PREFIX} --reads ${READS} --threads 12 --log test.log --resume --dv_max_auto --debug

clean:
	rm -rf ${OUTDIR} test.log

pclean:
	rm -f ${OUTDIR}/${PREFIX}_{ava.abc,ava.mci,ava.paf,ava_seq.tab,final.fasta,final_tophits.paf,mapped.fastq,mcl.out,report_dvs_hist.png,report.html,report.json,report.md} test.log
