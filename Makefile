ref.fa: hg38.fa HIVmini.fa PhiX.fa
	cat $^ > $@
	
mapped.sam: ref.fa 2023-12-05_DSC_231205_NB501055_0303_AHTWHJBGXC.fastq
	bwa index $<
	bwa mem -t 8 $^ > $@ 

positive.sam: mapped.sam
	samtools view -q20 -F16 $< > $@

negative.sam: mapped.sam
	samtools view -q20 -f16 $< > $@
#samtools view Reads.bam | gawk '(and(16, $2))' > reverseStrandReads.sam
#samtools view Reads.bam | gawk '(! and(16,$2))' > forwardStrandReads.sam
positive.txt: positive.sam
	echo "chromosome position sequence" > $@
	awk '{print $3, $4, $10}' $^ >> $@

negative.txt: negative.sam
	echo "chromosome position sequence" > $@
	awk '{if ($6 ~ /M/) $4 += substr($6, 1, 2); print $3, $4, $$10}' $^ >> $@

negative_complement.txt: negative.txt
	python reverse_complement_dna.py $^ > $@

data.txt: positive.txt negative_complement.txt
	cat $^ > $@

sites.txt: data.txt 
	python get_sites.py data.txt $@
