# The following requires the the SRA-toolkit v2.8.2 to be installed.
for i in `cat SRR_accessionsbatch01.txt`
do
	fastq-dump -I --split-files ${i}
done

# The downloaded files do not contain the index in a read of a particular orientation.
# It is, therefore, necessary to combine them in a single file per read and then demultiplex them using our own script.
cat *_1.fastq | gzip > forward.fastq.gz
cat *_2.fastq | gzip > reverse.fastq.gz
