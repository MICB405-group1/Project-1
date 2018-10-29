
#SETTINGS

base=~/"PData"
bwa_sub="BWA_out"
bwa_path=$base/$bwa_sub

bam_sub="BAM_out"
bam_path=$base/$bam_sub

bcf_sub="BCF_out"
bcf_path=$base/$bcf_sub

phy_sub="Phyl_out"
phy_path=$base/$phy_sub

vcf_file="VCF_out"


proj_path=/projects/micb405/resources/project_1
ref_file=ref_genome.fasta
remote_ref_file=$proj_path/$ref_file
threads=2 # can change the amount of threads being used here


echo "Please only say yes to the following if you have already created the directory in question"
#SETUP
#input, base -> create base directory || remove current and create base directory
echo ">>>> Creating a working root directory... Would you like to remove the existing directory? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
        if [ -d $base ]; then
			echo "  >>  Removing $base directory"
			rm -Rf $base
		else
			echo "  >>  The $base directory does not exist"
		fi
		echo "    --> Creating $base directory"
		mkdir $base
else
        if [ -d $base ]; then
			echo "  >>  Using the existed $base directory"
		else
			echo "  >>  The $base directory does not exist"
			echo "    --> Creating $base directory"
			mkdir $base
		fi
fi

#input, bwa_path, bwa_sub -> create BWA_out directory || remove current and create BWA_out
echo ">>>> Creating a BWA_out directory... Would you like to remove the existing directory? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
        if [ -d $bwa_path ]; then
			echo "  >> Removing $bwa_path directory"
			rm -Rf $bwa_path
		else
			echo "  >> The $bwa_path directory does not exist"
		fi
		echo "    --> Creating $bwa_sub directory"
		mkdir $bwa_path
else
        if [ -d $bwa_path ]; then
			echo "  >> Using the existed $bwa_path directory"
		else
			echo "  >> The $bwa_path directory does not exist"
			echo "    --> Creating $bwa_sub directory"
			mkdir $bwa_path
		fi
fi

#input, bam_path, bam_sub -> create BAM_out directory || remove current and create BAM_out
echo ">>>> Creating a BAM_out directory... Would you like to remove the existing directory? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
        if [ -d $bam_path ]; then
			echo "  >> Removing $bam_path directory"
			rm -Rf $bam_path
		else
			echo "  >> The $bam_path directory does not exist"
		fi
		echo "    --> Creating $bam_sub directory"
		mkdir $bam_path
else
        if [ -d $bam_path ]; then
			echo "  >> Using the existed $bam_path directory"
		else
			echo "  >> The $bam_path directory does not exist"
			echo "    --> Creating $bam_sub directory"
			mkdir $bam_path
		fi
fi

#input, bcf_path, bcf_sub -> create BCF_out directory || remove current and create BCF_out
echo ">>>> Creating a  BCF_out directory... Would you like to remove the existing directory? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
        if [ -d $bcf_path ]; then
			echo "  >> Removing $bcf_path directory"
			rm -Rf $bcf_path
		else
			echo "  >> The $bcf_path directory does not exist"
		fi
		echo "    --> Creating $bcf_sub directory"
		mkdir $bcf_path
else
        if [ -d $bcf_path ]; then
			echo "  >> Using the existed $bcf_path directory"
		else
			echo "  >> The $bcf_path directory does not exist"
			echo "    --> Creating $bcf_sub directory"
			mkdir $bcf_path
		fi
fi

#input, phy_path, phy_sub -> create Phyl_out directory || remove current and create Phyl_out
echo ">>>> Creating a Phyl_out directory... Would you like to remove the existing directory? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
        if [ -d $phy_path ]; then
			echo "  >> Removing $phy_path directory"
			rm -Rf $phy_path
		else
			echo "  >> The $phy_path directory does not exist"
		fi
		echo "    --> Creating $phy_sub directory"
		mkdir $phy_path
else
        if [ -d $phy_path ]; then
			echo "  >> Using the existed $phy_path directory"
		else
			echo "  >> The $phy_path directory does not exist"
			echo "    --> Creating $phy_sub directory"
			mkdir $phy_path
		fi
fi



#INDEXING
#lines 142-150 will first check to see if the user wants to index the reference  genome,
#then depending on their answer it will copy the file to the bwa_path which we created earlier
#and it will call bwa index on the file
echo ">>>> Index the reference ref_genome.fasta? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
	echo "  >> Indexing ref_genome.fasta"
	cp $remote_ref_file $bwa_path/
	bwa index $bwa_path/$ref_file
else
	echo "  >> Using the existing ref_genome.fasta index"
fi




#UNZIPPING, ALIGNING, AND BAM
#lines 159-164 loop through every pair of files in our prof_path (fwd and rev reads) and then
#align them with bwa mem and then pipe the output of that into samtools and create a BAM file
#output for each pair
for f in $proj_path/*_1.fastq.gz; do
    STEM=$(basename "${f}" | sed 's/_1.fastq.gz//g')

    bwa mem -t $threads $bwa_path/$ref_file $proj_path/$STEM\_1.fastq.gz $proj_path/$STEM\_2.fastq.gz | \
    samtools view --threads $threads -bS -F 4 - >$bam_path/$STEM.bam
done



#SORTING
#lines 172-177 loop through each BAM files we created in the previous steps (lines 159-164)
#and calls samtools sort on them and saves each one with the same basename but the suffix .sorted
echo "  >> Sorting the BAM files and storing them in $bam_path"
for f in $bam_path/*.bam; do
		echo "    --> Sorting: $f"
		STEM=$(basename "${f}" .bam)
		samtools sort --threads $threads "${f}" -o $bam_path/"${STEM}".sorted
done



#REMOVING DUPLICATES
#lines 184-189 loop through each of the sorted BAM files we created in the previous lines 172-177
#and remove duplicates from them using samtools rmdup paramter, they are them saved
#with the same basename as before with the new suffix .rmdup added after .sorted
echo "  >> Removing duplicates from the BAM files"
for f in $bam_path/*.sorted; do
		echo "    --> Removing duplicates: $f"
		STEM=$(basename "${f}" .sorted)
		samtools rmdup "${f}" $bam_path/"${STEM}".sorted.rmdup
done


#INDEXING BAM
#lines 196-200 loop through each of the duplicate removed, and sorted BAM files we have thus far
#and indexes each of them
echo "  >> Indexing the BAM files in $bam_path"
for f in $bam_path/*.sorted.rmdup; do
		echo "    --> Indexing: $f"
		samtools index -@ $threads "${f}"
done


#VARIANT CALLING
#lines 206-211 loop through each of the indexed bam files from the previous step and calls variants
#on them by using bcftools with the call bcftools mpileup, saving them to the ouput of the same
#basename but with the new suffix .vcf
echo "  >> Calling variants on the BAM files in $bam_path"
for f in $bam_path/*.sorted.rmdup; do
		echo "    --> Variant calling: $f"
		STEM=$(basename "${f}" .sorted.rmdup)
		bcftools mpileup --threads $threads -q 30 --fasta-ref $bwa_path/$ref_file "${f}" | bcftools call --threads $threads -O v -mv  | bcftools filter --exclude "QUAL < 200" -o $bcf_path/"${STEM}".vcf
done



#PARSE VCF
#the next two lines use the python code that was provided to us to parse the VCF files we just
#created
echo "  >> Parsing VCF files"
python vcf_to_fasta_het.py -x $bcf_path/ $vcf_file









#MSA WITH MUSCLE
#lines 232-240 prompt the user to see if they would like to create a new MSA with muscle. If they
#say no we warn them they will need an already created MSA in the $phy_path directory. If they
#say yes then we use the vcf_file we just created and use muscle and trimal to create the MSA
echo ">>>> Use MUSCLE? If no is selected then there must already be an MSA in $phy_path. [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
	echo "  >> Creating MSA with muscle and storing the result in $phy_sub directory"
	muscle -in $bcf_path/"${vcf_file}"_x.fasta -out $phy_path/"${vcf_file}"_muscle.mfa
    trimal -automated1 -in $phy_path/"${vcf_file}"_muscle.mfa -out $phy_path/"${vcf_file}"_muscle_t.mfa
else
	echo "  >> Using the existing MSA"
fi


#FASTTREE
echo ">>>> Build tree? [Y,n]"
read input
if [[ $input == "Y" || $input == "y" ]]; then
	echo "  >> Building tree from the $phy_sub directory..."
	FastTree $phy_path/"${vcf_file}"_muscle_t.mfa 1>$phy_path/"${vcf_file}"_muscle_t.nwk
else
	echo "  >> Using the existing tree from $phy_path"
fi
