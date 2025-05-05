

indir=${1? Input directory with files pos.txt sample_list.txt and mat_vcf.txt}

proj=${2? Output (project) prefix, e.g. ricerp }


#loadmatrix=./bin/loadmatrix_geno.v2.hdf1.10	
loadmatrix=./loadmatrix_geno

n_snp=`wc -l $indir/pos.txt| awk '{print $1}' `

n_sampl=`wc -l $indir/sample_list.txt| awk '{print $1}' `

m=1000

n=64

if [[ ! -f "Matrix.txt" ]] ; then
 ln -s $indir/mat_vcf.txt Matrix.txt
fi

# Transposed version: r=n_snp

r=$n_snp
$loadmatrix  -m $m -n $n -r $r \
    -t \
	 -o ${proj}_transp.h5 \
	# -i $indir/mat_vcf.txt

# Untransposed version - nope
# r=$n_snp
#
# $loadmatrix  -m $m -n $n -r $r \
#   -o ${proj}.h5 \
#   # -i $indir/mat_vcf.txt



