#!bash
#-- data used in this analysis
snp=./bac_spring_snp.vcf2plink_geno02maf001   # prefix of genotype data (PLINK binary format which has three files including: .bed, .bim and .fam)
phe=./bac_spring_biomarker.phe
outdir=./Spring_gaws_res

# GWAS analysis for all available OTU 
##-- total OTU number
total=`wc -l ${phe}.id | cut -d ' ' -f 1`
##-- run one by one
for i in $(seq 1 $total); do
    name=$(sed -n "${i}p" "${phe}.id" | awk '{print $1}')

# Make subfolder, ensuring quotes protect paths with spaces
spt_dir="${outdir}/${name}"
mkdir -p "$spt_dir"
    
    # kinship to correct the relatedness of indiviuals
    gemma --bfile $snp -p $phe -miss 0.1 -maf 0.01 -gk -outdir $spt_dir -o snp.kinship

    # first ten PC to represent the population structure/strafation.
    plink --bfile $snp --pca 10 --out $spt_dir/PC10
    cut -d ' ' -f 3- $spt_dir/PC10.eigenvec | sed 's/^/1 /' > $spt_dir/snp_PC10

    # Genome-wide association study (GWAS) by GEMMA
    gemma --bfile $snp -p $phe -miss 0.1 -n $i -k $spt_dir/snp.kinship.cXX.txt -c $spt_dir/snp_PC10 -lmm -outdir $spt_dir -o $name
    
    # plot GWAS result
    cat > $spt_dir/$name.CMplot.r <<EOF
library(CMplot)
df <- read.table("$spt_dir/$name.assoc.for_plot.txt", header=TRUE)
CMplot(df, plot.type = "m", LOG10=TRUE, cex=0.4, threshold=2.2E-5, threshold.lty=2, threshold.col="black",  amplify=TRUE, signal.col="red", signal.cex=1,width=14, bin.size=1e6, chr.den.col=c("darkgreen", "yellow", "red"), verbose=TRUE, height=5, file="pdf")
EOF
    Rscript $spt_dir/$name.CMplot.r
done
