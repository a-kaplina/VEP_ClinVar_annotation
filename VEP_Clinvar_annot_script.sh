#!/bin/bash
set -euo pipefail

# Check input argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <absolute/path/to/source_folder>"
    echo "Example: $0 /home/user/data"
    exit 1
fi

SOURCE_DIR="$1"

# Check if source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
    echo "ERROR: Source directory $SOURCE_DIR does not exist"
    exit 1
fi

# Define paths
REF=~/reference/hg38/hg38.fa
CLINVAR=~/clinvar38/clinvar.vcf.gz
VEPCACHEDIR=~/.vep/
VEPPLUGINS=~/.vep/Plugins
REF_VEP=~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create directory structure
INPUT_DIR="${SOURCE_DIR}/input"
OUTPUT_DIR="${SOURCE_DIR}/output"
LOG_FILE="${OUTPUT_DIR}/annotation_log.txt"

mkdir -p "$INPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Initialize log file
echo "=========================================" > "$LOG_FILE"
echo "VCF Annotation Pipeline Log" >> "$LOG_FILE"
echo "Start time: $(date)" >> "$LOG_FILE"
echo "Source directory: $SOURCE_DIR" >> "$LOG_FILE"
echo "=========================================" >> "$LOG_FILE"

# Function to log and print
log() {
    echo "$1" | tee -a "$LOG_FILE"
}

# Function to count variants in VCF
count_variants() {
    local vcf_file=$1
    if [ -f "$vcf_file" ]; then
        bcftools view -H "$vcf_file" 2>/dev/null | wc -l
    else
        echo "0"
    fi
}

log ""
log "========================================="
log "Source directory: $SOURCE_DIR"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"
log "========================================="

# Function for VEP annotation
anno() {
    local invcf=$1
    local outvcf=$2
    local tmpvcf=${outvcf%%.vcf.gz}.tmp.vcf.gz
    local threads=8
    
    log "  - Running VEP..."
    ~/ensembl-vep/vep --assembly GRCh38 --fork "$threads" \
        -i "$invcf" -o "$tmpvcf" --vcf --force_overwrite \
        --cache --offline --cache_version 115 --compress_output bgzip \
        --pick \
        --gencode_primary \
        --fasta "$REF_VEP" \
        --dir_plugins "$VEPPLUGINS" \
        --sift b --polyphen b \
        --af --af_1kg --max_af \
        --af_gnomade --af_gnomadg \
        --hgvs --hgvsg \
        --protein --canonical --mane \
        --domains \
        --pubmed --regulatory --check_existing \
        --tsl --uniprot \
        --gene_phenotype \
        --variant_class \
        --shift_hgvs 1 \
        --total_length
    
    if [ ! -f "$tmpvcf" ]; then
        log "ERROR: VEP failed for $invcf"
        return 1
    fi
    
    log "  - Indexing VEP output..."
    bcftools index "$tmpvcf"
    
    log "  - Adding ClinVar annotations..."
    bcftools annotate -a "$CLINVAR" \
        -c CLNDISDB,CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,CLNVC,GENEINFO,MC,ORIGIN,CLNSIGCONF \
        -Oz -o "$outvcf" "$tmpvcf"
    
    log "  - Cleaning temporary files..."
    rm "$tmpvcf" "$tmpvcf.csi" 2>/dev/null || true
    bcftools index -f "$outvcf"
}

# Find all VCF files in source directory
log "Looking for VCF files in $SOURCE_DIR..."
find "$SOURCE_DIR" -maxdepth 1 -name "*.vcf.gz" -type f | while read -r vcf; do
    BASENAME=$(basename "$vcf" .vcf.gz)
    log ""
    log "========================================="
    log "Processing: $BASENAME"
    log "========================================="
    
    # Define file paths for this sample
    SORTED_VCF="${INPUT_DIR}/${BASENAME}_sorted.vcf.gz"
    NORM_VCF="${INPUT_DIR}/${BASENAME}_norm.vcf.gz"
    NOCHR_VCF="${INPUT_DIR}/${BASENAME}_nochr.vcf.gz"
    FINAL_VCF="${OUTPUT_DIR}/${BASENAME}_annotated.vcf.gz"
    PHARMGKB_VCF="${OUTPUT_DIR}/${BASENAME}_annotated_Pharm.vcf.gz"
    REAL_VCF="${OUTPUT_DIR}/${BASENAME}_real_annotated_Pharm.vcf.gz"
    PATHOGENIC_VCF="${OUTPUT_DIR}/${BASENAME}_pathogenic_only.vcf.gz"
    GOOD_Q_VCF="${OUTPUT_DIR}/${BASENAME}_good_q_pathogenic_only.vcf.gz"
    VUS_VCF="${OUTPUT_DIR}/${BASENAME}_vus_only.vcf.gz"
    GOOD_Q_VUS_VCF="${OUTPUT_DIR}/${BASENAME}_good_q_vus_only.vcf.gz"
    
    # Track initial variant counts
    INITIAL_COUNT=$(count_variants "$vcf")
    log "Initial variant count in $BASENAME: $INITIAL_COUNT"
    
    # Step 1: Sort the VCF
    log "Step 1: Sorting VCF..."
    bcftools sort "$vcf" -o "$SORTED_VCF"
    bcftools index "$SORTED_VCF"
    log "  - Created: $SORTED_VCF ($(count_variants "$SORTED_VCF") variants)"
    
    # Step 2: Normalize VCF (split multiallelic sites)
    log "Step 2: Normalizing VCF..."
    bcftools norm -m -any -f "$REF" "$SORTED_VCF" -Oz -o "$NORM_VCF"
    bcftools index "$NORM_VCF"
    log "  - Created: $NORM_VCF ($(count_variants "$NORM_VCF") variants)"
    
    # Step 3: Create chromosome rename mapping
    log "Step 3: Preparing chromosome rename mapping..."
    cat > "${INPUT_DIR}/chr_rename.txt" << 'EOF'
chr1 1
chr2 2
chr3 3
chr4 4
chr5 5
chr6 6
chr7 7
chr8 8
chr9 9
chr10 10
chr11 11
chr12 12
chr13 13
chr14 14
chr15 15
chr16 16
chr17 17
chr18 18
chr19 19
chr20 20
chr21 21
chr22 22
chrX X
chrY Y
chrMT MT
EOF
    
    # Step 4: Remove 'chr' prefix from chromosome names
    log "Step 4: Removing 'chr' prefix from chromosome names..."
    bcftools annotate --rename-chrs "${INPUT_DIR}/chr_rename.txt" \
        -Oz -o "$NOCHR_VCF" "$NORM_VCF"
    bcftools index "$NOCHR_VCF"
    log "  - Created: $NOCHR_VCF ($(count_variants "$NOCHR_VCF") variants)"
    
    # Step 5: Remove intermediate files (sorted and normalized)
    log "Step 5: Removing intermediate files..."
    rm -f "$SORTED_VCF" "$SORTED_VCF.csi"
    rm -f "$NORM_VCF" "$NORM_VCF.csi"
    rm -f "${INPUT_DIR}/chr_rename.txt"
    log "  - Removed intermediate files"
    
    # Step 6: Run VEP and ClinVar annotation
    log "Step 6: Running VEP and ClinVar annotation..."
    anno "$NOCHR_VCF" "$FINAL_VCF"
    log "  - Created: $FINAL_VCF ($(count_variants "$FINAL_VCF") variants)"
    
    # Step 7: Remove NOCHR_VCF (no longer needed)
    log "Step 7: Removing temporary nochr VCF..."
    rm -f "$NOCHR_VCF" "$NOCHR_VCF.csi"
    log "  - Removed: $NOCHR_VCF"
    
    # Step 8: Add PharmGKB annotations
    log "Step 8: Adding PharmGKB annotations..."
    if [ -f "${SCRIPT_DIR}/pharmgkb_map_clean_CAT_all.txt" ]; then
        log "  - Found PharmGKB mapping file"
        
        # Create PharmGKB annotation map
        awk 'BEGIN {FS="\t"; OFS="\t"} 
             NR==FNR {id_info[$1]=$2; next} 
             $3 in id_info {print $1, $2, $4, $5, id_info[$3]}' \
             "${SCRIPT_DIR}/pharmgkb_map_clean_CAT_all.txt" \
             <(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$FINAL_VCF") \
             > "${OUTPUT_DIR}/${BASENAME}_pharmgkb_with_ref_alt.txt"
        
        PHARMGKB_COUNT=$(wc -l < "${OUTPUT_DIR}/${BASENAME}_pharmgkb_with_ref_alt.txt")
        log "  - Found $PHARMGKB_COUNT variants with PharmGKB annotations"
        
        # Add PharmGKB INFO field to VCF
        zcat "$FINAL_VCF" | \
        awk 'BEGIN {FS="\t"; OFS="\t"} 
             NR==FNR {key=$1":"$2":"$3":"$4; ann[key]=$5; next}
             {
                 if ($0 ~ /^#/) print
                 else {
                     key = $1":"$2":"$4":"$5
                     if (key in ann) {
                         if ($8 == ".") $8 = ann[key]
                         else $8 = $8 ";" ann[key]
                     }
                     print
                 }
             }' "${OUTPUT_DIR}/${BASENAME}_pharmgkb_with_ref_alt.txt" - | \
        sed '/^##fileformat=/a ##INFO=<ID=PHARMGKB,Number=1,Type=String,Description="PharmGKB clinical annotation (evidence 1-4 and PharmCAT)">' | \
        bgzip -c > "$PHARMGKB_VCF"
        
        bcftools index "$PHARMGKB_VCF"
        log "  - Created: $PHARMGKB_VCF ($(count_variants "$PHARMGKB_VCF") variants)"
        
        # Remove the temporary mapping file
        rm -f "${OUTPUT_DIR}/${BASENAME}_pharmgkb_with_ref_alt.txt"
        
        # Count PHARMGKB annotation levels
        PHARMGKB_TOTAL=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "PHARMGKB=" || echo "0")
        PHARMGKB_1A=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|1A" || echo "0")
        PHARMGKB_1B=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|1B" || echo "0")
        PHARMGKB_2A=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|2A" || echo "0")
        PHARMGKB_2B=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|2B" || echo "0")
        PHARMGKB_3=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|3" || echo "0")
        PHARMGKB_4=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|4" || echo "0")
        PHARMGKB_CAT=$(bcftools view "$PHARMGKB_VCF" 2>/dev/null | grep -c "|PharmCAT" || echo "0")
        
        log "    - PHARMGKB total: $PHARMGKB_TOTAL"
        log "    - Level 1A: $PHARMGKB_1A"
        log "    - Level 1B: $PHARMGKB_1B"
        log "    - Level 2A: $PHARMGKB_2A"
        log "    - Level 2B: $PHARMGKB_2B"
        log "    - Level 3: $PHARMGKB_3"
        log "    - Level 4: $PHARMGKB_4"
        log "    - PharmCAT: $PHARMGKB_CAT"
    else
        log "  - WARNING: PharmGKB mapping file not found in $SCRIPT_DIR"
        log "  - Skipping PharmGKB annotation"
        cp "$FINAL_VCF" "$PHARMGKB_VCF"
        bcftools index "$PHARMGKB_VCF"
    fi
    
    # Step 9: Remove FINAL_VCF (keep only PharmGKB version)
    log "Step 9: Removing intermediate annotated VCF..."
    rm -f "$FINAL_VCF" "$FINAL_VCF.csi"
    log "  - Removed: $FINAL_VCF"
    
    # Step 10: Filter for alternate alleles only (real variants)
    log "Step 10: Filtering for alternate alleles..."
    bcftools view -i 'GT="alt"' "$PHARMGKB_VCF" | bgzip -c > "$REAL_VCF"
    bcftools index "$REAL_VCF"
    REAL_COUNT=$(count_variants "$REAL_VCF")
    log "  - Created: $REAL_VCF ($REAL_COUNT variants with alternate alleles)"
    
    # Generate pathogenicity summaries
    log "Step 11: Generating pathogenicity summaries..."
    
    # Real variants pathogenicity table
    echo -e "COUNT\tCLNSIG\tCLNREVSTAT" > "${OUTPUT_DIR}/${BASENAME}_nonfilt_annot_clnsig_clnrevstat.txt"
    bcftools query -f '%CLNSIG\t%CLNREVSTAT\n' "$PHARMGKB_VCF" 2>/dev/null | \
        grep -v '^\.' | sort | uniq -c | sort -nr | \
        awk '{print $1"\t"$2"\t"$3}' >> "${OUTPUT_DIR}/${BASENAME}_nonfilt_annot_clnsig_clnrevstat.txt"
        
    # Real variants pathogenicity table
    echo -e "COUNT\tCLNSIG\tCLNREVSTAT" > "${OUTPUT_DIR}/${BASENAME}_real_clnsig_clnrevstat.txt"
    bcftools query -f '%CLNSIG\t%CLNREVSTAT\n' "$REAL_VCF" 2>/dev/null | \
        grep -v '^\.' | sort | uniq -c | sort -nr | \
        awk '{print $1"\t"$2"\t"$3}' >> "${OUTPUT_DIR}/${BASENAME}_real_clnsig_clnrevstat.txt"
    
    # Step 12: Extract pathogenic only
    log "Step 12: Extracting pathogenic/likely pathogenic variants..."
    bcftools view -i 'CLNSIG="Pathogenic" || CLNSIG="Likely_pathogenic" || CLNSIG="Pathogenic/Likely_pathogenic" || CLNSIG="Pathogenic/Pathogenic,_low_penetrance|other" || CLNSIG="Pathogenic|other" || CLNSIG="Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance/Established_risk_allele|risk_factor"' \
        "$REAL_VCF" -Oz -o "$PATHOGENIC_VCF"
    bcftools index "$PATHOGENIC_VCF"
    PATHO_COUNT=$(count_variants "$PATHOGENIC_VCF")
    log "  - Created: $PATHOGENIC_VCF ($PATHO_COUNT pathogenic/likely pathogenic variants)"
    
    # Step 13: Filter for good quality review status (Pathogenic)
    log "Step 13: Filtering for good quality review status (Pathogenic)..."
    bcftools view -i 'CLNREVSTAT="criteria_provided,_multiple_submitters,_no_conflicts" || CLNREVSTAT="reviewed_by_expert_panel" || CLNREVSTAT="practice_guideline" || CLNREVSTAT="criteria_provided,_single_submitter"' \
        "$PATHOGENIC_VCF" -Oz -o "$GOOD_Q_VCF"
    bcftools index "$GOOD_Q_VCF"
    GOOD_Q_COUNT=$(count_variants "$GOOD_Q_VCF")
    log "  - Created: $GOOD_Q_VCF ($GOOD_Q_COUNT high-quality pathogenic variants)"
    
    # Step 14: Extract VUS (Uncertain Significance) variants
    log "Step 14: Extracting VUS (Uncertain Significance) variants..."
    bcftools view -i 'CLNSIG="Uncertain_significance" || CLNSIG="Uncertain_significance|no_interpretation_for_the_single_submitter" || CLNSIG="Uncertain_significance|other"' \
        "$REAL_VCF" -Oz -o "$VUS_VCF"
    bcftools index "$VUS_VCF"
    VUS_COUNT=$(count_variants "$VUS_VCF")
    log "  - Created: $VUS_VCF ($VUS_COUNT VUS variants)"
    
    # Step 15: Filter VUS for good quality review status (2+ stars or single submitter with criteria)
    log "Step 15: Filtering VUS for good quality review status..."
    bcftools view -i 'CLNREVSTAT="criteria_provided,_multiple_submitters,_no_conflicts" || CLNREVSTAT="reviewed_by_expert_panel" || CLNREVSTAT="practice_guideline" || CLNREVSTAT="criteria_provided,_single_submitter"' \
        "$VUS_VCF" -Oz -o "$GOOD_Q_VUS_VCF"
    bcftools index "$GOOD_Q_VUS_VCF"
    GOOD_Q_VUS_COUNT=$(count_variants "$GOOD_Q_VUS_VCF")
    log "  - Created: $GOOD_Q_VUS_VCF ($GOOD_Q_VUS_COUNT high-quality VUS variants)"
    
    # Generate VUS pathogenicity summary
    echo -e "COUNT\tCLNSIG\tCLNREVSTAT" > "${OUTPUT_DIR}/${BASENAME}_vus_clnsig_clnrevstat.txt"
    bcftools query -f '%CLNSIG\t%CLNREVSTAT\n' "$GOOD_Q_VUS_VCF" 2>/dev/null | \
        grep -v '^\.' | sort | uniq -c | sort -nr | \
        awk '{print $1"\t"$2"\t"$3}' >> "${OUTPUT_DIR}/${BASENAME}_vus_clnsig_clnrevstat.txt"
    
    # Step 16: Run Python script for Pathogenic variants
    log "Step 16: Running Python annotation script for Pathogenic variants..."
    if [ -f "${SCRIPT_DIR}/make_long_table_script.py" ]; then
        python3 "${SCRIPT_DIR}/make_long_table_script.py" -i "$GOOD_Q_VCF" -o "${OUTPUT_DIR}/${BASENAME}_pathogenic_complete_annotation.tsv"
        log "  - Created: ${OUTPUT_DIR}/${BASENAME}_pathogenic_complete_annotation.tsv"
    else
        log "  - WARNING: make_long_table_script.py not found in $SCRIPT_DIR"
    fi
    
    # Step 17: Run Python script for VUS variants
    log "Step 17: Running Python annotation script for VUS variants..."
    if [ -f "${SCRIPT_DIR}/make_long_table_script.py" ]; then
        python3 "${SCRIPT_DIR}/make_long_table_script.py" -i "$GOOD_Q_VUS_VCF" -o "${OUTPUT_DIR}/${BASENAME}_vus_complete_annotation.tsv"
        log "  - Created: ${OUTPUT_DIR}/${BASENAME}_vus_complete_annotation.tsv"
    else
        log "  - WARNING: make_long_table_script.py not found in $SCRIPT_DIR"
    fi
    
    # Summary for this sample
    log ""
    log "=== SUMMARY FOR $BASENAME ==="
    log "Initial variants: $INITIAL_COUNT"
    log ""
    log "Pathogenic/LP variants:"
    log "  - All pathogenic: $PATHO_COUNT"
    log "  - High-quality pathogenic: $GOOD_Q_COUNT"
    log ""
    log "VUS (Uncertain Significance) variants:"
    log "  - All VUS: $VUS_COUNT"
    log "  - High-quality VUS: $GOOD_Q_VUS_COUNT"
    log ""
    log "Output files saved:"
    log "  - ${PHARMGKB_VCF} (all variants with PharmGKB)"
    log "  - ${PATHOGENIC_VCF} (all pathogenic/likely pathogenic)"
    log "  - ${GOOD_Q_VCF} (high-quality pathogenic)"
    log "  - ${VUS_VCF} (all VUS)"
    log "  - ${GOOD_Q_VUS_VCF} (high-quality VUS)"
    log "  - ${OUTPUT_DIR}/${BASENAME}_pathogenic_complete_annotation.tsv (pathogenic long table)"
    log "  - ${OUTPUT_DIR}/${BASENAME}_vus_complete_annotation.tsv (VUS long table)"
    log "================================"
    
done

log ""
log "========================================="
log "All files processed successfully!"
log ""
log "Final output files saved in: $OUTPUT_DIR"
log ""
log "For each sample, the following files are created:"
log "  [1] *_annotated_Pharm.vcf.gz - all variants with PharmGKB"
log "  [2] *_pathogenic_only.vcf.gz - all pathogenic/likely pathogenic variants"
log "  [3] *_good_q_pathogenic_only.vcf.gz - high-quality pathogenic variants"
log "  [4] *_vus_only.vcf.gz - all VUS (Uncertain Significance) variants"
log "  [5] *_good_q_vus_only.vcf.gz - high-quality VUS variants"
log "  [6] *_pathogenic_complete_annotation.tsv - long format table for pathogenic"
log "  [7] *_vus_complete_annotation.tsv - long format table for VUS"
log "  [8] *_real_clnsig_clnrevstat.txt - pathogenicity summary"
log "  [9] *_vus_clnsig_clnrevstat.txt - VUS summary"
log ""
log "Log file: $LOG_FILE"
log "End time: $(date)"
log "========================================="
