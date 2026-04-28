#!/usr/bin/env python3
"""
create_long_table_complete.py
"""

import gzip
import re
import pandas as pd
import argparse
import subprocess
from collections import defaultdict

def get_sample_names(vcf_file):
    """Get list of samples from VCF"""
    cmd = f"bcftools query -l {vcf_file}"
    result = subprocess.check_output(cmd, shell=True, text=True)
    return result.strip().split('\n')

def parse_pharmgkb_field(pharmgkb_value):
    """Parse of PharmGKB field"""
    if not pharmgkb_value or pharmgkb_value == '.':
        return None
    
    parts = pharmgkb_value.split('|')
    
    result = {}
    if len(parts) >= 1:
        result['PHARMGKB_Gene'] = parts[0]
    if len(parts) >= 2:
        result['PHARMGKB_Drug'] = parts[1]
    if len(parts) >= 3:
        result['PHARMGKB_Type'] = parts[2]
    if len(parts) >= 4:
        result['PHARMGKB_Effect'] = parts[3]
    if len(parts) >= 5:
        result['PHARMGKB_Level'] = parts[4]
    
    return result

def parse_clinvar_phenotypes(clndn_value, clndisdb_value):
    """Parse of phenotypes from ClinVar"""
    phenotypes = []
    if clndn_value and clndn_value != '.':
        # Split by | for multiple phenotypes
        for phen in clndn_value.split('|'):
            if phen and phen != 'not_provided':
                phenotypes.append(phen)
    
    diseases = []
    if clndisdb_value and clndisdb_value != '.':
        for dis in clndisdb_value.split('|'):
            if dis and ':' in dis:
                # Extract OMIM, MONDO, Orphanet IDs
                diseases.append(dis)
    
    return {
        'PHENOTYPES': '|'.join(phenotypes),
        'DISEASE_IDS': '|'.join(diseases)
    }

def extract_field_by_name(info, field_name):
    """Get field from INFO by name"""
    match = re.search(rf'{field_name}=([^;\t]+)', info)
    return match.group(1) if match else ''

def parse_vcf_complete(vcf_file, output_file):
    """Make full long format table with all fields"""
    print(f"Parsing {vcf_file}...")
    
    samples = get_sample_names(vcf_file)
    print(f"Samples found: {len(samples)}")
    
    # Complete CSQ indices from YOUR VCF header
    csq_fields = [
        'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 
        'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 
        'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 
        'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 
        'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE', 
        'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'TSL', 'ENSP', 'SWISSPROT', 
        'TREMBL', 'UNIPARC', 'UNIPROT_ISOFORM', 'GENE_PHENO', 'SIFT', 
        'PolyPhen', 'DOMAINS', 'HGVS_OFFSET', 'HGVSg', 'AF', 'AFR_AF', 
        'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 
        'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 
        'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_MID_AF', 
        'gnomADe_NFE_AF', 'gnomADe_REMAINING_AF', 'gnomADe_SAS_AF', 
        'gnomADg_AF', 'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 
        'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF', 
        'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF', 
        'gnomADg_REMAINING_AF', 'gnomADg_SAS_AF', 'MAX_AF', 
        'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 
        'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 
        'TRANSCRIPTION_FACTORS'
    ]
    
    data = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    with opener(vcf_file, 'rt') as f:
        for line_num, line in enumerate(f):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            pos = fields[1]
            rsid = fields[2] if fields[2] != '.' else ''
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_val = fields[6]
            info = fields[7]
            format_field = fields[8] if len(fields) > 8 else ''
            genotypes = fields[9:] if len(fields) > 9 else []
            
            # Extract CSQ - find transcript with PICK=1
            csq_match = re.search(r'CSQ=([^;\t]+)', info)
            csq_dict = {field: '' for field in csq_fields}
            
            if csq_match:
                csq_string = csq_match.group(1)
                transcripts = csq_string.split(',')
                
                # Try to find PICK=1 transcript
                pick_transcript = None
                for trans in transcripts:
                    values = trans.split('|')
                    if len(values) > 20 and values[20] and 'PICK' in values[20]:
                        pick_transcript = trans
                        break
                
                # If no PICK, use first
                if not pick_transcript:
                    pick_transcript = transcripts[0]
                
                values = pick_transcript.split('|')
                
                for i, field_name in enumerate(csq_fields):
                    if i < len(values):
                        csq_dict[field_name] = values[i]
            
            # Extract all INFO fields
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = 'True'
            
            # Extract PHARMGKB
            pharmgkb_raw = info_dict.get('PHARMGKB', '')
            pharmgkb_info = parse_pharmgkb_field(pharmgkb_raw) if pharmgkb_raw else {}
            
            # Parse ClinVar phenotypes
            clinvar_phenotypes = parse_clinvar_phenotypes(
                info_dict.get('CLNDN', ''),
                info_dict.get('CLNDISDB', '')
            )
            
            # Build complete variant record
            variant_record = {
                # Basic position info
                'CHROM': chrom,
                'POS': int(pos),
                'rsID': rsid,
                'REF': ref,
                'ALT': alt,
                'QUAL': qual,
                'FILTER': filter_val,
                
                # Quality scores (from your VCF example)
                'MCD': info_dict.get('MCD', ''),
                'GenCallScore': info_dict.get('GenCallScore', ''),
                
                # === CSQ FIELDS (all of them) ===
                'CSQ_Allele': csq_dict.get('Allele', ''),
                'Consequence': csq_dict.get('Consequence', ''),
                'IMPACT': csq_dict.get('IMPACT', ''),
                'Gene_Symbol': csq_dict.get('SYMBOL', ''),
                'Gene_ID': csq_dict.get('Gene', ''),
                'Feature_type': csq_dict.get('Feature_type', ''),
                'Feature': csq_dict.get('Feature', ''),
                'BIOTYPE': csq_dict.get('BIOTYPE', ''),
                'EXON': csq_dict.get('EXON', ''),
                'INTRON': csq_dict.get('INTRON', ''),
                'HGVSc': csq_dict.get('HGVSc', ''),
                'HGVSp': csq_dict.get('HGVSp', ''),
                'cDNA_position': csq_dict.get('cDNA_position', ''),
                'CDS_position': csq_dict.get('CDS_position', ''),
                'Protein_position': csq_dict.get('Protein_position', ''),
                'Amino_acids': csq_dict.get('Amino_acids', ''),
                'Codons': csq_dict.get('Codons', ''),
                'Existing_variation': csq_dict.get('Existing_variation', ''),
                'DISTANCE': csq_dict.get('DISTANCE', ''),
                'STRAND': csq_dict.get('STRAND', ''),
                'FLAGS': csq_dict.get('FLAGS', ''),
                'VARIANT_CLASS': csq_dict.get('VARIANT_CLASS', ''),
                'SYMBOL_SOURCE': csq_dict.get('SYMBOL_SOURCE', ''),
                'HGNC_ID': csq_dict.get('HGNC_ID', ''),
                'CANONICAL': csq_dict.get('CANONICAL', ''),
                'MANE': csq_dict.get('MANE', ''),
                'MANE_SELECT': csq_dict.get('MANE_SELECT', ''),
                'MANE_PLUS_CLINICAL': csq_dict.get('MANE_PLUS_CLINICAL', ''),
                'TSL': csq_dict.get('TSL', ''),
                'ENSP': csq_dict.get('ENSP', ''),
                'SWISSPROT': csq_dict.get('SWISSPROT', ''),
                'TREMBL': csq_dict.get('TREMBL', ''),
                'UNIPARC': csq_dict.get('UNIPARC', ''),
                'UNIPROT_ISOFORM': csq_dict.get('UNIPROT_ISOFORM', ''),
                'GENE_PHENO': csq_dict.get('GENE_PHENO', ''),
                
                # === FUNCTIONAL PREDICTIONS ===
                'SIFT': csq_dict.get('SIFT', ''),
                'PolyPhen': csq_dict.get('PolyPhen', ''),
                'DOMAINS': csq_dict.get('DOMAINS', ''),
                
                # === HGVS ===
                'HGVS_OFFSET': csq_dict.get('HGVS_OFFSET', ''),
                'HGVSg': csq_dict.get('HGVSg', ''),
                
                # === 1000 GENOMES FREQUENCIES ===
                'AF_1KG': csq_dict.get('AF', ''),
                'AFR_AF': csq_dict.get('AFR_AF', ''),
                'AMR_AF': csq_dict.get('AMR_AF', ''),
                'EAS_AF': csq_dict.get('EAS_AF', ''),
                'EUR_AF': csq_dict.get('EUR_AF', ''),
                'SAS_AF': csq_dict.get('SAS_AF', ''),
                
                # === GNOMAD EXOME FREQUENCIES ===
                'gnomADe_AF': csq_dict.get('gnomADe_AF', ''),
                'gnomADe_AFR_AF': csq_dict.get('gnomADe_AFR_AF', ''),
                'gnomADe_AMR_AF': csq_dict.get('gnomADe_AMR_AF', ''),
                'gnomADe_ASJ_AF': csq_dict.get('gnomADe_ASJ_AF', ''),
                'gnomADe_EAS_AF': csq_dict.get('gnomADe_EAS_AF', ''),
                'gnomADe_FIN_AF': csq_dict.get('gnomADe_FIN_AF', ''),
                'gnomADe_MID_AF': csq_dict.get('gnomADe_MID_AF', ''),
                'gnomADe_NFE_AF': csq_dict.get('gnomADe_NFE_AF', ''),
                'gnomADe_REMAINING_AF': csq_dict.get('gnomADe_REMAINING_AF', ''),
                'gnomADe_SAS_AF': csq_dict.get('gnomADe_SAS_AF', ''),
                
                # === GNOMAD GENOME FREQUENCIES ===
                'gnomADg_AF': csq_dict.get('gnomADg_AF', ''),
                'gnomADg_AFR_AF': csq_dict.get('gnomADg_AFR_AF', ''),
                'gnomADg_AMI_AF': csq_dict.get('gnomADg_AMI_AF', ''),
                'gnomADg_AMR_AF': csq_dict.get('gnomADg_AMR_AF', ''),
                'gnomADg_ASJ_AF': csq_dict.get('gnomADg_ASJ_AF', ''),
                'gnomADg_EAS_AF': csq_dict.get('gnomADg_EAS_AF', ''),
                'gnomADg_FIN_AF': csq_dict.get('gnomADg_FIN_AF', ''),
                'gnomADg_MID_AF': csq_dict.get('gnomADg_MID_AF', ''),
                'gnomADg_NFE_AF': csq_dict.get('gnomADg_NFE_AF', ''),
                'gnomADg_REMAINING_AF': csq_dict.get('gnomADg_REMAINING_AF', ''),
                'gnomADg_SAS_AF': csq_dict.get('gnomADg_SAS_AF', ''),
                
                # === MAX FREQUENCIES ===
                'MAX_AF': csq_dict.get('MAX_AF', ''),
                'MAX_AF_POPS': csq_dict.get('MAX_AF_POPS', ''),
                
                # === CLINVAR FIELDS (from CSQ) ===
                'CLIN_SIG_CSQ': csq_dict.get('CLIN_SIG', ''),
                'SOMATIC_CSQ': csq_dict.get('SOMATIC', ''),
                'PHENO_CSQ': csq_dict.get('PHENO', ''),
                
                # === PUBMED ===
                'PUBMED': csq_dict.get('PUBMED', ''),
                
                # === REGULATORY FEATURES ===
                'MOTIF_NAME': csq_dict.get('MOTIF_NAME', ''),
                'MOTIF_POS': csq_dict.get('MOTIF_POS', ''),
                'HIGH_INF_POS': csq_dict.get('HIGH_INF_POS', ''),
                'MOTIF_SCORE_CHANGE': csq_dict.get('MOTIF_SCORE_CHANGE', ''),
                'TRANSCRIPTION_FACTORS': csq_dict.get('TRANSCRIPTION_FACTORS', ''),
                
                # === CLINVAR INFO FIELDS (from bcftools annotate) ===
                'CLNSIG': info_dict.get('CLNSIG', ''),
                'CLNREVSTAT': info_dict.get('CLNREVSTAT', ''),
                'CLNDN': info_dict.get('CLNDN', ''),
                'CLNDISDB': info_dict.get('CLNDISDB', ''),
                'CLNHGVS': info_dict.get('CLNHGVS', ''),
                'CLNVC': info_dict.get('CLNVC', ''),
                'GENEINFO': info_dict.get('GENEINFO', ''),
                'MC': info_dict.get('MC', ''),
                'ORIGIN': info_dict.get('ORIGIN', ''),
                'CLNSIGCONF': info_dict.get('CLNSIGCONF', ''),
                
                # === PARSED CLINVAR PHENOTYPES ===
                'ClinVar_Phenotypes': clinvar_phenotypes.get('PHENOTYPES', ''),
                'ClinVar_Disease_IDS': clinvar_phenotypes.get('DISEASE_IDS', ''),
                
                # === PHARMGKB ===
                'PHARMGKB_Raw': pharmgkb_raw,
                'PHARMGKB_Gene': pharmgkb_info.get('PHARMGKB_Gene', ''),
                'PHARMGKB_Drug': pharmgkb_info.get('PHARMGKB_Drug', ''),
                'PHARMGKB_Type': pharmgkb_info.get('PHARMGKB_Type', ''),
                'PHARMGKB_Effect': pharmgkb_info.get('PHARMGKB_Effect', ''),
                'PHARMGKB_Level': pharmgkb_info.get('PHARMGKB_Level', ''),
                'Is_PharmCAT': 'Yes' if pharmgkb_info.get('PHARMGKB_Level') == 'PharmCAT' else 'No',
                'Is_High_Evidence': 'Yes' if pharmgkb_info.get('PHARMGKB_Level', '').startswith(('1', '2')) else 'No',
            }
            
            # Process each sample
            for i, sample in enumerate(samples):
                if i >= len(genotypes):
                    continue
                
                gt = genotypes[i]
                
                # Check if sample has alternate allele
                if gt in ['0/1', '1/0', '1/1', '0|1', '1|0', '1|1']:
                    row = variant_record.copy()
                    row['Sample_ID'] = sample
                    row['Genotype'] = gt
                    row['Genotype_Type'] = 'HET' if gt in ['0/1', '1/0', '0|1', '1|0'] else 'HOM_ALT'
                    
                    # Extract genotype fields if FORMAT present
                    if format_field and len(genotypes) > i:
                        gt_fields = dict(zip(format_field.split(':'), genotypes[i].split(':')))
                        for key, value in gt_fields.items():
                            if key not in ['GT']:  # GT already stored
                                row[f'GT_{key}'] = value
                    
                    data.append(row)
            
            if len(data) % 50000 == 0 and len(data) > 0:
                print(f"  Collected {len(data):,} records...")
    
    df = pd.DataFrame(data)
    print(f"\n Total records: {len(df):,}")
    print(f"   Unique variants: {df.groupby(['CHROM','POS','REF','ALT']).ngroups:,}")
    print(f"   Unique samples: {df['Sample_ID'].nunique():,}")
    
    if len(df) > 0:
        # Statistics
        print(f"\n Statistics:")
        
        # Pathogenic variants
        if 'CLNSIG' in df.columns:
            pathogenic = df[df['CLNSIG'].str.contains('Pathogenic', na=False) & ~df['CLNSIG'].str.contains('Benign', na=False)]
            print(f"   Pathogenic/LP variants: {pathogenic['CHROM'].nunique():,}")
        
        # PHARMGKB
        if 'PHARMGKB_Level' in df.columns:
            pharmgkb_records = df[df['PHARMGKB_Level'] != ''].shape[0]
            print(f"   PHARMGKB annotated records: {pharmgkb_records:,}")
        
        # Functional predictions
        if 'SIFT' in df.columns:
            sift_count = df[df['SIFT'] != ''].shape[0]
            print(f"   SIFT predictions: {sift_count:,}")
        
        if 'PolyPhen' in df.columns:
            polyphen_count = df[df['PolyPhen'] != ''].shape[0]
            print(f"   PolyPhen predictions: {polyphen_count:,}")
        
        # Frequencies
        if 'gnomADe_AF' in df.columns:
            gnomade_count = df[df['gnomADe_AF'] != ''].shape[0]
            print(f"   gnomAD exome frequencies: {gnomade_count:,}")
        
        if 'gnomADg_AF' in df.columns:
            gnomadg_count = df[df['gnomADg_AF'] != ''].shape[0]
            print(f"   gnomAD genome frequencies: {gnomadg_count:,}")
        
        # Phenotypes
        if 'ClinVar_Phenotypes' in df.columns:
            pheno_count = df[df['ClinVar_Phenotypes'] != ''].shape[0]
            print(f"   ClinVar phenotypes: {pheno_count:,}")
        
        # Available columns
        print(f"\n Total columns in output: {len(df.columns)}")
    
    # Save
    df.to_csv(output_file, index=False, sep='\t')
    print(f"\n Table saved to {output_file}")
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Create complete long format variant table with ALL fields')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()
    
    df = parse_vcf_complete(args.input, args.output)
    
    # Display sample of important columns
    print("\n First 5 rows (key columns):")
    display_cols = ['Sample_ID', 'Genotype', 'CHROM', 'POS', 'Gene_Symbol', 
                   'Consequence', 'IMPACT', 'HGVSc', 'HGVSp', 'SIFT', 'PolyPhen',
                   'CLNSIG', 'CLNREVSTAT', 'ClinVar_Phenotypes',
                   'PHARMGKB_Level', 'PHARMGKB_Drug', 
                   'AF_1KG', 'gnomADe_AF', 'gnomADg_AF', 'MAX_AF', 'MAX_AF_POPS']
    
    # Only show columns that exist
    display_cols = [col for col in display_cols if col in df.columns]
    if display_cols:
        print(df[display_cols].head(10).to_string())
    else:
        print(df.head().to_string())

if __name__ == "__main__":
    main()
