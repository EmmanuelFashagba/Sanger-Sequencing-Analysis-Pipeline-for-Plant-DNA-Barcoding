#!/usr/bin/env python3
"""
Sanger consensus + alignment + variant summary pipeline.
Version: 1.3 - Enhanced with strict sequence cleaning and tree visualization
"""

import os, re, sys, argparse, subprocess, textwrap, shutil, statistics
from collections import defaultdict, Counter
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Try to import visualization libraries
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    
try:
    from ete3 import Tree, TreeStyle
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False

# ---------- Configuration / Regex ----------

FILENAME_REGEX = re.compile(
    r"""^(?P<sample>.+?)_+(?P<gene>(?:RBCL|RBCLA|HSP70|WRKY5|SUT1|SUT|RBC[LA]?|MEHSP70)[A-Za-z0-9]*)_+(?P<dir>[FRfr])(?:[^/]*)\.ab1$""",
    re.IGNORECASE
)

# Quality trimming parameters
QUALITY_THRESHOLD = 20      # Minimum quality score
WINDOW_SIZE = 10           # Sliding window size
MIN_WINDOW_QUALITY = 15    # Minimum average quality in window
MAX_N_PERCENT = 5          # Maximum % of Ns allowed in window
MIN_FINAL_LENGTH = 100     # Minimum length after trimming
MAX_N_STRETCH = 5          # Maximum consecutive Ns allowed

# ---------- Utility Functions ----------

def which(cmd):
    return shutil.which(cmd) is not None

def log(msg):
    print(f"[INFO] {msg}")

def warn(msg):
    print(f"[WARN] {msg}")

def fail(msg):
    print(f"[ERROR] {msg}")
    sys.exit(1)

def parse_args():
    p = argparse.ArgumentParser(
        description="Batch process Sanger AB1 files → consensus → alignment → SNP summary → phylogenetic tree.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Example usage:
          %(prog)s --input ./ab1_files --out ./results
          %(prog)s --input ./ab1_files --out ./results --quality_threshold 25
          %(prog)s --input ./ab1_files --out ./results --aggressive_trim --strict_clean
          
        File naming convention:
          Forward: SampleName_GENE_F.ab1
          Reverse: SampleName_GENE_R.ab1
        """)
    )
    p.add_argument("--input", "-i", required=True, help="Input folder containing .ab1 traces")
    p.add_argument("--out", "-o", required=True, help="Output directory")
    p.add_argument("--genes", help="Comma-separated list of gene tokens to include (optional)")
    p.add_argument("--min_overlap", type=int, default=80, help="Min overlap (bp) for F/R merge (default: 80)")
    p.add_argument("--quality_threshold", type=int, default=20, help="Quality score threshold (default: 20)")
    p.add_argument("--window_size", type=int, default=10, help="Sliding window size for quality trimming (default: 10)")
    p.add_argument("--aggressive_trim", action="store_true", help="More aggressive N removal")
    p.add_argument("--strict_clean", action="store_true", help="Convert all non-ATCG to N (recommended)")
    p.add_argument("--no_align", action="store_true", help="Skip MAFFT & downstream steps")
    p.add_argument("--no_tree", action="store_true", help="Skip IQ-TREE even if installed")
    p.add_argument("--no_viz", action="store_true", help="Skip tree visualization")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    return p.parse_args()

def clean_sequence(seq, strict=True):
    """
    Clean sequence by handling ambiguous bases.
    
    If strict=True: Convert all ambiguous bases to N (recommended)
    If strict=False: Keep standard IUPAC codes
    """
    seq = seq.upper()
    
    if strict:
        # Only keep ATCG, convert everything else to N
        cleaned = ''.join([base if base in 'ATCG' else 'N' for base in seq])
    else:
        # Keep standard IUPAC codes
        valid_bases = 'ATCGRYSWKMBDHVN-'
        cleaned = ''.join([base if base in valid_bases else 'N' for base in seq])
    
    return cleaned

def read_ab1_with_quality(path):
    """
    Read AB1 file and return sequence with quality scores.
    Returns: (sequence, quality_scores)
    """
    try:
        record = SeqIO.read(path, "abi")
        seq = str(record.seq).upper()
        
        # Get quality scores - try different possible locations
        quality = None
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            quality = record.letter_annotations['phred_quality']
        elif hasattr(record, 'annotations'):
            # Look for quality data in annotations
            for key in ['abif_raw', 'peak_qualities', 'qual']:
                if key in record.annotations:
                    quality = record.annotations[key]
                    break
        
        # If no quality found, create fake quality scores
        if quality is None:
            warn(f"No quality scores found in {path}, using default")
            quality = [20] * len(seq)  # Default medium quality
            
        return seq, quality
    except Exception as e:
        warn(f"Failed to read {path}: {e}")
        return "", []

def quality_trim_sequence(seq, quality, args):
    """
    Advanced quality-based trimming using sliding window approach.
    """
    if not seq or not quality:
        return ""
    
    # First clean the sequence if strict mode
    if args.strict_clean:
        seq = clean_sequence(seq, strict=True)
    
    seq_array = np.array(list(seq))
    qual_array = np.array(quality)
    
    # Step 1: Find high-quality region using sliding window
    window = args.window_size
    threshold = args.quality_threshold
    
    # Calculate average quality for each window
    if len(seq) < window:
        return seq if np.mean(quality) >= threshold else ""
    
    # Find start position
    start_pos = 0
    for i in range(len(seq) - window + 1):
        window_qual = np.mean(qual_array[i:i+window])
        window_seq = seq[i:i+window]
        n_percent = (window_seq.count('N') / window) * 100
        
        if window_qual >= threshold and n_percent <= MAX_N_PERCENT:
            start_pos = i
            break
    
    # Find end position
    end_pos = len(seq)
    for i in range(len(seq) - window, -1, -1):
        window_qual = np.mean(qual_array[i:i+window])
        window_seq = seq[i:i+window]
        n_percent = (window_seq.count('N') / window) * 100
        
        if window_qual >= threshold and n_percent <= MAX_N_PERCENT:
            end_pos = i + window
            break
    
    # Trim sequence
    if start_pos < end_pos:
        trimmed_seq = seq[start_pos:end_pos]
    else:
        trimmed_seq = ""
    
    # Step 2: Additional N removal if aggressive trimming
    if args.aggressive_trim and trimmed_seq:
        trimmed_seq = remove_n_stretches(trimmed_seq)
    
    # Step 3: Final quality check
    if len(trimmed_seq) < MIN_FINAL_LENGTH:
        if args.verbose:
            warn(f"Sequence too short after trimming ({len(trimmed_seq)} bp)")
        return ""
    
    # Final cleaning to ensure no ambiguous bases
    if args.strict_clean:
        trimmed_seq = clean_sequence(trimmed_seq, strict=True)
    
    return trimmed_seq

def remove_n_stretches(seq, max_n_stretch=MAX_N_STRETCH):
    """
    Remove stretches of Ns longer than max_n_stretch.
    """
    import re
    
    # Replace runs of N longer than max_n_stretch with a single N
    pattern = f"N{{{max_n_stretch+1},}}"
    seq = re.sub(pattern, 'N' * max_n_stretch, seq)
    
    # Remove leading and trailing Ns
    seq = seq.strip('N')
    
    # If too many Ns remain, try to salvage good regions
    if seq and seq.count('N') / len(seq) > 0.1:  # More than 10% Ns
        # Split by N-stretches and keep longest good fragment
        fragments = [f for f in seq.split('N' * max_n_stretch) if len(f) > 50]
        if fragments:
            seq = max(fragments, key=len)
    
    return seq

def pairwise_overlap_quality(fwd, rev_rc, fwd_qual, rev_qual, min_overlap=80):
    """
    Enhanced overlap detection considering quality scores.
    """
    best = (0, 0, 0.0, 0.0)  # offset, overlap, identity, avg_quality
    max_shift = min(len(fwd), len(rev_rc))
    
    for offset in range(-max_shift + 10, max_shift - 10, 5):
        matches = 0
        total = 0
        qual_sum = 0
        
        for i, b in enumerate(fwd):
            j = i - offset
            if 0 <= j < len(rev_rc):
                total += 1
                # Consider quality when comparing
                if fwd[i] == rev_rc[j] and fwd[i] in "ACGT":
                    matches += 1
                    # Average quality at this position
                    if fwd_qual and rev_qual:
                        qual_sum += (fwd_qual[i] + rev_qual[j]) / 2
                
        if total >= min_overlap and total > 0:
            ident = matches / total
            avg_qual = qual_sum / total if matches > 0 else 0
            
            # Prefer high identity AND high quality overlaps
            score = ident * 0.7 + (avg_qual / 40) * 0.3  # Combined score
            if score > (best[2] * 0.7 + (best[3] / 40) * 0.3):
                best = (offset, total, ident, avg_qual)
                
    return best

def build_consensus_quality(fwd, rev_rc, fwd_qual, rev_qual, offset, strict_clean=True):
    """
    Build consensus using quality scores to resolve conflicts.
    """
    start = min(0, -offset)
    end = max(len(fwd), len(rev_rc) + offset)
    cons = []
    
    for pos in range(start, end):
        i = pos
        j = pos + offset
        
        # Get bases and qualities
        b1 = fwd[i] if 0 <= i < len(fwd) else "-"
        b2 = rev_rc[j] if 0 <= j < len(rev_rc) else "-"
        q1 = fwd_qual[i] if (fwd_qual and 0 <= i < len(fwd)) else 0
        q2 = rev_qual[j] if (rev_qual and 0 <= j < len(rev_rc)) else 0
        
        # Resolve base
        if b1 == b2 and b1 != "-":
            base = b1
        elif b1 == "-":
            base = b2
        elif b2 == "-":
            base = b1
        elif b1 in "ACGT" and b2 in "ACGT":
            # Use higher quality base
            base = b1 if q1 >= q2 else b2
        elif b1 in "ACGT":
            base = b1
        elif b2 in "ACGT":
            base = b2
        else:
            base = "N"
            
        if base and base != "-":
            # Clean ambiguous bases if strict mode
            if strict_clean and base not in "ACGTN":
                base = "N"
            cons.append(base)
            
    return "".join(cons)

def create_tree_visualization(treefile, output_dir, gene_name=""):
    """Create multiple tree visualizations"""
    visualizations_created = []
    
    # 1. Try matplotlib/Biopython visualization
    if MATPLOTLIB_AVAILABLE:
        try:
            tree = Phylo.read(treefile, "newick")
            
            # Create figure with better settings
            fig = plt.figure(figsize=(14, 10))
            ax = fig.add_subplot(1, 1, 1)
            
            # Draw tree with customization
            Phylo.draw(tree, ax=ax, 
                      show_confidence=True,
                      branch_labels=lambda c: c.confidence if c.confidence and c.confidence > 50 else "",
                      label_colors=lambda c: 'red' if c.confidence and c.confidence < 70 else 'black')
            
            # Add title
            plt.title(f"Phylogenetic Tree - {gene_name}" if gene_name else "Phylogenetic Tree", 
                     fontsize=16, pad=20)
            
            # Save
            output_file = output_dir / f"{Path(treefile).stem}_matplotlib.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            visualizations_created.append(str(output_file))
            log(f"Matplotlib tree visualization saved to {output_file}")
            
        except Exception as e:
            warn(f"Could not create matplotlib visualization: {e}")
    
    # 2. Try ETE3 visualization
    if ETE3_AVAILABLE:
        try:
            t = Tree(treefile)
            
            # Create ASCII representation
            ascii_file = output_dir / f"{Path(treefile).stem}_ascii.txt"
            with open(ascii_file, 'w') as f:
                f.write(f"ASCII Tree Representation - {gene_name}\n")
                f.write("=" * 50 + "\n\n")
                f.write(t.get_ascii(show_internal=True, attributes=["support"]))
            
            visualizations_created.append(str(ascii_file))
            log(f"ASCII tree saved to {ascii_file}")
            
            # Try to create image if display is available
            try:
                ts = TreeStyle()
                ts.show_leaf_name = True
                ts.show_branch_support = True
                ts.branch_vertical_margin = 10
                
                image_file = output_dir / f"{Path(treefile).stem}_ete3.png"
                t.render(str(image_file), tree_style=ts, dpi=300)
                visualizations_created.append(str(image_file))
                log(f"ETE3 tree image saved to {image_file}")
            except:
                pass  # Image rendering might fail on headless systems
                
        except Exception as e:
            warn(f"Could not create ETE3 visualization: {e}")
    
    # 3. Create an HTML viewer
    try:
        html_file = output_dir / f"{Path(treefile).stem}_viewer.html"
        with open(treefile, 'r') as f:
            newick_string = f.read().strip()
        
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Phylogenetic Tree - {gene_name}</title>
    <script src="https://d3js.org/d3.v3.min.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/veg/phylotree.js@master/phylotree.js"></script>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/veg/phylotree.js@master/phylotree.css">
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        #tree_display {{ border: 1px solid #ddd; }}
    </style>
</head>
<body>
    <h1>Phylogenetic Tree - {gene_name}</h1>
    <p>Tree file: {Path(treefile).name}</p>
    <svg id="tree_display" width="900" height="700"></svg>
    <script>
        var newick_string = `{newick_string}`;
        var tree = d3.layout.phylotree()
            .svg(d3.select("#tree_display"))
            .options({{
                'show-bootstrap': true,
                'left-right-spacing': 'fit-to-size',
                'top-bottom-spacing': 'fit-to-size',
                'show-scale': true
            }});
        tree(d3_phylotree_newick_parser(newick_string))
            .layout();
    </script>
</body>
</html>"""
        
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        visualizations_created.append(str(html_file))
        log(f"Interactive HTML tree viewer saved to {html_file}")
        
    except Exception as e:
        warn(f"Could not create HTML viewer: {e}")
    
    return visualizations_created

def codon_changes(aln_records):
    """
    Codon analysis with better handling of alignments.
    """
    out = []
    if not aln_records:
        return out
        
    length = len(aln_records[0].seq)
    
    for pos in range(0, length - 2, 3):
        codons = {}
        for rec in aln_records:
            codon = str(rec.seq[pos:pos+3]).upper()
            if len(codon) < 3:
                continue
            if "-" in codon:
                codons[rec.id] = "---"
            else:
                codons[rec.id] = codon
                
        unique_codons = {c for c in codons.values() if c != "---" and 'N' not in c}
        if len(unique_codons) <= 1:
            continue
            
        aa_map = {}
        for sid, c in codons.items():
            if c == "---":
                aa_map[sid] = "-"
            elif 'N' in c:
                aa_map[sid] = "?"
            else:
                try:
                    aa_map[sid] = str(Seq(c).translate())
                except:
                    aa_map[sid] = "?"
                    
        unique_aas = {a for a in aa_map.values() if a not in ["-", "?"]}
        effect = "synonymous" if len(unique_aas) == 1 else "nonsynonymous"
        
        out.append({
            "Codon_index_1based": pos//3 + 1,
            "Alignment_nt_positions": f"{pos+1}-{pos+3}",
            "Effect": effect,
            "Codons": ";".join([f"{sid}:{codons[sid]}" for sid in sorted(codons)]),
            "AAs": ";".join([f"{sid}:{aa_map[sid]}" for sid in sorted(aa_map)])
        })
    return out

def run_mafft(fasta_in, fasta_out):
    """Run MAFFT alignment"""
    if not which("mafft"):
        warn("MAFFT not found; skipping alignment.")
        shutil.copy(fasta_in, fasta_out)
        return False
        
    try:
        cmd = ["mafft", "--auto", "--quiet", fasta_in]
        with open(fasta_out, "w") as fout:
            result = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, check=True)
        log(f"MAFFT alignment completed: {fasta_out}")
        return True
    except subprocess.CalledProcessError as e:
        warn(f"MAFFT failed: {e}")
        shutil.copy(fasta_in, fasta_out)
        return False

def run_iqtree(aln_fasta, outdir):
    """Run IQ-TREE - compatible with both iqtree and iqtree2"""
    iqtree_cmd = None
    if which("iqtree2"):
        iqtree_cmd = "iqtree2"
    elif which("iqtree"):
        iqtree_cmd = "iqtree"
    else:
        warn("IQ-TREE not found; skipping tree.")
        return False
        
    try:
        if iqtree_cmd == "iqtree":
            cmd = [iqtree_cmd, "-s", aln_fasta, "-m", "TEST", "-bb", "1000", "-nt", "AUTO", "-quiet"]
        else:
            cmd = [iqtree_cmd, "-s", aln_fasta, "-m", "MFP", "-bb", "1000", "-nt", "AUTO", "-quiet"]
            
        log(f"Running {iqtree_cmd} for phylogenetic analysis...")
        result = subprocess.run(cmd, cwd=outdir, capture_output=True, text=True)
        
        if result.returncode != 0:
            warn(f"IQ-TREE failed: {result.stderr}")
            return False
            
        log("IQ-TREE analysis completed")
        return True
        
    except Exception as e:
        warn(f"IQ-TREE error: {e}")
        return False

# ---------- Main Pipeline ----------

def main():
    args = parse_args()
    in_dir = Path(args.input)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_dir.exists():
        fail(f"Input directory {in_dir} not found")

    log(f"Starting Sanger pipeline v1.3...")
    log(f"Input directory: {in_dir}")
    log(f"Output directory: {out_dir}")
    log(f"Quality threshold: {args.quality_threshold}")
    log(f"Window size: {args.window_size}")
    log(f"Strict cleaning: {'ENABLED' if args.strict_clean else 'DISABLED'}")
    if args.aggressive_trim:
        log("Aggressive trimming: ENABLED")

    # Parse AB1 files
    groups = defaultdict(lambda: defaultdict(dict))
    unmatched = []
    total_files = 0
    
    for path in in_dir.glob("*.ab1"):
        total_files += 1
        m = FILENAME_REGEX.match(path.name)
        if not m:
            unmatched.append(path.name)
            continue
            
        sample = m.group("sample")
        gene = m.group("gene").upper().replace("RBCLA","RBCL")
        direction = m.group("dir").upper()
        
        if args.verbose:
            log(f"Processing {path.name}: sample={sample}, gene={gene}, dir={direction}")
            
        # Read with quality scores
        raw_seq, quality = read_ab1_with_quality(str(path))
        if not raw_seq:
            warn(f"Failed to read sequence from {path.name}")
            continue
            
        # Quality-based trimming
        trimmed = quality_trim_sequence(raw_seq, quality, args)
        
        if args.verbose:
            log(f"  Raw length: {len(raw_seq)}, Trimmed length: {len(trimmed)}")
            if trimmed:
                n_count = trimmed.count('N')
                n_pct = (n_count/len(trimmed)*100)
                log(f"  N count: {n_count}, N%: {n_pct:.1f}%")
                # Check for non-ATCGN characters
                non_atcgn = set(trimmed) - set('ATCGN')
                if non_atcgn:
                    log(f"  Non-ATCGN characters found: {non_atcgn}")
            
        groups[gene][sample][direction] = {
            'seq': trimmed,
            'qual': quality[:len(trimmed)] if len(trimmed) < len(quality) else quality
        }

    log(f"Found {total_files} AB1 files")
    
    if unmatched:
        warn(f"{len(unmatched)} files did not match naming pattern (see unmatched_files.txt)")
        with open(out_dir/"unmatched_files.txt","w") as f:
            f.write("\n".join(unmatched))

    if args.genes:
        wanted = {g.strip().upper() for g in args.genes.split(",")}
        groups = {g: samples for g, samples in groups.items() if g in wanted}

    if not groups:
        fail("No valid AB1 files found matching the expected pattern")

    gene_reports = []
    all_tree_files = []

    for gene, samples in groups.items():
        log(f"\nProcessing gene: {gene} (samples: {len(samples)})")
        gene_dir = out_dir / gene
        gene_dir.mkdir(exist_ok=True)
        
        consensus_records = []
        per_sample_rows = []
        
        for sample, reads in samples.items():
            f_data = reads.get("F", {})
            r_data = reads.get("R", {})
            f_seq = f_data.get('seq', '')
            r_seq = r_data.get('seq', '')
            f_qual = f_data.get('qual', [])
            r_qual = r_data.get('qual', [])
            
            consensus = ""
            ov_len = ident = avg_qual = 0.0
            warn_flag = ""
            
            if f_seq and r_seq:
                # Reverse complement
                r_rc = str(Seq(r_seq).reverse_complement())
                r_qual_rc = r_qual[::-1] if r_qual else []
                
                # Find overlap with quality consideration
                offset, ov_len, ident, avg_qual = pairwise_overlap_quality(
                    f_seq, r_rc, f_qual, r_qual_rc, min_overlap=args.min_overlap
                )
                
                if ov_len >= args.min_overlap and ident >= 0.85:
                    consensus = build_consensus_quality(f_seq, r_rc, f_qual, r_qual_rc, offset, args.strict_clean)
                    if args.verbose:
                        log(f"  {sample}: F={len(f_seq)}bp, R={len(r_seq)}bp, overlap={ov_len}bp, identity={ident:.2%}, avg_qual={avg_qual:.1f}")
                else:
                    warn_flag = "LOW_OVERLAP" if ov_len < args.min_overlap else "LOW_IDENT"
                    # Use higher quality read
                    if f_qual and r_qual:
                        consensus = f_seq if np.mean(f_qual) >= np.mean(r_qual) else r_rc
                    else:
                        consensus = f_seq if len(f_seq) >= len(r_rc) else r_rc
            else:
                consensus = f_seq or (str(Seq(r_seq).reverse_complement()) if r_seq else "")
                warn_flag = "ONLY_ONE_READ"

            # Final cleaning and N removal
            if consensus:
                # Ensure strict cleaning
                if args.strict_clean:
                    consensus = clean_sequence(consensus, strict=True)
                
                # Final N removal
                if args.aggressive_trim:
                    consensus = remove_n_stretches(consensus)
                
                consensus = consensus.strip("N")
            
            if not consensus:
                warn(f"  {sample}: No consensus sequence generated")
                continue
                
            n_count = consensus.count("N")
            n_pct = (n_count/len(consensus)*100) if consensus else 0
            
            # Check for any remaining non-ATCGN
            non_atcgn = set(consensus) - set('ATCGN')
            if non_atcgn:
                warn(f"  {sample}: Non-ATCGN characters in consensus: {non_atcgn}")
            
            if len(consensus) < MIN_FINAL_LENGTH * 2:
                warn_flag = (warn_flag + ";SHORT") if warn_flag else "SHORT"

            rec_id = f"{sample}|{gene}"
            consensus_records.append(
                SeqRecord(Seq(consensus), 
                         id=rec_id, 
                         description=f"len={len(consensus)} Ns={n_count} N%={n_pct:.1f}%")
            )
            
            per_sample_rows.append({
                "Gene": gene,
                "Sample": sample,
                "Forward_len": len(f_seq),
                "Reverse_len": len(r_seq),
                "Consensus_len": len(consensus),
                "Overlap_bp": int(ov_len),
                "Overlap_identity": round(ident, 4),
                "Avg_overlap_quality": round(avg_qual, 1),
                "N_count": n_count,
                "%N_consensus": round(n_pct, 2),
                "Flags": warn_flag
            })

        if not consensus_records:
            warn(f"No consensus sequences generated for gene {gene}")
            continue

        # Write consensus FASTA
        fasta_path = gene_dir / f"{gene}_consensus.fasta"
        SeqIO.write(consensus_records, fasta_path, "fasta")
        log(f"Wrote {len(consensus_records)} consensus sequences to {fasta_path}")

        # Quality statistics
        n_percentages = [float(r.description.split('N%=')[1].split('%')[0]) for r in consensus_records]
        log(f"N% statistics: min={min(n_percentages):.1f}%, median={statistics.median(n_percentages):.1f}%, max={max(n_percentages):.1f}%")

        # Continue with alignment and analysis...
        if not args.no_align:
            aln_path = gene_dir / f"{gene}_aligned.fasta"
            alignment_success = run_mafft(str(fasta_path), str(aln_path))
            
            if alignment_success:
                aln_records = list(SeqIO.parse(aln_path, "fasta"))

                # SNP table
                snp_rows = []
                if aln_records and len(aln_records) > 1:
                    aln_len = len(aln_records[0].seq)
                    for pos in range(aln_len):
                        column = [str(r.seq[pos]) for r in aln_records]
                        # Only consider positions without N or gaps
                        uniq = {b for b in column if b not in ["N","-"]}
                        if len(uniq) > 1:
                            snp_rows.append({
                                "Gene": gene,
                                "Alignment_Position": pos+1,
                                **{rec.id: rec.seq[pos] for rec in aln_records}
                            })
                    
                    if snp_rows:
                        snp_df = pd.DataFrame(snp_rows)
                        snp_df.to_csv(gene_dir / f"{gene}_SNPs_nt.csv", index=False)
                        log(f"Found {len(snp_rows)} SNP positions")

                # Codon effects
                codon_rows = codon_changes(aln_records)
                if codon_rows:
                    codon_df = pd.DataFrame(codon_rows)
                    codon_df.to_csv(gene_dir / f"{gene}_codon_effects.csv", index=False)
                    log(f"Analyzed {len(codon_rows)} variable codon positions")

                # IQ-TREE and visualization
                if not args.no_tree and len(aln_records) > 2:
                    if run_iqtree(str(aln_path), gene_dir):
                        treefile = gene_dir / f"{gene}_aligned.fasta.treefile"
                        if treefile.exists():
                            all_tree_files.append(str(treefile))
                            # Create visualizations
                            if not args.no_viz:
                                viz_files = create_tree_visualization(str(treefile), gene_dir, gene)
                                if viz_files:
                                    log(f"Created {len(viz_files)} tree visualization(s)")
        else:
            snp_rows = []
            codon_rows = []

        # Sample QC table
        qc_df = pd.DataFrame(per_sample_rows)
        qc_df.to_csv(gene_dir / f"{gene}_read_report.csv", index=False)

        # Gene summary
        if consensus_records:
            lengths = [len(r.seq) for r in consensus_records]
            gene_reports.append({
                "Gene": gene,
                "N_samples": len(consensus_records),
                "Min_consensus_len": min(lengths),
                "Median_consensus_len": int(statistics.median(lengths)),
                "Max_consensus_len": max(lengths),
                "Mean_consensus_len": round(statistics.mean(lengths), 1),
                "Mean_N%": round(statistics.mean(n_percentages), 2),
                "Max_N%": round(max(n_percentages), 2),
                "SNP_sites": len(snp_rows),
                "Nonsynonymous_sites": sum(1 for r in codon_rows if r.get("Effect") == "nonsynonymous"),
                "Samples_with_flags": sum(1 for r in per_sample_rows if r.get("Flags"))
            })

    # Master summary
    if gene_reports:
        master_df = pd.DataFrame(gene_reports)
        master_df.to_csv(out_dir / "MASTER_summary.csv", index=False)
        log(f"\nWrote master summary to {out_dir / 'MASTER_summary.csv'}")

        # Enhanced HTML report
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Sanger Pipeline Report - Quality Enhanced v1.3</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 2px solid #4CAF50; padding-bottom: 10px; }}
        h2 {{ color: #555; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .timestamp {{ color: #666; font-size: 0.9em; }}
        .warning {{ color: #ff6600; font-weight: bold; }}
        .good {{ color: #009900; font-weight: bold; }}
        .info-box {{ background-color: #e7f3ff; border-left: 4px solid #2196F3; padding: 10px; margin: 20px 0; }}
        .tree-section {{ background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Sanger Pipeline Summary Report - Quality Enhanced v1.3</h1>
        <p class="timestamp">Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="info-box">
            <h3>Pipeline Parameters</h3>
            <ul>
                <li>Input directory: {in_dir}</li>
                <li>Quality threshold: {args.quality_threshold}</li>
                <li>Window size: {args.window_size}</li>
                <li>Strict sequence cleaning: <span class="{'good' if args.strict_clean else 'warning'}">{'ENABLED' if args.strict_clean else 'DISABLED'}</span></li>
                <li>Aggressive trimming: {'Yes' if args.aggressive_trim else 'No'}</li>
                <li>Total genes processed: {len(gene_reports)}</li>
            </ul>
        </div>
        
        <h2>Gene Summary</h2>
        {master_df.to_html(index=False, classes='summary-table')}
        
        <h2>Quality Metrics</h2>
        <ul>
            <li>Sequences with N% > 5%: <span class="{'warning' if any(df['Max_N%'] > 5 for _, df in master_df.iterrows()) else 'good'}">
                {sum(1 for _, df in master_df.iterrows() if df['Max_N%'] > 5)}
            </span></li>
            <li>Total SNPs detected: {sum(df['SNP_sites'] for _, df in master_df.iterrows())}</li>
            <li>Total nonsynonymous mutations: {sum(df['Nonsynonymous_sites'] for _, df in master_df.iterrows())}</li>
        </ul>
        
        {'<div class="tree-section"><h2>Phylogenetic Trees</h2><p>Tree files generated:</p><ul>' + ''.join(f'<li>{Path(f).name}</li>' for f in all_tree_files) + '</ul><p>View trees using FigTree or the included HTML viewers.</p></div>' if all_tree_files else ''}
        
        <h2>Output Files</h2>
        <ul>
            <li><strong>Per-gene directories</strong> containing:
                <ul>
                    <li>Consensus sequences (quality-trimmed, strict ATCG+N only)</li>
                    <li>Multiple sequence alignment</li>
                    <li>SNP positions (excluding N positions)</li>
                    <li>Codon effects (synonymous/nonsynonymous)</li>
                    <li>Read quality report (including overlap quality)</li>
                    {'<li>Phylogenetic tree files (.treefile, .iqtree, .log)</li>' if not args.no_tree else ''}
                    {'<li>Tree visualizations (PNG, HTML, ASCII)</li>' if all_tree_files and not args.no_viz else ''}
                </ul>
            </li>
        </ul>
        
        <h2>Notes</h2>
        <ul>
            <li>All sequences have been cleaned to contain only ATCG and N characters</li>
            <li>Quality trimming based on Phred scores with window size {args.window_size}</li>
            <li>Minimum quality threshold: {args.quality_threshold}</li>
            <li>N-stretches longer than {MAX_N_STRETCH} bases have been collapsed</li>
        </ul>
    </div>
</body>
</html>
"""
        with open(out_dir / "report.html", "w") as fh:
            fh.write(html_content)
        log(f"Wrote HTML report to {out_dir / 'report.html'}")

    log("\nPipeline completed successfully!")
    log(f"Results saved to: {out_dir}")
    
    # Final summary of visualizations
    if all_tree_files and not args.no_viz:
        log("\nTree visualizations created:")
        log("- PNG files: View with any image viewer")
        log("- HTML files: Open in web browser for interactive viewing")
        log("- ASCII files: View in any text editor")
        log("\nFor publication-quality trees, use FigTree with the .treefile files")

if __name__ == "__main__":
    main()