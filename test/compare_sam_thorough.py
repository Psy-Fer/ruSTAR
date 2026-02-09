#!/usr/bin/env python3
"""
Thorough comparison of ruSTAR vs STAR SAM output.

Usage:
    python3 compare_sam_thorough.py <rustar_dir> <star_dir>

Both directories should contain Aligned.out.sam and SJ.out.tab.
"""

import re
import sys
import os
from collections import defaultdict, Counter

def parse_sam(path):
    """Parse SAM file, return dict of read_name -> list of alignment records."""
    reads = defaultdict(list)
    header_lines = 0
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                header_lines += 1
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            reads[qname].append({
                "flag": flag,
                "rname": rname,
                "pos": pos,
                "mapq": mapq,
                "cigar": cigar,
            })
    return reads, header_lines


def classify_read(records):
    """Classify a read as unique/multi/unmapped based on its records."""
    if len(records) == 1 and (records[0]["flag"] & 4):
        return "unmapped"
    for r in records:
        if not (r["flag"] & 256):  # not secondary
            if r["flag"] & 4:
                return "unmapped"
            if r["mapq"] == 255:
                return "unique"
            return "multi"
    return "unmapped"


def get_primary(records):
    """Get primary alignment from records."""
    for r in records:
        if not (r["flag"] & 256) and not (r["flag"] & 4):
            return r
    return None


def cigar_category(cigar):
    """Categorize a CIGAR string."""
    cats = set()
    if cigar == "*":
        return {"unmapped"}
    ops = re.findall(r'\d+([MIDNSHP=X])', cigar)
    if "N" in ops:
        cats.add("spliced")
    if "I" in ops or "D" in ops:
        cats.add("indel")
    if "S" in ops:
        cats.add("softclip")
    if not cats:
        cats.add("pure_match")
    return cats


def parse_sj(path):
    """Parse SJ.out.tab, return dict of (chr, start, end) -> record."""
    junctions = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = int(fields[3])
            motif = int(fields[4])
            annotated = int(fields[5])
            uniq_count = int(fields[6])
            multi_count = int(fields[7])
            max_overhang = int(fields[8])
            key = (chrom, start, end)
            junctions[key] = {
                "strand": strand,
                "motif": motif,
                "annotated": annotated,
                "uniq_count": uniq_count,
                "multi_count": multi_count,
                "max_overhang": max_overhang,
            }
    return junctions


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <rustar_dir> <star_dir>")
        sys.exit(1)

    rustar_dir = sys.argv[1]
    star_dir = sys.argv[2]

    rustar_sam = os.path.join(rustar_dir, "Aligned.out.sam")
    star_sam = os.path.join(star_dir, "Aligned.out.sam")
    rustar_sj_path = os.path.join(rustar_dir, "SJ.out.tab")
    star_sj_path = os.path.join(star_dir, "SJ.out.tab")

    for f in [rustar_sam, star_sam]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            sys.exit(1)

    has_sj = os.path.exists(rustar_sj_path) and os.path.exists(star_sj_path)

    # ============================================================
    # PARSE FILES
    # ============================================================
    print("=" * 80)
    print("COMPREHENSIVE COMPARISON: ruSTAR vs STAR")
    print("=" * 80)
    print(f"\nruSTAR dir: {rustar_dir}")
    print(f"STAR dir:   {star_dir}")

    print("\nParsing SAM files...")
    rustar_reads, rustar_headers = parse_sam(rustar_sam)
    star_reads, star_headers = parse_sam(star_sam)

    print(f"  ruSTAR: {len(rustar_reads)} unique read names, {rustar_headers} header lines")
    print(f"  STAR:   {len(star_reads)} unique read names, {star_headers} header lines")

    # ============================================================
    # 1. MAPPING STATS COMPARISON
    # ============================================================
    print("\n" + "=" * 80)
    print("1. MAPPING STATS COMPARISON")
    print("=" * 80)

    rustar_class = Counter()
    star_class = Counter()

    for qname, records in rustar_reads.items():
        rustar_class[classify_read(records)] += 1

    for qname, records in star_reads.items():
        star_class[classify_read(records)] += 1

    rustar_total = sum(rustar_class.values())
    star_total = sum(star_class.values())

    print(f"\n{'Category':<15} {'ruSTAR':>10} {'%':>8} {'STAR':>10} {'%':>8} {'Diff':>8}")
    print("-" * 65)
    for cat in ["unique", "multi", "unmapped"]:
        rc = rustar_class[cat]
        sc = star_class[cat]
        rp = 100.0 * rc / rustar_total if rustar_total else 0
        sp = 100.0 * sc / star_total if star_total else 0
        print(f"{cat:<15} {rc:>10} {rp:>7.1f}% {sc:>10} {sp:>7.1f}% {rc - sc:>+8}")
    print(f"{'TOTAL':<15} {rustar_total:>10} {'':>8} {star_total:>10}")

    # ============================================================
    # 2. PER-READ AGREEMENT
    # ============================================================
    print("\n" + "=" * 80)
    print("2. PER-READ AGREEMENT")
    print("=" * 80)

    all_reads = set(rustar_reads.keys()) | set(star_reads.keys())

    both_mapped_agree_pos = 0
    both_mapped_agree_exact = 0
    both_mapped_disagree_pos = 0
    star_only_mapped = 0
    rustar_only_mapped = 0
    both_unmapped = 0
    both_mapped_agree_strand = 0
    both_mapped_agree_cigar = 0

    disagree_examples = []

    for qname in sorted(all_reads):
        r_records = rustar_reads.get(qname, [])
        s_records = star_reads.get(qname, [])

        r_class = classify_read(r_records) if r_records else "missing"
        s_class = classify_read(s_records) if s_records else "missing"

        r_mapped = r_class in ("unique", "multi")
        s_mapped = s_class in ("unique", "multi")

        if r_mapped and s_mapped:
            r_pri = get_primary(r_records)
            s_pri = get_primary(s_records)
            if r_pri and s_pri:
                same_chr = r_pri["rname"] == s_pri["rname"]
                same_pos = abs(r_pri["pos"] - s_pri["pos"]) <= 5
                same_strand = (r_pri["flag"] & 16) == (s_pri["flag"] & 16)
                same_cigar = r_pri["cigar"] == s_pri["cigar"]

                if same_chr and same_pos:
                    both_mapped_agree_pos += 1
                    if same_strand:
                        both_mapped_agree_strand += 1
                    if same_cigar:
                        both_mapped_agree_cigar += 1
                    if same_chr and same_pos and same_strand and same_cigar:
                        both_mapped_agree_exact += 1
                else:
                    both_mapped_disagree_pos += 1
                    if len(disagree_examples) < 30:
                        disagree_examples.append({
                            "qname": qname,
                            "rustar_chr": r_pri["rname"],
                            "rustar_pos": r_pri["pos"],
                            "rustar_strand": "-" if (r_pri["flag"] & 16) else "+",
                            "rustar_cigar": r_pri["cigar"],
                            "rustar_mapq": r_pri["mapq"],
                            "star_chr": s_pri["rname"],
                            "star_pos": s_pri["pos"],
                            "star_strand": "-" if (s_pri["flag"] & 16) else "+",
                            "star_cigar": s_pri["cigar"],
                            "star_mapq": s_pri["mapq"],
                        })
            else:
                both_mapped_disagree_pos += 1
        elif r_mapped and not s_mapped:
            rustar_only_mapped += 1
        elif not r_mapped and s_mapped:
            star_only_mapped += 1
        else:
            both_unmapped += 1

    total_both_mapped = both_mapped_agree_pos + both_mapped_disagree_pos

    print(f"\n{'Category':<48} {'Count':>8} {'%':>8}")
    print("-" * 68)
    print(f"{'Both mapped, agree position (<=5bp)':<48} {both_mapped_agree_pos:>8} {100.0 * both_mapped_agree_pos / len(all_reads):>7.1f}%")
    print(f"{'  ...of which agree strand':<48} {both_mapped_agree_strand:>8} {100.0 * both_mapped_agree_strand / max(both_mapped_agree_pos, 1):>7.1f}%")
    print(f"{'  ...of which agree CIGAR (exact)':<48} {both_mapped_agree_cigar:>8} {100.0 * both_mapped_agree_cigar / max(both_mapped_agree_pos, 1):>7.1f}%")
    print(f"{'  ...of which agree ALL (chr+pos+strand+cigar)':<48} {both_mapped_agree_exact:>8} {100.0 * both_mapped_agree_exact / max(both_mapped_agree_pos, 1):>7.1f}%")
    print(f"{'Both mapped, DISAGREE position':<48} {both_mapped_disagree_pos:>8} {100.0 * both_mapped_disagree_pos / len(all_reads):>7.1f}%")
    print(f"{'STAR only mapped':<48} {star_only_mapped:>8} {100.0 * star_only_mapped / len(all_reads):>7.1f}%")
    print(f"{'ruSTAR only mapped':<48} {rustar_only_mapped:>8} {100.0 * rustar_only_mapped / len(all_reads):>7.1f}%")
    print(f"{'Both unmapped':<48} {both_unmapped:>8} {100.0 * both_unmapped / len(all_reads):>7.1f}%")
    print("-" * 68)
    print(f"{'Total reads':<48} {len(all_reads):>8}")

    if total_both_mapped > 0:
        concordance = 100.0 * both_mapped_agree_pos / total_both_mapped
        print(f"\nPosition concordance (among both-mapped): {concordance:.1f}%")

    # ============================================================
    # 3. JUNCTION COMPARISON
    # ============================================================
    if has_sj:
        print("\n" + "=" * 80)
        print("3. JUNCTION COMPARISON (SJ.out.tab)")
        print("=" * 80)

        rustar_sj = parse_sj(rustar_sj_path)
        star_sj = parse_sj(star_sj_path)

        rustar_keys = set(rustar_sj.keys())
        star_keys = set(star_sj.keys())

        shared = rustar_keys & star_keys
        rustar_only_sj = rustar_keys - star_keys
        star_only_sj = star_keys - rustar_keys

        print(f"\n{'Category':<30} {'Count':>8}")
        print("-" * 40)
        print(f"{'Shared junctions':<30} {len(shared):>8}")
        print(f"{'STAR-only junctions':<30} {len(star_only_sj):>8}")
        print(f"{'ruSTAR-only junctions':<30} {len(rustar_only_sj):>8}")
        print(f"{'Total STAR junctions':<30} {len(star_keys):>8}")
        print(f"{'Total ruSTAR junctions':<30} {len(rustar_keys):>8}")

        # Motif comparison for shared junctions
        motif_names = {0: "non-canonical", 1: "GT/AG", 2: "CT/AC", 3: "GC/AG", 4: "CT/GC", 5: "AT/AC", 6: "GT/AT"}
        motif_agree = 0
        motif_disagree = 0
        motif_disagree_examples = []

        for key in shared:
            r_motif = rustar_sj[key]["motif"]
            s_motif = star_sj[key]["motif"]
            if r_motif == s_motif:
                motif_agree += 1
            else:
                motif_disagree += 1
                if len(motif_disagree_examples) < 5:
                    motif_disagree_examples.append((key, r_motif, s_motif))

        if shared:
            print(f"\nMotif agreement for shared junctions: {motif_agree}/{len(shared)} ({100.0*motif_agree/len(shared):.1f}%)")
            if motif_disagree_examples:
                print("\nMotif disagreements (first 5):")
                for key, rm, sm in motif_disagree_examples:
                    print(f"  {key[0]}:{key[1]}-{key[2]}  ruSTAR={motif_names.get(rm, str(rm))}  STAR={motif_names.get(sm, str(sm))}")

        # Coverage comparison for shared junctions
        if shared:
            print(f"\nCoverage comparison for shared junctions (top 20 by STAR unique count):")
            print(f"  {'Junction':<30} {'rSTAR_uniq':>11} {'STAR_uniq':>11} {'rSTAR_multi':>12} {'STAR_multi':>11}")
            print("  " + "-" * 78)
            shared_sorted = sorted(shared, key=lambda k: star_sj[k]["uniq_count"], reverse=True)
            for key in shared_sorted[:20]:
                r = rustar_sj[key]
                s = star_sj[key]
                jstr = f"{key[0]}:{key[1]}-{key[2]}"
                print(f"  {jstr:<30} {r['uniq_count']:>11} {s['uniq_count']:>11} {r['multi_count']:>12} {s['multi_count']:>11}")

        # Show STAR-only junctions
        if star_only_sj:
            print(f"\nSTAR-only junctions ({len(star_only_sj)}):")
            star_only_sorted = sorted(star_only_sj, key=lambda k: star_sj[k]["uniq_count"], reverse=True)
            for key in star_only_sorted[:10]:
                s = star_sj[key]
                jstr = f"{key[0]}:{key[1]}-{key[2]}"
                print(f"  {jstr:<30} motif={motif_names.get(s['motif'], str(s['motif'])):<15} uniq={s['uniq_count']:<5} multi={s['multi_count']:<5} annot={s['annotated']}")

        # Show ruSTAR-only junctions (first 20 by count)
        if rustar_only_sj:
            print(f"\nruSTAR-only junctions (top 20 of {len(rustar_only_sj)} by unique count):")
            rustar_only_sorted = sorted(rustar_only_sj, key=lambda k: rustar_sj[k]["uniq_count"], reverse=True)
            for key in rustar_only_sorted[:20]:
                r = rustar_sj[key]
                jstr = f"{key[0]}:{key[1]}-{key[2]}"
                print(f"  {jstr:<30} motif={motif_names.get(r['motif'], str(r['motif'])):<15} uniq={r['uniq_count']:<5} multi={r['multi_count']:<5} annot={r['annotated']}")

        # Summarize ruSTAR-only by motif
        if rustar_only_sj:
            motif_counter = Counter()
            for key in rustar_only_sj:
                motif_counter[rustar_sj[key]["motif"]] += 1
            print(f"\nruSTAR-only junctions by motif:")
            for motif_id, count in motif_counter.most_common():
                print(f"  {motif_names.get(motif_id, f'unknown({motif_id})'):<20} {count:>5}")
    else:
        print("\n(Skipping junction comparison - SJ.out.tab not found in both directories)")

    # ============================================================
    # 4. CIGAR PATTERN DISTRIBUTION
    # ============================================================
    print("\n" + "=" * 80)
    print("4. CIGAR PATTERN DISTRIBUTION")
    print("=" * 80)

    for label, reads in [("ruSTAR", rustar_reads), ("STAR", star_reads)]:
        cats = Counter()
        total_mapped = 0
        for qname, records in reads.items():
            pri = get_primary(records)
            if pri:
                total_mapped += 1
                for cat in cigar_category(pri["cigar"]):
                    cats[cat] += 1

        print(f"\n{label} (mapped reads: {total_mapped}):")
        print(f"  {'Category':<20} {'Count':>8} {'%':>8}")
        print("  " + "-" * 40)
        for cat in ["pure_match", "spliced", "indel", "softclip", "unmapped"]:
            c = cats.get(cat, 0)
            pct = 100.0 * c / total_mapped if total_mapped else 0
            print(f"  {cat:<20} {c:>8} {pct:>7.1f}%")
        print(f"  (Note: categories overlap - a read can be spliced + softclipped)")

    # ============================================================
    # 5. EXAMPLE DISAGREEMENTS
    # ============================================================
    print("\n" + "=" * 80)
    print("5. EXAMPLE READS WITH POSITION DISAGREEMENT (first 10)")
    print("=" * 80)

    if not disagree_examples:
        print("\nNo position disagreements found!")
    else:
        for i, ex in enumerate(disagree_examples[:10]):
            print(f"\n--- Read {i+1}: {ex['qname']} ---")
            print(f"  ruSTAR: {ex['rustar_chr']}:{ex['rustar_pos']} ({ex['rustar_strand']}) MAPQ={ex['rustar_mapq']} CIGAR={ex['rustar_cigar']}")
            print(f"  STAR:   {ex['star_chr']}:{ex['star_pos']} ({ex['star_strand']}) MAPQ={ex['star_mapq']} CIGAR={ex['star_cigar']}")

            if ex['rustar_chr'] != ex['star_chr']:
                print(f"  >> Different chromosome!")
            else:
                diff = abs(ex['rustar_pos'] - ex['star_pos'])
                print(f"  >> Same chr, position difference: {diff}bp")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    total = len(all_reads)
    print(f"""
Total reads:               {total}
Both mapped, agree:        {both_mapped_agree_pos} ({100.0*both_mapped_agree_pos/total:.1f}%)
Both mapped, disagree:     {both_mapped_disagree_pos} ({100.0*both_mapped_disagree_pos/total:.1f}%)
STAR only mapped:          {star_only_mapped} ({100.0*star_only_mapped/total:.1f}%)
ruSTAR only mapped:        {rustar_only_mapped} ({100.0*rustar_only_mapped/total:.1f}%)
Both unmapped:             {both_unmapped} ({100.0*both_unmapped/total:.1f}%)""")

    if has_sj:
        print(f"""
Junctions shared:          {len(shared)}/{len(star_keys)} STAR junctions ({100.0*len(shared)/max(len(star_keys),1):.1f}%)
ruSTAR extra junctions:    {len(rustar_only_sj)} (potential false positives)""")
    print()


if __name__ == "__main__":
    main()
