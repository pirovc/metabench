#!/usr/bin/env python3

import sys
import gzip
import json
import argparse
import numpy as np
from collections import defaultdict

from util import default_ranks, parse_tax, file_exists


def main():

    parser = argparse.ArgumentParser(
        prog="ganon benchmark evaluation", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input-results", metavar="", type=file_exists,
                        help="bioboxes profile TAXID <tab> RANK <tab> TAXPATH <tab> TAXPATHSN <tab> PERCENTAGE")
    parser.add_argument("-g", "--input-ground-truth", metavar="",
                        type=file_exists, help="bioboxes profile TAXID <tab> RANK <tab> TAXPATH <tab> TAXPATHSN <tab> PERCENTAGE")
    parser.add_argument("-d", "--input-database-profile", metavar="",
                        type=file_exists, help="taxid <tab> parent <tab> rank <tab> name")
    parser.add_argument("-t", "--taxonomy",       metavar="",
                        type=str, help="custom, ncbi or gtdb")
    parser.add_argument("-x", "--taxonomy-files", metavar="",
                        type=file_exists, nargs="*", help="")
    parser.add_argument("-r", "--ranks",          metavar="", type=str, nargs="*",
                        help="empty for default ranks (superkingdom phylum class order family genus species)")
    parser.add_argument("-t", "--thresholds", nargs="*", type=float, default=[0],
                        help="Thresholds to generate evaluations [0-100]")
    parser.add_argument("-j", "--output-json", type=str,
                        help="Output file for evaluation in json (cumulative mode)")
    args = parser.parse_args()

    tax = parse_tax(args.taxonomy, args.taxonomy_files)
    fixed_ranks = args.ranks if args.ranks else default_ranks[1:]

    gt_file = args.input_ground_truth
    db_file = args.input_database_profile
    res_file = args.input_results

    output_json = open(args.output_json, "w") if args.output_json else None

    # Ground truth: TAXID <tab> RANK <tab> TAXPATH <tab> TAXPATHSN <tab> PERCENTAGE
    gt = dict()
    for line in gzip.open(gt_file, "rt") if gt_file.endswith(".gz") else open(gt_file, "r"):
        if line[0] == "@":
            continue
        fields = line.rstrip().split("\t")
        taxid = tax.latest(fields[0])
        rank = tax.rank(taxid)
        if taxid == tax.undefined_node:
            sys.stderr.write(gt_file + ": " +
                             fields[0] + " not found in taxonomy\n")
            continue
        if rank not in fixed_ranks:
            sys.stderr.write(gt_file + ": " +
                             fields[0] + " not in defined ranks\n")
            continue
        if rank not in gt:
            gt[rank] = {}
        gt[rank][taxid] = float(fields[4])

    # Database profile: taxid <tab> parent <tab> rank <tab> name
    db = dict()
    if db_file:
        for line in gzip.open(db_file, "rt") if db_file.endswith(".gz") else open(db_file, "r"):
            if line[0] == "@":
                continue
            fields = line.rstrip().split("\t")
            taxid = tax.latest(fields[0])
            rank = tax.rank(taxid)
            if taxid == tax.undefined_node:
                sys.stderr.write(db_file + ": " +
                                 fields[0] + " not found in taxonomy\n")
                continue
            if rank not in fixed_ranks:
                sys.stderr.write(db_file + ": " +
                                 fields[0] + " not in defined ranks\n")
                continue
            if rank not in db:
                db[rank] = set()
            db[rank].add(taxid)

    # Results:  TAXID <tab> RANK <tab> TAXPATH <tab> TAXPATHSN <tab> PERCENTAGE
    res = dict()
    for line in gzip.open(res_file, "rt") if res_file.endswith(".gz") else open(res_file, "r"):
        if line[0] == "@":
            continue
        fields = line.rstrip().split("\t")
        taxid = tax.latest(fields[0])
        rank = tax.rank(taxid)
        if taxid == tax.undefined_node:
            sys.stderr.write(res_file + ": " +
                             fields[0] + " not found in taxonomy\n")
            continue
        if rank not in fixed_ranks:
            sys.stderr.write(res_file + ": " +
                             fields[0] + " not in defined ranks\n")
            continue
        if rank not in res:
            res[rank] = {}
        res[rank][taxid] = float(fields[4])

    profile_eval(res, gt, db, fixed_ranks, args.thresholds, output_json)

    if output_json:
        output_json.close()
    

def profile_eval(res, gt, db, fixed_ranks, thresholds, output_json):

    final_stats = {}
    y=defaultdict(list)
    x=defaultdict(list)
    
    fields = ["gt",
              "db",
              "classified",
              "tp",
              "fp",
              "fn",
              "sensitivity_max_db",
              "sensitivity",
              "precision",
              "f1_score",
              "l1n",
              "l2n"]

    for threshold in thresholds:
        db_ranks = defaultdict(int)
        tp_ranks = defaultdict(int)
        fp_ranks = defaultdict(int)
        fn_ranks = defaultdict(int)
        l1_ranks = defaultdict(float)
        l2_ranks = defaultdict(float)
        
        threshold_key = "threshold>=" + str(threshold)
        final_stats[threshold_key] = {f: defaultdict() for f in fields}

        res_thr = {}
        for rank in fixed_ranks:
            res_thr[rank] = {t:a for t, a in res[rank].items() if a >= threshold}

        for rank in fixed_ranks:
            # Keep only above threshold
            res_taxids = set(res_thr[rank].keys())
            gt_taxids = set(gt[rank].keys())
            db_taxids = set(db[rank]) if rank in db else set()

            # intersection is available on db
            db_ranks[rank] = len(db_taxids.intersection(gt_taxids))
            # intersection is true positive
            tp_ranks[rank] = len(res_taxids.intersection(gt_taxids))
            # everything on res not in gt is a false positive
            fp_ranks[rank] = len(res_taxids.difference(gt_taxids))
            # everything on gt not in res is a false negative
            fn_ranks[rank] = len(gt_taxids.difference(res_taxids))
            
            for res_taxid, res_abundance in res_thr[rank].items():
                gt_abundance = gt[rank][res_taxid] if res_taxid in gt[rank] else 0
                d = gt_abundance - res_abundance
                l1_ranks[rank] += abs(d)
                l2_ranks[rank] += pow(d,2)

        header = ["--" + threshold_key + "--"] + fields
                  
        print("\t".join(header), file=sys.stderr)

        for fr in fixed_ranks[::-1]:
            total_taxa_gt = len(gt[fr])
            total_classified_rank = len(res_thr[fr])
            tp = tp_ranks[fr]
            fp = fp_ranks[fr]
            fn = fn_ranks[fr]
            sens = tp/total_taxa_gt
            sens_max = tp/float(db_ranks[fr]) if db_ranks[fr] > 0 else 0
            prec = tp/total_classified_rank if total_classified_rank > 0 else 1
            f1s = (2*sens*prec)/float(sens+prec) if sens+prec > 0 else 0

            print(fr,
                  len(gt[fr]),
                  db_ranks[fr],
                  total_classified_rank,
                  tp,
                  fp,
                  fn,
                  "%.5f" % sens_max,
                  "%.5f" % sens,
                  "%.5f" % prec,
                  "%.5f" % f1s,
                  l1_ranks[fr],
                  l2_ranks[fr],
                  sep="\t", file=sys.stderr)
        
            if output_json:
                final_stats[threshold_key]["db"][fr] = db_ranks[fr]
                final_stats[threshold_key]["gt"][fr] = len(gt[fr])
                final_stats[threshold_key]["classified"][fr] = total_classified_rank
                final_stats[threshold_key]["tp"][fr] = tp
                final_stats[threshold_key]["fp"][fr] = fp
                final_stats[threshold_key]["fn"][fr] = fn
                final_stats[threshold_key]["sensitivity_max_db"][fr] = sens
                final_stats[threshold_key]["sensitivity"][fr] = sens
                final_stats[threshold_key]["precision"][fr] = prec
                final_stats[threshold_key]["f1_score"][fr] = f1s
                final_stats[threshold_key]["l1n"][fr] = l1_ranks[fr]
                final_stats[threshold_key]["l2n"][fr] = l2_ranks[fr]

            y[fr].append(prec)
            x[fr].append(sens)

    final_stats["summary"] = {"aupr": defaultdict()}
    for fr in fixed_ranks[::-1]:
        # add limits to calculate proper aupr
        aupr = np.trapz([1] + y[fr] + [0], x=[0] + x[fr] + [1])
        print("aupr", fr, aupr, sep="\t", file=sys.stderr)
        final_stats["summary"]["aupr"][fr] = aupr

    if output_json:
        json.dump(final_stats, output_json, indent=4)

if __name__ == "__main__":
    main()
