#!/usr/bin/env python3

import sys
import gzip
import json
import argparse
from pylca.pylca import LCA
from collections import defaultdict

from util import default_ranks, parse_tax, file_exists


def main():

    parser = argparse.ArgumentParser(
        prog="ganon benchmark evaluation", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input-results", metavar="", type=file_exists,
                        help="readid <tab> assembly id (or 0) <tab> taxid")
    parser.add_argument("-g", "--input-ground-truth", metavar="",
                        type=file_exists, help="readid <tab> assembly id (or 0) <tab> taxid")
    parser.add_argument("-d", "--input-database-profile", metavar="",
                        type=file_exists, help="taxid <tab> parent <tab> rank <tab> name")
    parser.add_argument("-t", "--taxonomy",       metavar="",
                        type=str, help="custom, ncbi or gtdb")
    parser.add_argument("-x", "--taxonomy-files", metavar="",
                        type=file_exists, nargs="*", help="")
    parser.add_argument("-r", "--ranks",          metavar="", type=str, nargs="*",
                        help="empty for default ranks (superkingdom phylum class order family genus species)")
    parser.add_argument("-c", "--output-cumulative", type=str,
                        help="Output file for evaluation in json (cumulative mode)")
    parser.add_argument("-r", "--output-rank", type=str,
                        help="Output file for evaluation in json (rank mode)")
    args = parser.parse_args()

    tax = parse_tax(args.taxonomy, args.taxonomy_files)
    fixed_ranks = ["root"] + args.ranks if args.ranks else default_ranks

    gt_file = args.input_ground_truth
    db_file = args.input_database_profile
    res_file = args.input_results

    output_cumu = open(args.output_cumulative,
                       "w") if args.output_cumulative else None
    output_rank = open(args.output_rank, "w") if args.output_rank else None

    # Ground truth: readid <tab> assembly <tab> taxid
    gt = dict()
    gt_leaf_taxids = set()
    for line in gzip.open(gt_file, "rt") if gt_file.endswith(".gz") else open(gt_file, "r"):
        if line[0] == "@":
            continue
        fields = line.rstrip().split("\t")
        taxid = tax.latest(fields[2])
        if taxid == tax.undefined_node:
            sys.stderr.write(gt_file + ": " +
                             fields[2] + " not found in taxonomy\n")
            continue
        readid = fields[0]
        assembly = fields[1]
        gt_leaf_taxids.add(taxid)
        gt[readid] = (assembly, taxid)

    # Database profile: taxid <tab> parent <tab> rank <tab> name
    db_taxids = set()
    db_assembly = set()
    if db_file:
        for line in gzip.open(db_file, "rt") if db_file.endswith(".gz") else open(db_file, "r"):
            if line[0] == "@":
                continue
            fields = line.rstrip().split("\t")
            taxid = tax.latest(fields[0])
            if taxid == tax.undefined_node:
                if fields[2] == "assembly":
                    db_assembly.add(fields[0])
                else:
                    sys.stderr.write(db_file + ": " +
                                     fields[0] + " not found in taxonomy\n")
                    continue
            else:
                db_taxids.add(taxid)

    # Results: readid <tab> assembly <tab> taxid
    # if no assembly readid <tab> 0 <tab> taxid
    res = dict()
    res_leaf_taxids = set()
    for line in gzip.open(res_file, "rt") if res_file.endswith(".gz") else open(res_file, "r"):
        if line[0] == "@":
            continue
        fields = line.rstrip().split("\t")
        taxid = tax.latest(fields[2])
        if taxid == tax.undefined_node:
            sys.stderr.write(res_file + ": " +
                             fields[2] + " not found in taxonomy\n")
            continue
        readid = fields[0]
        assembly = fields[1]
        res_leaf_taxids.add(taxid)
        res[readid] = (assembly, taxid)

    # Cumulative evaluation
    # filter tax for faster LCA (keep root out)
    valid_taxids = gt_leaf_taxids | res_leaf_taxids
    valid_taxids.discard(tax.root_node)
    tax.filter(valid_taxids)
    L = LCA(tax._nodes)

    # pre-build lineages for faster access
    tax.build_lineages(ranks=fixed_ranks)

    # check lineage of the gt taxids to check at which rank it could be classified given the database
    rank_gttaxid = {}
    for leaf_gttaxid in gt_leaf_taxids:
        if not rank_gttaxid.get(leaf_gttaxid):  # if not yet calculated
            for t in tax.lineage(leaf_gttaxid)[::-1]:
                if t in db_taxids:
                    rank_gttaxid[leaf_gttaxid] = tax.rank(
                        tax.closest_parent(t, fixed_ranks))
                    break

    # 1 - check if there is an accession (uniq assignment) and if it"s correct
    # 2 - check if taxid matches
    #   if taxid tool = taxid gt -> correct identification at taxonomic level
    #   if lca = tool results -> TP, meaning it got it right in a lower taxonomic level, read was more specific
    #   if lca = gt -> FP, meaning the classification was too specific
    cumulative_eval(res, gt, tax, fixed_ranks, L, rank_gttaxid, output_cumu)

    # Rank evaluation
    rank_eval(res, gt, tax, fixed_ranks, db_assembly, db_taxids, output_rank)

    if output_cumu:
        output_cumu.close()
    if output_rank:
        output_rank.close()


def print_log(text):
    sys.stderr.write(text+"\n")


def cumulative_eval(res, gt, tax, fixed_ranks, L, rank_gttaxid, output_cumu):

    stats = {"classified": 0, "unclassified": 0, "tp": 0, "fp": 0}
    classified_ranks = defaultdict(int)
    db_ranks = defaultdict(int)
    gt_ranks = defaultdict(int)
    tp_direct_ranks = defaultdict(int)
    tp_indirect_ranks = defaultdict(int)
    fp_lower_ranks = defaultdict(int)
    fp_higher_ranks = defaultdict(int)

    for readid, (_, gt_taxid) in gt.items():

        leaf_rank_gt = tax.rank(tax.closest_parent(gt_taxid, fixed_ranks))
        gt_ranks[leaf_rank_gt] += 1

        # if rank level taxid is present in the database (=could be classified)
        # define at which level a certain taxid could be classified given this db
        if gt_taxid in rank_gttaxid:
            db_ranks[rank_gttaxid[gt_taxid]] += 1

        if readid in res.keys():  # if read is classified
            res_taxid = res[readid][1]

            # taxonomic clasification
            r = tax.rank(tax.closest_parent(res_taxid, fixed_ranks))

            classified_ranks[r] += 1
            if r == tax.root_rank:  # root classification is equal to false
                fp_lower_ranks[r] += 1
            else:
                if res_taxid == gt_taxid:  # tp -> perfect classification
                    tp_direct_ranks[r] += 1
                else:
                    lca = L(gt_taxid, res_taxid)
                    # tp -> conservative classification (gt is lower on tree)
                    if lca == res_taxid:
                        tp_indirect_ranks[r] += 1
                    # fp -> classification to specific (gt is higher on tree)
                    elif lca == gt_taxid:
                        fp_higher_ranks[r] += 1
                    else:  # fp -> lca is higher than gt and res
                        fp_lower_ranks[r] += 1
        else:
            stats["unclassified"] += 1

    total_reads_gt = len(gt)
    stats["classified"] = total_reads_gt - stats["unclassified"]
    stats["tp"] = sum(tp_direct_ranks.values()) + \
        sum(tp_indirect_ranks.values())
    stats["fp"] = stats["classified"] - stats["tp"]

    final_stats = defaultdict(dict)
    header = ["--cumu_eval--",
              "gt",
              "db",
              "class",
              "tp",
              "fp",
              "sens_max_db",
              "sens",
              "prec",
              "f1s",
              "cs_db",
              "cs_class",
              "cs_tp",
              "cs_fp",
              "tp_direct",
              "tp_indirect",
              "fp_lower",
              "fp_higher"]

    # taxonomic stats
    print("\t".join(header), file=sys.stderr)
    cs_db = 0
    cs_class = 0
    cs_tp = 0
    cs_fp = 0

    for fr in fixed_ranks[::-1]:
        tp = tp_direct_ranks[fr] + tp_indirect_ranks[fr]
        fp = fp_lower_ranks[fr] + fp_higher_ranks[fr]

        cs_db += db_ranks[fr]  # make it cumulative
        cs_class += classified_ranks[fr]
        cs_tp += tp
        cs_fp += fp

        # if root, all available
        if fr == "root":
            cs_db = total_reads_gt

        sens_max_db = cs_tp/float(cs_db) if cs_db > 0 else 0
        sens = cs_tp/total_reads_gt
        prec = cs_tp/float(cs_class) if cs_class > 0 else 0
        f1s = (2*sens*prec)/float(sens+prec) if sens+prec > 0 else 0

        print(fr,
              gt_ranks[fr],
              db_ranks[fr],
              classified_ranks[fr],
              tp,
              fp,
              "%.5f" % sens_max_db,
              "%.5f" % sens,
              "%.5f" % prec,
              "%.5f" % f1s,
              cs_db,
              cs_class,
              cs_tp,
              cs_fp,
              tp_direct_ranks[fr],
              tp_indirect_ranks[fr],
              fp_lower_ranks[fr],
              fp_higher_ranks[fr], sep="\t", file=sys.stderr)

        if output_cumu:
            final_stats["gt"][fr] = gt_ranks[fr]
            final_stats["db"][fr] = db_ranks[fr]
            final_stats["classified"][fr] = classified_ranks[fr]
            final_stats["tp"][fr] = tp
            final_stats["fp"][fr] = fp
            final_stats["sensitivity_max_db"][fr] = sens_max_db
            final_stats["sensitivity"][fr] = sens
            final_stats["precision"][fr] = prec
            final_stats["f1_score"][fr] = f1s
            final_stats["cs_db"][fr] = cs_db
            final_stats["cs_classified"][fr] = cs_class
            final_stats["tp_direct"][fr] = tp_direct_ranks[fr]
            final_stats["tp_indirect"][fr] = tp_indirect_ranks[fr]
            final_stats["fp_lower"][fr] = fp_lower_ranks[fr]
            final_stats["fp_higher"][fr] = fp_higher_ranks[fr]

    if output_cumu:
        json.dump(final_stats, output_cumu, indent=4)


def rank_eval(res, gt, tax, fixed_ranks, db_assembly, db_taxids, output_rank):

    stats = {"classified": 0, "unclassified": 0, "tp": 0, "fp": 0}
    classified_ranks = defaultdict(int)
    classified_ranks_assembly = 0
    db_ranks = defaultdict(int)
    db_ranks_assembly = 0
    gt_ranks = defaultdict(int)
    gt_ranks_assembly = 0
    tp_ranks = defaultdict(int)
    tp_ranks_assembly = 0
    fp_ranks = defaultdict(int)
    fp_ranks_assembly = 0

    for readid, (gt_assembly, gt_taxid) in gt.items():

        # Check if there"s assembly id on ground truth and account for it
        if gt_assembly != "0":
            gt_ranks_assembly += 1
            # if assembly is present in the database (=could be classified)
            if gt_assembly in db_assembly:
                db_ranks_assembly += 1

        # already pre-built build_lineages() with fixed_ranks
        gt_lin = tax.lineage(gt_taxid)
        # Check if the taxid is on ground truth and account for it
        for idx, fr in enumerate(fixed_ranks):
            if gt_lin[idx]:
                gt_ranks[fr] += 1
                if gt_lin[idx] in db_taxids:
                    # if the is present in the database (=could be classified)
                    db_ranks[fr] += 1

        if readid in res.keys():  # if read is classified
            res_assembly = res[readid][0]
            res_taxid = res[readid][1]

            # has a unique assembly classification
            if res_assembly != "0":
                classified_ranks_assembly += 1
                if res_assembly == gt_assembly:  # it is correct
                    tp_ranks_assembly += 1
                else:
                    fp_ranks_assembly += 1

            # already pre-built build_lineages() with fixed_ranks
            res_lin = tax.lineage(res_taxid)
            # compare every taxonomic rank
            for idx, fr in enumerate(fixed_ranks):
                if res_lin[idx]:  # if there"s a classification for such rank
                    classified_ranks[fr] += 1
                    if gt_lin[idx] == res_lin[idx]:
                        tp_ranks[fr] += 1
                    else:
                        fp_ranks[fr] += 1
        else:
            stats["unclassified"] += 1

    total_reads_gt = len(gt)
    stats["classified"] = total_reads_gt - stats["unclassified"]
    classified_ranks["root"] = stats["classified"]  # all have root

    final_stats = defaultdict(dict)
    header = ["--rank_eval--",
              "gt",
              "db",
              "class",
              "tp",
              "fp",
              "sens_max_db",
              "sens",
              "prec",
              "f1s"]
    print("\t".join(header), file=sys.stderr)

    sens_assembly = tp_ranks_assembly/total_reads_gt
    sens_max_assembly = tp_ranks_assembly / \
        float(db_ranks_assembly) if db_ranks_assembly > 0 else 0
    prec_assembly = tp_ranks_assembly / \
        classified_ranks_assembly if classified_ranks_assembly > 0 else 0
    f1s_assembly = (2*sens_assembly*prec_assembly)/float(sens_assembly +
                                                         prec_assembly) if sens_assembly + prec_assembly > 0 else 0
    print("assembly",
          gt_ranks_assembly,
          db_ranks_assembly,
          classified_ranks_assembly,
          tp_ranks_assembly,
          fp_ranks_assembly,
          "%.5f" % sens_max_assembly,
          "%.5f" % sens_assembly,
          "%.5f" % prec_assembly,
          "%.5f" % f1s_assembly, sep="\t", file=sys.stderr)

    if output_rank:
        final_stats["db"]["assembly"] = db_ranks_assembly
        final_stats["gt"]["assembly"] = gt_ranks_assembly
        final_stats["classified"]["assembly"] = classified_ranks_assembly
        final_stats["tp"]["assembly"] = tp_ranks_assembly
        final_stats["fp"]["assembly"] = fp_ranks_assembly
        final_stats["sensitivity"]["assembly"] = sens_assembly
        final_stats["precision"]["assembly"] = prec_assembly
        final_stats["f1_score"]["assembly"] = f1s_assembly

    for fr in fixed_ranks[::-1]:
        # if root, all classified are true
        tp = tp_ranks[fr] if fr != "root" else classified_ranks[fr]
        fp = fp_ranks[fr]
        sens = tp/total_reads_gt
        sens_max = tp/float(db_ranks[fr]) if db_ranks[fr] > 0 else 0
        prec = tp/classified_ranks[fr] if classified_ranks[fr] > 0 else 0
        f1s = (2*sens*prec)/float(sens+prec) if sens+prec > 0 else 0

        print(fr, gt_ranks[fr], db_ranks[fr], classified_ranks[fr], tp, fp, "%.5f" %
              sens_max, "%.5f" % sens, "%.5f" % prec, "%.5f" % f1s, sep="\t", file=sys.stderr)

        if output_rank:
            final_stats["db"][fr] = db_ranks[fr]
            final_stats["gt"][fr] = gt_ranks[fr]
            final_stats["classified"][fr] = classified_ranks[fr]
            final_stats["tp"][fr] = tp
            final_stats["fp"][fr] = fp
            final_stats["sensitivity"][fr] = sens
            final_stats["precision"][fr] = prec
            final_stats["f1_score"][fr] = f1s

    if output_rank:
        json.dump(final_stats, output_rank, indent=4)


if __name__ == "__main__":
    main()
