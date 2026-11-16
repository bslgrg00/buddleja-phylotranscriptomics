"""
restore_bootstrap_support.py
Restores bootstrap support to TreeShrink-trimmed trees: tips removed by TreeShrink
are pruned from the original RAxML bootstrap tree.

Usage: python restore_bootstrap_support.py tree_dir treeshrink_dir [raxml_bs.tre] [tt] [tt_bs.tre]
Output: treeshrink_dir/<cluster>.tt_bs.tre
"""

import newick3, phylo3, os, sys
from tree_utils import remove_kink


def get_tip_labels(root):
    return set([n.label for n in root.iternodes() if n.istip])


def prune_tips(curroot, tips_to_remove):
    """Prune the named tips and remove resulting kinks."""
    going = True
    while going and len(curroot.leaves()) > 2:
        going = False
        for node in curroot.iternodes():
            if node.istip and node.label in tips_to_remove:
                node2 = node.prune()
                if len(curroot.leaves()) >= 2:
                    if (node2 == curroot and node2.nchildren == 2) or \
                       (node2 != curroot and node2.nchildren == 1):
                        node2, curroot = remove_kink(node2, curroot)
                going = True
                break
    return curroot


def main(tree_dir, ts_dir, tree_ending="raxml_bs.tre", tt_ending="tt",
         out_ending="tt_bs.tre"):
    if tree_dir[-1] != "/":
        tree_dir += "/"
    if ts_dir[-1] != "/":
        ts_dir += "/"

    filecount = 0
    n_restored = 0
    n_unchanged = 0
    n_skipped = 0

    for i in os.listdir(ts_dir):
        if not i.endswith("." + tt_ending):
            continue
        clusterID = i[:-(len(tt_ending) + 1)]
        orig_path = tree_dir + clusterID + "." + tree_ending
        tt_path = ts_dir + i
        out_path = ts_dir + clusterID + "." + out_ending

        if not os.path.exists(orig_path):
            print "SKIP:", clusterID, "-- original tree not found:", orig_path
            n_skipped += 1
            continue

        filecount += 1

        with open(orig_path, "r") as f:
            orig_tree = newick3.parse(f.readline())
        with open(tt_path, "r") as f:
            tt_tree = newick3.parse(f.readline())

        orig_tips = get_tip_labels(orig_tree)
        tt_tips = get_tip_labels(tt_tree)
        removed = orig_tips - tt_tips

        if len(removed) == 0:
            # Nothing removed by TreeShrink -- original tree already correct
            with open(out_path, "w") as outfile:
                outfile.write(newick3.tostring(orig_tree) + ";\n")
            n_unchanged += 1
        else:
            pruned = prune_tips(orig_tree, removed)
            with open(out_path, "w") as outfile:
                outfile.write(newick3.tostring(pruned) + ";\n")
            n_restored += 1
            print clusterID, "-- removed", len(removed), "tip(s), bootstrap preserved"

    print ""
    print "========================================"
    print "Total .tt files scanned  :", filecount
    print "Restored (tips removed)  :", n_restored
    print "Unchanged (no removal)   :", n_unchanged
    print "Skipped (no orig tree)   :", n_skipped
    print "========================================"

    assert filecount > 0, "No ." + tt_ending + " files found in " + ts_dir


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage:"
        print "python restore_bootstrap_support.py <tree_dir> <treeshrink_dir> " \
              "[tree_ending=raxml_bs.tre] [tt_ending=tt] [out_ending=tt_bs.tre]"
        sys.exit(0)

    tree_dir = sys.argv[1]
    ts_dir = sys.argv[2]
    tree_ending = sys.argv[3] if len(sys.argv) > 3 else "raxml_bs.tre"
    tt_ending = sys.argv[4] if len(sys.argv) > 4 else "tt"
    out_ending = sys.argv[5] if len(sys.argv) > 5 else "tt_bs.tre"

    main(tree_dir, ts_dir, tree_ending, tt_ending, out_ending)
