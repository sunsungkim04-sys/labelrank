"""Stand-in for `hyphy CPU=1 RELAX_seeded.bf --alignment A --tree T --test X --reference R --output J`.

Reads the tree it was given, prints the branch sets the way RELAX does, and writes a relax.json whose LRT is a
fixed function of the test tips: LRT = 9 if the test set holds a tip listed in FAKE_HIGH (comma-separated
environment variable), else a small number from the tip names. No statistics are involved; this exercises the
pipeline (setup -> run -> report) end to end.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from labelrank import newick  # noqa: E402


def main(argv):
    if argv[:1] == ["--version"]:
        print("HYPHY %s(MP) stand-in for tests" % os.environ.get("FAKE_VERSION", "2.5.101"))
        return
    a = {}
    i = 0
    while i < len(argv):
        if argv[i].startswith("--"):
            a[argv[i]] = argv[i + 1]
            i += 2
        else:
            i += 1
    if os.environ.get("FAKE_COUNT"):   # one line per engine start, to count double launches
        with open(os.environ["FAKE_COUNT"], "a") as fh:
            fh.write(os.getcwd() + "\n")
    if os.environ.get("FAKE_SLEEP"):   # a slow engine, for interrupt and timeout tests
        import time
        time.sleep(float(os.environ["FAKE_SLEEP"]))
    with open(a["--tree"]) as fh:
        root = newick.parse(fh.read())
    br = newick.branches(root)
    test = sorted(n.name for n in br if n.label == a["--test"])
    if a["--reference"] == "Unlabeled branches":
        ref = sorted(n.name for n in br if not n.label)
    else:
        ref = sorted(n.name for n in br if n.label == a["--reference"])
    nui = sorted(n.name for n in br if n.name not in set(test) | set(ref))
    print("* Selected %d branches as the _Reference_ set: `%s`" % (len(ref), ", ".join(ref)))
    print("* Selected %d branches as the _Test_ set: `%s`" % (len(test), ", ".join(test)))
    if nui:
        print("* %d branches are in the unclassified (nuisance) set: `%s`" % (len(nui), ", ".join(nui)))
    high = set(x for x in os.environ.get("FAKE_HIGH", "").split(",") if x)
    fixed = dict(x.split("=") for x in os.environ.get("FAKE_LRT", "").split(",") if x)   # tip=LRT, *=default
    if fixed:
        lrt = float(next((fixed[x] for x in test if x in fixed), fixed.get("*", "0")))
    elif high & set(test):
        lrt = 9.0
    else:
        lrt = (sum(ord(c) for c in "".join(test)) % 17) / 10.0
    alt, null, k = -1000.0 + lrt / 2, -1000.0, 0.5
    print("### Fitting the alternative model to test K != 1")
    print("* Log(L) = %.2f, AIC-c = 1.00 (2 estimated parameters)" % alt)
    print("* Relaxation/intensification parameter (K) =     %.2f" % k)
    print("### Fitting the null (K := 1) model")
    print("* Log(L) = %.2f, AIC-c = 1.00 (1 estimated parameters)" % null)
    import math
    with open(a["--output"], "w") as fh:
        fh.write(json.dumps({"fits": {"RELAX alternative": {"Log Likelihood": alt},
                                      "RELAX null": {"Log Likelihood": null}},
                             "test results": {"LRT": lrt, "p-value": math.erfc(math.sqrt(max(lrt, 0) / 2)),
                                              "relaxation or intensification parameter": k},
                             "branch attributes": {"x": "__INF__"},   # HyPhy writes bare inf
                             "analysis": {"settings": {}}}).replace('"__INF__"', "inf"))

if __name__ == "__main__":
    main(sys.argv[1:])
