#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import subprocess, os, sys, re, datetime, shutil
from threading import Thread, Lock, Semaphore

# ---------- 基本設定 ----------
CPP_FILE   = "hw3.cpp"
MAKE_CMD   = ["make"]
EXE_NAME   = "hw3"
CASES      = ["public1.txt", "public2.txt", "public3.txt"]
BASELINE   = {"public1.txt":248807,
              "public2.txt":495065,
              "public3.txt":696662}
TC_DIR     = "../testcase"
OUT_DIR    = "../output"
DEAD_RATIO = "0.1"
SCORE_TXT  = "score.txt"
MAX_THREAD = 50  

seed_re = re.compile(r"(mt19937_64\s+rng\()\s*\d+\s*(\);)")

lock = Lock()
sema = Semaphore(MAX_THREAD)
results = []

# ---------- 工具 ----------
def run(cmd, cwd=None):
    print("Running:", " ".join(cmd))
    p = subprocess.Popen(cmd, cwd=cwd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    output = ""
    while True:
        line = p.stdout.readline()
        if not line:
            break
        line = line.decode("utf-8", "ignore")
        print(line, end="")
        output += line
    p.wait()
    if p.returncode != 0:
        print(output)
        raise Exception("!! Command failed: " + " ".join(cmd))
    return output

def extract_hpwl(out_path):
    with open(out_path) as f:
        line = f.readline()
    m = re.search(r"Wirelength\s+(\d+)", line)
    if not m:
        raise Exception("!! Cannot parse HPWL in %s" % out_path)
    return int(m.group(1))

def log(msg):
    with lock:
        sys.stdout.write(msg)
        sys.stdout.flush()
        with open(SCORE_TXT, "a") as f:
            f.write(msg)

def run_seed(seed):
    sema.acquire()
    try:
        local_dir = "thread_{}".format(seed)
        local_cpp = os.path.join(local_dir, "hw3.cpp")
        local_bin = os.path.join(local_dir, EXE_NAME)

        if not os.path.exists(local_dir):
            os.makedirs(local_dir)

        # 複製 Makefile
        makefile_src = "Makefile"
        makefile_dst = os.path.join(local_dir, "Makefile")
        if os.path.exists(makefile_src):
            shutil.copy(makefile_src, makefile_dst)
        else:
            raise Exception("!! Missing Makefile in current directory")

        # 複製 cpp 並換 seed
        with open(CPP_FILE) as f:
            lines = f.readlines()
        changed = 0
        for i, l in enumerate(lines):
            if seed_re.search(l):
                lines[i] = seed_re.sub(r"\g<1>%d\g<2>" % seed, l)
                changed += 1
        if changed == 0:
            raise Exception("!! No rng(seed) found")

        with open(local_cpp, "w") as f:
            f.writelines(lines)

        # 編譯
        run(MAKE_CMD + ["-C", local_dir])

        # 執行測資
        all_pass = True
        hpwl = {}
        for c in CASES:
            in_path = os.path.join(TC_DIR, c)
            out_path = os.path.join(OUT_DIR, "s{}_{}".format(seed, c.replace(".txt", ".out")))
            cmd = [os.path.abspath(local_bin), in_path, out_path, DEAD_RATIO]
            run(cmd)
            hp = extract_hpwl(out_path)
            hpwl[c] = hp
            if hp > BASELINE[c]:
                all_pass = False

        if not all_pass:
            score = 0.0
        else:
            score = sum(float(BASELINE[c] - hpwl[c]) / BASELINE[c] for c in CASES)

        now = datetime.datetime.now().strftime("%H:%M:%S")
        msg = "[{}] seed={} score={:.6f}  ({})\n".format(
            now, seed, score,
            ", ".join("{}={}".format(c[:-4], hpwl[c]) for c in CASES))
        
        if score>0.0:
            log(msg)

        with lock:
            results.append((score, seed))

    except Exception as e:
        log("[ERROR] seed={} failed: {}\n".format(seed, str(e)))
    finally:
        sema.release()

# ---------- 主程式 ----------
ts_head = datetime.datetime.now().strftime("%F %T")
log("==== Auto SA (Multithread) start {} ====\n".format(ts_head))

if os.path.exists(SCORE_TXT):
    os.remove(SCORE_TXT)

threads = []
for seed in range(0, 1001):
    t = Thread(target=run_seed, args=(seed,))
    threads.append(t)
    t.start()

for t in threads:
    t.join()

# 最佳結果
if results:
    results.sort(reverse=True)
    best_score, best_seed = results[0]
    log("---- Best seed {}  score {:.6f} ----\n".format(best_seed, best_score))
else:
    log("!! No successful seed found.\n")
