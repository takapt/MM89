#! /usr/bin/env python3

import os
import shutil
import sys
import random
import time
from threading import Lock
from multiprocessing.pool import ThreadPool

HOME = os.environ['HOME']
MM = os.path.join(HOME, 'mm89')
exe = 'a.out' if len(sys.argv) < 2 else sys.argv[1]

exe_path = os.path.join(MM, exe)
copied_exe_path = os.path.join(MM, 'copied_' + str(random.randint(0, 10**5)))
shutil.copy(exe_path, copied_exe_path)

def make_command(seed):
    command = "java -jar tester.jar -exec '{}' -novis -seed {}".format(copied_exe_path, seed)
    return command

def get_score(seed):
    c = make_command(seed)

    start = time.time()
    try:
        output = os.popen(c).read()
    except:
        raise "exit"

    try:
        score = float(output.split()[-1])
    except Exception as e:
        sys.stderr.write('seed: {}'.format(seed) + '\n')
        sys.stderr.write(str(e))
        raise e

    if score == -1:
        score = 0
    exe_time = time.time() - start

    return {'seed': seed, 'score': score, 'time': exe_time}

def single(seeds):
    for seed in seeds:
        result = get_score(seed)
        print('{:5d} {:10.2f} {:.3f}'.format(seed, result['score'], result['time']))
        sys.stdout.flush()

def get_s(result):
    seed = result['seed']
    return '{:10d} {:10.4f} {:.3f}'.format(seed, result['score'], result['time'])



n = 200
# random.seed(810)
# seeds = [random.randint(100, 10**9) for _ in range(n)]
# seeds = list(range(1, n + 1))
seeds = list(range(100, 100 + n + 1))

lock = Lock()
done_seeds = 0
start_time = time.time()
def run(seed):
    r = get_score(seed)

    lock.acquire()

    global done_seeds
    done_seeds += 1
    sys.stderr.write('{:6.2f}% ({:5.4f}): {}'.format(done_seeds / n * 100, (time.time() - start_time) / 60, get_s(r)) + '\n')

#     print(get_s(r))
    sys.stdout.flush()

    lock.release()
    return r

def total(results):
    sum_score = 0
    sum_time = 0
    for r in results:
        sum_score += r['score']
        sum_time += r['time']
    return {'seed': -1, 'score': sum_score / len(results), 'time': sum_time}

def multi(seeds):
    pool = ThreadPool(5)
    results = pool.map(run, seeds)
    for result in results:
        print(get_s(result))
        sys.stdout.flush()

    print(get_s(total(results)))

try:
#     single(range(1, 10))
    multi(seeds)
finally:
    os.remove(copied_exe_path)


