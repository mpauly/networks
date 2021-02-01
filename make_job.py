import argparse
import os
import random
import string
import subprocess
import tempfile


def random_id(length=8):
    return "".join(random.sample(string.ascii_letters + string.digits, length))


BASEDIR = "/remote/pi312a/pauly/networks"
TEMPLATE_SERIAL = """
#!/bin/bash
#PBS -l nodes=1:ppn=4:medium_buster
#PBS -q medium_buster
#PBS -l mem=8gb,vmem=8gb
#PBS -l walltime=72:00:00
#PBS -N {name}_{batch_id}
#PBS -o {base_dir}/logs/{name}_{batch_id}.run.log
#PBS -e {base_dir}/logs/{name}_{batch_id}.err.log
OMP_NUM_THREADS=$PBS_NUM_PPN
export OMP_NUM_THREADS
cd {base_dir}
log_file={base_dir}/logs/{name}_{batch_id}.log
echo "------------------------------------------------------------------------" >> $log_file
echo "Job started on" `date`  >> $log_file
echo "------------------------------------------------------------------------"  >> $log_file
nice -19 ./queue_walk.x -i {interval} -W {start_walk} -w {batch_size} -l {walk_length} graphs/{name}.dat  >> $log_file
echo "------------------------------------------------------------------------"  >> $log_file
echo "Job ended on" `date`  >> $log_file
echo "------------------------------------------------------------------------"  >> $log_file
"""


parser = argparse.ArgumentParser(
    description="A helper script that allows to qsub a number of jobs"
)
parser.add_argument(
    "-i", "--interval", help="Interval for intermediate saves", type=int
)
parser.add_argument("graph", help="the graph to walk on in format roadnet_pa")
parser.add_argument("total_walkers", help="the total number of walkers", type=int)
parser.add_argument("batch_size", help="the batch size", type=int)
parser.add_argument("walk_length", help="the number of steps to take", type=int)
args = parser.parse_args()

batch_size = args.batch_size
total_walkers = args.total_walkers

name = args.graph

if batch_size > total_walkers:
    print("Error: Cannot have larger batches than total walkers")
    quit(1)

for batch_id in range(total_walkers // batch_size):
    start_walk = batch_size * batch_id
    base = "submit_{0}".format(random_id())
    open(base + ".qsub", "w").write(
        TEMPLATE_SERIAL.format(
            name=name,
            base_dir=BASEDIR,
            batch_id=batch_id,
            walk_length=args.walk_length,
            interval=args.interval,
            start_walk=start_walk,
            batch_size=batch_size,
        )
    )
    try:
        # subprocess.call("qsub " + base + ".qsub", shell=True)
        subprocess.call("cat " + base + ".qsub", shell=True)
    finally:
        os.remove(base + ".qsub")
