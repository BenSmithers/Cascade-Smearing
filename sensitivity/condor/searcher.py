"""
I need this script to look at the files that failed, figure out which ones that was, and then make a new submit script to re-run the failureis

But not just that! It takes the failed jobs and splits them into two sub-jobs
"""

from glob import glob
import os

from cascade.utils import SterileParams

all_logs = glob(os.path.join("/scratch/bsmithers", "genflux.*.log"))
failed = []

def was_held(filename):
    with open(filename, 'rt') as f:
        if "held" in f.read():
            return True
    return False

def get_fl(job_id, filename):
    with open(filename, 'rt') as f:
        if job_id in f.readline()

def set_line_jobno(line, job_no):
    prefix = "genflux"
    split = line.split()
    split[1] = prefix+"{:06d}".format()
    return " ".join(split)+"\n"


job_count = 0
new_dag_n = "resubmit.dag"
new_dag = open(new_dag_n)

# Now, we use the above function to filter out just the jobs that failed! 
for entry in filter(was_held, all_logs):
    f = open(entry, 'rt')
    f.readline() # skip the first line
    rel = f.readline()
    f.close()
    job_id = rel.split()[2]

    # find the job in the original dag file
    dag_file_name = "/scratch/bsmithers/sensy/genflux.submit"
    og_dag = open(dag_file_name, 'rt')
    full_line = ""
    matches = []
    while True:
        full_line = og_dag.readline()
        if job_id in full_line:
            matches.append(full_line)
        if len(matches)==2:
            break
        if full_line == "":
            break
    if len(matches)!=2:
        raise ValueError("Didn't find right number of matches for jobID {} in {}".format(job_id, dag_file_name))

    # matches should have two entries: one defining the job, and one defining its variables 
    # we want to split the job into two, so for each failed job we split it into two mini-jobs
    for line in matches:
        if "JOB" in line:
            new_dag.write(set_line_jobno(line, job_count))
            new_dag.write(set_line_jobno(line, job_count+1))
        elif "VARS" in line:
            ## ugggg 
            # the line will go like """ [...] 2.4"\n """"
            # so we do some trimming 
            trim = line[:-2]
            # add extra arg to set a little switch to only do half of the mass splittings 
            v1 = trim + " 0\"\n"
            v2 = trim + " 1\"\n"
            
            new_dag.write(set_line_jobno(v1, job_count))
            new_dag.write(set_line_jobno(v2, job_count+1))

        else:
            raise ValueError("Found unfamiliar line: {}".format(line))
    job_count += 2
new_dag.close()
