#!/usr/bin/env python3

"""
Submit this clustering script for sbatch to snakemake with:
    snakemake --rerun-incomplete -j 99 --debug --immediate-submit --cluster-config cluster.json --cluster 'grun_cluster.py {dependencies}'
"""

## In order to submit all the jobs to the moab queuing system, one needs to write a wrapper.
## This wrapper is inspired by Daniel Park https://github.com/broadinstitute/viral-ngs/blob/master/pipes/Broad_LSF/cluster-submitter.py
## I asked him questions on the snakemake google group and he kindly answered: https://groups.google.com/forum/#!topic/snakemake/1QelazgzilY

import sys
import re
import os
import errno
from snakemake.utils import read_job_properties

## snakemake will generate a jobscript containing all the (shell) commands from your Snakefile. 
## I think that's something baked into snakemake's code itself. It passes the jobscript as the last parameter.
## https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-job-properties

## make a directory for the logs from the cluster 
try: 
	os.makedirs("grun_log")
except OSError as exception:
	if exception.errno != errno.EEXIST:
		raise


jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

## the jobscript is something like snakejob.index_bam.23.sh
mo = re.match(r'(\S+)/snakejob\.\S+\.(\d+)\.sh', jobscript)
assert mo
sm_tmpdir, sm_jobid = mo.groups()

## set up jobname. 
jobname = "{rule}-{jobid}".format(rule = job_properties["rule"], jobid = sm_jobid)

## it is safer to use get method in case the key is not present
# the job_properties is a dictionary of dictonary. I set up job name in the Snake file under the params directive and associate the sample name with the 
# job

jobname_tag_sample = job_properties.get('params', {}). get('jobname')


if jobname_tag_sample:
	jobname = jobname + "-" + jobname_tag_sample


# access property defined in the cluster configuration file (Snakemake >=3.6.0), cluster.json
time = job_properties["cluster"]["time"]
mem = job_properties["cluster"]["MaxMem"]
queue = job_properties["cluster"]["queue"]
EmailNotice = job_properties["cluster"]["EmailNotice"]
email = job_properties["cluster"]["email"]

cmdline = 'grun -cput {time} -M {email} -q {queue} -n {jobname} -o grun_log/{out}.out -e grun_log/{err}.err'.format(time = time, email = email, queue = queue, jobname = jobname, out = jobname, err = jobname)
cmdline += ' -M {} '
#cmdline += ' -M {} -R rusage[mem={}] '.format(int(mem), int(mem))

# figure out job dependencies, the last argument is the jobscript which is baked in snakemake
# man bsub to see -w documentation
# exit(job_ID | "job_name"
#  The job state is EXIT, and the job's exit code
#  satisfies the comparison test.

#  If you specify an exit code with no operator,
#   the test is for equality (== is assumed).

#   If you specify only the job, any exit code
#   satisfies the test.

#Enclose the dependency expression in single quotes (')
#to prevent the shell from interpreting special
#                characters (space, any logic operator, or parentheses).
#                If you use single quotes for the dependency expression,
#                use double quotes (") for quoted items within it, such
#                as job names.

# e.g.
#  dependencies = {'1234', '2345', '122'}
#
# " -w '{}' ".format(" && ".join(["done({})".format(dependency) for dependency in dependencies]) )
# " -w 'done(1234) && done(2345) && done(122)' 

# -w done(1234) is the same as -w 1234

# single quote is important '{}'
dependencies = set(sys.argv[1:-1])
if dependencies:
	cmdline += " -w '{}' ".format(" && ".join(["done({})".format(dependency) for dependency in dependencies]))

# the actual job
cmdline += jobscript

# the part that strips bsub's output to just the job id
cmdline += r" | tail -1 | cut -f 2 -d \< | cut -f 1 -d \>"

# call the command
os.system(cmdline)
