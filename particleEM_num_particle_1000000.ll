# @ job_type = bluegene
# @ class = normal
#
# @ job_name = particleEM_num_particle_1000000
#
# The executable that will run your parallel application should always be specified as per the next line.
# @ executable = /usr/bin/mpirun
#
#
# Run job on 512 compute nodes. LoadLeveler will dynamically allocate the partition.
# @ bg_size = 128
#
# initialdir will be the initial directory. LoadLeveler changes to this
# directory before running the job. If you don't specify this, it defaults to your current working directory
# at the time the job was submitted.
# File names mentioned in the batch script that do not begin with a slash ( / ) are relative to the initial
# directory.
# The initial directory must exist on both the fen and the compute nodes.
# @ initialdir = /gpfs/home3/k/kyu/projects/EM/original
#
# If for example your jobid is 82, your output and error will be written in
# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.
# @ input = /dev/null
# @ output = $(job_name).$(jobid).out
# @ error = $(job_name).$(jobid).err
# 
# Maximum wall clock time for job will be 5 minutes.
# @ wall_clock_limit = 00:20:00
#
# Send email to yukwangmin@gmail.com when job has completed.
# @ notification = complete
# @ notify_user = yukwangmin@gmail.com
#
# Specify executable for your parallel application, and arguments to that executable.
# Note that the arguments to specify for the executable will vary depending upon the executable.
# Specify any special environment variables for your application that need to be in the environment 
# presented to the job on the compute nodes, they will vary
# depending upon the application -  some applications will not require any - so delete or modify the
# -env specification below.
# @ arguments =  -exe $(initialdir)/BGPEM \ 
-cwd $(initialdir) \
-mode VN\
#-env "OMP_NUM_THREADS=4" \
-args " --procNumX=8 --procNumY=8 --procNumZ=8 --setting=setting_MPI.xml" 
#
# The next statement marks the end of the job step. This example is a one-job-step batch job,
# so this is equivalent to saying that the next statement marks the end of the batch job.
# @ queue

