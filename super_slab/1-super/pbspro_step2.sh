#! /bin/bash                                                                                           
#PBS -N SuperTBH
#PBS -q workq
##---------------------------------------------------------------------------------------------------
## -l select=<# of chunks>[:<ncpus=#>], how many chunks. ncpus, how many cpus for each chunk.
## chunks is the set of resources allocated, typically, there is one chunk per MPI process.
## recommend select=2:ncpus=1 for 2-rank MPI job, select=1:ncpus=2 for 2-threads openMP job.
## select=<...>:n1=6130 or n=E5 or n=E5fat for request the 6130 CPU nodes, E5 CPU nodes, or E5fat node
##---------------------------------------------------------------------------------------------------
#PBS -l select=1:ncpus=4
##---------------------------------------------------------------------------------------------------
## -l place= <type>[:<sharing>][:<group>]
## place the chunks on which way.
## vnodes is a virtual node, a abstract object representing a sub set of resources of a machine(host).
###     free,     on any vnode(s)
###     pack,     all chunks will be taken from one host
###     scatter,  only one chunk is taken from any host(means try to distribute your MPI rank to nodes)
###     vscatter, only one chunk is taken from any vnode(same as above but to vnodes)
## recommned place=pack for avoiding the hosts-crossing jobs
## place=pack:excl(exclhost) for exclusively using the vnode(node)
## place=<...>:group=n or group=n1 , use place chunks to n or n1 group, will use others if ocuupied
##---------------------------------------------------------------------------------------------------
#PBS -l place=free
##---------------------------------------------------------------------------------------------------

if ! $(type module &>/dev/null) || [[ -z "$LMOD_CMD" ]]; then
    echo "The \`module' or lua mod command failed, pleae contact the administrator"
    exit
fi

module load matlab

cd $PBS_O_WORKDIR
echo -e "you are acquiring the resoures\n$(cat $PBS_NODEFILE)"

run="matlab -nodisplay -r main < SuperTBHmftn.m > out"  # type waht you mainly want to run here
eval $run 1> >(tee ${PBS_JOBID}.${PBS_JOBNAME}.out) 2> >(tee -a ${PBS_JOBID}.${PBS_JOBNAME}.out)
