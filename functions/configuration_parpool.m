
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
mkdir(scratch_dir);
pc = parcluster('local');
pc.JobStorageLocation = scratch_dir;
parpool(pc,pc.NumWorkers);
