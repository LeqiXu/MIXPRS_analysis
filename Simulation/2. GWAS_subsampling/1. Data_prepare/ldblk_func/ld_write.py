import numpy as np
import os
import h5py
import sys
import getopt

# parse params
long_opts_list = ['pop=', 'split=', 'geno=', 'fold=', 'help']
param_dict = {'pop': None, 'split': None, 'geno': None, 'fold': None}

if len(sys.argv) > 1:
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
    except:
        print('* Option not recognized.')
        print('* Use --help for usage information.\n')
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h" or opt == "--help":
            print(__doc__)
            sys.exit(0)
        elif opt == "--pop": param_dict['pop'] = arg
        elif opt == "--split": param_dict['split'] = arg
        elif opt == "--geno": param_dict['geno'] = arg
        elif opt == "--fold": param_dict['fold'] = arg
else:
    print(__doc__)
    sys.exit(0)

if param_dict['pop'] == None:
    print('* Please specify the population of the data using --pop\n')
    sys.exit(2)
elif param_dict['split'] == None:
    print('* Please specify the split of the data (discover, validate, discover_validate) using --split\n')
    sys.exit(2)
elif param_dict['geno'] == None:
    print('* Please specify the type of the geno data using --geno\n')
    sys.exit(2)
elif param_dict['fold'] == None:
    print('* Please specify the fold in the cross-validation of the data using --fold\n')
    sys.exit(2)

for key in param_dict:
    print('--%s=%s' % (key, param_dict[key]))

print('\n')

pop = param_dict['pop']
split = param_dict['split']
geno = param_dict['geno']
fold = param_dict['fold']

if geno == "1kg" or geno == "1kg1" or geno == "1kg2":
    geno_type = "1kg"
else:
    geno_type = geno

if fold == "0":
    cv = "all"
else:
    cv = "fold" + fold

# write LD files into h5py file
with open("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/"+pop+"_"+geno_type+"_blk_chr.txt") as ff:
    blk_chr = [int(line.strip()) for line in ff]
with open("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/"+pop+"_"+geno_type+"_blk_size.txt") as ff:
    blk_size = [int(line.strip()) for line in ff]

n_blk = len(blk_chr)
n_chr = max(blk_chr)

for chrom in range(1,n_chr+1):
    print("... parse chomosome %d ..." % chrom)
    if split == "discover_validate":
        hdf_chr = h5py.File("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/"+split+"/"+geno+"/ldblk_"+geno_type+"_"+pop.lower()+"/ldblk_"+geno_type+"_chr"+str(chrom)+".hdf5", "w")
    else:
        hdf_chr = h5py.File("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/"+split+"/"+geno+"_"+cv+"/ldblk_"+geno_type+"_"+pop.lower()+"/ldblk_"+geno_type+"_chr"+str(chrom)+".hdf5", "w")
    blk_cnt = 0
    for blk in range(n_blk):
        if blk_chr[blk] == chrom:
            if blk_size[blk] > 0:
                with open("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk/"+geno_type+"/"+pop+"/ldblk"+str(blk+1)+".ld") as ff:
                    ld = [[float(val) for val in (line.strip()).split()] for line in ff]
                print('blk %d size %s' % (blk+1, np.shape(ld)))
                with open("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_snplist/"+geno_type+"/"+pop+"/"+str(blk+1)) as ff:
                    snplist = [line.strip() for line in ff]
            else:
                ld = []; snplist = []
                print("blk %d size (0, 0)" % (blk+1))

            blk_cnt += 1
            hdf_blk = hdf_chr.create_group("blk_%d" % blk_cnt)
            hdf_blk.create_dataset("ldblk", data=np.array(ld), compression="gzip", compression_opts=9)
            hdf_blk.create_dataset("snplist", data=snplist, compression="gzip", compression_opts=9)
