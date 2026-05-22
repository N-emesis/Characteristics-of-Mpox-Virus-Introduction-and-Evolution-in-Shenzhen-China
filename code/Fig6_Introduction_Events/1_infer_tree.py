import subprocess
import shlex
import argparse

def infer_phylogeny(aln, outdir):
    iqtree_cmd = f'/root/miniforge3/envs/bio/bin/iqtree -redo -nt AUTO -ninit 10 -me 0.05 -bb 1000 -wbtl -czb -m GTR --prefix {outdir}/aligned -s {aln}'
    print("Running:", iqtree_cmd)
    subprocess.run(shlex.split(iqtree_cmd))
    return(f'{outdir}/aligned.treefile', f'{outdir}/aligned.ufboot')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--aln_file', default='data/aligned.fasta')
    parser.add_argument('--outdir', default='results')
    args = parser.parse_args()
    
    subprocess.run(shlex.split(f'mkdir -p {args.outdir}'))
    infer_phylogeny(args.aln_file, args.outdir)
    print('Tree inferred.')

if __name__ == "__main__":
    main()
