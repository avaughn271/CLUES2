import cli
def main():
    cli.arginfer_main()

if __name__ == "__main__":
    main()

'''
simulation: 
python3 __main__.py infer -I 1000 --thin 20 --burn 0 -n 5 -L 1e3  --Ne 5000 -r 1e-8 -mu 1e-8 \
        --tsfull /Users/amahmoudi/Ali/phd/github_projects/mcmc/test1/ts_sim/sim_r1/n5Ne5K_L1K_iter0.args \
        -O /Users/amahmoudi/Ali/phd/github_projects/mcmc/arginfer/output \
        --random-seed 5 -p --verbose --verify        
real data:
python3 __main__.py infer -I 30 --thin 0 --burn 0 -n 10 -L 1e5  --Ne 5000 -r 1e-8 -mu 1e-8 \
        --input_path /Users/amahmoudi/Ali/phd/github_projects/real_data/qced_data \
        --haplotype_name "haplotype_ready.txt" \
        --ancAllele_name "ancestral_allele_ready.txt" \
        --snpPos_name "SNP_pos_ready.txt" \
        -O /Users/amahmoudi/Ali/phd/github_projects/mcmc/ARGinfer/output \
        --random-seed 5 -p -v --verify
'''
