python prepare_workflow.py config/config.yaml
snakemake --configfile config/config.yaml \
--use-singularity --singularity-args  \
"--nv -B $(pwd)/output:/root/af_output -B $(pwd)/workflow/scripts/match_template_and_target_chain_ids.py:/app/utils/chain_pair_generator.py -B $(realpath ../../../../../af2db/AlphaFold):/AF_data" \
-j 300 -c all \
-p -k -w 30 --rerun-triggers mtime --workflow-profile profile --workflow-profile profile
