## Call 2mG and 6mmA on 16S rRNA

1. Convert predicted models to SA models
    * `/Users/andrewbailey/CLionProjects/functional_model_analysis/src/convert_csv_to_signalalign.py --dir /Users/andrewbailey/Desktop/rna_GSE124309_r9.4_mods_model --output_dir /Users/andrewbailey/Desktop/sa_rna_GSE124309_r9.4_mods_model --base_model /Users/andrewbailey/CLionProjects/functional_model_analysis/basesline_model/baseline_rna_ACGTXYZ_RNA.model --num_threads 4 --rna`
   
2. Generate Canonical Positions
   ```
   from py3helpers.seq_tools import ReferenceHandler
   from py3helpers.utils import merge_lists
   ecoli = ReferenceHandler("/Users/andrewbailey/CLionProjects/functional_model_analysis/basesline_model/temp/J01859.1.fa")
   miss_me_with_that = [515, 526, 965, 966, 1206, 1401, 1406, 1497, 1515, 1517, 1518]
   must_miss = merge_lists([list(range(x-20, x+20)) for x in miss_me_with_that])
   seq = ecoli.get_sequence("J01859.1", 0, 1518)
   pos = []
   min_gap = 41
   wait = False
   curr_gap = 0
   for i, x in enumerate(seq):
       if wait:
           curr_gap += 1
           if curr_gap == min_gap:
               curr_gap = 0
               wait = False
       elif i > 200 and i not in must_miss and i < 1480:
           if x == "A":
               pos.append((i, x))
               wait = True
           if x == "G":
               pos.append((i, x))
               wait = True
   
   print(len(pos))
   my_map = {"G": "GZ",
            "A": "AXY"}
   with open("/Users/andrewbailey/CLionProjects/functional_model_analysis/basesline_model/temp/16S_final.positions", "a") as fh:
       for x in pos:
           print("\t".join(["J01859.1", str(x[0]), "+", x[1], my_map[x[1]]]), file=fh)
   
   ```
   
3. re-run signalAlign
   * sudo docker run -v /home/ubuntu/ecoli_16s:/data --entrypoint /bin/bash -it adbailey4/signalalign@sha256:a350ce89a00e23b96f2224a0ca8fc84e53ba7d44fde2c75331218c73b4833b1a
   * bash run_16s_rrna.sh bailey-ding-model-prediction/rRNA_mod_analysis/run_1/sa_rna_GSE124309_r9.4_mods_model/ 8 bailey-ding-model-prediction/rRNA_mod_analysis/run_1/sa_rna_GSE124309_r9.4_mods_output/
