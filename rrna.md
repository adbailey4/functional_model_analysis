## Call 2mG and 6mmA on 16S rRNA

1. Convert predicted models to SA models
   * `/Users/andrewbailey/CLionProjects/functional_model_analysis/src/convert_csv_to_signalalign.py --dir /Users/andrewbailey/Desktop/rna_GSE124309_r9.4_mods_model --output_dir /Users/andrewbailey/Desktop/sa_rna_GSE124309_r9.4_mods_model --base_model /Users/andrewbailey/CLionProjects/functional_model_analysis/basesline_model/baseline_rna_ACGTXYZ_RNA.model --num_threads 4 --rna`
   
2. Generate Tombo Plots for 16S rRNA
   * `tombo resquiggle /home/ubuntu/ecoli_16s/tombo_analysis/temp_fast5 /home/ubuntu/ecoli_16s/reference/J01859.1.fa --processes 8 --num-most-common-errors 5 --rna`
   * `tombo detect_modifications de_novo --fast5-basedirs /home/ubuntu/ecoli_16s/tombo_analysis/temp_fast5 --statistics-file-basename 16S_de_novo --rna --processes 8`


3. Generate Canonical Positions
   * `run_16s_rrna/ecoli_rRNA_site_selection.ipynb`
    
   ```
   from py3helpers.seq_tools import ReferenceHandler
   from py3helpers.utils import merge_lists
   from tombo import tombo_helper, tombo_stats, resquiggle
   rrna_16S_stats = "/Users/andrewbailey/CLionProjects/functional_model_analysis/run_16s_rrna/16S_de_novo.tombo.stats"
   save_fig_path = None
   assert os.path.exists(rrna_16S_stats)
   ts = tombo_stats.TomboStats(rrna_16S_stats)
   for contig in ts:
      all_data = contig[4]
      position_fraction_modified = {x[2]: x[0] for x in all_data}
      break

   ecoli = ReferenceHandler("/Users/andrewbailey/CLionProjects/functional_model_analysis/run_16s_rrna/baseline_model/temp/J01859.1.fa")
   # https://www.ncbi.nlm.nih.gov/nuccore/J01859
   variations = [76, 80, 88, 89, 90, 92, 179, 182, 193, 194, 267, 272, 283, 285, 348, 630, 632, 640, 853, 915, 965, 1027, 1071, 1074, 1099, 1206, 1280, 1321, 1401, 1402, 1405, 1406, 1490, 1493, 1497, 1517, 1518]                
   mods = [515, 526, 965, 966, 1206, 1401, 1406, 1497, 1515, 1517, 1518]
   miss_pos = variations + mods
   must_miss = merge_lists([list(range(x-20, x+20)) for x in miss_pos])
   seq = ecoli.get_sequence("J01859.1", 0, 1518)
   pos = []
   min_gap = 18
   wait = False
   curr_gap = 0
   for i, x in enumerate(seq):
      if wait:
          curr_gap += 1
          if curr_gap == min_gap:
              curr_gap = 0
              wait = False
      elif i > 200 and i not in must_miss and i < 1480 and position_fraction_modified[i] < .05:
          if x == "G":
              pos.append((i, x))
              wait = True
   
   print(len(pos))
   my_map = {"G": "GZ",
           "A": "AXY"}
   with open("/Users/andrewbailey/CLionProjects/functional_model_analysis/run_16s_rrna/baseline_model/temp/16S_final.positions", "a") as fh:
      for x in pos:
          print("\t".join(["J01859.1", str(x[0]), "+", x[1], my_map[x[1]]]), file=fh)

   ```
   
4. re-run signalAlign
   * sudo docker run -v /home/ubuntu/ecoli_16s:/data --entrypoint /bin/bash -it adbailey4/signalalign@sha256:a350ce89a00e23b96f2224a0ca8fc84e53ba7d44fde2c75331218c73b4833b1a
   * bash run_16s_rrna.sh bailey-ding-model-prediction/rRNA_mod_analysis/models/sa_rna_GSE124309_r9.4_mods_model/ 96 bailey-ding-model-prediction/rRNA_mod_analysis/run_3/output/
