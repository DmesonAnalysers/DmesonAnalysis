# The workflow of the systematics

- fill the flow configuration like flow/config/config_flow_D0_3050_mother.yml
    - - include all pt bin and all bdt cut

- fill the flow/systematics/modifications_config_fit.yml
    - - the flow_file is a list, that contains all AnRes.root splite by pT. if you don't have, please run it again, refer to this configuraion file tool/config_pre_sys.yml
    - - Rebin: only 2 values will be used
    - - inv_mass_bins: steps * 2 MassMin * 2 MassMax
    - - use_inv_mass_bins: not enable
    - - sigma: only considered edge, to be add central value
    - - ignore the remains

- bash test.sh
  - - fill the varibles
  - - uncomment the part based on your needs, like run the central value fisrtly to obtain the sigma
  - - since running pt bin by pt bin, if the printed info is the next pt bin, you can just go to the next step to have a look at the syst_multitrial reuslts

- collect the results by hand 
    - - `python3 compute_syst_multitrial.py $output_dir/trails/trails_pt_30_35 $output_dir/pre -o $output_dir/sys -p`

    - - modify the `trails_pt_30_35`
