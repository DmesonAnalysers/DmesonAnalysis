# The workflow of the systematics

- prepare a config_flow.yml
    bdt cut: correlated and uncorrelated
    flow files: make sure the thnsparse has `mass` `pt` `sp` `cent` `bkg` `sig` 
    max_worker will be automatically set; corr:12 uncorr:25
    skip_cut: make sure it is correct
- fill the test.sh
- bash test.sh

