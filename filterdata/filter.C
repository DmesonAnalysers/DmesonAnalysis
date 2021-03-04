#include "FilterTree4ML.cc"
#include "TSystem.h"
// #include "/home/alidock/DmesonAnalysis/printTree.C"

void filter13TeV();

void filter() {
 gSystem->CompileMacro("FilterTree4ML.cc");

   //  FilterTree4ML("config_LctopK0s_Data.yml");
   //  FilterTree4ML("config_LctopK0s_MC_train.yml");
    FilterTree4ML("config_LctopK0s_MC_eff.yml");



 // filter data
//  FilterTree4ML("config_LctopK0s_pt_2_4_Data_13TeV.yml");
//  FilterTree4ML("config_LctopK0s_pt_4_6_Data_13TeV.yml");
//  FilterTree4ML("config_LctopK0s_pt_6_50_Data_13TeV.yml");

//  // filter MC_Lcdedicated (training)
//  FilterTree4ML("config_LctopK0s_pt_2_4_MC_13TeV.yml");
//  FilterTree4ML("config_LctopK0s_pt_4_6_MC_13TeV.yml");
//  FilterTree4ML("config_LctopK0s_pt_6_50_MC_13TeV.yml");

//  // filter MC Efficiency (apply)
//  FilterTree4ML("config_LctopK0s_pt_2_4_MC_13TeV_eff.yml");
//  FilterTree4ML("config_LctopK0s_pt_4_6_MC_13TeV_eff.yml");
//  FilterTree4ML("config_LctopK0s_pt_6_50_MC_13TeV_eff.yml");

 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/Data_filter_pT_2_4.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/Data_filter_pT_4_6.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/Data_filter_pT_6_50.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/Prompt_filter_pT_2_4.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/Prompt_filter_pT_4_6.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/Prompt_filter_pT_6_50.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/FD_filter_pT_2_4.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/FD_filter_pT_4_6.root");
 // tree2pdf("/home/alidock/DmesonAnalysis/filterdata/filtered/13TeV/FD_filter_pT_6_50.root");
}

