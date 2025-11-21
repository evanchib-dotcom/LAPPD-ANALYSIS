#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric> // for std::accumulate
#include <cmath>

std::vector<std::vector<int>> n_samples_above_thresh(const std::string& filepath){
    TFile *file = TFile::Open(filepath.c_str(), "READ");
        if(!file || file->isZombie()) {
            std::cerr << "Error opening file: " << file_path << std::endl;
            return {};
    }

    TTree *outputTree = (*TTree)file->Get("output");

    std::vector<std:vector<float>> *waveforms = nullptr;
    ouputTree->SetBranchAddress("waveforms", &waveforms);

     Long64_t nEntries = outputTree->GetEntries(); 
     std::vector<std::vector<int>> n_above_thresh ; 
        for(i = 0; i < nEntries; i++) {
            outputTree->GetEntry(i);
            std::vector<int> n_above_event;
            for(j = 0; j < 30; j++) {
                std::vector<float> &currWfFull = waveforms->at(j) ;
                std::vector<float> wf(currWfFull.begin() + 10, currWfFull.begin() + 1014);
                int pedStart = 500;
                int pedEnd = 900;
                int pedLen = pedEnd - pedStart; 
                
            }



     }



}                float ped = std::accumulate(wf.begin() + pedStart, wf.begin() + pedEnd, 0.0)/pedLen; 
