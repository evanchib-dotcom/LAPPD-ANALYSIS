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


HighFive::File hdf5(filePath, HighFive::File::ReadOnly);
 
const char* dig = "digitizer0"
char digtype[64] = ""
 sprintf(digtype, "%s/", dig);

HighFive::Group digGroup = hdf5.getGroup(digtype);
int bits = digGroup.getAttribute("bits").read<int_32>.();


int ngroups = 
int nchannels = 


//loop over groups and channels for digitizer
char channelPath[64] = ""
char groupPath[8] ="" //we will need this for triggerwaveforms
for(size_t igr = 0; igr < ngroups; igr ++) {
    for(size_t ich = 0; ich < nchannels ; ich++ ) {
        // create filepath for channel 
        sprintf(channelPath, "%s/gr%lu/ch%lu/", dig, igr, ich);
        HighFive::Group chanData = hdf5.getGroup(channelPath);
        HighFive:: DataSet wfWrapper = chanData.getDataSet("samples");
        std::vector<size_t> dim = wfWrapper.getSpace().getDimension();
        size_t nEvents = dim[0];
        size_t nSamples = dim[1];
        //get dimension and createArray in memory 
        std::vector<std::vector<float>> waveform = wfWrapper.read(waveform);
        getTime(waveform);
        sprintf(groupPath, "%s/gr%lu/", dig, igr);
        HighFive::Group gr = hdf5.getGroup(groupPath);
        HighFive::DataSet trigWrap = gr.getDataSet("tr") //btw assuming that all these key words unique to format used same ones as in tanners
        std::vector<std::vector<float>> trWaveform = trigWrap.read(trWaveform);






    }
}


 std::vector<std::vector<float>> calcTime(std::vector<std::vector<float>> waveform) {



    
 }








