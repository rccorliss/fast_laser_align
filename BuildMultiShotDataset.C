void BuildMultiShotDataset(const char* originalRootFile,  const int numShots, const int iterations=1){
    //const int numShots = inShots;
       // Open the original ROOT file
    TFile* originalFile = TFile::Open(originalRootFile);
    if (!originalFile || originalFile->IsZombie()) {
        std::cerr << "Error opening original ROOT file: " << originalRootFile << std::endl;
        return;
    }

    //assume the original file ends in ".root", and we want to create a new file with the same name but with "_multiShot" appended before the ".root"
    std::string outputRootFileStr(originalRootFile);
    size_t pos = outputRootFileStr.find_last_of(".");
    if (pos == std::string::npos) {
        std::cerr << "Error: Original ROOT file name does not contain an extension." << std::endl;
        originalFile->Close();
        return;
    }
    outputRootFileStr.insert(pos, "_multiShot" + std::to_string(numShots));
    //outputRootFileStr += ".root";
    const char* outputRootFile = outputRootFileStr.c_str();

    // Create a new ROOT file for the output
    TFile* outputFile = new TFile(outputRootFile, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output ROOT file: " << outputRootFile << std::endl;
        return;
    }

    // Create a new tree in the output file
    TTree* tree = new TTree("Scan4D", "Multi-shot dataset");

    //load the same tree from the original file
    TTree* inTree = (TTree*)originalFile->Get("Scan4D");
    if (!inTree) {
        std::cerr << "Error: Tree 'Scan4D' not found in original ROOT file." << std::endl;
        outputFile->Close();
        return;
    }

    // load the branches from the original tree
    float  fpos[4],pd1[4],pd2[4];
    for (int i = 0; i < 4; ++i) {
        inTree->SetBranchAddress(Form("fpos%d", i), &fpos[i]);
        inTree->SetBranchAddress(Form("pd1_%d", i), &pd1[i]);
        inTree->SetBranchAddress(Form("pd2_%d", i), &pd2[i]);
    }

    //create the variables for the new tree.  Here we will combine multiple entries from the original tree, so we need an additional index in the variables
    float fposNew[4][numShots*2], pd1New[4][numShots*2], pd2New[4][numShots*2];

    // Create branches in the new tree
    int n=numShots;
    tree->Branch("n", &n, "n/I");
    for (int i = 0; i < 4; ++i) {
        tree->Branch(Form("fpos%d", i), &(fposNew[i][0]), Form("fpos%d[n]/F", i));
        tree->Branch(Form("pd1_%d", i), &(pd1New[i][0]), Form("pd1_%d[n]/F", i));
        tree->Branch(Form("pd2_%d", i), &(pd2New[i][0]), Form("pd2_%d[n]/F", i));
        tree->Branch(Form("fpos_r%d", i), &(fposNew[i][numShots]), Form("fpos_r%d[n]/F", i));
        tree->Branch(Form("pd1_r%d", i), &(pd1New[i][numShots]), Form("pd1_r%d[n]/F", i));
        tree->Branch(Form("pd2_r%d", i), &(pd2New[i][numShots]), Form("pd2_r%d[n]/F", i));
    }

    //also add a distance-from-center value for network training, and other new derived variables.
    float d2[numShots*2], d3[numShots*2], dist[numShots*2];
    float x1[numShots*2], y1[numShots*2];
    float x2[numShots*2], y2[numShots*2];
    float s1[numShots*2], s2[numShots*2];
    tree->Branch("dist", &dist, "dist[n]/F");
    tree->Branch("d2", &d2, "d2[n]/F");
    tree->Branch("d3", &d3, "d3[n]/F");
    tree->Branch("x1", &x1, "x1[n]/F");
    tree->Branch("y1", &y1, "y1[n]/F");
    tree->Branch("x2", &x2, "x2[n]/F");
    tree->Branch("y2", &y2, "y2[n]/F");
    tree->Branch("s1", &s1, "s1[n]/F");
    tree->Branch("s2", &s2, "s2[n]/F");
    // Add branches for the random entries
    tree->Branch("distr", &(dist[numShots]), "distr[n]/F");
    tree->Branch("dr2", &(d2[numShots]), "dr2[n]/F");
    tree->Branch("dr3", &(d3[numShots]), "dr3[n]/F");
    tree->Branch("xr1", &(x1[numShots]), "xr1[n]/F");
    tree->Branch("yr1", &(y1[numShots]), "yr1[n]/F");
    tree->Branch("xr2", &(x2[numShots]), "xr2[n]/F");
    tree->Branch("yr2", &(y2[numShots]), "yr2[n]/F");
    tree->Branch("sr1", &(s1[numShots]), "sr1[n]/F");
    tree->Branch("sr2", &(s2[numShots]), "sr2[n]/F");

    Long64_t nEntries = inTree->GetEntries();
    for (int it=0;it<iterations; it++) {
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            inTree->GetEntry(entry);
            // Fill the new tree with the data from the original tree
            
            for (int shot = 0; shot < numShots*2; ++shot) {
                //for the first half we pick consecutive hits, then we draw the same number of random ones (I just didn't want to write this twice)
                //zero is the actual shot in order.
                s1[shot] = 0;
                s2[shot] = 0;
                for (int i = 0; i < 4; ++i) {
                    fposNew[i][shot] = fpos[i];
                    pd1New[i][shot] = pd1[i];
                    pd2New[i][shot] = pd2[i];
                    s1[shot] += pd1[i]; 
                    s2[shot] += pd2[i]; 
                }
                x1[shot]=(pd1[0]+pd1[3]-pd1[2]-pd1[1])/s1[shot];
                y1[shot]=(pd1[0]+pd1[2]-pd1[3]-pd1[1])/s1[shot];
                
                x2[shot]=(pd2[0]+pd2[3]-pd2[2]-pd2[1])/s1[shot];
                y2[shot]=(pd2[0]+pd2[2]-pd2[3]-pd2[1])/s2[shot];

                // Calculate the distance from the center for each shot
                d2[shot] = (fpos[2]-1);
                d3[shot] = (fpos[3]+3);
                dist[shot] = sqrt(d2[shot]*d2[shot] + d3[shot]*d3[shot]);
                //printf("center[%d]=%f\n",shot,center[shot]);
                //draw normal entries for the first numShots, then random entries for the rest
                int rndEntry = (entry+1)% nEntries; // This will cycle through the entries
                if (shot>numShots){
                    rndEntry = gRandom->Integer(nEntries);
                }
                inTree->GetEntry(rndEntry);
            }
            // Fill the tree with the new data for this entry
            tree->Fill();
        }
    }

    tree->Write();
    outputFile->Close();
    
    originalFile->Close();

    printf("Multi-shot dataset created successfully in: %s\n", outputRootFile);
    return;
}