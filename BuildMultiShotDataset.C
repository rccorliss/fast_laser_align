void LearnParams(TTree* inTree, float* fpos, float* start, float* stepsize, float* tolerance, int* nsteps, int* major, int* minor);

void BuildMultiShotDataset(const char* originalRootFile,  const int numShots, const int iterations=1){
    //const int numShots = inShots;
       // Open the original ROOT file
    TFile* originalFile = TFile::Open(originalRootFile);
    if (!originalFile || originalFile->IsZombie()) {
        std::cerr << "Error opening original ROOT file: " << originalRootFile << std::endl;
        return;
    }

    //assume the original file ends in ".root", and we want to create a new file with the same name but with "_multiShot" appended before the ".root"

    //get the original file name without the .root at the end:
    std::string inputFileBase(originalRootFile);
    inputFileBase = inputFileBase.substr(0, inputFileBase.find_last_of("."));

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
    const int numFixedEntries = 5; // Number of fixed entries for each shot (u,d,l,r))
    float fposNew[4][numShots+numFixedEntries], pd1New[4][numShots+numFixedEntries], pd2New[4][numShots+numFixedEntries];
    float dx1df[4], dx2df[4], dy1df[4], dy2df[4]; // Derivatives in terms of motor motion

    // Create branches in the new tree
    int n=numShots;
    int nf= numFixedEntries; // Number of fixed entries
    tree->Branch("n", &n, "n/I");
    tree->Branch("nf", &nf, "nf/I");
    for (int i = 0; i < 4; ++i) {
        tree->Branch(Form("fpos%d", i), &(fposNew[i][0]), Form("fpos%d[nf]/F", i));
        tree->Branch(Form("pd1_%d", i), &(pd1New[i][0]), Form("pd1_%d[nf]/F", i));
        tree->Branch(Form("pd2_%d", i), &(pd2New[i][0]), Form("pd2_%d[nf]/F", i));
        tree->Branch(Form("fpos_r%d", i), &(fposNew[i][numFixedEntries]), Form("fpos_r%d[n]/F", i));
        tree->Branch(Form("pd1_r%d", i), &(pd1New[i][numFixedEntries]), Form("pd1_r%d[n]/F", i));
        tree->Branch(Form("pd2_r%d", i), &(pd2New[i][numFixedEntries]), Form("pd2_r%d[n]/F", i));
        tree->Branch(Form("dx1df%d", i), &dx1df[i], Form("dx1df%d/F", i));
        tree->Branch(Form("dx2df%d", i), &dx2df[i], Form("dx2df%d/F", i));
        tree->Branch(Form("dy1df%d", i), &dy1df[i], Form("dy1df%d/F", i));
        tree->Branch(Form("dy2df%d", i), &dy2df[i], Form("dy2df%d/F", i));
    }

    //also add a distance-from-center value for network training, and other new derived variables.
    float d2[numShots+numFixedEntries], d3[numShots+numFixedEntries], dist[numShots+numFixedEntries];
    float x1[numShots+numFixedEntries], y1[numShots+numFixedEntries];
    float x2[numShots+numFixedEntries], y2[numShots+numFixedEntries];
    float s1[numShots+numFixedEntries], s2[numShots+numFixedEntries];
    float Df2, Df3; // Best guess delta in fpos to get to zero,zero in diode space.
    int neighbors; //boolean packing of whether we have neighbors in each direction (up, down, left, right)
    tree->Branch("neighbors", &neighbors, "neighbors/I");
    tree->Branch("dist", &dist, "dist[nf]/F");
    tree->Branch("d2", &d2, "d2[nf]/F");
    tree->Branch("d3", &d3, "d3[nf]/F");
    tree->Branch("x1", &x1, "x1[nf]/F");
    tree->Branch("y1", &y1, "y1[nf]/F");
    tree->Branch("x2", &x2, "x2[nf]/F");
    tree->Branch("y2", &y2, "y2[nf]/F");
    tree->Branch("s1", &s1, "s1[nf]/F");
    tree->Branch("s2", &s2, "s2[nf]/F");
    tree->Branch("Df2",&Df2, "Df2/F");
    tree->Branch("Df3",&Df3, "Df3/F");
    // Add branches for the step entries
    tree->Branch("dists", &(dist[numFixedEntries]), "dists[n]/F");
    tree->Branch("ds2", &(d2[numFixedEntries]), "ds2[n]/F");
    tree->Branch("ds3", &(d3[numFixedEntries]), "ds3[n]/F");
    tree->Branch("xs1", &(x1[numFixedEntries]), "xs1[n]/F");
    tree->Branch("ys1", &(y1[numFixedEntries]), "ys1[n]/F");
    tree->Branch("xs2", &(x2[numFixedEntries]), "xs2[n]/F");
    tree->Branch("ys2", &(y2[numFixedEntries]), "ys2[n]/F");
    tree->Branch("ss1", &(s1[numFixedEntries]), "ss1[n]/F");
    tree->Branch("ss2", &(s2[numFixedEntries]), "ss2[n]/F");

    //prep the 2D grid so we can read about gradients etc.
    float start[4], stepsize[4], tolerance[4];
    int nsteps[4];
    int major, minor; //major and minor axis indices
    LearnParams(inTree,fpos, start,stepsize,tolerance,nsteps,&major,&minor);
//we will look for how many multiples of stepsize[i] we are from start[i], and if the residual is less than tolerance[i], we can confidently put the data in that bin.
    // Print the learned parameters using printf
    printf("Learned parameters:\n");
    for (int i = 0; i < 4; ++i) {
        printf("fpos[%d]: start = %.2f, step size = %.2f, tolerance = %.2f, nsteps = %d\n", i, start[i], stepsize[i], tolerance[i], nsteps[i]);
    }
    float fpos_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    float pd1_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    float pd2_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    float x1_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float x2_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float y1_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float y2_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float d2_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float d3_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float dist_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float s1_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];
    float s2_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1];


    Long64_t nEntries = inTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        inTree->GetEntry(entry);
        // Fill the new tree with the data from the original tree
        //find where we are in the grid:
        int bin[4];
        for (int i = 0; i < 4; ++i) {
            if (nsteps[i] == 0) {
                bin[i] = 0; //if there are no steps, we are at the start.
            }else{
                bin[i] = (fpos[i] - start[i]+tolerance[i]) / stepsize[i];    
                if (bin[i] < 0) {
                    bin[i] = 0;
                    printf("Warning: fpos[%d] = %.2f +tolerance is below start[%d] = %.2f, setting bin[%d] to 0.\n", i, fpos[i], i, start[i], i);
                }
                if (bin[i] > nsteps[i]) {
                    bin[i] = nsteps[i];
                    printf("Warning: fpos[%d] = %.2f +tolerance is above start[%d] = %.2f + %.2f*nsteps[%d] = %d, setting bin[%d] to %d.\n", i, fpos[i], i, start[i], stepsize[i],i, nsteps[i], i, nsteps[i]);
                }
            }
        }        
        

        // Fill this grid point with the data from the original tree:
        int shot=0; //zero is 'where we are'.
        s1[shot] = 0;
        s2[shot] = 0;
        for (int i = 0; i < 4; ++i) {
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


        for (int i = 0; i < 4; ++i) {
            fpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][i] = fpos[i];
            pd1_grid[bin[0]][bin[1]][bin[2]][bin[3]][i] = pd1[i];
            pd2_grid[bin[0]][bin[1]][bin[2]][bin[3]][i] = pd2[i];
        }
        x1_grid[bin[0]][bin[1]][bin[2]][bin[3]] = x1[shot];
        x2_grid[bin[0]][bin[1]][bin[2]][bin[3]] = x2[shot];
        y1_grid[bin[0]][bin[1]][bin[2]][bin[3]] = y1[shot];
        y2_grid[bin[0]][bin[1]][bin[2]][bin[3]] = y2[shot];
        d2_grid[bin[0]][bin[1]][bin[2]][bin[3]] = d2[shot];
        d3_grid[bin[0]][bin[1]][bin[2]][bin[3]] = d3[shot];
        dist_grid[bin[0]][bin[1]][bin[2]][bin[3]] = dist[shot];
        s1_grid[bin[0]][bin[1]][bin[2]][bin[3]] = s1[shot];
        s2_grid[bin[0]][bin[1]][bin[2]][bin[3]] = s2[shot];
    }
    //we have now filled the grid with the data from the original tree.

    //create an additional grid for the derivatives, which we will use to fill the new tree:
    float dx1_dfpos_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    float dx2_dfpos_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    float dy1_dfpos_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    float dy2_dfpos_grid[nsteps[0]+1][nsteps[1]+1][nsteps[2]+1][nsteps[3]+1][4];
    //coord, coord, coord, coord, derivatives in each fpos direction.
    //loop over the grid and calculate the derivatives in terms of motor motion.
    for (int i = 0; i < nsteps[0]+1; ++i) {
        for (int j = 0; j < nsteps[1]+1; ++j) {
            for (int k = 0; k < nsteps[2]+1; ++k) {
                for (int l = 0; l < nsteps[3]+1; ++l) {
                    //calculate the derivatives in terms of motor motion.
                    for(int ii=0;ii<4;ii++){
                        dx1_dfpos_grid[i][j][k][l][ii] = 0;
                        dx2_dfpos_grid[i][j][k][l][ii] = 0;
                        dy1_dfpos_grid[i][j][k][l][ii] = 0;
                        dy2_dfpos_grid[i][j][k][l][ii] = 0;

                        //get the adjacent bin in the ith direction, and calculate the derivatives. 
                        int adjBin[4] = {i,j,k,l}; //start with the current bin
                        adjBin[ii] += 1; //move one step in the i-th direction
                        if (adjBin[ii] <= nsteps[ii]) { //check if the adjacent bin is within bounds
                            dx1_dfpos_grid[i][j][k][l][ii] = (x1_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]] - x1_grid[i][j][k][l]) / (fpos_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][ii] - fpos_grid[i][j][k][l][ii]);
                            dx2_dfpos_grid[i][j][k][l][ii] = (x2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]] - x2_grid[i][j][k][l]) / (fpos_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][ii] - fpos_grid[i][j][k][l][ii]);
                            dy1_dfpos_grid[i][j][k][l][ii] = (y1_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]] - y1_grid[i][j][k][l]) / (fpos_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][ii] - fpos_grid[i][j][k][l][ii]);
                            dy2_dfpos_grid[i][j][k][l][ii] = (y2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]] - y2_grid[i][j][k][l]) / (fpos_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][ii] - fpos_grid[i][j][k][l][ii]);
                        } 
                    }
                }
            }
        }
    }
  


    // Now fill the output tree with the spatially correlated data from the grids:



    for (int iteration = 0; iteration < iterations; ++iteration) {
       //loop over every coordinate in the grid, but linearize to avoid nested loops.
        printf("Iteration %d/%d\n", iteration+1, iterations);
       for (int i = 0; i < (nsteps[0]+1)*(nsteps[1]+1)*(nsteps[2]+1)*(nsteps[3]+1); ++i) {
            int bin[4];
            //linearize the index to the 4D grid:
            bin[0] = i / ((nsteps[1]+1) * (nsteps[2]+1) * (nsteps[3]+1));
            bin[1] = (i / ((nsteps[2]+1) * (nsteps[3]+1))) % (nsteps[1]+1);
            bin[2] = (i / (nsteps[3]+1)) % (nsteps[2]+1);
            bin[3] = i % (nsteps[3]+1);

            // Fill the tree with the data from this grid point
            for (int i=0; i < 4; ++i) {//the four fpos, pd1, pd2 values
                fposNew[i][0] = fpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][i];
                pd1New[i][0] = pd1_grid[bin[0]][bin[1]][bin[2]][bin[3]][i];
                pd2New[i][0] = pd2_grid[bin[0]][bin[1]][bin[2]][bin[3]][i];
            }
            x1[0] = x1_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            x2[0] = x2_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            y1[0] = y1_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            y2[0] = y2_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            d2[0] = d2_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            d3[0] = d3_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            dist[0] = dist_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            s1[0] = s1_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            s2[0] = s2_grid[bin[0]][bin[1]][bin[2]][bin[3]];
            for (int j = 0; j < 4; j++) {
                dx1df[j] = dx1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][j];
                dx2df[j] = dx2_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][j];
                dy1df[j] = dy1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][j];
                dy2df[j] = dy2_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][j];
            }

            //and assume we know the derivatives, so we can calculate the size of the step we need to take to get to zero, zero:
/*
    (df2/dx*Dx-df2/dy*Dy)
-----------------------------  =Df3
(dx/df3*df2/dx-dy/df3*df2/dy)

( Dx-dx/df3*Df3 ) / ( dx/df2 ) =Df2
*/            
            Df3=0;
            Df2=0;
            Df3= x1[0]/(dx1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][2])
                    - y1[0]/(dy1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][2]);
            Df3/= dx1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][3]/
                dx1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][2]
                -dy1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][3]/
                dy1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][2];
            Df2=x1[0]-dx1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][3]*Df3/
                dx1_dfpos_grid[bin[0]][bin[1]][bin[2]][bin[3]][2];

            //and if we make those changes, we can look at where we actually are in the grid:
            int binNew[4];
            binNew[0]=bin[0];
            binNew[1]=bin[1];
            binNew[2]=((fposNew[2][0]-Df2 - start[2] + tolerance[2]) / stepsize[2]);
            binNew[3]=((fposNew[3][0]-Df3 - start[3] + tolerance[3]) / stepsize[3]);

            //fill the entry in the new tree with the data from the grid at this point:
            //check if the new bin is within bounds
            for (int j = 0; j < 4; ++j) {
                if (binNew[j] < 0 || binNew[j] > nsteps[j]) {
              printf("Warning: binNew[%d] = %d is out of bounds, setting it to 0.\n", j, binNew[j]);
                    binNew[j] = 0; //if out of bounds, set to 0
                }
            }
            // Fill the tree with the data from this grid point
            for (int i=0; i < 4; i++) {//the four fpos, pd1, pd2 values
                    fposNew[i][numFixedEntries] = fpos_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]][i];
            pd1New[i][numFixedEntries] = pd1_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]][i];
            pd2New[i][numFixedEntries] = pd2_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]][i];
                }
                x1[numFixedEntries] = x1_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                x2[numFixedEntries] = x2_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                y1[numFixedEntries] = y1_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                y2[numFixedEntries] = y2_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                d2[numFixedEntries] = d2_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                d3[numFixedEntries] = d3_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                dist[numFixedEntries] = dist_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                s1[numFixedEntries] = s1_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];
                s2[numFixedEntries] = s2_grid[binNew[0]][binNew[1]][binNew[2]][binNew[3]];




            //now fill the next four entries with the ones in the adjacent bins, if they exist.
            neighbors = 0; //reset the neighbors counter
            for (int j = 0; j < 4; ++j) {//calculate and fill adjacent bins
                //calculate the adjacent coordinate, by moving one step in each of the two directions that have a step size.
                int adjBin[4] = {bin[0], bin[1], bin[2], bin[3]};
                if (j < 2) {
                    adjBin[major] += j%2?-1:1;
                } else {
                    adjBin[minor] += j%2?-1:1;
                }
                //check if the adjacent bin is within bounds
                bool inBounds = true;
                for (int k = 0; k < 4; ++k) {
                    if (adjBin[k] < 0 || adjBin[k] > nsteps[k]) {
                        inBounds = false;
                        break;
                    }
                }
                if (inBounds) {
                    neighbors |= (1 << j); //set the j-th bit in the neighbors variable
                    // Fill the tree with the data from this adjacent grid point
                    for (int i=0; i < 4; ++i) {
                        fposNew[i][j+1] = fpos_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][i];
                        pd1New[i][j+1] = pd1_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][i];
                        pd2New[i][j+1] = pd2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]][i];
                    }
                    x1[j+1] = x1_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    x2[j+1] = x2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    y1[j+1] = y1_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    y2[j+1] = y2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    d2[j+1] = d2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    d3[j+1] = d3_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    dist[j+1] = dist_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    s1[j+1] = s1_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                    s2[j+1] = s2_grid[adjBin[0]][adjBin[1]][adjBin[2]][adjBin[3]];
                } else {
                    //if the adjacent bin is out of bounds, we just fill it with zeros.
                    for (int i=0; i < 4; ++i) {
                        fposNew[i][j+1] = -99;
                        pd1New[i][j+1] = -99;
                        pd2New[i][j+1] = -99;
                    }   
                    x1[j+1] = -99;
                    x2[j+1] = -99;
                    y1[j+1] = -99;
                    y2[j+1] = -99;
                    d2[j+1] = -99;
                    d3[j+1] = -99;
                    dist[j+1] = -99;
                    s1[j+1] = -99;
                    s2[j+1] = -99;
                }
            }
            //now fill with the next step in the pattern, then random entries, which are just random entries from the original tree, or, better, random entries from the grid, since we already calculated some values
            for (int shot = 1; shot < numShots; ++shot) {
                int randBin[4];
                //generate random indices for the grid, ensuring they are within bounds
                for (int i = 0; i < 4; ++i) {
                    randBin[i] = gRandom->Uniform(0, nsteps[i]+1);
                }
                // Fill the tree with the data from this random grid point
                for (int i=0; i < 4; ++i) {
                    fposNew[i][numFixedEntries+shot] = fpos_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]][i];
                    pd1New[i][numFixedEntries+shot] = pd1_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]][i];
                    pd2New[i][numFixedEntries+shot] = pd2_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]][i];
                }
                x1[numFixedEntries+shot] = x1_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                x2[numFixedEntries+shot] = x2_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                y1[numFixedEntries+shot] = y1_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                y2[numFixedEntries+shot] = y2_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                d2[numFixedEntries+shot] = d2_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                d3[numFixedEntries+shot] = d3_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                dist[numFixedEntries+shot] = dist_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                s1[numFixedEntries+shot] = s1_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];
                s2[numFixedEntries+shot] = s2_grid[randBin[0]][randBin[1]][randBin[2]][randBin[3]];


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



void LearnParams(TTree* inTree, float* fpos, float* start, float* stepsize, float* tolerance, int* nsteps, int* major, int* minor) {
    //read through the input tree to figure out the starting coordinates, step sizes, and number of steps in each direction.


    float fposold[4];
    int nEntries = inTree->GetEntries();
    inTree->GetEntry(0);
    for (int i = 0; i < 4; ++i) {
        fposold[i] = fpos[i];
        start[i] = fpos[i];
    }
    vector<float> stepSizes[4];
    for (int i = 0; i < 4; ++i) {
        stepSizes[i].clear();
        stepsize[i] = 0;
        tolerance[i] = 0;
        nsteps[i] = 0;
    }
    int nLocalSteps[4]= {0, 0, 0, 0};

    for (int entry = 1; entry < nEntries; ++entry) {
        inTree->GetEntry(entry);
        //figure out which direction the laser moved in:
        int nMoved = 0;
        for (int i = 0; i < 4; ++i) {
            float delta=fpos[i] - fposold[i];
            if (nLocalSteps[i] == 0) {
                nMoved++;
            }
            if (fabs(delta) > 0.0005) { 
                stepSizes[i].push_back(delta);
                nLocalSteps[i]++;
            }
        }
        for (int i = 0; i < 4; ++i) {
            fposold[i] = fpos[i];
        }
        assert(nMoved > 0 && nMoved<3);
    }
  
    //okay, now we have the step size patterns, and we can look at what sort of motion we have:
    //1)First, figure out the hierarchy of the step sizes.
    vector<int> most(2, -1);
    vector<int> nmost(2, 0); //the number of steps in the two largest step sizes
    for( int i=0; i<4; ++i) {
        if(nLocalSteps[i] > nmost[0]) {
            nmost[1] = nmost[0];
            most[1] = most[0];
            nmost[0] = nLocalSteps[i];
            most[0] = i;
        } else if(nLocalSteps[i] > nmost[1]) {
            nmost[1] = nLocalSteps[i];
            most[1] = i;
        }
    } //most[0] is the largest, most[1] is the second largest.
    //2)the largest number of steps is the main direction, and the second largest is the switchback direction.

    //we can look at the pattern in the main direction.  If it has a consistent step size, differing by a sign, then we are snaking.  If it has two different step sizes, then we are zigzagging.
    //fit the step size data to a pair of gaussians, one constrained to be positive, the other negative.  If the fit is good, then we are zigzagging, otherwise we are snaking.
    float mean[3] = {0, 0, 0};//neg and pos mean position, then the 2nd most's mean etc.
    float numdir[3] = {0, 0, 0};//neg and pos number of entries
    float sigma[3] = {0, 0, 0};//neg and pos sigma
    for (int i=0;i<nmost[0];i++){
        if (stepSizes[most[0]][i] < 0) {
            mean[0] += stepSizes[most[0]][i];
            numdir[0]++;
        } else {
            mean[1] += stepSizes[most[0]][i];
            numdir[1]++;
        }
    }
    mean[0] /= numdir[0];
    mean[1] /= numdir[1];
    for (int i=0;i<nmost[0];i++){
        if (stepSizes[most[0]][i] < 0) {
            sigma[0] += (stepSizes[most[0]][i] - mean[0]) * (stepSizes[most[0]][i] - mean[0]);
        } else {
            sigma[1] += (stepSizes[most[0]][i] - mean[1]) * (stepSizes[most[0]][i] - mean[1]);
        }
    }
    sigma[0] = std::sqrt(sigma[0] / numdir[0])+0.01;
    sigma[1] = std::sqrt(sigma[1] / numdir[1])+0.01;

    //check if the two means are consistent with each other, using proper gaussian error propagation.  If so, we are switching back (=snaking):
    float meandiff=mean[0]+mean[1];
    float sigmadiff = std::sqrt(sigma[0]*sigma[0]+sigma[1]*sigma[1]);
    bool zigzag = false; //default is snaking
    if (fabs(meandiff) > 3*sigmadiff) {
        //if we're not consistent with zero, then we are zigzagging.
        zigzag = 1; //we are zigzagging.
        printf("Zigzag pattern detected: mean[0]=%.2f, mean[1]=%.2f, sigmadiff=%.2f\n", mean[0], mean[1], sigmadiff);
        printf("major axis: %.2f, minor axis: %.2f\n", most[0], most[1]);
        printf("nmost[0]=%d, nmost[1]=%d\n", nmost[0], nmost[1]);
        printf("mean[0]=%.2f, mean[1]=%.2f, sigma[0]=%.2f, sigma[1]=%.2f\n", mean[0], mean[1], sigma[0], sigma[1]);

    } else {
        printf("Snaking pattern detected.\n");
        printf("major axis: %.2f, minor axis: %.2f\n", most[0], most[1]);
        printf("nmost[0]=%d, nmost[1]=%d\n", nmost[0], nmost[1]);
        printf("mean[0]=%.2f, mean[1]=%.2f, sigma[0]=%.2f, sigma[1]=%.2f\n", mean[0], mean[1], sigma[0], sigma[1]); 
    }
            for (int i=0;i<4;i++){
            printf("Step sizes in direction %d:\n", i);
            for (int j=0;j<stepSizes[i].size();j++){
                printf("%.2f ", stepSizes[i][j]);
            }
            if (i== most[0]) {
                printf(" <-- main direction");
            } else if (i == most[1]) {
                printf(" <-- secondary direction");
            }
            printf("\n");
        }

    //the second direction is easier, since we only have one step size.  We can just take the mean and sigma of the step sizes in that direction.
    float mean2 = 0;
    float sigma2 = 0;
    for (int i=0;i<nmost[1];i++){
        mean2 += stepSizes[most[1]][i];
    }
    mean2 /= nmost[1];
    for (int i=0;i<nmost[1];i++){
        sigma2 += (stepSizes[most[1]][i] - mean2) * (stepSizes[most[1]][i] - mean2);
    }
    sigma2 = std::sqrt(sigma2 / nmost[1])+0.01;

    if(zigzag){
        //Extract the stepsize, number of steps, and how tightly we can cut (toleraance) for the main direction.  In zigzag, the stepsize is the one that's more common:
        if (numdir[0] > numdir[1]) {
            stepsize[most[0]] = mean[0];
            tolerance[most[0]] = 3*sigma[0];
         } else {
            stepsize[most[0]] = mean[1];
            tolerance[most[0]] = 3*sigma[1];
         }
    
        //the number of steps we counted will be nsteps*nrows+ncarriagereturns, since we do not have the final 'carriage return'.  hence the number we should report back is:
        nsteps[most[0]] = (nmost[0]-nmost[1])/(nmost[1]+1); 

        //and the second direction is just a simple count:
        nsteps[most[1]] = nmost[1];
        stepsize[most[1]] = mean2;
        tolerance[most[1]] = 3*sigma2;

        //and the remainder are zero already.

    }else{
        //we're snaking, so we have a single abs(step size) in the main direction, and a single step size in the second direction.
        stepsize[most[0]] = (abs(mean[1])+abs(mean[0]))/2;//average the two.
        tolerance[most[0]] =1.5*sigmadiff;//technically correct for a 3 sigma range, but we could be more permissive.

        //the number of steps we counted will be nsteps*nrows, since the carriage returns do not register on the main direction list. hence the number we should report back is:
        nsteps[most[0]] = nmost[0]/(nmost[1]+1); 
        //and the second direction is still just a simple count:
        nsteps[most[1]] = nmost[1];
        stepsize[most[1]] = mean2;
        tolerance[most[1]] = 3*sigma2;
    }
*major= most[0];
*minor= most[1];
    return;
}
