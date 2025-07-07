#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TMVA/Factory.h>
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TRandom3.h> // For simulating data

void TMVASimpleNNRegression() {


    // --------------------------------------------------------------------------------------------------
    // 2. Prepare TMVA Factory
    // --------------------------------------------------------------------------------------------------

    // Initialize TMVA
    TMVA::Tools::Instance();

    // Create a TMVA Factory instance.
    // The first argument is the output file name for TMVA's results.
    TFile *outputTMVAFile = TFile::Open("TMVANNRG.root", "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("TMVANNRegression", outputTMVAFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D");
    // Options:
    // !V: Verbose off (default is on, can be noisy)
    // !Silent: Suppress TMVA output to stdout (can be useful for scripting)
    // Color: Use color output for TMVA log
    // DrawProgressBar: Show progress bar during training
    // Transformations: Apply transformations to input variables.
    //   I: Identity (no transformation)
    //   D: Decorrelation
    //   P: Principal Component Analysis (PCA)
    //   G: Gaussianization (transforms to Gaussian distribution)
    // We use 'I' as a baseline, but 'G' or 'D' can often improve NN performance.

    // --------------------------------------------------------------------------------------------------
    // 3. Define Variables and Target
    // --------------------------------------------------------------------------------------------------

    // Add the input variables
    factory->AddVariable("pd1_0", "pd1_0", "", 'D');
    factory->AddVariable("pd1_1", "pd1_1", "", 'D');
    factory->AddVariable("pd1_2", "pd1_2", "", 'D');
    factory->AddVariable("pd1_3", "pd1_3", "", 'D');

    // Specify the target variable for regression
    factory->AddTarget("fpos", "x_target", "", 'D'); // The 'x' is our target

    // --------------------------------------------------------------------------------------------------
    // 4. Register the data for training and testing
    // --------------------------------------------------------------------------------------------------

    TFile *inputDataFile = TFile::Open("TMVANNRG_data.root");
    TTree *dataTreeRead = (TTree*)inputDataFile->Get("DataTree");

    // TMVA handles splitting the data internally.
    // Add the entire tree as signal (no background concept in single-target regression for now).
    factory->AddRegressionTree(dataTreeRead, 1.0); // Weight 1.0 for all events

    // Set the training and testing trees (TMVA will split from the added data)
    // The "SplitMode=Random" means random splitting. "SplitMode=Block" can be useful for time series.
    // "NormMode=NumEvents" normalizes event weights by total number of events.
    // "nTrain_Regression=0" means use all events not used for test for training.
    // "nTest_Regression=0" means use all events not used for training for testing.
    // Here we will use 50% for training and 50% for testing.
    factory->PrepareTrainingAndTestTree("SplitMode=Random:V=F:nTrain_Regression=5000:nTest_Regression=5000:NormMode=NumEvents");

    // --------------------------------------------------------------------------------------------------
    // 5. Book the Multi-Layer Perceptron (Neural Network) Method
    // --------------------------------------------------------------------------------------------------

    // Method options for MLP (Neural Network)
    // "H" specifies the hidden layers architecture. Here "N+N" means two hidden layers,
    // each with 'N' neurons where N is the number of input variables.
    // "H:N+5" would mean one hidden layer with (number of inputs + 5) neurons.
    // "H:N+N+5" means two hidden layers, first with N neurons, second with N+5 neurons.
    // "NeuronType=tanh": Activation function for hidden layers (tanh, sigmoid, relu).
    // "VarTransform=N": Normalizes input variables to range [0, 1] for NN. Recommended.
    // "BP_TrainingMethod=BFGS": Backpropagation training algorithm (BFGS, GD, GDL). BFGS is often good.
    // "TrainingCycles=1000": Number of training epochs.
    // "LearningRate=1e-2": Learning rate for backpropagation.
    // "TestRate=10": Test network performance every 10 cycles.
    // "Sampling=MC": Monte Carlo sampling for training events.
    // "Architecture=CPU": Use CPU for training (GPU option exists with CUDA/OpenCL, but requires ROOT compiled with it).
    TString mlpOptions = "H:N+N:"
                         "NeuronType=tanh:"
                         "VarTransform=N:" // Normalize input variables
                         "BP_TrainingMethod=BFGS:"
                         "TrainingCycles=2000:" // Increase cycles for better learning
                         "LearningRate=1e-2:"
                         "TestRate=50:"
                         "Sampling=MCS:"; // Use Monte Carlo Sampling of events

    factory->BookMethod(TMVA::Types::kMLP, "MLP", mlpOptions);

    // --------------------------------------------------------------------------------------------------
    // 6. Train, Test, and Evaluate
    // --------------------------------------------------------------------------------------------------

    // Train all booked methods
    factory->TrainAllMethods();

    // Test all trained methods
    factory->TestAllMethods();

    // Evaluate all trained methods
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------------------------------------------
    // 7. Apply the Trained Model (Example of how to use the Reader for prediction)
    // --------------------------------------------------------------------------------------------------

    std::cout << "\n--- Using the trained MLP model for prediction ---" << std::endl;

    // Create a TMVA Reader object
    TMVA::Reader *reader = new TMVA::Reader("!V:!Silent");

    // Declare the variables for the reader. These must be float (or Double_t if compiled with that option).
    // It's crucial that these variable names and types match those used in factory->AddVariable exactly.
    Float_t r_pd1_0, r_pd1_1, r_pd1_2, r_pd1_3;
    reader->AddVariable("pd1_0", &r_pd1_0);
    reader->AddVariable("pd1_1", &r_pd1_1);
    reader->AddVariable("pd1_2", &r_pd1_2);
    reader->AddVariable("pd1_3", &r_pd1_3);

    // Book the trained MLP method
    // The second argument ("MLP") must match the name given in factory->BookMethod
    // The third argument is the TMVA weights file produced during training.
    reader->BookMVA("MLP", "TMVANNRegression/weights/TMVANNRegression_MLP.weights.xml");

    std::cout << "Testing predictions on a few example data points:" << std::endl;

    // Simulate new data points for prediction
    for (int i = 0; i < 5; ++i) {
        r_pd1_0 = rand.Uniform(-5, 5);
        r_pd1_1 = rand.Uniform(0, 10);
        r_pd1_2 = rand.Gaus(10, 2);
        r_pd1_3 = rand.Exp(1);

        // Calculate the true x value based on our simulated function
        double true_x = 2.0 * r_pd1_0 * r_pd1_0 +
                        0.5 * std::sin(r_pd1_1) +
                        3.0 * std::log(r_pd1_2 + 1) +
                        -1.0 * r_pd1_3 * r_pd1_0; // Exclude noise for true_x comparison

        // Get the predicted value from the MLP model
        double predicted_x = reader->EvaluateRegression("MLP")[0]; // [0] because it's a single output regression

        std::cout << "Input: (";
        printf("%.2f, %.2f, %.2f, %.2f) -> ", r_pd1_0, r_pd1_1, r_pd1_2, r_pd1_3);
        printf("True X: %.4f, Predicted X: %.4f\n", true_x, predicted_x);
    }

    // Clean up
    delete reader;
    factory->DeleteAllMethods();
    delete factory;
    outputTMVAFile->Close();
    delete outputTMVAFile;
    inputDataFile->Close();
    delete inputDataFile;

    std::cout << "\nTMVASimpleNNRegression finished." << std::endl;
}