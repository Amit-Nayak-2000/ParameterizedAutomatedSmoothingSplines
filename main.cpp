//Automated Smoothing Spline Generator
//Amit Nayak & Dr. David Kopriva

#include "autosspline.hpp"

int main(int argc, char** argv){
    
    std::string controlfilename = argv[1];
    std::string filename, outputfilename;
    std::vector<double> x, y;
    double userdef, eps;
    int evalpoints;
    bool header, periodic;

    readControlFile(controlfilename, filename, outputfilename, userdef, eps, evalpoints, header, periodic);

    readDataset(filename, x, y, header);
    
    if(evalpoints == -1)
        evalpoints = x.size();
    
    PSpline smoothed(userdef, eps, x, y, evalpoints, periodic);

    smoothed.RunBisecant();

    smoothed.writeResults(outputfilename);
    
    return 0;
}