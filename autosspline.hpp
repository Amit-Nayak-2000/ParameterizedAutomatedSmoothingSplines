/**
* @brief Automated data smoothing to a desired length scale using smoothing splines.
* @author Amit Nayak 
* @author Dr. David Kopriva
* @version 1.0
*/
#ifndef AUTOSMOOTHINGSPLINE_HPP
#define AUTOSMOOTHINGSPLINE_HPP

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <limits>


//! @brief Smoothing Spline Class
class SSpline{
    public:
        //! Number of Points
        int npoint; 
        //! Smoothing Parameter
        double p; 
        //! x coordinates
        std::vector<double> x; 
        //! y coordinates
        std::vector<double> y; 
        //! 1/weights array
        std::vector<double> dy; 
        //! Spline Coefficient Multidimensional Array
        std::vector<std::vector<double>> A; 
        //! Coefficient Multidimensional Array for Spline Evaluation
        std::vector<std::vector<double>> CO; 

        //! Smoothing Spline Constructor.
        SSpline(const std::vector<double> &x = std::vector<double>(0), const std::vector<double> &y = std::vector<double>(0), 
                double p = 0.0, const std::vector<double> &weights = std::vector<double>(0)){
            //Check that x & y arrays are consistent.
            assert(x.size() == y.size());

            //if x & y data not given.
            if(x.size() == 0){
                this->npoint = 10;
                this->x = std::vector<double>(npoint);
                this->y = std::vector<double>(npoint);
                this->dy = std::vector<double>(npoint);
                this->A = std::vector<std::vector<double>>(this->npoint, std::vector<double>(4));
                this->CO = std::vector<std::vector<double>>(4, std::vector<double>(this->npoint - 1));

                return;
            }

            //If data is actually given, construct as normal.
            this->npoint = x.size();
            this->p = p;
            this->x = std::vector<double>(npoint);
            this->y = std::vector<double>(npoint);
            this->dy = std::vector<double>(npoint);
            this->A = std::vector<std::vector<double>>(this->npoint, std::vector<double>(4));
            this->CO = std::vector<std::vector<double>>(4, std::vector<double>(this->npoint - 1));
            
            //If weights are not given, set dy = 1, otherwise dy = 1/weights
            if(weights.size() == 0){
                for(int i = 0; i < npoint; i++){
                    this->x[i] = x[i];
                    this->y[i] = y[i];
                    this->dy[i] = 1;
                }
            }
            else{
                for(int i = 0; i < npoint; i++){
                    this->x[i] = x[i];
                    this->y[i] = y[i];
                    this->dy[i] = 1 / weights[i];
                }
            }
            
            smooth();

            transposeForEval();

        }

        //Smoothing Spline Algorithms:

        /** @brief Allocates work arrays, calls setupq, then chol1d. Then finishes constructing coefficient matrix. 
        * For full details on this function, see A Practical Guide to Splines by Carl DeBoor.
        */
        void smooth(){
            //Work arrays
            std::vector<std::vector<double>> V(npoint, std::vector<double>(7));
            std::vector<double> A3(npoint);  
            std::vector<double> A2(npoint);
            std::vector<double> A0(npoint);

            this->setupq(dy, V, A3); 
            this->chol1d(V, A3, A2, A0); 

            for(int i = 0; i < npoint; i++){    
                A[i][0] = y[i] - 6 * (1 - p) * pow(dy[i], 2) * A0[i];
            }

            for(int i = 0; i < npoint; i++){
                A[i][2] = A2[i] * 6 * p;
            }

            for(int i = 0; i < npoint - 1; i++){
                A[i][3] = (A[i+1][2] - A[i][2]) / V[i][3];
                A[i][1] = (A[i+1][0] - A[i][0]) / V[i][3] - (A[i][2] + A[i][3] / 3*V[i][3]) / 2*V[i][3];
            }

        }

        /** @brief Helper function for smooth(). For full details on this function, see A Practical Guide to Splines by Carl DeBoor.
        * @param dx array of (1/weights)
        * @param V Work Array
        * @return QTY Array to be used in coefficient calculation
        */
        void setupq(std::vector<double> &dx, std::vector<std::vector<double>> &V, std::vector<double> &QTY){
            

            for(int i = 0; i < npoint - 1; i++){
                V[i][3] = x[i+1] - x[i];
            }

            for(int i = 1; i < npoint - 1; i++){
                V[i][0] = dx[i-1] / V[i-1][3];
            }

            V[npoint - 1][0] = 0.0;

            for(int i = 1; i < npoint - 1; i++){
                V[i][1] = -dx[i] / V[i][3] - dx[i] / V[i-1][3];
            }

            for(int i = 1; i < npoint - 1; i++){
                V[i][2] = dx[i+1] / V[i][3];
            }

            for(int i = 1; i < npoint - 1; i++){
                V[i][4] = pow(V[i][0], 2) + pow(V[i][1], 2) + pow(V[i][2], 2);
            }

            for(int i = 1; i < npoint - 2; i++){
                V[i][5] = V[i][1] * V[i+1][0] + V[i][2] * V[i+1][1];
            }

            V[npoint - 2][5] = 0.0;

            for(int i = 1; i < npoint - 3; i++){
                V[i][6] = V[i][2] * V[i+2][0];
            }

            V[npoint - 3][6] = 0.0;
            V[npoint - 2][6] = 0.0;

            double prev = (y[1] - y[0]) / V[0][3];
            double diff;

            for(int i = 1; i < npoint - 1; i++){
                diff = (y[i+1] - y[i]) / V[i][3];
                QTY[i] = diff - prev;
                prev = diff;
            }

        }

        /** @brief Helper function for smooth(). For full details on this function, see A Practical Guide to Splines by Carl DeBoor.
        * @param V Work Array
        * @param QTY Array used in calculation for coefficients
        * @return U Array used in calculation for coefficients
        * @return QU Array used in calculation for coefficients
        */
        void chol1d(std::vector<std::vector<double>> &V, std::vector<double> &QTY, std::vector<double> &U, std::vector<double> &QU){

            double six1mp = 6*(1-p);
            double twop = 2*p;

            for(int i = 1; i < npoint - 1; i++){
                V[i][0] = six1mp * V[i][4] + twop*(V[i-1][3] + V[i][3]);
            }

            for(int i = 1; i < npoint - 1; i++){
                V[i][1] = six1mp * V[i][5] + p*V[i][3];
            }

            for(int i = 1; i < npoint - 1; i++){
                V[i][2] = six1mp * V[i][6];
            }

            if(npoint < 4){
                U[0] = 0.0;
                U[1] = QTY[1] / V[1][0];
                U[2] = 0.0;
            }
            else{

                double ratio;

                for(int i = 1; i < npoint - 2; i++){
                    ratio = V[i][1] / V[i][0];
                    V[i+1][0] = V[i+1][0] - ratio * V[i][1];
                    V[i+1][1] = V[i+1][1] - ratio * V[i][2];
                    V[i][1] = ratio;
                    ratio = V[i][2] / V[i][0];
                    V[i+2][0] = V[i+2][0] - ratio*V[i][2];
                    V[i][2] = ratio;
                }

                U[0] = 0.0;
                V[0][2] = 0.0; 
                U[1] = QTY[1];

                
                for(int i = 1; i < npoint - 2; i++){
                    U[i+1] = QTY[i+1] - V[i][1] * U[i] - V[i-1][2] * U[i-1];
                }

                U[npoint - 1] = 0.0;
                U[npoint - 2] = U[npoint - 2] / V[npoint - 2][0];

                for(int i = npoint - 3; i >= 1; i--){
                    U[i] = U[i] / V[i][0] - U[i+1]*V[i][1] - U[i+2]*V[i][2];
                }
            }

            double prev = 0.0;

            for(int i = 1; i < npoint; i++){
                QU[i] = (U[i] - U[i-1]) / V[i-1][3];
                QU[i-1] = QU[i] - prev;
                prev = QU[i];
            }

            QU[npoint - 1] = -QU[npoint - 1];
            

        }

        /** @brief Function for finding index for spline evaluation (Recursive). 
        * @param breaks Knots for spline
        * @param x Evaluation point
        * @param lower lower bound index
        * @param upper upper bound index
        * @return Array index
        */
        int findbreakindex(const std::vector<double> &breaks, const double &x, int lower, int upper){
            assert(lower < upper);

            //if the last (largest) element of breaks is smaller than or equal to x.
            if(breaks[upper] <= x){
                return upper;
            }
            
            int middle = lower + (upper - lower) * 0.5;

            if(breaks[middle] <= x && breaks[middle + 1] > x){
                return middle;
            }
            //if the middle element is smaller than x, the appropriate index lies between middle + 1 and upper
            else if(breaks[middle + 1] < x){
                return findbreakindex(breaks, x, middle + 1, upper);
            }
            //middle element is greater than x, appropriate index lies between lower and middle -1
            else{
                return findbreakindex(breaks, x, lower, middle);
            }
            
            //if could not find index
            return -1;

        }

        /** @brief Evaluate spline (piecewise polynomial). For full details on this function, see A Practical Guide to Splines by Carl DeBoor.
        * @param breaks Knots of spline
        * @param coeff Spline coefficients
        * @param eval Evaluation point.
        * @param jderiv Derivative order
        * @param K Spline order + 1 (this will be 4 for cubic spline)
        * @return Evaluated spline value
        */
        double ppvalu(const std::vector<double> &breaks, std::vector<std::vector<double>> &coeff, double eval, int jderiv, int K){
            int L = breaks.size() - 1;
            double result = 0.0;

            double fmmjdr = K - jderiv;

            if(fmmjdr <= 0){
                return result;
            }

            //returns c style index
            int i = findbreakindex(breaks, eval, 0, L);

            assert(i != -1);

            double h = eval - breaks[i];

            int m = K - 1;
        

            bool flag = true;
            

            while(flag){
                result = ((result / fmmjdr) * h) + coeff[m][i];
                m--;
                fmmjdr -= 1;

                if(fmmjdr <= 0){
                    flag = false;
                    break;
                }
            }


            return result;

        }

        /** @brief Transposes A and decrements a dimension appropriate for ppvalu.
        */
        void transposeForEval(){
            for(int i = 0; i < npoint-1; i++){
                CO[0][i] = A[i][0];
                CO[1][i] = A[i][1];
                CO[2][i] = A[i][2];
                CO[3][i] = A[i][3];
            }
        }

        /** @brief Change smoothing param and recalculate smoothing spline.
        * @param pnew New smoothing param
        */      
        void changeSmooth(double pnew){
            this->p = pnew;
            smooth();
            transposeForEval();
        }

        /** @brief Evaluate smoothing spline.
        * @param xeval Array of evaluation points
        * @return yeval Array of evaluated points
        */
        void evalSpline(const std::vector<double> &xeval, std::vector<double> &yeval){
            for(int i = 0; i < xeval.size(); i++){
                yeval[i] = ppvalu(x, CO, xeval[i], 0, 4);
            }
        }

};


//Forward Declarations of Helper Functions
std::vector<double> parametrize(const std::vector<double> &x, const std::vector<double> &y);
std::vector<double> linarray(double min, double max, int npoint);
void numfirstder(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &yprime);
void num2der(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &ydoubleprime);
void radiusofcurve(const std::vector<double> &t, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &result);
void readDataset(std::string filename, std::vector<double> &x, std::vector<double> &y);
void readControlFile(const std::string &controlfile, std::string &filename, std::string &outputfilename, double &userdef, double &tolerance, int &evalpoints, bool &header, bool &periodic);


//Helper Functions:

/** @brief Parametrizing Data (according to arclength).
* @param x Coordinate point array
* @param y Coordinate point array
* @return Parameterized array
*/
std::vector<double> parametrize(const std::vector<double> &x, const std::vector<double> &y){
    int npoint = x.size();
    double totalarclength = 0;

    for(int i = 1; i < npoint; i++){
        totalarclength += pow((pow((x[i] - x[i-1]), 2) + pow((y[i] - y[i-1]), 2)), 0.5);
    }
    
    std::vector<double> result(npoint);
    result[0] = 0;

    for(int i = 1; i < npoint; i++){
        result[i] = (result[i-1] + pow((pow((x[i] - x[i-1]), 2) + pow((y[i] - y[i-1]), 2)), 0.5));
    }

    for(int i = 0; i < npoint; i++){
        result[i] = result[i] / totalarclength;
    }

    return result;
}


/** @brief Generate a linearly spaced array.
* @param npoint Number of points
* @param min Lowest value in desired array
* @param max Maximmum value in desired array
* @return Linearly spaced array with desired number of points between min & max
*/
std::vector<double> linarray(double min, double max, int npoint){
    std::vector<double> result(npoint);

    double spacing = (max-min)/double(npoint-1);

    for(int i = 0; i < npoint; i++) {
        result[i] = min + i * spacing;
    }

    return result;
}


/** @brief Numerical first derivative.
* @param x Data array to calculate derivative of y WRT x
* @param y Data array to calculate derivative of y WRT x
* @return yprime Array that contains dy/dx
*/
void numfirstder(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &yprime){
    int N = x.size();

    for(int i = 0; i < N; i++){
        //beginning edge case
        //forward difference approximation
        if(i == 0){ 
            yprime[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
        }
        //ending edge case
        //backward difference approximation
        else if (i == (N - 1)){
            yprime[i] = (y[i] - y[i-1]) / (x[i] - x[i-1]);
        }
        //centered difference approximation
        else{
            yprime[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1]);
        }
    }

} 


/** @brief Numerical second derivative.
* @param x Data array to calculate 2nd derivative of y WRT x
* @param y Data array to calculate 2nd derivative of y WRT x
* @return ydoubleprime Array that contains d2y/dx2
*/
void num2der(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &ydoubleprime){
    int N = x.size();

    for(int i = 0; i < N; i++){
        //beginning edge case
        //forward difference approximation
        if(i == 0){
            ydoubleprime[i] = (2*y[i] - 5*y[i+1] + 4*y[i+2] - y[i+3]) / pow((x[i+1] - x[i]), 2);
        }
        //ending edge case
        //backward difference approximation
        else if (i == (N - 1)){
            ydoubleprime[i] = (2*y[i] - 5*y[i-1] + 4*y[i-2] - y[i-3]) / pow((x[i] - x[i-1]), 2);
        }
        //centered difference approximation
        else{
            ydoubleprime[i] = (y[i+1] - 2*y[i] + y[i-1]) / pow((x[i+1] - x[i]), 2);
        }
    }

} 


/** @brief Parametrized radius of curvature.
* @param x Coordinate point array
* @param y Coordinate point array
* @param t Parameterized array
* @return Result array that contains parameterized radius of curvature
*/
void radiusofcurve(const std::vector<double> &t, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &result){
    int N = t.size();

    std::vector<double> xprime(N);
    std::vector<double> yprime(N);
    std::vector<double> xdoubleprime(N);
    std::vector<double> ydoubleprime(N);

    numfirstder(t, x, xprime);
    numfirstder(t, y, yprime);
    num2der(t, x, xdoubleprime);
    num2der(t, y, ydoubleprime);


    double denominator, zerolim;
    zerolim = 10*std::numeric_limits<double>::epsilon();

    for(int i = 0; i < N; i++){
        denominator = std::abs(xprime[i]*ydoubleprime[i] - yprime[i]*xdoubleprime[i]);

        //if denominator is close to 0, set radius of curvature to max value.
        if(denominator <= zerolim){
            result[i] = std::numeric_limits<double>::max();
            continue;
        }

        result[i] = std::abs(pow((pow(xprime[i], 2) + pow(yprime[i], 2)),1.5)) / denominator;
    }
    
}


/** @brief Reading in data.
* @note Assuming array inputs are empty.
* @param filename Filename of data .csv file
* @param header boolean to keep track if there is a header
* @return x x-coordinate points
* @return y y-coordinate points
*/
void readDataset(std::string filename, std::vector<double> &x, std::vector<double> &y, bool header){
    std::ifstream DataFile(filename);
    std::string line;
    std::string tempstring;
    
    int i = 0;
    double xparam;
    double yparam;

    while(DataFile.good()){
        //skip header
        if(header){
            //pass by header line
            getline(DataFile, line);
            std::stringstream parsestring(line);
            header = false;
        }

        //read in line
        getline(DataFile, line);
        std::stringstream parsestring(line);

        //read x coord
        getline(parsestring, tempstring, ',');
        xparam = std::stod(tempstring);
        
        //read y coord
        getline(parsestring, tempstring, ',');
        yparam = std::stod(tempstring);

        if(i > 0 && xparam == x[i-1] && yparam == y[i-1]){
            continue;
        }

        x.push_back(xparam);
        y.push_back(yparam);
        i++;

        
        
    }

    DataFile.close();
}


/** @brief Reading in control file.
* @param controlfile Filename of control file
* @return filename Filename of .csv data file
* @return outputfilename Filename of output file
* @return userdef User defined radius of curvature
* @return tolerance User defined tolerance of radius of curvature
* @return evalpoints Number of evaluation points for smoothed data
* @return header If there is a header in the input file
* @return periodic If the dataset is periodic
*/
void readControlFile(const std::string &controlfile, std::string &filename, std::string &outputfilename, double &userdef, double &tolerance, int &evalpoints, bool &isHeader, bool &periodic){

    // std::string controlfile = "userdef.control";
    std::ifstream DataFile(controlfile);
    //bulk string to read line.
    std::string line;
    //string to hold key
    std::string key;
    //string to hold value
    std::string value;
    //map to hold key value pairs.
    std::map<std::string, std::string> ParseMap;
    
    bool header = false;
    while(DataFile.good()){
        if(header == false){
            //pass by header line
            getline(DataFile, line);
            std::stringstream parsestring(line);
            header = true;
        }
        else{
            //read in line
            getline(DataFile, line);
            std::stringstream parsestring(line);
            //read text before =
            getline(parsestring, key, '=');
            //read text after =
            getline(parsestring, value, '=');
        }

        ParseMap[key] = value;
    }

    //Fill in mandatory params:
    filename = ParseMap.find("input_filename ")->second;
    outputfilename = ParseMap.find("output_filename ")->second;

    //remove any whitespaces from filenames:
    filename.erase(remove_if(filename.begin(), filename.end(), isspace), filename.end());
    outputfilename.erase(remove_if(outputfilename.begin(), outputfilename.end(), isspace), outputfilename.end());

    //Fill in radius of curvature, tolerance, and if there is a header:
    userdef = std::stod(ParseMap.find("radius_of_curvature ")->second);
    tolerance = std::stod(ParseMap.find("tolerance ")->second);
    isHeader = std::stoi(ParseMap.find("header_in_input_file? ")->second);

    //Remove any possible whitespaces from optional params as its likely that people may keep a space:
    //Evaluation points: (Optional Param)
    auto evaliter = ParseMap.find("number_of_evaluation_points ");
    evaliter->second.erase(remove_if(evaliter->second.begin(), evaliter->second.end(), isspace), evaliter->second.end());

    if(evaliter->second.empty() == true){
        evalpoints = -1;
    }
    else{
        evalpoints = std::stoi(evaliter->second);
    }

    //Periodic Datatset: (Optional Param)
    auto piter = ParseMap.find("periodic_dataset? ");
    piter->second.erase(remove_if(piter->second.begin(), piter->second.end(), isspace), piter->second.end());

    if(piter->second.empty() == true){
        periodic = 0;
    }
    else{
        periodic = std::stoi(piter->second);
    }

    DataFile.close();
}


//! @brief Paramterized Spline Wrapper Class
class PSpline{
    public:
        //! Smoothing spline for x
        SSpline Xspline;
        //! Smoothing spline for y
        SSpline Yspline;

        //! Upperbound for bisection
        double u;
        //! Lowerbound for bisection
        double l; 
        //! Userdefined radius of curvature
        double userdef; 
        //! Iteration tolerance
        double eps; 
        //! Final smoothing parameter
        double fp; 
        //! Final radius of curvature
        double result;
        //! Boolean to keep track if dataset is periodic
        bool periodic;
        //! # of evaluation points for smoothed data
        int evalsize;

        //! Parameter array
        std::vector<double> t;
        //! Evaluation array in parameterized space
        std::vector<double> tvals;
        //! xdata
        std::vector<double> x;
        //! ydata
        std::vector<double> y;
        //! smoothed x
        std::vector<double> xnew;
        //! smoothed y
        std::vector<double> ynew;

        //! @brief Parameterized Smoothing Spline Constructor
        PSpline(double userdef, double eps, std::vector<double> x, std::vector<double> y, int evalsize, bool periodic = false, double u = 0.9999999, double l = 0.0000001){

            this->u = u;
            this->l = l;
            this->userdef = userdef;
            this->eps = eps;
            this->evalsize = evalsize;
            this->periodic = periodic;
            

            this->x = std::vector<double>(x.size());
            this->y = std::vector<double>(y.size());
            this->xnew = std::vector<double>(evalsize);
            this->ynew = std::vector<double>(evalsize);

            for(int i = 0; i < x.size(); i++){
                this->x[i] = x[i];
                this->y[i] = y[i];
            }

            this->t = parametrize(this->x, this->y);
            this->tvals = linarray(0.0, 1.0, evalsize);
        }

        /** @brief Bisection-Secant Algorithm Hybrid
        * @details Construct an initial guess for the smoothing parameter (0.5), increase tolerance, and evaluate smoothing splines.
        * Iteration: Calculate radius of curvature at every point and obtain the minimum, and check if it satsifies the users requirement.
        * If the min radius of curvature is satisfactory, call secant function to finish smoothing the data.
        * If it is not satisfactory, adjust the bounds for bisection and keep iterating.
        */
        void RunBisecant(){
            //middle point for bisection
            double p = (u + l) / 2;

            eps *= 100;

            int evalsize = tvals.size();

            //construct inital splines for x and y
            Xspline = SSpline(t, x, p);
            Yspline = SSpline(t, y, p);

            //evaluate initial spline points
            Xspline.evalSpline(tvals, xnew);
            Yspline.evalSpline(tvals, ynew);

            //vector to contain radius of curvature
            std::vector<double> radius(evalsize);

            double M;

            //Bisection begins here
            for(int i = 0; i < 150; i++){
                //get radius of curvature
                radiusofcurve(tvals, xnew, ynew, radius);

                //get minimum radius of curvature
                M = *std::min_element(radius.begin(), radius.end());

                //check if obtained min radius of curvature satisfies the user requirement.
                if(std::abs(M - userdef) < eps){
                    std::cout << "Bisection finished with " << i+1 << " iterations." << std::endl;
                    std::cout << std::endl;
                    fp = p;
                    result = M;
                    
                    //call secant to finish iteration to find optimal smoothing parameter and smooth data.
                    secant();
                    
                    return;
                }
                else if (M < userdef){
                    //recalculate smoothing param
                    u = p;
                    p += 0.5*(l - p);

                    //recalculate splines with new smoothing param
                    Xspline.changeSmooth(p);
                    Yspline.changeSmooth(p);
                    Xspline.evalSpline(tvals, xnew);
                    Yspline.evalSpline(tvals, ynew);
                }
                else if (M > userdef){
                    //recalculate smoothing param
                    l = p;
                    p += 0.5*(u - p);

                    //recalculate splines with new smoothing param
                    Xspline.changeSmooth(p);
                    Yspline.changeSmooth(p);
                    Xspline.evalSpline(tvals, xnew);
                    Yspline.evalSpline(tvals, ynew);
                }

                fp = p;
                result = M;

            }
        

        }

        /** @brief Secant portion of Bisection-Secant Algorithm Hybrid
        * @details Decrease tolerance back to user defined tolerance, and proceed to use the secant method to obtain desired smoothness.
        * Iteration: compute new potential smoothing parameter according to secant method and check if the smoothing parameter
        * produces the desired minimum radius of curvature. If desired smoothness is reached, iteration stops. Otherwise, 
        * keep iterating till desired smoothness is reached.
        */
        void secant(){
            double a = l;
            double b = u;
            eps /= 100;

            //construct inital splines for a
            Xspline.changeSmooth(a);
            Yspline.changeSmooth(a);


            //evaluate initial spline points
            Xspline.evalSpline(tvals, xnew);
            Yspline.evalSpline(tvals, ynew);
            
            double fa = secantFunc();

            //change spline param to B
            Xspline.changeSmooth(b);
            Yspline.changeSmooth(b);
            Xspline.evalSpline(tvals, xnew);
            Yspline.evalSpline(tvals, ynew);

            double fb = secantFunc();

            double d;

            //Secant method begins here
            for(int i = 0; i <= 150; i++){
                d = ((b - a) / (fb - fa));
                b = a;
                fb = fa;
                d = d*fa;
                a = a - d;

                //evaluate new a
                Xspline.changeSmooth(a);
                Yspline.changeSmooth(a);
                Xspline.evalSpline(tvals, xnew);
                Yspline.evalSpline(tvals, ynew);

                //call helper function which provides the difference of the obtained and user defined min radius of curvature.
                fa = secantFunc();

                if(std::abs(fa) <= eps){
                    std::cout << "Secant finished with " << i+1 << " iterations." << std::endl;
                    std::cout << std::endl;
                    fp = a;
                    result = fa + userdef;
                    std::cout << "Final radius of curvature is " << result << "." << std::endl;
                    return;
                }
                
            }

            fp = a;
            result = fa + userdef;

        }

        /** @brief Helper function called by secant().
        * @return Absolute difference between obtained minimum radius of curvature & user defined radius of curvature 
        * @details Calculates radius of curvature at every point. 
        * After, find the minimum of that array and subtract the user defined radius of curvature.
        * return the absolute value of the difference.
        */
        double secantFunc(){
            std::vector<double> radius(evalsize);

            radiusofcurve(tvals, xnew, ynew, radius);
            
            return (*std::min_element(radius.begin(), radius.end()) - userdef);
        }

        /** @brief Writing results to output file (csv format).
        * @param OutFileName string that contains the name of desired output file.
        * @details Creates an output file and writes the smoothed data to it. 
        */
        void writeResults(std::string OutFileName){
            
            assert(xnew.size() == ynew.size());

            std::ofstream DataFile(OutFileName);

            DataFile << "x-smoothed, y-smoothed\n";

            for(int i = 0; i < xnew.size() - 1; i++){
                DataFile << xnew[i] << ", " << ynew[i] << "\n";
            }
            
            //If data is periodic, average the first and last point
            if(periodic){
                DataFile << 0.5*(xnew[0] + xnew[xnew.size() - 2]) << ", " << 0.5*(ynew[0] + ynew[ynew.size() - 2]) << "\n";
            }

            

            DataFile.close();
        }


};

#endif