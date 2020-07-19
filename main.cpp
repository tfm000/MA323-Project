//Main.cpp
// *** Main file for MA323 Project Q3 & Q4 -- Candidate Number: 40946 *** //


// Header Files
#include "Normals.hpp"
#include "Question3.hpp"
#include "Question4.hpp"
#include "RandomNumberGenerators.hpp"
using namespace std;

int main(){
    srand(time(NULL));

    // Setting Global Variables //
    int n = 10000, m = 1000;                                                     // Simulation Variables.
    double C = 15, H = 10, S0 = 10, r = 0.01, T = 1, sigma = 0.2, gamma = 0.75; // Option Variables.


    // *** QUESTION 3 CODE *** //


    // Calculating Analytical Black-Scholes Solution //
    double ANSolution = analytical(C,H,S0,r,T,sigma);

    // Calculating Monte-Carlo estimate and Confidence Intervals //
    MonteCarlo MC(C,H,S0,r,T,sigma,n);                          // Our Monte-Carlo Object.
    double MCEstimate = MC.getEstimate(), MCvar = MC.myanalyticalVariance;

    // Calculating 95% CI
    double width95 = MC.ConfidenceInterval(0.05); // epsilon = 0.05
    double left95 = MC.myleftCI, right95 = MC.myrightCI;

    // Calculating 99% CI
    double width99 = MC.ConfidenceInterval(0.01); // epsilon = 0.01
    double left99 = MC.myleftCI, right99 = MC.myrightCI;

    // Calculating Control-Variates estimate //
    ControlVariates CV(C,H,S0,r,T,sigma,n);
    double CVEstimate = CV.getEstimate(), CVvar = CV.myvariance;

    // QUESTION 3 OUTPUT STREAM //
    cout<<"Question 3\n\n";
    cout<<"Pricing Method\t\tPrice\t\tVariance\t95% CI\t\t\t99% CI";
    cout<<"\n\nAnalytical\t\t"<<ANSolution<<"\t\t"<<0.0;
    cout<<"\n\nMonte-Carlo*\t\t"<<MCEstimate<<"\t\t"<<MCvar<<"\t("<<left95<<","<<right95<<")\t("<<left99<<","<<right99<<")";
    cout<<"\n\nControl Variates*\t"<<CVEstimate<<"\t\t"<<CVvar;
    cout<<"\n\n95% Monte-Carlo CI width: "<<width95;
    cout<<"\n99% Monte-Carlo CI width: "<<width99;
    cout<<"\nControl Variates variance reduction (wrt Monte-Carlo Estimation) of: "<<(1.0 - CVvar/MCvar)*100<<"%";
    cout<<"\n\n*Using "<<n<<" iterations";

    // *** QUESTION 4 CODE *** //


    cout<<"\n\n\n\nQuestion 4 (May take a while to compute)\n\n";

    // Calculating Monte-Carlo estimate and Confidence Intervals using both CEV asset pricing schemes //
    CEVMonteCarlo CEV_MC(C,H,S0,r,T,sigma,n,gamma,m); // Our CEV Monte-Carlo Object.
    double MilsteinSchemeMCEstimate = CEV_MC.getEstimate();
    double EulerSchemeMCEstimate = CEV_MC.myEulerEstimate;
    double EulerVar = CEV_MC.myEulerVariance, MilsteinVar = CEV_MC.myMilsteinVariance;

    // Calculating 95% CIs
    double MilsteinWidth95 = CEV_MC.ConfidenceInterval(0.05), EulerWidth95 = CEV_MC.EulerCIWidth;
    double EulerLeft95 = CEV_MC.EulerLeftCI, EulerRight95 = CEV_MC.EulerRightCI;
    double MilsteinLeft95 = CEV_MC.MilsteinLeftCI, MilsteinRight95 = CEV_MC.MilsteinRightCI;

    // Calculating 99% CIs
    double MilsteinWidth99 = CEV_MC.ConfidenceInterval(0.01), EulerWidth99 = CEV_MC.EulerCIWidth;
    double EulerLeft99 = CEV_MC.EulerLeftCI, EulerRight99 = CEV_MC.EulerRightCI;
    double MilsteinLeft99 = CEV_MC.MilsteinLeftCI, MilsteinRight99 = CEV_MC.MilsteinRightCI;

    // Calculating CEV Milstein Control-Variates estimate //
    CEVControlVariates CEV_CV(C,H,S0,r,T,sigma,n,gamma,m); // Our CEV Monte-Carlo Object.
    double bstarhatCEV = CEV_CV.bstarhat(); // our b*hat i.e. the optimal choice of b.
    double b_SubOptimal = 1.0; // a sub-optimal choice for b.
    double CEVMilstein_bstarhat = CEV_CV.getEstimate(bstarhatCEV);
    double optimalVar = CEV_CV.myvariance;
    double CEVMilstein_b_SubOptimal = CEV_CV.getEstimate(b_SubOptimal);
    double subOptimalVar = CEV_CV.myvariance;

    // QUESTION 4 OUTPUT STREAM //
    cout<<"Pricing Method\t\t\t\tPrice\t\tVariance\t95% CI\t\t\t99% CI";
    cout<<"\n\nMonte-Carlo**\t\t\t\t"<<EulerSchemeMCEstimate<<"\t\t"<<EulerVar<<"\t("<<EulerLeft95<<","<<EulerRight95<<")\t("<<EulerLeft99<<","<<EulerRight99<<")";
    cout<<"\n(CEV & Euler approx)";
    cout<<"\n\nMonte-Carlo**\t\t\t\t"<<MilsteinSchemeMCEstimate<<"\t\t"<<MilsteinVar<<"\t("<<MilsteinLeft95<<","<<MilsteinRight95<<")\t("<<MilsteinLeft99<<","<<MilsteinRight99<<")";
    cout<<"\n(CEV & Milstein approx)";

    cout<<"\n\nOptimal b Control Variates**\t\t"<<CEVMilstein_bstarhat<<"\t\t"<<optimalVar;
    cout<<"\n(CEV & Milstein approx)";
    cout<<"\n\nSub-optimal b Control Variates**\t"<<CEVMilstein_b_SubOptimal<<"\t\t"<<subOptimalVar;
    cout<<"\n(CEV & Milstein approx)";


    cout<<"\n\n95% Monte-Carlo (CEV Euler approx) CI width: "<<EulerWidth95;
    cout<<"\n99% Monte-Carlo (CEV Euler approx) CI width: "<<EulerWidth99;
    cout<<"\n95% Monte-Carlo (CEV Milstein approx) CI width: "<<MilsteinWidth95;
    cout<<"\n99% Monte-Carlo (CEV Milstein approx) CI width: "<<MilsteinWidth99;
    cout<<"\n\nOptimal b = b*hat = "<<bstarhatCEV;
    cout<<"\nSub-optimal b = "<<b_SubOptimal;
    cout<<"\nControl Variates Optimal b variance reduction (wrt CEV Monte-Carlo Estimation) of: "<<(1.0 - optimalVar/MilsteinVar)*100<<"%";
    cout<<"\nControl Variates Sub-optimal b variance reduction (wrt CEV Monte-Carlo Estimation) of: "<<(1.0 - subOptimalVar/MilsteinVar)*100<<"%";
    cout<<"\nControl Variates Optimal b variance reduction wrt Sub-optimal b of: "<<(1.0 - optimalVar/subOptimalVar)*100<<"%";

    cout<<"\n\n**Using "<<n<<" iterations and "<<m<<" grid partitions\n\n";
    return 0;
}


